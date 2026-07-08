#!/usr/bin/env python3
"""Self-contained primate mtDNA coordinate liftover workflow.

The workflow prepares raw chrM FASTAs, rotates species and human references,
aligns each species to human, builds coordinate maps, lifts raw VCF/COV files to
canonical human chrM coordinates, and writes liftover QC reports.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import gzip
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DNA = set("ACGTNacgtn")
OUTDIRS = [
    "prepared_fastas",
    "rotated_fastas",
    "alignments",
    "maps",
    "vcf_lifted_raw",
    "cov_lifted",
    "reports",
]


@dataclass
class FastaRecord:
    name: str
    seq: str


@dataclass
class Sample:
    name: str
    species_fasta: Path
    vcf: Path
    cov: Path
    species_chrom: str = "chrM"
    rotate_anchor: Optional[int] = None
    target_sequence: Optional[str] = None
    species: Optional[str] = None


@dataclass
class SampleStats:
    sample: str
    species_bases: int = 0
    mapped_positions: int = 0
    unmapped_positions: int = 0
    gaps: int = 0
    ambiguous_mappings: int = 0
    vcf_variants_lifted: int = 0
    vcf_variants_failed_liftover: int = 0
    cov_positions_lifted: int = 0
    cov_positions_failed_liftover: int = 0
    ref_match_count: int = 0
    ref_mismatch_count: int = 0
    notes: List[str] = field(default_factory=list)


def read_fasta(path: Path, target: Optional[str] = None) -> FastaRecord:
    records: List[FastaRecord] = []
    name: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(FastaRecord(name, "".join(seq_chunks).upper()))
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        records.append(FastaRecord(name, "".join(seq_chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    if target:
        matches = [r for r in records if r.name == target]
        if len(matches) != 1:
            raise ValueError(f"Expected one target sequence {target!r} in {path}, found {len(matches)}")
        rec = matches[0]
    elif len(records) == 1:
        rec = records[0]
    else:
        mt = [r for r in records if r.name.lower() in {"chrm", "mt", "m", "mitochondrion", "mitochondrial"}]
        if len(mt) != 1:
            raise ValueError(f"{path} contains multiple records; set target_sequence/human_target_sequence")
        rec = mt[0]
    bad = set(rec.seq) - DNA
    if bad:
        raise ValueError(f"Unexpected FASTA bases in {path}: {''.join(sorted(bad))}")
    return rec


def write_fasta(path: Path, name: str, seq: str, width: int = 80) -> None:
    with path.open("w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i : i + width] + "\n")


def rotate_sequence(seq: str, anchor_1based: int) -> str:
    if anchor_1based < 1 or anchor_1based > len(seq):
        raise ValueError(f"Rotation anchor {anchor_1based} outside 1..{len(seq)}")
    i = anchor_1based - 1
    return seq[i:] + seq[:i]


def rotated_to_original(pos: int, anchor: int, length: int) -> int:
    return ((pos + anchor - 2) % length) + 1


def original_to_rotated(pos: int, anchor: int, length: int) -> int:
    return ((pos - anchor) % length) + 1


def restore_human_pos(rotated_pos: int, human_len: int, restore_offset: int) -> int:
    return ((rotated_pos + restore_offset - 1) % human_len) + 1


def infer_anchor(species_seq: str, human_seq: str) -> int:
    """Infer a species rotation anchor from a shared k-mer; fall back to 1."""
    upper_human = human_seq.upper()
    upper_species = species_seq.upper()
    for k in (31, 25, 21, 15):
        for start in range(0, max(1, len(upper_human) - k + 1), max(1, k // 3)):
            kmer = upper_human[start : start + k]
            if len(kmer) == k and "N" not in kmer:
                hit = upper_species.find(kmer)
                if hit >= 0:
                    return hit + 1
    return 1


def write_simple_alignment(species_fa: Path, human_fa: Path, out_fa: Path) -> None:
    """Write a deterministic ungapped two-sequence alignment for smoke tests."""
    s = read_fasta(species_fa)
    h = read_fasta(human_fa)
    max_len = max(len(s.seq), len(h.seq))
    write_fasta(out_fa, s.name, s.seq.ljust(max_len, "-"))
    with out_fa.open("a") as out:
        out.write(f">{h.name}\n")
        padded = h.seq.ljust(max_len, "-")
        for i in range(0, len(padded), 80):
            out.write(padded[i : i + 80] + "\n")


def run_alignment(species_fa: Path, human_fa: Path, out_fa: Path, cfg: configparser.ConfigParser) -> None:
    aligner = cfg.get("alignment", "aligner", fallback="mafft")
    opts = cfg.get("alignment", "aligner_options", fallback="--auto --quiet").split()
    fallback = cfg.getboolean("alignment", "allow_simple_alignment_fallback", fallback=True)
    use_conda_env = cfg.getboolean("alignment", "use_conda_env", fallback=True)
    module_load = cfg.get("alignment", "module_load", fallback="miniconda/24.11.3").strip()
    conda_env = cfg.get("alignment", "conda_env", fallback="mafft_env").strip()
    tmp = out_fa.with_suffix(".input.fa")
    tmp.write_text(species_fa.read_text() + human_fa.read_text())
    try:
        if use_conda_env:
            quoted_cmd = " ".join([shlex.quote(aligner), *[shlex.quote(o) for o in opts], shlex.quote(str(tmp))])
            shell_lines = ["source /etc/profile >/dev/null 2>&1 || true"]
            if module_load:
                shell_lines.append(f"module load {shlex.quote(module_load)}")
            if conda_env:
                shell_lines.append('source "$(conda info --base)/etc/profile.d/conda.sh"')
                shell_lines.append(f"conda activate {shlex.quote(conda_env)}")
            shell_lines.append(quoted_cmd)
            with out_fa.open("w") as out:
                subprocess.run("\n".join(shell_lines), check=True, stdout=out, shell=True, executable="/bin/bash")
            return
        if shutil.which(aligner):
            with out_fa.open("w") as out:
                subprocess.run([aligner, *opts, str(tmp)], check=True, stdout=out)
            return
        if not fallback:
            raise RuntimeError(f"Aligner {aligner!r} not found")
    except (subprocess.CalledProcessError, RuntimeError):
        if not fallback:
            raise
    finally:
        tmp.unlink(missing_ok=True)
    write_simple_alignment(species_fa, human_fa, out_fa)


def load_alignment(path: Path) -> Tuple[FastaRecord, FastaRecord]:
    recs = []
    name = None
    chunks: List[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append(FastaRecord(name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        elif line:
            chunks.append(line.strip())
    if name is not None:
        recs.append(FastaRecord(name, "".join(chunks).upper()))
    if len(recs) != 2:
        raise ValueError(f"Expected two aligned FASTA records in {path}, found {len(recs)}")
    if len(recs[0].seq) != len(recs[1].seq):
        raise ValueError(f"Aligned sequences have different lengths in {path}")
    return recs[0], recs[1]


def build_map(sample: Sample, aln: Path, map_path: Path, species_anchor: int, human_anchor: int, human_len: int, restore_offset: int) -> Tuple[Dict[int, dict], SampleStats]:
    species, human = load_alignment(aln)
    stats = SampleStats(sample.name)
    species_rot = 0
    human_rot = 0
    pos_map: Dict[int, dict] = {}
    rows: List[dict] = []
    for s_base, h_base in zip(species.seq, human.seq):
        if s_base != "-":
            species_rot += 1
        if h_base != "-":
            human_rot += 1
        if s_base == "-" or h_base == "-":
            stats.gaps += 1
        if s_base == "-":
            continue
        stats.species_bases += 1
        original = rotated_to_original(species_rot, species_anchor, len(species.seq.replace("-", "")))
        human_canonical = ""
        status = "unmapped_gap" if h_base == "-" else "mapped"
        if h_base != "-":
            human_canonical_i = restore_human_pos(human_rot, human_len, restore_offset)
            if not (1 <= human_canonical_i <= human_len):
                raise ValueError(f"Restored human position out of bounds: {human_canonical_i}")
            human_canonical = str(human_canonical_i)
            stats.mapped_positions += 1
        else:
            stats.unmapped_positions += 1
        row = {
            "sample": sample.name,
            "species_chrom": sample.species_chrom,
            "species_pos_original": str(original),
            "species_pos_rotated": str(species_rot),
            "human_pos_rotated": str(human_rot) if h_base != "-" else "",
            "human_pos_canonical": human_canonical,
            "species_base": s_base,
            "human_base": h_base,
            "map_status": status,
        }
        rows.append(row)
        pos_map[original] = row
    with map_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()) if rows else ["sample"] , delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return pos_map, stats


def require_existing_file(path: Path, label: str) -> Path:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def require_gzipped_vcf(path: Path) -> Path:
    if not path.name.endswith(".vcf.gz"):
        raise ValueError(f"Input VCF must end with .vcf.gz: {path}")
    return require_existing_file(path, "VCF")


def find_species_fasta(species: str, cfg: configparser.ConfigParser) -> Path:
    fasta_dir = Path(cfg.get("paths", "species_fasta_dir", fallback="").strip())
    if not fasta_dir:
        raise ValueError("sample_ref_file uses a species column; configure [paths] species_fasta_dir")
    extensions = [e.strip() for e in cfg.get("paths", "species_fasta_extensions", fallback=".fa,.fasta,.fna").split(",") if e.strip()]
    candidates = [fasta_dir / f"{species}{ext}" for ext in extensions]
    existing = [p for p in candidates if p.exists() and p.is_file()]
    if len(existing) == 1:
        return existing[0]
    if not existing:
        tried = ", ".join(str(p) for p in candidates)
        raise FileNotFoundError(f"No species FASTA found for species {species!r}; tried: {tried}")
    raise ValueError(f"Multiple species FASTAs found for species {species!r}: {existing}")


def find_sample_file(sample: str, cfg: configparser.ConfigParser, dir_key: str, pattern_key: str, label: str) -> Path:
    base_dir = Path(cfg.get("paths", dir_key, fallback="").strip())
    if not base_dir:
        raise ValueError(f"Configure [paths] {dir_key} to resolve {label} files from sample names")
    patterns = [p.strip() for p in cfg.get("paths", pattern_key, fallback=f"{{sample}}*").split(",") if p.strip()]
    matches: List[Path] = []
    for pattern in patterns:
        matches.extend(sorted(base_dir.glob(pattern.format(sample=sample))))
    unique = sorted({p for p in matches if p.is_file()})
    if len(unique) == 1:
        return unique[0]
    if not unique:
        rendered = ", ".join(str(base_dir / p.format(sample=sample)) for p in patterns)
        raise FileNotFoundError(f"No {label} file found for sample {sample!r}; tried: {rendered}")
    raise ValueError(f"Multiple {label} files found for sample {sample!r}: {unique}")


def _looks_like_header(fields: Sequence[str]) -> bool:
    normalized = {f.strip().lower() for f in fields}
    return "sample" in normalized or "species" in normalized or "species_fasta" in normalized


def iter_sample_ref_rows(path: Path) -> Iterable[dict]:
    """Yield sample rows from headered or legacy headerless TSV manifests.

    Headered manifests should include at least a ``sample`` column and either a
    ``species`` or ``species_fasta`` column. Headerless manifests are interpreted
    as two-column TSVs where the first column is ``sample`` and the second column
    is ``species``. Extra columns are ignored.
    """
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return
        if _looks_like_header(first):
            dict_reader = csv.DictReader(handle, fieldnames=first, delimiter="\t")
            for row in dict_reader:
                if row.get("sample"):
                    yield row
            return
        yield from _iter_headerless_sample_ref_rows([first, *reader], path)


def _iter_headerless_sample_ref_rows(rows: Sequence[Sequence[str]], path: Path) -> Iterable[dict]:
    for line_number, fields in enumerate(rows, start=1):
        if not fields or not fields[0].strip() or fields[0].lstrip().startswith("#"):
            continue
        if len(fields) < 2:
            raise ValueError(
                f"Headerless sample_ref_file {path} line {line_number} must contain "
                "at least sample and species columns"
            )
        yield {
            "sample": fields[0].strip(),
            "species": fields[1].strip(),
        }


def sample_from_row(row: dict, cfg: configparser.ConfigParser) -> Sample:
    sample_name = row["sample"]
    species = row.get("species") or None
    species_fasta = Path(row["species_fasta"]) if row.get("species_fasta") else find_species_fasta(species or sample_name, cfg)
    vcf = Path(row["vcf"]) if row.get("vcf") else find_sample_file(sample_name, cfg, "vcf_dir", "vcf_pattern", "VCF")
    cov = Path(row["cov"]) if row.get("cov") else find_sample_file(sample_name, cfg, "cov_dir", "cov_pattern", "COV")
    return Sample(
        sample_name,
        require_existing_file(species_fasta, "species FASTA"),
        require_gzipped_vcf(vcf),
        require_existing_file(cov, "COV"),
        row.get("species_chrom") or "chrM",
        int(row["rotate_anchor"]) if row.get("rotate_anchor") else None,
        row.get("target_sequence") or None,
        species,
    )


def load_rotate_overrides(cfg: configparser.ConfigParser) -> Dict[str, int]:
    path = cfg.get("paths", "rotate_pos_file", fallback="").strip()
    if not path:
        return {}
    override_path = Path(path)
    if not override_path.exists():
        raise FileNotFoundError(f"rotate_pos_file not found: {override_path}")
    overrides: Dict[str, int] = {}
    with override_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            overrides[row["sample"]] = int(row["rotate_anchor"])
    return overrides


def load_samples(cfg: configparser.ConfigParser, sample_filter: Optional[str]) -> List[Sample]:
    samples: List[Sample] = []
    ref = cfg.get("paths", "sample_ref_file", fallback="").strip()
    if ref:
        ref_path = Path(ref)
        if not ref_path.exists():
            raise FileNotFoundError(
                f"sample_ref_file not found: {ref_path}. Set [paths] sample_ref_file "
                "to an existing TSV or leave it blank and configure [paths] samples."
            )
        for row in iter_sample_ref_rows(ref_path):
            samples.append(sample_from_row(row, cfg))
    else:
        for name in [s.strip() for s in cfg.get("paths", "samples", fallback="").split(",") if s.strip()]:
            sec = f"sample:{name}"
            row = {
                "sample": name,
                "species": cfg.get(sec, "species", fallback=""),
                "species_fasta": cfg.get(sec, "species_fasta", fallback=""),
                "vcf": cfg.get(sec, "vcf", fallback=""),
                "cov": cfg.get(sec, "cov", fallback=""),
                "species_chrom": cfg.get(sec, "species_chrom", fallback="chrM"),
                "rotate_anchor": cfg.get(sec, "rotate_anchor", fallback=""),
                "target_sequence": cfg.get(sec, "target_sequence", fallback=""),
            }
            samples.append(sample_from_row(row, cfg))
    rotate_overrides = load_rotate_overrides(cfg)
    for s in samples:
        if s.name in rotate_overrides:
            s.rotate_anchor = rotate_overrides[s.name]
    if sample_filter:
        samples = [s for s in samples if s.name == sample_filter]
    if not samples:
        raise ValueError("No samples selected; configure sample_ref_file or [paths] samples")
    return samples


def lift_vcf(sample: Sample, pos_map: Dict[int, dict], out_vcf: Path, human_seq: str, target_chrom: str, check_ref: bool, fail_on_mismatch: bool, stats: SampleStats) -> None:
    with gzip.open(sample.vcf, "rt") as inp, out_vcf.open("w") as out:
        for line in inp:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write('##INFO=<ID=SRC_CHROM,Number=1,Type=String,Description="Original species chromosome">\n')
                out.write('##INFO=<ID=SRC_POS,Number=1,Type=Integer,Description="Original species position">\n')
                out.write('##INFO=<ID=SRC_REF,Number=1,Type=String,Description="Original species REF allele">\n')
                out.write('##INFO=<ID=SRC_ALT,Number=.,Type=String,Description="Original species ALT allele(s)">\n')
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            src_chrom, src_pos, _id, ref, alt = parts[:5]
            row = pos_map.get(int(src_pos))
            if not row or row["map_status"] != "mapped":
                stats.vcf_variants_failed_liftover += 1
                continue
            pos = int(row["human_pos_canonical"])
            human_ref = human_seq[pos - 1 : pos - 1 + len(ref)]
            mismatch = check_ref and human_ref.upper() != ref.upper()
            stats.ref_mismatch_count += int(mismatch)
            stats.ref_match_count += int(not mismatch)
            if mismatch and fail_on_mismatch:
                stats.vcf_variants_failed_liftover += 1
                continue
            info = parts[7] if parts[7] not in {"", "."} else ""
            src_info = f"SRC_CHROM={src_chrom};SRC_POS={src_pos};SRC_REF={ref};SRC_ALT={alt}"
            parts[0] = target_chrom
            parts[1] = str(pos)
            parts[7] = src_info if not info else info + ";" + src_info
            out.write("\t".join(parts) + "\n")
            stats.vcf_variants_lifted += 1


def lift_cov(sample: Sample, pos_map: Dict[int, dict], out_cov: Path, target_chrom: str, stats: SampleStats) -> None:
    with sample.cov.open() as inp, out_cov.open("w") as out:
        header = inp.readline().rstrip("\n")
        cols = header.split("\t")
        has_header = any(c.lower() in {"chrom", "pos", "position"} for c in cols)
        if has_header:
            out.write(header + "\tSRC_CHROM\tSRC_POS\n")
        else:
            inp.seek(0)
        for line in inp:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pos_idx = 1 if not has_header else next((i for i, c in enumerate(cols) if c.lower() in {"pos", "position"}), 1)
            chrom_idx = 0
            src_pos = int(parts[pos_idx])
            row = pos_map.get(src_pos)
            if not row or row["map_status"] != "mapped":
                stats.cov_positions_failed_liftover += 1
                continue
            src_chrom = parts[chrom_idx]
            parts[chrom_idx] = target_chrom
            parts[pos_idx] = row["human_pos_canonical"]
            out.write("\t".join(parts + [src_chrom, str(src_pos)]) + "\n")
            stats.cov_positions_lifted += 1


def write_report(stats: SampleStats, path: Path) -> None:
    fields = [f for f in stats.__dataclass_fields__ if f != "notes"]
    with path.open("w") as out:
        for f in fields:
            out.write(f"{f}\t{getattr(stats, f)}\n")
        for note in stats.notes:
            out.write(f"note\t{note}\n")


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sample")
    args = ap.parse_args(argv)

    cfg = configparser.ConfigParser()
    cfg.read(args.config)
    outdir = Path(cfg.get("paths", "output_dir"))
    dirs = {d: outdir / d for d in OUTDIRS}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    human_raw = read_fasta(Path(cfg["paths"]["human_fasta"]), cfg.get("fasta", "human_target_sequence", fallback="").strip() or None)
    human_name = cfg.get("fasta", "human_sequence_name", fallback="human_chrM")
    human_prepared = dirs["prepared_fastas"] / "human_chrM.prepared.fa"
    write_fasta(human_prepared, human_name, human_raw.seq)

    human_len = cfg.getint("coordinates", "human_len")
    restore_offset = cfg.getint("coordinates", "human_restore_offset", fallback=0)
    target_chrom = cfg.get("coordinates", "target_chrom", fallback="chrM")
    check_ref = cfg.getboolean("liftover", "check_ref_against_human_fasta", fallback=True)
    fail_ref = cfg.getboolean("liftover", "fail_on_ref_mismatch", fallback=False)
    human_anchor_cfg = cfg.get("coordinates", "human_rotate_anchor", fallback="").strip()

    summaries: List[SampleStats] = []
    for sample in load_samples(cfg, args.sample):
        species_raw = read_fasta(sample.species_fasta, sample.target_sequence)
        species_name = cfg.get("fasta", "species_sequence_name_template", fallback="{sample}_chrM").format(sample=sample.name)
        species_prepared = dirs["prepared_fastas"] / f"{sample.name}.prepared.fa"
        write_fasta(species_prepared, species_name, species_raw.seq)

        species_anchor = sample.rotate_anchor or infer_anchor(species_raw.seq, human_raw.seq)
        human_anchor = int(human_anchor_cfg) if human_anchor_cfg else infer_anchor(human_raw.seq, species_raw.seq)
        species_rotated = dirs["rotated_fastas"] / f"{sample.name}.rotated.fa"
        human_rotated = dirs["rotated_fastas"] / f"{sample.name}.human.rotated.fa"
        write_fasta(species_rotated, species_name, rotate_sequence(species_raw.seq, species_anchor))
        write_fasta(human_rotated, human_name, rotate_sequence(human_raw.seq, human_anchor))

        aln = dirs["alignments"] / f"{sample.name}.species_to_human.aligned.fa"
        run_alignment(species_rotated, human_rotated, aln, cfg)
        map_path = dirs["maps"] / f"{sample.name}.coordinate_map.tsv"
        pos_map, stats = build_map(sample, aln, map_path, species_anchor, human_anchor, human_len, restore_offset)
        stats.notes.append(f"species_rotate_anchor={species_anchor}")
        stats.notes.append(f"human_rotate_anchor={human_anchor}")

        lift_vcf(sample, pos_map, dirs["vcf_lifted_raw"] / f"{sample.name}.lifted.raw.vcf", human_raw.seq, target_chrom, check_ref, fail_ref, stats)
        lift_cov(sample, pos_map, dirs["cov_lifted"] / f"{sample.name}.lifted.cov", target_chrom, stats)
        write_report(stats, dirs["reports"] / f"{sample.name}.coordinate_liftover_qc.tsv")
        summaries.append(stats)

    summary_path = dirs["reports"] / "all_samples.coordinate_liftover_summary.tsv"
    fields = [f for f in SampleStats.__dataclass_fields__ if f != "notes"]
    with summary_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow({f: getattr(s, f) for f in fields})
    return 0


if __name__ == "__main__":
    sys.exit(main())
