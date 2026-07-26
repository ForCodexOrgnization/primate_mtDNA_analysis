"""Reusable tRNAscan-SE parsing, execution, and position-index utilities.

The ``.ss`` sequence is always interpreted in mature tRNA (5' to 3')
orientation.  Genomic bases are independently read from the input FASTA.
"""
from __future__ import annotations

import csv
import gzip
import shlex
import subprocess
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

INDEX_FORMAT_VERSION = "2"
INDEX_COLUMNS = ["index_format_version", "base_orientation", "pair_type_orientation",
                 "coordinate_space", "reference_key", "chrom", "pos", "trna_id",
                 "trna_begin", "trna_end", "strand", "aa", "anticodon", "score",
                 "local_pos", "base_genomic", "base_rna", "struct_char", "struct_class",
                 "struct_element", "paired_local_pos", "paired_genomic_pos",
                 "paired_base_genomic", "paired_base_rna", "pair_bases_rna", "pair_type",
                 "pair_state", "base", "paired_base"]


@dataclass
class TRNARecord:
    chrom: str
    number: str
    begin: int
    end: int
    aa: str
    anticodon: str
    score: float | str = ""
    sequence: str = ""
    structure: str = ""
    pairs: dict[int, int] = field(default_factory=dict)

    @property
    def strand(self) -> str:
        return "+" if self.begin <= self.end else "-"

    @property
    def trna_id(self) -> str:
        return f"{self.chrom}.trna{self.number}"

    def genomic_pos(self, local_pos: int) -> int:
        return self.begin + local_pos - 1 if self.strand == "+" else self.begin - local_pos + 1


def _open(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def parse_trnascan_out(path: str | Path) -> list[TRNARecord]:
    """Parse standard tRNAscan-SE tabular ``-o`` output (v1/v2 headings)."""
    records = []
    with _open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("Sequence") or set(line.strip()) <= {"-", " ", "\t"}:
                continue
            cols = line.split()
            if len(cols) < 9:
                continue
            try:
                begin, end = int(cols[2]), int(cols[3])
                score = float(cols[8])
            except ValueError:
                continue
            records.append(TRNARecord(cols[0], cols[1], begin, end, cols[4], cols[5], score))
    return records


def parse_trnascan_ss(path: str | Path) -> list[dict]:
    """Parse tRNAscan-SE ``-f`` records, accepting wrapped Seq/Str values."""
    result, current = [], None
    with _open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if not line.startswith(("Seq:", "Str:", "Type:")) and re.search(r"\.trna?\d+\s+\([^)]*\)", line, re.I):
                if current: result.append(current)
                current={"id":line.split()[0],"sequence":"","structure":""}
            elif line.startswith("Seq:"):
                value = line[4:].strip()
                # A record header contains coordinates; sequence lines contain bases.
                if "(" in value and ")" in value:
                    if current:
                        result.append(current)
                    current = {"id": value.split()[0], "sequence": "", "structure": ""}
                elif current:
                    current["sequence"] += "".join(value.split()).upper()
            elif line.startswith("Str:") and current:
                current["structure"] += "".join(line[4:].split())
            elif line.startswith("Type:") and current:
                bits = line.replace(":", " ").split()
                current["aa"] = bits[1] if len(bits) > 1 else ""
                current["anticodon"] = bits[bits.index("Anticodon") + 1] if "Anticodon" in bits else ""
    if current:
        result.append(current)
    if not result:
        raise ValueError(f"Missing or malformed tRNAscan secondary-structure records: {path}")
    for record in result:
        if not record["sequence"] or not record["structure"]:
            raise ValueError(f"Incomplete tRNAscan secondary-structure record {record['id']} in {path}")
        if len(record["sequence"]) != len(record["structure"]):
            raise ValueError(f"Sequence/structure length mismatch for {record['id']} in {path}")
    return result


def build_pairs(structure: str) -> dict[int, int]:
    """Return 1-based local pairing positions for bracket and angle notation."""
    pairs, stacks = {}, {">": [], "(": [], "[": [], "{": []}
    closing = {"<": ">", ")": "(", "]": "[", "}": "{"}
    for pos, char in enumerate(structure, 1):
        if char in stacks:
            stacks[char].append(pos)
        elif char in closing:
            stack = stacks[closing[char]]
            if not stack:
                raise ValueError(f"Unbalanced structure at position {pos}")
            mate = stack.pop()
            pairs[pos], pairs[mate] = mate, pos
    if any(stacks.values()):
        raise ValueError("Unbalanced tRNA secondary structure")
    return pairs


def infer_structural_elements(structure: str, pairs: dict[int, int] | None = None) -> dict[int, str]:
    """Assign conventional tRNA elements by local position and stem topology."""
    pairs = pairs or build_pairs(structure)
    n = len(structure)
    out = {}
    # Canonical mature-tRNA ranges; proportional fallback keeps truncated calls useful.
    ranges = [(1, 7, "acceptor_stem"), (8, 9, "d_arm_connector"), (10, 25, "d_arm"),
              (26, 27, "anticodon_connector"), (28, 42, "anticodon_arm"),
              (43, 48, "variable_loop"), (49, 65, "t_arm"), (66, 72, "acceptor_stem"),
              (73, 10**9, "discriminator_or_CCA")]
    for pos in range(1, n + 1):
        out[pos] = next(name for lo, hi, name in ranges if lo <= pos <= hi)
    return out


def merge_trnascan_records(out_records: list[TRNARecord], ss_records: list[dict]) -> list[TRNARecord]:
    """Merge output and structure records without relying on record ordering alone."""
    unused = list(ss_records)
    for record in out_records:
        aliases = {record.trna_id, f"{record.chrom}.trna{record.number}", f"{record.chrom}.tRNA{record.number}"}
        match = next((x for x in unused if x["id"] in aliases), None)
        if match is None and len(unused) == len(out_records):
            match = unused[0]
        if match is None:
            raise ValueError(f"No .ss record found for {record.trna_id}")
        unused.remove(match)
        record.sequence, record.structure = match["sequence"], match["structure"]
        record.pairs = build_pairs(record.structure)
    if not out_records:
        raise ValueError("No tRNAs detected in tRNAscan output")
    return out_records


def trnascan_mode_args(mode: str) -> list[str]:
    modes = {"euk": ["-E"], "bact": ["-B"], "arch": ["-A"], "general": ["-G"],
             "mito_mammal": ["-M", "mammal"], "mito_vert": ["-M", "vert"],
             "organellar": ["-O"]}
    if mode not in modes:
        raise ValueError(f"Unsupported tRNAscan mode {mode!r}; choose {', '.join(modes)}")
    return modes[mode]


def run_trnascan(fasta: str | Path, prefix: str | Path, trnascan_bin="tRNAscan-SE",
                 trnascan_mode="mito_mammal", trnascan_threads=1,
                 trnascan_extra_args: str | Iterable[str] = "") -> dict[str, Path]:
    prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {k: Path(str(prefix) + suffix) for k, suffix in
               {"out": ".trnascan.out", "ss": ".trnascan.ss", "stats": ".trnascan.stats",
                "bed": ".trnascan.bed", "fasta": ".trnascan.fa"}.items()}
    extra = shlex.split(trnascan_extra_args) if isinstance(trnascan_extra_args, str) else list(trnascan_extra_args)
    command = [str(trnascan_bin), *trnascan_mode_args(trnascan_mode), "--thread", str(trnascan_threads),
               "-o", str(outputs["out"]), "-f", str(outputs["ss"]), "-m", str(outputs["stats"]),
               "-b", str(outputs["bed"]), "-a", str(outputs["fasta"]), *extra, str(fasta)]
    subprocess.run(command, check=True)
    for key in ("out", "ss"):
        if not outputs[key].is_file() or outputs[key].stat().st_size == 0:
            raise RuntimeError(f"tRNAscan-SE did not create a nonempty {outputs[key]}")
    return outputs


def read_fasta(path: str | Path) -> dict[str, str]:
    records, name, chunks = {}, None, []
    with _open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None: records[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else: chunks.append(line.strip())
    if name is not None: records[name] = "".join(chunks).upper()
    if not records: raise ValueError(f"No FASTA records in {path}")
    return records


def _rna(genomic: str, strand: str) -> str:
    base = genomic.upper().replace("T", "U")
    return base if strand == "+" else {"A": "U", "U": "A", "C": "G", "G": "C"}.get(base, base)


def _pair_type(a: str, b: str) -> str:
    pair = (a, b)
    if pair in {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}: return "WC"
    if pair in {("G", "U"), ("U", "G")}: return "GU_wobble"
    return "non_WC"


def build_trna_position_index(reference_key: str, fasta: str | Path, trnascan_out: str | Path,
                              trnascan_ss: str | Path, output: str | Path,
                              chrom_normalization="none", mismatch_rate_threshold=0.0) -> dict:
    """Build the v2 reference-coordinate position index and return build metrics."""
    seqs = read_fasta(fasta)
    records = merge_trnascan_records(parse_trnascan_out(trnascan_out), parse_trnascan_ss(trnascan_ss))
    rows, mismatches, compared = [], 0, 0
    for rec in records:
        chrom = rec.chrom
        if chrom not in seqs and len(seqs) == 1: chrom = next(iter(seqs))
        if chrom not in seqs: raise ValueError(f"tRNAscan sequence {rec.chrom!r} not found in FASTA")
        out_chrom = chrom
        if chrom_normalization == "strip_chr" and out_chrom.lower().startswith("chr"): out_chrom = out_chrom[3:]
        elif chrom_normalization == "add_chr" and not out_chrom.lower().startswith("chr"): out_chrom = "chr" + out_chrom
        elif chrom_normalization == "mitochondrial_alias": out_chrom = "MT" if out_chrom.lower().removeprefix("chr") in {"m", "mt"} else out_chrom
        elif chrom_normalization != "none": raise ValueError(f"Unsupported chromosome normalization: {chrom_normalization}")
        elements = infer_structural_elements(rec.structure, rec.pairs)
        for local, ss_base in enumerate(rec.sequence, 1):
            pos = rec.genomic_pos(local)
            if pos < 1 or pos > len(seqs[chrom]): raise ValueError(f"{rec.trna_id} position {pos} outside FASTA")
            genomic = seqs[chrom][pos - 1]; rna = _rna(genomic, rec.strand)
            if ss_base.upper().replace("T", "U") not in {rna, "N"}: mismatches += 1
            compared += 1
            mate = rec.pairs.get(local); mate_pos = rec.genomic_pos(mate) if mate else None
            mate_genomic = seqs[chrom][mate_pos - 1] if mate_pos else ""
            mate_rna = _rna(mate_genomic, rec.strand) if mate else ""
            ptype = _pair_type(rna, mate_rna) if mate else "."
            struct_class = "stem" if mate else "loop"
            rows.append({"index_format_version": INDEX_FORMAT_VERSION, "base_orientation": "genomic_and_rna",
                         "pair_type_orientation": "transcript_rna", "coordinate_space": "original_reference",
                         "reference_key": reference_key, "chrom": out_chrom, "pos": pos, "trna_id": rec.trna_id,
                         "trna_begin": rec.begin, "trna_end": rec.end, "strand": rec.strand, "aa": rec.aa,
                         "anticodon": rec.anticodon, "score": rec.score, "local_pos": local,
                         "base_genomic": genomic, "base_rna": rna, "struct_char": rec.structure[local-1],
                         "struct_class": struct_class, "struct_element": elements[local],
                         "paired_local_pos": mate or ".", "paired_genomic_pos": mate_pos or ".",
                         "paired_base_genomic": mate_genomic or ".", "paired_base_rna": mate_rna or ".",
                         "pair_bases_rna": f"{rna}-{mate_rna}" if mate else ".", "pair_type": ptype,
                         "pair_state": "paired" if mate else "unpaired", "base": genomic,
                         "paired_base": mate_genomic or "."})
    rate = mismatches / compared if compared else 0
    if rate > mismatch_rate_threshold:
        raise ValueError(f"FASTA/tRNAscan sequence mismatch rate {rate:.6g} exceeds threshold {mismatch_rate_threshold:.6g} ({mismatches}/{compared})")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with _open(output, "wt") as handle:
        writer = csv.DictWriter(handle, INDEX_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return {"records": records, "rows": rows, "fasta_length": sum(map(len, seqs.values())),
            "n_fasta_sequence_mismatch": mismatches}
