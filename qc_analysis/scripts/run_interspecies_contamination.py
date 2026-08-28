#!/usr/bin/env python3
"""Report cohort-level cross-species contamination from lifted mtDNA VCFs."""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

FIELDS = ("sample species interspecies_status classification reason recipient_species_n "
          "n_lowA n_lowA_after_species_background best_source_species best_source_sample "
          "overlap_count overlap_fraction matched_low_vaf_median vaf_coherence "
          "source_species_count source_sample_count").split()


def resolve(value: object) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def number(value: str) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def read_metadata(path: Path, sample_col: str, species_col: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not rows.fieldnames or sample_col not in rows.fieldnames or species_col not in rows.fieldnames:
            raise ValueError(f"metadata must contain {sample_col!r} and {species_col!r}: {path}")
        result = {}
        for row in rows:
            sample, species = row[sample_col].strip(), row[species_col].strip()
            if sample and species:
                if sample in result and result[sample] != species:
                    raise ValueError(f"conflicting species for sample {sample}")
                result[sample] = species
        return result


def discover(directory: Path, pattern: str) -> dict[str, Path]:
    if "{sample}" not in pattern:
        raise ValueError("input_vcf_pattern must contain {sample}")
    prefix, suffix = pattern.split("{sample}", 1)
    found = {}
    for ending in dict.fromkeys((suffix, suffix[:-3] if suffix.endswith(".gz") else suffix + ".gz")):
        for path in directory.glob(f"{prefix}*{ending}"):
            sample = path.name[len(prefix):len(path.name)-len(ending)]
            if sample and path.is_file():
                if sample in found:
                    raise ValueError(f"ambiguous lifted VCF for {sample}: {found[sample]}, {path}")
                found[sample] = path
    return found


def alleles(path: Path, dp_min: float) -> list[tuple[tuple[str, int, str, str], float]]:
    opener = gzip.open if path.suffix == ".gz" else open
    result = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] != "PASS" or "," in f[4]:
                continue
            ref, alt = f[3].upper(), f[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            dp = number(fmt.get("DP", "")) or number(info.get("DP", ""))
            af = number(fmt.get("AF", "").split(",")[0])
            ad = [number(x) for x in fmt.get("AD", "").split(",")]
            if af is None and len(ad) == 2 and None not in ad and sum(ad) > 0:
                af = ad[1] / sum(ad)
            try:
                key = (f[0], int(f[1]), ref, alt)
            except ValueError:
                continue
            if dp is not None and dp >= dp_min and af is not None:
                result.append((key, af))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    sec = read_simple_yaml(args.config).get("interspecies_contamination") or {}
    if sec.get("enabled", True) is False:
        print("[interspecies_contamination] disabled; skipping.")
        return 0
    paths, settings = sec.get("paths", {}) or {}, sec.get("settings", {}) or {}
    vcf_dir = resolve(paths.get("input_vcf_dir", "results/qc/coordinate_liftover/vcf_lifted_raw"))
    metadata_path = resolve(paths.get("sample_ref_file", "config/sample_ref_file.tsv"))
    output_dir = resolve(paths.get("output_dir", "results/qc/interspecies_contamination"))
    output = output_dir / "reports/interspecies_contamination_report.tsv"
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"output exists (use --overwrite): {output}")
    vcfs = discover(vcf_dir, str(settings.get("input_vcf_pattern", "{sample}.lifted.raw.vcf")))
    if not vcfs:
        raise ValueError(f"no post-liftover VCFs found in {vcf_dir}")
    metadata = read_metadata(metadata_path, str(paths.get("metadata_sample_column", "sample")),
                             str(paths.get("metadata_species_column", "species")))
    missing = sorted(set(vcfs) - set(metadata))
    if missing:
        raise ValueError("lifted VCF samples missing species metadata: " + ", ".join(missing))
    dp_min = float(settings.get("dp_min", 100)); low_min = float(settings.get("low_vaf_min", .01))
    low_max = float(settings.get("low_vaf_max", .20)); high_min = float(settings.get("high_vaf_min", .99))
    min_overlap = int(settings.get("min_overlap", 3)); min_fraction = float(settings.get("min_overlap_fraction", .5))
    tolerance = float(settings.get("vaf_coherence_tolerance", .03)); min_coherence = float(settings.get("min_vaf_coherence", .7))
    calls = {sample: alleles(path, dp_min) for sample, path in vcfs.items()}
    species_samples = defaultdict(set)
    high_index = defaultdict(list)  # allele inverted index: never all-vs-all intersections
    for sample, rows in calls.items():
        species_samples[metadata[sample]].add(sample)
        for allele, af in rows:
            if af >= high_min:
                high_index[allele].append(sample)
    report = []
    for recipient in sorted(vcfs):
        species = metadata[recipient]
        raw_low = [(key, af) for key, af in calls[recipient] if low_min <= af <= low_max]
        retained = [(key, af) for key, af in raw_low if not any(
            other != recipient and metadata[other] == species for other in high_index.get(key, ()))]
        by_species, by_sample = defaultdict(dict), defaultdict(dict)
        for key, af in retained:
            for source in high_index.get(key, ()):
                if source != recipient and metadata[source] != species:
                    by_species[metadata[source]][key] = af
                    by_sample[source][key] = af
        ranked_species = sorted(by_species, key=lambda x: (-len(by_species[x]), x))
        best_species = ranked_species[0] if ranked_species else ""
        overlap = len(by_species.get(best_species, {})); denominator = len(retained)
        fraction = overlap / denominator if denominator else 0.0
        eligible_samples = [s for s in by_sample if metadata[s] == best_species]
        ranked_samples = sorted(eligible_samples, key=lambda x: (-len(by_sample[x]), x))
        best_sample = ranked_samples[0] if ranked_samples else ""
        values = list(by_species.get(best_species, {}).values())
        median = statistics.median(values) if values else None
        coherence = sum(abs(v - median) <= tolerance for v in values) / len(values) if values else 0.0
        tied_species = len(ranked_species) > 1 and len(by_species[ranked_species[0]]) == len(by_species[ranked_species[1]])
        strong = overlap >= min_overlap and fraction >= min_fraction
        if not retained:
            status, classification, reason = "PASS", "NO_INFORMATIVE_LOW_VAF", "no low-VAF alleles remain after recipient-species background removal"
        elif not strong:
            status, classification, reason = "PASS", "NO_CROSS_SPECIES_SIGNAL", "cross-species overlap is below configured thresholds"
        elif len(species_samples[species]) == 1:
            status, classification, reason = "WARN", "SINGLETON_RECIPIENT_SPECIES", "recipient-species homoplasmic background cannot be established"
        elif tied_species:
            status, classification, reason = "WARN", "AMBIGUOUS_SOURCE_SPECIES", "multiple source species have equal best overlap"
        elif coherence < min_coherence:
            status, classification, reason = "WARN", "VAF_INCOHERENT", "matched low-VAF alleles do not meet coherence threshold"
        else:
            status, classification, reason = "FAIL", "INTERSPECIES_CONTAMINATION", "coherent low-VAF alleles match a different-species homoplasmic source"
        report.append(dict(sample=recipient, species=species, interspecies_status=status,
                           classification=classification, reason=reason,
                           recipient_species_n=len(species_samples[species]), n_lowA=len(raw_low),
                           n_lowA_after_species_background=len(retained), best_source_species=best_species,
                           best_source_sample=best_sample, overlap_count=overlap,
                           overlap_fraction=f"{fraction:.6f}", matched_low_vaf_median="NA" if median is None else f"{median:.6f}",
                           vaf_coherence=f"{coherence:.6f}", source_species_count=len(by_species),
                           source_sample_count=len(by_sample)))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(report)
    temporary.replace(output)
    print(f"Wrote {output} ({len(report)} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
