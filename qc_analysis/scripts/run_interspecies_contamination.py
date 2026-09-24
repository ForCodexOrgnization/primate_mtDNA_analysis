#!/usr/bin/env python3
"""Report cohort-level cross-species contamination from current lifted mtDNA VCFs.

Production PASS/WARN/FAIL classification intentionally preserves the historical
hard-threshold logic. In parallel, this module reports a 0-1 evidence score that
mirrors the intraspecies framework while accounting for cross-species allele
specificity, genome-wide dispersion, AF coherence, source-sample concentration,
source-species separation, source-genus background correction, and project/cohort provenance.

Project/cohort provenance modifies only the source-matching component. It never
creates contamination evidence by itself.
"""
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

FIELDS = """sample species interspecies_status classification reason recipient_species_n
n_lowA n_lowA_after_species_background best_source_species best_source_sample
overlap_count overlap_fraction best_source_sample_overlap best_source_sample_fraction
best_source_species_overlap best_source_species_fraction matched_low_vaf_median
matched_low_vaf_mad vaf_coherence source_species_count source_sample_count
target_project target_cohort best_source_project best_source_cohort
target_source_same_project target_source_same_cohort provenance_relationship_basis
provenance_source_factor source_specificity_assessable
target_genus best_source_genus source_genus_specificity_assessable best_source_genus_species_n
best_overlap_background_adjusted best_fraction_background_adjusted
best_overlap_mean_cross_species_frequency best_overlap_mean_source_specificity
best_overlap_mean_source_genus_frequency best_overlap_mean_source_genus_specificity
best_overlap_mean_effective_source_specificity
best_overlap_effective_specific_overlap best_source_sample_effective_overlap
best_source_sample_concentration best_source_species_adjusted_overlap
second_source_species second_source_species_adjusted_overlap best_source_species_separation
overlap_positions overlap_occupied_bins overlap_bin_entropy_normalized
overlap_circular_span_bp overlap_max_local_fraction
contamination_score_gate_pass contamination_score_version
contamination_score_source_basis contamination_score_source_composite_index
contamination_score_source_total_pre_provenance contamination_score_source_total
contamination_score_dispersion contamination_score_af_coherence
contamination_score_source_concentration contamination_score_species_separation
contamination_score_raw_10 contamination_score contamination_score_interpretation""".split()


def resolve(value: object) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def number(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def clean_meta(value) -> str:
    value = str(value or "").strip()
    return "" if not value or value.upper() == "NA" else value


def read_metadata(path: Path, sample_col: str, species_col: str) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Read sample/species plus optional project/cohort from the same metadata table."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}, {}

    first = rows[0]
    headered = sample_col in first and species_col in first
    species_by_sample: dict[str, str] = {}
    provenance: dict[str, dict[str, str]] = {}

    if headered:
        sample_index = first.index(sample_col)
        species_index = first.index(species_col)
        project_aliases = ("project", "bioproject", "bioproject_accession", "study_accession",
                           "project_accession", "study", "ena_study", "sra_study")
        cohort_aliases = ("cohort", "cohort_id", "cohort_name", "study_cohort",
                          "dataset", "dataset_id", "source_cohort")
        lower_to_name = {x.strip().lower(): x for x in first}
        project_name = next((lower_to_name[x] for x in project_aliases if x in lower_to_name), None)
        cohort_name = next((lower_to_name[x] for x in cohort_aliases if x in lower_to_name), None)
        for values in rows[1:]:
            if len(values) <= max(sample_index, species_index):
                raise ValueError(f"metadata row must contain sample and species columns: {path}")
            row = dict(zip(first, values))
            sample = values[sample_index].strip()
            species = values[species_index].strip()
            if not sample or not species:
                continue
            if sample in species_by_sample and species_by_sample[sample] != species:
                raise ValueError(f"conflicting species for sample {sample}")
            species_by_sample[sample] = species
            provenance[sample] = {
                "project": clean_meta(row.get(project_name)) if project_name else "",
                "cohort": clean_meta(row.get(cohort_name)) if cohort_name else "",
            }
    else:
        for values in rows:
            if len(values) < 2:
                continue
            sample, species = values[0].strip(), values[1].strip()
            if sample and species:
                species_by_sample[sample] = species
                provenance[sample] = {"project": "", "cohort": ""}
    return species_by_sample, provenance


def read_optional_provenance(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and "sample" not in rows[0]:
        return {}
    return {
        str(row.get("sample", "")).strip(): {
            "project": clean_meta(row.get("project")),
            "cohort": clean_meta(row.get("cohort")),
        }
        for row in rows if str(row.get("sample", "")).strip()
    }


def merge_provenance(primary, fallback):
    out = dict(fallback)
    for sample, values in primary.items():
        merged = dict(out.get(sample, {}))
        for key in ("project", "cohort"):
            if values.get(key):
                merged[key] = values[key]
        out[sample] = merged
    return out


def provenance_relationship(target_sample: str, source_sample: str, provenance: dict, settings: dict) -> dict:
    target = provenance.get(target_sample, {})
    source = provenance.get(source_sample, {}) if source_sample else {}
    tp, tc = target.get("project", ""), target.get("cohort", "")
    sp, sc = source.get("project", ""), source.get("cohort", "")
    same_project = None if not tp or not sp else tp == sp
    same_cohort = None if not tc or not sc else tc == sc

    if same_project is False:
        basis = "different_project"
        factor = float(settings.get("provenance_different_project_factor", .60))
    elif same_cohort is False:
        basis = "same_project_different_cohort" if same_project is True else "different_cohort_project_unknown"
        factor = float(settings.get("provenance_different_cohort_factor", .75))
    elif same_cohort is True:
        basis = "same_cohort"
        factor = float(settings.get("provenance_same_cohort_factor", 1.00))
    elif same_project is True:
        basis = "same_project_cohort_unknown"
        factor = float(settings.get("provenance_same_project_factor", .90))
    else:
        basis = "provenance_unknown_neutral"
        factor = float(settings.get("provenance_unknown_factor", 1.00))

    return {
        "target_project": tp or "NA",
        "target_cohort": tc or "NA",
        "best_source_project": sp or "NA",
        "best_source_cohort": sc or "NA",
        "target_source_same_project": "YES" if same_project is True else "NO" if same_project is False else "NA",
        "target_source_same_cohort": "YES" if same_cohort is True else "NO" if same_cohort is False else "NA",
        "provenance_relationship_basis": basis,
        "provenance_source_factor": factor,
    }


def collection_ok_samples(cfg: dict) -> set[str]:
    collect = cfg.get("collect_variant_calling") or {}
    report = resolve(collect.get("outdir", "results/qc/collected_variant_calling_results")) / "reports/variant_calling_collection_summary.tsv"
    if not report.is_file():
        raise ValueError(f"current collection summary is required for interspecies QC: {report}")
    result = set()
    with report.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("sample") or row.get("Sample") or "").strip()
            status = (row.get("status") or row.get("collection_status") or "").strip()
            if sample and status == "OK":
                result.add(sample)
    return result


def liftover_ok_samples(cfg: dict) -> set[str]:
    section = cfg.get("coordinate_liftover") or {}
    paths = section.get("paths") or {}
    reports = resolve(paths.get("output_dir", "results/qc/coordinate_liftover")) / "reports"
    result = set()
    if not reports.is_dir():
        return result
    for report in reports.glob("*.coordinate_liftover_qc.tsv"):
        status = ""
        with report.open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) >= 2 and row[0].strip() == "status":
                    status = row[1].strip()
                    break
        if status == "completed":
            result.add(report.name.removesuffix(".coordinate_liftover_qc.tsv"))
    return result


def discover(directory: Path, pattern: str) -> dict[str, Path]:
    if "{sample}" not in pattern:
        raise ValueError("input_vcf_pattern must contain {sample}")
    prefix, suffix = pattern.split("{sample}", 1)
    found = {}
    endings = dict.fromkeys((suffix, suffix[:-3] if suffix.endswith(".gz") else suffix + ".gz"))
    for ending in endings:
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
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[6] != "PASS" or "," in fields[4]:
                continue
            ref, alt = fields[3].upper(), fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            info = dict(x.split("=", 1) for x in fields[7].split(";") if "=" in x)
            dp = number(fmt.get("DP", ""))
            if dp is None:
                dp = number(info.get("DP", ""))
            af = number(fmt.get("AF", "").split(",")[0])
            ad = [number(x) for x in fmt.get("AD", "").split(",")]
            if af is None and len(ad) == 2 and None not in ad and sum(ad) > 0:
                af = ad[1] / sum(ad)
            try:
                key = (fields[0], int(fields[1]), ref, alt)
            except ValueError:
                continue
            if dp is not None and dp >= dp_min and af is not None:
                result.append((key, af))
    return result


def specificity_weight(freq: float, settings: dict) -> float:
    full = float(settings.get("source_bg_freq_full_weight_max", .05))
    three_quarter = float(settings.get("source_bg_freq_three_quarter_max", .10))
    half = float(settings.get("source_bg_freq_half_weight_max", .25))
    zero = float(settings.get("source_bg_freq_zero_weight_min", .50))
    if freq <= full:
        return 1.0
    if freq <= three_quarter:
        return .75
    if freq <= half:
        return .50
    if freq < zero:
        return .25
    return 0.0


def genus_name(species: str) -> str:
    """Return the genus token from normalized or whitespace-delimited species labels."""
    text = str(species or "").strip().replace(" ", "_")
    return text.split("_", 1)[0] if text else ""


def genus_specificity(carrier_n: int, eligible_n: int) -> float:
    """Return 1 for source-species-unique alleles and 0 for genus-wide alleles.

    The normalization makes an allele present in exactly one of the eligible
    source-genus species retain full specificity even when the genus contains
    only a few sampled species:
        (eligible_n - carrier_n) / (eligible_n - 1)
    """
    if eligible_n < 2:
        return 1.0
    carrier_n = max(1, min(int(carrier_n), int(eligible_n)))
    return max(0.0, min(1.0, (eligible_n - carrier_n) / (eligible_n - 1)))


def median_abs_deviation(values: list[float]) -> float | None:
    if not values:
        return None
    med = statistics.median(values)
    return statistics.median(abs(x - med) for x in values)


def overlap_dispersion(keys, mt_length: int, bin_bp: int, window_bp: int) -> dict:
    positions = sorted({int(key[1]) for key in keys})
    if not positions:
        return {
            "overlap_positions": "",
            "overlap_occupied_bins": 0,
            "overlap_bin_entropy_normalized": 0.0,
            "overlap_circular_span_bp": 0,
            "overlap_max_local_fraction": 0.0,
        }
    total_bins = max(1, math.ceil(mt_length / bin_bp))
    counts = defaultdict(int)
    for pos in positions:
        counts[min((pos - 1) // bin_bp, total_bins - 1)] += 1
    probs = [count / len(positions) for count in counts.values()]
    entropy = -sum(p * math.log(p) for p in probs if p > 0)
    entropy_norm = entropy / math.log(total_bins) if total_bins > 1 else 0.0

    if len(positions) == 1:
        circular_span = 0
    else:
        extended = positions + [positions[0] + mt_length]
        gaps = [extended[i + 1] - extended[i] for i in range(len(positions))]
        circular_span = mt_length - max(gaps)

    max_in_window = 0
    for start in positions:
        n = sum(min((pos - start) % mt_length, (start - pos) % mt_length) <= window_bp / 2 for pos in positions)
        max_in_window = max(max_in_window, n)

    return {
        "overlap_positions": ",".join(str(x) for x in positions),
        "overlap_occupied_bins": len(counts),
        "overlap_bin_entropy_normalized": entropy_norm,
        "overlap_circular_span_bp": circular_span,
        "overlap_max_local_fraction": max_in_window / len(positions),
    }


def dispersion_points(entropy: float) -> float:
    if entropy >= .90:
        return 2.0
    if entropy >= .85:
        return 1.5
    if entropy >= .75:
        return 1.0
    if entropy >= .65:
        return .5
    return 0.0


def af_coherence_points(mad: float | None, overlap: int) -> float:
    if mad is None:
        return 0.0
    if mad <= .005:
        score = 1.5
    elif mad <= .01:
        score = 1.0
    elif mad <= .02:
        score = .5
    else:
        score = 0.0
    return min(score, 1.0) if overlap < 5 else score


def concentration_points(value: float) -> float:
    if value >= .80:
        return 1.5
    if value >= .60:
        return 1.0
    if value >= .40:
        return .5
    return 0.0


def species_separation_points(value: float) -> float:
    """Score how clearly the best source species exceeds the runner-up."""
    if value >= .60:
        return 1.0
    if value >= .40:
        return .67
    if value >= .20:
        return .33
    return 0.0


def score_interpretation(score: float | None) -> str:
    if score is None:
        return "not_scored_insufficient_evidence"
    if score >= .70:
        return "strong_evidence"
    if score >= .50:
        return "candidate_evidence"
    if score >= .30:
        return "weak_ambiguous_evidence"
    return "little_evidence"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("interspecies_contamination") or {}
    if sec.get("enabled", True) is False:
        print("[interspecies_contamination] disabled; skipping.")
        return 0
    paths, settings = sec.get("paths", {}) or {}, sec.get("settings", {}) or {}

    vcf_dir = resolve(paths.get("input_vcf_dir", "results/qc/coordinate_liftover/vcf_lifted_raw"))
    metadata_path = resolve(paths.get("sample_ref_file", "config/sample_ref_file.tsv"))
    output_dir = resolve(paths.get("output_dir", "results/qc/interspecies_contamination"))
    output = output_dir / "reports/interspecies_contamination_report.tsv"
    if output.exists():
        print(f"[interspecies_contamination] replacing existing report: {output}", file=sys.stderr)

    metadata, metadata_provenance = read_metadata(
        metadata_path,
        str(paths.get("metadata_sample_column", "sample")),
        str(paths.get("metadata_species_column", "species")),
    )
    provenance_path = resolve(paths.get(
        "sample_provenance_file",
        "results/qc/sample_deduplication/reports/deduplicated_sample_ref_file.tsv",
    ))
    provenance = merge_provenance(metadata_provenance, read_optional_provenance(provenance_path))

    current_ok = collection_ok_samples(cfg)
    liftover_ok = liftover_ok_samples(cfg)
    current_samples = set(metadata) & current_ok & liftover_ok
    discovered = discover(vcf_dir, str(paths.get("input_vcf_pattern", "{sample}.lifted.raw.vcf")))
    ignored = sorted(set(discovered) - current_samples)
    if ignored:
        print("[interspecies_contamination] ignoring stale/non-current lifted VCFs: " + ", ".join(ignored[:20]), file=sys.stderr)
    vcfs = {sample: path for sample, path in discovered.items() if sample in current_samples}
    if not vcfs:
        raise ValueError(f"no current collection-OK, liftover-completed VCFs found in {vcf_dir}")

    dp_min = float(settings.get("dp_min", 100))
    low_min = float(settings.get("low_vaf_min", .01))
    low_max = float(settings.get("low_vaf_max", .20))
    high_min = float(settings.get("high_vaf_min", .99))
    min_overlap = int(settings.get("min_overlap", 3))
    min_fraction = float(settings.get("min_overlap_fraction", .5))
    min_informative = int(settings.get("min_informative_low_variants", 5))
    min_sample_overlap = int(settings.get("min_source_sample_overlap", 3))
    min_sample_fraction = float(settings.get("min_source_sample_fraction", .5))
    tolerance = float(settings.get("vaf_coherence_tolerance", .03))
    min_coherence = float(settings.get("min_vaf_coherence", .7))
    mt_length = int(settings.get("human_mt_length", 16569))
    bin_bp = int(settings.get("dispersion_bin_bp", 1000))
    window_bp = int(settings.get("dispersion_window_bp", 1000))
    source_overlap_saturation = float(settings.get("score_source_overlap_saturation", 5.0))
    source_fraction_saturation = float(settings.get("score_source_fraction_saturation", .70))
    min_genus_species = int(settings.get("min_source_genus_species_for_specificity", 3))

    calls = {sample: alleles(path, dp_min) for sample, path in vcfs.items()}
    species_samples = defaultdict(set)
    high_index = defaultdict(list)
    for sample, rows in calls.items():
        species_samples[metadata[sample]].add(sample)
        for allele, af in rows:
            if af >= high_min:
                high_index[allele].append((sample, af))

    all_species = set(species_samples)
    report = []
    for recipient in sorted(vcfs):
        species = metadata[recipient]
        raw_low = [(key, af) for key, af in calls[recipient] if low_min <= af <= low_max]

        retained = [
            (key, af) for key, af in raw_low
            if not any(source != recipient and metadata[source] == species for source, _ in high_index.get(key, ()))
        ]

        by_species, by_sample = defaultdict(dict), defaultdict(dict)
        source_af_by_sample = defaultdict(dict)
        for key, af in retained:
            for source, source_af in high_index.get(key, ()):
                if source != recipient and metadata[source] != species:
                    by_species[metadata[source]][key] = af
                    by_sample[source][key] = af
                    source_af_by_sample[source][key] = source_af

        ranked_species = sorted(by_species, key=lambda x: (-len(by_species[x]), x))
        best_species = ranked_species[0] if ranked_species else ""
        overlap_keys = set(by_species.get(best_species, {}))
        overlap = len(overlap_keys)
        denominator = len(retained)
        fraction = overlap / denominator if denominator else 0.0

        eligible_samples = [sample for sample in by_sample if metadata[sample] == best_species]
        ranked_samples = sorted(eligible_samples, key=lambda x: (-len(by_sample[x]), x))
        best_sample = ranked_samples[0] if ranked_samples else ""
        sample_keys = set(by_sample.get(best_sample, {}))
        sample_overlap = len(sample_keys)
        sample_fraction = sample_overlap / denominator if denominator else 0.0

        values = list(by_species.get(best_species, {}).values())
        median = statistics.median(values) if values else None
        mad = median_abs_deviation(values)
        coherence = sum(abs(value - median) <= tolerance for value in values) / len(values) if values else 0.0

        # Preserve historical production classification.
        tied_species = len(ranked_species) > 1 and len(by_species[ranked_species[0]]) == len(by_species[ranked_species[1]])
        strong = overlap >= min_overlap and fraction >= min_fraction
        sample_supported = sample_overlap >= min_sample_overlap and sample_fraction >= min_sample_fraction

        if not retained:
            status, classification, reason = "PASS", "NO_INFORMATIVE_LOW_VAF", "no low-VAF alleles remain after recipient-species background removal"
        elif len(retained) < min_informative:
            status, classification, reason = "WARN", "INSUFFICIENT_INFORMATIVE_LOW_VAF", "too few low-VAF alleles remain after recipient-species background removal"
        elif not strong:
            status, classification, reason = "PASS", "NO_CROSS_SPECIES_SIGNAL", "cross-species overlap is below configured thresholds"
        elif len(species_samples[species]) == 1:
            status, classification, reason = "WARN", "SINGLETON_RECIPIENT_SPECIES", "recipient-species homoplasmic background cannot be established"
        elif tied_species:
            status, classification, reason = "WARN", "AMBIGUOUS_SOURCE_SPECIES", "multiple source species have equal best overlap"
        elif coherence < min_coherence:
            status, classification, reason = "WARN", "VAF_INCOHERENT", "matched low-VAF alleles do not meet coherence threshold"
        elif not sample_supported:
            status, classification, reason = "WARN", "INSUFFICIENT_SOURCE_SAMPLE_SUPPORT", "species-level signal is not sufficiently supported by one source sample"
        else:
            status, classification, reason = "FAIL", "INTERSPECIES_CONTAMINATION", "coherent low-VAF alleles match a different-species homoplasmic source"

        # Global specificity is measured at the species level to avoid
        # overweighting taxa with many samples. V3 additionally applies a
        # source-genus background correction when that genus contains enough
        # independently represented species in the current cohort.
        eligible_other_species = all_species - {species}
        specificity_assessable = len(eligible_other_species) >= 2
        target_genus = genus_name(species)

        def source_genus_context(source_species):
            source_genus = genus_name(source_species)
            eligible_genus_species = {
                sp for sp in eligible_other_species if genus_name(sp) == source_genus
            }
            assessable = len(eligible_genus_species) >= min_genus_species
            return source_genus, eligible_genus_species, assessable

        def effective_key_weight(key, source_species):
            carrying_species = {
                metadata[source] for source, _ in high_index.get(key, ())
                if source != recipient and metadata[source] != species
            }
            global_freq = (
                len(carrying_species) / len(eligible_other_species)
                if eligible_other_species else 1.0
            )
            global_weight = (
                specificity_weight(global_freq, settings)
                if specificity_assessable else 1.0
            )

            source_genus, eligible_genus_species, genus_assessable = source_genus_context(source_species)
            genus_carriers = carrying_species & eligible_genus_species
            genus_freq = (
                len(genus_carriers) / len(eligible_genus_species)
                if genus_assessable else None
            )
            genus_weight = (
                genus_specificity(len(genus_carriers), len(eligible_genus_species))
                if genus_assessable else 1.0
            )
            effective_weight = min(global_weight, genus_weight) if genus_assessable else global_weight
            return {
                "global_freq": global_freq,
                "global_weight": global_weight,
                "source_genus": source_genus,
                "eligible_genus_species_n": len(eligible_genus_species),
                "genus_assessable": genus_assessable,
                "genus_freq": genus_freq,
                "genus_weight": genus_weight,
                "effective_weight": effective_weight,
            }

        best_source_genus, best_source_genus_species, genus_assessable = source_genus_context(best_species)
        best_key_metrics = {
            key: effective_key_weight(key, best_species) for key in overlap_keys
        }
        key_freq = {key: m["global_freq"] for key, m in best_key_metrics.items()}
        key_global_weight = {key: m["global_weight"] for key, m in best_key_metrics.items()}
        key_weight = {key: m["effective_weight"] for key, m in best_key_metrics.items()}
        genus_freq_values = [
            m["genus_freq"] for m in best_key_metrics.values()
            if m["genus_freq"] is not None
        ]
        genus_specificity_values = [
            m["genus_weight"] for m in best_key_metrics.values()
            if m["genus_assessable"]
        ]

        adjusted_overlap = sum(key_weight.get(key, 1.0) for key in overlap_keys)
        adjusted_fraction = adjusted_overlap / denominator if denominator else 0.0
        mean_cross_species_frequency = statistics.mean(key_freq.values()) if key_freq else 0.0
        mean_specificity = statistics.mean(key_global_weight.values()) if key_global_weight else 0.0
        mean_genus_frequency = statistics.mean(genus_freq_values) if genus_freq_values else None
        mean_genus_specificity = statistics.mean(genus_specificity_values) if genus_specificity_values else None
        mean_effective_specificity = statistics.mean(key_weight.values()) if key_weight else 0.0

        best_sample_effective = sum(key_weight.get(key, 1.0) for key in sample_keys)
        sample_concentration = best_sample_effective / adjusted_overlap if adjusted_overlap > 0 else 0.0

        adjusted_by_species = {}
        for source_species, source_rows in by_species.items():
            adjusted_by_species[source_species] = sum(
                effective_key_weight(key, source_species)["effective_weight"]
                for key in source_rows
            )

        ranked_adjusted_species = sorted(
            adjusted_by_species,
            key=lambda source_species: (-adjusted_by_species[source_species], source_species),
        )
        best_adjusted_overlap = adjusted_by_species.get(best_species, 0.0)
        second_species = next(
            (source_species for source_species in ranked_adjusted_species if source_species != best_species),
            "",
        )
        second_adjusted_overlap = adjusted_by_species.get(second_species, 0.0) if second_species else 0.0
        species_separation = (
            max(0.0, best_adjusted_overlap - second_adjusted_overlap) / best_adjusted_overlap
            if best_adjusted_overlap > 0 else 0.0
        )

        dispersion = overlap_dispersion(overlap_keys, mt_length, bin_bp, window_bp)
        provenance_fields = provenance_relationship(recipient, best_sample, provenance, settings)

        score_gate = denominator >= min_informative and overlap >= min_overlap and bool(best_species)
        if score_gate:
            overlap_strength = min(1.0, adjusted_overlap / source_overlap_saturation) if source_overlap_saturation > 0 else 0.0
            fraction_strength = min(1.0, adjusted_fraction / source_fraction_saturation) if source_fraction_saturation > 0 else 0.0
            source_index = math.sqrt(max(0.0, overlap_strength * fraction_strength))
            source_pre = 4.0 * source_index
            provenance_factor = max(0.0, min(1.0, float(provenance_fields["provenance_source_factor"])))
            source_total = source_pre * provenance_factor
            score_dispersion = dispersion_points(float(dispersion["overlap_bin_entropy_normalized"]))
            score_af = af_coherence_points(mad, overlap)
            score_concentration = concentration_points(sample_concentration)
            score_separation = species_separation_points(species_separation)
            raw_10 = source_total + score_dispersion + score_af + score_concentration + score_separation
            score = raw_10 / 10.0
            score_basis = ("global_plus_source_genus_specificity_adjusted" if genus_assessable else
                           "global_cross_species_specificity_adjusted" if specificity_assessable else
                           "raw_overlap_specificity_unassessable")
        else:
            source_index = source_pre = source_total = 0.0
            score_dispersion = score_af = score_concentration = score_separation = 0.0
            raw_10 = 0.0
            score = None
            score_basis = "not_scored_gate_failed"

        row = dict(
            sample=recipient,
            species=species,
            interspecies_status=status,
            classification=classification,
            reason=reason,
            recipient_species_n=len(species_samples[species]),
            n_lowA=len(raw_low),
            n_lowA_after_species_background=denominator,
            best_source_species=best_species,
            best_source_sample=best_sample,
            overlap_count=overlap,
            overlap_fraction=f"{fraction:.6f}",
            best_source_sample_overlap=sample_overlap,
            best_source_sample_fraction=f"{sample_fraction:.6f}",
            best_source_species_overlap=overlap,
            best_source_species_fraction=f"{fraction:.6f}",
            matched_low_vaf_median="NA" if median is None else f"{median:.6f}",
            matched_low_vaf_mad="NA" if mad is None else f"{mad:.6f}",
            vaf_coherence=f"{coherence:.6f}",
            source_species_count=len(by_species),
            source_sample_count=len(by_sample),
            source_specificity_assessable="YES" if specificity_assessable else "NO",
            target_genus=target_genus,
            best_source_genus=best_source_genus,
            source_genus_specificity_assessable="YES" if genus_assessable else "NO",
            best_source_genus_species_n=len(best_source_genus_species),
            best_overlap_background_adjusted=f"{adjusted_overlap:.6f}",
            best_fraction_background_adjusted=f"{adjusted_fraction:.6f}",
            best_overlap_mean_cross_species_frequency=f"{mean_cross_species_frequency:.6f}",
            best_overlap_mean_source_specificity=f"{mean_specificity:.6f}",
            best_overlap_mean_source_genus_frequency="NA" if mean_genus_frequency is None else f"{mean_genus_frequency:.6f}",
            best_overlap_mean_source_genus_specificity="NA" if mean_genus_specificity is None else f"{mean_genus_specificity:.6f}",
            best_overlap_mean_effective_source_specificity=f"{mean_effective_specificity:.6f}",
            best_overlap_effective_specific_overlap=f"{adjusted_overlap:.6f}",
            best_source_sample_effective_overlap=f"{best_sample_effective:.6f}",
            best_source_sample_concentration=f"{sample_concentration:.6f}",
            best_source_species_adjusted_overlap=f"{best_adjusted_overlap:.6f}",
            second_source_species=second_species,
            second_source_species_adjusted_overlap=f"{second_adjusted_overlap:.6f}",
            best_source_species_separation=f"{species_separation:.6f}",
            contamination_score_gate_pass="YES" if score_gate else "NO",
            contamination_score_version="v3_genus_corrected_cross_species_project_cohort",
            contamination_score_source_basis=score_basis,
            contamination_score_source_composite_index=f"{source_index:.6f}",
            contamination_score_source_total_pre_provenance=f"{source_pre:.6f}",
            contamination_score_source_total=f"{source_total:.6f}",
            contamination_score_dispersion=f"{score_dispersion:.6f}",
            contamination_score_af_coherence=f"{score_af:.6f}",
            contamination_score_source_concentration=f"{score_concentration:.6f}",
            contamination_score_species_separation=f"{score_separation:.6f}",
            contamination_score_raw_10=f"{raw_10:.6f}",
            contamination_score="NA" if score is None else f"{score:.6f}",
            contamination_score_interpretation=score_interpretation(score),
        )
        row.update(provenance_fields)
        row.update({
            key: f"{value:.6f}" if isinstance(value, float) else value
            for key, value in dispersion.items()
        })
        report.append(row)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report)
    temporary.replace(output)
    print(f"Wrote {output} ({len(report)} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
