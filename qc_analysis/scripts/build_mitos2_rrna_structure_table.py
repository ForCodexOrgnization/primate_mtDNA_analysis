#!/usr/bin/env python3
"""Build reference-level RNA intervals and structures from final MITOS2 output.

Production structures come from ``result.mitos`` and are reconciled against
the final GFF before they are mapped to the exact coordinate FASTA.  SVG
geometry is deliberately never parsed.
"""
from __future__ import annotations

import argparse, csv, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import (
    IUPAC_DNA_BASES, IUPAC_RNA_BASES, orient_dna_base_to_rna, rrna_pair_state,
    rrna_pair_type, yaml,
)
from qc_analysis.lib.mitos_rna import (
    assignment_to_bases,
    normalize_rrna_gene as normalize_result_rrna_gene,
    parse_result_mitos,
    per_base_assignments,
    reconcile_result_mitos_record,
)

RRNA_STRUCTURE_FIELDS = (
    "reference_key reference_species coordinate_reference_accession "
    "coordinate_reference_fasta coordinate_reference_sequence_sha256 "
    "rrna_gene genomic_pos local_pos base struct_class paired_genomic_pos "
    "paired_local_pos paired_base pair_type pair_state annotation_source "
    "structure_source struct_element model_name model_position strand "
    "confidence mitos2_raw_dir"
).split()

TRNA_STRUCTURE_FIELDS = (
    "reference_key reference_species coordinate_reference_accession "
    "coordinate_reference_fasta coordinate_reference_sequence_sha256 "
    "trna_gene gene_raw genomic_pos local_pos base struct_class "
    "paired_genomic_pos paired_local_pos paired_base pair_type pair_state "
    "strand annotation_source structure_source confidence"
).split()

RRNA_REGION_FIELDS = (
    "reference_key reference_species coordinate_reference_accession "
    "coordinate_reference_fasta coordinate_reference_sequence_sha256 "
    "rrna_gene start end strand length annotation_source source_file"
).split()

PAIR_OPEN = {"<": ">", "(": ")", "[": "]", "{": "}"}
PAIR_CLOSE = {v: k for k, v in PAIR_OPEN.items()}
GAP_CHARS = {"-", "."}


def normalize_rrna_gene(gene):
    key = re.sub(r"[^A-Z0-9]", "", str(gene or "").upper())
    if key.startswith("MT"):
        key = key[2:]
    return {
        "RRNS": "MT-RNR1", "RNR1": "MT-RNR1", "12S": "MT-RNR1",
        "SRRNA": "MT-RNR1", "SSU": "MT-RNR1",
        "RRNL": "MT-RNR2", "RNR2": "MT-RNR2", "16S": "MT-RNR2",
        "LRRNA": "MT-RNR2", "LSU": "MT-RNR2",
    }.get(key, gene if gene in {"MT-RNR1", "MT-RNR2"} else str(gene or ""))


def attrs(text):
    result = {}
    for item in str(text or "").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key.lower()] = value.strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            result[key.lower()] = value.strip(' "')
    return result


def infer_rrna_gene(raw, declared=""):
    text = " ".join([str(raw or ""), str(declared or "")]).lower()
    compact = re.sub(r"[^a-z0-9]", "", text)
    if any(token in compact for token in ("rrns", "rnr1", "12s", "srna", "smallsubunit")):
        return "MT-RNR1"
    if any(token in compact for token in ("rrnl", "rnr2", "16s", "lrna", "largesubunit")):
        return "MT-RNR2"
    return normalize_rrna_gene(raw)


def parse_rrna_features_from_gff(raw_dir, ref):
    path = Path(raw_dir) / "result.gff"
    if not path.exists():
        return []
    features = []
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 9 or not columns[3].isdigit() or not columns[4].isdigit():
            continue
        declared = columns[2].lower()
        at = attrs(columns[8])
        raw = at.get("name") or at.get("gene") or at.get("gene_id") or at.get("id") or ""
        gene = infer_rrna_gene(raw, declared)
        if declared == "rrna" or gene in {"MT-RNR1", "MT-RNR2"}:
            features.append({
                **ref, "feature_type": "rRNA", "gff_seqid": columns[0],
                "gene": gene, "gene_raw": raw, "start": columns[3],
                "end": columns[4], "strand": columns[6] or "+",
            })
    return features


def read_fasta_records(path):
    records = []
    name = None
    seq = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq).upper()))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def circular_interval_coordinates(start, end, sequence_length):
    if sequence_length <= 0:
        raise ValueError("sequence_length must be positive")
    if start < 1:
        raise ValueError(f"start must be >= 1, observed {start}")
    if end < start:
        raise ValueError(f"end must be >= start, observed start={start}, end={end}")
    if end - start + 1 > sequence_length:
        raise ValueError("Circular feature length exceeds one full circular genome")
    return [(position - 1) % sequence_length for position in range(start, end + 1)]


def reference_record_for_feature(records, feature):
    seqid = feature.get("gff_seqid", "")
    if seqid:
        for name, seq in records:
            if name == seqid:
                return name, seq
    if len(records) == 1:
        return records[0]
    raise ValueError(f"No FASTA record matches rRNA GFF seqid {seqid!r}")


def feature_model(feature, records):
    name, seq = reference_record_for_feature(records, feature)
    start, end = int(feature["start"]), int(feature["end"])
    coords = circular_interval_coordinates(start, end, len(seq))
    strand = feature.get("strand", "+") or "+"
    if strand == "-":
        coords = list(reversed(coords))
    genomic_positions = [coord + 1 for coord in coords]
    genomic_bases = [seq[coord] for coord in coords]
    rna_sequence = "".join(orient_dna_base_to_rna(base, strand) or "N" for base in genomic_bases)
    return {
        "gene": normalize_rrna_gene(feature.get("gene")),
        "strand": strand,
        "record_name": name,
        "genomic_positions": genomic_positions,
        "genomic_bases": genomic_bases,
        "rna_sequence": rna_sequence,
    }


def text_file(path):
    if path.suffix.lower() in {".fa", ".fasta", ".fna", ".gz", ".bam", ".png", ".pdf", ".svg"}:
        return False
    try:
        path.read_text(errors="strict")
        return True
    except (UnicodeDecodeError, OSError):
        return False


def stockholm_candidates(raw_dir):
    raw = Path(raw_dir)
    if not raw.exists():
        return []
    candidates = []
    for path in sorted(raw.rglob("*")):
        if not path.is_file() or not text_file(path):
            continue
        try:
            snippet = path.read_text(errors="replace")[:200000]
        except OSError:
            continue
        if "# STOCKHOLM" in snippet or "#=GC SS_cons" in snippet or re.search(r"^#=GR\s+\S+\s+SS\s+", snippet, re.M):
            candidates.append(path)
    return candidates


def parse_stockholm(path):
    blocks = []
    seqs, gr_ss, ss_cons = {}, {}, []

    def finish():
        if not seqs:
            return
        consensus = "".join(ss_cons)
        for name, parts in seqs.items():
            structure = "".join(gr_ss.get(name, [])) or consensus
            if structure:
                blocks.append({
                    "path": str(path), "sequence_name": name,
                    "sequence": "".join(parts), "structure": structure,
                })

    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("# STOCKHOLM"):
            if seqs:
                finish()
            seqs, gr_ss, ss_cons = {}, {}, []
            continue
        if line == "//":
            finish()
            seqs, gr_ss, ss_cons = {}, {}, []
            continue
        if line.startswith("#=GC SS_cons"):
            parts = line.split(None, 2)
            if len(parts) == 3:
                ss_cons.append(parts[2].strip())
            continue
        if line.startswith("#=GR"):
            parts = line.split(None, 3)
            if len(parts) == 4 and parts[2] == "SS":
                gr_ss.setdefault(parts[1], []).append(parts[3].strip())
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            seqs.setdefault(parts[0], []).append(parts[1].strip())
    if seqs:
        finish()
    return blocks


def structure_pairs(structure):
    pairs = {}
    stacks = {char: [] for char in PAIR_OPEN}
    for index, char in enumerate(structure):
        if char in PAIR_OPEN:
            stacks[char].append(index)
        elif char in PAIR_CLOSE:
            opener = PAIR_CLOSE[char]
            if stacks[opener]:
                left = stacks[opener].pop()
                pairs[left] = index
                pairs[index] = left
    return pairs


def ungapped_sequence(sequence):
    return "".join(char.upper().replace("T", "U") for char in sequence if char not in GAP_CHARS)


def alignment_column_to_local(sequence):
    mapping = {}
    local = 0
    for column, char in enumerate(sequence):
        if char in GAP_CHARS:
            continue
        if char.upper().replace("T", "U") in IUPAC_RNA_BASES or char.upper() in IUPAC_DNA_BASES:
            local += 1
            mapping[column] = local
    return mapping, local


def infer_gene_for_structure(record, models):
    text = f"{Path(record['path']).name} {record['sequence_name']}"
    gene = infer_rrna_gene(text)
    if gene in models:
        return gene
    seq = ungapped_sequence(record["sequence"])
    exact = [gene for gene, model in models.items() if seq == model["rna_sequence"]]
    if len(exact) == 1:
        return exact[0]
    by_length = [gene for gene, model in models.items() if len(seq) == len(model["rna_sequence"])]
    return by_length[0] if len(by_length) == 1 else ""


def assignments_from_record(record, model):
    sequence = record["sequence"]
    structure = record["structure"]
    if len(sequence) != len(structure):
        return {}, f"structure length {len(structure)} != sequence length {len(sequence)}"
    col_to_local, n_residues = alignment_column_to_local(sequence)
    if n_residues != len(model["genomic_positions"]):
        return {}, f"structure residues {n_residues} != rRNA feature length {len(model['genomic_positions'])}"
    pairs = structure_pairs(structure)
    assignments = {}
    for column, local_pos in col_to_local.items():
        char = structure[column]
        partner_column = pairs.get(column)
        if partner_column is not None and partner_column in col_to_local:
            assignments[local_pos] = {
                "struct_class": "stem",
                "paired_local_pos": col_to_local[partner_column],
                "structure_source": record["path"],
                "model_name": record["sequence_name"],
                "model_position": local_pos,
            }
        elif char in PAIR_OPEN or char in PAIR_CLOSE:
            assignments[local_pos] = {
                "struct_class": "unknown",
                "structure_source": record["path"],
                "model_name": record["sequence_name"],
                "model_position": local_pos,
            }
        else:
            assignments[local_pos] = {
                "struct_class": "loop",
                "structure_source": record["path"],
                "model_name": record["sequence_name"],
                "model_position": local_pos,
            }
    return assignments, ""


def base_row(ref, raw_dir, model, local_pos, assignment):
    index = local_pos - 1
    genomic_pos = model["genomic_positions"][index]
    base = model["genomic_bases"][index]
    strand = model["strand"]
    klass = assignment.get("struct_class", "unknown")
    paired_local = assignment.get("paired_local_pos", ".")
    paired_genomic = paired_base = pair_kind = "."
    state = rrna_pair_state(".", klass)
    if klass == "stem" and isinstance(paired_local, int):
        partner_index = paired_local - 1
        paired_genomic = model["genomic_positions"][partner_index]
        paired_base = model["genomic_bases"][partner_index]
        pair_kind = rrna_pair_type(orient_dna_base_to_rna(base, strand), orient_dna_base_to_rna(paired_base, strand))
        state = rrna_pair_state(pair_kind, klass)
    elif klass == "loop":
        state = "unpaired"
    return {
        "reference_key": ref.get("reference_key", ""),
        "reference_species": ref.get("reference_species", ""),
        "coordinate_reference_accession": ref.get("coordinate_reference_accession", ""),
        "coordinate_reference_fasta": ref.get("coordinate_reference_fasta", ""),
        "coordinate_reference_sequence_sha256": ref.get("coordinate_reference_sequence_sha256", ""),
        "rrna_gene": model["gene"],
        "genomic_pos": genomic_pos,
        "local_pos": local_pos,
        "base": base,
        "struct_class": klass,
        "paired_genomic_pos": paired_genomic,
        "paired_local_pos": paired_local if paired_local != "." else ".",
        "paired_base": paired_base,
        "pair_type": pair_kind,
        "pair_state": state,
        "annotation_source": "MITOS2",
        "structure_source": assignment.get("structure_source", "."),
        "struct_element": ".",
        "model_name": assignment.get("model_name", "."),
        "model_position": assignment.get("model_position", "."),
        "strand": strand,
        "confidence": ".",
        "mitos2_raw_dir": str(raw_dir),
    }


def build_reference_rrna_structure_rows(ref, features, fasta, raw_dir):
    """Build final rRNA rows, preferring validated ``result.mitos`` records."""
    result_path = Path(raw_dir) / "result.mitos"
    if result_path.is_file():
        return build_reference_rna_structure_rows(ref, features, fasta, raw_dir, "rRNA")

    # Explicitly legacy-only support for old retained Stockholm fixtures.  It
    # is unreachable when a final result.mitos exists and can never override it.
    rrna_features = [f for f in features if f.get("feature_type") == "rRNA" or normalize_rrna_gene(f.get("gene")) in {"MT-RNR1", "MT-RNR2"}]
    if not rrna_features:
        rrna_features = parse_rrna_features_from_gff(raw_dir, ref)
    if not rrna_features:
        return [], "no_rrna_features", "No MITOS2 rRNA feature intervals were available."
    records = read_fasta_records(fasta)
    models = {}
    notes = []
    for feature in rrna_features:
        try:
            model = feature_model(feature, records)
            if model["gene"] in {"MT-RNR1", "MT-RNR2"}:
                models[model["gene"]] = model
        except Exception as exc:
            notes.append(f"{feature.get('gene', 'rRNA')}: {type(exc).__name__}: {exc}")
    if not models:
        return [], "failed_rrna_feature_mapping", "; ".join(notes) or "No rRNA features mapped to the coordinate FASTA."

    assignments = {gene: {} for gene in models}
    candidates = stockholm_candidates(raw_dir)
    parse_notes = []
    for path in candidates:
        for record in parse_stockholm(path):
            gene = infer_gene_for_structure(record, models)
            if gene not in models or assignments[gene]:
                continue
            mapped, note = assignments_from_record(record, models[gene])
            if mapped:
                assignments[gene] = mapped
            elif note:
                parse_notes.append(f"{Path(path).name}: {note}")

    rows = []
    for gene, model in sorted(models.items()):
        for local_pos in range(1, len(model["genomic_positions"]) + 1):
            assignment = assignments[gene].get(local_pos, {"struct_class": "unknown"})
            rows.append(base_row(ref, raw_dir, model, local_pos, assignment))

    parsed_genes = [gene for gene, mapped in assignments.items() if mapped]
    if parsed_genes and len(parsed_genes) == len(models):
        return rows, "parsed_legacy_stockholm_structure", "Legacy fallback: parsed Stockholm/Infernal secondary structure for all MITOS2 rRNA features."
    if parsed_genes:
        missing = sorted(set(models) - set(parsed_genes))
        return rows, "partial_legacy_stockholm_structure", "Legacy fallback: parsed Stockholm/Infernal secondary structure for " + ",".join(sorted(parsed_genes)) + "; missing " + ",".join(missing)
    if candidates:
        return rows, "no_matched_legacy_stockholm_structure", "; ".join(parse_notes) or "Legacy Stockholm structure files were present but did not match rRNA features."
    return rows, "no_result_mitos_rna_structure", "No final result.mitos RNA structure was found; SVG plots were not parsed."


def _rna_model(feature, records):
    model = feature_model(feature, records)
    if feature.get("feature_type") == "tRNA":
        model["gene"] = str(feature.get("gene_raw") or feature.get("gene") or "")
    return model


def _rna_base_row(ref, feature, model, local_pos, assignment, feature_type, confidence="."):
    row = {
        "reference_key": ref.get("reference_key", ""),
        "reference_species": ref.get("reference_species", ""),
        "coordinate_reference_accession": ref.get("coordinate_reference_accession", ""),
        "coordinate_reference_fasta": ref.get("coordinate_reference_fasta", ""),
        "coordinate_reference_sequence_sha256": ref.get("coordinate_reference_sequence_sha256", ""),
        **assignment_to_bases(model, local_pos, assignment),
        "annotation_source": "MITOS2",
        "structure_source": assignment.get("structure_source", "."),
        "confidence": confidence or ".",
    }
    if feature_type == "rRNA":
        row.update({
            "rrna_gene": normalize_result_rrna_gene(feature.get("gene_raw") or feature.get("gene")),
            "struct_element": ".", "model_name": ".", "model_position": local_pos,
            "mitos2_raw_dir": str(feature.get("raw_dir", "")),
        })
    else:
        raw_gene = str(feature.get("mitos_gene_raw") or feature.get("gene_raw") or feature.get("gene") or "")
        row.update({"trna_gene": re.sub(r"\([^)]*\)$", "", raw_gene), "gene_raw": raw_gene})
    return row


def build_reference_rna_structure_rows(ref, features, fasta, raw_dir, feature_type):
    """Expand final MITOS RNA structures in RNA 5'->3' local orientation."""
    wanted = feature_type.lower()
    final_features = [f for f in features if str(f.get("feature_type", "")).lower() == wanted]
    if feature_type == "rRNA" and not final_features:
        final_features = parse_rrna_features_from_gff(raw_dir, ref)
    noun = "rrna" if feature_type == "rRNA" else "trna"
    if not final_features:
        return [], f"no_mitos2_{noun}_features", f"No final MITOS2 {feature_type} GFF features were available."

    records = read_fasta_records(fasta)
    mitos_records = [r for r in parse_result_mitos(Path(raw_dir) / "result.mitos") if r["feature_type"] == feature_type]
    reconciled = []
    interval_mismatches = []
    for record in mitos_records:
        match, status, note = reconcile_result_mitos_record(record, final_features)
        reconciled.append((record, match, status, note))
        if status == "result_mitos_gff_interval_mismatch":
            interval_mismatches.append(f"{record['gene_raw']}: {note}")

    rows, notes = [], []
    parsed_features = 0
    length_mismatches = 0
    for feature in final_features:
        try:
            model = _rna_model(feature, records)
        except Exception as exc:
            notes.append(f"{feature.get('gene_raw', feature.get('gene', feature_type))}: FASTA mapping failed: {exc}")
            continue
        matches = [(record, status) for record, match, status, _note in reconciled if match is feature]
        assignment_map = {}
        confidence = "."
        mitos_gene_raw = ""
        if len(matches) == 1 and matches[0][0].get("structure"):
            record = matches[0][0]
            mitos_gene_raw = record.get("gene_raw", "")
            try:
                assignment_map = per_base_assignments(
                    record["structure"], len(model["genomic_positions"]), record["source_file"]
                )
                parsed_features += 1
                confidence = record.get("score") or "."
            except ValueError as exc:
                if "structure length" in str(exc):
                    length_mismatches += 1
                notes.append(f"{record['gene_raw']}: {exc}")
        elif len(matches) > 1:
            notes.append(f"{feature.get('gene_raw')}: multiple result.mitos records match the final GFF feature")
        else:
            notes.append(f"{feature.get('gene_raw')}: no usable final result.mitos structure")
        for local_pos in range(1, len(model["genomic_positions"]) + 1):
            assignment = assignment_map.get(local_pos, {"struct_class": "unknown", "structure_source": "."})
            feature_with_raw = {**feature, "raw_dir": str(raw_dir), "mitos_gene_raw": mitos_gene_raw}
            rows.append(_rna_base_row(ref, feature_with_raw, model, local_pos, assignment, feature_type, confidence))

    total = len(final_features)
    if parsed_features == total:
        status = "parsed_result_mitos_structure"
        note = f"Parsed and GFF-validated final result.mitos structure for all {total} {feature_type} features."
    elif parsed_features:
        status = "partial_result_mitos_structure"
        note = f"Parsed {parsed_features} of {total} final {feature_type} structures; " + "; ".join(notes)
    elif interval_mismatches:
        status = "result_mitos_gff_interval_mismatch"
        note = "; ".join(interval_mismatches + notes)
    elif length_mismatches:
        status = "result_mitos_structure_length_mismatch"
        note = "; ".join(notes)
    else:
        status = f"no_result_mitos_{noun}_structure"
        note = "; ".join(notes) or f"No usable final result.mitos {feature_type} structure was found."
    return rows, status, note


def build_reference_trna_structure_rows(ref, features, fasta, raw_dir):
    return build_reference_rna_structure_rows(ref, features, fasta, raw_dir, "tRNA")


def build_reference_rrna_region_rows(ref, features):
    rows = []
    for feature in features:
        if str(feature.get("feature_type", "")).lower() != "rrna":
            continue
        start, end = int(feature["start"]), int(feature["end"])
        rows.append({
            "reference_key": ref.get("reference_key", ""),
            "reference_species": ref.get("reference_species", ""),
            "coordinate_reference_accession": ref.get("coordinate_reference_accession", ""),
            "coordinate_reference_fasta": ref.get("coordinate_reference_fasta", ""),
            "coordinate_reference_sequence_sha256": ref.get("coordinate_reference_sequence_sha256", ""),
            "rrna_gene": normalize_result_rrna_gene(feature.get("gene_raw") or feature.get("gene")),
            "start": start, "end": end, "strand": feature.get("strand", "+"),
            "length": end - start + 1, "annotation_source": "MITOS2",
            "source_file": feature.get("source_file", ""),
        })
    return rows


def write_tsv(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_from_config(config_path, sample=None):
    from qc_analysis.scripts import run_mitos2_annotation as mitos2

    cfg = yaml(config_path)
    section = cfg.get("mitos2_annotation") or {}
    paths = section.get("paths", {})
    refs = mitos2.references(paths, sample)
    output = paths.get(
        "mitos2_reference_rrna_structure_table",
        str(Path(paths.get("output_dir", "results/qc/mitos2_annotation")) / "all_mitos2_reference_rrna_structure.tsv"),
    )
    all_rows = []
    for ref, _linked in refs:
        raw_dir = Path(paths["mitos2_raw_dir"]) / ref.get("task_key", ref["reference_key"])
        features, _diag = mitos2.parse_outputs(raw_dir, ref)
        rows, _status, _note = build_reference_rrna_structure_rows(ref, features, ref.get("coordinate_reference_fasta", ""), raw_dir)
        all_rows.extend(rows)
    write_tsv(output, RRNA_STRUCTURE_FIELDS, all_rows)
    return output, len(all_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample")
    args = parser.parse_args()
    output, n_rows = build_from_config(args.config, args.sample)
    print(f"Wrote {n_rows} reference rRNA structure rows to {output}.")


if __name__ == "__main__":
    main()
