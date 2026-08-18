"""Parsing and coordinate-safe expansion of final MITOS2 RNA structures."""
from __future__ import annotations

import re
from pathlib import Path

from qc_analysis.lib.match_utils import orient_dna_base_to_rna, rrna_pair_state, rrna_pair_type


PAIR_OPEN = {"(": ")", "[": "]", "{": "}", "<": ">"}
PAIR_CLOSE = {closer: opener for opener, closer in PAIR_OPEN.items()}
STRUCTURE_RE = re.compile(r"^[.()\[\]{}<>]+$")


def normalize_strand(value):
    text = str(value or "").strip().lower()
    if text in {"1", "+", "+1", "plus", "forward"}:
        return "+"
    if text in {"-1", "-", "minus", "reverse"}:
        return "-"
    return str(value or "").strip()


def normalize_rrna_gene(gene):
    key = re.sub(r"[^A-Z0-9]", "", str(gene or "").upper())
    if key.startswith("MT"):
        key = key[2:]
    return {
        "RRNS": "MT-RNR1", "RNR1": "MT-RNR1", "12S": "MT-RNR1",
        "SRRNA": "MT-RNR1", "SSU": "MT-RNR1",
        "RRNL": "MT-RNR2", "RNR2": "MT-RNR2", "16S": "MT-RNR2",
        "LRRNA": "MT-RNR2", "LSU": "MT-RNR2",
    }.get(key, str(gene or ""))


def comparable_gene(gene, feature_type):
    if str(feature_type).lower() == "rrna":
        return normalize_rrna_gene(gene)
    # Parenthesized anticodons are presentation details, but L1/L2 and S1/S2
    # remain part of the identity and therefore must not be collapsed.
    return re.sub(r"\([^)]*\)$", "", str(gene or "").strip(), flags=re.I).lower()


def _structure_field(columns):
    """Return the RNA structure from the MITOS-specific tail columns.

    MITOS2's custom format has fixed core fields (seqid, type, name, method,
    start, end, strand, score), followed by feature-specific fields.  The
    structure is identified only inside that tail and must use MITOS-supported
    dot/bracket symbols.  Choosing the longest candidate handles the optional
    tRNA anticodon/model columns without pretending this is BED or GFF.
    """
    candidates = [
        value.strip() for value in columns[8:]
        if len(value.strip()) > 1 and STRUCTURE_RE.fullmatch(value.strip())
    ]
    return max(candidates, key=len) if candidates else ""


def parse_result_mitos(path):
    """Parse final tRNA/rRNA records from MITOS2 ``result.mitos``.

    Raw MITOS coordinates are retained.  The documented MITOS convention is
    normalized explicitly from zero-based inclusive to one-based inclusive;
    callers must still validate that interval against final ``result.gff``.
    """
    path = Path(path)
    if not path.is_file():
        return []
    records = []
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) < 8 or columns[1].strip().lower() not in {"trna", "rrna"}:
            continue
        try:
            start, end = int(columns[4]), int(columns[5])
        except ValueError:
            continue
        feature_type = "tRNA" if columns[1].strip().lower() == "trna" else "rRNA"
        gene_raw = columns[2].strip()
        records.append({
            "feature_type": feature_type,
            "gene": normalize_rrna_gene(gene_raw) if feature_type == "rRNA" else gene_raw,
            "gene_raw": gene_raw,
            "mitos_start": start,
            "mitos_end": end,
            "normalized_start": start + 1,
            "normalized_end": end + 1,
            "strand": normalize_strand(columns[6]),
            "score": columns[7].strip(),
            "structure": _structure_field(columns),
            "source_file": str(path),
            "line_number": line_number,
        })
    return records


def reconcile_result_mitos_record(record, gff_features):
    """Match one MITOS record to the same final GFF RNA feature."""
    target_gene = comparable_gene(record.get("gene_raw"), record.get("feature_type"))
    candidates = [
        feature for feature in gff_features
        if str(feature.get("feature_type", "")).lower() == str(record.get("feature_type", "")).lower()
        and comparable_gene(feature.get("gene_raw") or feature.get("gene"), feature.get("feature_type")) == target_gene
    ]
    strand_candidates = [f for f in candidates if normalize_strand(f.get("strand")) == record.get("strand")]
    exact = [
        f for f in strand_candidates
        if int(f["start"]) == record["normalized_start"]
        and int(f["end"]) == record["normalized_end"]
    ]
    if len(exact) == 1:
        return exact[0], "matched_final_gff", ""
    if len(exact) > 1:
        return None, "ambiguous_final_gff_match", f"{len(exact)} identical final GFF matches"
    if candidates:
        observed = ",".join(f"{f.get('start')}-{f.get('end')}({normalize_strand(f.get('strand'))})" for f in candidates)
        expected = f"{record['normalized_start']}-{record['normalized_end']}({record['strand']})"
        return None, "result_mitos_gff_interval_mismatch", f"normalized result.mitos={expected}; final GFF={observed}"
    return None, "no_final_gff_feature_match", f"no final GFF match for {record.get('feature_type')} {record.get('gene_raw')}"


def structure_pairs(structure):
    """Parse balanced dot-bracket notation and return reciprocal zero-based pairs."""
    if not STRUCTURE_RE.fullmatch(structure or ""):
        raise ValueError("structure contains unsupported symbols")
    pairs = {}
    stacks = {opener: [] for opener in PAIR_OPEN}
    for index, symbol in enumerate(structure):
        if symbol in PAIR_OPEN:
            stacks[symbol].append(index)
        elif symbol in PAIR_CLOSE:
            opener = PAIR_CLOSE[symbol]
            if not stacks[opener]:
                raise ValueError(f"unmatched closing bracket {symbol!r} at local position {index + 1}")
            partner = stacks[opener].pop()
            pairs[index] = partner
            pairs[partner] = index
    unmatched = [(symbol, index + 1) for symbol, stack in stacks.items() for index in stack]
    if unmatched:
        symbol, position = unmatched[0]
        raise ValueError(f"unmatched opening bracket {symbol!r} at local position {position}")
    if any(pairs.get(partner) != position for position, partner in pairs.items()):
        raise ValueError("dot-bracket pair relationships are not reciprocal")
    return pairs


def per_base_assignments(structure, feature_length, source_file):
    if len(structure) != feature_length:
        raise ValueError(f"structure length {len(structure)} != RNA feature length {feature_length}")
    pairs = structure_pairs(structure)
    assignments = {}
    for index, symbol in enumerate(structure):
        partner = pairs.get(index)
        assignments[index + 1] = {
            "struct_class": "stem" if partner is not None else "loop",
            "paired_local_pos": partner + 1 if partner is not None else ".",
            "structure_source": str(source_file),
        }
    return assignments


def assignment_to_bases(model, local_pos, assignment):
    index = local_pos - 1
    genomic_pos = model["genomic_positions"][index]
    base = model["genomic_bases"][index]
    klass = assignment.get("struct_class", "unknown")
    paired_local = assignment.get("paired_local_pos", ".")
    paired_genomic = paired_base = pair_kind = "."
    state = rrna_pair_state(".", klass)
    if klass == "stem" and isinstance(paired_local, int):
        partner_index = paired_local - 1
        if partner_index < 0 or partner_index >= len(model["genomic_positions"]):
            raise ValueError(f"paired local position {paired_local} is outside the RNA feature")
        paired_genomic = model["genomic_positions"][partner_index]
        paired_base = model["genomic_bases"][partner_index]
        pair_kind = rrna_pair_type(
            orient_dna_base_to_rna(base, model["strand"]),
            orient_dna_base_to_rna(paired_base, model["strand"]),
        )
        state = rrna_pair_state(pair_kind, klass)
    elif klass == "loop":
        state = "unpaired"
    return {
        "genomic_pos": genomic_pos, "local_pos": local_pos, "base": base,
        "struct_class": klass, "paired_genomic_pos": paired_genomic,
        "paired_local_pos": paired_local, "paired_base": paired_base,
        "pair_type": pair_kind, "pair_state": state,
        "strand": model["strand"],
    }
