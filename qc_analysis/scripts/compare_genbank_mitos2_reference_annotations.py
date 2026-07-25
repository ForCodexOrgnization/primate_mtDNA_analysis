#!/usr/bin/env python3
"""Compare GenBank and MITOS2 annotations by biological sequence identity."""
from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256

GENES = ('MT-ND1 MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-ND5 MT-ND6 MT-CO1 '
         'MT-CO2 MT-CO3 MT-CYB MT-ATP6 MT-ATP8').split()
GROUP_FIELDS = '''comparison_sequence_sha256 genbank_reference_keys mitos2_reference_keys genbank_accessions mitos2_accessions mitos2_target_species n_genbank_groups n_mitos2_groups sequence_usage_category'''.split()
GENE_FIELDS = GROUP_FIELDS + '''genbank_reference_key mitos2_reference_key gene genbank_present mitos2_present genbank_n_rows mitos2_n_rows genbank_n_unique_positions mitos2_n_unique_positions genbank_start genbank_end mitos2_start mitos2_end genbank_strand mitos2_strand position_overlap position_union position_jaccard genbank_only_positions mitos2_only_positions same_position_set same_strand same_codon_triplets same_codon_position_mapping same_reference_bases start_delta end_delta length_delta wraps_origin_genbank wraps_origin_mitos2 ordered_coordinate_match sequence_compatibility_category gene_comparison_category'''.split()
SUMMARY_FIELDS = GROUP_FIELDS + '''genbank_reference_key mitos2_reference_key target_species reference_species genbank_accession mitos2_accession genbank_record_sequence_sha256 mitos2_input_sequence_sha256 sequence_compatibility_category sequence_match coordinate_comparison_performed n_genes_genbank n_genes_mitos2 n_exact_genes n_minor_difference_genes n_moderate_difference_genes n_major_difference_genes n_missing_genes n_strand_mismatches all_13_genes_present_genbank all_13_genes_present_mitos2 all_13_exact_match reference_comparison_category note'''.split()
DIAG_FIELDS = GROUP_FIELDS + '''genbank_reference_key mitos2_reference_key target_species reference_species genbank_accession mitos2_accession genbank_record_sequence_sha256 mitos2_input_sequence_sha256 sequence_compatibility_category candidate_matching_basis reason_comparison_skipped'''.split()
SIGNATURE_FIELDS = ('gene', 'pos', 'strand', 'codon_index', 'codon_pos_in_triplet',
                    'codon_seq', 'ref_base_genome')


def v(row, key):
    return (row.get(key) or '').strip()


def yes(value):
    return 'yes' if value else 'no'


def read_groups(path, mitos=False):
    """Read annotation groups, deriving only a display key for legacy MITOS2."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Required reference codon table is missing: {path}')
    groups = defaultdict(list)
    with path.open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = v(row, 'reference_key')
            if not key and mitos:
                key = v(row, 'target_species')
            if not key:
                key = v(row, 'coordinate_reference_accession') or v(row, 'accession_version')
            row['reference_key'] = key
            groups[key].append(row)
    return groups


def group_hash(rows, mitos=False):
    """Return the comparison hash; GenBank deliberately has no FASTA fallback."""
    if not rows:
        return ''
    if not mitos:
        return v(rows[0], 'genbank_record_sequence_sha256')
    existing = v(rows[0], 'mitos2_input_sequence_sha256')
    if existing:
        return existing
    path = v(rows[0], 'mitos2_input_fasta') or v(rows[0], 'coordinate_reference_fasta')
    if not path:
        return ''
    try:
        return normalized_fasta_sequence_sha256(path)['sequence_sha256']
    except (OSError, ValueError):
        return ''


def annotation_signature(rows):
    return tuple(sorted(tuple(v(row, field) for field in SIGNATURE_FIELDS) for row in rows))


def accession(rows):
    row = rows[0]
    return v(row, 'coordinate_reference_accession') or v(row, 'accession_version') or v(row, 'accession')


def joined(values):
    return ';'.join(sorted({value for value in values if value}))


def ordered(rows):
    return [(v(r, 'pos'), v(r, 'codon_index'), v(r, 'codon_pos_in_triplet'), v(r, 'strand'))
            for r in sorted(rows, key=lambda x: (int(v(x, 'codon_index') or 0),
                                                  int(v(x, 'codon_pos_in_triplet') or 0),
                                                  int(v(x, 'pos') or 0)))]


def group_metadata(sequence_hash, genbank_groups, mitos_groups):
    gkeys = list(genbank_groups)
    mkeys = list(mitos_groups)
    targets = [v(rows[0], 'target_species') for rows in mitos_groups.values()]
    return {
        'comparison_sequence_sha256': sequence_hash,
        'genbank_reference_keys': joined(gkeys),
        'mitos2_reference_keys': joined(mkeys),
        'genbank_accessions': joined(accession(rows) for rows in genbank_groups.values()),
        'mitos2_accessions': joined(accession(rows) for rows in mitos_groups.values()),
        'mitos2_target_species': joined(targets),
        'n_genbank_groups': len(genbank_groups),
        'n_mitos2_groups': len(mitos_groups),
        'sequence_usage_category': ('shared_sequence_multiple_targets'
                                    if len(set(filter(None, targets))) > 1
                                    else 'shared_sequence_single_target'),
    }


def gene_row(g_rows, m_rows, gene, compatibility, minor, moderate, metadata):
    a = [r for r in g_rows if v(r, 'gene') == gene]
    b = [r for r in m_rows if v(r, 'gene') == gene]
    ap, bp = {v(r, 'pos') for r in a}, {v(r, 'pos') for r in b}
    strand_a, strand_b = {v(r, 'strand') for r in a}, {v(r, 'strand') for r in b}
    overlap, union = len(ap & bp), len(ap | bp)
    jac = overlap / union if union else 1.0
    wrapa = bool(a and ordered(a)[0][0] != min(ap, key=int))
    wrapb = bool(b and ordered(b)[0][0] != min(bp, key=int))
    bases = lambda xs: {v(x, 'pos'): v(x, 'ref_base_genome') for x in xs}
    triplets = lambda xs: [(v(x, 'codon_index'), v(x, 'codon_pos_in_triplet'), v(x, 'codon_seq')) for x in sorted(xs, key=lambda x: (int(v(x, 'codon_index') or 0), int(v(x, 'codon_pos_in_triplet') or 0)))]
    mapping = lambda xs: [(v(x, 'pos'), v(x, 'codon_index'), v(x, 'codon_pos_in_triplet')) for x in xs]
    samepos = ap == bp
    samestrand = strand_a == strand_b and len(strand_a) == 1
    ordermatch = ordered(a) == ordered(b)
    if not a:
        kind = 'missing_in_genbank'
    elif not b:
        kind = 'missing_in_mitos2'
    elif not samestrand:
        kind = 'strand_mismatch'
    elif samepos and ordermatch and triplets(a) == triplets(b) and bases(a) == bases(b):
        kind = 'exact_match'
    elif samepos:
        kind = 'coordinate_match_codon_difference'
    elif jac >= minor:
        kind = 'minor_boundary_difference'
    elif jac >= moderate:
        kind = 'moderate_boundary_difference'
    else:
        kind = 'major_coordinate_difference'
    bounds = lambda xs: (min((int(v(r, 'pos')) for r in xs), default=''), max((int(v(r, 'pos')) for r in xs), default=''))
    sa, ea = bounds(a)
    sb, eb = bounds(b)
    values = [v(g_rows[0], 'reference_key'), v(m_rows[0], 'reference_key'), gene, yes(a), yes(b), len(a), len(b), len(ap), len(bp), sa, ea, sb, eb, joined(strand_a), joined(strand_b), overlap, union, f'{jac:.12g}', joined(sorted(ap-bp, key=int)), joined(sorted(bp-ap, key=int)), yes(samepos), yes(samestrand), yes(triplets(a) == triplets(b)), yes(mapping(a) == mapping(b)), yes(bases(a) == bases(b)), (sb-sa if a and b else ''), (eb-ea if a and b else ''), len(bp)-len(ap), yes(wrapa), yes(wrapb), yes(ordermatch), compatibility, kind]
    return {**metadata, **dict(zip(GENE_FIELDS[len(GROUP_FIELDS):], values))}


def atomic(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name + '.', text=True)
    with os.fdopen(fd, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def compare(genbank_table, mitos_table, output, summary_output, mismatch_output,
            strict=True, allow_rotation=False, fail_no_shared=True, minor=.99,
            moderate=.90, sample_reference_map=None, species_fasta_dir='',
            fasta_extensions=''):
    """Compare once per shared normalized sequence hash.

    Deprecated resolution and rotation arguments remain in the API so existing
    configurations continue to run, but biological pairing is hash-only.
    """
    del strict, allow_rotation, sample_reference_map, species_fasta_dir, fasta_extensions
    gb = read_groups(genbank_table)
    mi = read_groups(mitos_table, mitos=True)
    gb_by_hash, mi_by_hash = defaultdict(dict), defaultdict(dict)
    for key, rows in gb.items():
        gb_by_hash[group_hash(rows)][key] = rows
    for key, rows in mi.items():
        mi_by_hash[group_hash(rows, True)][key] = rows
    gb_hashes, mi_hashes = set(gb_by_hash) - {''}, set(mi_by_hash) - {''}
    shared = gb_hashes & mi_hashes
    diagnostics, summaries, gene_rows = [], [], []

    def diagnostic(sequence_hash, ggroups, mgroups, reason):
        metadata = group_metadata(sequence_hash, ggroups, mgroups)
        if not ggroups or not mgroups:
            metadata['sequence_usage_category'] = 'unmatched_sequence'
        grow = next(iter(ggroups.values()), [{}])
        mrow = next(iter(mgroups.values()), [{}])
        diagnostics.append({**metadata,
            'genbank_reference_key': v(grow[0], 'reference_key'),
            'mitos2_reference_key': v(mrow[0], 'reference_key'),
            'target_species': v(mrow[0], 'target_species'),
            'reference_species': v(grow[0], 'reference_species'),
            'genbank_accession': accession(grow) if ggroups else '',
            'mitos2_accession': accession(mrow) if mgroups else '',
            'genbank_record_sequence_sha256': sequence_hash if ggroups else '',
            'mitos2_input_sequence_sha256': sequence_hash if mgroups else '',
            'sequence_compatibility_category': reason,
            'candidate_matching_basis': 'normalized_sequence_sha256',
            'reason_comparison_skipped': reason})

    for sequence_hash in sorted(gb_hashes - mi_hashes):
        diagnostic(sequence_hash, gb_by_hash[sequence_hash], {}, 'no_matching_mitos2_sequence_hash')
    for sequence_hash in sorted(mi_hashes - gb_hashes):
        diagnostic(sequence_hash, {}, mi_by_hash[sequence_hash], 'no_matching_genbank_sequence_hash')

    inconsistent = 0
    single_target = multiple_target = 0
    for sequence_hash in sorted(shared):
        ggroups, mgroups = gb_by_hash[sequence_hash], mi_by_hash[sequence_hash]
        if len({annotation_signature(rows) for rows in ggroups.values()}) > 1:
            diagnostic(sequence_hash, ggroups, mgroups, 'same_sequence_inconsistent_genbank_annotation')
            inconsistent += 1
            continue
        if len({annotation_signature(rows) for rows in mgroups.values()}) > 1:
            diagnostic(sequence_hash, ggroups, mgroups, 'same_sequence_inconsistent_mitos2_annotation')
            inconsistent += 1
            continue
        metadata = group_metadata(sequence_hash, ggroups, mgroups)
        if metadata['sequence_usage_category'] == 'shared_sequence_multiple_targets':
            multiple_target += 1
        else:
            single_target += 1
        g_rows, m_rows = next(iter(ggroups.values())), next(iter(mgroups.values()))
        rows = [gene_row(g_rows, m_rows, gene, 'exact_sequence_match', minor, moderate, metadata) for gene in GENES]
        gene_rows.extend(rows)
        kinds = [row['gene_comparison_category'] for row in rows]
        category = ('all_13_exact' if all(k == 'exact_match' for k in kinds) else
                    'strand_mismatch' if 'strand_mismatch' in kinds else
                    'gene_missing' if any(k.startswith('missing') for k in kinds) else
                    'only_minor_boundary_differences' if all(k in ('exact_match', 'minor_boundary_difference') for k in kinds) else
                    'moderate_or_major_differences')
        summaries.append({**metadata,
            'genbank_reference_key': v(g_rows[0], 'reference_key'), 'mitos2_reference_key': v(m_rows[0], 'reference_key'),
            'target_species': v(m_rows[0], 'target_species'), 'reference_species': v(g_rows[0], 'reference_species'),
            'genbank_accession': accession(g_rows), 'mitos2_accession': accession(m_rows),
            'genbank_record_sequence_sha256': sequence_hash, 'mitos2_input_sequence_sha256': sequence_hash,
            'sequence_compatibility_category': 'exact_sequence_match', 'sequence_match': 'yes', 'coordinate_comparison_performed': 'yes',
            'n_genes_genbank': len({v(r, 'gene') for r in g_rows}), 'n_genes_mitos2': len({v(r, 'gene') for r in m_rows}),
            'n_exact_genes': kinds.count('exact_match'), 'n_minor_difference_genes': kinds.count('minor_boundary_difference'),
            'n_moderate_difference_genes': kinds.count('moderate_boundary_difference'),
            'n_major_difference_genes': kinds.count('major_coordinate_difference') + kinds.count('coordinate_match_codon_difference'),
            'n_missing_genes': sum(k.startswith('missing') for k in kinds), 'n_strand_mismatches': kinds.count('strand_mismatch'),
            'all_13_genes_present_genbank': yes(all(r['genbank_present'] == 'yes' for r in rows)),
            'all_13_genes_present_mitos2': yes(all(r['mitos2_present'] == 'yes' for r in rows)),
            'all_13_exact_match': yes(category == 'all_13_exact'), 'reference_comparison_category': category, 'note': ''})

    # Missing hashes are diagnostics too, but never candidates for pairing.
    if '' in gb_by_hash:
        diagnostic('', gb_by_hash[''], {}, 'no_matching_mitos2_sequence_hash')
    if '' in mi_by_hash:
        diagnostic('', {}, mi_by_hash[''], 'no_matching_genbank_sequence_hash')
    atomic(output, GENE_FIELDS, gene_rows)
    atomic(summary_output, SUMMARY_FIELDS, summaries)
    atomic(mismatch_output, DIAG_FIELDS, diagnostics)
    print(f'unique GenBank sequence hashes: {len(gb_hashes)}')
    print(f'unique MITOS2 sequence hashes: {len(mi_hashes)}')
    print(f'shared sequence hashes: {len(shared)}')
    print(f'single-target shared groups: {single_target}')
    print(f'multiple-target shared groups: {multiple_target}')
    print(f'inconsistent same-sequence annotation groups: {inconsistent}')
    print(f'unmatched GenBank hashes: {len(gb_hashes-mi_hashes)}')
    print(f'unmatched MITOS2 hashes: {len(mi_hashes-gb_hashes)}')
    print(f'reference comparison rows written: {len(summaries)}')
    print(f'gene comparison rows written: {len(gene_rows)}')
    if fail_no_shared and not gene_rows:
        raise RuntimeError('No consistent shared sequence hashes after diagnostics were written.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    for name in ('genbank-table', 'mitos2-table', 'sample-reference-map', 'output', 'reference-summary-output', 'mismatch-output'):
        parser.add_argument('--' + name)
    parser.add_argument('--strict-sequence-match', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--allow-rotation-equivalent', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--fail-on-no-shared-references', action=argparse.BooleanOptionalAction, default=None)
    args = parser.parse_args()
    config = yaml(args.config)
    section = config.get('genbank_mitos2_comparison', {})
    paths, settings = section.get('paths', {}), section.get('settings', {})
    try:
        compare(args.genbank_table or paths.get('genbank_reference_codon_table'), args.mitos2_table or paths.get('mitos2_reference_codon_table'), args.output or paths.get('gene_comparison_table'), args.reference_summary_output or paths.get('reference_summary_table'), args.mismatch_output or paths.get('sequence_mismatch_table'), fail_no_shared=settings.get('fail_on_no_shared_references', True) if args.fail_on_no_shared_references is None else args.fail_on_no_shared_references)
    except Exception as error:
        raise SystemExit(f'GenBank-versus-MITOS2 comparison failed: {error}')


if __name__ == '__main__':
    main()
