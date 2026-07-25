#!/usr/bin/env python3
"""Annotate coordinate-lifted VCF records with source and human codon matches."""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import (  # noqa: E402
    human_pos, info_format, info_parse, inject_headers, open_text, rows,
    sample_names, source, write_summary, yaml,
)

FIELDS = [
    ('MTCODON_STATUS', 'Codon match status'),
    ('MTCODON_MATCH', 'A gene/phase-matched human codon matches the source reference or alternate codon'),
    ('MTCODON_STRICT_PHASE', 'Strict phase matching enabled'),
    ('MTCODON_GENE_MATCH', 'At least one source-human pair has matching genes'),
    ('MTCODON_PHASE_MATCH', 'At least one gene-matched pair has matching codon phase'),
    ('MTCODON_PRIMATE_GENE', 'Source gene for the selected best pair'),
    ('MTCODON_PRIMATE_CODON', 'Source reference codon for the selected best pair'),
    ('MTCODON_PRIMATE_ALT_CODON', 'Source alternate codon constructed from strand-aware SRC_ALT'),
    ('MTCODON_PRIMATE_PHASE', 'Source codon phase for the selected best pair'),
    ('MTCODON_HUMAN_GENE', 'Human gene for the selected best pair'),
    ('MTCODON_HUMAN_CODON', 'Human codon for the selected best pair'),
    ('MTCODON_HUMAN_PHASE', 'Human codon phase for the selected best pair'),
    ('MTCODON_N_PRIMATE_ANNOTATIONS', 'Number of source CDS annotations at this position'),
    ('MTCODON_N_HUMAN_ANNOTATIONS', 'Number of human CDS annotations at this position'),
    ('MTCODON_N_PAIR_CANDIDATES', 'Number of evaluated source-human annotation pairs'),
    ('MTCODON_OVERLAPPING_CDS', 'Whether either position has multiple CDS annotations'),
    ('MTCODON_PRIMATE_GENES', 'Sorted unique source genes at this position'),
    ('MTCODON_HUMAN_GENES', 'Sorted unique human genes at this position'),
    ('MTCODON_MATCHING_GENES', 'Sorted unique genes with gene, phase, and codon matches'),
    ('MTCODON_AMBIGUOUS_BEST_MATCH', 'Whether distinct annotations tie for the best score'),
]


def complement_base(base):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(str(base).upper(), str(base).upper())


def mutate_codon(codon, phase, alt_base):
    if codon in {'', None, '.', 'NA'} or phase in {'', None, '.', 'NA'} or alt_base in {'', None, '.', 'NA'}:
        return '.'
    try:
        phase = int(phase)
    except (TypeError, ValueError):
        return '.'
    bases, alt = list(str(codon).upper()), str(alt_base).upper()
    if len(bases) != 3 or not 1 <= phase <= 3 or not alt:
        return '.'
    bases[phase - 1] = alt[0]
    return ''.join(bases)


def load_codon_index(path, key_column=None):
    """Stream a plain or gzipped codon table into a one-to-many position index."""
    index = defaultdict(list)
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            try:
                pos = int(row['pos'])
            except (KeyError, TypeError, ValueError):
                continue
            reference_key = (row.get(key_column) or '').strip() if key_column else ''
            index[(reference_key, pos)].append(row)
    return index


# Historical public helper retained, now returning the one-to-many index.
def load(path, key_column=None):
    return load_codon_index(path, key_column), key_column


def find_overlapping_annotations(path, table, key_column=None):
    """Return a non-mutating report of positions annotated to multiple genes."""
    report = []
    for (reference_key, pos), annotations in load_codon_index(path, key_column).items():
        genes = sorted({(row.get('gene') or '').strip() for row in annotations if row.get('gene')})
        if len(genes) > 1:
            report.append({'table': table, 'reference_key': reference_key, 'position': pos,
                           'genes': ','.join(genes), 'number_of_annotations': len(annotations)})
    return sorted(report, key=lambda row: (row['table'], row['reference_key'], row['position']))


def load_sample_reference_map(path):
    mapping = {}
    for row in rows(path):
        sample, reference_key = (row.get('sample') or '').strip(), (row.get('reference_key') or '').strip()
        if sample and reference_key:
            mapping[sample] = reference_key
    return mapping


def _value(row, name, default='.'):
    return row.get(name, default) if row else default


def _row_tie_break(row):
    return tuple(str(row.get(name, '')) for name in
                 ('gene', 'codon_seq', 'codon_pos_in_triplet', 'strand'))


def evaluate_candidates(source_rows, human_rows, source_alt):
    candidates = []
    for source_row in source_rows:
        strand_alt = complement_base(source_alt) if source_row.get('strand', '+') == '-' else source_alt
        source_gene = source_row.get('gene', '')
        source_phase = str(source_row.get('codon_pos_in_triplet', ''))
        source_codon = source_row.get('codon_seq', '.')
        alternate_codon = mutate_codon(source_codon, source_phase, strand_alt)
        for human_row in human_rows:
            human_gene = human_row.get('gene', '')
            human_phase = str(human_row.get('codon_pos_in_triplet', ''))
            human_codon = human_row.get('codon_seq', '')
            gene_match = source_gene == human_gene
            phase_match = source_phase == human_phase
            codon_match = human_codon in {source_codon, alternate_codon}
            candidates.append({
                'source': source_row, 'human': human_row, 'source_gene': source_gene,
                'human_gene': human_gene, 'source_phase': source_phase, 'human_phase': human_phase,
                'source_codon': source_codon, 'human_codon': human_codon,
                'alternate_codon': alternate_codon, 'gene_match': gene_match,
                'phase_match': phase_match, 'codon_match': codon_match,
            })
    return candidates


def _primary_score(candidate):
    # Codon equality is relevant only after gene and phase equality.
    return (int(candidate['gene_match']),
            int(candidate['gene_match'] and candidate['phase_match']),
            int(candidate['gene_match'] and candidate['phase_match'] and candidate['codon_match']))


def _tie_break(candidate):
    return tuple(str(candidate[name]) for name in (
        'source_gene', 'human_gene', 'source_codon', 'human_codon', 'source_phase', 'human_phase',
        'alternate_codon'))


def annotate(source_rows, human_rows, source_alt, strict):
    """Build INFO values for all annotation pairs and a deterministic representative."""
    candidates = evaluate_candidates(source_rows, human_rows, source_alt)
    source_genes = sorted({row.get('gene', '') for row in source_rows if row.get('gene')})
    human_genes = sorted({row.get('gene', '') for row in human_rows if row.get('gene')})
    overlap = len(source_rows) > 1 or len(human_rows) > 1
    vals = {
        'MTCODON_STRICT_PHASE': 'yes' if strict else 'no',
        'MTCODON_N_PRIMATE_ANNOTATIONS': str(len(source_rows)),
        'MTCODON_N_HUMAN_ANNOTATIONS': str(len(human_rows)),
        'MTCODON_N_PAIR_CANDIDATES': str(len(candidates)),
        'MTCODON_OVERLAPPING_CDS': 'yes' if overlap else 'no',
        'MTCODON_PRIMATE_GENES': ','.join(source_genes) or '.',
        'MTCODON_HUMAN_GENES': ','.join(human_genes) or '.',
    }
    any_gene = any(c['gene_match'] for c in candidates)
    any_phase = any(c['gene_match'] and c['phase_match'] for c in candidates)
    valid = [c for c in candidates if c['gene_match'] and c['phase_match'] and c['codon_match']]
    vals.update(MTCODON_GENE_MATCH='yes' if any_gene else 'no',
                MTCODON_PHASE_MATCH='yes' if any_phase else 'no',
                MTCODON_MATCH='yes' if valid else 'no',
                MTCODON_MATCHING_GENES=','.join(sorted({c['source_gene'] for c in valid})) or '.')
    best = None
    ambiguous = False
    if candidates:
        best_score = max(map(_primary_score, candidates))
        tied = [c for c in candidates if _primary_score(c) == best_score]
        ambiguous = len({_tie_break(c) for c in tied}) > 1
        best = min(tied, key=_tie_break)
    vals['MTCODON_AMBIGUOUS_BEST_MATCH'] = 'yes' if ambiguous else 'no'
    # Preserve useful legacy context when only one side is annotated, while making
    # its representative independent of table order too.
    selected_source = best['source'] if best else min(source_rows, key=_row_tie_break, default=None)
    selected_human = best['human'] if best else min(human_rows, key=_row_tie_break, default=None)
    selected_alt = '.'
    if selected_source:
        strand_alt = complement_base(source_alt) if selected_source.get('strand', '+') == '-' else source_alt
        selected_alt = mutate_codon(selected_source.get('codon_seq'),
                                    selected_source.get('codon_pos_in_triplet'), strand_alt)
    vals.update(
        MTCODON_PRIMATE_GENE=_value(selected_source, 'gene'),
        MTCODON_PRIMATE_CODON=_value(selected_source, 'codon_seq'),
        MTCODON_PRIMATE_ALT_CODON=selected_alt,
        MTCODON_PRIMATE_PHASE=_value(selected_source, 'codon_pos_in_triplet'),
        MTCODON_HUMAN_GENE=_value(selected_human, 'gene'),
        MTCODON_HUMAN_CODON=_value(selected_human, 'codon_seq'),
        MTCODON_HUMAN_PHASE=_value(selected_human, 'codon_pos_in_triplet'),
    )
    return vals, candidates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True); parser.add_argument('--sample')
    parser.add_argument('--input'); parser.add_argument('--output')
    args = parser.parse_args(); config = yaml(args.config); section = config['codon_match']
    paths, settings = section['paths'], section['settings']; strict = bool(settings.get('strict_phase_match', True))
    reference_table, map_table = paths.get('reference_codon_table'), paths.get('sample_reference_map')
    if reference_table and map_table:
        primate = load_codon_index(reference_table, 'reference_key')
        sample_references = load_sample_reference_map(map_table)
    else:  # compatibility with historical sample-level tables
        primate = load_codon_index(paths['all_primate_position_codon_table'], 'sample')
        sample_references = {}
    human = load_codon_index(paths['human_codon_table'])
    samples = [args.sample] if args.sample else sample_names(config)
    if args.input:
        samples = [args.sample or Path(args.input).name.split('.')[0]]
    if not samples:
        raise SystemExit('No samples found; supply --sample or --input.')
    allrows = []
    for sample in samples:
        inp = Path(args.input) if args.input else Path(paths['input_vcf_dir']) / str(settings['input_vcf_pattern']).format(sample=sample)
        out = Path(args.output) if args.output else Path(paths['output_dir']) / 'vcf_codon' / f"{sample}{settings['output_suffix']}"
        if not inp.exists():
            raise SystemExit(f'Missing input VCF for {sample}: {inp}')
        out.parent.mkdir(parents=True, exist_ok=True); header = []; body = []; counts = Counter(); metrics = Counter()
        with open_text(inp) as handle:
            for line in handle:
                if line.startswith('#'):
                    header.append(line); continue
                fields = line.rstrip('\n').split('\t'); info = info_parse(fields[7])
                _, source_pos, _, source_alt = source(info); lifted_pos = human_pos(fields, info)
                reference_key = sample_references.get(sample, sample)
                source_rows = primate.get((reference_key, source_pos), []) if source_pos else []
                human_rows = human.get(('', lifted_pos), []) if lifted_pos else []
                vals, candidates = annotate(source_rows, human_rows, source_alt, strict)
                if not source_pos or not lifted_pos: status = 'MISSING_COORD'
                elif not source_rows: status = 'SKIPPED_NONCODING'
                elif not human_rows: status = 'NO_HUMAN_CODON'
                elif vals['MTCODON_MATCH'] == 'yes': status = 'PASS'
                elif strict and vals['MTCODON_GENE_MATCH'] == 'no': status = 'GENE_MISMATCH'
                elif strict and vals['MTCODON_PHASE_MATCH'] == 'no': status = 'PHASE_MISMATCH'
                else: status = 'MISMATCH'
                vals['MTCODON_STATUS'] = status
                info.update(vals); fields[7] = info_format(info); body.append('\t'.join(fields) + '\n'); counts[status] += 1
                metrics['records_with_overlapping_source_cds'] += len(source_rows) > 1
                metrics['records_with_overlapping_human_cds'] += len(human_rows) > 1
                metrics['records_with_overlapping_cds'] += vals['MTCODON_OVERLAPPING_CDS'] == 'yes'
                metrics['records_with_ambiguous_best_match'] += vals['MTCODON_AMBIGUOUS_BEST_MATCH'] == 'yes'
                metrics['records_with_multiple_pair_candidates'] += len(candidates) > 1
        with out.open('w') as handle:
            handle.writelines(inject_headers(header, FIELDS, 'MTCODON')); handle.writelines(body)
        row = {'sample': sample, 'input_vcf': str(inp), 'output_vcf': str(out), 'total_records': len(body),
               **{f'status_{name}': counts[name] for name in ['PASS', 'SKIPPED_NONCODING', 'NO_HUMAN_CODON', 'GENE_MISMATCH', 'PHASE_MISMATCH', 'MISMATCH', 'MISSING_COORD']},
               **{name: metrics[name] for name in ['records_with_overlapping_source_cds', 'records_with_overlapping_human_cds', 'records_with_overlapping_cds', 'records_with_ambiguous_best_match', 'records_with_multiple_pair_candidates']},
               'strict_phase_match': strict, 'status': 'completed'}
        write_summary(Path(paths['reports_dir']) / f'{sample}.codon_match_summary.tsv', row); allrows.append(row)
    if allrows:
        path = Path(paths['reports_dir']) / 'all_samples.codon_match_summary.tsv'; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(allrows[0]), delimiter='\t'); writer.writeheader(); writer.writerows(allrows)


if __name__ == '__main__':
    main()
