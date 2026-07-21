#!/usr/bin/env python3
"""Download primate GenBank records and build sample-level CDS codon annotations."""
import argparse
import csv
import re
import sys
import time
import gzip
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml

try:
    from Bio import Entrez, SeqIO
except ImportError:  # checked in main so importing helpers remains possible in tests
    Entrez = SeqIO = None

OUTPUT_FIELDS = "file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source annotation_fallback_used coordinate_reference_fasta coordinate_reference_accession".split()
SUMMARY_FIELDS = "sample species accession_query accession_source accession_note manifest_file matched_manifest_species species_fasta_path accession_record genbank_file n_cds_features n_coding_position_rows n_genes min_pos max_pos status note".split()
FAIL_FIELDS = "sample species accession_query reason".split()
GENES = {'ND1':'MT-ND1','ND2':'MT-ND2','ND3':'MT-ND3','ND4':'MT-ND4','ND4L':'MT-ND4L','ND5':'MT-ND5','ND6':'MT-ND6','COX1':'MT-CO1','COI':'MT-CO1','COX2':'MT-CO2','COII':'MT-CO2','COX3':'MT-CO3','COIII':'MT-CO3','CYTB':'MT-CYB','ATP6':'MT-ATP6','ATP8':'MT-ATP8'}

def value(row, key):
    return (row.get(key) or '').strip()

def normalize_gene(raw):
    key = re.sub(r'[^A-Z0-9]', '', raw.upper().replace('MT', '', 1))
    return GENES.get(key, raw)

def species_key(species):
    return re.sub(r'_+', '_', re.sub(r'\s+', '_', species.strip().lower())).strip('_')

def read_samples(path, sample_column, species_column):
    """Read normal header TSVs, while retaining legacy two-column sample/species files."""
    with Path(path).open(newline='') as handle:
        raw = list(csv.reader(handle, delimiter='\t'))
    if not raw: return []
    header = raw[0]
    if sample_column in header:
        return [dict(zip(header, row)) for row in raw[1:] if any(row)]
    # Historical config/sample_ref_file.tsv has no header and is sample, species.
    return [{sample_column: row[0], species_column: row[1] if len(row) > 1 else ''}
            for row in raw if row and row[0].strip()]

def accession_for(row, columns):
    return next((value(row, col) for col in columns if value(row, col)), '')

def configured_columns(settings, name, defaults):
    return [item.strip() for item in str(settings.get(name, defaults)).split(',') if item.strip()]

def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return []
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def manifest_rows(paths, settings):
    """Load the resolved manifest first, followed by optional fallback manifests."""
    candidates = [paths.get('reference_summary_file', '')]
    candidates.extend(str(paths.get('reference_summary_fallback_files', '')).split(','))
    loaded = []
    for candidate in candidates:
        candidate = candidate.strip()
        if candidate:
            loaded.append((candidate, read_tsv(candidate)))
    return loaded

def choose_manifest_match(matches, accession_columns):
    """Choose the most review-ready chrM record without excluding usable rows."""
    def score(row):
        review = value(row, 'manual_review_required').lower()
        pairing = value(row, 'reference_pairing_status').lower()
        return (
            review in ('false', 'no', '0'),
            bool(value(row, accession_columns[0])),
            any(word in pairing for word in ('usable', 'final', 'same')),
            bool(value(row, 'final_reference_strategy')),
        )
    best = max(matches, key=score)
    # Record ambiguity only when none of the stated preferences selected a row.
    note = 'multiple_manifest_matches' if len(matches) > 1 and all(score(row) == score(best) for row in matches) else ''
    return best, note

def find_species_fasta(species, fasta_dir, extensions):
    directory = Path(fasta_dir)
    if not directory.is_dir():
        return None
    normalized_extensions = sorted((item.strip() for item in str(extensions).split(',') if item.strip()), key=len, reverse=True)
    target = species_key(species)
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        name = path.name
        for extension in normalized_extensions:
            if name.endswith(extension):
                name = name[:-len(extension)]
                break
        if species_key(name) == target:
            return path
    return None

def fasta_accession(path, regex):
    pattern = re.compile(regex)
    match = pattern.search(path.name)
    if match:
        return match.group(1), 'fasta_filename'
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                match = pattern.search(line)
                return (match.group(1), 'fasta_header') if match else ('', '')
    return '', ''

def resolve_accession(metadata, direct_columns, manifests, settings, paths):
    """Resolve direct sample metadata, then chrM manifests, then a species FASTA."""
    direct = accession_for(metadata, direct_columns)
    if direct:
        return direct, 'sample_ref_file', '', '', '', ''
    manifest_species_columns = configured_columns(settings, 'reference_summary_species_columns', 'target_species,final_chrM_species,preprint_REFERENCE_SPECIES,final_wg_ref_species')
    manifest_accession_columns = configured_columns(settings, 'reference_summary_accession_columns', 'final_chrM_genbank_accn,final_chrM_refseq_accn,final_chrM_accession,chrM_source_accession')
    target = species_key(metadata['species'])
    for index, (manifest_file, rows) in enumerate(manifests):
        matches = [row for row in rows if any(species_key(value(row, column)) == target for column in manifest_species_columns)]
        if matches:
            # An incomplete preferred row must not obscure a less-preferred row
            # that actually supplies a chrM accession.
            usable_matches = [row for row in matches if accession_for(row, manifest_accession_columns)]
            if usable_matches:
                selected, note = choose_manifest_match(usable_matches, manifest_accession_columns)
                accession = accession_for(selected, manifest_accession_columns)
                matched = next((value(selected, column) for column in manifest_species_columns if species_key(value(selected, column)) == target), '')
                return accession, 'reference_manifest' if index == 0 else 'reference_manifest_fallback', note, manifest_file, matched, ''
    if settings.get('infer_accession_from_fasta', False):
        fasta_path = find_species_fasta(metadata['species'], paths.get('species_fasta_dir', ''), paths.get('species_fasta_extensions', '.fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz'))
        if fasta_path:
            accession, source = fasta_accession(fasta_path, settings.get('accession_regex', r'([A-Z]{1,3}_[0-9]+(?:\\.[0-9]+)?)'))
            if accession:
                return accession, source, '', '', '', str(fasta_path)
            return '', 'unresolved', '', '', '', str(fasta_path)
    return '', 'unresolved', '', '', '', ''

def manifest_coordinate_reference(metadata, manifest_file, manifests, settings):
    """Return the materialized coordinate FASTA/accession for a manifest-resolved sample."""
    if not manifest_file:
        return '', '', ''
    target = species_key(metadata['species'])
    columns = configured_columns(settings, 'reference_summary_species_columns', 'target_species,final_chrM_species,preprint_REFERENCE_SPECIES,final_wg_ref_species')
    for path, rows in manifests:
        if path != manifest_file:
            continue
        matches = [r for r in rows if any(species_key(value(r, c)) == target for c in columns)]
        if matches:
            selected, _ = choose_manifest_match(matches, configured_columns(settings, 'reference_summary_accession_columns', 'final_chrM_genbank_accn,final_chrM_refseq_accn,final_chrM_accession,chrM_source_accession'))
            return (value(selected, 'chrM_expected_output_fasta'), value(selected, 'final_chrM_accession'),
                    value(selected, 'final_chrM_species'))
    return '', '', ''

def safe_filename(accession):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', accession) + '.gb'

def write_tsv(path, fields, records):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        out = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        out.writeheader(); out.writerows(records)

def coding_coordinates(feature):
    """Return genomic 0-based positions in coding orientation, including joined CDSs."""
    parts = list(feature.location.parts)
    strand = feature.location.strand or 1
    if strand == -1: parts.reverse()
    positions = []
    for part in parts:
        start, end = int(part.start), int(part.end)
        positions.extend(range(end - 1, start - 1, -1) if strand == -1 else range(start, end))
    return positions

def qualifier(feature, name):
    values = feature.qualifiers.get(name, [])
    return str(values[0]) if values else ''

def parse_record(record, metadata, filename):
    rows, n_cds = [], 0
    for feature in record.features:
        if feature.type != 'CDS': continue
        n_cds += 1
        coords = coding_coordinates(feature)
        sequence = str(feature.extract(record.seq)).upper()
        # codon_start is a one-based offset into the annotated CDS.
        try: offset = max(0, int(qualifier(feature, 'codon_start') or '1') - 1)
        except ValueError: offset = 0
        coords, sequence = coords[offset:], sequence[offset:]
        usable = len(sequence) - (len(sequence) % 3)
        tail = len(sequence) - usable
        raw = qualifier(feature, 'gene') or qualifier(feature, 'label') or qualifier(feature, 'product')
        codon_start = qualifier(feature, 'codon_start') or '1'
        for index in range(0, usable, 3):
            codon, triplet = sequence[index:index + 3], coords[index:index + 3]
            if len(triplet) != 3: continue
            for phase, coordinate in enumerate(triplet, 1):
                rows.append({
                    'file_name': Path(filename).name, 'seq_name': record.id,
                    'sample': value(metadata, 'sample'), 'species': value(metadata, 'species'),
                    'species_key': species_key(value(metadata, 'species')), 'accession': value(metadata, 'accession_query'),
                    'accession_version': record.id, 'reference_id': value(metadata, 'reference_id'),
                    'family': value(metadata, 'family'), 'pos': coordinate + 1,
                    'ref_base_genome': str(record.seq[coordinate]).upper(), 'gene': normalize_gene(raw),
                    'gene_raw': raw, 'product': qualifier(feature, 'product'), 'protein_id': qualifier(feature, 'protein_id'),
                    'strand': '+' if (feature.location.strand or 1) == 1 else '-', 'codon_index': index // 3 + 1,
                    'codon_pos_in_triplet': phase, 'codon_seq': codon,
                    'codon_pos1_genomic': triplet[0] + 1, 'codon_pos2_genomic': triplet[1] + 1,
                    'codon_pos3_genomic': triplet[2] + 1, 'codon_start_qualifier': codon_start,
                    'transl_table': qualifier(feature, 'transl_table'), 'cds_tail_incomplete_bases': tail,
                    'annotation_source': 'GenBank', 'annotation_fallback_used': 'no',
                    'coordinate_reference_fasta': value(metadata, 'coordinate_reference_fasta'),
                    'coordinate_reference_accession': value(metadata, 'coordinate_reference_accession') or value(metadata, 'accession_query'),
                })
    return rows, n_cds

def download(accession, destination, settings):
    if not settings.get('email'): print('WARNING: NCBI Entrez email is unset; set build_primate_codon_table.settings.email to identify requests.', file=sys.stderr)
    else: Entrez.email = str(settings['email'])
    with Entrez.efetch(db='nuccore', id=accession, rettype=settings.get('rettype', 'gb'), retmode=settings.get('retmode', 'text')) as handle:
        text = handle.read()
    if not text.strip(): raise RuntimeError('NCBI returned an empty GenBank record')
    destination.parent.mkdir(parents=True, exist_ok=True); destination.write_text(text)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True); parser.add_argument('--sample'); parser.add_argument('--force-download', action='store_true'); parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    section = yaml(args.config).get('build_primate_codon_table')
    if not section: raise SystemExit('Missing build_primate_codon_table section in config.')
    paths, settings = section['paths'], section.get('settings', {})
    sample_col, species_col = settings.get('sample_column', 'sample'), settings.get('species_column', 'species')
    sample_rows = read_samples(paths['sample_ref_file'], sample_col, species_col)
    if args.sample: sample_rows = [r for r in sample_rows if value(r, sample_col) == args.sample]
    if not sample_rows: raise SystemExit('No samples selected.')
    columns = configured_columns(settings, 'accession_columns', 'accession,accession_version,reference_id,seq_name')
    manifests = manifest_rows(paths, settings)
    if SeqIO is None and not args.dry_run:
        raise SystemExit('Biopython is required for GenBank download/parsing.')
    failures, output, summary, cache = [], [], [], {}
    for raw_metadata in sample_rows:
        metadata = dict(raw_metadata); metadata['sample'] = value(metadata, sample_col); metadata['species'] = value(metadata, species_col)
        accession, source, accession_note, manifest_file, matched_species, fasta_path = resolve_accession(metadata, columns, manifests, settings, paths)
        manifest_fasta, manifest_accession, manifest_species = manifest_coordinate_reference(metadata, manifest_file, manifests, settings)
        # A manifest is coordinate authority: never fall back to original-species FASTA.
        if manifest_fasta:
            fasta_path = manifest_fasta
        if manifest_accession:
            accession = manifest_accession
        metadata['accession_query'] = accession
        coordinate_fasta = fasta_path or (str(find_species_fasta(metadata['species'], paths.get('species_fasta_dir', ''), paths.get('species_fasta_extensions', '.fa,.fasta,.fna'))) or '')
        metadata['coordinate_reference_fasta'] = coordinate_fasta
        metadata['coordinate_reference_accession'] = accession
        base = {'sample':metadata['sample'], 'species':metadata['species'], 'accession_query':accession,
                'accession_source':source, 'accession_note':accession_note,
                'manifest_file':manifest_file, 'matched_manifest_species':matched_species,
                'species_fasta_path':fasta_path, 'accession_record':'', 'genbank_file':'',
                'n_cds_features':0, 'n_coding_position_rows':0, 'n_genes':0, 'min_pos':'',
                'max_pos':'', 'status':'failed', 'note':''}
        if not accession:
            reason = 'No accession found from sample_ref_file, reference manifest, or FASTA.'
            failures.append({**base, 'reason':reason}); base['note'] = reason
            base['status'] = 'dry_run_unresolved' if args.dry_run else 'failed'
            summary.append(base); continue
        gb = Path(paths['genbank_dir']) / safe_filename(accession); base['genbank_file'] = str(gb)
        try:
            force = args.force_download or bool(settings.get('force_download', False))
            if not args.dry_run and (force or not (gb.exists() and settings.get('skip_existing_genbank', True))):
                download(accession, gb, settings); time.sleep(float(settings.get('sleep_seconds', 0.34)))
            if args.dry_run and not gb.exists():
                base.update(status='dry_run_resolved', note='Would download GenBank record.'); summary.append(base); continue
            if args.dry_run:
                base.update(status='dry_run_resolved', note='Would reuse cached GenBank record.'); summary.append(base); continue
            if not gb.exists(): raise FileNotFoundError(f'GenBank file is missing: {gb}')
            if accession not in cache: cache[accession] = SeqIO.read(str(gb), 'genbank')
            parsed, n_cds = parse_record(cache[accession], metadata, gb); output.extend(parsed)
            notes = []
            if len(parsed) < 5000 or len(parsed) > 13000: notes.append(f'Coding rows outside expected mammalian range: {len(parsed)}')
            base.update(accession_record=cache[accession].id, n_cds_features=n_cds, n_coding_position_rows=len(parsed), n_genes=len({r['gene'] for r in parsed}), min_pos=min((r['pos'] for r in parsed), default=''), max_pos=max((r['pos'] for r in parsed), default=''), status='completed' if parsed else 'failed', note='; '.join(notes or ([] if parsed else ['No CDS codon rows parsed.'])))
            if not parsed: failures.append({'sample':metadata['sample'], 'species':metadata['species'], 'accession_query':accession, 'reason':'No CDS codon rows parsed.'})
        except Exception as exc:
            reason = f'{type(exc).__name__}: {exc}'; failures.append({'sample':metadata['sample'], 'species':metadata['species'], 'accession_query':accession, 'reason':reason}); base['note'] = reason
        summary.append(base)
    # Select exactly one source per sample: valid GenBank first, then MITOS2 fallback.
    mitos_paths = yaml(args.config).get('mitos2_annotation', {}).get('paths', {})
    mitos_path = mitos_paths.get('mitos2_cds_table', '')
    mitos_reference_path = mitos_paths.get('mitos2_reference_cds_table', '')
    mitos_rows = read_tsv(mitos_path) if mitos_path else []
    mitos_reference_rows = read_tsv(mitos_reference_path) if mitos_reference_path else []
    genbank_rows = list(output)
    selected = []
    for metadata in sample_rows:
        sample = value(metadata, sample_col)
        gb_rows = [row for row in output if row['sample'] == sample]
        if gb_rows:
            selected.extend(gb_rows)
            continue
        fallback = [dict(row) for row in mitos_rows if value(row, 'sample') == sample]
        if not fallback:
            summary_row = next((r for r in summary if r['sample'] == sample), {})
            accession = value(summary_row, 'accession_query')
            coordinate_fasta = value(summary_row, 'species_fasta_path')
            fallback = [dict(row) for row in mitos_reference_rows if value(row, 'coordinate_reference_accession') == accession or (coordinate_fasta and value(row, 'coordinate_reference_fasta') == coordinate_fasta)]
            for row in fallback:
                row.update(sample=sample, species=value(metadata, species_col), species_key=species_key(value(metadata, species_col)))
        if fallback and settings.get('use_mitos2_if_genbank_fails', True):
            for row in fallback:
                row['annotation_source'] = 'MITOS2'; row['annotation_fallback_used'] = 'yes'
            selected.extend(fallback)
            for row in summary:
                if row['sample'] == sample:
                    row['status'] = 'completed_mitos2_fallback'; row['n_coding_position_rows'] = len(fallback); row['n_genes'] = len({x['gene'] for x in fallback}); row['note'] = '; '.join(filter(None, [row['note'], 'GenBank CDS unavailable; MITOS2 fallback used.']))
    output = selected
    # Compare raw GenBank and MITOS2 CDS boundaries wherever both sources exist.
    comparison_path = yaml(args.config).get('mitos2_annotation', {}).get('paths', {}).get('genbank_mitos2_comparison_table', 'results/qc/codon_table_build/genbank_vs_mitos2_cds_comparison.tsv')
    comparison_fields = 'record_type sample original_species reference_species coordinate_reference_accession gene genbank_start genbank_end genbank_strand genbank_length mitos2_start mitos2_end mitos2_strand mitos2_length start_delta end_delta length_delta strand_match coordinate_match_status n_genbank_cds_genes n_mitos2_cds_genes n_gene_matches n_exact_boundary_matches n_near_boundary_matches n_strand_mismatches missing_in_genbank missing_in_mitos2 comparison_status'.split()
    comparisons = []
    for sample in sorted({value(r, 'sample') for r in genbank_rows + mitos_rows}):
        g = {}; m = {}
        for row in [x for x in genbank_rows if value(x, 'sample') == sample]: g.setdefault(normalize_gene(value(row, 'gene')), []).append(row)
        for row in [x for x in mitos_rows if value(x, 'sample') == sample]: m.setdefault(normalize_gene(value(row, 'gene')), []).append(row)
        exact = near = mismatch = matches = 0
        for gene in sorted(set(g) | set(m)):
            gr, mr = g.get(gene, []), m.get(gene, [])
            if not gr or not mr:
                status = 'missing_mitos2' if gr else 'missing_genbank'
                comparisons.append({'record_type':'gene','sample':sample,'original_species':value((gr or mr)[0], 'species'),'reference_species':'','coordinate_reference_accession':value((gr or mr)[0], 'coordinate_reference_accession'),'gene':gene,'coordinate_match_status':status})
                continue
            gg, mm = gr[0], mr[0]; gs,ge=min(int(x['pos']) for x in gr),max(int(x['pos']) for x in gr); ms,me=min(int(x['pos']) for x in mr),max(int(x['pos']) for x in mr)
            sd,ed=ms-gs,me-ge; sm=value(gg,'strand') == value(mm,'strand'); status='strand_mismatch' if not sm else ('exact' if not sd and not ed else ('near_boundary_delta_le_3' if abs(sd)<=3 and abs(ed)<=3 else 'different'))
            matches += 1; exact += status == 'exact'; near += status == 'near_boundary_delta_le_3'; mismatch += status == 'strand_mismatch'
            comparisons.append({'record_type':'gene','sample':sample,'original_species':value(gg,'species'),'reference_species':'','coordinate_reference_accession':value(gg,'coordinate_reference_accession'),'gene':gene,'genbank_start':gs,'genbank_end':ge,'genbank_strand':value(gg,'strand'),'genbank_length':len(gr),'mitos2_start':ms,'mitos2_end':me,'mitos2_strand':value(mm,'strand'),'mitos2_length':len(mr),'start_delta':sd,'end_delta':ed,'length_delta':len(mr)-len(gr),'strand_match':'yes' if sm else 'no','coordinate_match_status':status})
        comparisons.append({'record_type':'summary','sample':sample,'n_genbank_cds_genes':len(g),'n_mitos2_cds_genes':len(m),'n_gene_matches':matches,'n_exact_boundary_matches':exact,'n_near_boundary_matches':near,'n_strand_mismatches':mismatch,'missing_in_genbank':','.join(sorted(set(m)-set(g))),'missing_in_mitos2':','.join(sorted(set(g)-set(m))),'comparison_status':'compared' if g and m else 'source_unavailable'})
    if settings.get('compare_genbank_and_mitos2', True): write_tsv(comparison_path, comparison_fields, comparisons)
    # Global checks are warnings; they intentionally do not discard useful records.
    required = {'sample','pos','gene','codon_seq','codon_pos_in_triplet'}
    global_notes = []
    if not required.issubset(OUTPUT_FIELDS): global_notes.append('Internal error: output schema lacks required columns.')
    if any(r['codon_pos_in_triplet'] not in (1,2,3) for r in output): global_notes.append('Invalid codon_pos_in_triplet values detected.')
    if any(r['codon_seq'] and len(r['codon_seq']) != 3 for r in output): global_notes.append('Non-triplet codon_seq values detected.')
    if global_notes:
        for row in summary: row['note'] = '; '.join(filter(None, [row['note'], *global_notes]))
    if not args.dry_run: write_tsv(paths['output_table'], OUTPUT_FIELDS, output)
    write_tsv(paths['failed_downloads_table'], FAIL_FIELDS, failures); write_tsv(paths['summary_table'], SUMMARY_FIELDS, summary)
    if args.dry_run:
        resolved = sum(row['status'] == 'dry_run_resolved' for row in summary)
        unresolved = sum(row['status'] == 'dry_run_unresolved' for row in summary)
        print(f'Resolved accessions for {resolved} samples.')
        print(f'Unresolved accessions for {unresolved} samples.')
    else:
        print(f'Built {len(output)} coding-position rows for {sum(r["status"] == "completed" for r in summary)} samples.')

if __name__ == '__main__': main()
