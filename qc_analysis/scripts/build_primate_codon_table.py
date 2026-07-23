#!/usr/bin/env python3
"""Download primate GenBank records and build sample-level CDS codon annotations."""
import argparse
import csv
import re
import sys
import time
import gzip
import os
import tempfile
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml

repo_root = Path(__file__).resolve().parents[2]

try:
    from Bio import Entrez, SeqIO
except ImportError:  # checked in main so importing helpers remains possible in tests
    Entrez = SeqIO = None

OUTPUT_FIELDS = "file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source annotation_fallback_used coordinate_reference_fasta coordinate_reference_accession".split()
SUMMARY_FIELDS = "sample species accession_query accession_source accession_note manifest_file matched_manifest_species species_fasta_path accession_record genbank_file n_cds_features n_coding_position_rows n_genes min_pos max_pos status note".split()
FAIL_FIELDS = "sample species accession_query reason".split()
MITOS2_FALLBACK_SELECTION_FIELDS = "sample species accession coordinate_reference_fasta coordinate_reference_accession original_sample_fasta canonical_sample_fasta original_mitos2_fasta canonical_mitos2_fasta fasta_match group_row_count normalized_gene_count has_all_13_protein_coding_genes coding_row_count_in_expected_range fallback_match_mode n_candidate_rows n_selected_rows_before_dedup n_selected_rows_after_dedup n_duplicate_rows_collapsed n_candidate_reference_groups selected_reference_group selection_status rejection_reason note".split()
MITOS2_NUMERIC_FIELDS = ('pos', 'codon_index', 'codon_pos_in_triplet', 'codon_pos1_genomic', 'codon_pos2_genomic', 'codon_pos3_genomic')
EXPECTED_MAMMALIAN_CODING_TOTAL = 11400
SUCCESS_STATUSES = {
    'completed',
    'completed_mitos2_fallback',
}
GENES = {'ND1':'MT-ND1','ND2':'MT-ND2','ND3':'MT-ND3','ND4':'MT-ND4','ND4L':'MT-ND4L','ND5':'MT-ND5','ND6':'MT-ND6','COX1':'MT-CO1','COI':'MT-CO1','COX2':'MT-CO2','COII':'MT-CO2','COX3':'MT-CO3','COIII':'MT-CO3','CYTB':'MT-CYB','ATP6':'MT-ATP6','ATP8':'MT-ATP8'}
PROTEIN_CODING_GENES = frozenset(GENES.values())

def value(row, key):
    return (row.get(key) or '').strip()

def canonical_path(raw, repo_root):
    raw = (raw or "").strip()
    if not raw:
        return ""

    path = Path(raw)
    if not path.is_absolute():
        path = repo_root / path

    try:
        return str(path.resolve())
    except OSError:
        return str(path.absolute())

def summarize_build_status(summary):
    """Count final sample outcomes, including both successful annotation routes."""
    counts = {
        'total': len(summary), 'completed': 0, 'completed_genbank': 0,
        'completed_mitos2_fallback': 0, 'failed': 0, 'other': 0,
    }
    for row in summary:
        status = value(row, 'status')
        if status == 'completed':
            counts['completed_genbank'] += 1
            counts['completed'] += 1
        elif status == 'completed_mitos2_fallback':
            counts['completed_mitos2_fallback'] += 1
            counts['completed'] += 1
        elif status == 'failed':
            counts['failed'] += 1
        else:
            counts['other'] += 1
    return counts

def warn_if_output_summary_disagree(summary, output):
    """Warn, without failing, if final output and successful status sets differ."""
    successful_samples = {value(row, 'sample') for row in summary
                          if value(row, 'status') in SUCCESS_STATUSES}
    output_samples = {value(row, 'sample') for row in output if value(row, 'sample')}
    if successful_samples != output_samples:
        print('WARNING: build summary and final codon table sample sets disagree.', file=sys.stderr)
        print('Successful-without-output: ' + ', '.join(sorted(successful_samples - output_samples)), file=sys.stderr)
        print('Output-without-success-status: ' + ', '.join(sorted(output_samples - successful_samples)), file=sys.stderr)
    return successful_samples == output_samples

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

def build_fasta_index(fasta_dir, extensions):
    """Index the first deterministic FASTA for each normalized species name."""
    directory = Path(fasta_dir)
    if not directory.is_dir():
        return {}
    # Preserve legacy longest-extension matching; lexical filenames break ties.
    normalized_extensions = sorted((item.strip() for item in str(extensions).split(',') if item.strip()), key=len, reverse=True)
    index = {}
    # Match the legacy directory scan: lexical file order first, longest matching
    # extension second, then retain the first matching species.
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        for extension in normalized_extensions:
            if path.name.endswith(extension):
                index.setdefault(species_key(path.name[:-len(extension)]), path)
                break
    return index

def find_species_fasta(species, fasta_dir, extensions, fasta_index=None):
    if fasta_index is not None:
        return fasta_index.get(species_key(species))
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

def resolve_accession(metadata, direct_columns, manifests, settings, paths, fasta_index=None):
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
        fasta_path = find_species_fasta(metadata['species'], paths.get('species_fasta_dir', ''), paths.get('species_fasta_extensions', '.fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz'), fasta_index)
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

def valid_phase(value):
    try:
        return int(value) in (1, 2, 3)
    except Exception:
        return False

def normalize_mitos2_numeric_fields(row):
    """Make MITOS2 TSV values comparable to GenBank-derived numeric values."""
    for field in MITOS2_NUMERIC_FIELDS:
        raw = value(row, field)
        if raw:
            try:
                row[field] = str(int(raw))
            except (TypeError, ValueError):
                row[field] = raw
    return row

def fallback_duplicate_key(row, include_sample=False):
    fields = ('coordinate_reference_fasta', 'coordinate_reference_accession', 'gene',
              'pos', 'codon_index', 'codon_pos_in_triplet')
    if include_sample:
        fields = ('sample', 'pos', 'gene', 'codon_index', 'codon_pos_in_triplet')
    return tuple(value(row, field) for field in fields)

def deduplicate_rows(rows, key):
    """Keep the first deterministic occurrence of each annotation position."""
    unique, seen = [], set()
    for row in rows:
        row_key = key(row)
        if row_key not in seen:
            seen.add(row_key)
            unique.append(row)
    return unique, len(rows) - len(unique)

def select_reference_fallback(rows, coordinate_fasta, accession, mitos_reference_rows_by_canonical_fasta=None):
    """Choose one complete reference group before assigning a fallback sample.

    Row-level deduplication is deliberately not part of selection: a duplicated
    annotation set must not make several reference groups appear to be one.
    """
    exact_fasta_row_ids = ({id(row) for row in
                            mitos_reference_rows_by_canonical_fasta.get(coordinate_fasta, [])}
                           if coordinate_fasta and mitos_reference_rows_by_canonical_fasta is not None else None)
    groups = defaultdict(list)
    for row in rows:
        group = (value(row, 'coordinate_reference_fasta'), value(row, 'coordinate_reference_accession'))
        groups[group].append(row)
    profiles = []
    for group, group_rows in groups.items():
        canonical_group_fasta = (value(group_rows[0], '_canonical_coordinate_reference_fasta') or
                                 canonical_path(group[0], repo_root))
        genes = {normalize_gene(value(row, 'gene')) for row in group_rows}
        n_rows = len(group_rows)
        profile = {
            'group': group, 'rows': group_rows, 'n_rows': n_rows,
            'n_genes': len(genes), 'has_13_genes': PROTEIN_CODING_GENES.issubset(genes),
            'in_expected_range': 5000 <= n_rows <= 13000,
            'canonical_fasta': canonical_group_fasta,
            'fasta_match': bool(coordinate_fasta and canonical_group_fasta == coordinate_fasta and
                                (exact_fasta_row_ids is None or id(group_rows[0]) in exact_fasta_row_ids)),
            'accession_match': bool(accession and group[1] == accession),
            'distance': abs(n_rows - EXPECTED_MAMMALIAN_CODING_TOTAL),
        }
        # The lexical group key is intentionally last: it makes selection stable
        # while preserving whether the biological criteria were tied.
        profile['rank'] = (-profile['fasta_match'], -profile['accession_match'],
                           -profile['has_13_genes'], -profile['in_expected_range'],
                           profile['distance'])
        profiles.append(profile)
    if not profiles:
        return [], 'none', [], False
    profiles.sort(key=lambda p: (*p['rank'], p['group'][0], p['group'][1]))
    selected = profiles[0]
    ambiguous = sum(profile['rank'] == selected['rank'] for profile in profiles) > 1
    for profile in profiles:
        if profile is selected:
            profile['rejection_reason'] = ''
            continue
        reason = []
        for label, key in (('no_exact_coordinate_reference_fasta', 'fasta_match'),
                           ('no_exact_coordinate_reference_accession', 'accession_match'),
                           ('does_not_contain_all_13_normalized_protein_coding_genes', 'has_13_genes'),
                           ('coding_row_count_outside_5000_13000', 'in_expected_range')):
            if selected[key] and not profile[key]:
                reason.append(label)
                break
        if not reason and profile['distance'] > selected['distance']:
            reason.append('coding_row_count_farther_from_expected_mammalian_total')
        if not reason:
            reason.append('deterministic_lexical_tiebreaker')
        profile['rejection_reason'] = ';'.join(reason)
    match_mode = ('coordinate_fasta' if selected['fasta_match'] else
                  'accession' if selected['accession_match'] else 'gene_count_row_count')
    return [dict(row) for row in selected['rows']], match_mode, profiles, ambiguous

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

class EntrezRateLimiter:
    """Thread-safe limiter that spaces Entrez request starts globally."""
    def __init__(self, requests_per_second):
        if requests_per_second <= 0:
            raise ValueError('requests per second must be a positive number')
        self.interval = 1.0 / requests_per_second
        self.lock = threading.Lock()
        self.next_request = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            scheduled = max(now, self.next_request)
            self.next_request = scheduled + self.interval
        time.sleep(max(0, scheduled - now))


def positive_integer(raw, option):
    try:
        number = int(raw)
    except (TypeError, ValueError):
        raise SystemExit(f'{option} must be a positive integer.')
    if number < 1:
        raise SystemExit(f'{option} must be a positive integer.')
    return number


def resolve_workers(cli_value, settings, setting_name, fallback=1):
    if cli_value is not None:
        return positive_integer(cli_value, '--' + setting_name.replace('_', '-'))
    configured = settings.get(setting_name)
    if configured not in (None, ''):
        return positive_integer(configured, f'build_primate_codon_table.settings.{setting_name}')
    if setting_name == 'workers' and os.environ.get('SLURM_CPUS_PER_TASK'):
        return positive_integer(os.environ['SLURM_CPUS_PER_TASK'], 'SLURM_CPUS_PER_TASK')
    return fallback


def download(accession, destination, settings, limiter=None):
    if not settings.get('email'): print('WARNING: NCBI Entrez email is unset; set build_primate_codon_table.settings.email to identify requests.', file=sys.stderr)
    else: Entrez.email = str(settings['email'])
    if settings.get('entrez_api_key'):
        Entrez.api_key = str(settings['entrez_api_key'])
    if limiter:
        limiter.wait()
    with Entrez.efetch(db='nuccore', id=accession, rettype=settings.get('rettype', 'gb'), retmode=settings.get('retmode', 'text')) as handle:
        text = handle.read()
    if not text.strip(): raise RuntimeError('NCBI returned an empty GenBank record')
    destination.parent.mkdir(parents=True, exist_ok=True)
    # A same-directory temporary file and replace prevent readers seeing partial cache files.
    fd, temporary = tempfile.mkstemp(prefix=destination.name + '.', suffix='.tmp', dir=destination.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            handle.write(text)
        if not Path(temporary).stat().st_size:
            raise RuntimeError('NCBI returned an empty GenBank record')
        Path(temporary).replace(destination)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


@dataclass
class SampleResult:
    sample: str
    order: int
    metadata: dict
    rows: list
    summary: dict
    failure: dict | None
    accession: str
    genbank_file: str
    status: str
    error: str = ''


def prepare_sample(order_metadata, context):
    order, raw_metadata = order_metadata
    settings, paths, columns, manifests, fasta_index, dry_run = context
    sample_col, species_col = settings.get('sample_column', 'sample'), settings.get('species_column', 'species')
    metadata = dict(raw_metadata); metadata['sample'] = value(metadata, sample_col); metadata['species'] = value(metadata, species_col)
    accession, source, accession_note, manifest_file, matched_species, fasta_path = resolve_accession(metadata, columns, manifests, settings, paths, fasta_index)
    manifest_fasta, manifest_accession, _ = manifest_coordinate_reference(metadata, manifest_file, manifests, settings)
    if manifest_fasta: fasta_path = manifest_fasta
    if manifest_accession: accession = manifest_accession
    metadata['accession_query'] = accession
    coordinate = fasta_path or str(find_species_fasta(metadata['species'], paths.get('species_fasta_dir', ''), paths.get('species_fasta_extensions', '.fa,.fasta,.fna'), fasta_index) or '')
    metadata['coordinate_reference_fasta'] = coordinate; metadata['coordinate_reference_accession'] = accession
    metadata['_canonical_coordinate_reference_fasta'] = canonical_path(coordinate, repo_root)
    original_species_fasta = fasta_path or coordinate
    base = {'sample':metadata['sample'], 'species':metadata['species'], 'accession_query':accession, 'accession_source':source, 'accession_note':accession_note, 'manifest_file':manifest_file, 'matched_manifest_species':matched_species, 'species_fasta_path':original_species_fasta, '_canonical_species_fasta_path':canonical_path(original_species_fasta, repo_root), 'accession_record':'', 'genbank_file':'', 'n_cds_features':0, 'n_coding_position_rows':0, 'n_genes':0, 'min_pos':'', 'max_pos':'', 'status':'failed', 'note':''}
    if not accession:
        reason = 'No accession found from sample_ref_file, reference manifest, or FASTA.'
        base['note'] = reason; base['status'] = 'dry_run_unresolved' if dry_run else 'failed'
        return SampleResult(metadata['sample'], order, metadata, [], base, {**base, 'reason':reason}, accession, '', base['status'], reason)
    gb = Path(paths['genbank_dir']) / safe_filename(accession); base['genbank_file'] = str(gb)
    return SampleResult(metadata['sample'], order, metadata, [], base, None, accession, str(gb), 'prepared')


def parse_sample(result):
    if result.status != 'prepared': return result
    try:
        record = SeqIO.read(result.genbank_file, 'genbank')
        rows, n_cds = parse_record(record, result.metadata, result.genbank_file)
        notes = []
        if len(rows) < 5000 or len(rows) > 13000: notes.append(f'Coding rows outside expected mammalian range: {len(rows)}')
        result.rows = rows
        result.summary.update(accession_record=record.id, n_cds_features=n_cds, n_coding_position_rows=len(rows), n_genes=len({r['gene'] for r in rows}), min_pos=min((r['pos'] for r in rows), default=''), max_pos=max((r['pos'] for r in rows), default=''), status='completed' if rows else 'failed', note='; '.join(notes or ([] if rows else ['No CDS codon rows parsed.'])))
        result.status = result.summary['status']
        if not rows: result.failure = {'sample':result.sample, 'species':result.metadata['species'], 'accession_query':result.accession, 'reason':'No CDS codon rows parsed.'}
    except Exception as exc:
        result.error = f'{type(exc).__name__}: {exc}'; result.summary['note'] = result.error; result.failure = {'sample':result.sample, 'species':result.metadata['species'], 'accession_query':result.accession, 'reason':result.error}; result.status = 'failed'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True); parser.add_argument('--sample'); parser.add_argument('--force-download', action='store_true'); parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--workers'); parser.add_argument('--download-workers')
    args = parser.parse_args()
    section = yaml(args.config).get('build_primate_codon_table')
    if not section: raise SystemExit('Missing build_primate_codon_table section in config.')
    paths, settings = section['paths'], section.get('settings', {})
    workers = resolve_workers(args.workers, settings, 'workers')
    api_key = bool(settings.get('entrez_api_key'))
    download_workers = resolve_workers(args.download_workers, settings, 'download_workers', 1) if args.download_workers is not None or settings.get('download_workers') not in (None, '') else 1
    sample_col, species_col = settings.get('sample_column', 'sample'), settings.get('species_column', 'species')
    sample_rows = read_samples(paths['sample_ref_file'], sample_col, species_col)
    if args.sample: sample_rows = [r for r in sample_rows if value(r, sample_col) == args.sample]
    if not sample_rows: raise SystemExit('No samples selected.')
    columns = configured_columns(settings, 'accession_columns', 'accession,accession_version,reference_id,seq_name'); manifests = manifest_rows(paths, settings)
    if SeqIO is None and not args.dry_run: raise SystemExit('Biopython is required for GenBank download/parsing.')
    fasta_index = build_fasta_index(paths.get('species_fasta_dir', ''), paths.get('species_fasta_extensions', '.fa,.fasta,.fna,.fa.gz,.fasta.gz,.fna.gz'))
    context = (settings, paths, columns, manifests, fasta_index, args.dry_run)
    # Resolution is serial (manifest/FASTA coordination); parsing is parallel and immutable per sample.
    prepared = [prepare_sample(item, context) for item in enumerate(sample_rows)]
    resolved_accessions = sum(result.status == 'prepared' for result in prepared)
    print(f'Resolved {len(prepared)} samples:', file=sys.stderr)
    print(f'  {resolved_accessions} with GenBank accession', file=sys.stderr)
    print(f'  {len(prepared) - resolved_accessions} requiring fallback/unresolved', file=sys.stderr)
    force = args.force_download or bool(settings.get('force_download', False))
    if args.dry_run:
        for result in prepared:
            if result.status == 'prepared':
                result.summary.update(status='dry_run_resolved', note='Would reuse cached GenBank record.' if Path(result.genbank_file).exists() else 'Would download GenBank record.')
                result.status = 'dry_run_resolved'
    else:
        unique = {r.accession: Path(r.genbank_file) for r in prepared if r.status == 'prepared'}
        missing = [(a, p) for a, p in unique.items() if force or not (p.exists() and p.stat().st_size > 0 and settings.get('skip_existing_genbank', True))]
        print('GenBank cache:', file=sys.stderr)
        print(f'  {len(unique)} unique accessions', file=sys.stderr)
        print(f'  {len(missing)} missing downloads', file=sys.stderr)
        rate = float(settings.get('requests_per_second_with_api_key' if api_key else 'requests_per_second_without_api_key', 10 if api_key else 3))
        limiter = EntrezRateLimiter(rate)
        def fetch(item):
            accession, path = item
            # Another invocation may finish safely before this worker begins.
            if not force and path.exists() and path.stat().st_size > 0: return accession, None
            try: download(accession, path, settings, limiter); return accession, None
            except Exception as exc: return accession, f'{type(exc).__name__}: {exc}'
        download_errors = {}
        if missing:
            with ThreadPoolExecutor(max_workers=min(download_workers, len(missing))) as executor:
                download_errors = dict(executor.map(fetch, missing))
        for result in prepared:
            if result.status == 'prepared' and result.accession in download_errors and download_errors[result.accession]:
                result.error = download_errors[result.accession]; result.summary['note'] = result.error; result.failure = {'sample':result.sample, 'species':result.metadata['species'], 'accession_query':result.accession, 'reason':result.error}; result.status = 'failed'
        parse_targets = [result for result in prepared if result.status == 'prepared']
        print(f'  {len(parse_targets)} parse targets', file=sys.stderr)
        if parse_targets:
            with ThreadPoolExecutor(max_workers=min(workers, len(parse_targets))) as executor:
                parsed_by_order = {result.order: result for result in executor.map(parse_sample, parse_targets)}
            prepared = [parsed_by_order.get(result.order, result) for result in prepared]
    prepared.sort(key=lambda r: r.order)
    failures = [r.failure for r in prepared if r.failure]; summary = [r.summary for r in prepared]; output = [row for r in prepared for row in r.rows]
    # Select exactly one source per sample: valid GenBank first, then MITOS2 fallback.
    mitos_paths = yaml(args.config).get('mitos2_annotation', {}).get('paths', {})
    mitos_path = mitos_paths.get('mitos2_cds_table', '')
    mitos_reference_path = mitos_paths.get('mitos2_reference_cds_table', '')
    mitos_rows = read_tsv(mitos_path) if mitos_path else []
    mitos_reference_rows = read_tsv(mitos_reference_path) if mitos_reference_path else []
    for mitos_row in mitos_rows + mitos_reference_rows:
        # Keep the source string for output, while comparing the physical file.
        mitos_row['_canonical_coordinate_reference_fasta'] = canonical_path(
            value(mitos_row, 'coordinate_reference_fasta'), repo_root)
    mitos_reference_rows_by_canonical_fasta = defaultdict(list)
    for mitos_row in mitos_reference_rows:
        mitos_reference_rows_by_canonical_fasta[
            value(mitos_row, '_canonical_coordinate_reference_fasta')].append(mitos_row)
    genbank_rows = list(output)
    selected, fallback_selection_summary = [], []
    genbank_rows_by_sample = defaultdict(list)
    for row in output: genbank_rows_by_sample[row['sample']].append(row)
    mitos_rows_by_sample = defaultdict(list)
    for row in mitos_rows: mitos_rows_by_sample[value(row, 'sample')].append(row)
    summary_by_sample = {row['sample']: row for row in summary}
    for metadata in sample_rows:
        sample = value(metadata, sample_col)
        gb_rows = genbank_rows_by_sample[sample]
        if gb_rows:
            selected.extend(gb_rows)
            continue
        candidate_rows = [dict(row) for row in mitos_rows_by_sample[sample]]
        match_mode = 'sample_level'
        profiles = []
        ambiguous = False
        if candidate_rows:
            summary_row = summary_by_sample.get(sample, {})
            original_sample_fasta = (value(summary_row, 'species_fasta_path') or
                                     value(summary_row, 'coordinate_reference_fasta'))
            fallback, match_mode, profiles, ambiguous = select_reference_fallback(
                candidate_rows,
                value(summary_row, '_canonical_species_fasta_path') or
                value(summary_row, '_canonical_coordinate_reference_fasta'),
                value(summary_row, 'accession_query'))
        else:
            summary_row = summary_by_sample.get(sample, {})
            accession = value(summary_row, 'accession_query')
            original_sample_fasta = (value(summary_row, 'species_fasta_path') or
                                     value(summary_row, 'coordinate_reference_fasta'))
            coordinate_fasta = (value(summary_row, '_canonical_species_fasta_path') or
                                value(summary_row, '_canonical_coordinate_reference_fasta'))
            # Group every reference candidate before choosing; restricting the input
            # to the first matching FASTA/accession would hide competing groups.
            fallback, match_mode, profiles, ambiguous = select_reference_fallback(
                mitos_reference_rows, coordinate_fasta, accession,
                mitos_reference_rows_by_canonical_fasta)
            for row in fallback:
                row.update(sample=sample, species=value(metadata, species_col), species_key=species_key(value(metadata, species_col)))
        candidate_groups = len(profiles)
        selected_group_tuple = (value(fallback[0], 'coordinate_reference_fasta'),
                                value(fallback[0], 'coordinate_reference_accession')) if fallback else ('', '')
        selected_group = '|'.join(selected_group_tuple) if fallback else ''
        n_candidate_rows = len(fallback)
        if fallback and settings.get('use_mitos2_if_genbank_fails', True):
            for row in fallback:
                normalize_mitos2_numeric_fields(row)
                row['annotation_source'] = 'MITOS2'; row['annotation_fallback_used'] = 'yes'
                row['gene'] = normalize_gene(value(row, 'gene'))
            # Reference-level input tables can contain repeated annotation sets;
            # remove those before and after assigning the target sample identity.
            fallback, reference_duplicates = deduplicate_rows(fallback, fallback_duplicate_key)
            fallback, sample_duplicates = deduplicate_rows(fallback, lambda row: fallback_duplicate_key(row, include_sample=True))
            duplicates_collapsed = reference_duplicates + sample_duplicates
            assert len({fallback_duplicate_key(row, include_sample=True) for row in fallback}) == len(fallback), \
                f'MITOS2 fallback duplicate keys remain for {sample}'
            assert len(fallback) <= 13000, f'MITOS2 fallback coding rows exceed 13000 for {sample}: {len(fallback)}'
            selected.extend(fallback)
            notes = ['GenBank CDS unavailable; MITOS2 fallback used.']
            if candidate_groups > 1:
                notes.append(f'MITOS2 fallback selected one reference group from {candidate_groups} candidates: {selected_group}')
            if ambiguous:
                notes.append('MITOS2 fallback reference-group selection ambiguous before deterministic lexical tie-breaker.')
            if duplicates_collapsed:
                notes.append(f'MITOS2 fallback duplicate rows collapsed: {duplicates_collapsed}')
            if len(fallback) < 5000 or len(fallback) > 13000:
                notes.append(f'Coding rows outside expected mammalian range after MITOS2 fallback: {len(fallback)}')
            row = summary_by_sample.get(sample)
            if row:
                row['status'] = 'completed_mitos2_fallback'; row['n_coding_position_rows'] = len(fallback); row['n_genes'] = len({normalize_gene(value(x, 'gene')) for x in fallback}); row['note'] = '; '.join(filter(None, [row['note'], *notes]))
            for profile in profiles:
                group = profile['group']
                is_selected = group == selected_group_tuple
                fallback_selection_summary.append({
                    'sample': sample, 'species': value(metadata, species_col),
                    'accession': value(summary_by_sample.get(sample, {}), 'accession_query'),
                    'coordinate_reference_fasta': group[0], 'coordinate_reference_accession': group[1],
                    'original_sample_fasta': original_sample_fasta,
                    'canonical_sample_fasta': (value(summary_row, '_canonical_species_fasta_path') or
                                               value(summary_row, '_canonical_coordinate_reference_fasta')),
                    'original_mitos2_fasta': group[0],
                    'canonical_mitos2_fasta': profile['canonical_fasta'],
                    'fasta_match': 'yes' if profile['fasta_match'] else 'no',
                    'group_row_count': profile['n_rows'], 'normalized_gene_count': profile['n_genes'],
                    'has_all_13_protein_coding_genes': 'yes' if profile['has_13_genes'] else 'no',
                    'coding_row_count_in_expected_range': 'yes' if profile['in_expected_range'] else 'no',
                    'fallback_match_mode': match_mode, 'n_candidate_rows': n_candidate_rows,
                    'n_selected_rows_before_dedup': n_candidate_rows,
                    'n_selected_rows_after_dedup': len(fallback) if is_selected else '',
                    'n_duplicate_rows_collapsed': duplicates_collapsed if is_selected else '',
                    'n_candidate_reference_groups': candidate_groups, 'selected_reference_group': selected_group,
                    'selection_status': ('ambiguous_selected' if ambiguous else 'selected') if is_selected else 'rejected',
                    'rejection_reason': profile['rejection_reason'], 'note': '; '.join(notes[1:]),
                })
    output = selected
    # A final successful fallback supersedes an earlier GenBank-resolution failure.
    final_status_by_sample = {value(row, 'sample'): value(row, 'status') for row in summary}
    failures = [row for row in failures
                if final_status_by_sample.get(value(row, 'sample')) not in SUCCESS_STATUSES]
    # Compare raw GenBank and MITOS2 CDS boundaries wherever both sources exist.
    comparison_path = yaml(args.config).get('mitos2_annotation', {}).get('paths', {}).get('genbank_mitos2_comparison_table', 'results/qc/codon_table_build/genbank_vs_mitos2_cds_comparison.tsv')
    comparison_fields = 'record_type sample original_species reference_species coordinate_reference_accession gene genbank_start genbank_end genbank_strand genbank_length mitos2_start mitos2_end mitos2_strand mitos2_length start_delta end_delta length_delta strand_match coordinate_match_status n_genbank_cds_genes n_mitos2_cds_genes n_gene_matches n_exact_boundary_matches n_near_boundary_matches n_strand_mismatches missing_in_genbank missing_in_mitos2 comparison_status'.split()
    comparisons = []
    for sample in sorted(set(genbank_rows_by_sample) | set(mitos_rows_by_sample)):
        g = {}; m = {}
        for row in genbank_rows_by_sample[sample]: g.setdefault(normalize_gene(value(row, 'gene')), []).append(row)
        for row in mitos_rows_by_sample[sample]: m.setdefault(normalize_gene(value(row, 'gene')), []).append(row)
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
    if any(not valid_phase(r.get('codon_pos_in_triplet')) for r in output): global_notes.append('Invalid codon_pos_in_triplet values detected.')
    if any(r['codon_seq'] and len(r['codon_seq']) != 3 for r in output): global_notes.append('Non-triplet codon_seq values detected.')
    if global_notes:
        for row in summary: row['note'] = '; '.join(filter(None, [row['note'], *global_notes]))
    if not args.dry_run: write_tsv(paths['output_table'], OUTPUT_FIELDS, output)
    fallback_summary_path = mitos_paths.get('mitos2_fallback_selection_summary_table', 'results/qc/codon_table_build/mitos2_fallback_selection_summary.tsv')
    write_tsv(fallback_summary_path, MITOS2_FALLBACK_SELECTION_FIELDS, fallback_selection_summary)
    write_tsv(paths['failed_downloads_table'], FAIL_FIELDS, failures); write_tsv(paths['summary_table'], SUMMARY_FIELDS, summary)
    if args.dry_run:
        resolved = sum(row['status'] == 'dry_run_resolved' for row in summary)
        unresolved = sum(row['status'] == 'dry_run_unresolved' for row in summary)
        print(f'Resolved accessions for {resolved} samples.')
        print(f'Unresolved accessions for {unresolved} samples.')
    else:
        counts = summarize_build_status(summary)
        assert counts['completed'] == counts['completed_genbank'] + counts['completed_mitos2_fallback']
        warn_if_output_summary_disagree(summary, output)
        print('Final annotation:', file=sys.stderr)
        print(f"  {counts['completed_genbank']} GenBank", file=sys.stderr)
        print(f"  {counts['completed_mitos2_fallback']} MITOS2 fallback", file=sys.stderr)
        print(f"  {counts['failed']} failed", file=sys.stderr)
        print(f"Built {len(output)} coding-position rows for {counts['completed']} samples "
              f"({counts['completed_genbank']} GenBank, {counts['completed_mitos2_fallback']} MITOS2 fallback; "
              f"{counts['failed']} failed, {counts['other']} other).")

if __name__ == '__main__': main()
