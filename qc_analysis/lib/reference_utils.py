"""Reference sequence identity helpers shared by QC annotation tools."""
from __future__ import annotations
import gzip
import hashlib
from pathlib import Path

_FASTA_HASH_CACHE = {}

def normalized_sequence_sha256(sequence):
    """Return SHA256 for uppercase, whitespace-free sequence bases."""
    normalized = ''.join(str(sequence).split()).upper()
    if not normalized:
        raise ValueError('Sequence is empty after normalization.')
    return hashlib.sha256(normalized.encode('ascii')).hexdigest()

def normalized_fasta_sequence_sha256(path):
    """Hash FASTA bases only and return provenance; supports gzipped FASTA.

    Results are cached by resolved path for the lifetime of the process.
    """
    raw = str(path or '')
    if not raw:
        raise FileNotFoundError('Required FASTA path is empty.')
    resolved = Path(raw).expanduser().resolve()
    if resolved in _FASTA_HASH_CACHE:
        return dict(_FASTA_HASH_CACHE[resolved])
    if not resolved.is_file():
        raise FileNotFoundError(f'Required FASTA file is missing: {resolved}')
    opener = gzip.open if resolved.name.endswith('.gz') else open
    pieces = []
    with opener(resolved, 'rt') as handle:
        for line in handle:
            if not line.startswith('>'):
                pieces.append(''.join(line.split()))
    sequence = ''.join(pieces).upper()
    if not sequence:
        raise ValueError(f'FASTA contains no sequence bases: {resolved}')
    result = {'fasta_path': str(resolved), 'sequence_length': len(sequence),
              'sequence_sha256': normalized_sequence_sha256(sequence),
              'hash_source': 'normalized_fasta_sequence', 'hash_status': 'ok'}
    _FASTA_HASH_CACHE[resolved] = result
    return dict(result)
