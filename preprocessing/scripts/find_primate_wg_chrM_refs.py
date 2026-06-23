#!/usr/bin/env python3
"""
Find whole-genome references and mitochondrial references for a list of species.

Version 5 changes:
  - Whole-genome assembly chrM hits are accepted only if their contig length is within
    the expected complete mitochondrial genome range. Short MT/mitochondrial fragments
    are recorded as invalid_length and the workflow continues to RefSeq mito FASTA,
    then nuccore, then fallback species.
  - WG-embedded chrM accessions ignore placeholder values such as na/none/- and prefer
    RefSeq accession, then GenBank accession, then contig name.
  - Output includes same_species_wg_chrM_status so invalid assembly chrM hits are visible.
  - Assembly-report chrM candidate selection now prioritizes complete-length records before
    RefSeq/GenBank accession preference. A short MT record can no longer mask a valid
    complete chrM record in the same assembly report.
  - A final safety check rejects any final chrM outside the expected complete mitogenome
    length range and continues to same-species RefSeq mito FASTA, nuccore, then fallback
    species.

Inputs:
  1) species list TSV/CSV with at least a species column, e.g. Chlorocebus_sabaeus
     Optional columns: sample_count, preprint_REFERENCE_SPECIES
  2) RefSeq mitochondrion genomic FASTA, e.g. mitochondrion.1.1.genomic.fna or .gz
  3) phylogeny tree in Newick format, with tip labels preferably as Genus_species

Reference search logic:
  For each target species:
    1. Search same-species whole-genome reference from NCBI RefSeq, NCBI GenBank, DNA Zoo.
       If found, check whether that whole-genome assembly contains chrM.
       If chrM exists in any same-species NCBI assembly, report the assembly/accession/contig where chrM was found.
       If no chrM in same-species assembly, search same-species chrM from local RefSeq mitochondrion FASTA, then NCBI nuccore.
       If still no chrM, search nearest species in the phylogeny until chrM is found.
       Final WG reference remains the same-species WG reference.

    2. If no same-species whole-genome reference exists:
       Search same-species chrM from local RefSeq mitochondrion FASTA, then NCBI nuccore.
       If no same-species chrM, search nearest species in the phylogeny until chrM is found.
       Then search nearest whole-genome reference from the phylogeny until a WG reference is found.

Outputs:
  outdir/species_reference_chrM_summary.tsv
  outdir/species_reference_chrM_summary.status_counts.tsv
  outdir/all_candidate_wg_refs.tsv
  outdir/nuccore_mito_hits.tsv

Dependencies:
  Python 3 standard library only.
  Internet is required for NCBI/DNA Zoo lookup unless cached files are already present.
"""

import argparse
import csv
import gzip
import json
import os
import re
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict, deque
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any

NCBI_REFSEQ_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
NCBI_GENBANK_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
DNAZOO_BASE = "https://dnazoo.s3.wasabisys.com/"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

BAD_ORG_RE = re.compile(
    r"virus|phage|plasmid|bacterium|bacteria|fungus|fungal|synthetic|construct|vector|uncultured|environmental|metagenome|associated",
    re.IGNORECASE,
)

CHR_M_RE = re.compile(r"mitochondr|\bchrM\b|\bMT\b|mitogenome", re.IGNORECASE)

# Nuccore mitochondrial reference hits must be full mitochondrial genomes.
# This intentionally rejects partial cytochrome-b/COX1/12S records and nuclear genes
# with mitochondrial-related names.
NUCCORE_BAD_TITLE_RE = re.compile(
    r"\b(partial|fragment|gene|cds|mRNA|nuclear|antiviral|pseudogene)\b|"
    r"control region|D-loop|cytochrome b|\bcytb\b|\bcox[0-9a-z]*\b|\bcoi\b|\bco1\b|\b12S\b|\b16S\b|\brrna\b|\btrna\b",
    re.IGNORECASE,
)
NUCCORE_COMPLETE_TITLE_RE = re.compile(
    r"complete mitochondrial genome|complete mitogenome|mitochondr(?:ion|ial).*complete genome|complete genome.*mitochondr(?:ion|ial)",
    re.IGNORECASE,
)
MIN_MITO_LEN = 14000
MAX_MITO_LEN = 25000

ASSEMBLY_LEVEL_SCORE = {
    "chromosome": 50,
    "complete genome": 45,
    "scaffold": 30,
    "contig": 20,
}

REFSEQ_CATEGORY_SCORE = {
    "reference genome": 20,
    "representative genome": 15,
    "na": 0,
    "": 0,
}


def log(msg: str) -> None:
    sys.stderr.write(f"[INFO] {msg}\n")
    sys.stderr.flush()


def warn(msg: str) -> None:
    sys.stderr.write(f"[WARN] {msg}\n")
    sys.stderr.flush()


def sanitize(x: Any) -> str:
    if x is None:
        return ""
    s = str(x)
    s = s.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def clean_field(x: Any) -> str:
    """Return a cleaned field, treating common assembly-report placeholders as empty."""
    s = sanitize(x)
    if s.lower() in {"", "na", "n/a", "none", "null", "-"}:
        return ""
    return s


def first_nonempty(*vals: Any) -> str:
    """First non-empty value after placeholder cleanup."""
    for v in vals:
        s = clean_field(v)
        if s:
            return s
    return ""


def valid_mito_length(x: Any) -> bool:
    """Complete primate mitochondrial genomes should normally be about 14-25 kb."""
    s = clean_field(x)
    if not s:
        return False
    try:
        # Some NCBI fields may contain commas.
        n = int(str(s).replace(",", ""))
    except Exception:
        return False
    return MIN_MITO_LEN <= n <= MAX_MITO_LEN


def norm_species_name(x: str) -> str:
    """Normalize species names for comparison: Genus_species lowercase."""
    if x is None:
        return ""
    s = str(x).strip().strip("'\"")
    s = s.replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s.lower()


def species_to_space(x: str) -> str:
    return re.sub(r"_+", " ", str(x).strip()).strip()


def first_two_words(x: str) -> str:
    toks = re.findall(r"[A-Za-z]+", str(x))
    if len(toks) >= 2:
        return f"{toks[0]}_{toks[1]}".lower()
    return norm_species_name(x)


def organism_matches_species(organism_name: str, target_species: str) -> bool:
    """Exact species-level match based on first two latin words."""
    if not organism_name or BAD_ORG_RE.search(organism_name):
        return False
    return first_two_words(organism_name) == norm_species_name(target_species)


def safe_urlretrieve(url: str, path: str, force: bool = False, timeout: int = 60) -> None:
    if os.path.exists(path) and os.path.getsize(path) > 0 and not force:
        return
    log(f"Downloading {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "primate-ref-check/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(path, "wb") as out:
        out.write(r.read())


def urlopen_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "primate-ref-check/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


@dataclass
class AssemblyHit:
    query_species: str
    matched_species: str
    source: str  # NCBI_RefSeq, NCBI_GenBank, DNAZoo
    assembly_accession: str
    bioproject: str = ""
    biosample: str = ""
    taxid: str = ""
    species_taxid: str = ""
    organism_name: str = ""
    refseq_category: str = ""
    assembly_level: str = ""
    genome_rep: str = ""
    seq_rel_date: str = ""
    asm_name: str = ""
    submitter: str = ""
    ftp_path: str = ""
    dna_zoo_prefix: str = ""
    score: int = 0
    has_chrM: str = "unknown"
    chrM_contig_name: str = ""
    chrM_sequence_name: str = ""
    chrM_genbank_accn: str = ""
    chrM_refseq_accn: str = ""
    chrM_ucsc_name: str = ""
    chrM_length: str = ""
    chrM_record: str = ""


@dataclass
class MitoHit:
    query_species: str
    matched_species: str
    source: str  # same_species_refseq_mito_fasta, same_species_nuccore, nearest_species_refseq_mito_fasta, nearest_species_nuccore, whole_genome_assembly
    similarity_rank: str = ""
    phylo_distance: str = ""
    accession: str = ""
    contig_name: str = ""
    title_or_header: str = ""
    length: str = ""
    assembly_accession: str = ""
    assembly_source: str = ""
    genbank_accn: str = ""
    refseq_accn: str = ""
    ucsc_name: str = ""


class MitoFastaIndex:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.by_species: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        self._build_index()

    def _open(self):
        if self.fasta_path.endswith(".gz"):
            return gzip.open(self.fasta_path, "rt", encoding="utf-8", errors="replace")
        return open(self.fasta_path, "rt", encoding="utf-8", errors="replace")

    def _extract_species_from_header(self, header: str) -> str:
        # Examples:
        # >NC_023958.1 Chlorocebus sabaeus mitochondrion, complete genome
        # >YP_... not expected for genomic fna
        text = header[1:] if header.startswith(">") else header
        toks = text.split()
        if len(toks) >= 3:
            # First token is accession; next two are Genus species.
            return f"{toks[1]}_{toks[2]}"
        return ""

    def _build_index(self):
        log(f"Indexing mitochondrial FASTA: {self.fasta_path}")
        current_header = None
        current_len = 0
        with self._open() as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if current_header is not None:
                        self._add_record(current_header, current_len)
                    current_header = line
                    current_len = 0
                else:
                    current_len += len(line.strip())
            if current_header is not None:
                self._add_record(current_header, current_len)
        log(f"Indexed mitochondrial FASTA species: {len(self.by_species)}")

    def _add_record(self, header: str, length: int):
        text = header[1:] if header.startswith(">") else header
        accession = text.split()[0] if text.split() else ""
        sp = self._extract_species_from_header(header)
        if sp:
            self.by_species[norm_species_name(sp)].append({
                "accession": accession,
                "header": sanitize(header),
                "length": str(length),
                "species": sp,
            })

    def find(self, species: str) -> Optional[MitoHit]:
        hits = self.by_species.get(norm_species_name(species), [])
        if not hits:
            return None
        # Prefer complete-looking and length near 16 kb.
        def score(h):
            title = h.get("header", "")
            length = int(h.get("length", "0") or 0)
            s = 0
            if re.search(r"complete", title, re.I):
                s += 20
            if re.search(r"mitochond", title, re.I):
                s += 10
            if 14000 <= length <= 25000:
                s += 10
            if h.get("accession", "").startswith("NC_"):
                s += 5
            return s
        valid_hits = [h for h in hits if valid_mito_length(h.get("length", ""))]
        pool = valid_hits if valid_hits else hits
        best = sorted(pool, key=score, reverse=True)[0]
        if not valid_mito_length(best.get("length", "")):
            warn(f"RefSeq mitochondrion FASTA hit for {species} has unusual length {best.get('length', '')}: {best.get('accession', '')}")
        return MitoHit(
            query_species=species,
            matched_species=norm_species_name(species),
            source="refseq_mito_fasta",
            accession=best.get("accession", ""),
            contig_name=best.get("accession", ""),
            title_or_header=best.get("header", ""),
            length=best.get("length", ""),
        )


class NewickTree:
    class Node:
        def __init__(self, name="", length=0.0):
            self.name = name
            self.length = length
            self.children = []
            self.parent = None

    def __init__(self, newick_path: str):
        self.newick_path = newick_path
        self.root = None
        self.tip_nodes = {}
        if newick_path:
            self._parse_file(newick_path)

    def _parse_file(self, path: str):
        text = open(path, "rt", encoding="utf-8", errors="replace").read().strip()
        self.root, pos = self._parse_subtree(text, 0)
        self._collect_tips(self.root)
        log(f"Parsed Newick tree tips: {len(self.tip_nodes)}")

    def _skip_ws(self, s, i):
        while i < len(s) and s[i].isspace():
            i += 1
        return i

    def _parse_label(self, s, i):
        i = self._skip_ws(s, i)
        if i < len(s) and s[i] in "'\"":
            quote = s[i]
            i += 1
            start = i
            while i < len(s) and s[i] != quote:
                i += 1
            label = s[start:i]
            i += 1 if i < len(s) else 0
            return label, i
        start = i
        while i < len(s) and s[i] not in ":,();":
            i += 1
        return s[start:i].strip(), i

    def _parse_length(self, s, i):
        i = self._skip_ws(s, i)
        if i < len(s) and s[i] == ":":
            i += 1
            start = i
            while i < len(s) and s[i] not in ",();":
                i += 1
            try:
                return float(s[start:i]), i
            except Exception:
                return 1.0, i
        return 1.0, i

    def _parse_subtree(self, s, i):
        i = self._skip_ws(s, i)
        if i < len(s) and s[i] == "(":
            i += 1
            node = self.Node()
            while True:
                child, i = self._parse_subtree(s, i)
                child.parent = node
                node.children.append(child)
                i = self._skip_ws(s, i)
                if i < len(s) and s[i] == ",":
                    i += 1
                    continue
                if i < len(s) and s[i] == ")":
                    i += 1
                    break
                break
            label, i = self._parse_label(s, i)
            node.name = label
            length, i = self._parse_length(s, i)
            node.length = length
            return node, i
        else:
            label, i = self._parse_label(s, i)
            length, i = self._parse_length(s, i)
            return self.Node(label, length), i

    def _collect_tips(self, node):
        if not node.children:
            raw_key = norm_species_name(node.name)
            species_key = first_two_words(node.name)
            # Store both exact tip label and genus_species alias. This fixes trees whose
            # tips contain extra labels, subspecies, or suffixes. Prefer the first tip
            # for a genus_species key to keep behavior stable.
            if raw_key and raw_key not in self.tip_nodes:
                self.tip_nodes[raw_key] = node
            if species_key and species_key not in self.tip_nodes:
                self.tip_nodes[species_key] = node
            return
        for c in node.children:
            self._collect_tips(c)

    def _canonical_tree_key(self, target_species: str) -> str:
        candidates = [norm_species_name(target_species), first_two_words(target_species)]
        for k in candidates:
            if k in self.tip_nodes:
                return k
        return ""

    def nearest_species(self, target_species: str) -> List[Tuple[str, int, float]]:
        """Return list of (tip_name_norm, rank, distance), ordered by patristic distance.

        If the exact target species is absent from the tree, fall back to species in
        the same genus. Those rows get a large synthetic distance but still allow
        nearest-reference search to continue instead of silently failing.
        """
        if not self.tip_nodes:
            return []
        target_key = self._canonical_tree_key(target_species)
        if not target_key:
            genus = norm_species_name(target_species).split("_")[0]
            genus_hits = sorted([k for k in self.tip_nodes if k.split("_")[0] == genus])
            if genus_hits:
                warn(f"Target species not found in tree tips: {target_species}; using same-genus tree fallback")
                maxn = getattr(self, "max_nearest", 200)
                return [(sp, i + 1, 999000.0 + i) for i, sp in enumerate(genus_hits[:maxn])]
            warn(f"Target species not found in tree tips and no same-genus fallback: {target_species}")
            return []
        start = self.tip_nodes[target_key]
        adj = defaultdict(list)
        def add_edges(n):
            for c in n.children:
                w = c.length if c.length is not None else 1.0
                adj[n].append((c, w))
                adj[c].append((n, w))
                add_edges(c)
        root = start
        while root.parent is not None:
            root = root.parent
        add_edges(root)
        import heapq
        dist = {start: 0.0}
        heap = [(0.0, id(start), start)]
        while heap:
            d, _, n = heapq.heappop(heap)
            if d != dist[n]:
                continue
            for nb, w in adj[n]:
                nd = d + (w if w is not None else 1.0)
                if nb not in dist or nd < dist[nb]:
                    dist[nb] = nd
                    heapq.heappush(heap, (nd, id(nb), nb))
        # De-duplicate aliases that point to the same Node.
        node_to_key = {}
        for key, node in self.tip_nodes.items():
            if key == target_key:
                continue
            if node in dist and node not in node_to_key:
                node_to_key[node] = key
        tips = [(key, dist[node]) for node, key in node_to_key.items()]
        tips.sort(key=lambda x: (x[1], x[0]))
        return [(sp, i + 1, d) for i, (sp, d) in enumerate(tips)]


class ReferenceFinder:
    def __init__(self, outdir: str, mito_index: MitoFastaIndex, tree: NewickTree,
                 email: str = "", api_key: str = "", force_download: bool = False,
                 max_nearest: int = 200, delay: float = 0.34):
        self.outdir = outdir
        self.cache_dir = os.path.join(outdir, "cache")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.mito_index = mito_index
        self.tree = tree
        self.email = email
        self.api_key = api_key
        self.max_nearest = max_nearest
        # Expose max_nearest to NewickTree same-genus fallback.
        setattr(self.tree, "max_nearest", max_nearest)
        self.delay = delay
        self.assemblies_by_species: Dict[str, List[AssemblyHit]] = defaultdict(list)
        self.assembly_report_cache: Dict[str, Optional[AssemblyHit]] = {}
        self.nuccore_cache: Dict[str, Optional[MitoHit]] = {}
        self.all_candidate_hits: List[AssemblyHit] = []
        self.nuccore_hits: List[MitoHit] = []
        self._cache_lock = threading.RLock()
        self._results_lock = threading.Lock()
        self._eutils_lock = threading.Lock()
        self._last_eutils_request = 0.0
        self._prepare_assembly_summaries(force_download)

    def _prepare_assembly_summaries(self, force: bool):
        refseq_path = os.path.join(self.cache_dir, "assembly_summary_refseq.txt")
        genbank_path = os.path.join(self.cache_dir, "assembly_summary_genbank.txt")
        safe_urlretrieve(NCBI_REFSEQ_SUMMARY, refseq_path, force=force)
        safe_urlretrieve(NCBI_GENBANK_SUMMARY, genbank_path, force=force)
        self._parse_assembly_summary(refseq_path, "NCBI_RefSeq")
        self._parse_assembly_summary(genbank_path, "NCBI_GenBank")
        log(f"Indexed NCBI species with assemblies: {len(self.assemblies_by_species)}")

    def _parse_assembly_summary(self, path: str, source: str):
        with open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 20:
                    continue
                organism = f[7]
                if BAD_ORG_RE.search(organism):
                    continue
                sp_key = first_two_words(organism)
                hit = AssemblyHit(
                    query_species="",
                    matched_species=sp_key,
                    source=source,
                    assembly_accession=f[0],
                    bioproject=f[1],
                    biosample=f[2],
                    taxid=f[5],
                    species_taxid=f[6],
                    organism_name=organism,
                    refseq_category=f[4],
                    assembly_level=f[11],
                    genome_rep=f[13],
                    seq_rel_date=f[14],
                    asm_name=f[15],
                    submitter=f[17],
                    ftp_path=f[19],
                )
                hit.score = self._assembly_score(hit)
                self.assemblies_by_species[sp_key].append(hit)

    def _assembly_score(self, h: AssemblyHit) -> int:
        s = 0
        if h.source == "NCBI_RefSeq":
            s += 1000
        elif h.source == "NCBI_GenBank":
            s += 800
        elif h.source == "DNAZoo":
            s += 600
        s += ASSEMBLY_LEVEL_SCORE.get(h.assembly_level.lower(), 0)
        s += REFSEQ_CATEGORY_SCORE.get(h.refseq_category.lower(), 0)
        if h.genome_rep.lower() == "full":
            s += 5
        return s

    def find_same_species_wg(self, species: str) -> List[AssemblyHit]:
        sp_key = norm_species_name(species)
        hits = []
        for h in self.assemblies_by_species.get(sp_key, []):
            hh = AssemblyHit(**asdict(h))
            hh.query_species = species
            hits.append(hh)
        # DNA Zoo only as supplementary; exact prefix check.
        dz = self.check_dnazoo(species)
        if dz:
            hits.append(dz)
        hits.sort(key=lambda x: x.score, reverse=True)
        return hits

    def best_wg(self, species: str) -> Optional[AssemblyHit]:
        hits = self.find_same_species_wg(species)
        if not hits:
            return None
        return hits[0]

    def check_dnazoo(self, species: str) -> Optional[AssemblyHit]:
        # DNA Zoo prefixes are not perfectly standardized. Try common variants only.
        sp = norm_species_name(species)
        genus, _, epithet = sp.partition("_")
        guesses = [
            f"{genus}_{epithet}/",
            f"{genus.capitalize()}_{epithet}/",
            f"{genus}.{epithet}/",
            f"{genus.capitalize()}.{epithet}/",
            f"{genus[:1].capitalize()}{epithet}/",
        ]
        for prefix in guesses:
            url = DNAZOO_BASE + "?" + urllib.parse.urlencode({"list-type": "2", "prefix": prefix, "max-keys": "5"})
            try:
                xml = urlopen_text(url, timeout=20)
                root = ET.fromstring(xml)
                keys = [e.text for e in root.iter() if e.tag.endswith("Key") and e.text]
                if keys:
                    return AssemblyHit(
                        query_species=species,
                        matched_species=sp,
                        source="DNAZoo",
                        assembly_accession=f"DNAZoo:{prefix.rstrip('/')}",
                        organism_name=species_to_space(species),
                        assembly_level="unknown",
                        ftp_path=DNAZOO_BASE + prefix,
                        dna_zoo_prefix=prefix,
                        score=600,
                        has_chrM="unknown",
                    )
            except Exception:
                continue
        return None

    def check_assembly_chrM(self, hit: AssemblyHit) -> AssemblyHit:
        if hit.source not in ("NCBI_RefSeq", "NCBI_GenBank") or not hit.ftp_path:
            hit.has_chrM = "unknown"
            return hit
        with self._cache_lock:
            cached_marker = self.assembly_report_cache.get(hit.assembly_accession, "__missing__")
        if cached_marker != "__missing__":
            cached = cached_marker
            if cached is not None:
                # copy chrM fields
                for k, v in asdict(cached).items():
                    if k.startswith("chrM") or k == "has_chrM":
                        setattr(hit, k, v)
            return hit
        base = os.path.basename(hit.ftp_path.rstrip("/"))
        report_url = hit.ftp_path.rstrip("/") + f"/{base}_assembly_report.txt"
        report_path = os.path.join(self.cache_dir, f"{hit.assembly_accession}_assembly_report.txt")
        try:
            safe_urlretrieve(report_url, report_path, force=False, timeout=45)
            chr_hit = self._parse_assembly_report_for_chrM(report_path)
            if chr_hit:
                # Store the best mitochondrial-looking record from the assembly report,
                # but accept it as a usable chrM only if the length is compatible with
                # a complete mitochondrial genome. This prevents short MT fragments
                # (for example ~200 bp records) from blocking fallback to RefSeq mito
                # FASTA or NCBI nuccore complete mitochondrial genomes.
                hit.chrM_sequence_name = clean_field(chr_hit.get("Sequence-Name", ""))
                hit.chrM_contig_name = first_nonempty(
                    chr_hit.get("Sequence-Name", ""),
                    chr_hit.get("RefSeq-Accn", ""),
                    chr_hit.get("GenBank-Accn", ""),
                )
                hit.chrM_genbank_accn = clean_field(chr_hit.get("GenBank-Accn", ""))
                hit.chrM_refseq_accn = clean_field(chr_hit.get("RefSeq-Accn", ""))
                hit.chrM_ucsc_name = clean_field(chr_hit.get("UCSC-style-name", ""))
                hit.chrM_length = clean_field(chr_hit.get("Sequence-Length", ""))
                hit.chrM_record = clean_field(chr_hit.get("raw", ""))
                if valid_mito_length(hit.chrM_length):
                    hit.has_chrM = "yes"
                else:
                    hit.has_chrM = "invalid_length"
                    warn(
                        f"Assembly {hit.assembly_accession} has mitochondrial-looking record "
                        f"{hit.chrM_contig_name} with length {hit.chrM_length}; "
                        "not accepting as complete chrM"
                    )
            else:
                hit.has_chrM = "no"
            with self._cache_lock:
                self.assembly_report_cache[hit.assembly_accession] = AssemblyHit(**asdict(hit))
        except Exception as e:
            hit.has_chrM = "assembly_report_unavailable"
            warn(f"Failed to check assembly report for {hit.assembly_accession}: {e}")
            with self._cache_lock:
                self.assembly_report_cache[hit.assembly_accession] = None
        return hit

    def _parse_assembly_report_for_chrM(self, path: str) -> Optional[Dict[str, str]]:
        cols = [
            "Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type",
            "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name"
        ]
        candidates = []
        with open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                raw = line.rstrip("\n")
                f = raw.split("\t")
                if len(f) < 8:
                    continue
                rec = {cols[i]: f[i] if i < len(f) else "" for i in range(len(cols))}
                rec["raw"] = sanitize(raw)
                joined = " ".join(f)
                role = rec.get("Sequence-Role", "")
                assigned = rec.get("Assigned-Molecule", "")
                asm_unit = rec.get("Assembly-Unit", "")
                seq_name = rec.get("Sequence-Name", "")
                if CHR_M_RE.search(joined) or assigned.upper() in ("MT", "M", "CHRM") or seq_name.lower() in ("mt", "chrm", "chrmt"):
                    candidates.append(rec)
                elif "non-nuclear" in role.lower() and "mitochond" in asm_unit.lower():
                    candidates.append(rec)
        if not candidates:
            return None
        # Choose candidate in two stages:
        #   1) If any mitochondrial-looking record has complete-mitogenome length,
        #      only consider those complete-length records.
        #   2) If no complete-length record exists, return the best invalid candidate
        #      only for reporting/status; check_assembly_chrM() will mark it as
        #      invalid_length and the workflow will continue to RefSeq mito FASTA /
        #      nuccore fallback.
        valid_candidates = [r for r in candidates if valid_mito_length(r.get("Sequence-Length", ""))]
        pool = valid_candidates if valid_candidates else candidates

        def score(r):
            s = 0
            # Complete-length records should dominate accession/name preference.
            if valid_mito_length(r.get("Sequence-Length", "")):
                s += 1000
            if r.get("RefSeq-Accn", "").startswith("NC_"):
                s += 100
            elif clean_field(r.get("RefSeq-Accn", "")):
                s += 50
            if clean_field(r.get("GenBank-Accn", "")):
                s += 30
            if "mitochond" in r.get("Assigned-Molecule-Location/Type", "").lower():
                s += 20
            if r.get("Assigned-Molecule", "").upper() in ("MT", "M", "CHRM"):
                s += 10
            return s
        return sorted(pool, key=score, reverse=True)[0]

    def find_wg_with_chrM_among_hits(self, hits: List[AssemblyHit]) -> Optional[AssemblyHit]:
        checked = []
        for h in hits:
            hc = self.check_assembly_chrM(h)
            checked.append(hc)
        # Which source has chrM? Prefer RefSeq with chrM, then GenBank with chrM.
        with_chrM = [h for h in checked if h.has_chrM == "yes"]
        if not with_chrM:
            return None
        with_chrM.sort(key=lambda x: x.score, reverse=True)
        return with_chrM[0]

    def is_usable_mito_hit(self, mh: Optional[MitoHit]) -> bool:
        """A final chrM must have complete-mitogenome length unless length is unavailable.

        For sources produced by this script, length should normally be available. If a
        whole-genome assembly record has a short mitochondrial-looking contig, reject it
        and keep searching. This is the key guard against 200 bp MT fragments becoming
        final chrM references.
        """
        if mh is None:
            return False
        if valid_mito_length(mh.length):
            return True
        # Do not accept known-length short/long records. Missing length is not treated
        # as usable here because all accepted FASTA/nuccore/WG candidates should have length.
        return False

    def make_mito_hit_from_wg_chrM(self, target: str, wg_chrM: Optional[AssemblyHit],
                                  source_label: str, rank: str = "0", dist: str = "0") -> Optional[MitoHit]:
        """Convert a WG assembly chrM record to MitoHit only if it has valid length."""
        if wg_chrM is None:
            return None
        if wg_chrM.has_chrM != "yes" or not valid_mito_length(wg_chrM.chrM_length):
            warn(
                f"Rejecting WG-embedded chrM for {target}: "
                f"assembly={wg_chrM.assembly_accession}, contig={wg_chrM.chrM_contig_name}, "
                f"length={wg_chrM.chrM_length}, status={wg_chrM.has_chrM}"
            )
            return None
        return MitoHit(
            query_species=target,
            matched_species=norm_species_name(target) if rank == "0" else wg_chrM.matched_species,
            source=source_label,
            similarity_rank=rank,
            phylo_distance=dist,
            accession=first_nonempty(wg_chrM.chrM_refseq_accn, wg_chrM.chrM_genbank_accn, wg_chrM.chrM_contig_name),
            contig_name=wg_chrM.chrM_contig_name,
            title_or_header=wg_chrM.chrM_record,
            length=wg_chrM.chrM_length,
            assembly_accession=wg_chrM.assembly_accession,
            assembly_source=wg_chrM.source,
            genbank_accn=wg_chrM.chrM_genbank_accn,
            refseq_accn=wg_chrM.chrM_refseq_accn,
            ucsc_name=wg_chrM.chrM_ucsc_name,
        )

    def ensure_final_chrM_usable(self, target: str, current: Optional[MitoHit],
                                 preprint_ref: str = "") -> Tuple[Optional[MitoHit], Dict[str, str], List[str]]:
        """Reject invalid final chrM and continue fallback search.

        Returns the usable MitoHit, nearest-search metadata if used, and notes.
        """
        info = {"target_in_tree": "", "nearest_search_mode": "", "nearest_fallback_reason": ""}
        notes = []
        if self.is_usable_mito_hit(current):
            return current, info, notes
        if current is not None:
            notes.append(
                f"rejected_invalid_final_chrM:{current.source}:{current.accession}:{current.length}"
            )
            warn(
                f"Final chrM candidate for {target} rejected because length is not complete-mitogenome range: "
                f"source={current.source}, accession={current.accession}, length={current.length}"
            )

        mh = self.search_same_species_chrM(target)
        if self.is_usable_mito_hit(mh):
            notes.append("final_chrM_recovered_from_same_species_refseq_or_nuccore")
            return mh, info, notes

        mh, info = self.search_nearest_chrM(target, preprint_ref)
        if self.is_usable_mito_hit(mh):
            notes.append("final_chrM_recovered_from_cross_species_fallback")
            return mh, info, notes
        if mh is not None:
            notes.append(f"cross_species_chrM_candidate_rejected_invalid_length:{mh.source}:{mh.accession}:{mh.length}")
        return None, info, notes

    def search_nuccore_mito(self, species: str, source_label: str) -> Optional[MitoHit]:
        key = norm_species_name(species)
        with self._cache_lock:
            cached_marker = self.nuccore_cache.get(key, "__missing__")
        if cached_marker != "__missing__":
            cached = cached_marker
            if cached:
                mh = MitoHit(**asdict(cached))
                mh.source = source_label
                mh.query_species = species
                return mh
            return None

        sp_space = species_to_space(species)
        # Strict query: complete mitochondrial genomes only. Do not fall back to
        # broad mitochondrial gene searches, because those often return partial
        # cytochrome-b/COX1/12S fragments or nuclear genes.
        terms = [
            f'"{sp_space}"[Organism] AND (mitochondrion[Title] OR mitochondrial[Title] OR mitogenome[Title]) AND ("complete genome"[Title] OR "complete mitochondrial genome"[Title] OR "complete mitogenome"[Title]) NOT partial[Title] NOT fragment[Title] NOT gene[Title] NOT cds[Title] NOT cytb[Title] NOT COX[Title] NOT 12S[Title] NOT 16S[Title] NOT rRNA[Title] NOT tRNA[Title]',
        ]
        best = None
        rejected = []
        for term in terms:
            ids = self._eutils_esearch(term, retmax=50)
            if not ids:
                continue
            summaries = self._eutils_esummary(ids)
            candidates = []
            for s in summaries:
                title = s.get("Title") or s.get("title") or ""
                acc = s.get("AccessionVersion") or s.get("accessionversion") or s.get("Caption") or s.get("caption") or ""
                slen = str(s.get("Length") or s.get("slen") or s.get("Length", ""))
                cand = {"accession": acc, "title": title, "length": slen}
                ok, reason = self._is_complete_mito_nuccore(cand)
                if ok:
                    candidates.append(cand)
                else:
                    rejected.append((cand, reason))
            if candidates:
                best = self._choose_best_nuccore(candidates)
                break
        if not best:
            if rejected:
                r0, reason = rejected[0]
                warn(f"nuccore rejected non-complete mito hit for {species}: {r0.get('accession')} {reason} {r0.get('title')}")
            with self._cache_lock:
                self.nuccore_cache[key] = None
            return None
        mh_cache = MitoHit(
            query_species=species,
            matched_species=key,
            source="nuccore",
            accession=best["accession"],
            contig_name=best["accession"],
            title_or_header=best["title"],
            length=best["length"],
        )
        with self._cache_lock:
            self.nuccore_cache[key] = mh_cache
        mh = MitoHit(**asdict(mh_cache))
        mh.source = source_label
        with self._results_lock:
            self.nuccore_hits.append(mh)
        return mh

    def _is_complete_mito_nuccore(self, h: Dict[str, str]) -> Tuple[bool, str]:
        title = h.get("title", "")
        acc = h.get("accession", "")
        try:
            length = int(h.get("length") or 0)
        except Exception:
            length = 0
        if BAD_ORG_RE.search(title):
            return False, "bad_organism_title"
        if NUCCORE_BAD_TITLE_RE.search(title):
            return False, "partial_or_gene_like_title"
        if not NUCCORE_COMPLETE_TITLE_RE.search(title):
            return False, "not_complete_mito_title"
        if not (MIN_MITO_LEN <= length <= MAX_MITO_LEN):
            return False, f"length_outside_{MIN_MITO_LEN}_{MAX_MITO_LEN}"
        if not acc:
            return False, "missing_accession"
        return True, "ok"

    def _wait_for_eutils_slot(self) -> None:
        """Rate-limit NCBI E-utility requests across worker threads."""
        if self.delay <= 0:
            return
        with self._eutils_lock:
            now = time.monotonic()
            wait = self.delay - (now - self._last_eutils_request)
            if wait > 0:
                time.sleep(wait)
            self._last_eutils_request = time.monotonic()

    def _eutils_esearch(self, term: str, retmax: int = 20) -> List[str]:
        params = {"db": "nuccore", "term": term, "retmode": "json", "retmax": str(retmax)}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
        try:
            self._wait_for_eutils_slot()
            data = json.loads(urlopen_text(url, timeout=45))
            return data.get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            warn(f"nuccore esearch failed for term={term}: {e}")
            return []

    def _eutils_esummary(self, ids: List[str]) -> List[Dict[str, Any]]:
        if not ids:
            return []
        params = {"db": "nuccore", "id": ",".join(ids), "retmode": "json"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS_BASE}/esummary.fcgi?" + urllib.parse.urlencode(params)
        try:
            self._wait_for_eutils_slot()
            data = json.loads(urlopen_text(url, timeout=45))
            result = data.get("result", {})
            return [result[i] for i in result.get("uids", []) if i in result]
        except Exception as e:
            warn(f"nuccore esummary failed for ids={ids[:3]}...: {e}")
            return []

    def _choose_best_nuccore(self, hits: List[Dict[str, str]]) -> Dict[str, str]:
        # Called only after strict complete-mito filtering. Prefer RefSeq NC_, then
        # length closest to a typical primate mitogenome (~16.5 kb).
        def score(h):
            acc = h.get("accession", "")
            try:
                length = int(h.get("length") or 0)
            except Exception:
                length = 0
            s = 0.0
            if acc.startswith("NC_"):
                s += 100.0
            s -= abs(length - 16500) / 1000.0
            return s
        return sorted(hits, key=score, reverse=True)[0]

    def search_same_species_chrM(self, species: str, preferred_prefix: str = "same_species") -> Optional[MitoHit]:
        mh = self.mito_index.find(species)
        if mh:
            mh.source = f"{preferred_prefix}_refseq_mito_fasta"
            mh.query_species = species
            mh.matched_species = norm_species_name(species)
            mh.similarity_rank = "0"
            mh.phylo_distance = "0"
            return mh
        mh = self.search_nuccore_mito(species, f"{preferred_prefix}_nuccore")
        if mh:
            mh.similarity_rank = "0"
            mh.phylo_distance = "0"
            return mh
        return None

    def _nearest_candidates_with_fallback(self, target_species: str, reference_species: str = "") -> Tuple[List[Tuple[str, str, str]], Dict[str, str]]:
        """Return candidate species for cross-species fallback with explicit mode.

        Priority when target is absent from tree:
          1. same-genus proxy species already present in the tree
          2. REFERENCE_SPECIES / preprint_REFERENCE_SPECIES fallback

        Returns:
          candidates: list of (species_key, rank, distance)
          info: metadata for summary output
        """
        target_key = norm_species_name(target_species)
        canonical = self.tree._canonical_tree_key(target_species) if self.tree else ""
        info = {
            "target_in_tree": "yes" if canonical else "no",
            "nearest_search_mode": "none",
            "nearest_fallback_reason": "",
        }

        ref_key = norm_species_name(reference_species)

        if canonical:
            nearest = self.tree.nearest_species(target_species)[: self.max_nearest]
            candidates = [(sp, str(rank), f"{dist:.8g}") for sp, rank, dist in nearest]
            if ref_key and ref_key not in {sp for sp, _, _ in candidates} and ref_key != target_key:
                candidates.append((ref_key, "REFERENCE_SPECIES", "NA"))
            info["nearest_search_mode"] = "phylogeny_distance"
            info["nearest_fallback_reason"] = "target_species_found_in_tree"
            return candidates, info

        genus = target_key.split("_")[0] if target_key else ""
        if genus and getattr(self.tree, "tip_nodes", None):
            genus_hits = sorted([k for k in self.tree.tip_nodes if k.split("_")[0] == genus])
            # De-duplicate aliases pointing to the same tree node.
            seen_nodes = set()
            dedup = []
            for k in genus_hits:
                node = self.tree.tip_nodes[k]
                if id(node) in seen_nodes:
                    continue
                seen_nodes.add(id(node))
                dedup.append(k)
            if dedup:
                candidates = [(sp, str(i + 1), f"same_genus_proxy_{i + 1}") for i, sp in enumerate(dedup[: self.max_nearest])]
                if ref_key and ref_key not in {sp for sp, _, _ in candidates} and ref_key != target_key:
                    candidates.append((ref_key, "REFERENCE_SPECIES", "NA"))
                info["nearest_search_mode"] = "same_genus_proxy"
                info["nearest_fallback_reason"] = "target_species_absent_from_tree_using_same_genus_tree_tips"
                return candidates, info

        if ref_key:
            info["nearest_search_mode"] = "REFERENCE_SPECIES_fallback"
            info["nearest_fallback_reason"] = "target_species_absent_from_tree_no_same_genus_using_REFERENCE_SPECIES"
            return [(ref_key, "REFERENCE_SPECIES", "NA")], info

        info["nearest_fallback_reason"] = "target_species_absent_from_tree_no_same_genus_no_REFERENCE_SPECIES"
        return [], info

    def search_nearest_chrM(self, target_species: str, reference_species: str = "") -> Tuple[Optional[MitoHit], Dict[str, str]]:
        candidates, info = self._nearest_candidates_with_fallback(target_species, reference_species)
        for sp_key, rank, dist in candidates:
            mh = self.mito_index.find(sp_key)
            if mh:
                mh.query_species = target_species
                mh.matched_species = sp_key
                if rank == "REFERENCE_SPECIES" or info.get("nearest_search_mode") == "REFERENCE_SPECIES_fallback":
                    mh.source = "REFERENCE_SPECIES_refseq_mito_fasta"
                elif info.get("nearest_search_mode") == "same_genus_proxy":
                    mh.source = "same_genus_proxy_refseq_mito_fasta"
                else:
                    mh.source = "nearest_species_refseq_mito_fasta"
                mh.similarity_rank = str(rank)
                mh.phylo_distance = str(dist)
                return mh, info
            source_prefix = "nearest_species_nuccore"
            if rank == "REFERENCE_SPECIES" or info.get("nearest_search_mode") == "REFERENCE_SPECIES_fallback":
                source_prefix = "REFERENCE_SPECIES_nuccore"
            elif info.get("nearest_search_mode") == "same_genus_proxy":
                source_prefix = "same_genus_proxy_nuccore"
            mh = self.search_nuccore_mito(sp_key, source_prefix)
            if mh:
                mh.query_species = target_species
                mh.matched_species = sp_key
                mh.similarity_rank = str(rank)
                mh.phylo_distance = str(dist)
                return mh, info
        return None, info

    def search_nearest_wg(self, target_species: str, reference_species: str = "") -> Tuple[Optional[AssemblyHit], str, str, Dict[str, str]]:
        candidates, info = self._nearest_candidates_with_fallback(target_species, reference_species)
        for sp_key, rank, dist in candidates:
            hits = self.find_same_species_wg(sp_key)
            if hits:
                best = hits[0]
                best.query_species = target_species
                best.matched_species = sp_key
                return best, str(rank), str(dist), info
        return None, "", "", info

    def analyze_species(self, row: Dict[str, str]) -> Dict[str, str]:
        target = row["species"]
        preprint_ref = row.get("preprint_REFERENCE_SPECIES", "") or row.get("REFERENCE_SPECIES", "")
        sample_count = row.get("sample_count", "") or row.get("samples", "") or row.get("n_samples", "")

        same_hits = self.find_same_species_wg(target)
        if same_hits:
            with self._results_lock:
                self.all_candidate_hits.extend(same_hits)

        same_wg_best = same_hits[0] if same_hits else None
        same_wg_chrM = self.find_wg_with_chrM_among_hits(same_hits) if same_hits else None

        final_wg = None
        final_wg_rank = ""
        final_wg_dist = ""
        final_chrM = None
        nearest_chrM_info = {"target_in_tree": "", "nearest_search_mode": "", "nearest_fallback_reason": ""}
        nearest_wg_info = {"target_in_tree": "", "nearest_search_mode": "", "nearest_fallback_reason": ""}
        notes = []

        if same_wg_best:
            final_wg = same_wg_best
            final_wg_rank = "0"
            final_wg_dist = "0"
            if same_wg_chrM:
                final_chrM = self.make_mito_hit_from_wg_chrM(
                    target,
                    same_wg_chrM,
                    source_label="same_species_whole_genome_assembly",
                    rank="0",
                    dist="0",
                )
            else:
                final_chrM = self.search_same_species_chrM(target)
                if not final_chrM:
                    final_chrM, nearest_chrM_info = self.search_nearest_chrM(target, preprint_ref)
                    if final_chrM:
                        notes.append("same_species_wg_found_but_chrM_from_cross_species_fallback")
        else:
            final_chrM = self.search_same_species_chrM(target)
            if not final_chrM:
                final_chrM, nearest_chrM_info = self.search_nearest_chrM(target, preprint_ref)
            final_wg, final_wg_rank, final_wg_dist, nearest_wg_info = self.search_nearest_wg(target, preprint_ref)
            if not final_wg:
                notes.append("no_same_or_cross_species_wg_ref_found")

        # Optional: preprint reference fallback only if no same/nearest WG found.
        # This helps when the phylogeny tree does not include the preprint reference species.
        preprint_wg = None
        preprint_wg_chrM = None
        if preprint_ref:
            pre_hits = self.find_same_species_wg(preprint_ref)
            if pre_hits:
                preprint_wg = pre_hits[0]
                preprint_wg_chrM = self.find_wg_with_chrM_among_hits(pre_hits)
                if final_wg is None:
                    final_wg = preprint_wg
                    final_wg_rank = "preprint_REFERENCE_SPECIES"
                    final_wg_dist = "NA"
                    notes.append("final_wg_used_preprint_REFERENCE_SPECIES_fallback")
                if final_chrM is None and preprint_wg_chrM:
                    final_chrM = self.make_mito_hit_from_wg_chrM(
                        target,
                        preprint_wg_chrM,
                        source_label="preprint_reference_species_whole_genome_assembly",
                        rank="preprint_REFERENCE_SPECIES",
                        dist="NA",
                    )
                    if final_chrM is not None:
                        final_chrM.matched_species = norm_species_name(preprint_ref)
                        notes.append("final_chrM_used_preprint_REFERENCE_SPECIES_fallback")
                elif final_chrM is None:
                    # Still try preprint reference species local mt then nuccore.
                    mh = self.search_same_species_chrM(preprint_ref, preferred_prefix="preprint_reference_species")
                    if mh:
                        mh.query_species = target
                        mh.matched_species = norm_species_name(preprint_ref)
                        mh.similarity_rank = "preprint_REFERENCE_SPECIES"
                        mh.phylo_distance = "NA"
                        final_chrM = mh
                        notes.append("final_chrM_used_preprint_REFERENCE_SPECIES_fallback")

        # Final safety net: no record outside the complete mitochondrial genome length
        # range can be emitted as final_chrM. If a short WG MT fragment slipped through
        # any previous step, reject it here and continue fallback searching.
        final_chrM, recovered_chrM_info, recovered_notes = self.ensure_final_chrM_usable(target, final_chrM, preprint_ref)
        if recovered_chrM_info.get("nearest_search_mode"):
            nearest_chrM_info = recovered_chrM_info
        notes.extend(recovered_notes)

        final_strategy = self._strategy(same_wg_best, same_wg_chrM, final_wg, final_chrM)

        final_chrM_accession = ""
        if final_chrM:
            final_chrM_accession = first_nonempty(final_chrM.refseq_accn, final_chrM.genbank_accn, final_chrM.accession, final_chrM.contig_name)

        out = {
            "target_species": target,
            "sample_count": sample_count,
            "preprint_REFERENCE_SPECIES": preprint_ref,
            "target_in_tree_for_chrM_search": nearest_chrM_info.get("target_in_tree", ""),
            "chrM_nearest_search_mode": nearest_chrM_info.get("nearest_search_mode", ""),
            "chrM_nearest_fallback_reason": nearest_chrM_info.get("nearest_fallback_reason", ""),
            "target_in_tree_for_wg_search": nearest_wg_info.get("target_in_tree", ""),
            "wg_nearest_search_mode": nearest_wg_info.get("nearest_search_mode", ""),
            "wg_nearest_fallback_reason": nearest_wg_info.get("nearest_fallback_reason", ""),

            "same_species_wg_ref_found": "yes" if same_wg_best else "no",
            "same_species_wg_ref_source": same_wg_best.source if same_wg_best else "",
            "same_species_wg_assembly_accession": same_wg_best.assembly_accession if same_wg_best else "",
            "same_species_wg_assembly_level": same_wg_best.assembly_level if same_wg_best else "",
            "same_species_wg_organism_name": same_wg_best.organism_name if same_wg_best else "",
            "same_species_wg_ftp_path": same_wg_best.ftp_path if same_wg_best else "",
            "same_species_wg_has_chrM": "yes" if same_wg_chrM else ("no" if same_wg_best else ""),
            "same_species_wg_chrM_status": same_wg_chrM.has_chrM if same_wg_chrM else (same_wg_best.has_chrM if same_wg_best else ""),
            "same_species_wg_chrM_assembly_source": same_wg_chrM.source if same_wg_chrM else "",
            "same_species_wg_chrM_assembly_accession": same_wg_chrM.assembly_accession if same_wg_chrM else "",
            "same_species_wg_chrM_contig_name": same_wg_chrM.chrM_contig_name if same_wg_chrM else "",
            "same_species_wg_chrM_genbank_accn": same_wg_chrM.chrM_genbank_accn if same_wg_chrM else "",
            "same_species_wg_chrM_refseq_accn": same_wg_chrM.chrM_refseq_accn if same_wg_chrM else "",
            "same_species_wg_chrM_length": same_wg_chrM.chrM_length if same_wg_chrM else "",

            "final_wg_ref_species": final_wg.matched_species if final_wg else "",
            "final_wg_ref_source": final_wg.source if final_wg else "",
            "final_wg_ref_similarity_rank": final_wg_rank,
            "final_wg_ref_phylo_distance": final_wg_dist,
            "final_wg_assembly_accession": final_wg.assembly_accession if final_wg else "",
            "final_wg_assembly_level": final_wg.assembly_level if final_wg else "",
            "final_wg_organism_name": final_wg.organism_name if final_wg else "",
            "final_wg_ftp_path": final_wg.ftp_path if final_wg else "",

            "final_chrM_species": final_chrM.matched_species if final_chrM else "",
            "final_chrM_source": final_chrM.source if final_chrM else "",
            "final_chrM_similarity_rank": final_chrM.similarity_rank if final_chrM else "",
            "final_chrM_phylo_distance": final_chrM.phylo_distance if final_chrM else "",
            "final_chrM_accession": final_chrM_accession,
            "final_chrM_contig_name": final_chrM.contig_name if final_chrM else "",
            "final_chrM_title_or_header": final_chrM.title_or_header if final_chrM else "",
            "final_chrM_length": final_chrM.length if final_chrM else "",
            "final_chrM_assembly_accession": final_chrM.assembly_accession if final_chrM else "",
            "final_chrM_assembly_source": final_chrM.assembly_source if final_chrM else "",
            "final_chrM_genbank_accn": final_chrM.genbank_accn if final_chrM else "",
            "final_chrM_refseq_accn": final_chrM.refseq_accn if final_chrM else "",
            "final_chrM_ucsc_name": final_chrM.ucsc_name if final_chrM else "",

            "preprint_ref_wg_found": "yes" if preprint_wg else ("no" if preprint_ref else ""),
            "preprint_ref_wg_source": preprint_wg.source if preprint_wg else "",
            "preprint_ref_wg_assembly_accession": preprint_wg.assembly_accession if preprint_wg else "",
            "preprint_ref_wg_assembly_level": preprint_wg.assembly_level if preprint_wg else "",
            "preprint_ref_wg_has_chrM": "yes" if preprint_wg_chrM else ("no" if preprint_wg else ""),
            "preprint_ref_wg_chrM_status": preprint_wg_chrM.has_chrM if preprint_wg_chrM else (preprint_wg.has_chrM if preprint_wg else ""),
            "preprint_ref_wg_chrM_contig_name": preprint_wg_chrM.chrM_contig_name if preprint_wg_chrM else "",
            "preprint_ref_wg_chrM_refseq_accn": preprint_wg_chrM.chrM_refseq_accn if preprint_wg_chrM else "",
            "preprint_ref_wg_chrM_genbank_accn": preprint_wg_chrM.chrM_genbank_accn if preprint_wg_chrM else "",

            "final_reference_strategy": final_strategy,
            "notes": ";".join(notes),
        }
        return {k: clean_field(v) for k, v in out.items()}

    def _strategy(self, same_wg_best, same_wg_chrM, final_wg, final_chrM) -> str:
        if same_wg_best and same_wg_chrM:
            return "same_species_wg_with_chrM"
        if same_wg_best and final_chrM:
            if final_chrM.source == "same_species_refseq_mito_fasta":
                return "same_species_wg_plus_same_species_refseq_chrM"
            if final_chrM.source == "same_species_nuccore":
                return "same_species_wg_plus_same_species_nuccore_chrM"
            if final_chrM.source.startswith("nearest_species"):
                return "same_species_wg_plus_nearest_species_chrM"
            if final_chrM.source.startswith("preprint_reference_species"):
                return "same_species_wg_plus_preprint_reference_species_chrM"
            return "same_species_wg_plus_other_chrM"
        if (not same_wg_best) and final_wg and final_chrM:
            if final_chrM.similarity_rank == "0":
                return "nearest_species_wg_plus_same_species_chrM"
            return "nearest_species_wg_plus_nearest_species_chrM"
        if final_chrM and not final_wg:
            return "no_wg_ref_but_chrM_found"
        if final_wg and not final_chrM:
            return "wg_ref_found_but_no_chrM_found"
        return "no_reference_found"


def read_species_table(path: str) -> List[Dict[str, str]]:
    # Auto-detect delimiter.
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(4096)
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    rows = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise ValueError("Input species table has no header")
        fields = reader.fieldnames
        species_col = None
        for c in ["species", "GENUS_SPECIES", "target_species", "FINAL_PRIMATE_NAME"]:
            if c in fields:
                species_col = c
                break
        if species_col is None:
            raise ValueError(f"Cannot find species column. Expected one of: species, GENUS_SPECIES, target_species, FINAL_PRIMATE_NAME. Found: {fields}")
        for r in reader:
            sp = r.get(species_col, "").strip()
            if not sp:
                continue
            r["species"] = sp
            rows.append(r)
    # De-duplicate by species, keeping first row.
    seen = set()
    out = []
    for r in rows:
        k = norm_species_name(r["species"])
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out


def write_tsv(path: str, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None):
    if fieldnames is None:
        fieldnames = []
        for r in rows:
            for k in r.keys():
                if k not in fieldnames:
                    fieldnames.append(k)
    with open(path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: sanitize(r.get(k, "")) for k in fieldnames})


def main():
    ap = argparse.ArgumentParser(description="Find WG references and chrM references for species list")
    ap.add_argument("--species", required=True, help="Species list TSV/CSV; must contain species or GENUS_SPECIES column")
    ap.add_argument("--mito-fasta", required=True, help="RefSeq mitochondrion genomic FASTA, can be .gz")
    ap.add_argument("--tree", required=True, help="Newick tree with species tip names")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--email", default="", help="Email for NCBI E-utilities")
    ap.add_argument("--api-key", default="", help="NCBI API key for E-utilities")
    ap.add_argument("--max-nearest", type=int, default=200, help="Maximum nearest species to try in phylogeny")
    ap.add_argument("--force-download", action="store_true", help="Re-download NCBI assembly summaries")
    ap.add_argument("--delay", type=float, default=0.34, help="Delay between NCBI E-utility requests")
    ap.add_argument("--threads", type=int, default=1, help="Number of species to analyze concurrently; NCBI E-utility requests remain rate-limited by --delay")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    species_rows = read_species_table(args.species)
    log(f"Loaded target species: {len(species_rows)}")
    mito_index = MitoFastaIndex(args.mito_fasta)
    tree = NewickTree(args.tree)
    finder = ReferenceFinder(
        outdir=args.outdir,
        mito_index=mito_index,
        tree=tree,
        email=args.email,
        api_key=args.api_key,
        force_download=args.force_download,
        max_nearest=args.max_nearest,
        delay=args.delay,
    )

    def analyze_one(i: int, row: Dict[str, str]) -> Tuple[int, Dict[str, str]]:
        log(f"[{i}/{len(species_rows)}] {row['species']}")
        try:
            return i - 1, finder.analyze_species(row)
        except Exception as e:
            warn(f"Failed species {row.get('species')}: {e}")
            return i - 1, {
                "target_species": row.get("species", ""),
                "final_reference_strategy": "error",
                "notes": sanitize(str(e)),
            }

    threads = max(1, args.threads)
    summary_rows = [None] * len(species_rows)
    if threads == 1 or len(species_rows) <= 1:
        for i, row in enumerate(species_rows, 1):
            idx, result = analyze_one(i, row)
            summary_rows[idx] = result
    else:
        log(f"Analyzing species with {threads} worker threads")
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(analyze_one, i, row) for i, row in enumerate(species_rows, 1)]
            for future in as_completed(futures):
                idx, result = future.result()
                summary_rows[idx] = result
    summary_rows = [r for r in summary_rows if r is not None]

    summary_path = os.path.join(args.outdir, "species_reference_chrM_summary.tsv")
    write_tsv(summary_path, summary_rows)

    # Counts by final strategy.
    counts = defaultdict(int)
    for r in summary_rows:
        counts[r.get("final_reference_strategy", "") or "NA"] += 1
    count_rows = [{"final_reference_strategy": k, "n_species": v} for k, v in sorted(counts.items())]
    write_tsv(os.path.join(args.outdir, "species_reference_chrM_summary.status_counts.tsv"), count_rows)

    # Candidate WG references.
    cand_rows = [asdict(h) for h in finder.all_candidate_hits]
    if cand_rows:
        write_tsv(os.path.join(args.outdir, "all_candidate_wg_refs.tsv"), cand_rows)

    # Nuccore hits.
    nuc_rows = [asdict(h) for h in finder.nuccore_hits]
    if nuc_rows:
        write_tsv(os.path.join(args.outdir, "nuccore_mito_hits.tsv"), nuc_rows)

    log(f"Done. Main output: {summary_path}")


if __name__ == "__main__":
    main()
