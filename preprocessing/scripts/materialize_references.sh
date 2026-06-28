#!/usr/bin/env bash
set -euo pipefail
MANIFEST=${1:-references/manifests/reference_materialization_manifest.tsv}
OUTDIR=${OUTDIR:-results/preprocessing/reference_materialization}
LOCAL_MITO_FASTA=${LOCAL_MITO_FASTA:-}
mkdir -p "$OUTDIR" references/wg references/chrM/embedded_from_wg references/chrM/independent references/chrM/dnazoo references/manifests
DL="$OUTDIR/reference_download_manifest.tsv"; EX="$OUTDIR/chrM_extraction_manifest.tsv"; IH="$OUTDIR/in_house_score_reference_inputs.tsv"; CC="$OUTDIR/chrM_candidate_check.tsv"; RESOLVED_RESULTS="$OUTDIR/reference_materialization_manifest.resolved.tsv"; RESOLVED_REFS="references/manifests/reference_materialization_manifest.resolved.tsv"
PYTHON_COMMAND=${PYTHON_COMMAND:-python3}
WGET_COMMAND=${WGET_COMMAND:-wget}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
CURL_COMMAND=${CURL_COMMAND:-curl}
EFETCH_COMMAND=${EFETCH_COMMAND:-efetch}
export WGET_COMMAND SAMTOOLS_COMMAND CURL_COMMAND EFETCH_COMMAND
echo -e "target_species\tassembly_accession\tstatus\twg_fasta_path\twg_fai_path\twg_assembly_report_path\tmessage" > "$DL"
echo -e "target_species\tchrM_reference_context\tstatus\tchrM_fasta_path\tchrM_fai_path\tmessage" > "$EX"
echo -e "target_species\tattempted_assembly_accession\tattempt_type\twg_fasta_path\tassembly_report_path\tcandidate_names\tcandidates_found_in_fai\tcandidates_missing_from_fai\tstatus\tfallback_action" > "$CC"
"$PYTHON_COMMAND" - "$MANIFEST" "$DL" "$EX" "$CC" "$RESOLVED_RESULTS" "$RESOLVED_REFS" <<'PY'
import csv, gzip, os, re, shutil, subprocess, sys
man, dl, ex, cc, resolved_results, resolved_refs = sys.argv[1:]


def run(cmd):
    subprocess.check_call(cmd)


def run_with_output(cmd):
    try:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        output = (e.output or "").strip().replace("\t", " ").replace("\n", " | ")
        raise RuntimeError(f"{' '.join(cmd)} failed with exit {e.returncode}: {output}") from e


def download_file(url, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".download"
    if os.path.exists(tmp):
        os.remove(tmp)
    run([os.environ.get("WGET_COMMAND", "wget"), "-O", tmp, url])
    os.replace(tmp, out_path)


def gunzip_to_fasta(gz_path, fasta_path):
    tmp = fasta_path + ".tmp"
    with gzip.open(gz_path, "rb") as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.replace(tmp, fasta_path)


def normalize_species(value):
    return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")


def reference_pairing_status(target, wg_species, chrM_species):
    t = normalize_species(target)
    wg = normalize_species(wg_species)
    mt = normalize_species(chrM_species)
    has_wg = bool(wg)
    has_mt = bool(mt)
    if has_wg and has_mt and wg == t and mt == t:
        return "same_species_wg_same_species_chrM"
    if has_wg and has_mt and wg == t and mt != t:
        return "same_species_wg_cross_species_chrM"
    if has_wg and has_mt and wg != t and mt == t:
        return "cross_species_wg_same_species_chrM"
    if has_wg and has_mt and wg != t and mt != t:
        return "cross_species_wg_cross_species_chrM"
    if has_wg and not has_mt:
        return "wg_only_no_chrM"
    if not has_wg and has_mt:
        return "chrM_only_no_wg"
    return "no_reference_found"


def dna_zoo_species_token(value):
    raw = str(value or "").strip()
    token = re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_")
    return token or normalize_species(raw)


def is_dnazoo_source(value):
    return bool(re.search(r"dna\s*zoo|dnazoo", str(value or ""), re.I))


def materialize_dnazoo_wg(row, asm, fallback_species):
    token = dna_zoo_species_token(row.get("final_wg_ref_species") or fallback_species)
    stable_id = asm or f"DNAZOO_{token}"
    wg_dir = os.path.join("references", "wg", stable_id)
    gz_path = os.path.join(wg_dir, f"{stable_id}.genome.fa.gz")
    fasta_path = os.path.join(wg_dir, f"{stable_id}.genome.fa")
    url = f"https://dnazoo.s3.wasabisys.com/{token}/{token}.fasta.gz"
    print(f"[DNA Zoo] attempting WG URL: {url}", file=sys.stderr)
    print(f"[DNA Zoo] WG output: {fasta_path}", file=sys.stderr)
    download_file(url, gz_path)
    gunzip_to_fasta(gz_path, fasta_path)
    run([os.environ.get("SAMTOOLS_COMMAND", "samtools"), "faidx", fasta_path])
    return stable_id, fasta_path, "", url


def materialize_dnazoo_chrM(row, fallback_species):
    token = dna_zoo_species_token(row.get("final_chrM_species") or fallback_species)
    chrout = os.path.join("references", "chrM", "dnazoo", f"{token}_MT.fasta")
    url = f"https://dnazoo.s3.wasabisys.com/{token}/{token}_MT.fasta"
    print(f"[DNA Zoo] attempting chrM URL: {url}", file=sys.stderr)
    print(f"[DNA Zoo] chrM output: {chrout}", file=sys.stderr)
    download_file(url, chrout)
    run([os.environ.get("SAMTOOLS_COMMAND", "samtools"), "faidx", chrout])
    return chrout, url


def ensure_indexable_wg_fasta(wg_path, src_url):
    # NCBI genomic.fna.gz files are standard gzip, not bgzip. samtools faidx cannot
    # index ordinary gzip, so store the working WG FASTA decompressed as .fa.
    if wg_path.endswith(".gz"):
        fasta_path = wg_path[:-3]
        gz_path = wg_path
    else:
        fasta_path = wg_path
        gz_path = wg_path + ".gz"
    if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
        try:
            run([os.environ.get("SAMTOOLS_COMMAND", "samtools"), "faidx", fasta_path])
            return fasta_path
        except Exception:
            os.remove(fasta_path)
    download_file(src_url, gz_path)
    gunzip_to_fasta(gz_path, fasta_path)
    run([os.environ.get("SAMTOOLS_COMMAND", "samtools"), "faidx", fasta_path])
    return fasta_path


def clean_path(path):
    if not path:
        return
    for candidate in [path, path + ".fai", path + ".gz"]:
        try:
            if os.path.isdir(candidate):
                shutil.rmtree(candidate)
            elif os.path.exists(candidate):
                os.remove(candidate)
        except FileNotFoundError:
            pass


def clean_wg_materialization(asm):
    if asm:
        clean_path(os.path.join("references", "wg", asm))


def ncbi_summary_cache_path(kind):
    return os.path.join(os.path.dirname(dl), f"assembly_summary_{kind}.txt")


def ensure_ncbi_summary(kind):
    path = ncbi_summary_cache_path(kind)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/{kind}/assembly_summary_{kind}.txt"
    download_file(url, path)
    return path


def find_ncbi_assembly_row(accession):
    if not accession:
        return None
    fields = None
    for kind in ["refseq", "genbank"]:
        path = ensure_ncbi_summary(kind)
        with open(path, newline="") as handle:
            for line in handle:
                if line.startswith("# assembly_accession"):
                    fields = line[2:].rstrip("\n").split("\t")
                    continue
                if line.startswith("#"):
                    continue
                vals = line.rstrip("\n").split("\t")
                if fields and len(vals) >= len(fields) and vals[0] == accession:
                    return dict(zip(fields, vals))
    return None


def find_gca_gcf_partner(accession):
    row = find_ncbi_assembly_row(accession)
    if not row:
        return None
    partner = row.get("gbrs_paired_asm", "")
    if not partner or partner in ("na", "-", accession):
        return None
    return find_ncbi_assembly_row(partner)


def parse_assembly_report_chrM_candidates(report_path):
    candidates = []
    if not report_path or not os.path.exists(report_path):
        return candidates
    with open(report_path, newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sequence_role = fields[1].lower()
            assigned_molecule = fields[2].lower()
            genbank_accn = fields[4]
            refseq_accn = fields[6]
            ucsc_name = fields[9]
            if (
                assigned_molecule in ("mt", "m", "mitochondrion", "mitochondria", "chrm")
                or "mitochond" in assigned_molecule
                or "mitochond" in sequence_role
                or ucsc_name == "chrM"
            ):
                candidates.extend([ucsc_name, refseq_accn, genbank_accn, fields[0]])
    return [c for c in dict.fromkeys(candidates) if c and c not in ("na", "-")]



def swap_gca_gcf_accession(accession):
    if accession.startswith("GCF_"):
        return "GCA_" + accession[4:]
    if accession.startswith("GCA_"):
        return "GCF_" + accession[4:]
    return ""


def find_swapped_gca_gcf_row(accession):
    swapped = swap_gca_gcf_accession(accession)
    if not swapped:
        return None
    return find_ncbi_assembly_row(swapped)


def fasta_index_names(wg_path):
    fai = wg_path + ".fai" if wg_path else ""
    if not fai or not os.path.exists(fai):
        return set()
    names = set()
    with open(fai, newline="") as handle:
        for line in handle:
            if line.strip():
                names.add(line.split("\t", 1)[0])
    return names


def candidate_presence(candidates, wg_path):
    names = fasta_index_names(wg_path)
    found = []
    missing = []
    for cand in [c for c in dict.fromkeys(candidates) if c]:
        cand_found = cand in names or any(n.startswith(cand + ".") for n in names)
        (found if cand_found else missing).append(cand)
    return found, missing


def append_candidate_check(target, attempted_assembly, attempt_type, wg_path, report_path, candidates, status, fallback_action):
    found, missing = candidate_presence(candidates, wg_path)
    with open(cc, "a", newline="") as handle:
        handle.write("\t".join([
            target,
            attempted_assembly,
            attempt_type,
            wg_path or "",
            report_path or "",
            ",".join(candidates),
            ",".join(found),
            ",".join(missing),
            status,
            fallback_action,
        ]) + "\n")
    return found, missing


def append_manual_reason(row, reason):
    current = row.get("manual_review_reason", "")
    parts = [p for p in current.split(";") if p] + [p for p in reason.split(";") if p]
    row["manual_review_reason"] = ";".join(dict.fromkeys(parts))


def manifest_assembly_accessions(row):
    vals = [row.get(k, "") for k in [
        "chrM_source_assembly_accession",
        "final_chrM_assembly_accession",
        "final_chrM_accession",
    ]]
    return [v for v in dict.fromkeys(vals) if v.startswith(("GCA_", "GCF_"))]


def independent_chrM_accessions(row):
    vals = [row.get(k, "") for k in [
        "final_chrM_accession",
        "final_chrM_refseq_accn",
        "final_chrM_genbank_accn",
        "chrM_source_accession",
    ]]
    return [v for v in dict.fromkeys(vals) if v and v not in ("na", "-")]


def try_embedded_assembly_chrM(row, target, accession, attempt_type, candidate_sets):
    rowinfo = find_ncbi_assembly_row(accession)
    if not rowinfo:
        return None, f"{attempt_type}_{accession}_not_found_in_ncbi_summary"
    asm2 = rowinfo.get("assembly_accession", accession)
    chrout2 = os.path.join("references", "chrM", "embedded_from_wg", f"{asm2}.chrM.fa")
    try:
        asm2, wg2, report2 = materialize_wg_from_row(rowinfo)
        report_cands = parse_assembly_report_chrM_candidates(report2)
        all_cands = []
        for cset in [report_cands] + candidate_sets:
            all_cands.extend(cset)
        run_with_output(["bash", "preprocessing/scripts/extract_chrM_from_wg.sh", wg2, chrout2] + all_cands)
        append_candidate_check(target, asm2, attempt_type, wg2, report2, all_cands, "success", "extracted_embedded_chrM")
        return (asm2, wg2, report2, chrout2), "success"
    except Exception as err:
        clean_path(chrout2)
        append_candidate_check(target, asm2, attempt_type, locals().get("wg2", ""), locals().get("report2", ""), locals().get("all_cands", []), "failure", str(err).replace("\t", " "))
        return None, str(err).replace("\t", " ")


def try_independent_chrM(row):
    errors = []
    for acc in independent_chrM_accessions(row):
        chrout = os.path.join("references", "chrM", "independent", f"{acc}.fa")
        try:
            subprocess.check_call(["bash", "preprocessing/scripts/download_independent_chrM.sh", acc, chrout, os.environ.get("LOCAL_MITO_FASTA", "")])
            return acc, chrout, "success"
        except Exception as e:
            clean_path(chrout)
            errors.append(f"{acc}:{str(e).replace(chr(9), ' ')}")
    return "", "", "; ".join(errors) if errors else "no_independent_chrM_accession"

def materialize_wg_from_row(row):
    asm = row.get("assembly_accession", "")
    ftp = row.get("ftp_path", "")
    if not asm or not ftp or ftp in ("na", "-"):
        raise RuntimeError("paired_assembly_missing_ftp_or_accession")
    wg_path = f"references/wg/{asm}/{asm}.genome.fa"
    report_path = f"references/wg/{asm}/{asm}.assembly_report.txt"
    os.makedirs(os.path.dirname(wg_path), exist_ok=True)
    base = ftp.rstrip("/").split("/")[-1]
    src = f"{ftp.rstrip('/')}/{base}_genomic.fna.gz"
    rep = f"{ftp.rstrip('/')}/{base}_assembly_report.txt"
    wg_path = ensure_indexable_wg_fasta(wg_path, src)
    download_file(rep, report_path)
    return asm, wg_path, report_path


rows = list(csv.DictReader(open(man), delimiter="\t"))
for r in rows:
    asm = r.get("final_wg_assembly_accession", "")
    ftp = r.get("final_wg_ftp_path", "")
    target = r.get("target_species", "")
    manifest_wg = r.get("wg_expected_output_fasta", "")
    report = f"references/wg/{asm}/{asm}.assembly_report.txt" if asm else ""
    wg = manifest_wg[:-3] if manifest_wg.endswith(".gz") else manifest_wg
    status = "skipped"
    msg = "missing_ftp_or_assembly"
    if is_dnazoo_source(r.get("final_wg_ref_source", "")):
        try:
            asm, wg, report, attempted_url = materialize_dnazoo_wg(r, asm, target)
            if not r.get("final_wg_assembly_accession", ""):
                r["final_wg_assembly_accession"] = asm
            r["wg_expected_output_fasta"] = wg
            status = "success"
            msg = "downloaded_dnazoo_decompressed_and_indexed"
        except Exception as e:
            status = "failure"
            msg = f"dnazoo_wg_download_failed:{str(e).replace(chr(9), ' ')}"
    elif asm and ftp:
        os.makedirs(os.path.dirname(wg), exist_ok=True)
        base = ftp.rstrip("/").split("/")[-1]
        src = f"{ftp.rstrip('/')}/{base}_genomic.fna.gz"
        rep = f"{ftp.rstrip('/')}/{base}_assembly_report.txt"
        try:
            wg = ensure_indexable_wg_fasta(wg, src)
            download_file(rep, report)
            status = "success"
            msg = "downloaded_decompressed_and_indexed"
        except Exception as e:
            status = "failure"
            msg = str(e).replace("\t", " ")
    chrout = r.get("chrM_expected_output_fasta", "")
    ctx = r.get("chrM_reference_context", "")
    estatus = "skipped"
    emsg = "missing_chrM_ref"
    if is_dnazoo_source(r.get("final_chrM_source", "")) or is_dnazoo_source(r.get("final_wg_ref_source", "")) and not r.get("final_chrM_accession", "").startswith(("NC_", "CM_", "J")):
        try:
            chrout, attempted_chrM_url = materialize_dnazoo_chrM(r, target)
            r["chrM_reference_context"] = "independent_chrM_ref"
            r["chrM_extraction_strategy"] = "download_dnazoo_mt_fasta"
            r["chrM_expected_output_fasta"] = chrout
            ctx = r["chrM_reference_context"]
            estatus = "success"
            emsg = "downloaded_dnazoo_mt_fasta"
        except Exception as e:
            estatus = "failure"
            emsg = f"dnazoo_chrM_download_failed:{str(e).replace(chr(9), ' ')}"
    elif chrout and ctx == "embedded_in_wg_ref" and wg:
        cands = [r.get(k, "") for k in ["final_chrM_contig_name", "final_chrM_refseq_accn", "final_chrM_genbank_accn", "final_chrM_accession", "final_chrM_ucsc_name"]]
        report_cands = parse_assembly_report_chrM_candidates(report)
        current_cands = cands + report_cands
        try:
            run_with_output(["bash", "preprocessing/scripts/extract_chrM_from_wg.sh", wg, chrout] + current_cands)
            estatus = "success"
            emsg = "extracted_from_wg"
            append_candidate_check(target, asm, "current_final_wg", wg, report, current_cands, "success", "none")
        except Exception as e:
            first_error = str(e).replace("\t", " ")
            found, missing = append_candidate_check(target, asm, "current_final_wg", wg, report, report_cands, "failure", "evaluate_report_mismatch")
            report_mismatch = bool(report_cands) and not found
            fallback_notes = []
            recovered = False
            if report_mismatch:
                fallback_notes.append("wg_fasta_report_chrM_mismatch")
                paired_row = None
                paired_lookup_error = ""
                try:
                    paired_row = find_swapped_gca_gcf_row(asm)
                except Exception as err:
                    paired_lookup_error = str(err).replace("\t", " ")
                if paired_row:
                    paired_asm = paired_row.get("assembly_accession", "")
                    result, perr = try_embedded_assembly_chrM(
                        r, target, paired_asm, "paired_gca_gcf_after_report_mismatch", [cands, report_cands]
                    )
                    if result:
                        paired_asm, paired_wg, paired_report, paired_chrout = result
                        r["chrM_reference_context"] = "embedded_in_paired_wg_ref"
                        r["chrM_extraction_strategy"] = "extracted_from_paired_gca_gcf_after_report_mismatch"
                        r["chrM_source_assembly_accession"] = paired_asm
                        r["final_chrM_assembly_accession"] = paired_asm
                        r["chrM_expected_output_fasta"] = paired_chrout
                        ctx = r["chrM_reference_context"]
                        chrout = paired_chrout
                        estatus = "success"
                        emsg = f"wg_fasta_report_chrM_mismatch; recovered_with_paired_gca_gcf:{paired_asm}"
                        recovered = True
                    else:
                        fallback_notes.append(f"paired_gca_gcf_chrM_failed:{paired_asm}:{perr}")
                else:
                    fallback_notes.append("no_swapped_gca_gcf_partner_found" + (f":{paired_lookup_error}" if paired_lookup_error else ""))

            if not recovered:
                for src_asm in manifest_assembly_accessions(r):
                    if src_asm == asm or src_asm == swap_gca_gcf_accession(asm):
                        continue
                    result, merr = try_embedded_assembly_chrM(r, target, src_asm, "manifest_defined_chrM_assembly", [cands, report_cands])
                    if result:
                        src_asm, src_wg, src_report, src_chrout = result
                        r["chrM_reference_context"] = "embedded_in_wg_ref"
                        r["chrM_extraction_strategy"] = "extracted_from_manifest_chrM_assembly_after_report_mismatch" if report_mismatch else "extracted_from_manifest_chrM_assembly_after_initial_failure"
                        r["chrM_source_assembly_accession"] = src_asm
                        r["final_chrM_assembly_accession"] = src_asm
                        r["chrM_expected_output_fasta"] = src_chrout
                        chrout = src_chrout
                        estatus = "success"
                        emsg = "; ".join(fallback_notes + [f"recovered_with_manifest_assembly:{src_asm}"])
                        recovered = True
                        break
                    fallback_notes.append(f"manifest_assembly_chrM_failed:{src_asm}:{merr}")

            if not recovered:
                acc, independent_chrout, ierr = try_independent_chrM(r)
                if acc:
                    r["chrM_reference_context"] = "independent_chrM_ref"
                    r["chrM_extraction_strategy"] = "paired_gca_gcf_failed_use_independent_chrM"
                    r["chrM_source_accession"] = acc
                    r["chrM_expected_output_fasta"] = independent_chrout
                    ctx = r["chrM_reference_context"]
                    chrout = independent_chrout
                    estatus = "success"
                    emsg = "; ".join(fallback_notes + [f"used_independent_chrM:{acc}"])
                    if report_mismatch:
                        append_manual_reason(r, "wg_fasta_report_chrM_mismatch;paired_gca_gcf_chrM_failed;used_independent_chrM")
                else:
                    estatus = "failure"
                    emsg = "; ".join(fallback_notes + [first_error, ierr])
    elif chrout and ctx == "independent_chrM_ref":
        acc = r.get("final_chrM_accession", "")
        try:
            subprocess.check_call(["bash", "preprocessing/scripts/download_independent_chrM.sh", acc, chrout, os.environ.get("LOCAL_MITO_FASTA", "")])
            estatus = "success"
            emsg = "downloaded_or_extracted"
        except Exception as e:
            estatus = "failure"
            emsg = str(e).replace("\t", " ")
    with open(dl, "a") as handle:
        handle.write("\t".join([target, asm, status, wg, wg + ".fai" if wg else "", report, msg]) + "\n")
    with open(ex, "a") as handle:
        handle.write("\t".join([target, ctx, estatus, chrout, chrout + ".fai" if chrout else "", emsg]) + "\n")

if rows:
    fieldnames = list(rows[0].keys())
    for row in rows:
        row["reference_pairing_status"] = reference_pairing_status(
            row.get("target_species", ""),
            row.get("final_wg_ref_species", ""),
            row.get("final_chrM_species", ""),
        )
    for extra_col in ["manual_review_reason", "chrM_reference_context", "chrM_extraction_strategy", "chrM_source_assembly_accession", "chrM_source_accession", "final_chrM_assembly_accession", "chrM_expected_output_fasta", "wg_expected_output_fasta", "reference_pairing_status"]:
        if extra_col not in fieldnames:
            fieldnames.append(extra_col)
    for resolved_manifest in [resolved_results, resolved_refs]:
        os.makedirs(os.path.dirname(resolved_manifest), exist_ok=True)
        with open(resolved_manifest, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
PY
"$PYTHON_COMMAND" - "$RESOLVED_RESULTS" "$DL" "$EX" "$IH" <<'PY'
import csv, re, sys
from collections import Counter
man,dl,ex,ih=sys.argv[1:]
def normalize_species(value):
  return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")
def reference_pairing_status(target, wg_species, chrM_species):
  t=normalize_species(target); wg=normalize_species(wg_species); mt=normalize_species(chrM_species)
  has_wg=bool(wg); has_mt=bool(mt)
  if has_wg and has_mt and wg == t and mt == t: return "same_species_wg_same_species_chrM"
  if has_wg and has_mt and wg == t and mt != t: return "same_species_wg_cross_species_chrM"
  if has_wg and has_mt and wg != t and mt == t: return "cross_species_wg_same_species_chrM"
  if has_wg and has_mt and wg != t and mt != t: return "cross_species_wg_cross_species_chrM"
  if has_wg and not has_mt: return "wg_only_no_chrM"
  if not has_wg and has_mt: return "chrM_only_no_wg"
  return "no_reference_found"
dlmap={r['target_species']:r for r in csv.DictReader(open(dl), delimiter='\t')}
exmap={r['target_species']:r for r in csv.DictReader(open(ex), delimiter='\t')}
cols='target_species reference_pairing_status chrM_reference_context final_wg_ref_species final_wg_assembly_accession wg_fasta_path wg_fai_path wg_assembly_report_path final_chrM_species final_chrM_accession chrM_fasta_path chrM_fai_path chrM_extraction_strategy final_reference_strategy manual_review_required manual_review_reason'.split()
counts=Counter()
with open(ih,'w',newline='') as out_handle:
  w=csv.DictWriter(out_handle, fieldnames=cols, delimiter='\t'); w.writeheader()
  for r in csv.DictReader(open(man), delimiter='\t'):
    d=dlmap.get(r['target_species'],{}); e=exmap.get(r['target_species'],{})
    reasons=[x for x in [r.get('manual_review_reason','')] if x]
    if d.get('status') not in ('success','skipped'): reasons.append('wg_download_'+d.get('status','failure'))
    if e.get('status') not in ('success','skipped'): reasons.append('chrM_materialization_'+e.get('status','failure'))
    status=reference_pairing_status(r.get('target_species',''), r.get('final_wg_ref_species',''), r.get('final_chrM_species',''))
    out={c:r.get(c,'') for c in cols}; out.update({'reference_pairing_status':status,'wg_fasta_path':d.get('wg_fasta_path',r.get('wg_expected_output_fasta','')),'wg_fai_path':d.get('wg_fai_path',''),'wg_assembly_report_path':d.get('wg_assembly_report_path',''),'chrM_fasta_path':e.get('chrM_fasta_path',r.get('chrM_expected_output_fasta','')),'chrM_fai_path':e.get('chrM_fai_path',''),'manual_review_required':'yes' if reasons else 'no','manual_review_reason':';'.join(dict.fromkeys(reasons))})
    counts[status]+=1
    w.writerow(out)
print("in_house_score_reference_inputs.tsv reference_pairing_status counts:", file=sys.stderr)
for k,v in counts.most_common():
  print(f"  {k}\t{v}", file=sys.stderr)
non_missing=sum(v for k,v in counts.items() if k)
if non_missing and counts["cross_species_wg_cross_species_chrM"] / non_missing > 0.9:
  print("WARNING: >90% of non-missing rows are cross_species_wg_cross_species_chrM", file=sys.stderr)
PY
cp "$IH" references/manifests/in_house_score_reference_inputs.tsv
cmp -s "$IH" references/manifests/in_house_score_reference_inputs.tsv || { echo "ERROR: in-house score input copies differ" >&2; exit 1; }
# Post-run validation examples:
# md5sum results/preprocessing/reference_materialization/in_house_score_reference_inputs.tsv references/manifests/in_house_score_reference_inputs.tsv
# awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next} {count[$h["reference_pairing_status"]]++} END {for(k in count) print count[k], k}' references/manifests/in_house_score_reference_inputs.tsv | sort -nr
"$PYTHON_COMMAND" - "$RESOLVED_RESULTS" <<'PY'
import csv, sys
from collections import Counter
counts=Counter(r.get("reference_pairing_status","") for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"))
print("reference_materialization_manifest.resolved.tsv reference_pairing_status counts:", file=sys.stderr)
for k,v in counts.most_common():
    print(f"  {k}\t{v}", file=sys.stderr)
non_missing=sum(v for k,v in counts.items() if k)
if non_missing and counts["cross_species_wg_cross_species_chrM"] / non_missing > 0.9:
    print("WARNING: >90% of non-missing rows are cross_species_wg_cross_species_chrM", file=sys.stderr)
PY
cmp -s "$RESOLVED_RESULTS" "$RESOLVED_REFS" || { echo "ERROR: resolved manifest copies differ" >&2; exit 1; }
