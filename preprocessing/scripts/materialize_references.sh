#!/usr/bin/env bash
set -euo pipefail
MANIFEST=${1:-references/manifests/reference_materialization_manifest.tsv}
OUTDIR=${OUTDIR:-results/preprocessing/reference_materialization}
LOCAL_MITO_FASTA=${LOCAL_MITO_FASTA:-}
mkdir -p "$OUTDIR" references/wg references/chrM/embedded_from_wg references/chrM/independent references/manifests
DL="$OUTDIR/reference_download_manifest.tsv"; EX="$OUTDIR/chrM_extraction_manifest.tsv"; IH="$OUTDIR/in_house_score_reference_inputs.tsv"
PYTHON_COMMAND=${PYTHON_COMMAND:-python3}
WGET_COMMAND=${WGET_COMMAND:-wget}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
CURL_COMMAND=${CURL_COMMAND:-curl}
EFETCH_COMMAND=${EFETCH_COMMAND:-efetch}
export WGET_COMMAND SAMTOOLS_COMMAND CURL_COMMAND EFETCH_COMMAND
echo -e "target_species\tassembly_accession\tstatus\twg_fasta_path\twg_fai_path\twg_assembly_report_path\tmessage" > "$DL"
echo -e "target_species\tchrM_reference_context\tstatus\tchrM_fasta_path\tchrM_fai_path\tmessage" > "$EX"
"$PYTHON_COMMAND" - "$MANIFEST" "$DL" "$EX" <<'PY'
import csv, os, subprocess, sys
man, dl, ex = sys.argv[1:]
rows=list(csv.DictReader(open(man), delimiter='\t'))
for r in rows:
  asm=r.get('final_wg_assembly_accession',''); ftp=r.get('final_wg_ftp_path',''); target=r.get('target_species','')
  wg=r.get('wg_expected_output_fasta',''); report=f"references/wg/{asm}/{asm}.assembly_report.txt" if asm else ''
  status='skipped'; msg='missing_ftp_or_assembly'
  if asm and ftp and 'dnazoo' not in r.get('final_wg_ref_source','').lower():
    os.makedirs(os.path.dirname(wg), exist_ok=True); base=ftp.rstrip('/').split('/')[-1]
    src=f"{ftp.rstrip('/')}/{base}_genomic.fna.gz"; rep=f"{ftp.rstrip('/')}/{base}_assembly_report.txt"
    try:
      subprocess.check_call([os.environ.get('WGET_COMMAND','wget'),'-c','-O',wg,src]); subprocess.call([os.environ.get('WGET_COMMAND','wget'),'-c','-O',report,rep]); subprocess.check_call([os.environ.get('SAMTOOLS_COMMAND','samtools'),'faidx',wg]); status='success'; msg='downloaded'
    except Exception as e: status='failure'; msg=str(e).replace('\t',' ')
  elif 'dnazoo' in r.get('final_wg_ref_source','').lower(): status='manual_review'; msg='dnazoo_download_not_implemented'
  open(dl,'a').write('\t'.join([target,asm,status,wg,wg+'.fai' if wg else '',report,msg])+'\n')
  chrout=r.get('chrM_expected_output_fasta',''); ctx=r.get('chrM_reference_context',''); estatus='skipped'; emsg='missing_chrM_ref'
  if chrout and ctx=='embedded_in_wg_ref' and wg:
    cands=[r.get(k,'') for k in ['final_chrM_contig_name','final_chrM_refseq_accn','final_chrM_genbank_accn','final_chrM_accession','final_chrM_ucsc_name']]
    try: subprocess.check_call(['bash','preprocessing/scripts/extract_chrM_from_wg.sh',wg,chrout]+cands); estatus='success'; emsg='extracted_from_wg'
    except Exception as e: estatus='failure'; emsg=str(e).replace('\t',' ')
  elif chrout and ctx=='independent_chrM_ref':
    acc=r.get('final_chrM_accession','')
    try: subprocess.check_call(['bash','preprocessing/scripts/download_independent_chrM.sh',acc,chrout,os.environ.get('LOCAL_MITO_FASTA','')]); estatus='success'; emsg='downloaded_or_extracted'
    except Exception as e: estatus='failure'; emsg=str(e).replace('\t',' ')
  open(ex,'a').write('\t'.join([target,ctx,estatus,chrout,chrout+'.fai' if chrout else '',emsg])+'\n')
PY
"$PYTHON_COMMAND" - "$MANIFEST" "$DL" "$EX" "$IH" <<'PY'
import csv, sys
man,dl,ex,ih=sys.argv[1:]
dlmap={r['target_species']:r for r in csv.DictReader(open(dl), delimiter='\t')}
exmap={r['target_species']:r for r in csv.DictReader(open(ex), delimiter='\t')}
cols='target_species reference_pairing_status chrM_reference_context final_wg_ref_species final_wg_assembly_accession wg_fasta_path wg_fai_path wg_assembly_report_path final_chrM_species final_chrM_accession chrM_fasta_path chrM_fai_path chrM_extraction_strategy final_reference_strategy manual_review_required manual_review_reason'.split()
w=csv.DictWriter(open(ih,'w',newline=''), fieldnames=cols, delimiter='\t'); w.writeheader()
for r in csv.DictReader(open(man), delimiter='\t'):
  d=dlmap.get(r['target_species'],{}); e=exmap.get(r['target_species'],{})
  reasons=[x for x in [r.get('manual_review_reason','')] if x]
  if d.get('status') not in ('success','skipped'): reasons.append('wg_download_'+d.get('status','failure'))
  if e.get('status') not in ('success','skipped'): reasons.append('chrM_materialization_'+e.get('status','failure'))
  out={c:r.get(c,'') for c in cols}; out.update({'wg_fasta_path':d.get('wg_fasta_path',r.get('wg_expected_output_fasta','')),'wg_fai_path':d.get('wg_fai_path',''),'wg_assembly_report_path':d.get('wg_assembly_report_path',''),'chrM_fasta_path':e.get('chrM_fasta_path',r.get('chrM_expected_output_fasta','')),'chrM_fai_path':e.get('chrM_fai_path',''),'manual_review_required':'yes' if reasons else 'no','manual_review_reason':';'.join(dict.fromkeys(reasons))})
  w.writerow(out)
PY
cp "$IH" references/manifests/in_house_score_reference_inputs.tsv
