#!/usr/bin/env python3
"""Build final variant-calling reference FASTA packages from preprocessing outputs."""
from __future__ import annotations
import argparse, csv, os, re, shutil, subprocess, sys
from pathlib import Path

MANIFEST_COLUMNS = [
    "target_species","safe_species_id","REF_TYPE","MTLIKE_PATTERN","ValidAnnotatedMitoContig","HasValidAnnotatedMito",
    "MaskPriority","numt_mask_bed","numt_mask_applied_to_whole_ref","whole_fasta","whole_fai","whole_dict",
    "chrM_fasta","chrM_fai","chrM_dict","chrM_shift_fasta","chrM_shift_fai","chrM_shift_dict",
    "non_control_interval","control_region_shifted_interval","shift_back_chain","wg_fasta_source","chrM_fasta_source",
    "chrM_reference_context","reference_pairing_status","final_reference_strategy","build_status","build_message",
]

def safe_id(value: str) -> str:
    x = re.sub(r"[:/\s]+", "_", (value or "unknown").strip())
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x).strip("_")
    return x or "unknown"

def read_tsv(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def fasta_iter(path: Path):
    name = None; seq = []
    with path.open() as h:
        for line in h:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name=line[1:].split()[0]; seq=[]
            else:
                seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

def write_record(path: Path, name: str, seq: str):
    with path.open("w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80): out.write(seq[i:i+80] + "\n")

def select_chrm_record(chrm_fa: Path):
    recs = list(fasta_iter(chrm_fa))
    if len(recs) == 1: return recs[0][1]
    candidates = [(n,s) for n,s in recs if re.search(r"mitochond|chrM|MT\b|mtdna", n, re.I)]
    if len(candidates) == 1: return candidates[0][1]
    raise ValueError(f"{chrm_fa} has {len(recs)} records and chrM cannot be selected unambiguously")

def load_bed(path: Path):
    intervals = []
    if not path or not path.exists() or path.stat().st_size == 0: return intervals
    with path.open() as h:
        for line in h:
            if not line.strip() or line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f) >= 3: intervals.append((f[0], int(f[1]), int(f[2])))
    return intervals

def mask_seq(seq: str, ivals):
    arr=list(seq)
    for a,b in ivals:
        a=max(0,a); b=min(len(arr),b)
        if a < b: arr[a:b] = "N"*(b-a)
    return "".join(arr)

def build_whole(wg: Path, out: Path, chrm_seq: str, valid_contig: str, bed: Path|None):
    by_contig = {}
    for c,a,b in load_bed(bed): by_contig.setdefault(c,[]).append((a,b))
    chrm_count = 0
    with out.open("w") as fh:
        for name, seq in fasta_iter(wg):
            out_name = "chrM" if valid_contig and valid_contig != "NA" and name == valid_contig else name
            if out_name == "chrM": chrm_count += 1
            if name in by_contig and out_name != "chrM": seq = mask_seq(seq, by_contig[name])
            fh.write(f">{out_name}\n")
            for i in range(0, len(seq), 80): fh.write(seq[i:i+80]+"\n")
        if chrm_count == 0:
            chrm_count = 1
            fh.write(">chrM\n")
            for i in range(0, len(chrm_seq), 80): fh.write(chrm_seq[i:i+80]+"\n")
    if chrm_count != 1:
        raise ValueError(f"final whole reference {out} has {chrm_count} chrM records; expected exactly one")

def write_intervals(sid: str, chrm_len: int, interval_dir: Path):
    non = interval_dir / f"{sid}_non_control_region.interval_list"
    ctrl = interval_dir / f"{sid}_control_region_shifted.interval_list"
    end = max(1, chrm_len)
    non.write_text(f"@SQ\tSN:chrM\tLN:{end}\nchrM\t1\t{end}\t+\tnon_control_region\n")
    ctrl.write_text(f"@SQ\tSN:chrM\tLN:{end}\nchrM\t1\t{min(8000,end)}\t+\tcontrol_region_shifted\n")
    return non, ctrl

def write_chain(sid: str, chrm_len: int, path: Path):
    path.write_text(f"chain 1 chrM {chrm_len} + 0 {chrm_len} chrM {chrm_len} + 0 {chrm_len} 1\n{chrm_len}\n\n")

def run_index(fa: Path, samtools: str, bwa: str, gatk: str):
    dict_path = fa.with_suffix(".dict")
    subprocess.run([samtools,"faidx",str(fa)], check=True)
    subprocess.run([bwa,"index",str(fa)], check=True)
    subprocess.run([gatk,"CreateSequenceDictionary","-R",str(fa),"-O",str(dict_path)], check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ref-inputs", default="references/manifests/in_house_score_reference_inputs.tsv")
    ap.add_argument("--score-file", default="results/preprocessing/in_house_score/merged_in_house_score.tsv")
    ap.add_argument("--outdir", default="references/variant_calling")
    ap.add_argument("--mask-ref-types", default=os.environ.get("MASK_REF_TYPES", "#C-likely_comp,#C-Ambiguous"))
    ap.add_argument("--skip-index", action="store_true")
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS_COMMAND","samtools"))
    ap.add_argument("--bwa", default=os.environ.get("BWA_COMMAND","bwa"))
    ap.add_argument("--gatk", default=os.environ.get("GATK_COMMAND","gatk"))
    args=ap.parse_args()
    outdir=Path(args.outdir); dirs={n:outdir/n for n in ["Ref_whole","Ref_chrM","Ref_chrM_shift","interval","shift_back_chain"]}
    for d in dirs.values(): d.mkdir(parents=True, exist_ok=True)
    scores = {r.get("Species") or r.get("target_species"): r for r in read_tsv(Path(args.score_file))}
    mask_types={x.strip() for x in args.mask_ref_types.split(",") if x.strip()}
    rows=[]; seen=set()
    for r in read_tsv(Path(args.ref_inputs)):
        key=(r.get("target_species",""), r.get("wg_fasta_path",""), r.get("chrM_fasta_path",""))
        if key in seen: continue
        seen.add(key)
        sp=key[0]; sid=safe_id(sp); score=scores.get(sp,{})
        base={c:"" for c in MANIFEST_COLUMNS}; base.update({"target_species":sp,"safe_species_id":sid,"wg_fasta_source":key[1],"chrM_fasta_source":key[2],
            "REF_TYPE":score.get("REF_TYPE",""),"MTLIKE_PATTERN":score.get("MTLIKE_PATTERN",""),"ValidAnnotatedMitoContig":score.get("ValidAnnotatedMitoContig",""),
            "HasValidAnnotatedMito":score.get("HasValidAnnotatedMito",""),"MaskPriority":score.get("MaskPriority",""),
            "chrM_reference_context":r.get("chrM_reference_context",""),"reference_pairing_status":r.get("reference_pairing_status",""),"final_reference_strategy":r.get("final_reference_strategy","")})
        try:
            wg=Path(key[1]); chrm=Path(key[2])
            if not wg.exists() or not chrm.exists(): raise FileNotFoundError("missing WG or chrM FASTA source")
            chrm_seq=select_chrm_record(chrm)
            chrm_fa=dirs["Ref_chrM"]/f"{sid}.chrM.fa"; write_record(chrm_fa,"chrM",chrm_seq)
            shift_seq=chrm_seq[8000:]+chrm_seq[:8000] if len(chrm_seq)>8000 else chrm_seq
            shift_fa=dirs["Ref_chrM_shift"]/f"{sid}.chrM_shift.fa"; write_record(shift_fa,"chrM",shift_seq)
            mask_bed = None
            if score.get("MinimalMaskBED"):
                minimal_bed = Path(score["MinimalMaskBED"])
                mask_bed = minimal_bed.with_name(minimal_bed.name.replace(".minimal_chrMcover.bed", ".FINAL_numt_mask.bed"))
            # Prefer sibling FINAL bed created by in_house_score; do not mask chrM or chrM_shift.
            apply_mask=base["REF_TYPE"] in mask_types and mask_bed and mask_bed.exists() and sum(1 for l in mask_bed.open() if l.strip() and not l.startswith("#"))>0
            whole=dirs["Ref_whole"]/f"{sid}.whole.fa"; build_whole(wg, whole, chrm_seq, base["ValidAnnotatedMitoContig"], mask_bed if apply_mask else None)
            non, ctrl=write_intervals(sid,len(chrm_seq),dirs["interval"]); chain=dirs["shift_back_chain"]/f"{sid}_ShiftBack.chain"; write_chain(sid,len(chrm_seq),chain)
            if not args.skip_index:
                for fa in (whole,chrm_fa,shift_fa): run_index(fa,args.samtools,args.bwa,args.gatk)
            base.update({"numt_mask_bed":str(mask_bed or ""),"numt_mask_applied_to_whole_ref":"yes" if apply_mask else "no","whole_fasta":str(whole),"whole_fai":str(whole)+".fai","whole_dict":str(whole.with_suffix(".dict")),"chrM_fasta":str(chrm_fa),"chrM_fai":str(chrm_fa)+".fai","chrM_dict":str(chrm_fa.with_suffix(".dict")),"chrM_shift_fasta":str(shift_fa),"chrM_shift_fai":str(shift_fa)+".fai","chrM_shift_dict":str(shift_fa.with_suffix(".dict")),"non_control_interval":str(non),"control_region_shifted_interval":str(ctrl),"shift_back_chain":str(chain),"build_status":"success","build_message":"ok"})
        except Exception as e:
            base.update({"build_status":"failed","build_message":str(e)})
        rows.append(base)
    manifest=outdir/"variant_calling_reference_manifest.tsv"
    with manifest.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=MANIFEST_COLUMNS,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"Wrote {manifest}")
if __name__ == "__main__": main()
