#!/usr/bin/env python3
"""Build a version-2 tRNA position index from tRNAscan-SE output."""
import argparse, csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.trnascan_utils import build_trna_position_index, run_trnascan

SUMMARY_COLUMNS = ["reference_key", "fasta", "fasta_length", "trnascan_out", "trnascan_ss",
                   "output_index", "n_trna_records", "n_index_rows", "n_positive_strand_trna",
                   "n_negative_strand_trna", "n_stem_positions", "n_loop_positions",
                   "n_missing_structure", "n_fasta_sequence_mismatch", "status", "notes"]

def write_summary(path, row):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer=csv.DictWriter(handle, SUMMARY_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerow(row)

def build(args):
    output=Path(args.output)
    if output.exists() and not args.overwrite: raise FileExistsError(f"Output exists (use --overwrite): {output}")
    out, ss = args.trnascan_out, args.trnascan_ss
    if args.run_trnascan:
        prefix=args.trnascan_prefix or str(output).removesuffix(".trna_position_index.tsv.gz").removesuffix(".tsv.gz")
        made=run_trnascan(args.fasta,prefix,args.trnascan_bin,args.trnascan_mode,args.threads,args.trnascan_extra_args)
        out,ss=str(made["out"]),str(made["ss"])
    if not out or not ss: raise ValueError("Provide --trnascan-out and --trnascan-ss, or --run-trnascan")
    result=build_trna_position_index(args.reference_key,args.fasta,out,ss,output,args.chrom_normalization,args.max_sequence_mismatch_rate)
    records,rows=result["records"],result["rows"]
    return {"reference_key":args.reference_key,"fasta":args.fasta,"fasta_length":result["fasta_length"],
            "trnascan_out":out,"trnascan_ss":ss,"output_index":str(output),"n_trna_records":len(records),
            "n_index_rows":len(rows),"n_positive_strand_trna":sum(r.strand=="+" for r in records),
            "n_negative_strand_trna":sum(r.strand=="-" for r in records),
            "n_stem_positions":sum(r["struct_class"]=="stem" for r in rows),
            "n_loop_positions":sum(r["struct_class"]=="loop" for r in rows),
            "n_missing_structure":sum(not r["struct_char"] for r in rows),
            "n_fasta_sequence_mismatch":result["n_fasta_sequence_mismatch"],"status":"completed","notes":""}

def parser():
    p=argparse.ArgumentParser(); p.add_argument("--reference-key",required=True); p.add_argument("--fasta",required=True)
    p.add_argument("--trnascan-out"); p.add_argument("--trnascan-ss"); p.add_argument("--run-trnascan",action="store_true")
    p.add_argument("--trnascan-bin",default="tRNAscan-SE"); p.add_argument("--trnascan-mode",default="mito_mammal")
    p.add_argument("--threads",type=int,default=1); p.add_argument("--trnascan-extra-args",default=""); p.add_argument("--trnascan-prefix")
    p.add_argument("--output",required=True); p.add_argument("--overwrite",action="store_true"); p.add_argument("--summary")
    p.add_argument("--chrom-normalization",default="none",choices=["none","strip_chr","add_chr","mitochondrial_alias"])
    p.add_argument("--max-sequence-mismatch-rate",type=float,default=0.0); return p

def main():
    args=parser().parse_args(); summary=args.summary or str(args.output)+".summary.tsv"
    try: row=build(args)
    except Exception as exc:
        row={"reference_key":args.reference_key,"fasta":args.fasta,"trnascan_out":args.trnascan_out or "",
             "trnascan_ss":args.trnascan_ss or "","output_index":args.output,"status":"failed","notes":str(exc)}
        write_summary(summary,row); raise
    write_summary(summary,row)
if __name__ == "__main__": main()
