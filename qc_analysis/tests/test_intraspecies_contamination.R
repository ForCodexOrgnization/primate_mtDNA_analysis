# Smoke/integration test for installations with R available. Run with:
# Rscript qc_analysis/tests/test_intraspecies_contamination.R
tmp <- tempfile("intraspecies-"); dir.create(tmp)
input <- file.path(tmp,"variants.tsv")
write.table(data.frame(Sample=c(rep("A",8),rep("B",8),"C"),Species=c(rep("sp",16),"singleton"),CHROM="chrM",POS=seq_len(17),REF="A",ALT="G",Type="SNV",FILTER="PASS",DP=100,VAF=c(rep(.1,5),rep(.9,3),rep(.99,8),.1),check.names=FALSE),input,sep="\t",row.names=FALSE,quote=FALSE)
status <- system2("Rscript",c("qc_analysis/scripts/run_intraspecies_contamination.R","--variant-table",input,"--outdir",file.path(tmp,"out")))
stopifnot(status==0)
x <- read.delim(file.path(tmp,"out/tables/final_contamination_summary.tsv"),check.names=FALSE)
stopifnot(nrow(x)==3, x$contamination_status[x$Sample=="C"]=="insufficient_singleton_species", all(x$mirror_calibration_status=="not_calibrated_no_file"))
