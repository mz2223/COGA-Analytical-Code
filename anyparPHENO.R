library(GWAF) 
library(glue)
library(parallel)

num_cores <- 15
cat("num_cores", num_cores, "\n")
cl <- makeCluster(num_cores)

args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

args_full <- commandArgs(trailingOnly = FALSE)
print(args_full)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_full[grep(file_arg, args_full)])
print(script_path)
PHENO_DIR <- dirname(normalizePath(script_path))
print(PHENO_DIR)
PHENO_SEX <- basename(PHENO_DIR)
print(PHENO_SEX)
parts <- strsplit(PHENO_SEX, "_")[[1]]
SEX <- tail(parts, 1)
PHENO <- paste(head(parts, -1), collapse = "_")
#PHENO <- parts[1]
#SEX <- parts[2]
print(PHENO)
print(SEX)

cat("chr", i, "\n")
chrdir <- glue("/N/project/mmge_audprs/coga/ea_fams/plink_files/chr{i}")
csv_files <- list.files(path = chrdir, pattern = "\\.raw.csv$", full.names = TRUE)

process_file <- function(csvfile) {
    print(csvfile)
    number <- sub(".*/([0-9]+)\\.raw\\.csv$", "\\1", csvfile)
    print(number)

    geepack.lgst.batch(phenfile=glue("/N/project/mmge_audprs/coga/ea_fams/analyses/withdrawal/EA_withdrawal_symptom_{SEX}.csv"),
            pedfile="/N/project/mmge_audprs/coga/ea_fams/analyses/pedfile_ea_gwaf.csv",
            genfile=glue("{csvfile}"),
            # genfile <- glue("/N/project/mmge_audprs/coga/ea_fams/plink_files/chr{i}/{csvfile}"),
            phen=glue("{PHENO}"),
            model="a",
            covars=c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
            outfile=glue("/N/project/mmge_audprs/coga/ea_fams/analyses/{PHENO_SEX}/chr{i}/{number}.{PHENO}.out"),
            col.names=T,
            sep.ped=,
            sep.phe=,
            sep.gen=)
}

clusterExport(cl, list("process_file", "csv_files", "i","PHENO", "PHENO_SEX", "SEX"))

clusterEvalQ(cl, library(geepack))
clusterEvalQ(cl, library(GWAF))
clusterEvalQ(cl, library(glue))
clusterEvalQ(cl, library(coxme))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(bdsmatrix))
clusterEvalQ(cl, library(lme4))
clusterEvalQ(cl, library(Matrix))

results <- parLapply(cl, csv_files, process_file)

stopCluster(cl)

