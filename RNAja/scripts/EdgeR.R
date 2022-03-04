suppressMessages(library('knitr', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('edgeR', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('limma', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('RColorBrewer', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('mixOmics', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('baySeq', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('Cairo', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('plotly', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('Rsubread', warn.conflict = FALSE, quietly = TRUE))

# input file

gcsv_hisat <- snakemake@input[["gcsv_hisat"]]
sample_info <- snakemake@input[["sample_info"]]
# <- snakemake@input[["de_comparisons_file"]]
gtf_file <- snakemake@input[["gtf_file"]]

outDir <- snakemake@params[["out_dir"]]
invisible(dir.exists(outDir) || dir.create(path = outDir, mode = "770"))

#normCount_file <- snakemake@params[["normCount"]]
#counts_file <- snakemake@params[["Counts"]]
#resDE_file <- snakemake@params[["resDE"]]
diane_counts_file <- snakemake@output[["dianeCounts"]]

# ---------------------------------------------------------------------------------

pheno_data <- read.csv(sample_info)

sampleInfo <- data.frame(fileName= pheno_data$SampleName)
for (i in 3:ncol(pheno_data)){
  sampleInfo <- cbind(sampleInfo, Condition = pheno_data[,i])
}
N_groups = ""
for (i in 2:ncol(sampleInfo)){
  N_groups = paste0(N_groups,sampleInfo[,i])
}
N = length(unique(N_groups))

rawCountTable <- as.matrix(read.csv(gcsv_hisat,row.names="gene_id"))

dgList <- DGEList(rawCountTable, group=rep(1:N,each=3), genes=rownames(rawCountTable))

#pheno_data <- read.csv(sample_info)
#sampleInfo <- data.frame(fileName= pheno_data$SampleName, condition = pheno_data$Condition1, replicate= pheno_data$Condition2)
#rawCountTable <- as.matrix(read.csv(gcsv_hisat,row.names="gene_id"))
#dgList <- DGEList(rawCountTable, group=sampleInfo$condition, genes=rownames(rawCountTable))

# Output for Diane

df = as.data.frame(rawCountTable)
df <- data.frame(Gene = row.names(df), df)
write.table(df, file=diane_counts_file, sep=";", row.names = FALSE)



# filter
keep <- filterByExpr(dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# normalisation
dgList <- calcNormFactors(dgList, method="TMM")

eff.lib.size <- dgList$samples$lib.size*dgList$samples$norm.factors
normCounts <- cpm(dgList)

# output for diffexdb

#write.table(dgList$counts,file=counts_file,sep=";")
#write.table(normCounts, file=normCount_file, sep=";")

