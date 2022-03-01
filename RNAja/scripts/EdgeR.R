suppressMessages(library('knitr', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('edgeR', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('limma', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('RColorBrewer', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('mixOmics', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('baySeq', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('Cairo', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('plotly', warn.conflict = FALSE, quietly = TRUE))
suppressMessages(library('Rsubread', warn.conflict = FALSE, quietly = TRUE))

sample_info <- snakemake@input[["sample_info"]]
comparaisons_file <- snakemake@input[["de_comparisons_file"]]
gtf_file <- snakemake@input[["gtf_file"]]

outDir <- snakemake@params[["out_dir"]]
invisible(dir.exists(outDir) || dir.create(path = outDir, mode = "770"))

normCount_file <- snakemake@params[["normCount"]]
counts_file <- snakemake@params[["Counts"]]
#resDE_file <- snakemake@params[["resDE"]]
diane_counts_file <- snakemake@params[["dianeCounts"]]

# ---------------------------------------------------------------------------------

pheno_data <- read.csv(sample_info)
sampleInfo <- data.frame(fileName= pheno_data$SampleName, condition = pheno_data$Treatment, replicate= pheno_data$Experiment)
rawCountTable <- as.matrix(read.csv(gcsv_hisat,row.names="gene_id"))

# Output for Diane

df = as.data.frame(rawCountTable)
df <- data.frame(Gene = row.names(df), df)
write.table(df, file=diane_counts_file, sep=";", row.names = FALSE)

dgList <- DGEList(rawCountTable, group=sampleInfo$condition, genes=rownames(rawCountTable))

# filter
keep <- filterByExpr(dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# normalisation
dgList <- calcNormFactors(dgList, method="TMM")

eff.lib.size <- dgList$samples$lib.size*dgList$samples$norm.factors
normCounts <- cpm(dgList)

# output for diffexdb

write.table(dgList$counts,file=counts_file,sep=";")
write.table(normCounts, file=normCount_file, sep=";")

