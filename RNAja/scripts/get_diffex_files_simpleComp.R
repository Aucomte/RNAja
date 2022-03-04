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
gcsv_hisat <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_sev/3_count/STRINGTIE/HISAT_gene_count_matrix.csv"
sample_info <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/DATA2/sample_info.txt"

normCount_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/test/normCounts_sev"
counts_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/test/Counts_sev"
resDE_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/test/resde"

setwd("/home/comtea/Documents/PHIM/BURKADAPT/RNAja/test/")

pheno_data <- read.csv(sample_info)

sampleInfo <- data.frame(fileName= pheno_data$SampleName, Condition = pheno_data[,3])
rawCountTable <- as.matrix(read.csv(gcsv_hisat,row.names="gene_id"))
dgList <- DGEList(rawCountTable, group=sampleInfo$Condition, genes=rownames(rawCountTable))

# filter
keep <- filterByExpr(dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# normalisation
dgList <- calcNormFactors(dgList, method="TMM")

eff.lib.size <- dgList$samples$lib.size*dgList$samples$norm.factors
normCounts <- cpm(dgList)

# output for diffexdb

write.table(dgList$counts,file="counts_file",sep=";")
write.table(normCounts, file="normCount_file", sep=";")

# design

Condition = factor(sampleInfo[,2])

design = model.matrix(~Condition)

dgList <- estimateDisp(dgList, design, robust=TRUE)
fit <- glmQLFit(dgList, design)

## -------------

qlf_MOCKvsRYMV <- glmQLFTest(fit)
tab = as.data.frame(qlf_MOCKvsRYMV$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "MOCK_vs_RYMV.txt",sep="\t", quote = FALSE, row.names = FALSE)

