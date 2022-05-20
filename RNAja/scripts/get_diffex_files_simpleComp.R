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
gcsv_hisat <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all_R/3_count/STRINGTIE/HISAT_gene_count_matrix.csv"
sample_info <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/DATA/sample_info_all_R.txt"

normCount_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all_R/4_DE_analysis/normcount.csv"
counts_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all_R/4_DE_analysis/Counts.csv"

setwd("/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all_R/4_DE_analysis/")

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

write.table(dgList$counts,file="counts_file",sep=";",quote=FALSE)
write.table(normCounts, file="normCount_file", sep=";",quote=FALSE)

# dispersion

dgList <- estimateDisp(dgList)

## -------------

w = exactTest(dgList, pair=c("Ct1R","PkM130R"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-AspB510L.txt",sep="\t", quote = FALSE, row.names = FALSE)


## -------------

w = exactTest(dgList, pair=c("Ct_2","Asp_B510"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-AspB510L.txt",sep="\t", quote = FALSE, row.names = FALSE)

w = exactTest(dgList, pair=c("Ct_2","BP_CfBP3573"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-BpCfBP3573L.txt",sep="\t", quote = FALSE, row.names = FALSE)

w = exactTest(dgList, pair=c("Ct_2","Ec_DH5a"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-EcDH5aL.txt",sep="\t", quote = FALSE, row.names = FALSE)

w = exactTest(dgList, pair=c("Ct_2","Pd_WCS374r"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-PdWCS374rL.txt",sep="\t", quote = FALSE, row.names = FALSE)

w = exactTest(dgList, pair=c("Ct_2","Hs_SmR1"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-HsSmR1L.txt",sep="\t", quote = FALSE, row.names = FALSE)

w = exactTest(dgList, pair=c("Ct_2","MIX"))
tab= as.data.frame(topTags(w, n=nrow(w$table)))
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$FDR
write.table(mat, file = "Ct2L-MIXL.txt",sep="\t", quote = FALSE, row.names = FALSE)

## -------------

my.contrasts <- makeContrasts(
CvsK = Ct_1_L-Pk_M130_L,
CvsP = Ct_1_L-Pp_psJN_L,
CvsV = Ct_1_L-Bv_TVV75_L,
levels=design)

qlf_CvsK <- glmQLFTest(fit, contrast=my.contrasts[,"CvsK"])
tab = as.data.frame(qlf_CvsK$table)
FDR <- p.adjust(qlf_CvsK$table$PValue, method="fdr")
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=FDR
write.table(mat, file = "Ct1L-PkM130L.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_CvsP <- glmQLFTest(fit, contrast=my.contrasts[,"CvsP"])
tab = as.data.frame(qlf_CvsP$table)
FDR <- p.adjust(qlf_CvsP$table$PValue, method="fdr")
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=FDR
write.table(mat, file = "Ct1L-PppsJNL.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_CvsV <- glmQLFTest(fit, contrast=my.contrasts[,"CvsV"])
tab = as.data.frame(qlf_CvsV$table)
FDR <- p.adjust(qlf_CvsV$table$PValue, method="fdr")
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=FDR
write.table(mat, file = "Ct1L-BvTVV75L.txt",sep="\t", quote = FALSE, row.names = FALSE)