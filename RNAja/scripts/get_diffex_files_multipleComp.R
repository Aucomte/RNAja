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
gcsv_hisat <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all/3_count/STRINGTIE/HISAT_gene_count_matrix.csv"
sample_info <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/DATA/sample_info_all.txt"

normCount_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all/4_DE_analysis/normCounts.csv"
counts_file <- "/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all/4_DE_analysis/Counts.csv"

setwd("/home/comtea/Documents/PHIM/BURKADAPT/RNAja/output_all/4_DE_analysis/")

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

# filter
keep <- filterByExpr(dgList)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]

# normalisation
dgList <- calcNormFactors(dgList, method="TMM")

eff.lib.size <- dgList$samples$lib.size*dgList$samples$norm.factors
normCounts <- cpm(dgList)

# output for diffexdb

write.table(dgList$counts,file=counts_file,sep=";", quote = FALSE)
write.table(normCounts, file=normCount_file, sep=";", quote = FALSE)

# design

Condition1 = factor(sampleInfo[,2])
Condition1 = relevel(Condition1, ref="F")

Condition2 = factor(sampleInfo[,3])
Condition2 = relevel(Condition2, ref="C")

Group <- factor(paste(sampleInfo[,2], sampleInfo[,3],sep=""))
cbind(sampleInfo, Group)

design = model.matrix(~0+Group)
colnames(design) <- levels(Group)

dgList <- estimateDisp(dgList, design, robust=TRUE)
fit <- glmQLFit(dgList, design)

## -------------

my.contrasts <- makeContrasts(
R.CTRLvs121 = RCTRL-R121,
R.CTRLvs1858 = RCTRL-R1858,
R.CTRLvs2428 = RCTRL-R2428,
R.CTRLvs3267 = RCTRL-R3267,
R.CTRLvs3269 = RCTRL-R3269,
R.CTRLvs44 = RCTRL-R44,
R.CTRLvsMIX = RCTRL-RMIX,
F.CTRLvs121 = FCTRL-F121,
F.CTRLvs1858 = FCTRL-F1858,
F.CTRLvs2428 = FCTRL-F2428,
F.CTRLvs3267 = FCTRL-F3267,
F.CTRLvs3269 = FCTRL-F3269,
F.CTRLvs44 = FCTRL-F44,
F.CTRLvsMIX = FCTRL-FMIX,
levels=design)

qlf_R.CTRLvs121 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs121"])
tab = as.data.frame(qlf_R.CTRLvs121$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R121.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvs1858 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs1858"])
tab = as.data.frame(qlf_R.CTRLvs1858$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R1858.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvs2428 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs2428"])
tab = as.data.frame(qlf_R.CTRLvs2428$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R2428.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvs3267 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs3267"])
tab = as.data.frame(qlf_R.CTRLvs3267$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R3267.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvs3269 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs3269"])
tab = as.data.frame(qlf_R.CTRLvs3269$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R3269.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvs44 <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvs44"])
tab = as.data.frame(qlf_R.CTRLvs44$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_R44.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CTRLvsMIX <- glmQLFTest(fit, contrast=my.contrasts[,"R.CTRLvsMIX"])
tab = as.data.frame(qlf_R.CTRLvsMIX$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RCTRL_vs_RMIX.txt",sep="\t", quote = FALSE, row.names = FALSE)

#-

qlf_F.CTRLvs121 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs121"])
tab = as.data.frame(qlf_F.CTRLvs121$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F121.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvs1858 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs1858"])
tab = as.data.frame(qlf_F.CTRLvs1858$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F1858.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvs2428 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs2428"])
tab = as.data.frame(qlf_F.CTRLvs2428$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F2428.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvs3267 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs3267"])
tab = as.data.frame(qlf_F.CTRLvs3267$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F3267.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvs3269 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs3269"])
tab = as.data.frame(qlf_F.CTRLvs3269$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F3269.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvs44 <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvs44"])
tab = as.data.frame(qlf_F.CTRLvs44$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_F44.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CTRLvsMIX <- glmQLFTest(fit, contrast=my.contrasts[,"F.CTRLvsMIX"])
tab = as.data.frame(qlf_F.CTRLvsMIX$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FCTRL_vs_FMIX.txt",sep="\t", quote = FALSE, row.names = FALSE)


## -------------

my.contrasts <- makeContrasts(
R.CvsK = RC-RK,
R.CvsP = RC-RP,
R.CvsV = RC-RV,
F.CvsK = FC-FK,
F.CvsP = FC-FP,
F.CvsV = FC-FV,
RvsF.C = RC-FC,
RvsF.K = RK-FK,
RvsF.P = RP-FP,
RvsF.V = RV-FV,
levels=design)

qlf_R.CvsK <- glmQLFTest(fit, contrast=my.contrasts[,"R.CvsK"])
FDR <- p.adjust(qlf_R.CvsK$table$PValue, method="BH")
tab = as.data.frame(qlf_R.CvsK$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=FDR
write.table(mat, file = "RC_vs_RK.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CvsP <- glmQLFTest(fit, contrast=my.contrasts[,"R.CvsP"])
tab = as.data.frame(qlf_R.CvsP$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RC_vs_RP.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_R.CvsV <- glmQLFTest(fit, contrast=my.contrasts[,"R.CvsV"])
tab = as.data.frame(qlf_R.CvsV$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RC_vs_RV.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CvsK <- glmQLFTest(fit, contrast=my.contrasts[,"F.CvsK"])
tab = as.data.frame(qlf_F.CvsK$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FC_vs_FK.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CvsP <- glmQLFTest(fit, contrast=my.contrasts[,"F.CvsP"])
tab = as.data.frame(qlf_F.CvsP$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FC_vs_FP.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_F.CvsV <- glmQLFTest(fit, contrast=my.contrasts[,"F.CvsV"])
tab = as.data.frame(qlf_F.CvsV$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "FC_vs_FV.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_RvsF.C <- glmQLFTest(fit, contrast=my.contrasts[,"RvsF.C"])
tab = as.data.frame(qlf_RvsF.C$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RC_vs_FC.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_RvsF.P <- glmQLFTest(fit, contrast=my.contrasts[,"RvsF.P"])
tab = as.data.frame(qlf_RvsF.P$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RP_vs_FP.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_RvsF.V <- glmQLFTest(fit, contrast=my.contrasts[,"RvsF.V"])
tab = as.data.frame(qlf_RvsF.V$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RV_vs_FV.txt",sep="\t", quote = FALSE, row.names = FALSE)

qlf_RvsF.K <- glmQLFTest(fit, contrast=my.contrasts[,"RvsF.K"])
tab = as.data.frame(qlf_RvsF.K$table)
mat = matrix(ncol = 6, nrow= nrow(tab))
colnames(mat) = c("genes","logFC","logCPM","NA","PValue","FDR")
mat[,1]=rownames(tab)
mat[,2]=tab$logFC
mat[,3]=tab$logCPM
mat[,4]="NA"
mat[,5]=tab$PValue
mat[,6]=tab$F
write.table(mat, file = "RK_vs_FK.txt",sep="\t", quote = FALSE, row.names = FALSE)