rm(list = ls())

#library(devtools)
#devtools::install_github("boxiangliu/locuscomparer")


#引用包
library(dplyr)
library(data.table)
library(coloc)
library(vroom)
library(locuscomparer)

#输入文件
qtlFile="18840_205_SH2D1B_SH21B.txt.gz"       #基因或者蛋白数据
gwasFile='finngen_R10_AB1_OTHER_SEPSIS.gz'     #结局数据

geneChr="chr1"             #基因所在的染色体(https://www.ncbi.nlm.nih.gov/gene)
geneStart=162395268        #染色体起始位置
geneEnd=162412136          #染色体终止位置

#设置工作目录
setwd("C:\\Users\\lexb\\Desktop\\Draco")

##读取整理eQTL
data1=vroom(qtlFile)
#挑选列名
data2 <- data1 %>% dplyr::select("rsids","Chrom","Pos","effectAllele","otherAllele",
                                 "ImpMAF","Beta","SE","Pval","N")
#更改列名
colnames(data2) <- c('SNP','chrom',"Pos",'effect_allele','other_allele',
                     "MAF","beta","se","P","samplesize")
#整理共定位所需数据
data3 <- as.data.frame(data2)
data3$varbeta <- data3$se^2
data3$z = data3$beta/data3$se
data3 <- subset(data3, !duplicated(SNP))
#此处需要在NCBI查询38版本启动子位点上下1MB提取cis-PQTL
geneData <- data3 %>% filter(chrom==geneChr, Pos>geneStart-1000000, Pos<geneEnd+1000000) %>% na.omit()
lead <- geneData %>% dplyr::arrange(P)
leadPos <- lead$Pos[1]
QTLdata <-  geneData %>% filter(Pos>leadPos-500000,Pos<leadPos+500000) %>% na.omit()
##读取整理GWAS，以下同上
data0 <- vroom(gwasFile)
finn_info <- fread('R10_manifest.xls',data.table = F)
trait_row <- finn_info[grepl(gwasFile, finn_info$path_https),]
data0$ncase.outcome <- trait_row$num_cases
data0$ncontrol.outcome <- trait_row$num_controls
data0$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

data1 <- data0 %>% dplyr::select("rsids","#chrom","pos","alt","ref",
                                 "af_alt","beta","sebeta","pval","samplesize.outcome","ncase.outcome")
colnames(data1) <- c('SNP','chrom',"pos",'effect_allele','other_allele',
                     "eaf","beta","se","P","samplesize","number_cases")
data2 <- as.data.frame(data1)
data2$varbeta <- data2$se^2
data2$MAF <- ifelse(data2$eaf<0.5,data2$eaf,1-data2$eaf)
data3 <- subset(data2, !duplicated(SNP))
data3$s <- data3$number_cases/data3$samplesize
data3$z = data3$beta/data3$se
GWASdata <- data3 %>% na.omit()

sameSNP <- intersect(QTLdata$SNP,GWASdata$SNP)
QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
GWASdata <- GWASdata[GWASdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
#共定位分析
coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,
                                 N=QTLdata$samplesize,MAF=QTLdata$MAF,z = QTLdata$z,
                                 pvalues=QTLdata$P,type="quant"), 
                   dataset2=list(snp=GWASdata$SNP,beta=GWASdata$beta,varbeta=GWASdata$varbeta,
                                 N=GWASdata$samplesize,MAF=GWASdata$MAF,z = GWASdata$z,
                                 pvalues=GWASdata$P,type="cc"))

result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
#一般用于GWAS 多出一项coloc_data[["dataset2"]][["s"]] <- NULL
result$results %>% dplyr::arrange(desc(SNP.PP.H4))
result

#绘图
GWAS_fn <- GWASdata[,c('SNP','P')] %>%
  dplyr::rename(rsid = SNP, pval = P)
pqtl_fn <- QTLdata[,c('SNP','P')] %>%
  dplyr::rename(rsid = SNP, pval = P)

#输出图形
pdf(file="coloc.pdf", width=8, height=6)
print(locuscompare(in_fn1 = GWAS_fn, 
                   in_fn2 = pqtl_fn,
                   snp = NULL,
                   title1 = 'GWAS', 
                   title2 = 'pqtl'))
dev.off()




