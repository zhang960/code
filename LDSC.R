setwd(dir="D:/mendelian/mendelian")  #设置你自己的工作路径
rm(list = ls())

# -加载需要用到的包，没安装的需要提前安装一下
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(dplyr)


# 暴露，需要用的文件放在你设置的工作路径下
FileNames <-list.files(paste0(getwd()),pattern="txt.gz$")

# 结局
out<-fread("convert_ieu-b-4980.vcf.gz.csv")#读取本地结局


out <- data.frame(out)
##对数据进行格式转换，注意列名
outcome<-format_data(out,type="outcome",
                     snp_col = "SNP",
                     effect_allele_col = "effect.allele",
                     other_allele_col = "other.allele",
                     beta_col = "beta",
                     se_col = "se",
                     samplesize_col = "samplesize",
                     pval_col = "pval",
                     eaf_col = "eaf",
                     chr_col = "chr",
                     pos_col = "pos"
)
rm(out)
outcome$outcome <- "Sepsis"



#######循环代码######################

#FileNames
#1:length(FileNames)中的1可以改，看你前面做了几个，从而接着做

for (i in 1:length(FileNames)) {
  
  
  #自己改i的值，就是做第几个。把循环给你取消了
  # i <- 1  
  
  cat("正在运行第",i, "个\n")
  exp_dat_id <- FileNames[i]
  exp <- sub("\\.txt\\.gz", "", exp_dat_id)
  
  ##逐个读入暴露数据
  cat("读入暴露数据\n")
  exposure<- fread(FileNames[i])
  
  exposure$Phenotype <- exp
  exposure <- data.frame(exposure)
  cat("格式转换\n")
  # head(exposure)
    ###数据格式转换，注意列名#
    exposure<-format_data(exposure,type="exposure",
                          snp_col = "rsids",
                          phenotype_col = "Phenotype",
                          effect_allele_col = "effectAllele",
                          other_allele_col = "otherAllele",
                          beta_col = "Beta",
                          se_col = "SE",
                          samplesize_col = "N",
                          pval_col = "Pval",
                          eaf_col = "ImpMAF",
                          chr_col = "Chrom",
                          pos_col = "Pos"
    )
    

  
  ###依赖文件放到工作路径内
  ld="d:\\mendelian\\mendelian\\reference\\1000Gv3EUR\\"
  wld="d:\\mendelian\\mendelian\\reference\\LDSC\\reference\\eur_w_ld_chr\\"
  


  LDSC_rg<-function(expo,outcome,an,sample_prev=NA,
                    population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=200){
    id.o<-outcome$outcome[1]
    id.e<-expo$exposure[1]
    
    expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
    expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                        ,A1=effect_allele.exposure
                        ,A2=other_allele.exposure)
    expo<-as_tibble(expo)
    
    out<-outcome%>%mutate(Z=beta.outcome/se.outcome)
    out<-out%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                              ,A1=effect_allele.outcome
                              ,A2=other_allele.outcome)
    out<-as_tibble(out)
    
    
    dat<-list(expo,out)
    names(dat)<-c(id.e,id.o)
    
    rm(expo)
    
    
    res<-try(ldscr::ldsc_rg(dat,ancestry = an,sample_prev=sample_prev,
                            population_prev=population_prev,ld=ld,wld=wld,
                            n_blocks=n_blocks,chr_filter=chr_filter))
    
    return(res)
    
  }
  cat("LDSC分析\n")
  LDSC_res<-LDSC_rg(exposure,outcome,ld=ld,wld=wld,an="EUR")
  
  h2 <- LDSC_res[["h2"]]
  rg <- LDSC_res[["rg"]]
  
  
  
  fwrite(h2,file = paste0("results/",exp,"_h2.csv"))
  fwrite(rg,file = paste0("results/",exp,"_rg.csv"))
  
}

# ##结果合并
# fs=list.files("D:/mendelian randomization/test/test/mendelian", pattern = "rg.csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
# df = map_dfr(fs, read.csv)
# write.csv(df,"rg_res.csv")
# 
# fs=list.files("D:/mendelian randomization/test/test/mendelian", pattern = "h2.csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
# df = map_dfr(fs, read.csv)
# write.csv(df,"h2_res.csv")



