### R codes for mutation timing
rm(list = ls())

# Install dependency packages and load in data
#devtools::install_github("gerstung-lab/MutationTimeR")
library(MutationTimeR)
library(tidyverse)

download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-SCNA-segment.csv")
download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-table.csv")
download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-WES.maf")
OCCC.seg <- read.csv("CGMH-OCCC-SCNA-segment.csv")
OCCC.mut <- read.delim("CGMH-OCCC-WES.maf")
OCCC.table <- read.csv("CGMH-OCCC-table.csv")

### time the somatic mutation alone with the SCNA status using MutationTimeR package
mut.timing <- function(OCCC.ID){
  cat(paste("Timing gene mutation for", OCCC.ID, "\n"))
  tmp.path <- paste0(OCCC.ID,".temp.vcf")
  
  mut.cat <- subset(OCCC.mut, Tumor_Sample_Barcode == OCCC.ID) %>% 
    subset(Chromosome %in% 1:22) %>% 
    separate(col = "Otherinfo.13", into = c("GT","AD","AF","DP"), sep = ":", remove = F, convert = T) %>% 
    mutate(Alt_c = as.numeric(str_extract(AD, "\\d*$"))) %>% 
    mutate(Ref_c = DP - Alt_c)
  
  purity <- OCCC.table$CCF[OCCC.table$ID == OCCC.ID]
  ploidy <- OCCC.table$Tumor.ploidy[OCCC.table$ID == OCCC.ID]
  WGD <- ifelse(OCCC.table$GD.status[OCCC.table$ID == OCCC.ID] == "GD", TRUE, FALSE)
  
  seg <- subset(OCCC.seg, ID == OCCC.ID)
  seg$chromosome <- sub("chr","",seg$chromosome)
  
  ## for gene mutations with no available CN estimation, use a average CN at arm level
  mut.cat$seq.range <- sapply(1:nrow(mut.cat), function(x){
    df <- subset(seg, chromosome == mut.cat$Chr & start.pos < mut.cat$Start & end.pos > mut.cat$End)
    ifelse(nrow(df)>0, "yes", "no")
  }) 
  mut.cat.arm <- subset(mut.cat, seq.range == "no")
  mut.cat <- subset(mut.cat, seq.range == "yes")
  
  cat(paste("Number of gene mut included:", nrow(mut.cat),"\n"))
  cat(paste("Number of gene mut not in sequenza range:", nrow(mut.cat.arm),"\n"))
  cat(paste("Estimated purity:", purity, "\n"))
  cat(paste("Estimated ploidy:", round(ploidy, digits = 2),"\n"))
  cat(paste("Set WGD to", WGD, "\n\n"))
  
  if(nrow(mut.cat)>0){
    vcf.tmp <- data.frame(
      Chr=mut.cat$Chromosome, Start=mut.cat$Start_Position, ID = ".", 
      Ref=mut.cat$Reference_Allele, Alt=mut.cat$Tumor_Seq_Allele2, 
      Qual = ".", FILTER = "PASSED",
      info = paste0("t_alt_count=",mut.cat$Alt_c,";t_ref_count=",mut.cat$Ref_c),
      Format = "AD:DP", Sample = paste0(".,",mut.cat$Alt_c,":",mut.cat$DP)
    )
    
    header <- readLines("Data_vcf_header.txt") 
    cat(header, file = tmp.path, sep = "\n")
    write.table(vcf.tmp, file = tmp.path, append = T, quote = F, sep = "\t", 
                row.names = F, col.names = F)
    
    vcf <- MutationTimeR::readVcf(tmp.path) 
    unlink(tmp.path)
    
    gr <- GRanges(Rle(c(seg$chromosome)), IRanges(start = seg$start.pos, end = seg$end.pos),
                  major_cn=as.integer(seg$CNa), minor_cn=as.integer(seg$CNb), 
                  clonal_frequency=purity)
    
    ### running mutationTimer package
    clusters <- data.frame(cluster = 1:2, proportion = purity*c(1, 0.5), n_ssms=c(90,10))
    
    mt <- mutationTime(vcf = vcf, cn = gr, isWgd = WGD, n.boot=1000, clusters = clusters)
    vcf <- addMutTime(vcf, mt$V)
    vcf.mut.timing <- info(vcf) %>% data.frame()
    vcf.mut.timing <- mutate(vcf.mut.timing, mutation = rownames(vcf.mut.timing), ID = OCCC.ID, 
                             WGD = WGD, purity = purity, ploidy = ploidy)
    
    mut.cat <- mutate(
      mut.cat, 
      mutation = str_glue("{Chromosome}:{Start_Position}_{Reference_Allele}/{Tumor_Seq_Allele2}")
    )
    vcf.mut.timing <- left_join(
      vcf.mut.timing, 
      dplyr::select(mut.cat,mutation,Hugo_Symbol,Variant_Classification,AF,DP,seq.range),
      by = "mutation") %>% arrange(CLS, desc(AF))
  }
  
  ### summarize SCNA seg to arm levels if some gene mutations have no CN estimation
  if(nrow(mut.cat.arm)>0){
    CNV.region <- read.delim("Data_GRCh38p7_region.tsv")
    seg.arm <- subset(OCCC.seg, ID == OCCC.ID)
    seg.arm$chromosome <- sapply(1:nrow(seg.arm), function(x){
      arm <- ifelse(seg.arm$start.pos[x] < CNV.region$centromere_start[CNV.region$chr == seg.arm$chromosome[x]],
                    "", ".5")
      return(paste0(seg.arm$chromosome[x],arm))
    })
    
    tmp <- lapply(unique(seg.arm$chromosome), function(x){
      df <- subset(seg.arm, chromosome == x) 
      seg.chr <- data.frame(chromosome = x, start.pos = min(df$start.pos), end.pos = max(df$end.pos),
                            major_CN = round(sum(df$CNa*(df$end.pos-df$start.pos))/sum(df$end.pos-df$start.pos, digits = 0)),
                            minor_CN = round(sum(df$CNb*(df$end.pos-df$start.pos))/sum(df$end.pos-df$start.pos, digits = 0)))
      return(seg.chr)
    })
    tmp <- do.call("rbind", tmp)
    tmp$chromosome <- sub("\\.5","",tmp$chromosome) %>% sub(pattern = "chr", replacement = "")
    
    gr.arm <- GRanges(Rle(c(tmp$chromosome)), 
                      IRanges(start = tmp$start.pos, end = tmp$end.pos),
                      major_cn=as.integer(tmp$major_CN), minor_cn=as.integer(tmp$minor_CN), 
                      clonal_frequency=purity)
    
    vcf.tmp <- data.frame(Chr=mut.cat.arm$Chromosome, Start=mut.cat.arm$Start_Position, ID = ".", 
                          Ref=mut.cat.arm$Reference_Allele, Alt=mut.cat.arm$Tumor_Seq_Allele2, 
                          Qual = ".", FILTER = "PASSED",
                          info = paste0("t_alt_count=",mut.cat.arm$Alt_c,";t_ref_count=",mut.cat.arm$Ref_c),
                          Format = "AD:DP", Sample = paste0(".,",mut.cat.arm$Alt_c,":",mut.cat.arm$DP))

    header <- readLines("Data_vcf_header.txt") 
    cat(header, file = tmp.path, sep = "\n")
    write.table(vcf.tmp, file = tmp.path, append = T, quote = F, sep = "\t", 
                row.names = F, col.names = F)
    
    vcf.arm <- readVcf(tmp.path) 
    unlink(tmp.path)
    
    clusters.arm <- data.frame(cluster = 1:2, proportion = purity*c(1, 0.5), n_ssms=c(90,10))
    
    mt.arm <- mutationTime(vcf = vcf.arm, cn = gr.arm, isWgd = WGD, n.boot=1000, 
                           clusters = clusters.arm)
    vcf.arm <- addMutTime(vcf.arm, mt.arm$V)
    vcf.arm.timing <- info(vcf.arm) %>% data.frame()
    vcf.arm.timing <- mutate(vcf.arm.timing, mutation = rownames(vcf.arm.timing), ID = OCCC.ID, 
                             WGD = WGD, purity = purity, ploidy = ploidy)
    
    mut.cat.arm <- mutate(
      mut.cat.arm, 
      mutation = str_glue("{Chromosome}:{Start_Position}_{Reference_Allele}/{Tumor_Seq_Allele2}"))
    vcf.arm.timing <- left_join(vcf.arm.timing, 
                                dplyr::select(mut.cat.arm,mutation,Hugo_Symbol,
                                              Variant_Classification,AF,DP,seq.range),
                                by = "mutation") %>% arrange(CLS, desc(AF))
  }
  
  if(nrow(mut.cat)>0 & nrow(mut.cat.arm)>0){
    vcf.mut.timing <- rbind(vcf.mut.timing, vcf.arm.timing) %>% arrange(CLS, desc(AF))
    return(vcf.mut.timing)
  } else if(nrow(mut.cat)>0 & nrow(mut.cat.arm)==0) {
    return(vcf.mut.timing)
  } else if(nrow(annovar)==0 & nrow(annovar.arm)>0) {
    return(vcf.arm.timing)
  }
  
}

OCCC.mut.timing <- lapply(unique(OCCC.table$ID), mut.timing)

OCCC.mut.timing <- do.call("rbind", OCCC.mut.timing)
