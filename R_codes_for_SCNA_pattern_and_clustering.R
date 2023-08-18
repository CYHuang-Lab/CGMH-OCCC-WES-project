### R codes for SCNA analysis and consensus clustering
rm(list = ls())

# library dependency packages and load in data
library(tidyverse)
library(ggpubr)
library(ConsensusClusterPlus)
library(ComplexHeatmap)

download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-SCNA-segment.csv")
download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-table.csv")
download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/Data_GRCh38p7_region.tsv")
OCCC.seg <- read.csv("CGMH-OCCC-SCNA-segment.csv")

### summarize CNV at 5MB resolution using sequenza output
seg.reader <- function(OCCC.ID, res.mb = 5){
  cat(paste("Reading SCNA segments :", OCCC.ID, "\n"))
  res <- res.mb*1000000
  seg <- subset(OCCC.seg, ID == OCCC.ID)
  ctm <- read.delim("Data_GRCh38p7_region.tsv")
  
  ## This function will exclude centromere regions from downstream analysis
  seg.nc <- lapply(paste0("chr",1:22), function(x){
    seg.chr <- subset(seg, chromosome == x)
    df1 <- subset(seg.chr, end.pos <= ctm$centromere_start[ctm$chr==x])
    df2 <- subset(seg.chr, start.pos > ctm$centromere_end[ctm$chr==x])
    
    df3 <- subset(seg.chr, start.pos <= ctm$centromere_start[ctm$chr==x] & 
                    end.pos > ctm$centromere_start[ctm$chr==x] & 
                    end.pos <= ctm$centromere_end[ctm$chr==x])
    if(nrow(df3)>0){
      df3$end.pos <- ctm$centromere_start[ctm$chr==x]
    }
    
    df4 <- subset(seg.chr, end.pos > ctm$centromere_end[ctm$chr==x] & 
                    start.pos > ctm$centromere_start[ctm$chr==x] & 
                    start.pos <= ctm$centromere_end[ctm$chr==x])
    if(nrow(df4)>0){
      df4$start.pos <- ctm$centromere_end[ctm$chr==x]
    }
    
    df5 <- subset(seg.chr, start.pos <= ctm$centromere_start[ctm$chr==x] & 
                    end.pos > ctm$centromere_end[ctm$chr==x])
    if(nrow(df5)>0){
      df5 <- df5[c(1,1),]
      df5$end.pos[1] <- ctm$centromere_start[ctm$chr==x]
      df5$start.pos[2] <- ctm$centromere_end[ctm$chr==x]
    }
    
    df <- rbind(df1,df2,df3,df4,df5) %>% arrange(start.pos)
    return(df)
  }) 
  seg <- do.call("rbind", seg.nc)
  
  seg.cn <- lapply(paste0("chr",1:22), function(x){
    seg.chr <- subset(seg, chromosome == x)
    max.bin <- ceiling(max(seg.chr$end.pos)/res)
    seg.bin <- seq(0, res*max.bin, res)
    
    seg.cn <- lapply(1:max.bin, function(x){
      df1 <- subset(seg.chr, start.pos > seg.bin[x] & end.pos <= seg.bin[x+1])
      
      df2 <- subset(seg.chr, start.pos <= seg.bin[x] & end.pos > seg.bin[x+1])
      if(nrow(df2)>0){
        df2$start.pos <- seg.bin[x]
        df2$end.pos <- seg.bin[x+1]
      }
      
      df3 <- subset(seg.chr, start.pos <= seg.bin[x] & end.pos > seg.bin[x] & end.pos <= seg.bin[x+1])
      if(nrow(df3)>0){
        df3$start.pos <- seg.bin[x]
      }
      
      df4 <- subset(seg.chr, start.pos > seg.bin[x] & start.pos <= seg.bin[x+1] & end.pos > seg.bin[x+1])
      if(nrow(df4)>0){
        df4$end.pos <- seg.bin[x+1]
      }
      
      df <- rbind(df1, df2, df3, df4)
      if(nrow(df)>0){
        CNa.sum <- (sum(df$CNa*(df$end.pos-df$start.pos))/sum(df$end.pos-df$start.pos)) %>% 
          round()
        CNb.sum <- (sum(df$CNb*(df$end.pos-df$start.pos))/sum(df$end.pos-df$start.pos)) %>% 
          round()
        CNt.sum <- CNa.sum + CNb.sum
      } else{
        CNa.sum <- NA
        CNb.sum <- NA
        CNt.sum <- NA
      }
      df <- data.frame(START = seg.bin[x]+1, END = seg.bin[x+1], 
                       CNa.sum = CNa.sum, CNb.sum = CNb.sum, CNt.sum = CNt.sum)
      return(df)
    })
    seg.cn <- do.call("rbind", seg.cn)
    seg.cn$CHR <- x
    return(seg.cn)
  })
  seg.cn <- do.call("rbind",seg.cn)
  seg.cn$ID <- OCCC.ID
  return(seg.cn)
}
OCCC.seg.binned <- lapply(unique(OCCC.seg$ID)[1:10], seg.reader)
OCCC.seg.binned <- do.call("rbind", OCCC.seg.binned)
OCCC.table <- read.csv("CGMH-OCCC-table.csv")

OCCC.seg.binned <- left_join(
  OCCC.seg.binned, OCCC.table[,c("ID","BG.ploidy")], by = "ID") %>% 
  mutate(CNV = CNt.sum - BG.ploidy, 
         START = as.integer(START), END = as.integer(END)) %>% 
  mutate(SEG = str_glue("{CHR}:{START}-{END}")) %>% 
  dplyr::select(ID, CHR, START, END, SEG, starts_with("CN"), BG.ploidy)

OCCC.seg.binned$CNa.sum[is.nan(OCCC.seg.binned$CNa.sum)] <- NA
OCCC.seg.binned$CNb.sum[is.nan(OCCC.seg.binned$CNb.sum)] <- NA
OCCC.seg.binned$CNt.sum[is.nan(OCCC.seg.binned$CNt.sum)] <- NA
OCCC.seg.binned$CNV[is.nan(OCCC.seg.binned$CNV)] <- NA

OCCC.seg.binned <- mutate(OCCC.seg.binned, CNV.status = case_when(
  !is.na(CNt.sum) & CNa.sum == 0 & CNb.sum == 0 ~ "CN.del",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 1 & CNb.sum > 0 ~ "CN1.gain.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 1 & CNb.sum == 0 ~ "CN1.gain.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV >= 2 & CNb.sum > 0 ~ "CN2.gain.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV >= 2 & CNb.sum == 0 ~ "CN2.gain.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == -1 & CNb.sum > 0 ~ "CN1.loss.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == -1 & CNb.sum == 0 ~ "CN1.loss.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV <= -2 & CNb.sum > 0 ~ "CN2.loss.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV <= -2 & CNb.sum == 0 ~ "CN2.loss.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 0 & CNb.sum == 0 ~ "CN0.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 0 & CNb.sum > 0 ~ "CN.neutral",
  TRUE ~ "Unknown"
))

CNV.sum <- data.frame(SEG = unique(OCCC.seg.binned$SEG)) %>% 
  separate(col = "SEG", into = c("CHR","RANGE"), sep = ":", remove = F) %>% 
  separate(col = "RANGE", into = c("START","END"), sep = "-", remove = F, convert = T)

CNV.sum <- rbind(CNV.sum, data.frame(
  SEG = paste0("chr",1:22,".space1"), CHR = paste0("chr",1:22), RANGE = NA, 
  START = max(CNV.sum$START) + 1, END = NA)) %>% 
  rbind(data.frame(
    SEG = paste0("chr",1:22,".space2"), CHR = paste0("chr",1:22), RANGE = NA, 
    START = max(CNV.sum$START) + 2, END = NA)) %>% 
  mutate(CHR = factor(CHR, levels = paste0("chr",1:22))) %>% arrange(CHR, START)

CNV.sum$pos <- 1:nrow(CNV.sum)
plot.chr <- group_by(CNV.sum, CHR) %>% 
  summarize(width = length(CHR), break.start = min(pos)-1, break.end = max(pos)) %>% 
  mutate(pos = break.end-0.5*width)

CNV.gain <- count(OCCC.seg.binned, SEG, CNV.status) %>% 
  separate(col = "SEG", into = c("CHR","RANGE"), sep = ":", remove = F) %>% 
  separate(col = "RANGE", into = c("START","END"), sep = "-", remove = F, convert = T) %>% 
  mutate(CHR = factor(CHR, levels = paste0("chr",1:22))) %>% 
  subset(str_detect(CNV.status, "gain")) %>% 
  arrange(CHR, START) %>% 
  mutate(CNV.status = factor(CNV.status, levels = c("CN2.gain.loh","CN2.gain.het",
                                                    "CN1.gain.loh","CN1.gain.het")),
         SEG = factor(SEG, levels = unique(CNV.sum$SEG)))
CNV.gain <- left_join(CNV.gain, CNV.sum[,c("SEG","pos")], by = "SEG")

CNV.loss <- count(OCCC.seg.binned, SEG, CNV.status) %>% 
  separate(col = "SEG", into = c("CHR","RANGE"), sep = ":", remove = F) %>% 
  separate(col = "RANGE", into = c("START","END"), sep = "-", remove = F, convert = T) %>% 
  mutate(CHR = factor(CHR, levels = paste0("chr",1:22))) %>% 
  subset(str_detect(CNV.status, "loss|CN0|del")) %>% 
  arrange(CHR, START) %>% 
  mutate(CNV.status = factor(CNV.status, levels = c("CN.del","CN2.loss.loh","CN2.loss.het",
                                                    "CN1.loss.loh","CN1.loss.het","CN0.loh")),
         SEG = factor(SEG, levels = unique(CNV.sum$SEG)))
CNV.loss <- left_join(CNV.loss, CNV.sum[,c("SEG","pos")], by = "SEG")


p1 <- ggplot(CNV.gain) + 
  geom_rect(data=plot.chr, mapping=aes(xmin=break.start, xmax=break.end, ymin=0, ymax=80),
            fill = rep(c("gray90","white"),11)) + 
  geom_hline(yintercept = seq(0,60,20), color = "gray80", size=0.25) + 
  geom_col(aes(x=pos, y=100*n/n_distinct(OCCC.seg$ID), fill=CNV.status), width=1, show.legend = F) + 
  geom_segment(x=0, xend = 20, y = 70, yend = 70, color = "black", size = 1) + 
  geom_text(x=10, y=72, label="100MB", size = 4) + 
  scale_x_continuous(breaks = c(0,plot.chr$break.end), minor_breaks = NULL, labels = NULL) +
  theme_void() + theme(axis.title = element_blank()) + 
  scale_fill_manual(values = c("CN2.gain.loh"="red3","CN2.gain.het"="salmon1",
                               "CN1.gain.loh"="sienna","CN1.gain.het"="tan1"))

p2 <- ggplot(plot.chr) + 
  geom_rect(aes(xmin=break.start, xmax=break.end, ymin=0, ymax=1), 
            fill = rep(c("gray80","white"),11), color="black") + 
  geom_text(aes(x=pos, y=0.5, label=1:22), size = 4) + 
  theme_void()

p3 <- ggplot(CNV.loss) + 
  geom_rect(data=plot.chr, mapping=aes(xmin=break.start, xmax=break.end, ymin=0, ymax=-80),
            fill = rep(c("gray90","white"),11)) + 
  geom_hline(yintercept = c(0,-20,-40), color = "gray80", size=0.25) + 
  geom_col(aes(x=pos, y=-100*n/n_distinct(OCCC.seg$ID), fill=CNV.status), width=1, show.legend = F) + 
  scale_x_continuous(breaks = c(0,plot.chr$break.end), minor_breaks = NULL, labels = NULL) +
  theme_void() + theme(axis.title = element_blank()) + 
  scale_fill_manual(values = c("CN2.loss.loh"="steelblue4","CN2.loss.het"="steelblue1",
                               "CN1.loss.loh"="limegreen","CN1.loss.het"="palegreen",
                               "CN0.loh"="gray50", "CN.del"="darkviolet"))

ggpubr::ggarrange(p1,p2,p3,heights = c(12,1,12), ncol=1, align="v")


## perform consensus clustering
OCCC.SCNA <- OCCC.seg.binned[,c("ID","SEG","CNV")] %>% spread(SEG, CNV, fill=0)
rownames(OCCC.SCNA) <- OCCC.SCNA[,1]
OCCC.SCNA <- OCCC.SCNA[,-1]

cluster.results <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = t(OCCC.SCNA),    
  maxK = 6,            
  reps = 10000,        
  pItem = 0.8,         
  pFeature = 1,        
  title = "./cluster.km",
  clusterAlg = "km",   
  distance = "euclidean", 
  seed = 1220,         
  plot = "pdf"
)      

icl <- calcICL(cluster.results, title="./cluster.k.3", plot="pdf")

## choose k = 3 clustering results 
clus.km3 <- cluster.results[[3]]$consensusClass   

ID.order <- unlist(sapply(1:3,function(x){names(clus.km3)[clus.km3==x]}))

rownames(OCCC.table) <- OCCC.table$ID

ComplexHeatmap::Heatmap(
  as.matrix(OCCC.SCNA[ID.order,grep("space",CNV.sum$SEG,invert = T,value = T)]), 
  show_column_names = F, 
  cluster_columns = F, 
  cluster_rows = T, 
  column_split = rep(1:22, plot.chr$width-2), row_split = rep(paste("Cluster",1:3), table(clus.km3)), 
  row_names_gp = gpar(fontsize = 6), 
  col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "gray95", "red")), na_col = "white",
  right_annotation = rowAnnotation(
    Mut.sig.group = OCCC.table[ID.order,"MS.Group"],
    col = list(
      Mut.sig.group = c("Unclassified"="gray80","MMR"="blue","APOBEC"="green",
                        "HRD"="yellow", "NHEJ"="tomato")
    )
  ))
