### R codes for survival analysis 
rm(list = ls())

# Install dependency packages and load in data
library(dplyr)
library(survival)
library(survminer)

OCCC.clinical <- read.csv("./CGMH-OCCC-table.csv")

OCCC.clinical <- mutate(OCCC.clinical, 
                        stage.group = ifelse(Stage %in% c("I","IA","IB","IC","II","IIA","IIB"),
                                             "Early stage","Late stage")
)
rownames(OCCC.clinical) <- OCCC.clinical$ID

gene.panel <- c("ARID1A","PIK3CA","TERT","KRAS","PPP2R1A","TP53","ATM",
                "CHD4","PIK3R1","ARID1B","CREBBP","KMT2C","MUC16","NFE2L2","SPOP")

plot.gene.OS <- function(gene, plot.surv = c(TRUE, FALSE)){
  #browser()
  tmp <- OCCC.clinical[,c("ID","OS.month","OS.event",gene)]
  colnames(tmp)[4] <- "gene"
  
  if(gene == "TERT"){
    tmp$gene <- factor(tmp$gene, levels = c("Wild type","upstream","UTR5"))
    pal.mut <- c("gray50","orange","pink")
  } else{
    tmp$gene <- factor(ifelse(tmp$gene == "Wild type", "Wild type", "Mutant"), 
                              levels = c("Wild type","Mutant"))
    pal.mut <- c("gray50","maroon1")
  }
  
  y <- coxph(formula = Surv(OS.month, OS.event == "1") ~ gene, data = tmp)
  z <- summary(y)
  HR <- round(z$coefficients[2],2)
  CI.lower <- round(z$conf.int[3],2)
  CI.upper <- round(z$conf.int[4],2)
  p.value <- ifelse(z$coefficients[5]<0.0001, z$coefficients[5], round(z$coefficients[5],4))
  
  if(isTRUE(plot.surv)){
    p <- survminer::ggsurvplot(
      data = tmp, fit = survfit(Surv(OS.month, OS.event == "1") ~ gene, data = tmp), 
      pval = T, surv.median.line = "hv", size = 0.75,
      legend = c(0.8, 0.3), legend.title = gene, 
      font.legend = list(size = 10, color = "black"),
      legend.labs = paste0(names(table(tmp$gene)), ", n = ", table(tmp$gene)),
      palette = pal.mut, risk.table = F, xlab = "", ylab = "", break.x.by = 12)
    
    return(p)
  } else if(!isTRUE(plot.surv) & gene == "TERT"){
    HR.table <- data.frame(
      Gene = c("TERT.promoter","TERT.UTR5"), Outcome = rep("OS",2), 
      MT.n = as.numeric(table(tmp$gene)[c("upstream","UTR5")]), 
      WT.n = rep(table(tmp$gene)["Wild type"],2), 
      HR = round(z$coefficients[,2], 2), 
      CI.lower = round(z$conf.int[,3],2), CI.upper = round(z$conf.int[,4],2), 
      p.value = z$coefficients[,5])
  } else {
    HR.table <- data.frame(
      Gene = gene, Outcome = "OS", 
      MT.n = table(tmp$gene)["Mutant"], 
      WT.n = table(tmp$gene)["Wild type"], 
      HR = HR, CI.lower = CI.lower, CI.upper = CI.upper, p.value = z$coefficients[5])
  } 
  return(HR.table)
}
plot.gene.OS(gene = "TERT", plot.surv = F)

OS.table <- lapply(gene.panel, function(x){
  plot.gene.OS(gene=x, plot.surv = F)})

OS.table <- do.call("rbind", OS.table)

OS.table$q.value <- p.adjust(OS.table$p.value, method = "BH")
rownames(OS.table) <- OS.table$Gene
write.csv(OS.table, "gene.mut.OS.csv", row.names = F, quote = T)



## Cox hazard model uni/multivariate analysis
library(MASS)

OCCC.clinical$TERT <- factor(OCCC.clinical$TERT, levels = c("Wild type", "upstream", "UTR5"))

fit  <- coxph(Surv(OS.month, OS.event == "1") ~ Age + stage.group + TERT, data = OCCC.clinical)

fit0 <- coxph(Surv(OS.month, OS.event == "1") ~ 1, data = OCCC.clinical)

fitf <- stepAIC(fit0, scope = formula(fit), direction = "both", k=2)

fitb <- stepAIC(fit, direction = "backward", k=2)

coxph(Surv(OS.month, OS.event == "1") ~ stage.group + TERT, data = OCCC.clinical) %>% 
  summary()

ggforest(model = coxph(Surv(OS.month, OS.event == "1") ~ stage.group + TERT, data = OCCC.clinical), 
         data = OCCC.clinical, fontsize = 1, main = "Forest plot of hazard ratio for OS")
