##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(dplyr)
library(GEOquery)
library(readxl)
library(tidyverse)
library(cancerclass)

# 12 genes
PSS_Sig.genes <- c("GIMAP7","PSMB8","CD27","CCR7","TAGAP","UQCR10","HCLS1","LCK","TNFAIP3","ISG15","GIMAP4","HLA-DRB1")


rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}


# single ROC plot
rocplot.single.V2 <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC = ",signif(auc, 2), " (95% CI: ", signif(ci.lower, 2), "-", signif(ci.upper, 2), ")", sep=""))
  }
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = "")) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_manual(labels = annotation, values = "#000000") +
    theme(
      plot.title = element_text(face="bold", size=15, color = "darkgreen"), 
      plot.subtitle = element_text(face="bold.italic", size=12, color = "orange"),
      plot.caption = element_text(face="plain", size=12, color = "blue"),
      axis.text.x = element_text(face="plain", size=14),
      axis.text.y = element_text(face="plain", size=14),
      axis.title.x = element_text(face="bold", size=18),
      axis.title.y = element_text(face="bold", size=18, angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.justification=c(1,0), 
      legend.position=c(0.95,0.15),
      legend.text = element_text(size = 14),
      legend.title=element_blank(),
      legend.key = element_blank()
    )+
    labs(title=title)
  return(p)
}


#################### GSE143153 ####################
GSE143153_eSet = GEOquery::getGEO("GSE143153")

exp_GSE143153 <- read.table("GSE143153.txt", header = T, sep = "\t")
exp_GSE143153 = exp_GSE143153[!duplicated(exp_GSE143153$GeneName),]  # 按照symbol列去重

GSE143153_PSS_Sig.genes <- exp_GSE143153[exp_GSE143153$GeneName %in% PSS_Sig.genes,]
rownames(GSE143153_PSS_Sig.genes) <- GSE143153_PSS_Sig.genes$GeneName
GSE143153_PSS_Sig.genes <- GSE143153_PSS_Sig.genes[,-1]
GSE143153_PSS_Sig.genes <- as.matrix(GSE143153_PSS_Sig.genes)
GSE143153_PSS_Sig.genes <- GSE143153_PSS_Sig.genes[rowSums(GSE143153_PSS_Sig.genes) > 1,]

dim(GSE143153_PSS_Sig.genes)
sort(rowSums(GSE143153_PSS_Sig.genes))
pdata_GSE143153 <- pData(GSE143153_eSet[[1]])
names(pdata_GSE143153)[39] <- "class"
GSE143153_phenoData <- new("AnnotatedDataFrame", data = pdata_GSE143153)

GSE143153_ExpSet <- ExpressionSet(assayData = as.matrix(GSE143153_PSS_Sig.genes), phenoData = GSE143153_phenoData)
pData(GSE143153_ExpSet)
dim(pData(GSE143153_ExpSet))
dim(exprs(GSE143153_ExpSet))
nrow(GSE143153_PSS_Sig.genes)
exprs(GSE143153_ExpSet)

nvalidation_GSE143153 <- nvalidate(GSE143153_ExpSet, ngenes = 2:(nrow(GSE143153_PSS_Sig.genes) + 1), ntrain = "balanced", method = "welch.test", dist = "cor")
pData(GSE143153_ExpSet)[["class"]]
predictor_GSE143153 <- fit(GSE143153_ExpSet, method = "welch.test") # must remove 0 or low expression genes < 10, otherwise, fit function does not work!!!
prediction_GSE143153 <- predict(predictor_GSE143153, GSE143153_ExpSet, "Primary SS", ngenes=nrow(GSE143153_PSS_Sig.genes), dist = "cor")
str(prediction_GSE143153)
prediction_GSE143153@predictor

out_GSE143153 <- as.factor(c("2","2","1","1","1","2","2","2","2","2","2","1","1","1","1","1","2","2","2","2",
                             "2","2","2","2","2","1","1","1","1","1","1","1"))

z_GSE143153 <- as.numeric(prediction_GSE143153@prediction[,'z'])

Test_GSE143153 <- cbind(out_GSE143153, z_GSE143153)
colnames(Test_GSE143153) <- c('grp','res')
Test_GSE143153 <- as.data.frame(Test_GSE143153)

p.GSE143153 <- rocplot.single.V2(Test_GSE143153$grp, Test_GSE143153$res, title = "") 


#################### GSE84844 ####################
GSE84844_eSet = GEOquery::getGEO("GSE84844")

#通过exprs()函数获取表达矩阵并校正
exp_GSE84844 <- exprs(GSE84844_eSet[[1]])

# 通过pData函数获取分组信息
pdata_GSE84844 <- pData(GSE84844_eSet[[1]])
table(pdata_GSE84844$characteristics_ch1)
names(pdata_GSE84844)[10] <- "class"

GPL_GSE84844 <- GSE84844_eSet[[1]]@annotation  # 平台信息——提取芯片平台编号
library(hgu133plus2.db) # 参考：http://www.bio-info-trainee.com/1399.html
ls("package:hgu133plus2.db") 

ids_GSE84844 <- toTable(hgu133plus2SYMBOL)  # 提取 
colnames(ids_GSE84844) = c("probe_id", "symbol")
exp_GSE84844 = as.data.frame(exp_GSE84844) 
exp_GSE84844$probe_id = rownames(exp_GSE84844)  # 将行名变为列名为probe_id的一列

# exp是原来的表达矩阵
exp2_GSE84844 = merge(exp_GSE84844, ids_GSE84844, by.x="probe_id", by.y="probe_id")  # 合并数据
exp2_GSE84844 = exp2_GSE84844[!duplicated(exp2_GSE84844$symbol),]  # 按照symbol列去重

# 数据框probe_exp的行名变成symbol
exp2_GSE84844 = exp2_GSE84844[,c(62, 2:61)]

GSE84844_PSS_Sig.genes <- exp2_GSE84844[exp2_GSE84844$symbol %in% PSS_Sig.genes,]
rownames(GSE84844_PSS_Sig.genes) <- GSE84844_PSS_Sig.genes$symbol
GSE84844_PSS_Sig.genes <- GSE84844_PSS_Sig.genes[,-1]
GSE84844_PSS_Sig.genes <- as.matrix(GSE84844_PSS_Sig.genes)
GSE84844_PSS_Sig.genes <- GSE84844_PSS_Sig.genes[rowSums(GSE84844_PSS_Sig.genes) > 1,]

dim(GSE84844_PSS_Sig.genes)
sort(rowSums(GSE84844_PSS_Sig.genes))
GSE84844_phenoData <- new("AnnotatedDataFrame", data = pdata_GSE84844)

all(pdata_GSE84844$sample[1:30]==colnames(GSE84844_PSS_Sig.genes)[1:30])
all(pdata_GSE84844$sample[31:60]==colnames(GSE84844_PSS_Sig.genes)[31:60])

GSE84844_ExpSet <- ExpressionSet(assayData = as.matrix(GSE84844_PSS_Sig.genes), phenoData = GSE84844_phenoData)
pData(GSE84844_ExpSet)
dim(pData(GSE84844_ExpSet))
dim(exprs(GSE84844_ExpSet))
nrow(GSE84844_PSS_Sig.genes)
exprs(GSE84844_ExpSet)

nvalidation_GSE84844 <- nvalidate(GSE84844_ExpSet, ngenes = 2:(nrow(GSE84844_PSS_Sig.genes) + 1), ntrain = "balanced", method = "welch.test", dist = "cor")
pData(GSE84844_ExpSet)[["class"]]
predictor_GSE84844 <- fit(GSE84844_ExpSet, method = "welch.test") # must remove 0 or low expression genes < 10, otherwise, fit function does not work!!!
prediction_GSE84844 <- predict(predictor_GSE84844, GSE84844_ExpSet, "disease: primary Sjogren's syndrome", ngenes=nrow(GSE84844_PSS_Sig.genes), dist = "cor")
str(prediction_GSE84844)
prediction_GSE84844@predictor

out_GSE84844 <- as.factor(rep(c(1, 2),c(30, 30)))
out_GSE84844 

z_GSE84844 <- as.numeric(prediction_GSE84844@prediction[,'z'])
Test_GSE84844 <- cbind(out_GSE84844, z_GSE84844)
colnames(Test_GSE84844) <- c('grp','res')
Test_GSE84844 <- as.data.frame(Test_GSE84844)

p.GSE84844 <- rocplot.single.V2(Test_GSE84844$grp, Test_GSE84844$res, title = "") 

