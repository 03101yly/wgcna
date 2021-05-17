#考虑到已经有论文对GSE113253的数据进行分析，我更改为分析成骨成脂早期数据，
#分析原代MSC、成骨分化三天、成脂分化三天的数据，合并GSE113253和GSE80614的数据

#GSE113253数据获取、整理、差异分析、可视化
          setwd(dir = "D:/000/bulk-WGCNA/data")          #设置工作目录
          rm(list=ls())
          suppressPackageStartupMessages({
            library(GEOquery)
            library(tidyr)
            library(dplyr)
            library(limma)
            library(Glimma)
            library(edgeR)
            library(Homo.sapiens)
            library(RColorBrewer)
            library(EnhancedVolcano)
            library(pheatmap)
            library(ggplot2)
            library(WGCNA)
            library(ggbiplot)
            library(AnnotationHub)
          })
#0 数据准备（导入数据，整理列名，转换为TPM（可选），整理样本信息，提取所需列，保留编码蛋白基因，删除低表达基因，标准化，质量控制）
#导入数据GSE113253的数据
          GSE113253_seq <- read.delim("GSE113253_GeneE.txt",header = T,sep = "\t",row.names = 1,fill = T)
          GSE113253_seq <- tidyr::separate(GSE113253_seq,Annotation.Divergence,c("symbol"),sep = "\\|")
          rownames(GSE113253_seq) <- GSE113253_seq$symbol     #前提是symbol中没有重复
          efflen <- GSE113253_seq$Length    #转录本的长度，用于counts转TPM，本处好像用不到，因为WGCNA和edgeR都可以采用counts
     #整理列名
          GSE113253_seq <- GSE113253_seq[,3:35]
          col <- c(paste0("RNA_14dOb_BM_rep",1:3),paste0("RNA_7dOb_BM_rep",1:3),paste0("RNA_3dOb_BM_rep",1:3),
                   paste0("RNA_1dOb_BM_rep",1:3),paste0("RNA_4hOb_BM_rep",1:3),paste0("RNA_Msc_BM_rep",1:3),
                   paste0("RNA_4hAd_BM_rep",1:3),paste0("RNA_1dAd_BM_rep",1:3),paste0("RNA_3dAd_BM_rep",1:3),
                   paste0("RNA_7dAd_BM_rep",1:3),paste0("RNA_14dAd_BM_rep",1:3))
          GSE113253_seq <-GSE113253_seq[,col]   #按照GSE号进行排序
          colname <- paste0("GSM",3100994:3101026)
          colnames(GSE113253_seq) <- colname
      #counts转TPM
          effLen <- data.frame(a=efflen,b=efflen,c=efflen,d=efflen,e=efflen,f=efflen,g=efflen,h=efflen,i=efflen,j=efflen,
                            a1=efflen,b1=efflen,c1=efflen,d1=efflen,e1=efflen,f1=efflen,g1=efflen,h1=efflen,i1=efflen,j1=efflen,
                            a2=efflen,b2=efflen,c2=efflen,d2=efflen,e2=efflen,f2=efflen,g2=efflen,h2=efflen,i2=efflen,j2=efflen,
                            a3=efflen,b3=efflen,c3=efflen)
          Counts2TPM <- function(counts, effLen){
              rate <- log(counts) - log(effLen)
              denom <- log(sum(exp(rate)))
              exp(rate - denom + log(1e6))
          }
          GSE113253_TPM <- Counts2TPM(GSE113253_seq,effLen)
    
     #整理基因数据，根据ensembl中的数据，只保留编码蛋白的基因,删除symbol重复的基因
          edb <- AnnotationHub()[["AH89426"]] # human, Ensembl v103.
          GSE113253_rowData <- data.frame(symbol = rownames(GSE113253_TPM))
          library(HGNChelper)
          GSE113253_rowData[,2:4] <- checkGeneSymbols(GSE113253_rowData$symbol)
          GSE113253_rowData <- GSE113253_rowData[!is.na(GSE113253_rowData$Suggested.Symbol),]
          GSE113253_rowData$ensemblID <- mapIds(edb, keys=GSE113253_rowData$Suggested.Symbol, keytype="SYMBOL", column="GENEID")
          GSE113253_rowData$cod <- mapIds(edb, keys=GSE113253_rowData$Suggested.Symbol, keytype="SYMBOL", column="TXBIOTYPE")
          GSE113253_rowData <- GSE113253_rowData[!GSE113253_rowData$cod %in% c("lncRNA","misc_RNA","miRNA","scaRNA","scRNA","snoRNA","snRNA","TEC"),]
          GSE113253_rowData <- na.omit(GSE113253_rowData)
          
          table(duplicated(GSE113253_rowData$Suggested.Symbol))
          GSE113253_TPM <- GSE113253_TPM[rownames(GSE113253_TPM) %in% GSE113253_rowData$symbol,]    #将表达矩阵中不存在于对应关系中的probe_id行删除
          GSE113253_rowData <- GSE113253_rowData[match(rownames(GSE113253_TPM),GSE113253_rowData$symbol),]          #将对应关系中的probe_id按照表达矩阵中的顺序排列
          GSE113253_rowData$median <- apply(GSE113253_TPM,1,median)        #GSE113253_rowData新建median这一列，同时对GSE113253_TPM这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
          GSE113253_rowData <- GSE113253_rowData[order(GSE113253_rowData$Suggested.Symbol,GSE113253_rowData$median,decreasing = T),]   #对GSE113253_rowData$symbol按照GSE113253_rowData$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的GSE113253_rowData
          GSE113253_rowData <- GSE113253_rowData[!duplicated(GSE113253_rowData$Suggested.Symbol),]   #将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果
          GSE113253_TPM <- GSE113253_TPM[rownames(GSE113253_TPM) %in% GSE113253_rowData$symbol,] #新的GSE113253_rowData取出probe_id这一列，将GSE113253_TPM按照取出的这一列中的每一行组成一个新的GSE113253_TPM
          GSE113253_rowData <- GSE113253_rowData[match(rownames(GSE113253_TPM),GSE113253_rowData$symbol),]          #将GSE113253_rowData中的probe_id按照表达矩阵中的顺序排列
          identical(GSE113253_rowData$symbol,rownames(GSE113253_TPM))                #判断一下两个名字和顺序是否一样
          rownames(GSE113253_TPM) <- GSE113253_rowData[,4]                                #将rowData的symbol列更改为表达矩阵的行名

     #删除低表达基因
          palmieri_medians <- apply(GSE113253_TPM, 1, median)
          hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, 
                           main = "Histogram of the median intensities", 
                           border = "antiquewhite4",
                           xlab = "Median intensities")
          threshold <- apply(GSE113253_TPM, 1,function(x){sum(x > 0) >= 1})    #这里的阈值，我设置的极宽
          table(threshold)
          GSE113253_TPM <- subset(GSE113253_TPM, threshold)  # 提取过滤转录本
     #整理样本信息
          GSE113253_metadata <- data.frame(GSE=colname,group1=col,group2=col)
          GSE113253_metadata <- tidyr::separate(GSE113253_metadata,group2,c("RNA","group"),sep = "_")
          GSE113253_metadata$RNA <- "seq"
          colnames(GSE113253_metadata)[3] <- "batch"
     #提取所需列
          GSE113253_samples <- GSE113253_metadata[GSE113253_metadata$group %in% c('3dOb','Msc',"3dAd"),]
          GSE113253_counts <- GSE113253_TPM[,colnames(GSE113253_TPM) %in% GSE113253_samples$GSE]
          
          #至此，提取了GSE113253的数据，分别为GSE113253_counts、GSE113253_samples、GSE113253_rowData
          
#提取GSE80614的数据
          gse <- getGEO('GSE80614',destdir = ".")
          GSE80614_assay <- exprs(gse[[1]])
          GSE80614_colData <- pData(phenoData(gse[[1]]))
          GSE80614_rowData <- pData(featureData(gse[[1]]))
          
          #整理基因数据，根据ensembl中的数据，只保留编码蛋白的基因,删除symbol重复的基因
          GSE80614_rowData <- GSE80614_rowData[,c(1,7,8)]
          library(HGNChelper)
          GSE80614_rowData[,4:6] <- checkGeneSymbols(GSE80614_rowData$ILMN_Gene)
          GSE80614_rowData <- GSE80614_rowData[!is.na(GSE80614_rowData$Suggested.Symbol),]
          library(AnnotationHub)
          edb <- AnnotationHub()[["AH89426"]] # human, Ensembl v103.
          GSE80614_rowData$ensemblID <- mapIds(edb, keys=GSE80614_rowData$Suggested.Symbol, keytype="SYMBOL", column="GENEID")
          GSE80614_rowData$cod <- mapIds(edb, keys=GSE80614_rowData$Suggested.Symbol, keytype="SYMBOL", column="TXBIOTYPE")
          GSE80614_rowData <- GSE80614_rowData[!GSE80614_rowData$cod %in% c("lncRNA","misc_RNA","miRNA","scaRNA","scRNA","snoRNA","snRNA","TEC"),]
          GSE80614_rowData <- na.omit(GSE80614_rowData)
          
          table(duplicated(GSE80614_rowData$ILMN_Gene))
          GSE80614_assay <- GSE80614_assay[rownames(GSE80614_assay) %in% rownames(GSE80614_rowData),]    #将表达矩阵中不存在于对应关系中的probe_id行删除
          GSE80614_rowData <- GSE80614_rowData[match(rownames(GSE80614_assay),rownames(GSE80614_rowData)),]          #将对应关系中的probe_id按照表达矩阵中的顺序排列
          GSE80614_rowData$median <- apply(GSE80614_assay,1,median)        #GSE80614_rowData新建median这一列，同时对GSE80614_assay这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
          GSE80614_rowData <- GSE80614_rowData[order(GSE80614_rowData$ILMN_Gene,GSE80614_rowData$median,decreasing = T),]   #对GSE80614_rowData$symbol按照GSE80614_rowData$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的GSE80614_rowData
          GSE80614_rowData <- GSE80614_rowData[!duplicated(GSE80614_rowData$ILMN_Gene),]   #将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果
          GSE80614_assay <- GSE80614_assay[rownames(GSE80614_assay) %in% rownames(GSE80614_rowData),] #新的GSE80614_rowData取出probe_id这一列，将GSE80614_assay按照取出的这一列中的每一行组成一个新的GSE80614_assay
          GSE80614_rowData <- GSE80614_rowData[match(rownames(GSE80614_assay),rownames(GSE80614_rowData)),]          #将GSE80614_rowData中的probe_id按照表达矩阵中的顺序排列
          identical(GSE80614_rowData$ID,rownames(GSE80614_assay))                #判断一下两个名字和顺序是否一样
          rownames(GSE80614_assay) <- GSE80614_rowData[,6]                                #将rowData的symbol列更改为表达矩阵的行名
          
          #整理样本信息
          GSE80614_metadata <- GSE80614_colData[,c(1,8)]
          GSE80614_metadata <- GSE80614_metadata[match(colnames(GSE80614_assay),rownames(GSE80614_colData)),]
          GSE80614_metadata$GSE <- rownames(GSE80614_metadata)
          GSE80614_metadata$batch <- "array"
          colnames(GSE80614_metadata)[2] <-"group1" 
          GSE80614_metadata <- separate(GSE80614_metadata,title,c("group"),sep="-")
          GSE80614_metadata <- GSE80614_metadata[,c(3,2,4,1)]
          
          ##提取所需列
          GSE80614_samples <- GSE80614_metadata[GSE80614_metadata$group %in% c('AD_0hr','AD_3d',"OS_0hr","OS_3d"),]
          GSE80614_counts <- GSE80614_assay[,colnames(GSE80614_assay) %in% GSE80614_samples$GSE]
          
          #至此，提取了GSE80614的数据，分别为GSE80614_counts、GSE80614_samples、GSE80614_rowData
          save(GSE113253_counts,GSE113253_samples,GSE113253_rowData,GSE80614_counts,GSE80614_samples,GSE80614_rowData,file = "01整理后counts及samples.RData")
          load("01整理后counts及samples.RData")
#合并两个数据集
          #对于测序的counts数据，批次矫正采用sva包的ComBat_seq,其他的采用ComBat
          library(limma)
          library(sva)
          library(EDASeq)
          library(RColorBrewer)
          #获取交集基因
          geneList <- list(seq=rownames(GSE113253_counts),array=rownames(GSE80614_counts))
          intersectGenes=Reduce(intersect, geneList)
          #合并数据集
          GSE113253_counts <- normalizeBetweenArrays(GSE113253_counts) #不知道是否需要
          allTab <- cbind(GSE113253_counts[intersectGenes,],GSE80614_counts[intersectGenes,])
          #标注批次
          batchType <- c(rep("seq",ncol(GSE113253_counts)), rep("array",ncol(GSE80614_counts)))
          #矫正前查看数据
          #创建SeqExpressionSet对象
          set <- newSeqExpressionSet(as.matrix(allTab),
                                     phenoData = data.frame(batchType,row.names = colnames(allTab)))
          set@phenoData@data
          set
          # 绘制PCA及箱线图检查数据
          col <- brewer.pal(3,"Set2")
          plotRLE(set,outline=FALSE,ylim=c(-4,4),col=col[1:2])
          plotPCA(set,col=col[1:2],cex=0.6)
          
          #数据矫正
          counts <- ComBat(allTab, batchType, par.prior = TRUE)
          samples <- rbind(GSE113253_samples,GSE80614_samples)
          
# 绘制PCA及箱线图检查数据
          set <- newSeqExpressionSet(as.matrix(counts),
                                     phenoData = data.frame(batchType,row.names = colnames(allTab)))
          set@phenoData@data
          set
          col <- brewer.pal(3,"Set2")
          plotRLE(set,outline=FALSE,ylim=c(-0.1,0.1),col=c(rep(col[1],9),rep(col[2],12)))
          
          #绘制PCA图
          samples$group2 <- c(rep("OBs",3),rep("MSCs",3),rep("ADs",3),rep("MSCs",3),rep("ADs",3),
                              rep("MSCs",3),rep("OBs",3))
          samples$group2 <- factor(samples$group2,
                                   levels = c("MSCs","OBs","ADs"))
          counts_PCA <- t(counts)
          #画图
          result <- prcomp(counts_PCA,scale = T)
          ggscreeplot(result)    #碎石图
          ggbiplot(result,obs.scale = 1,var.scale = 1,
                   groups = samples$group2,ellipse = T,circle = T,var.axes = F)+
            scale_color_brewer(palette = "Set1")+
            theme(legend.direction = "horizontal",legend.position = "top")
          
          save.image(file="数据合并.RData")
          save(counts,samples,file = "02合并后counts及samples.RData")
          load("02合并后counts及samples.RData")
          