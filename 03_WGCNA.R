#WGCNA
#0 数据准备
          setwd(dir = "D:/000/bulk-WGCNA/data")
          load("02合并后counts及samples.RData")
          load("03DEGs.RData")
          library(WGCNA)
          datExpr <- counts
          datTraits = samples
          rownames(datTraits) <- datTraits$GSE
          
          gene_WGCNA <- c(rownames(OBvsMsc_gene),rownames(ADvsMsc_gene))
          datExpr <- datExpr[rownames(datExpr) %in% gene_WGCNA,]
          datExpr = t(datExpr)   ## 转置
          gsg = goodSamplesGenes(datExpr, verbose = 3)    #检查缺失值
          gsg$allOK
          datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
          
#1 确定最佳beta值
          powers = c(c(1:10), seq(from = 12, to=20, by=2))       #幂指数范围
          sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)  #选取1到20为候选的power值，使用pickSoftThreshold()函数对相应的power值进行计算
          sft$powerEstimate     
          
          #只是为了看看，可做可不做的
          par(mfrow = c(1,2))
          cex1 = 0.9
          ###拟合指数(Scale-free topology fit index)与power值(the soft-thresholding power)散点图
          plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
               main = paste("Scale independence"));
          text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               labels=powers,cex=cex1,col="red");
          abline(h=0.50,col="red") #可以修改
          ###平均连通性（Mean connectivity）与power值散点图
          plot(sft$fitIndices[,1], sft$fitIndices[,5],
               xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
               main = paste("Mean connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
          #根据上图把power值设为10吧
#2 一步法构建共表达矩阵
          #把输入的表达矩阵的几千个基因归类成了几十个模块。
          #大体思路：
          #①计算基因间的邻接性，根据邻接性计算基因间的相似性，然后推出基因间的相异性系数，并据此得到基因间的系统聚类树。
          #②然后按照混合动态剪切树的标准，设置每个基因模块最少的基因数目为30。
          #③依次计算每个模块的特征向量值，然后对模块进行聚类分析，将距离较近的模块合并为新的模块。
          cor <- WGCNA::cor    #https://blog.csdn.net/liyunfan00/article/details/91686840
          net = blockwiseModules(
            datExpr,
            power = 10,   #即sft$powerEstimate
            maxBlockSize = 5000,
            TOMType = "unsigned", minModuleSize = 30,
            reassignThreshold = 0, mergeCutHeight = 0.25,
            numericLabels = TRUE, pamRespectsDendro = FALSE,
            saveTOMs = F, 
            verbose = 3
          )
          cor<-stats::cor
          table(net$colors) 
          team <- data.frame(names(net$colors),net$colors);write.csv(team,file = "team.csv")
          
#3 模块可视化
          # Convert labels to colors for plotting
          mergedColors = labels2colors(net$colors)     #灰色默认是无法归类于任何模块的那些基因
          table(mergedColors)
          # Plot the dendrogram(树状图) and the module colors underneath
          plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                              "Module colors",
                              dendroLabels = FALSE, hang = 0.03,
                              addGuide = TRUE, guideHang = 0.05)
          #plotDendroAndColors函数接受一个聚类的对象，以及该对象里面包含的所有个体所对应的颜色。
          
#4 模块与性状的关系
          ## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
          table(datTraits$group2)
          
          nGenes = ncol(datExpr)
          nSamples = nrow(datExpr)
          design=model.matrix(~0+ datTraits$group2)
          colnames(design)=levels(datTraits$group2)
          head(design)
          moduleColors <- labels2colors(net$colors)
          
          # Recalculate MEs with color labels
          MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
          MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
          moduleTraitCor = cor(MEs, design , use = "p");
          moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
          
          sizeGrWindow(10,6)
          # Will display correlations and their p-values
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                             signif(moduleTraitPvalue, 1), ")", sep = "");
          dim(textMatrix) = dim(moduleTraitCor)
          par(mar = c(6, 8.5, 3, 3));
          # Display the correlation values within a heatmap plot,这个图就是把moduleTraitCor这个矩阵给用热图可视化一下
          labeledHeatmap(Matrix = moduleTraitCor,
                         xLabels = colnames(design),
                         yLabels = names(MEs),
                         ySymbols = names(MEs),
                         colorLabels = FALSE,
                         colors = greenWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 1.5,
                         zlim = c(-1,1),
                         main = paste("Module-trait relationships"))
        
#5 感兴趣性状的模块的具体基因分析
          ##模块与基因的相关性矩阵 
          #names (colors) of the modules
          modNames = substring(names(MEs), 3)
          geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
          ## 算出每个模块跟基因的皮尔森相关系数矩阵
          ## MEs是每个模块在每个样本里面的值
          ## datExpr是每个基因在每个样本的表达量
          MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
          names(geneModuleMembership) = paste("MM", modNames, sep="");
          names(MMPvalue) = paste("p.MM", modNames, sep="");
          
          #性状与基因的相关性矩阵
          
          ## OB_3d与brown块的关系
          ## 这里把是否属于 OB_3d 表型这个变量用0,1进行数值化。
          OB_3d = as.data.frame(design[,2]);
          names(OB_3d) = "OB_3d"
          geneTraitSignificance = as.data.frame(cor(datExpr, OB_3d, use = "p"));
          GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
          names(geneTraitSignificance) = paste("GS.", names(OB_3d), sep="");
          names(GSPvalue) = paste("p.GS.", names(OB_3d), sep="");
          #把两个相关性矩阵联合起来,指定感兴趣模块进行分析
          module = "brown"
          column = match(module, modNames);
          moduleGenes = moduleColors==module;
          sizeGrWindow(7, 7);
          par(mfrow = c(1,1));
          verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                             abs(geneTraitSignificance[moduleGenes, 1]),
                             xlab = paste("Module Membership in", module, "module"),
                             ylab = "Gene significance for AD_terminal",
                             main = paste("Module membership vs. gene significance\n"),
                             cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
          
          ## AD_3d与blue块的关系
          ## 这里把是否属于 AD_terminal 表型这个变量用0,1进行数值化。
          AD_3d = as.data.frame(design[,3]);
          names(AD_3d) = "AD_3d"
          geneTraitSignificance = as.data.frame(cor(datExpr, AD_3d, use = "p"));
          GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
          names(geneTraitSignificance) = paste("GS.", names(AD_3d), sep="");
          names(GSPvalue) = paste("p.GS.", names(AD_3d), sep="");
          #把两个相关性矩阵联合起来,指定感兴趣模块进行分析
          module = "blue"
          column = match(module, modNames);
          moduleGenes = moduleColors==module;
          sizeGrWindow(7, 7);
          par(mfrow = c(1,1));
          verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                             abs(geneTraitSignificance[moduleGenes, 1]),
                             xlab = paste("Module Membership in", module, "module"),
                             ylab = "Gene significance for OS_terminal",
                             main = paste("Module membership vs. gene significance\n"),
                             cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#6 感兴趣模块的hub基因箱式图
          gene_blue <- names(net$colors[net$colors==2])
          gene_brown <- names(net$colors[net$colors==3])
          write.csv(gene_blue,file = "gene_blue.csv")
          write.csv(gene_brown,file = "gene_brown.csv")
          
#7 富集分析
       #GO分析
          library(clusterProfiler)
          library(enrichplot)
          #OBvsMsc blue模块
          gene_go <- AnnotationDbi::mapIds(Homo.sapiens,gene_blue, column= "ENTREZID",keytype = "SYMBOL")
          ego_MF_OB <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          ego_BP_OB <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          ego_CC_OB <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          write.csv(ego_CC_OB,file = "ego_CC_OB.csv")
          write.csv(ego_BP_OB,file = "ego_BP_OB.csv")
          write.csv(ego_MF_OB,file = "ego_MF_OB.csv")
          #每组挑选10个感兴趣的通路
          slimgo <- data.frame(Ontology=c(rep("Biological process",10),
                                          rep("Molecular function",10),
                                          rep("Cellular component",10)),
                               rbind(ego_BP_OB@result[1:10,],ego_MF_OB[1:10,],ego_CC_OB[1:10,]))
          #首先需要将Description转换为因子
          slimgo$Description=factor(slimgo$Description,levels=slimgo$Description)
          colnames(slimgo)
          pp <- ggplot(data = slimgo, mapping = aes(x=Description,y=Count,fill=Ontology))
          pp
          pp+geom_bar(stat="identity")
          pp+geom_bar(stat="identity")+coord_flip()
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))+guides(fill=FALSE)
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))+guides(fill=FALSE)+theme_bw()
      #ADvsMsc brown
          gene_go <- AnnotationDbi::mapIds(Homo.sapiens,gene_brown, column= "ENTREZID",keytype = "SYMBOL")
          ego_MF_AD <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          ego_BP_AD <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          ego_CC_AD <- enrichGO(gene          = gene_go,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
          write.csv(ego_CC_AD,file = "ego_CC_AD.csv")
          write.csv(ego_BP_AD,file = "ego_BP_AD.csv")
          write.csv(ego_MF_AD,file = "ego_MF_AD.csv")
          #每组挑选10个感兴趣的通路
          slimgo <- data.frame(Ontology=c(rep("Biological process",10),
                                          rep("Molecular function",10),
                                          rep("Cellular component",10)),
                               rbind(ego_BP_AD@result[1:10,],ego_MF_AD[1:10,],ego_CC_AD[1:10,]))
          #首先需要将Description转换为因子
          slimgo$Description=factor(slimgo$Description,levels=slimgo$Description)
          colnames(slimgo)
          pp <- ggplot(data = slimgo, mapping = aes(x=Description,y=Count,fill=Ontology))
          pp
          pp+geom_bar(stat="identity")
          pp+geom_bar(stat="identity")+coord_flip()
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))+guides(fill=FALSE)
          pp+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Description)))+guides(fill=FALSE)+theme_bw()
          
          
#8 hub基因
          # Recalculate topological overlap
          TOM = TOMsimilarityFromExpr(datExpr, power = 8); 
          # Select module
          module = "brown";
          # Select module probes
          probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
          inModule = (moduleColors==module);
          modProbes = probes[inModule]; 
          ## 也是提取指定模块的基因名
          # Select the corresponding Topological Overlap
          modTOM = TOM[inModule, inModule];
          dimnames(modTOM) = list(modProbes, modProbes)
          ## 模块对应的基因关系矩阵 
          
          cyt = exportNetworkToCytoscape(
            modTOM,
            edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
            nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
            weighted = TRUE,
            threshold = 0.02,
            nodeNames = modProbes, 
            nodeAttr = moduleColors[inModule]
          );
          
          
          #通过string，MCODE找blue模块、brown模块中的hub基因
          
          
          
          #hub基因表达箱式图
          gene_os_hub <- c("UBE2C","TOP2A","NUSAP1","MELK","KIF4A","CDC20","CCNB2","CCNA2","BIRC5")
          exp_os_hub <- counts[gene_os_hub,]
          score <- samples
          
          library(reshape2)
          library(ggpubr)
          
          sameGene=intersect(row.names(exp), gene)
          geneExp=as.data.frame(t(exp[sameGene,]))
          sameSample=intersect(rownames(geneExp), row.names(score))
          geneExp=geneExp[sameSample,]
          score=score[sameSample,]
          rt=cbind(geneExp, ICIscore=score[,"diff"])
          
          data=melt(rt, id.vars=c("ICIscore"))
          colnames(data)=c("ICIscore", "Gene", "Expression")
          data$ICIscore=factor(data$ICIscore, levels=c("AD", "OS"))
          
          bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
          bioCol=bioCol[1:length(levels(factor(rt[,"ICIscore"])))]
          p=ggboxplot(data, x="Gene", y="Expression", color="ICIscore", 
                      ylab="10",
                      xlab="",
                      legend.title="ICI score",
                      palette=bioCol)
          p=p+rotate_x_text(50)
          pdf(file="boxplot.pdf",width=30,height=6)
          p+stat_compare_means(aes(group=ICIscore),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
          dev.off()
          