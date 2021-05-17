#差异分析及GO、KEGG、GSEA，提取差异基因进行WGCNA分析
         setwd(dir = "D:/000/bulk-WGCNA/data")          #设置工作目录
         rm(list=ls())
         suppressPackageStartupMessages({
           library(edgeR)
           library(EnhancedVolcano)
           library(clusterProfiler)
           library(enrichplot)
           library(org.Hs.eg.db)
           library(Homo.sapiens)
         })
         load("02合并后counts及samples.RData")
#1 差异分析
    #创建DGElist
          x <- DGEList(counts = counts,samples = samples,genes = rownames(counts))
    #归一化基因表达    
          x <- calcNormFactors(x, method = "TMM")
          x$samples$norm.factors
    #创建设计矩阵和对比
          design <- model.matrix(~0+x$samples$group2) 
          colnames(design) <- c("MSC","OB","AD")
          design
          #成对比较
          contr.matrix <- makeContrasts(OBvsMSC = OB-MSC, 
                                        ADvsMSC = AD-MSC,
                                        levels = colnames(design))
          contr.matrix
    #从表达计数数据中删除异方差
          v <- voom(x, design, plot=TRUE)
          v
    #拟合线性模型
          vfit <- lmFit(v, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit)
          plotSA(efit, main="Final model: Mean-variance trend")
          
    #检查DE基因数量
          #在多个对比中皆差异表达的基因可以从decideTests的结果中提取，
          #其中的0代表不差异表达的基因，1代表上调的基因，-1代表下调的基因。
          dt <- decideTests(efit,p.value=0.05,lfc = 0.2)
          summary(dt)        #显著性的判断使用默认的校正p值阈值，即5%。
          vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
          
          OBvsMsc <- topTable(efit, coef=1, n=Inf)
          ADvsMsc <- topTable(efit, coef=2, n=Inf)
          
          OBvsMsc_UP <- OBvsMsc[OBvsMsc$adj.P.Val<0.05 & OBvsMsc$logFC>0.2,]
          OBvsMsc_DOWN <- OBvsMsc[OBvsMsc$adj.P.Val<0.05 & OBvsMsc$logFC<(-0.2),]
          ADvsMsc_UP <- ADvsMsc[ADvsMsc$adj.P.Val<0.05 & ADvsMsc$logFC>0.2,]
          ADvsMsc_DOWN <- ADvsMsc[ADvsMsc$adj.P.Val<0.05 & ADvsMsc$logFC<(-0.2),]
          
          OBvsMsc_gene <- rbind(OBvsMsc_UP,OBvsMsc_DOWN)
          ADvsMsc_gene <- rbind(ADvsMsc_UP,ADvsMsc_DOWN)
          
          gene_intersect_1 <- intersect(OBvsMsc_UP$genes,ADvsMsc_DOWN$genes)
          gene_intersect_2 <- intersect(OBvsMsc_DOWN$genes,ADvsMsc_UP$genes)
          gene_intersect <- c(gene_intersect_1,gene_intersect_2)
          gene_intersect
          
          save(dt,OBvsMsc,ADvsMsc,OBvsMsc_UP,OBvsMsc_DOWN,ADvsMsc_UP,ADvsMsc_DOWN,OBvsMsc_gene,
               ADvsMsc_gene,gene_intersect_1,gene_intersect_2,gene_intersect,file = "03DEGs.RData")
          
#火山图
          EnhancedVolcano(OBvsMsc,
                          lab = rownames(OBvsMsc),
                          x = 'logFC',
                          y = 'P.Value',
                          pCutoff = 10e-10,
                          FCcutoff = 0.5,
                          xlim = c(-1, 1.5),
                          ylim = c(0, -log10(10e-20)),
                          pointSize = 1,
                          labSize = 3.5,
                          title = 'OBvsMsc',
                          subtitle = 'Differential expression',
                          caption = 'FC cutoff, 0.5; p-value cutoff, 10e-10',
                          legendPosition = "top",
                          legendLabSize = 14,
                          col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                          colAlpha = 0.9,
                          drawConnectors = TRUE,
                          hline = c(10e-10),
                          widthConnectors = 0.5)
          
          EnhancedVolcano(ADvsMsc,
                          lab = rownames(ADvsMsc),
                          x = 'logFC',
                          y = 'P.Value',
                          pCutoff = 10e-12,
                          FCcutoff = 0.5,
                          xlim = c(-1, 2),
                          ylim = c(0, -log10(10e-20)),
                          pointSize = 2,
                          labSize = 3.5,
                          title = 'ADvsMsc',
                          subtitle = 'Differential expression',
                          caption = 'FC cutoff, 0.6; p-value cutoff, 10e-12',
                          legendPosition = "top",
                          legendLabSize = 14,
                          col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                          colAlpha = 0.9,
                          drawConnectors = TRUE,
                          hline = c(10e-12),
                          widthConnectors = 0.5)
#热图      
      #OBvsMsc
          data1 <- rownames(head(OBvsMsc[order(OBvsMsc$logFC, decreasing = T),],10))
          data2 <- rownames(head(OBvsMsc[order(OBvsMsc$logFC, decreasing = F),],10))
          data3 <- c(data1,data2)
          counts_heatmap_1 <- counts[data3,samples$group2 %in% c("OBs","MSCs")]
          counts_heatmap_1 <-counts_heatmap_1[,c(13:15,1:12)]
          annotation_col <- data.frame(treatment=c(rep("OBs",6),rep("MSCs",9)))
          rownames(annotation_col) <- colnames(counts_heatmap_1)
          pheatmap(counts_heatmap_1,
                   scale = "row",
                   border="white", # 设置边框为白色
                   show_rownames = T, #去掉横、纵坐标id
                   show_colnames = T,
                   annotation_col = annotation_col,
                   treeheight_col = 20,
                   treeheight_row = 20,
                   fontsize_row = 10,
                   fontsize_col = 12)
      #ADvsMsc    
          data1 <- rownames(head(ADvsMsc[order(ADvsMsc$logFC, decreasing = T),],10))
          data2 <- rownames(head(ADvsMsc[order(ADvsMsc$logFC, decreasing = F),],10))
          data3 <- c(data1,data2)
          counts_heatmap_1 <- counts[data3,samples$group2 %in% c("ADs","MSCs")]
          annotation_col <- samples[samples$group2 %in% c("ADs","MSCs"),]
          annotation_col <- data.frame(annotation_col$group2) 
          rownames(annotation_col) <- colnames(counts_heatmap_1)
          colnames(annotation_col) <- "treatment"
          pheatmap(counts_heatmap_1,
                   scale = "row",
                   border="white", # 设置边框为白色
                   show_rownames = T, #去掉横、纵坐标id
                   show_colnames = T,
                   annotation_col = annotation_col,
                   fontsize_row = 10,
                   fontsize_col = 12,
                   cluster_cols = T, treeheight_col = 20, # 分别设置横、纵向聚类树高
                   cluster_rows = T, treeheight_row = 20)
#基因交集
          library(ggsci)
          mycol <- pal_nejm()(8)
          venn.diagram(x = list(OBvsMsc = c(OBvsMsc_UP$genes,OBvsMsc_DOWN$genes),
                                ADvsMsc = c(ADvsMsc_UP$genes,ADvsMsc_DOWN$genes)),
                       filename = "inter1.tiff",
                       scaled=F,lwd = 2,fill = mycol[1:2],alpha = 0.75,
                       label.col = "white",cex = 3,fontfamily = "serif",fontface = "bold",
                       cat.col = mycol[1:2],cat.cex = 3,cat.fontfamily = "serif",cat.fontface = "bold",cat.dist = c(0.2, 0.2));
          venn.diagram(x = list(OBvsMsc_UP = OBvsMsc_UP$genes,
                                ADvsMsc_DOWN = ADvsMsc_DOWN$genes),
                       filename = "inter2.tiff",
                       scaled=F,lwd = 2,fill = mycol[3:4],alpha = 0.75,
                       label.col = "white",cex = 3,fontfamily = "serif",fontface = "bold",
                       cat.col = mycol[3:4],cat.cex = 3,cat.fontfamily = "serif",cat.fontface = "bold",cat.dist = c(0.2, 0.2));
          venn.diagram(x = list(OBvsMsc_DOWN = OBvsMsc_DOWN$genes,
                                ADvsMsc_UP = ADvsMsc_UP$genes),
                       filename = "inter3.tiff",
                       scaled=F,lwd = 2,fill = mycol[5:6],alpha = 0.75,
                       label.col = "white",cex = 3,fontfamily = "serif",fontface = "bold",
                       cat.col = mycol[5:6],cat.cex = 3,cat.fontfamily = "serif",cat.fontface = "bold",cat.dist = c(0.2, 0.2));
#GSEA分析
      #OBvsMsc    
          #把symbol转换为ensembleid
          gene_GSEA <- mapIds(Homo.sapiens,rownames(OBvsMsc),column = "GENEID",keytype = "SYMBOL")   #发现有好多symbol没有转换成功
          GSEA_value <- OBvsMsc$logFC
          names(GSEA_value) <- gene_GSEA
          GSEA_value <- GSEA_value[order(GSEA_value,decreasing = T)]
          head(GSEA_value)
          gsea_OB_BP <- gseGO(geneList     = GSEA_value,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       nPerm        = 1000,   #下面都是默认的，是不是可以省略
                       minGSSize    = 100,
                       maxGSSize    = 500,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
          write.csv(gsea_OB_BP,file = "gsea_OB_BP.csv")
          gseaplot2(gsea_OB_BP, c("GO:0001503","GO:0006631","GO:0006869","GO:0031098","GO:0072593"),
                    color = colorspace::rainbow_hcl(5))
     #ADvsMsc     
          #把symbol转换为ensembleid
          gene_GSEA <- mapIds(Homo.sapiens,rownames(ADvsMsc),column = "GENEID",keytype = "SYMBOL")   #发现有好多symbol没有转换成功
          GSEA_value <- ADvsMsc$logFC
          names(GSEA_value) <- gene_GSEA
          GSEA_value <- GSEA_value[order(GSEA_value,decreasing = T)]
          head(GSEA_value)
          gsea_AD_BP <- gseGO(geneList     = GSEA_value,
                              OrgDb        = org.Hs.eg.db,
                              ont          = "BP",
                              nPerm        = 1000,   #下面都是默认的，是不是可以省略
                              minGSSize    = 100,
                              maxGSSize    = 500,
                              pvalueCutoff = 0.05,
                              verbose      = FALSE)
          write.csv(gsea_AD_BP,file = "gsea_AD_BP.csv")
          gseaplot2(gsea_AD_BP, c("GO:0006631","GO:0045444","GO:0045047","GO:0048705","GO:0030278"),color = colorspace::rainbow_hcl(5))
           
#GO分析
          library(clusterProfiler)
          library(enrichplot)
      #OBvsMsc
          gene_go <- AnnotationDbi::mapIds(Homo.sapiens,OBvsMsc_gene$genes, column= "ENTREZID",keytype = "SYMBOL")
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
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
          write.csv(ego_CC_OB,file = "ego_CC_OB.csv")
          write.csv(ego_BP_OB,file = "ego_BP_OB.csv")
          write.csv(ego_MF_OB,file = "ego_MF_OB.csv")
          #每组挑选10个感兴趣的通路
          term_OB <- read.csv("GO_OBs.csv",header = T)
          slimgo <- data.frame(Ontology=c(rep("Biological process",5),
                                          rep("Molecular function",5),
                                          rep("Cellular component",5)),
                               term_OB)
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
      #ADvsMsc
          gene_go <- AnnotationDbi::mapIds(Homo.sapiens,ADvsMsc_gene$genes, column= "ENTREZID",keytype = "SYMBOL")
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
          term_AD <- read.csv("GO_ADs.csv",header = T)
          slimgo <- data.frame(Ontology=c(rep("Biological process",5),
                                          rep("Molecular function",5),
                                          rep("Cellular component",5)),
                               term_AD)
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
          