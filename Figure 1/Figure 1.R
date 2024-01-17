rm(list = ls())
load('Figure 1/LUAD_FQ.RData')
load('Figure 1/TCGA_LUAD.Rdata')
load('Figure 1/Immune.Rdata')
load('Figure 1/genelist.Rdata')


library(scRNAtoolVis)
library(Seurat)
library(jjPlot)
library(grid)
#----
reduc <-
  data.frame(Seurat::Embeddings(LUAD_FQ, reduction = "umap"))
#metadata
meta <- LUAD_FQ@meta.data
#combine
combine <- cbind(reduc,meta)

label_id <- combine %>% group_by(seurat_clusters) %>%
  summarise(x_m = median(UMAP_1),y_m = median(UMAP_2))
label_id$facet <- "umap"

cell_num <- combine %>% group_by(seurat_clusters,celltype) %>%
  summarise(num = n())
cell_num$y <- 1:nrow(cell_num)

combine$facet <- "umap"
cell_num$facet <- "numbers"

#
draw_number_circle <- function(data,params,size){
  grobTree(pointsGrob(x = 0.5,y = 0.5,
                      size = unit(2,"char"),
                      pch = 16,
                      gp = gpar(
                        col = alpha(data$colour %||% "grey50",data$alpha),
                        fill = alpha(data$fill %||% "grey50",data$alpha),
                        lwd = (data$linewidth %||% 0.5)* .pt,
                        lty = data$linetype %||% 1)),
           textGrob(label = data$label,
                    x = rep(0.5,3),y = rep(0.5,3),
                    gp = gpar(col = "black"))
  )
  
}

color <- c('#5DC2CC','#70ACC2','#F5ECE6','#E3B8B8','#F9A647','#D7AA23',
       '#959434','#F6B4A6','#B5CFD4')
ggplot(data = combine,aes(x = UMAP_1,y = UMAP_2,
                          color = celltype,label = celltype)) +
  geom_point(key_glyph = draw_number_circle) +
  #geom_markArrow(rel.pos = 0.05) +
  coord_cartesian(xlim = c(-15,15),ylim = c(-15,15)) +
  guides(color = guide_legend(override.aes = list(label = 0:8)),
         size = 50)+
  scale_color_manual(values = color)+
  #geom_point(size=0.05,shape=1,color='brown')+
  labs(x='UMAP1',y='UMAP2')+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'darkred',face='bold'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'))
ggsave(filename = 'Figure 1/Figure 1A.pdf',width =7.5 ,height = 4)


#----
reduc <-
  data.frame(Seurat::Embeddings(Immune, reduction = "umap"))
#metadata
meta <- Immune@meta.data
#combine
combine <- cbind(reduc,meta)

label_id <- combine %>% group_by(seurat_clusters) %>%
  summarise(x_m = median(umap_1),y_m = median(umap_2))
label_id$facet <- "umap"

cell_num <- combine %>% group_by(seurat_clusters,celltype) %>%
  summarise(num = n())
cell_num$y <- 1:nrow(cell_num)

combine$facet <- "umap"
cell_num$facet <- "numbers"

color <- c('#5DC2CC','#70ACC2','#F5ECE6','#959434','#F6B4A6',
           '#B5CFD4')
ggplot(data = combine,aes(x = umap_1,y = umap_2,
                          color = celltype,label = celltype)) +
  geom_point(key_glyph = draw_number_circle) +
  #geom_markArrow(rel.pos = 0.05) +
  coord_cartesian(xlim = c(-15,10),ylim = c(-15,10)) +
  guides(color = guide_legend(override.aes = list(label = 0:5)),
         size = 50)+
  scale_color_manual(values = color)+
  #geom_point(size=0.05,shape=1,color='brown')+
  labs(x='UMAP1',y='UMAP2')+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'darkred',face='bold'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'))
ggsave(filename = 'Figure 1/Figure 1B.pdf',width =5.5 ,height = 4)
#----
library(MySeuratWrappers)
my36colors <- c('#5DC2CC','#70ACC2','#F5ECE6','#E3B8B8','#F9A647','#D7AA23',
                '#959434','#F6B4A6','#B5CFD4')
Idents(LUAD_FQ) <- LUAD_FQ@meta.data$celltype 

pdf(file='Figure 1/Figure 1C.pdf',width=7.5,height=4)
VlnPlot(LUAD_FQ, features = markerGenes,  
        stacked=T,pt.size=0,  
        cols = my36colors,
        direction = "horizontal",
        x.lab = NULL, y.lab = NULL)+
  theme(axis.text.x = element_blank(),   
        axis.ticks.x =element_blank())+
  scale_x_discrete(limits = rev(c( "B cells", "CD4 +T cells","dendritic cells","endothelial cells and fibroblasts",
                                   "epithelial cells","mast cells",'monocytic cells','neutrophil cells','NK T cells')))
dev.off()

#----
library(clusterProfiler)
library(org.Hs.eg.db)
logtpm <- logtpm[genelist$ENSEMBL,]
logtpm <- na.omit(logtpm)

library(tibble)
data <- logtpm
data <- t(data)%>%as.data.frame()
riskdata <- merge(TCGA_LUAD_clin[1:3],data,by.x = 1,by.y = 0)
library(survival)
library(survminer)
riskdata <- riskdata[riskdata$OS.time>0,]
riskdata <- riskdata[!is.na(riskdata$OS.time),]  

outTab=data.frame()
riskdata$OS <- as.numeric(riskdata$OS)

for(i in colnames(riskdata[,4:ncol(riskdata)])){
  cox <- coxph(Surv(OS.time,OS) ~ get(i), data = riskdata)
  coxSummary = summary(cox)
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

outTab[,2:5] <- apply(outTab[,2:5],2,as.numeric)
gene <- subset(outTab,pvalue< 0.05)

NMF_expr <- logtpm[gene$id,]
library(ConsensusClusterPlus) ## 
## CC 
results = ConsensusClusterPlus(as.matrix(NMF_expr),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='Figure 1/ConsensusCluster/new3',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="png")
icl <- calcICL(results,title = 'Figure 1/ConsensusCluster/new',plot = 'png')

## PAC
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK

PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
library(ggthemr)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(size=4,shape=21,color='darkred',fill='skyblue')+
  ylab('Proportion of ambiguous clustering')+
  xlab('Cluster number K')
dev.off()
#4.5*6
## 


clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]
sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)##
table(sub$Cluster) #C1 289 C2 208
group <- substr(sub$Cluster,2,2)
nmf.input <- NMF_expr

cc <- sub$Cluster
names(cc) <- sub$Sample
cc2 <- sort(cc)

my <- results[[2]][["ml"]]
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
my2 <- my[names(cc2),names(cc2)]
library(pheatmap)
library(tidyverse)
library(ggsci)
library(scales)
pdf(file='Figure 1/pheat-cluster.pdf',width=5.5,height=4.5)
pheatmap(1-my2,show_colnames = F,show_rownames = F,
         treeheight_row = 20,treeheight_col = 20,
         cluster_rows = F,cluster_cols = F,
         clustering_method = 'complete',border_color = 'black',
         color = colorRampPalette(c("#DC8C6B","white"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Subtype=cc2),
         annotation_col = data.frame(Subtype=cc2),
         annotation_colors = list(Subtype=c('C1'='#F6B4A6','C2'='#70ACC2')))
dev.off()


load('Figure 1/Cluster.rda')
load('Figure 1/TCGA_mRNA.Rdata')
library(IOBR)
library(ComplexHeatmap)
TCGA_mRNA <- TCGA_mRNA[,Cluster$Sample]
estimate <- deconvo_tme(eset = TCGA_mRNA, method = "estimate")
rownames(estimate) <- rownames(Cluster)
estimate$group <- Cluster$Cluster
estimate <- estimate %>% rownames_to_column("sample")
estimate <- estimate[,c('ImmuneScore_estimate','group')]
b <- gather(estimate,key= Estimate,value = Proportion,-c(group))

ggplot(estimate,aes(group,ImmuneScore_estimate)) +
  geom_boxplot(aes(color=group), outlier.colour = NA, size=1, fill=NA) +
  geom_jitter(aes(fill=group,color=group), width = 0.2, shape=21, size=3, alpha=0.7) +
  theme_bw(base_rect_size = 1.5) +
  labs(x=NULL,y='Immune score', title = 'Immune score') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 2,size = 15),
        axis.title.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank()) +
  scale_color_manual(values = c('#E9B074','#b5D6DF')) +
  scale_fill_manual(values = c('#E9B074','#b5D6DF'))+  
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
ggsave(filename = 'Figure 1/Figure 1E.pdf',width = 4,height = 7)

load('Figure 1/module12.Rdata')
load('Figure 1/module21.Rdata')
pdf(file='Figure 1/Figure 1G.pdf',width=6,height=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes12, column12]),
                   
                   abs(geneTraitSignificance[moduleGenes12, 1]),
                   
                   xlab = paste("Module Membership in module 12"),
                   
                   ylab = "Gene significance for Cluster",
                   
                   main = paste("Module membership vs. gene significance\n"),
                   
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = '#b5D6DF',
                   pch = 21,bg= alpha('#b5D6DF',0.8) ,cex = 1.7)
dev.off()

pdf(file='Figure 1/Figure 1F.pdf',width=6,height=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes21, column21]),
                   
                   abs(geneTraitSignificance[moduleGenes21, 1]),
                   
                   xlab = paste("Module Membership in module 21"),
                   
                   ylab = "Gene significance for Cluster",
                   
                   main = paste("Module membership vs. gene significance\n"),
                   
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = '#E9B074',
                   pch = 21,bg= alpha('#E9B074',0.8) ,cex = 1.7)
dev.off()







