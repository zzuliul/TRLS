load('Figure 6/Group.Rdata')
load('Figure 6/TCGA_mRNA.Rdata')
library(data.table)
library(tibble)
library(dplyr)
library(maftools)
library(ggsci)
maf <-  read.maf('Figure 6/TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf.gz',isTCGA = T)
Cluster <- Group[,c(1,5)]
Cluster <- Cluster%>%column_to_rownames('ID')
Cluster$sample <- substr(rownames(Cluster),1,12)
gg <- maf@variant.type.summary%>%as.data.frame()
gg$INDEL <- gg$INS+gg$DEL
colnames(gg)
gg <- merge(gg,Cluster[,1:2],by.x=1,by.y = 2)
maf2 <-  read.maf('Figure 6/TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf.gz',isTCGA = T)
c1 <- as.character(gg[gg$Group=='High',][,1])
c2 <- as.character(gg[gg$Group=='Low',][,1]) 
maf3 <- subsetMaf(maf2,tsb = gg$Tumor_Sample_Barcode)
oncoplot(maf = maf3,top = 20,writeMatrix = T,bgCol = 'WhiteSmoke',
         #clinicalFeatures = 'TMB',
         removeNonMutated = F,
         #sortByAnnotation = TRUE,
         sampleOrder = c(c1,c2),
         #annotationOrder = sort(gg$TMB,decreasing = T),
         showTitle =F,drawColBar = F)

dd <- read.table('onco_matrix.txt',h=T,sep = '\t',check.names = F)
dd[dd=='0'] <- ''
dd2 <- as.matrix(dd)
dim(dd2) <- c(ncol(dd2)*nrow(dd2),1)
unique(dd2)

if(T){
  mycol <- pal_npg("nrc")(10)[-1] 
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = 'WhiteSmoke', col = NA)) #不要背景色
    },
    Nonsense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[1], col = NA)) 
    },
    Missense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[2], col = NA)) 
    },
    Frame_Shift_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[3], col = NA)) 
    },
    Frame_Shift_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[4], col = NA)) 
    },
    Splice_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[5], col = NA)) 
    },
    In_Frame_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[6], col = NA)) 
    },
    In_Frame_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[7], col = NA)) 
    },
    Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[9], col = NA)) 
    },
    Multi_Hit = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '805736', col = NA)) 
    }
  )
  
  col = c("Nonsense_Mutation" = mycol[1], 
          "Missense_Mutation" = mycol[2], 
          "Frame_Shift_Del" = mycol[3], 
          "Frame_Shift_Ins" = mycol[4], 
          "Splice_Site" = mycol[5], 
          "In_Frame_Del" = mycol[6],
          "In_Frame_Ins" = mycol[7],
          'Translation_Start_Site' = mycol[9],
          "Multi_Hit" = '#805736')
}
library('ComplexHeatmap')
library('circlize')
oncoPrint(dd,col = col,alter_fun = alter_fun,alter_fun_is_vectorized = F)
#6*9   5*7
sub <- gg
xx <- sub[sub$Tumor_Sample_Barcode%in%colnames(dd),]
table(xx$Group)## 
xx <- xx[order(xx$Group),]
rownames(xx) <- xx[,1]
xx <- xx[,-(1:3)]
xx <- xx[,c(4,1,2,3)]
colnames(xx)[3] <- 'TMB'
top <- HeatmapAnnotation(df = xx,border = T,
                         annotation_name_side = 'right',
                         annotation_name_gp = gpar(fontsize=10),
                         show_annotation_name = T,
                         na_col = NA,
                         show_legend = T,
                         gap = unit(1, "mm"),
                         col = list(TMB = colorRamp2(c(0,10),colors = c('white','#B3CDD4')),
                                    INDEL=colorRamp2(c(0,8),colors = c('white','#BEC8D8')),
                                    SNP=colorRamp2(c(0,9),colors = c('white','#F3CC91')),
                                    Group=c('High'='#DC7B58','Low'='#4088A6')))
pdf("Figure 6/Figure6A.pdf")
oncoPrint(dd,col = col,
          alter_fun = alter_fun,
          alter_fun_is_vectorized = F,
          top_annotation = top,
          column_split = c(rep(1,149),rep(2,338)),
          column_title = NULL,
          pct_side = 'right',
          row_names_side = 'left',
          #row_names_gp = gpar(fontface='italic'),
          border=T)
dev.off()

#------------------------------------------------------------------
oncoplot(maf = maf3,top=20,writeMatrix = T,removeNonMutated = F)

pp <- t(read.table('onco_matrix.txt',h=T,row.names = 1,sep = '\t',check.names = F))%>%as.data.frame()
pp[pp=='0'] <- 'x'
pp[pp==''] <- 'x'

ID <- maf3@clinical.data
gg1 <- gg
rownames(gg1) <- gg1$Tumor_Sample_Barcode
gg1 <- gg1[ID$Tumor_Sample_Barcode,]
identical(ID$Tumor_Sample_Barcode,rownames(gg1))
ID$Group <-gg1$Group

table(ID$Group)
C1 <- C2  <- c()
for (i in colnames(pp)) {
  C1[i] <- sum(pp[ID$Tumor_Sample_Barcode[ID$Group=='High'],i]!='x')/length(ID$Tumor_Sample_Barcode[ID$Group=='High'])
  C2[i] <- sum(pp[ID$Tumor_Sample_Barcode[ID$Group=='Low'],i]!='x')/length(ID$Tumor_Sample_Barcode[ID$Group=='Low'])
}
my <- data.frame(C1=C1,C2=C2)
p <- c()
for (i in rownames(my)) {
  
  dd <- data.frame(J1=c(my[i,1]*149,(1-my[i,1])*149),
                   J2=c(my[i,2]*338,(1-my[i,2])*338))
  p <- c(p,fisher.test(dd)$p.value)
}
ll <- data.frame(p=p,sig = ifelse(p<0.0001,'****',ifelse(p<0.001,'***',ifelse(p<0.01,'**',ifelse(p<0.05,'*','')))),
                 row.names = rownames(my))
my$ID <- rownames(my)
ll$Name <- rownames(ll)
mydata <- cbind(my,ll)

mydata$C1 <- sprintf('%0.2f',mydata$C1)%>%as.numeric()
mydata$C2 <- sprintf('%0.2f',mydata$C2)%>%as.numeric()
hh <- as.matrix(mydata[,1:2])
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", hh[i, j]), x, y, gp = gpar(fontsize = 10.5))
}
right <- HeatmapAnnotation(dd = anno_text(x = mydata$sig,which = 'row',just = -0.1),which = 'row')
top <- HeatmapAnnotation(Cluster=c('High','Low'),
                         border = T,
                         annotation_name_side = 'right',
                         annotation_name_gp = gpar(fontsize=8),
                         show_annotation_name = F,
                         na_col = NA,
                         show_legend = T,
                         gap = unit(1, "mm"),
                         col = list(Cluster=c('High'='#DC7B58','Low'='#4088A6')))
plo <- Heatmap(hh,cell_fun = cell_fun,
               border = T,
               top_annotation = top,
               
               right_annotation = right,
               col = colorRamp2(c(0.09,0.45,0.8),colors = c('white','#F3CC91',pal_npg()(10)[1])),
               column_split = 1:2,
               rect_gp = gpar(color='balck'),
               column_title = NULL,show_column_names = F,
               column_title_rot = 90,
               cluster_rows =F,cluster_columns = F,
               row_names_side = 'left',
               row_names_gp = gpar(fontface='italic',color='balck'),
               heatmap_legend_param = list(title= "Frequency",
                                           at=c(0,0.4,0.8),
                                           title_position = "topcenter",
                                           legend_direction="vertical")
)
plo

pdf("Figure 6/Figures 6B.pdf")
draw(plo,heatmap_legend_side = "bottom", annotation_legend_side="right",)
dev.off()

#---------
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Rmisc)
library(dplyr)

Segment1 <- read.table('Figure 6/focal_input.seg.txt',header = T)
Segment <- Segment1
Segment$bases=Segment$End.bp-Segment$Start.bp

Segment$Sample <- substr(Segment$Sample,1,16)

data=data.frame()
for (i in 1:length(table(Segment$Sample))) {
  tmp=Segment[Segment$Sample==names(table(Segment$Sample))[i],]
  FGA=sum(tmp[abs(tmp$Seg.CN)>0.2,"bases"])/ sum(tmp[,"bases"])
  FGG=sum(tmp[tmp$Seg.CN>0.2,"bases"])/sum(tmp[,"bases"])
  FGL=sum(tmp[tmp$Seg.CN< -0.2,"bases"])/sum(tmp[,"bases"])
  tmp=data.frame(Patient=names(table(Segment$Sample))[i],FGA=FGA,FGG=FGG,FGL=FGL)
  data=rbind(data,tmp)
}

sub <- Group[,c(1,5)]
data2 <- merge(sub,data,by=1)

# -------------------------------------------------------------------------

ggdata <- pivot_longer(data2,3:5,names_to = 'FF',values_to = 'VV')
p1 <- ggplot(ggdata,aes(Group,VV,color=Group))+
  geom_boxplot(fill=NA)+
  geom_jitter(width = 0.2,size=1.5,shape=21,aes(fill=Group))+
  stat_compare_means(label = 'p.signif',method = "t.test",label.x = 1.5)+
  #scale_y_log10()+
  facet_wrap(~FF,scales = 'free')+
  scale_color_manual(values = c('#DC7B58','#4088A6'))+  #'#21b6af','#eeba4d'
  scale_fill_manual(values = c('#DC7B58','#4088A6'))+
  theme_bw(base_rect_size = 1.5)+
  labs(y='Percent Genome Altered')+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14,colour = 'black'),
        axis.title.x = element_blank(),
        strip.text = element_text(size=15),
        legend.title = element_blank(),
        #legend.position = 'none',
        legend.text = element_text(size=12,colour = 'black'),
        axis.ticks.x = element_blank())
p1
ggsave(filename = 'Figure 6/Figure 6F.pdf',width = 6,height = 4)

# -------------------------------------------------------------------------

# 
focalload <- read.table('Figure 6/focal_data_by_genes1.txt', header = T, sep = "\t", check.names = F)
rownames(focalload) <- focalload$`Gene Symbol`
head(focalload)[, 1:8]

broadload <- read.table("Figure 6/broad_values_by_arm.txt", header = T, sep = "\t", check.names = F)
rownames(broadload) <- broadload$`Chromosome Arm`


focalload <- focalload[, -c(1:3)]
focalgainload <- focalload
focallossload <- focalload
# gain 
focalgainload[focalgainload > 0]  <- 1
focalgainload[focalgainload < 0]  <- 0
focalgainload <- data.frame(colSums(focalgainload))
colnames(focalgainload) <- "focal_gain_load"
# loss
focallossload[focallossload > 0]  <- 0
focallossload[focallossload < 0]  <- 1
focallossload <- data.frame(colSums(focallossload))
colnames(focallossload) <- "focal_loss_load"

## broad
broadload <- broadload[, -c(1)]
broadload_gain <- broadload
broadload_loss <- broadload
# gain 
broadload_gain[broadload_gain > 0] <- 1
broadload_gain[broadload_gain < 0] <- 0
broadload_gain <- data.frame(colSums(broadload_gain))
colnames(broadload_gain) <- "broad_gain_load"
# loss 
broadload_loss[broadload_loss > 0] <- 0
broadload_loss[broadload_loss < 0] <- 1
broadload_loss <- data.frame(colSums(broadload_loss))
colnames(broadload_loss) <- "broad_loss_load"

##
copyloadlist <- list(focalgainload = focalgainload, focallossload = focallossload,
                     broadload_gain = broadload_gain, broadload_loss = broadload_loss)
for (i in 1:length(copyloadlist)){
  tmpdata <- copyloadlist[[i]]
  copyloadlist[[i]] <- data.frame(barcode = rownames(tmpdata), tmpdata)
}
copyload <- Reduce(function(x, y) merge(x = x, y = y, by = "barcode"), 
                   copyloadlist)
copyload$barcode <- substr(copyload$barcode,1,16)

# -------------------------------------------------------------------------
library(tibble)
dd <- merge(sub,copyload,by=1)
table(sub$Group)

colnames(dd)[3:6] <- c("Focal Gain","Focal Loss",'Arm Gain','Arm Loss')

ggdata2 <- pivot_longer(dd,3:6,names_to = 'FF',values_to = 'VV')

tmp <- rbind(ggdata,ggdata2)

p2 <- ggplot(ggdata2,aes(Group,VV,color=Group))+
  geom_boxplot(fill=NA)+
  geom_jitter(width = 0.2,size=1.5,shape=21,aes(fill=Group))+
  stat_compare_means(label = 'p.signif',method = "t.test",label.x = 1.5)+
  #scale_y_log10()+
  facet_wrap(~FF,nrow = 1,scales = 'free')+
  scale_color_manual(values = c('#DC7B58','#4088A6'))+
  scale_fill_manual(values = c('#DC7B58','#4088A6'))+
  theme_bw(base_rect_size = 1.5)+
  labs(y='Copy Number Load')+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14,colour = 'black'),
        axis.title.x = element_blank(),
        strip.text = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=12,colour = 'black'),
        axis.ticks.x = element_blank())
p2
ggsave(filename = "Figures 6/Figures 6G",width = 9.5,height = 4)

#
maf <-  read.maf('Figure 6/TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf.gz',isTCGA = T)
TMB <- tmb(maf= maf)

rownames(Group) <- substr(Group$ID,1,12)
TMB <- merge(Group,TMB,by.x = 0,by.y=1)
h <- TMB

library(ggstatsplot)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
h <- h[-340,]
myColors<-c(brewer.pal(5,"Set2"))
cor.test(h$TRLS, h$total_perMB_log,method=c("pearson", "kendall", "spearman"))

ggplot(h,aes(TRLS,total_perMB_log)) +
  geom_point(col='darkblue',alpha=0.4,size=1.5) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#EA6B66") +
  geom_rug(col="#62A398",size=0.1) +
  ggtitle(label = '')+
  theme_bw(base_rect_size = 1.5) +
  annotate('text',x =1,y=1,label='r = 0.480\np < 0.001',
           hjust=0,size=5,fontface = 'italic',color='grey30')+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=15,hjust=0.5))
ggsave(filename = 'Figure 6/Figure 6D.pdf',width = 4.5,height = 4.4)

#
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(tibble)
library(ggplot2)
library(cowplot)

tmp <- tmp1[5]$TCGA_LUAD
tmp <- tmp[,c(1,4)]
testset <- TCGA_mRNA

fnSig <- "Figure 6/pcbc-stemsig.tsv" 
w <- read.delim(fnSig, header = FALSE, row.names = 1 ) %>% as.matrix() %>% drop()
w[1:10]
#
exprSet_processed <- TCGA_mRNA
X <- exprSet_processed %>%
  rownames_to_column(var="gene_id") %>%
  filter( gene_id %in% names(w) ) %>%
  column_to_rownames( "gene_id" ) %>% as.matrix()

# Reduce the signature to the common set of genes.
stopifnot( all( rownames(X) %in% names(w)))
w <- w[ rownames(X) ]
w[1:5]

# Score the Matrix `X` using Spearman correlation. 
s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)
s[1:5]

ss <- data.frame(Sample.ID = names(s), mRNAsi = s)

input <- ss
input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$mRNAsi <- as.numeric(input$mRNAsi)
input$Group <- Group$Group
input$TRLS <- Group$TRLS

darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"

My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))


cor.test(input$TRLS, input$mRNAsi,method=c("pearson", "kendall", "spearman"))

ggplot(input,aes(TRLS,mRNAsi)) +
  geom_point(col='darkblue',alpha=0.4,size=1.5) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#EA6B66") +
  geom_rug(col="#62A398",size=0.1) +
  ggtitle(label = '')+
  theme_bw(base_rect_size = 1.5) +
  annotate('text',x =1,y=1,label='r = 0.480\np < 0.001',
           hjust=0,size=5,fontface = 'italic',color='grey30')+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=15,hjust=0.5))

ggsave(filename = 'Figure 6/Figure 6C.pdf',width = 4.5,height = 4.4)

#
library(data.table)
expr <- TCGA_mRNA
colnames(expr) <- substr(colnames(expr),1,12)

nal <- fread('Figure 6/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv',data.table = F)
nal <- data.frame(ID=nal$sample,NAL=nal$neoantigen_num)
nal <- nal[nal$ID %in% colnames(expr),]

Group$ID <- substr(Group$ID,1,12)
nal_df <- merge(nal,Group,by=1)
nal_df <- nal_df[,c(1,2,6)]
nal_df$NAL <- log(nal_df$NAL+1)
nal_df$Group <- factor(nal_df$Group)

colnames(nal_df) <- c('ID','NAL_log','Group')

#

nal_df %>%
  ggplot(aes(Group,NAL_log)) +
  geom_boxplot(aes(color=Group), outlier.colour = NA, size=1, fill=NA) +
  geom_jitter(aes(fill=Group,color=Group), width = 0.25, shape=21, size=2, alpha=0.7) +
  stat_compare_means(method = 't.test', label = 'p.signif', label.x = 1.5,size=5) +
  theme_bw(base_rect_size = 1.5) +
  labs(x=NULL,y='Relative Expression', title = 'Neoantigen') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 12),
        axis.title.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'black'),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank()) +
  scale_color_manual(values = c('#DC7B58','#4088A6')) +
  scale_fill_manual(values = c('#DC7B58','#4088A6'))

ggsave(filename = 'Figure 6/Figure 6E.pdf', width = 2.3, height = 4)
