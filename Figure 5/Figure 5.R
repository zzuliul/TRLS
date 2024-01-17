rm(list = ls())
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(GSVA)
library(ggsci)
library(tidyr)
library(ggpubr)

load('Figure 5/TCGA_mRNA.Rdata')
load('Figure 5/Group.Rdata')
load("Figure 5/immune_landscape.RData")
load('Figure 5/CIC_Marker.rda')

cellMarker <- data.table::fread("Figure 5/cellMarker.csv")
colnames(cellMarker)[2] <- "celltype"
type <- split(cellMarker,cellMarker$celltype)
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})


Group <- my[,c(1:4,156)]
colnames(Group)[4] <- 'TRLS'
expr <- as.matrix(TCGA_mRNA)
#
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
a <- gsva_data %>% t() %>% as.data.frame()
a$group <- Group$Group
b <- gather(a,key=ssGSEA,value = Expression,-c(group))
ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = c("#E2AC72","#9FBAE0"),outlier.size = 0.005,outlier.alpha = 0.6)+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

ggsave(filename = 'Figure 5/Figure 5A.pdf', width = 12, height = 6.2)
dev.off()

#----
library(ggplot2)
library(data.table)
library(tidyverse)
ABSOLUTE_scores$ID <- substr(rownames(ABSOLUTE_scores),1,12)
ABSOLUTE_scores <- ABSOLUTE_scores[!duplicated(ABSOLUTE_scores$ID),-5]
rownames(ABSOLUTE_scores) <- substr(rownames(ABSOLUTE_scores),1,12)

LUAD_im <- immune_landscape[immune_landscape$TCGA.Study == "LUAD",]
LUAD_ab <- ABSOLUTE_scores[rownames(LUAD_im),]
LUAD_il<- merge(LUAD_im,LUAD_ab,by=0) %>%  column_to_rownames("Row.names")
rm(immune_landscape,ABSOLUTE_scores,LUAD_im,LUAD_ab)

rownames(Group) <- substr(Group$ID,1,12)
metadata <- merge(Group,LUAD_il,by = 0) %>% 
  column_to_rownames("Row.names")
ggdata <- metadata[,-c(1:4,6:8)] 
ggdata <- apply(ggdata[,-1],2,as.numeric) %>% 
  as.data.frame()
ggdata$Group <- factor(metadata$Group)
#--------------------------------------------------------------
colnames(ggdata)
img <- c("Group",
         "SNV.Neoantigens","Indel.Neoantigens","CTA.Score",
         "Intratumor.Heterogeneity",  
         "Number.of.Segments","Fraction.Altered",
         "LOH_n_seg","LOH_frac_altered",
         "Homologous.Recombination.Defects","Aneuploidy.Score",
         "TCR.Richness","TCR.Shannon")
dat <- ggdata[,img] 
rownames(dat) <- rownames(metadata)
df <- data.frame()
for (i in 2:13) {
  m1=  mean(dat[dat$Group=="High",i][!is.na(dat[dat$Group=="High",i])])/mean(dat[,i][!is.na(dat[,i])]) ;m1=log2(m1)
  m2= mean(dat[dat$Group=="Low",i][!is.na(dat[dat$Group=="Low",i])])/mean(dat[,i][!is.na(dat[,i])]);m2=log2(m2)
  df <- rbind(df,data.frame("High"=m1,"Low"=m2))
}
df <- df %>% as.matrix()

rownames(df) <- c("SNV Neoantigens","Indel Neoantigens","CTA Score",
                  
                  "Intratumor Heterogeneity",  
                  "Number of Segments","Fraction Altered",
                  "LOH_n_seg","LOH_frac_altered",
                  
                  "Homologous Recombination Defects","Aneuploidy Score",
                  
                  "TCR Richness","TCR Shannon")

colnames(df) <- c("High","Low")

library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(export)
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", df[i, j]),
            x, y, gp = gpar(fontsize = 10.5))
}
top <- HeatmapAnnotation(Cluster=c("High","Low"),
                         border = T,
                         annotation_name_side = 'right',
                         annotation_name_gp = gpar(fontsize=10),
                         show_annotation_name = F,
                         na_col = NA,
                         show_legend = F,
                         gap = unit(1, "mm"),
                         col = list(Cluster=c('MSC1'="#EFC000FF",'MSC2'="#197EC099")))
plo <- Heatmap(df,
               #cell_fun = cell_fun,
               border = T,
               heatmap_width = unit(10, "cm"),
               heatmap_height = unit(1, "cm")*nrow(df),
               #top_annotation = top,
               #right_annotation = right,
               col = colorRamp2(c(-1,0,1),colors = c("#009698","white","#cd5c5c")),
               #column_split = 1:2,
               row_split = rep(c(1,2,3,4),time= c(3,5,2,2)),
               gap = unit(2, "mm"),
               rect_gp = gpar(color='grey',lwd=1.2,size=0.4),
               
               show_column_names = F,
               cluster_rows = F,cluster_columns = F,
               
               row_names_side = 'right',
               row_names_gp = gpar(fontface='bold',color='balck',fontsize = 12),
               
               heatmap_legend_param = list(title= "log2(Ratio)",
                                           at=c(-1,0,1),
                                           title_position = "topcenter",
                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                           legend_direction="horizontal"))
plo

draw(plo,heatmap_legend_side = "bottom")
graph2pdf(file='Figure 5/Figure 5B.pdf',height=6.2,width=5)


#----
library(ggcor)
library(dplyr)
library(tibble)
library(readxl)
library(GSVA)
library(ggplot2)
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}



#
load('Figure 5/CIC_Marker.rda')
expr <- TCGA_mRNA
names(CIC_Marker) <- c(paste("step1:","Release of cancer cell antigens"),
                       paste("step2:","Cancer antigen presentation"),
                       paste("step3:","Priming and activation"),
                       paste("step4:","Basophil recruiting"),
                       paste("step4:","B cell recruiting"),
                       paste("step4:","CD4 T cell recruiting"),
                       paste("step4:","CD8 T cell recruiting"),
                       paste("step4:","Dendritic cell recruiting"),
                       paste("step4:","Eosinophil recruiting"),
                       paste("step4:","Macrophage recruiting"),
                       paste("step4:","MDSC recruiting"),
                       paste("step4:","Monocyte recruiting"),
                       paste("step4:","Neutrophil recruiting"),
                       paste("step4:","NK cell recruiting"),
                       paste("step4:","T cell recruiting"),
                       paste("step4:","Th1 cell recruiting"),
                       paste("step4:","Th22 cell recruiting"),
                       paste("step4:","Th2 cell recruiting"),
                       paste("step4:","Treg cell recruiting"),
                       paste("step5:","Infiltration of immune cells into tumors"),
                       paste("step6:","Recognition of cancer cells by T cells"),
                       paste("step7:","Killing of cancer cells"))



#
immPath.score <- gsva(expr = as.matrix(expr),
                      CIC_Marker, 
                      method = "ssgsea")

CIC.score <- immPath.score
a <- as.data.frame(t(CIC.score))

Group <- Group[,c(1,4)]
clin <- merge(a,Group,by.x=0,by.y=1)
clin <- clin%>%column_to_rownames('Row.names')

clin <- as.data.frame(t(clin))


outTab=data.frame()
corFilter=0             
pvalueFilter=1    
for(i in row.names(clin)){
  for(j in row.names(clin)){
    x=as.numeric(clin[i,])
    y=as.numeric(clin[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="postive"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="negative"))
    }
  }
}

table(outTab$d1)
sub <- outTab[outTab$d1%in%'TRLS',]
sub <- sub[order(sub$d1,decreasing = F),]
sub <- sub[!sub$d2%in%'TRLS',]
colnames(sub)[3:4] <- c('r','p.value')
sub$p.value <- as.numeric(sub$p.value)
sub$pd <- ifelse(sub$p.value<0.05,'< 0.05','>= 0.05')
sub$r <- as.numeric(sub$r)
sub$rd <- cut(sub$r, breaks = c(-Inf,-0.4,-0.2,0.2,0.4, Inf),
              labels = c("<= -0.4",'-0.4 - -0.2','-0.2 - 0.2', "0.2 - 0.4", ">= 0.4"))
table(sub$rd)
table(sub$d1)
sub$d1[sub$d1=='TRLS'] <- 'TRLS'
varechem <- as.data.frame(a)
col <- c("#26B1D0","#DA6003","#A2A2A288","#E8288E","#65A818")
show_col(col)
quickcor(varechem, type = "upper") + 
  geom_colour() +
  anno_link(aes(colour = rd, size = pd,), 
            #,linetype = Regulation), 
            data = sub) +
  scale_size_manual(values = c(1, 0.5)) +
  #scale_linetype_manual(values = c("dashed","solid"))+
  scale_color_manual(values = c("#6FB7CCFF","#AFD7E2FF","#A2A2A288","#FFCACAFF","#FFA1A1FF")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  remove_axis("x")+
  guides(colour = guide_legend(title = "Subtype's r", 
                               override.aes = list(size = 3), 
                               order = 1),
         size = guide_legend(title = "Subtype's p",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         fill = guide_colorbar(title = "Pearson's r", order =3))

ggsave(filename = "Figures/Figure5/Figure5C in top right", width = 10,height = 8)





immPath <- read_excel("Figure 5/Table S9.xlsx",sheet = 1,skip = 1)
immPath <- as.data.frame(immPath)[1:18,]
rownames(immPath) <- immPath[,1]
immPath <- immPath[,-1]
immPath.list <- list()
for (i in rownames(immPath)) {
  tmp <- immPath[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp,",",fixed = T)))
  tmp <- gsub(" ","",tmp)
  immPath.list[[i]] <- tmp
}

#
immPath.score <- gsva(expr = as.matrix(expr),
                      immPath.list, 
                      method = "ssgsea")
immtherapyPath.score <- immPath.score

#


a <- as.data.frame(t(immtherapyPath.score))
colnames(a)[1] <- "IFN_Gamma_signiture"


clin <- merge(a,Group,by.x=0,by.y=1)
clin <- clin%>%column_to_rownames('Row.names')
clin <- as.data.frame(t(clin))

outTab=data.frame()
corFilter=0             
pvalueFilter=1    
for(i in row.names(clin)){
  
  for(j in row.names(clin)){
    x=as.numeric(clin[i,])
    y=as.numeric(clin[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="postive"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(d1=j,d2=i,cor,pvalue,Regulation="negative"))
    }
  }
}

table(outTab$d1)

sub <- outTab[outTab$d1%in%'TRLS',]
sub <- sub[order(sub$d1,decreasing = F),]
sub <- sub[!sub$d2%in%'TRLS',]
colnames(sub)[3:4] <- c('r','p.value')
sub$p.value <- as.numeric(sub$p.value)
sub$pd <- ifelse(sub$p.value<0.05,'< 0.05','>= 0.05')
sub$r <- as.numeric(sub$r)
sub$rd <- cut(sub$r, breaks = c(-Inf,-0.4,-0.2,0.2,0.4, Inf),
              labels = c("<= -0.4",'-0.4 - -0.2','-0.2 - 0.2', "0.2 - 0.4", ">= 0.4"))
table(sub$rd)
table(sub$d1)
sub$d1[sub$d1=='TRLS'] <- 'TRLS'
varechem <- as.data.frame(a)

quickcor(varechem, type = "lower") + 
  geom_colour() +
  anno_link(aes(colour = rd, size = pd), data = sub) +
  scale_size_manual(values = c(1, 0.5)) +
  scale_linetype_manual(values = c("dashed","solid"))+
  scale_color_manual(values = c("#6FB7CCFF","#AFD7E2FF","#A2A2A288","#FFCACAFF","#FFA1A1FF")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = '#E11953',midpoint=0) +
  remove_axis("x")+
  guides(colour = guide_legend(title = "Subtype's r", 
                               override.aes = list(size = 3), 
                               order = 1),
         size = guide_legend(title = "Subtype's p",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         fill = guide_colorbar(title = "Pearson's r", order =3))

ggsave(filename = "Figures/ggcor plot in bottom left.pdf", width = 10,height = 8)




TCPA_expr <- read.csv(file = 'Figure 5/TCGA-LUAD-L4.csv')
TCPA_expr$Sample_ID <- substr(TCPA_expr$Sample_ID,1,16)
TCPA_expr <- TCPA_expr[,-c(2,3,4)]
TCPA <- merge(Group,TCPA_expr,by=1)
gg <- c('PD1','PDL1','CTLA4')
expr <- TCPA[,c('ID','TRLS','PDL1')]
expr <- expr%>%column_to_rownames('ID')

expr1 <- expr[,1:2]

cor.test(expr$TRLS, expr$PDL1,method=c("pearson", "kendall", "spearman"))

ggplot(expr,aes(RS,PDL1)) +
  geom_point(col='#D8975E',alpha=0.6,size=1.5) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#BF6F87") +
  geom_rug(col="#515792",size=0.1) +
  ggtitle(label = '')+
  theme_bw(base_rect_size = 1.5) +
  annotate('text',x =0.6,y=1.9,label='r = 0.149\np < 0.001',
           hjust=0,size=5,fontface = 'italic',color='grey30')+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=15,hjust=0.5))
ggsave(filename = 'Figure 5/Figure 5E.pdf', width = 5.5, height = 5)

expr <- as.data.frame(t(TCGA_mRNA))
expr <- merge(Group,expr,by.x =1,by.y = 0)
expr <- expr%>%column_to_rownames('ID')
expr <- expr[,c('TRLS','CD274')]
cor.test(expr$TRLS, expr$CD274,method=c("pearson", "kendall", "spearman"))

ggplot(expr,aes(TRLS,CD274)) +
  geom_point(col='#F6B4A6',alpha=0.95,size=1.5) +
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#BF6F87") +
  geom_rug(col="#515792",size=0.1) +
  ggtitle(label = '')+
  theme_bw(base_rect_size = 1.5) +
  annotate('text',x =0.55,y=6.85,label='r = 0.149\np < 0.001',
           hjust=0,size=5,fontface = 'italic',color='grey30')+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=15,hjust=0.5))
ggsave(filename = 'Figure 5/Figure5 D.pdf', width = 5.5, height = 5)


#----

load('Figure 5/Group.Rdata')
result <- read.csv('Figure 5/TIDE_output.csv',h=T)
table(result$Responder)
result2 <- merge(result[,c(1,4)],Group[,c(1,5)],by=1)
library(ggplot2)
library(ggpubr) 

ggplot(result2,aes(Group,TIDE,fill=Group))+ 
  geom_boxplot(aes(color=Group), outlier.colour = NA, size=1, fill=NA) +
  geom_jitter(aes(fill=Group,color=Group), width = 0.2, shape=21, size=1.5, alpha=0.7)+
  theme_classic()+
  ylab('TIDE Score')+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_npg()+
  stat_compare_means(comparisons = list(c('High','Low')),
                     label = 'p.signif',
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")))+ 
  stat_compare_means(label.y = max(result2$TIDE)+1)+
  scale_color_manual(values = c('#DC7B58','#4088A6')) +
  scale_fill_manual(values = c('#DC7B58','#4088A6'))
ggsave(filename = "Figures/Figure5/Figure 5F.pdf", width = 10,height = 8)

#----------
#submap----
#GSE93157
tmp <- matrix(c( 1.0000000,  0.07592408 , 1 ,1.000000000,
                 0.5824176, 0.01898102, 0.5044955 ,0.95004995), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("High_p","Low_p","High_b","Low_b"),c("PDL1-noR","PDL1-R")))

pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[1:5],
         gaps_row = 2,
         display_numbers = matrix(ifelse(tmp>0.01&tmp < 0.05,paste0('p=',round(tmp,3)), ifelse(tmp<0.01,'p<0.01','')),nrow(tmp)),number_color = "black",
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "Figure 5/Figure5 G.pdf"
)

