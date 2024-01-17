rm(list = ls())
library(ggplot2)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggbreak)
library(tibble)
load('Figure 2/model.Rdata')
model$Model <- gsub('\\[',' [',model$Model)
model$Model <- gsub('α','alpha',model$Model)
model$Model <- gsub('(Enet \\[alpha=1])','Lasso',model$Model)
model$Model <- gsub('(Enet \\[alpha=0])','Ridge',model$Model)

dd <- model%>%group_by(Model)%>%summarise(Cindex=mean(Cindex))
range(dd$Cindex)
dd$ll <- sprintf("%.3f",dd$Cindex)

dd%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width=0.7,stat = 'identity',fill='#1993CF',position = position_dodge(width = 0.9))+
  scale_x_continuous(breaks = c(0,0.60,0.7),labels = c(0,0.60,0.7),limits = c(0,1))+
  scale_x_break(c(0.05,0.54),scales = 20)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  ggtitle('Mean Cindex',)+
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.15,vjust = -3),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())+
  geom_text(aes(label=ll),hjust=1.1,size=2.7,color='black')
ggsave(filename = 'Figure 2/Figure 2A right.pdf',width = 3.5,height = 17)


dd2 <- pivot_wider(model,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)
dd3 <- dd2
dd3$Mean <- rowMeans(dd2[,-1])
dd3 <- dd3[order(dd3$Mean,decreasing = T),]
rownames(dd3) <- NULL
dd3 <- dd3%>%column_to_rownames('Model')

Cohort <- c('#0F505E','#E8D0C4','#AFAD88','#D16474','#C9A8D4')
names(Cohort) <- colnames(dd3)[-6]

Top = HeatmapAnnotation(Cohort=colnames(dd3)[-6],
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 9,col='black'),
                                                     title_gp = gpar(fontsize = 9, fontface = "bold",col='black'),ncol=1),
                        gap=unit(1, "mm"),
                        border = T,
                        col=list(Cohort=Cohort),
                        show_annotation_name = TRUE,
                        annotation_name_side="right",
                        annotation_name_gp = gpar(fontsize = 10,col='black', fontface = "bold"))
dd3[,-6] <- round(dd3[-6],3)
range(dd3[,-6])

pdf('Figure 2/Figure 2A left.pdf',width = 6.5,height=17)
Heatmap(as.matrix(dd3[,-6]),col = colorRamp2(c(0.49,0.625,0.76),c('#b5D6DF','white','#E9B074')),
        rect_gp = gpar(color='grey30',lwd=0.2),name = 'Cindex',
        cluster_rows = F,cluster_columns = F,column_names_gp = gpar(fontsize=13),
        row_names_side = 'right',column_names_side = 'top',column_names_rot = 0,column_names_centered = T,
        cell_fun = function(j,i,x,y,width,height,fill){grid.text(dd3[i,j],x,y,gp = gpar(fontsize = 9))},
        column_split = 1:5,top_annotation = Top,
        show_column_names = F,column_title = NULL,
        row_title = NULL
)
dev.off()
