rm(list = ls())
load('Figure 3/ALL_Cohort_Score.Rdata')
load('Figure 3/TRLS.Rdata')
load('Figure 3/New.LUAD.Rdata')
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(ggsci)
library(compareC)
library(ggkm)
library(ggplot2)
library(timeROC)
library(export)

tmp1 <- lapply(tmp1,function(x){
  rownames(x) <- x[,1]    ###
  x <-x[,1:4]
  return(x)})
rr <- data.frame()
for (i in names(tmp1)) {
  #i <- names(tmp1)[1]
  dd <- tmp1[[i]]
  scox <- summary(coxph(Surv(OS.time,OS)~RS,dd))
  cc <- scox$concordance[1]%>%as.numeric()
  se <- scox$concordance[2]%>%as.numeric()
  tt <- timeROC(dd$OS.time,dd$OS,dd$RS,cause = 1,weighting = 'marginal',times = c(1,2,3),ROC = T)
  rr <- rbind(rr,data.frame(ID=i,Cindex=cc,Cse=se,
                            AUC1=as.numeric(tt$AUC[1]),
                            AUC2=as.numeric(tt$AUC[2]),
                            AUC3=as.numeric(tt$AUC[3])))
}

AUC <- data.frame()
for (i in 1:7) {
  tt <- timeROC(tmp1[[i]]$OS.time,tmp1[[i]]$OS,
                tmp1[[i]]$RS,cause = 1,weighting = 'marginal',
                times = c(1,2,3),ROC = T)
  aa <- as.data.frame(t(tt$AUC))
  rownames(aa) <- names(tmp1)[i]
  AUC <- rbind(AUC,aa)
}
library(ggsci)
library(scales)
library(RColorBrewer)

mycol <- c('#006687','#FFDDCF','#00B7DD')
# -----------------------------------------------------------------
for (i in names(tmp1)) {

  #i <- names(rs1)[1]
  tt <- timeROC(tmp1[[i]]$OS.time,tmp1[[i]]$OS,
                tmp1[[i]]$RS,cause = 1,weighting = 'marginal',
                times = c(1,2,3),ROC = T)
  
  plot(tt,time=1,title=FALSE,lwd=2,col=mycol[1])
  plot(tt,time=2,col=mycol[2],add=TRUE,title=FALSE,lwd=2)
  plot(tt,time=3,col=mycol[3],add=TRUE,title=FALSE,lwd=2)
  legend.paste <- c(paste0("1-year AUC = ",sprintf("%.3f",tt$AUC[1])),
                    paste0("2-year AUC = ",sprintf("%.3f",tt$AUC[2])),
                    paste0("3-yaer AUC = ",sprintf("%.3f",tt$AUC[3])))
  legend(legend.paste,fill=mycol[1:3],bty="n",cex=1,
         border = NA,y.intersp=1, x.intersp=0.2,x = 0.5,y = 0.2)
  title(main = i)
  graph2pdf(file= paste0('Figure 3/',i,'-ROC.pdf'),height = 5.5,width = 5)
  dev.off()
}


#----






sur <- list()
for (i in names(LUAD_clin)[1:6]) {
  v <- LUAD_clin[[i]]
  sur1 <- merge(rs[[i]],v[,-c(1:3)],by=0)
  sur1 <- list(sur1)
  names(sur1) <- i
  sur <- c(sur,sur1)
}
#TCGA_LUAD----
tt <- sur$TCGA_LUAD
tt <- tt[,-c(5:10,13)]
table(tt$Age)
tt$Age <- ifelse(tt$Age=='>65',1,0)# >65,1
table(tt$Gender)
tt$Gender <- ifelse(tt$Gender=='Female',0,1)
table(tt$Stage)

tt$Stage <- ifelse(tt$Stage=='I',0,ifelse(tt$Stage=='II',0,
                                          ifelse(tt$Stage=='III',1,ifelse(tt$Stage=='IV',1,NA))))
table(tt$T)
tt$T <- ifelse(tt$T=='T1',0,ifelse(tt$T=='T2',0,
                                   ifelse(tt$T=='T3',1,ifelse(tt$T=='T4',1,NA))))

table(tt$N)
tt$N <- ifelse(tt$N=='N0',0,ifelse(tt$N=='N1',0,
                                   ifelse(tt$N=='N2',1,ifelse(tt$N=='N3',1,NA))))

table(tt$M)
tt$M <- substr(tt$M,2,2)

dd <- data.frame()  
for (i in colnames(tt)[4:10]) {
  
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  tt2 <- na.omit(tt[,c('OS.time','OS','RS',i)])
  p <- compareC(tt2$OS.time,tt2$OS,tt2$RS,tt2[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}

dd$ID <- c('TRLS','Age','Gender','Stage','pT','pN',"pM")
dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*',''))))
dd$ID <- factor(dd$ID,levels = dd$ID[c(9,1:8)])
ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('TCGA-LUAD')+
  ylab('C-index (Compared with TRLS)')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 50,hjust = 1),
        axis.text.y = element_text(size=15),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values =c('#5DC2CC','#C7A4C4','#70ACC2','#F5ECE6',
                              '#F6B4A6','#B5CFD4',"#91D1C2FF",'#FBC2A6'))+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))+
  geom_text(aes(y=0.7128748,label=ll),size=5)
ggsave(filename = 'Figure 3/Figure 3H.pdf',width = 4.5,height = 4.5)

#GSE30219----
tt <- sur$GSE30219
tt <- tt[,-(5:6)]
tt <- tt[,-9]
table(tt$Age)
tt$Age <- ifelse(tt$Age=='>65',1,0)
table(tt$Gender)
tt$Gender <- ifelse(tt$Gender=='Female',0,1)
table(tt$T)
tt$T <- ifelse(tt$T=='T1',0,ifelse(tt$T=='T2',0,
                                   ifelse(tt$T=='T3',1,ifelse(tt$T=='T4',1,NA))))

table(tt$N)
tt$N <- substr(tt$N,2,2)


dd <- data.frame()  
for (i in colnames(tt)[4:8]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  tt2 <- na.omit(tt[,c('OS.time','OS','RS',i)])
  p <- compareC(tt2$OS.time,tt2$OS,tt2$RS,tt2[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}



dd$ID <- c('TRLS','Age','Gender','pT','pN')
dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*','*'))))
dd$ID <- factor(dd$ID,levels = dd$ID[c(1,2:5)])

ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('GSE30219')+
  ylab('C-index (Compared with TRLS)')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(size=12,angle = 50,hjust = 1),
        axis.text.y = element_text(size=11),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values =c('#5DC2CC','#C7A4C4','#70ACC2','#F5ECE6',
                              '#F6B4A6','#B5CFD4',"#91D1C2FF",'#FBC2A6'))+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.85))+
  geom_text(aes(y=0.818119,label=ll),size=5)
ggsave(filename = 'Figure 3/Figure 3L.pdf',width = 4.5,height = 4.5)
#GSE31210----


tt <- sur$GSE31210
tt <- tt[,-c(5,6,9)]
table(tt$Age)
tt$Age <- ifelse(tt$Age=='>65',1,0)##>65 1
table(tt$Gender)
tt$Gender <- ifelse(tt$Gender=='Female',0,1)
table(tt$Stage)
tt$Stage <- ifelse(tt$Stage=='I',0,1)
table(tt$Smoke)
tt$Smoke <- ifelse(tt$Smoke=='Yes',1,0)
table(tt$EGFR)
tt$EGFR <- ifelse(tt$EGFR=='WT',0,1)## Mut 1
table(tt$KRAS)
tt$KRAS <- ifelse(tt$KRAS=='WT',0,1)## Mut 1
table(tt$ALK)
tt$ALK <- ifelse(tt$ALK=='WT',0,1)
dd <- data.frame()  
for (i in colnames(tt)[4:11]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  tt2 <- na.omit(tt[,c('OS.time','OS','RS',i)])
  p <- compareC(tt2$OS.time,tt2$OS,tt2$RS,tt2[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}

dd$ID <- c('TRLS','Age','Gender','Stage','Smoke','EGFR','KRAS','ALK')
dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*',''))))
dd$ID <- factor(dd$ID,levels = dd$ID[c(1,2:8)])

ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('GSE31210')+
  ylab('C-index (Compared with TRLS)')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(size=12,angle = 50,hjust = 1),
        axis.text.y = element_text(size=11),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values =c('#5DC2CC','#C7A4C4','#70ACC2','#F5ECE6',
                              '#F6B4A6','#B5CFD4',"#91D1C2FF",'#FBC2A6'))+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.80))+
  geom_text(aes(y=0.7569849,label=ll),size=5)
ggsave(filename = 'Figure 3/Figure 3K.pdf',width = 4.5,height = 4.5)

#


# GSE50081 ---------------------------------------------------------------
tt <- sur$GSE50081
tt <- tt[,-c(5,6)]
tt <- tt[,-10]
table(tt$Age)
tt$Age <- ifelse(tt$Age=='>65',1,0)##>65 1
table(tt$Gender)
tt$Gender <- ifelse(tt$Gender=='Female',0,1)
table(tt$Smoke)
tt$Smoke <- ifelse(tt$Smoke=='Yes',1,0)
table(tt$Stage)
tt$Stage <- ifelse(tt$Stage=='I',0,1)
table(tt$T)
tt$T <- ifelse(tt$T=='T1',0,ifelse(tt$T=='T2',0,
                                   ifelse(tt$T=='T3',1,ifelse(tt$T=='T4',1,NA))))
table(tt$N)
tt$N <- ifelse(tt$N=='N0',0,1)

dd <- data.frame()  
for (i in colnames(tt)[4:10]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  tt2 <- na.omit(tt[,c('OS.time','OS','RS',i)])
  p <- compareC(tt2$OS.time,tt2$OS,tt2$RS,tt2[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}

dd$ID <- c('TRLS','Age','Gender','Stage','T','N','Smoke')
dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*',''))))
dd$ID <- factor(dd$ID,levels = dd$ID[c(1,2:7)])

#
ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('GSE50081')+
  ylab('C-index (Compared with TRLS)')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(size=12,angle = 50,hjust = 1),
        axis.text.y = element_text(size=11),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values =c('#5DC2CC','#C7A4C4','#70ACC2','#F5ECE6',
                              '#F6B4A6','#B5CFD4',"#91D1C2FF",'#FBC2A6'))+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))+
  geom_text(aes(y=0.7032752,label=ll),size=5)
ggsave(filename = 'Figure 3/Figure 3J.pdf',width = 4.5,height = 4.5)


# GSE72094 ---------------------------------------------------------------
tt <- sur$GSE72094
tt <- tt[,-13]
table(tt$Age)
tt$Age <- ifelse(tt$Age=='>65',1,0)
table(tt$Gender)
tt$Gender <- ifelse(tt$Gender=='Female',0,1)#Male, 1
table(tt$Smoke)
tt$Smoke <- ifelse(tt$Smoke=='Yes',1,0)
table(tt$Stage)
tt$Stage <- ifelse(tt$Stage=='I',0,ifelse(tt$Stage=='II',0,
                                          ifelse(tt$Stage=='III',1,ifelse(tt$Stage=='IV',1,NA))))
table(tt$EGFR)
tt$EGFR <- ifelse(tt$EGFR=='WT',0,1)## Mut 1
table(tt$KRAS)
tt$KRAS <- ifelse(tt$KRAS=='WT',0,1)## Mut 1
table(tt$TP53)
tt$TP53 <- ifelse(tt$TP53=='WT',0,1)
table(tt$STK11)
tt$STK11 <- ifelse(tt$STK11=='WT',0,1)

dd <- data.frame()  
for (i in colnames(tt)[4:12]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  tt2 <- na.omit(tt[,c('OS.time','OS','RS',i)])
  p <- compareC(tt2$OS.time,tt2$OS,tt2$RS,tt2[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}

dd$ID <- c('TRLS','Age','Gender','Stage','EGFR','KRAS','TP53','STK11','Smoke')
dd$ll <- ifelse(dd$P<0.0001,'****',ifelse(dd$P<0.001,'***',ifelse(dd$P<0.01,'**',ifelse(dd$P<0.05,'*',''))))
dd$ID <- factor(dd$ID,levels = dd$ID[c(1,2:9)])
#
ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle('GSE72094')+
  ylab('C-index (Compared with TRLS)')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(size=12,angle = 50,hjust = 1),
        axis.text.y = element_text(size=11),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values =c('#5DC2CC','#C7A4C4','#70ACC2','#F5ECE6',
                              '#F6B4A6','#B5CFD4',"#91D1C2FF",'#FBC2A6','#008594'))+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.79))+
  geom_text(aes(y=0.7376291,label=ll),size=5)
ggsave(filename = 'Figure 3/Figure 3I.pdf',width = 4.5,height = 4.5)
