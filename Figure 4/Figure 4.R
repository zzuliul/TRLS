rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(compareC)
library(ggsci)
load('Figure 4/ALL_Cohort_Score.Rdata')
uni <- data.frame()
for (i in names(tmp)) {
  dd <- tmp[[i]]
  for(j in colnames(dd)[4:ncol(dd)]){
    scox <- summary(coxph(Surv(OS.time,OS)~get(j),dd))
    p <- compareC(dd$OS.time,dd$OS,dd$RS,dd[,j])$pval
    uni <- rbind(uni,data.frame(ID=i,A=j,
                                HR=scox$conf.int[,1],
                                HR.95L=scox$conf.int[,3],
                                HR.95R=scox$conf.int[,4],
                                pvalue=scox$coefficients[,5],
                                cindex=scox$concordance[1]%>%as.numeric(),
                                cse=scox$concordance[2]%>%as.numeric(),
                                CP=p))
  }
}

uni$ID <- factor(uni$ID,levels = c("TCGA_LUAD",names(tmp)[-5]))

####C index####
####TCGA
unique(uni$ID)
dd <- uni[uni$ID=='TCGA_LUAD',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill='#990000')+
  ylab(NULL)+xlab(NULL)+
  labs(title ="TCGA")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-TCGA_LUAD-signture.pdf',width = 4,height =21)


####GSE30219
unique(uni$ID)
dd <- uni[uni$ID=='GSE30219',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_futurama(alpha = 0.7)(12)[3])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE30219")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-GSE30219-signture.pdf',width = 4,height =21)


########
unique(uni$ID)
dd <- uni[uni$ID=='GSE31210',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_futurama(alpha = 0.7)(12)[4])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE31210")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-GSE31210-signture.pdf',width = 4,height =21)


########
unique(uni$ID)
dd <- uni[uni$ID=='GSE50081',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_futurama(alpha = 0.7)(12)[5])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE50081")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-GSE50081-signture.pdf',width = 4,height =21)


########
unique(uni$ID)
dd <- uni[uni$ID=='GSE3141',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_futurama(alpha = 0.7)(12)[6])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE3141")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-GSE3141-signture.pdf',width = 4,height =21)


########
unique(uni$ID)
dd <- uni[uni$ID=='GSE72094',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_futurama(alpha = 0.7)(12)[7])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="GSE72094")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-GSE72094-signture.pdf',width = 4,height =21)


########
unique(uni$ID)
dd <- uni[uni$ID=='Meta-Cohort',]
dd$ll <- ifelse(dd$CP<0.0001,'****',ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*',''))))

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_jama(alpha = 0.7)(7)[5])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="Meta")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.74,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(limits=c(0.4,0.8),breaks = seq(0.4,0.8,0.2))
ggsave(filename = 'Figure 4/all-Meta-Cohort-signture.pdf',width = 4,height =21)

