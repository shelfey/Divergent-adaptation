####PhyloPCA analysis

library(phytools)

tree <- read.tree("RAxML.tree")
pheno<-read.csv("MEs_traits_use.csv",header = T,sep = ",",row.names = 1)
pca_result<- phyl.pca(tree,pheno, method="BM", mode="corr")
a<-scores(pca_result)
write.table(a,"phylPCA_pheno_raxml.txt",col.names = T,sep = "\t",row.names = T,quote = F)

pheno<-read.csv("948_hubgenes_tpm.txt",header = T,sep = "\t",row.names = 1)
pheno1<-t(pheno)
pca_result<- phyl.pca(tree,pheno1, method="BM", mode="corr")
a<-scores(pca_result)
write.table(a,"phylPCA_hubgene_raxml.txt",col.names = T,sep = "\t",row.names = T,quote = F)


#### vector analysis 
da1<-read.csv("phylPCA_pheno_raxml.txt",header = T,sep = "\t",row.names = 1)
#da1<-read.csv("phylPCA_hubgene_raxml.txt",header = T,sep = "\t",row.names = 1)
head(da1)

da1$group<-as.character(c(rep("pml",6),rep("pmh",8),rep("pr",7),rep("ot",9)))
exp2<-melt(da1,id.vars=c("group"))
md<-summarySE(exp2, measurevar="value", groupvars=c("variable","group"),na.rm=T)


des2 <- subset(md,md$group=="pml")
# array for storing means
n.sim<-5000
means <- sds <- array(NA,dim=c(n.sim,dim(des2)[1]))
colnames(means) <- colnames(sds) <- des2$variable
# start simulation
for(t in 1:n.sim){
  for(s in 1:dim(des2)[1]){
    draws <- rnorm(des2$N[s],des2$value[s],des2$sd[s])
    means[t,s] <- mean(draws)
    sds[t,s] <- sd(draws)
  }}
# store means in data object
observed_pminL <- list(means=des2$value,sds=des2$sd)
simulated_pminL <- list(means=means,sds=sds)



n.sim<-5000
means <- sds <- array(NA,dim=c(n.sim,dim(des2)[1]))
colnames(means) <- colnames(sds) <- des2$variable
# start simulation
for(t in 1:n.sim){
  for(s in 1:dim(des2)[1]){
    draws <- rnorm(des2$N[s],des2$value[s],des2$sd[s])
    means[t,s] <- mean(draws)
    sds[t,s] <- sd(draws)
  }}
# store means in data object
observed_pminL0 <- list(means=des2$value,sds=des2$sd)
simulated_pminL0 <- list(means=means,sds=sds)

des2 <- subset(md,md$group=="pmh")
# array for storing means
n.sim<-5000
means <- sds <- array(NA,dim=c(n.sim,dim(des2)[1]))
colnames(means) <- colnames(sds) <- des2$variable
# start simulation
for(t in 1:n.sim){
  for(s in 1:dim(des2)[1]){
    draws <- rnorm(des2$N[s],des2$value[s],des2$sd[s])
    means[t,s] <- mean(draws)
    sds[t,s] <- sd(draws)
  }}
# store means in data object
observed_pminH <- list(means=des2$value,sds=des2$sd)
simulated_pminH <- list(means=means,sds=sds)


des2 <- subset(md,md$group=="pr")
# array for storing means
n.sim<-5000
means <- sds <- array(NA,dim=c(n.sim,dim(des2)[1]))
colnames(means) <- colnames(sds) <- des2$variable
# start simulation
for(t in 1:n.sim){
  for(s in 1:dim(des2)[1]){
    draws <- rnorm(des2$N[s],des2$value[s],des2$sd[s])
    means[t,s] <- mean(draws)
    sds[t,s] <- sd(draws)
  }}
observed_prinH <- list(means=des2$value,sds=des2$sd)
simulated_prinH <- list(means=means,sds=sds)


des2 <- subset(md,md$group=="ot")
# array for storing means
n.sim<-5000
means <- sds <- array(NA,dim=c(n.sim,dim(des2)[1]))
colnames(means) <- colnames(sds) <- des2$variable
# start simulation
for(t in 1:n.sim){
  for(s in 1:dim(des2)[1]){
    draws <- rnorm(des2$N[s],des2$value[s],des2$sd[s])
    means[t,s] <- mean(draws)
    sds[t,s] <- sd(draws)
  }}
observed_otinH <- list(means=des2$value,sds=des2$sd)
simulated_otinH <- list(means=means,sds=sds) 



# get phenotypic means and the standard deviations
meanpminL <- simulated_pminL$means
meanpminH <- simulated_pminH$means
meanprinH <- simulated_prinH$means
meanotinH <- simulated_otinH$means
meanpminL0 <- simulated_pminL0$means

# Centralize means
zmeanpminL <- (meanpminL-meanpminL)
zmeanpminH <- (meanpminH-meanpminL)
zmeanprinH <- (meanprinH-meanpminL)
zmeanotinH <- (meanotinH-meanpminL)
zmeanpminL0 <- (meanpminL0-meanpminL)

#vector length
length_PBC<-
  length_PBD<-
  length_PCD<-
  length_PB<-
  length_PC<-
  length_PD<-
  length_PM0<-array(NA,dim=c(5000,1))

for (k in 1:1) {
  length_PBC[,k] <- apply(zmeanprinH,1,function(v){sqrt(sum(v^2))})
  #length of plasticity vector A
  length_PBD[,k] <- apply(zmeanotinH,1,function(v){sqrt(sum(v^2))})
  # length of plasticity vector A
  
  length_PB[,k] <- apply(zmeanpminH,1,function(v){sqrt(sum(v^2))})
  #length of plasticity vector A
  length_PC[,k] <- apply(zmeanpminL,1,function(v){sqrt(sum(v^2))})
}  

all<- data.frame(pr=length_PBC,ot=length_PBD,pmh=length_PB,pml=length_PC)

pheno<-melt(all)
colnames(pheno)<-c("group","va")
compaired <- list(c("pml","pr"),c("pmh","pml"),c("ot","pml"))
colnames(pheno)[1]<- "group"
pheno$group<-factor(pheno$group,levels = c("pml","pmh","pr","ot")) 
ggplot(pheno,aes(group,va,fill=group,color=group))+geom_boxplot(width=0.5)+
  stat_compare_means(method = "t.test",label = "P.signif",comparisons = compaired)+
  scale_fill_manual(values = c("#7eaf4a","#f2cc61","#ea8750","#7eaf4a2E"))+
  scale_color_manual(values = c("#f4f2f2","#f4f2f2","#f4f2f2","#f4f2f2"))+
  ylab("Distance")+
  theme(axis.text = element_text(color="black",size = rel(1.2)),axis.title=element_text(size=15))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(color = "black",size=0.5))

write.table(pheno,"pheno_vector_length.txt",col.names = F,sep = "\t",row.names = F,quote = F)

####turkey HSD test for vector length across species
head(pheno)
dim(pheno)
fit <- aov(va ~ group, data = pheno)
tukey_result <- TukeyHSD(fit)
summary(tukey_result)
print(tukey_result)

tukey_df <- as.data.frame(tukey_result$`group`)[c(1:3),]
tukey_df$group1<-rownames(tukey_df)
data<-tukey_df
data$group1<-factor(data$group1,levels = c("pmh-pml","pr-pml","ot-pml"))
ggplot(data, aes(x = group1, y = diff, color = factor(group1))) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +  # 鍒涘缓鐐瑰浘
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),  # 閿欒绾胯〃绀虹疆淇″尯闂?
    position = position_dodge(width = 0.5),  # 璁剧疆浣嶇疆
    width = 0.3  # 閿欒绾垮搴?
  ) +
  labs(
    x = "group",
    y = "Difference of mean",
    title = "Tukey's HSD"
  )+coord_flip()+ylim(0,1)+geom_hline(yintercept = 0,colour="#990000",lwd=1,linetype="dashed")+
  scale_color_manual(values =c("#49912C","#A5CDF8","#1F78B3"))+
  theme(axis.text = element_text(color="black",size = rel(1.2)),axis.title=element_text(size=15))+theme(legend.background = NULL)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = "black",size=0.5))


####sweed windows overlapped with regulatory (key_up2k.pos & key_down2k.pos) and genic regions (key_gb0.pos), these genomic regions can be achieved from gff file.
BEDTools.v2.17.0/bedtools-2.17.0/bin/bedtools intersect -a key_gb0.pos -b OT.sweed2k.result -wb > key.OT.sweed2k.result.gb0
BEDTools.v2.17.0/bedtools-2.17.0/bin/bedtools intersect -a key_up2k.pos -b OT.sweed2k.result -wb > key.OT.sweed2k.result.up2k
BEDTools.v2.17.0/bedtools-2.17.0/bin/bedtools intersect -a key_down2k.pos -b OT.sweed2k.result -wb > key.OT.sweed2k.result.down2k
key<-read.csv("key.OT.sweed2k.result.gb0",header = F,sep = "\t")
#key<-read.csv("key.OT.sweed2k.result.up2k",header = F,sep = "\t")
#key<-read.csv("key.OT.sweed2k.result.down2k",header = F,sep = "\t")
key$length<- abs(key$V2-key$V3)
key1<-subset(key,key$length>500)
write.table(key1,"key.OT.sweed2k.result1.up2k",col.names = F,row.names = F,sep = "\t",quote = F)

#####sweed sampling
fst<-read.csv("OT.sweed2k.result",header = F,sep = "\t")
gb0<-read.csv("key.OT.sweed2k.result1.gb0",header = F, sep = "\t")[,4:7]
colnames(gb0)<-c("gene","V1","V2","V3")
up2k<-read.csv("key.OT.sweed2k.result1.up2k",header = F, sep = "\t")[,c(4,7,8,9)]
colnames(up2k)<-c("gene","V1","V2","V3")
down2k<-read.csv("key.OT.sweed2k.result1.down2k",header = F, sep = "\t")[,c(4,7,8,9)]
colnames(down2k)<-c("gene","V1","V2","V3")


fst<-fst[,c(1,5,6,3)]#reorder results of sweed
fst<-na.omit(fst)
a<-fst[fst$V3 > quantile(fst$V3,0.99),] 

dxy_high<-fst[,3]
dxy_high[dxy_high< quantile(fst$V3,0.99)]<-0
dxy_high[dxy_high >= quantile(fst$V3,0.99)]<-"A"
dxy_high<-as.data.frame(dxy_high)
colnames(dxy_high)<-"fsth"
fst_pich<-cbind(c(1:nrow(fst)),fst,dxy_high) ##annotate the top1 CLR windows 
colnames(fst_pich)[1:5]<-c("pos","V1","V2","V3","V4")

gb0<-merge(gb0,fst_pich, by=c("V1","V2","V3"))
up2k<-merge(up2k,fst_pich, by=c("V1","V2","V3"))
down2k<-merge(down2k,fst_pich, by=c("V1","V2","V3"))


gb0_order<-c()
gb0_p<-c()
up2k_order<-c()
up2k_p<-c()
down2k_order<-c()
down2k_p<-c()
gb0_obj<-c()
up2k_obj<-c()
down2k_obj<-c()

for (i in c("PMO007418")) ###enter gene IDs 
{
  gb0_obj[i]<-mean(subset(gb0,gb0$gene==i)$V4)
  a<- nrow(subset(gb0,gb0$gene==i))
  up2k_obj[i]<-mean(subset(up2k,up2k$gene==i)$V4)
  b<-nrow(subset(up2k,up2k$gene==i))
  down2k_obj[i]<-mean(subset(down2k,down2k$gene==i)$V4)
  c<-nrow(subset(down2k,down2k$gene==i))
  
  gb0_back<-c()  
  up2k_back<-c()
  down2k_back<-c()
  for (j in 1:1000) {
    gb0_back[j]<-mean(fst_pich[sample(1:nrow(fst_pich),a,replace = F),]$V4)
    up2k_back[j]<-mean(fst_pich[sample(1:nrow(fst_pich),b,replace = F),]$V4)
    down2k_back[j]<-mean(fst_pich[sample(1:nrow(fst_pich),c,replace = F),]$V4)
  }
  
  write.csv(gb0_back,paste0(i,"PM_gb0.txt" ))
  write.csv(up2k_back,paste0(i,"PM_up2k.txt" ))
  write.csv(down2k_back,paste0(i,"PM_down2k.txt" ))
  
  gb0_order[i]<-gb0_back[order(gb0_back)][950]
  gb0_p[i]<- sum(gb0_back> gb0_obj[i])/1000
  
  up2k_order[i]<-up2k_back[order(up2k_back)][950]
  up2k_p[i]<- sum(up2k_back> up2k_obj[i])/1000
  
  down2k_order[i]<-down2k_back[order(down2k_back)][950]
  down2k_p[i]<- sum(down2k_back> down2k_obj[i])/1000
}

all<-cbind(gb0_order,gb0_obj,gb0_p,up2k_order,up2k_obj,up2k_p,down2k_order,down2k_obj,down2k_p)
all
write.table(all,"PMH_new_sweed_sampling.txt",col.names = T,sep = "\t", row.names = T,quote = F)




