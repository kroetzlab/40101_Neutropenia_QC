# 40101_Neutropenia_QC

##Code obtained from Chen then altered

setwd("~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/")

library(GenABEL)
library(gdata)
library(SNPassoc)
library(survival)
library(coin)
library(xtable)
options(width=150)


sourcedat="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/L14_Cutoff_Top_090126_FinalReport.csv"
tools::md5sum(sourcedat)

##Call pheno file, (not sure what md5sum means and what it does), and read it in
phenofile="ac_arm_anc_40101jan32014_fmEdits.csv"
tools::md5sum(phenofile)
pheno=read.csv(phenofile)

randompheno="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/pheno 40101 All Patients for Chen 29JUN11.xls"
tools::md5sum(randompheno)
randompheno=read.xls(randompheno)
randomphenoCA=randompheno[randompheno$arm=="CA",]

##Plate design file
designfile="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/40101-link-list.csv"
tools::md5sum(designfile)
link=read.csv(designfile)

cutlocilistfile="Human610-Quadv1_C_ListofCutLoci.txt"
tools::md5sum(cutlocilistfile)


 
##note the tfamfile="40101-all1.tfam" sex coded as 2 instead of "female"
### Had to figure this out but now the convert.snp.tped comand runs

convert.snp.tped(tpedfile="40101-all.tped",
                 tfamfile="40101-all1.tfam",
                 outfile="40101-all-top-geno.raw",
                 bcast = 100000)

GWA40101CA=load.gwaa.data(phe=paste("40101-phdata.txt",sep=""),gen=paste("40101-all-top-geno.raw",sep=""),force=T)

GWA40101CA@gtdata@nsnps
GWA40101CA@gtdata@nids

##this piece of code is to call the expIDs and return the PatIDs (link is the full list of patient samples).
re=link
row.names(re)= re$expid
re["batch1l_plate01_G01",]

##The original code had skip=2 in the snpexc=read.table command but I received an error because there were 4 lines before 
##the data began. This is why I changed it to skip=4.
snpexc=read.table(cutlocilistfile,skip=4)
snpexc=as.character(snpexc[,1])

length(snpexc)
length(GWA40101CA@gtdata@snpnames)
length(setdiff(snpexc,GWA40101CA@gtdata@snpnames)) ##snpexc and GWA40101CA@gtdata@snpnames are in the wrong order, see next line of code
length(setdiff(GWA40101CA@gtdata@snpnames,snpexc)) ##removes the SNPs in the list

GWA40101CA=GWA40101CA[,setdiff(GWA40101CA@gtdata@snpnames,snpexc)]
save(GWA40101CA, file=paste("40101CAALLTOP_942.RData",sep=""))

##check.marker's ibs.exclude arguement: "both", "lower" or "none" – whether both samples with 
##IBS>ibs.threshold should be excluded, the one with lower call rate, or no check 
##(equivalent to use of 'ibs.mrk = -1').
aa=check.marker(GWA40101CA,ibs.exclude="both")
bb=check.marker(GWA40101CA,ibs.exclude="lower")

##The raw data contains \Sexpr{GWA40101CA@gtdata@nids} samples. 
###check.marker function found length(aa$ibsfail)/2 duplicated. Among them, 11 pairs have been 
###reported as the technical replicates while one pair have been found as unintended replicates. 

##1. Twenty-two expids will be excluded as the replicates with the lower callrates.

##2. "batch1l_plate07_F12" will be excluded as the unintended replicates. As design, 
##"batch1l_plate07_F12"==111102, while in the outcome, "batch1l_plate07_F12"=="batch1l_plate01_G01".

##3. Two expids will be excluded as the intended replicates but did not turn out to be 
##duplicates ("batch1l_plate01_G01","batch1l_plate07_D12").

##4. Two pairs expids which are dupliciated as designed but not shown as the duplicated as the 
##outcome

##5. Three CEPH samples will be excluded.

##6. Two male samples will be excluded("batch1l\_plate02\_D11","batch2l\_plate08\_B06").

##7. Patid 111318("batch1l_plate07_G09") consented to 60202, the PG companion, but never 
##registered to 40101 the treatment study, so they are being excluded, as we have no treatment 
##info on them.

##8. Two expids (batch1l_plate04_A09,batch2l_plate08_E04) are excluded from our initinal analysis and these need to be validated.



##analysisgwa
temp=link[link$patid!=""&link$patid!="dH2O",] ##basically a plate list excluding dH2O & NoID
##(978-31=947)
length(unique(temp$patid)) ##Amount of samples that are not duplicates of the 947 (922)

duppatid=as.character(temp[duplicated(temp$patid),"patid"]) ##determines which elements are duplicates calling by patid (25)

dupexpid=as.character(temp$expid[is.element(as.character(temp$patid),duppatid)]) 
## The previous command determined which of the 947(temp$patid) patids are the same as the 
##25(duppatid) returning them as a character vector of expids

NN=temp[is.element(as.character(temp$patid),duppatid),]
## The previous command determined which of the 947(temp$patid) patids are the same as the 
##25(duppatid) returning them as a data frame of characters containing expids and patids.
##So NN is a data frame of duplicates

NN$INIT=ifelse(is.element(as.character(NN$expid),aa$ibsfail),1,0) 
##Creates an additional column INIT in the data frame NN that determined which of the elements
##of NN (by expid) are the same as aa$ibsfail (check.marker output)
##If yes=1 if no=0. No=0 being additional patients (non assigned duplicates) that fail IBS because they are duplicates

NN[order(NN$patid),] ## Orders NN by patid

##aa$ibsfail had 46 (23 pairs) of duplicates patids. 
###Among these 46 expids, one ("batch1l_plate07_F12") is not in the design duplicated files (NN$expid). 
length(aa$ibsfail) ##Amount of patients 
sum(is.element(aa$ibsfail,NN$expid)) ##
aa$ibsfail[!is.element(aa$ibsfail,NN$expid)] 
##Among the 46 (23 pairs) of duplicates patids, 46 expids, one ("batch1l_plate07_F12") is not in the designed duplicated files
##bb$ibsfail have 23 low callrate duplicated expids. All those need to be removed.

length(bb$ibsfail) ##
sum(is.element(bb$ibsfail,NN$expid)) ##

exclude1=bb$ibsfail ##23 dup samples that failed IBS

exclude2="batch1l_plate07_F12" #Among the 46 (23 pairs) of duplicates patids, is not in the design duplicated files

exclude3=c("batch1l_plate01_G01","batch1l_plate07_D12") ##Two expids will be excluded as the intended replicates but did not turn out to be duplicates

exclude4=as.character(NN$expid[!is.element(as.character(NN$expid),aa$ibsfail)]) ##samples that are designed as duplicates but not recognized as such.

exclude5=as.character(temp$expid[substring(temp$patid,1,2)=="CE"|substring(temp$patid,1,2)=="NA"]) #removing CEPH patients

exclude6=c("batch1l_plate02_D11","batch2l_plate08_B06") #male samples, also accounted for in CEPH samples

exclude7=c("batch1l_plate07_G09") ##consented but no treatment info

exclude8=c("batch2l_plate10_C11", "batch1l_plate02_A04") ##genotyped but did not give consent
exc=c(exclude1,exclude2,exclude3,exclude4,exclude5,exclude6,exclude7,exclude8)
GWA40101CAana=GWA40101CA[!is.element(GWA40101CA@phdata$id,exc),]
GWA40101CAana=GWA40101CAana[,setdiff(GWA40101CA@gtdata@snpnames,aa$nocall)] ##removing failed snps
GWA40101CAana@gtdata@nids
GWA40101CAana@gtdata@nsnps
save(GWA40101CAana,file=paste("GWA40101CAana907.RData",sep=""))

##SNP filter
CHECK907=check.marker(GWA40101CAana, callrate=0.99, extr.call=0.99, p.level=1e-08, het.fdr=0)
GWA40101CAana@gtdata@nsnps


##Genetic Euro
GDAT=GWA40101CAana

###Most of this section is handled by PLINK
nonautosnps=GDAT@gtdata@snpnames[is.element(GDAT@gtdata@chromosome,c("23","24","25","26"))]
write.table(nonautosnps,file=paste("nonautosnps.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table(GDAT@phdata[,c("expid","id")],file=paste("expidlist907",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

symsnps=GDAT@gtdata@snpnames[is.element(as.character(GDAT@gtdata@coding),c("CG","GC","AT","TA"))]
write.table(symsnps,file=paste("symetricsnps.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

highldfile="highLDrsidlist.txt"

system(paste("./plink", "--tped", paste("40101-all.tped",sep=""), "--tfam", paste("40101-all.tfam",sep=""),"--missing-genotype 0", 
             "--keep", paste("expidlist907",sep=""),"--exclude", "Human610-Quadv1_C_ListofCutLoci.txt","--make-bed",
             "--out", paste("GWACA907NoLC",sep="")),intern=TRUE)
##Not sure what is different between "Human610-Quadv1_C_ListofCutLoci_Chen.txt" and
##Human610-Quadv1_C_ListofCutLoci.txt" since I don't have the former I just used the latter


system(paste("./plink", "--bfile", paste("GWACA907NoLC",sep=""),"--missing-genotype 0", "--keep", paste("expidlist907",sep=""),
             "--exclude", paste("nonautosnps.txt",sep=""),"--make-bed", "--out", paste("GWACA907NoLCAUTO",sep="")),intern=TRUE)


system(paste("./plink","--bfile", paste("GWACA907NoLCAUTO",sep=""),"--exclude", paste("symetricsnps.txt",sep=""),"--make-bed",
             "--out", paste("GWACA907AUTO",sep="")),intern=TRUE)


system(paste("./plink","--bfile", paste("GWACA907AUTO",sep=""),"--exclude", highldfile,"--make-bed",
             "--out", paste("GWACA907AUTONOLD",sep="")),intern=TRUE)

system(paste("./plink","--bfile", paste("GWACA907AUTONOLD",sep=""),"--geno 0.01","--make-bed",
             "--out", paste("GWACA907AUTONOLDAgrFltr",sep="")),intern=TRUE)
##I had to remove the --maf 0.05 flag to enable the command to run or else I get the error
##Warning message:running command './plink --bfile GWACA907AUTONOLD --maf 0.05 --geno 0.01 
##--make-bed --out GWACA907AUTONOLDAgrFltr' had status 1 
##Because "553932 SNPs failed frequency test ( MAF < 0.05 )"  


system(paste("./plink","--bfile", paste("GWACA907AUTONOLDAgrFltr",sep=""),"--mind 0.021","--make-bed",
             "--out", paste("GWACA907AUTONOLDAgrFltr",sep="")),intern=TRUE)

 
##prune           
system(paste("./plink","--bfile", paste("GWACA907AUTONOLDAgrFltr",sep=""),"--indep-pairwise 1500 150 0.2","--out", 
                 paste("GWACA907AUTONOLDAgrFltrthinning|tee",sep=""), paste("indpairwise.log",sep="")),intern=TRUE)

system(paste("./plink","--bfile", paste("GWACA907AUTONOLDAgrFltr",sep=""),"--extract", paste("GWACA907AUTONOLDAgrFltrthinning.prune.in",sep=""),
             "--make-bed", "--out", paste("GWACA907AUTONOLDAgrFltrThinned",sep="")),intern=TRUE)

##"Reading list of SNPs to extract [ GWACA907AUTONOLDAgrFltrthinning.prune.in ] ... 0 read"
##There aren't any contents in the GWACA907AUTONOLDAgrFltrthinning.prune.in file so this command
##can not run to completion. I am not sure if Chen made this file herself, but I could not progress pass this part


##eigensoftcode
##eigensoft does not run on Mac only in Linux, therefore I did not run the following code myself and it is left unaltered from Chen.
stem=paste("GWACA907AUTONOLDAgrFltrThinned",sep="") #Sets the stemname for .bed, .bim and .fam files
temp=read.table(paste(stem,".fam",sep=""))
temp$V1="CA"
temp$V6=1
write.table(temp,file=paste(stem,".fam",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
altnormstyle="N0"
numoutevec=10 
numoutlieriter=0                #sets maximum number of outlier removal iterations (0 turns off, 5 = original default)
nsnpldregress=0 
noxdata="YES" 
numgamma=10                     #Sets the number of axes for which SNP loadings are calculated (strictly, will be min(numgamma,numoutevec)
numplot=10                      #Sets the number of axes for which plots calculated (strictly, will be min (numgamma,numoutevec,numplot) for Q-Q plot #and min(numoutevec,numplot) for histogram)
gamplot="NO"                    #If "yes" then plot SNP loadings for PC1 against physical position
heatplot="NO"                   #If "yes" then plot "heatmap" of genotypes for each individual
ESOFTdir="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/EIG5.0"   #Sets the location of the EIGENSOFT directory
numoutlierevec=10               #number of principal components along which to remove outliers during each outlier removal iteration
outliersigmathresh=6            #number of standard deviations which an individual must exceed, along one of the top (numoutlierevec) PC's, in order for that individual to be removed as an outlier
kmeans="NO"                     #Option to perform kmeans analysis (for "inversion" regions)
snpweightout="NO"               #Option to call EIGENSOFT's SNP loading option

#args=c("stem=test","numgamma=2")  #use for debuggging
#Override with args if set
t=commandArgs()
if (charmatch("--args",t,nomatch=-1)>=0) args = t[((1:length(t))[t=="--args"]+1):length(t)] else args=""
if (charmatch("stem=",args,nomatch=-1)>=0) stem = strsplit(args[charmatch("stem=",args)],split="=")[[1]][2]
if (charmatch("altnormstyle=",args,nomatch=-1)>=0) altnormstyle = strsplit(args[charmatch("altnormstyle=",args)],split="=")[[1]][2]
if (charmatch("numoutevec=",args,nomatch=-1)>=0) numoutevec = as.numeric(strsplit(args[charmatch("numoutevec=",args)],split="=")[[1]][2])
if (charmatch("numoutlieriter=",args,nomatch=-1)>=0) numoutlieriter = as.numeric(strsplit(args[charmatch("numoutlieriter=",args)],split="=")[[1]][2])
if (charmatch("nsnpldregress=",args,nomatch=-1)>=0) nsnpldregress = as.numeric(strsplit(args[charmatch("nsnpldregress=",args)],split="=")[[1]][2])
if (charmatch("noxdata=",args,nomatch=-1)>=0) noxdata = strsplit(args[charmatch("noxdata=",args)],split="=")[[1]][2]
if (charmatch("numgamma=",args,nomatch=-1)>=0) numgamma = as.numeric(strsplit(args[charmatch("numgamma=",args)],split="=")[[1]][2])
if (charmatch("numplot=",args,nomatch=-1)>=0) numplot = as.numeric(strsplit(args[charmatch("numplot=",args)],split="=")[[1]][2])
if (charmatch("gamplot=",args,nomatch=-1)>=0) gamplot = strsplit(args[charmatch("gamplot=",args)],split="=")[[1]][2]
if (charmatch("heatplot=",args,nomatch=-1)>=0) heatplot = strsplit(args[charmatch("heatplot=",args)],split="=")[[1]][2]
if (charmatch("ESOFTdir=",args,nomatch=-1)>=0) ESOFTdir = strsplit(args[charmatch("ESOFTdir=",args)],split="=")[[1]][2]
if (charmatch("numoutlierevec=",args,nomatch=-1)>=0) numoutlierevec = as.numeric(strsplit(args[charmatch("numoutlierevec=",args)],split="=")[[1]][2])
if (charmatch("numoutlierevec=",args,nomatch=-1)>=0) numoutlierevec = as.numeric(strsplit(args[charmatch("numoutlierevec=",args)],split="=")[[1]][2])
if (charmatch("outliersigmathresh=",args,nomatch=-1)>=0) outliersigmathresh = as.numeric(strsplit(args[charmatch("outliersigmathresh=",args)],split="=")[[1]][2])
if (charmatch("kmeans=",args,nomatch=-1)>=0) kmeans = strsplit(args[charmatch("kmeans=",args)],split="=")[[1]][2]
if (charmatch("snpweightout=",args,nomatch=-1)>=0) snpweightout = strsplit(args[charmatch("snpweightout=",args)],split="=")[[1]][2]


##Save output also to text file
sink(file=paste(stem,".Rout",sep=""),type="output",split=TRUE)


##Make copies of .bim and .fam files, Create .par file with same stem name
print("Reading arguments...")
system(paste("cp ",stem,".bim ",stem,".pedsnp",sep=""))
system(paste("cp ",stem,".fam ",stem,".pedind",sep=""))
FID = file(paste(stem,".par",sep=""),"w")
writeLines(paste("genotypename:    ",stem,".bed",sep=""),FID)
writeLines(paste("snpname:         ",stem,".pedsnp",sep=""),FID)
writeLines(paste("indivname:       ",stem,".pedind",sep=""),FID)
writeLines(paste("evecoutname:     ",stem,".evec",sep=""),FID)
writeLines(paste("evaloutname:     ",stem,".eval",sep=""),FID)
writeLines(paste("altnormstyle:    ",altnormstyle,sep=""),FID)
writeLines(paste("numoutevec:      ",as.character(numoutevec),sep=""),FID)
writeLines(paste("numoutlieriter:  ",as.character(numoutlieriter),sep=""),FID)
writeLines(paste("nsnpldregress:   ",as.character(nsnpldregress),sep=""),FID)
writeLines(paste("noxdata:         ",noxdata,sep=""),FID)
writeLines(paste("numoutlierevec:  ",as.character(numoutlierevec),sep=""),FID)
writeLines(paste("outliersigmathresh:  ",as.character(outliersigmathresh),sep=""),FID)
writeLines(paste("outlieroutname:  ",stem,".outliers",sep=""),FID)
writeLines(paste("phylipoutname:   ",stem,".fst",sep=""),FID)
if (snpweightout=="YES") writeLines(paste("snpweightoutname:   ",stem,".load",sep=""),FID)
close(FID)


##Call EIGENSOFT
##print("Calling EIGENSOFT...")
system(paste("nice ", ESOFTdir, "/bin/smartpca -p ", stem, ".par > ", stem, ".Sout", sep=""))


##Call TWstats: twstats program to compute Tracy-Widom statistics (statistical significance of each principal component
##print("Calling TWstats...")
##system(paste("nice ", ESOFTdir, "/bin/twstats -t ", ESOFTdir, "/POPGEN/twtable -i ", stem, ".eval >>  ", stem, ".Sout", sep=""))


##Read in values from .evec file
pcafile = paste(stem,".evec",sep="")
FIDpca = file(pcafile,"r")
lambda = strsplit(readLines(FIDpca,n=1),split=" +")[[1]] [c(-1,-2)]  #Note " +" means match any number of spaces
NNaxes = length(lambda)
eigvec = as.matrix( read.table(FIDpca,header=FALSE,row.names=1,comment.char="",colClasses=c("character",rep("numeric",NNaxes),"NULL")) )
close(FIDpca)


##Get Genetic Euro
dim(eigvec)
eigvec=as.data.frame(eigvec)
eigvecfile=data.frame(id=substring(rownames(eigvec),4,22),PC1=eigvec$V2,PC2=eigvec$V3,PC3=eigvec$V4)

phenopat=randomphenoCA[,c("patid","expid","race","ethnic")]
names(phenopat)[2]="expidDee"
link$patid=as.character(link$patid)
phenopat=merge(phenopat,link,by.x="patid",by.y="patid",all.x=TRUE)
phenopat=phenopat[!is.na(phenopat$expid),]
names(phenopat)[which(names(phenopat)=="expid")]="id"

GWA40101CAana=add.phdata(GWA40101CAana,phenopat)
GWA40101CAana=add.phdata(GWA40101CAana,eigvecfile)

GDAT=GWA40101CAana

pc1mean=mean(GDAT@phdata$PC1[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])
pc1std=sd(GDAT@phdata$PC1[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])

pc2mean=mean(GDAT@phdata$PC2[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])
pc2std=sd(GDAT@phdata$PC2[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])

pc3mean=mean(GDAT@phdata$PC3[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])
pc3std=sd(GDAT@phdata$PC3[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])

ii1=GDAT@phdata$PC1<=(pc1mean+2*pc1std)&GDAT@phdata$PC1>=(pc1mean-2*pc1std)
ii2=GDAT@phdata$PC2<=(pc2mean+2*pc2std)&GDAT@phdata$PC2>=(pc2mean-2*pc2std)
ii3=GDAT@phdata$PC3<=(pc3mean+2*pc3std)&GDAT@phdata$PC3>=(pc3mean-2*pc3std)
#ii4=gwa_auto@phdata$race=="White"|gwa_auto@phdata$race=="Unknown"

GDAT@phdata$genetic_euro=ifelse(ii1&ii2&ii3,1,0)
GDAT@phdata$race=as.character(GDAT@phdata$race)
GDAT@phdata$ethnic=as.character(GDAT@phdata$ethnic)
table(GDAT@phdata$genetic_euro)
GWA40101CACAUC=GDAT[GDAT@phdata$genetic_euro==1,]

save(GDAT,file=paste("GWA40101CACAUCana907.RData",sep=""))



##PCAplot
pdf("/data/CALGB/40101/CA/PROC/PCA/40101CAPC1vs2.pdf",width=8,height=8)
plot(GDAT@phdata$PC1,GDAT@phdata$PC2,col=c(1:6)[factor(GDAT@phdata$race)],pch=c(3,19)[factor(GDAT@phdata$genetic_euro)],xlab="",ylab="",cex=c(0.5,0.25)[factor(GDAT@phdata$genetic_euro)])
legend(0.08,0.2,legend=levels(factor(GDAT@phdata$race)),col=c(1:6),pch=19,cex=0.75)
legend(0.08,0.08,legend=c("Cauc","Non-Cauc"),col=1,pch=c(19,3),cex=0.75)
graphics.off()

pdf("/data/CALGB/40101/CA/PROC/PCA/40101CAPC2vs3.pdf",width=8,height=8)
plot(GDAT@phdata$PC2,GDAT@phdata$PC3,col=c(1:6)[factor(GDAT@phdata$race)],pch=c(3,19)[factor(GDAT@phdata$genetic_euro)],xlab="",ylab="",cex=c(0.5,0.25)[factor(GDAT@phdata$genetic_euro)])
legend(0.13,0.1,legend=levels(factor(GDAT@phdata$race)),col=c(1:6),pch=19,cex=0.75)
legend(0.13,0.0,legend=c("Cauc","Non-Cauc"),col=1,pch=c(19,3),cex=0.75)
graphics.off()

pdf("/data/CALGB/40101/CA/PROC/PCA/40101CAPC1vs3.pdf",width=8,height=8)
plot(GDAT@phdata$PC1,GDAT@phdata$PC3,col=c(1:6)[factor(GDAT@phdata$race)],pch=c(3,19)[factor(GDAT@phdata$genetic_euro)],xlab="",ylab="",cex=c(0.5,0.25)[factor(GDAT@phdata$genetic_euro)])
legend(0.08,0.2,legend=levels(factor(GDAT@phdata$race)),col=c(1:6),pch=19,cex=0.75)
legend(0.08,0.08,legend=c("Cauc","Non-Cauc"),col=1,pch=c(19,3),cex=0.75)
graphics.off()

save(GWA40101CACAUC,file=paste("GWA40101CACAUC757.RData",sep=""))


##PCAplot
    
GDAT@gtdata@nids
###907

GDAT@gtdata@nsnps
###588426
pruneinSNP=scan("GWACA907AUTONOLDAgrFltrthinning.prune.in", what="character")

##GDATauto=GDAT[,!is.element(GDAT@gtdata@chromosome,c("23","24","25","26"))]
GDATauto=GDAT[,pruneinSNP]
GDATauto@gtdata@nsnps
###573058

#Using R code to compute Eigenstrat
##eigenstrat.1<-function(geno){                 #ind x snp matrix of genotypes \in 0,1,2
##  geno<-t(geno)
##  nMis<-rowSums(is.na(geno))
##  geno<-geno[nMis==0,]                      #remove snps with missing data
##  avg<-rowSums(geno)/ncol(geno)             # get allele frequency times 2
##  keep<-avg!=0&avg!=2                       # remove sites with non-polymorphic data
##  avg<-avg[keep]
##  geno<-geno[keep,]
##  snp<-nrow(geno)                           #number of snps used in analysis
##  ind<-ncol(geno)                           #number of individuals used in analuysis
##  freq<-avg/2                               #frequency
##  M <- (geno-avg)/sqrt(freq*(1-freq))       #normalize the genotype matrix
##  X<-t(M)%*%M                               #get the (almost) covariance matrix
##  X<-X/(sum(diag(X))/(snp-1))
##  E<-eigen(X)
##  class(E)<-"eigenstrat"
##  return(E)


eigenstrat.1<-function(geno){                 #ind x snp matrix of genotypes \in 0,1,2
    geno<-t(geno)                               #snp x ind matrix of genotypes \in 0,1,2
    avg<-rowMeans(geno, na.rm=TRUE)
    #nMis<-rowSums(is.na(geno))
    #geno<-geno[nMis==0,]                      #remove snps with missing data
    #avg<-rowSums(geno)/ncol(geno)             # get allele frequency times 2
    keep<-avg!=0&avg!=2                       # remove sites with non-polymorphic data
    avg<-avg[keep]
    geno<-geno[keep,]
    snp<-nrow(geno)                           #number of snps used in analysis
    ind<-ncol(geno)                           #number of individuals used in analuysis
    freq<-avg/2                               #frequency
    M<-(geno-avg)/sqrt(freq*(1-freq))       #normalize the genotype matrix
    M[is.na(M)]<-0
    X<-t(M)%*%M                               #get the (almost) covariance matrix
    X<-X/(sum(diag(X))/(snp-1))
    E<-eigen(X)
    class(E)<-"eigenstrat"
    return(E)
}

tempdat=as.numeric(GDATauto@gtdata)
eiganval=eigenstrat.1(tempdat)
tt=eiganval[[2]]

GDATauto=add.phdata(GDATauto,data.frame(id=GDATauto@gtdata@idnames,IBSPC1=tt[,1],IBSPC2=tt[,2],IBSPC3=tt[,3]))

pc1mean=mean(GDATauto@phdata$IBSPC1[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])
pc1std=sd(GDATauto@phdata$IBSPC1[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])

pc2mean=mean(GDATauto@phdata$IBSPC2[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])
pc2std=sd(GDATauto@phdata$IBSPC2[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])

pc3mean=mean(GDATauto@phdata$IBSPC3[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])
pc3std=sd(GDATauto@phdata$IBSPC3[GDATauto@phdata$race=="White"&GDATauto@phdata$ethnic=="Non-Hispanic"])

ii1=GDATauto@phdata$IBSPC1<=(pc1mean+2*pc1std)&GDATauto@phdata$IBSPC1>=(pc1mean-2*pc1std)
ii2=GDATauto@phdata$IBSPC2<=(pc2mean+2*pc2std)&GDATauto@phdata$IBSPC2>=(pc2mean-2*pc2std)
ii3=GDATauto@phdata$IBSPC3<=(pc3mean+2*pc3std)&GDATauto@phdata$IBSPC3>=(pc3mean-2*pc3std)
#ii4=gwa_auto@phdata$race=="White"|gwa_auto@phdata$race=="Unknown"

GDATauto@phdata$genetic_euro=ifelse(ii1&ii2&ii3,1,0)
GDATauto@phdata$race=as.character(GDATauto@phdata$race)
GDATauto@phdata$ethnic=as.character(GDATauto@phdata$ethnic)
table(GDATauto@phdata$genetic_euro)
GWA40101CACAUC=GDATauto[GDATauto@phdata$genetic_euro==1,]

#Use IBS function
IBS=ibs(GDAT,weight="freq")
diag(IBS)<-hom(GDAT)$Var
MDS=cmdscale(as.dist(1-IBS))
GDAT=add.phdata(GDAT,data.frame(id=rownames(MDS),IBSPC1=MDS[,1],IBSPC2=MDS[,2]))

pc1mean=mean(GDAT@phdata$IBSPC1[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])
pc1std=sd(GDAT@phdata$IBSPC1[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])

pc2mean=mean(GDAT@phdata$IBSPC2[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])
pc2std=sd(GDAT@phdata$IBSPC2[GDAT@phdata$race=="White"&GDAT@phdata$ethnic=="Non-Hispanic"])

ii1=GDAT@phdata$IBSPC1<=(pc1mean+2*pc1std)&GDAT@phdata$PC1>=(pc1mean-2*pc1std)
ii2=GDAT@phdata$IBSPC2<=(pc2mean+2*pc2std)&GDAT@phdata$PC2>=(pc2mean-2*pc2std)

GDAT@phdata$IBSgenetic_euro=ifelse(ii1&ii2,1,0)

table(GDAT@phdata$genetic_euro,GDAT@phdata$IBSgenetic_euro,exclude=NULL)


##Cauc SNP filter
    
CHECK=check.marker(GWA40101CACAUC, callrate=0.99, extr.call=0.99, p.level=1e-08, het.fdr=0)
GWA40101CACAUC@gtdata@nsnps

GDATCAUCreduced=GWA40101CACAUC[,CHECK$snpok]
GDATCAUCreduced@gtdata@nsnps

GDATCAUCreducedauto=GDATCAUCreduced[,!is.element(GDATCAUCreduced@gtdata@chromosome,c("23","24","25","26"))]
GDATCAUCreducedauto@gtdata@nsnps

save(GDATCAUCreducedauto, file=paste(PROCSTEM,"GWA40101CACAUC757reducedauto.RData",sep=""))
