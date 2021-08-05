#code for analysis of SNP loci (from filtered VCF file)

library("adegenet")
library("hierfstat")
library("pegas")
library(vcfR)
library(reshape2)
library(ggplot2)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(dplyr)
library(diveRsity)
library("assignPOP")

#roseate tern color palette 
myCol=c("darkorange1","pink", "lightgrey","lightskyblue", "black")
#################################################################################################

#adegenet: https://grunwaldlab.github.io/Population_Genetics_in_R/reading_vcf.html

#input files and make genlight for analysis 

#convert vcf to genlight via https://knausb.github.io/vcfR_documentation/export_genlight_snpclone.html
#input vcf file
vcf<- read.vcfR("modern_SNPs_raw_filter_90_thinned.recode.vcf")

#input structure file with hwe loci removed (71 genotypes, 2043 loci)
rost_gen<-read.structure("modern_SNPs_raw_filter_90_thinned_structure.stru")

#convert to genlight format
#rost_snp <- vcfR2genlight(vcf)
#to genind format
rost_gen<-vcfR2genind(vcf)
#look at matrix
(as.matrix(rost_gen))[c(1:5), 1:12]

#read in population file 
snp_pop=read.csv("snp_pops.csv")
#assign population, ie put it in the population slot 
#pop(rost_snp)=snp_pop$Pop
pop(rost_gen)=snp_pop$Pop
#check 
popNames(rost_gen)
#make sure is diploid
ploidy(rost_gen) <- 2

###############################################################################
#HWE test

#chi squared test of HWE, 2 p-values: 1 analytical, 1 from permutations
(snps_hwe<-hw.test(rost_gen, B=1000))
write.csv(x=snps_hwe, file = "snps_hwe_thinned")

#check for each population
(snps_hwe_pop<-seppop(rost_gen) %>% lapply (hw.test, B = 1000))

#just extract p-values 
(snps_hwe_pop_mat <- sapply(snps_hwe_pop, "[", i = TRUE, j = 3))
write.csv(x=snps_hwe_pop_mat, file = "snps_hwe_thinned")

###############################################################################

#population genetics
#https://popgen.nescent.org/StartSNP.html

#create hierfstat object
rost2 <- genind2hierfstat(rost_gen)

#genetic diversity 
div <- summary(rost_gen)
div
#Number of alleles per group: 3354 4609 4669 4850
#Percentage of missing data: 3.62 %
#plot observed het
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
#plot observed vs expected het
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

#bartlett test of obv vs expected
bartlett.test(list(div$Hexp, div$Hobs))
#p=0.15, no diff btwn obvs, exp

#W&C FST and FIS 
wc(rost_gen)

#pairwise Fst
genet.dist(rost_gen, method = "WC84")

##################################################################

rostpop="modern_SNPs_raw_filter_90_thinned_genepop.txt"

#diveRsity
#basic statistics (allelic richness etc using either rarefaction or bootstrapping)
#fis: inbreeding coefficient
#ar: allelic richness

basicStats(infile = "modern_SNPs_raw_filter_90_thinned_genepop.txt", outfile = "modern_SNPs_raw_filter_90_thinned_genepop_output_rare", fis_ci = TRUE,
           ar_ci = TRUE, fis_boots = 1000, ar_boots = 1000,
           mc_reps = 1000, rarefaction = TRUE, ar_alpha = 0.05,fis_alpha = 0.05)

#differentiation
fastDivPart(infile = "modern_SNPs_raw_filter_90_thinned_genepop.txt", outfile = "ROST_divPart", gp = 3,pairwise = TRUE, fst = TRUE,bs_locus = FALSE,
        bs_pairwise = TRUE,boots = 1000, plot = FALSE,para = FALSE)

##################################################################
#read in stats

#anova to compare pop genetic parameters
compare=read.csv("anova.csv")

#allelic richness
m = aov(ar~pop, data=compare)
summary(m)
AR=TukeyHSD(m)
# Tuckey test representation :
plot(AR, las=1)

#observed heterozygosity
m1 = aov(obs_het~pop, data=compare)
summary(m1)
ho=TukeyHSD(m1)
# Tuckey test representation :
plot(ho, las=1)

#expected heterozygosity
m2 = aov(exp_het~pop, data=compare)
summary(m2)
he=TukeyHSD(m2)
# Tuckey test representation :
plot(he, las=1)


#################################################################

#DACP

#identify number of groups with k means 
grp=find.clusters(rost_gen, max.n.clust=10)
#biggest increase after n=2
#look at group assignment 
table(pop(rost_gen), grp$grp)
#MASS separate from other 2 

#cross-validation to determine number of principal components to use
#k-means: use as many PCs as need, DAPC: want to minimize PCs, overfitting bad 
set.seed(999) #15 pcs
rostx <- xvalDapc(tab(rost_gen, NA.method = "mean"), pop(rost_gen))
#check results 
rostx[-1]

set.seed(999)
system.time(rostx <- xvalDapc(tab(rost_snp, NA.method = "mean"), pop(rost_snp),
                              n.pca = 10:20, n.rep = 1000,
                              parallel = "multicore", ncpus = 4L))

#thinned set: 15 PCAs 
#n.da is number of populations - 1 
#group membership without a priori populations: 
dapc1 <- dapc(rost_gen, var.contrib = TRUE, n.pca=15, n.da=1, grp$grp)

scatter(dapc1,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

#print contents of the object
print.dapc(dapc1)
#summary/useful info 
summary.dapc(dapc1)
#predict indiviudal assignment to groups 
predict.dapc(dapc1)

#no prior, populations from clusters
dapc2 <- dapc(rost_gen, var.contrib = TRUE, n.pca=15, n.da=3)

#print contents of the object
print.dapc(dapc2)
#summary/useful info 
summary.dapc(dapc2)
#predict indiviudal assignment to groups 
predict.dapc(dapc2)

scatter(dapc2,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

#admixture plot
compoplot(dapc2, col=myCol,lab="", ncol=2)

#plot of population assignment 
assignplot(dapc2, subset=1:71)

loadingplot(dapc1$var.contr, threshold=quantile(dapc1$var.contr,0.75))
############################################################################################
#cross-validate population assignmnet: DAPC with training population 
#training pop: 25% of each population but 50% for AZ

#subset object: 

train=rost_gen[c(2, 3, 5:8, 18:21, 41:45),]
#check samples
(as.matrix(train))[c(1:12), 1:3]
set.seed(999) #2
trainx <- xvalDapc(tab(train, NA.method = "mean"), pop(train))
#check results 
trainx[-1]
#2 PCs

new=rost_gen[c(3:6, 9:17, 22:40, 46:71)]

dapc3 <- dapc(train, var.contrib = TRUE, n.pca=2, n.da=3)
pred.assign<-predict.dapc(dapc3, newdata=new)
print.dapc(dapc2)
#summary/useful info 
summary.dapc(dapc2)

#posterior membership probabilities 
pred.assign$posterior

#mean assigned accurately
mean(as.character(pred.assign$assign)==as.character(pop(new)))

table.value(table(pred.assign$assign, pop(new)), col.lab=levels(pop(new)))

col <- rainbow(length(levels(pop(train))))
col.points <- transp(col[as.integer(pop(train))],.2)
scatter(dapc3, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, legend=TRUE)
par(xpd=TRUE)
points(dapc3$ind.coord[,1], dapc3$ind.coord[,2], pch=20,
       col=col.points, cex=5)
col.sup <- col[as.integer(pop(new))]
points(pred.assign$ind.scores[,1], pred.assign$ind.scores[,2], 
       pch=15,col=transp(col.sup,.7), cex=2)
#################################################################

#assignpop

pop <- read.Structure( "modern_SNPs_raw_filter_90_thinned_structure.txt")

#MCMC validation: 
assign.MC( pop, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.5, 0.75, 1),
           loci.sample="fst", iterations=100, model="svm", dir="Result-folder_assignpop/")
#calculate MCMC accuracy 
accuMC <- accuracy.MC(dir = "Result-folder_assignpop/")

#assignment matrix
assign.matrix( dir="Result-folder_assignpop/")

#plot accuracy 
accuracy.plot(accuMC, pop = "all")

#####################

#k-fold validation 
assign.kfold(pop, k.fold=c(2, 3, 4), train.loci=c(0.1, 0.25, 0.5, 1), 
              loci.sample="random", model="lda", dir="Result-folder2/")

#k-folds accuracy:
accuKF <- accuracy.kfold(dir = "Result-folder2/")

#1: AZ, 2: FL. 3: MASS, 4: USVI
assign.matrix( dir="Result-folder2/")

#membership plot 
membership.plot(dir = "Result-folder2/")

#################################################################
#write genind to structure
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)

genind2structure <- function(obj, file="", pops=FALSE){
        if(!"genind" %in% class(obj)){
                warning("Function was designed for genind objects.")
        }
        
        # get the max ploidy of the dataset
        pl <- max(obj@ploidy)
        # get the number of individuals
        S <- adegenet::nInd(obj)
        # column of individual names to write; set up data.frame
        tab <- data.frame(ind=rep(indNames(obj), each=pl))
        # column of pop ids to write
        if(pops){
                popnums <- 1:adegenet::nPop(obj)
                names(popnums) <- as.character(unique(adegenet::pop(obj)))
                popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
                tab <- cbind(tab, data.frame(pop=popcol))
        }
        loci <- adegenet::locNames(obj) 
        # add columns for genotypes
        tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                                 dimnames=list(NULL,loci)))
        
        # begin going through loci
        for(L in loci){
                thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                                          dimnames(obj@tab)[[2]]), 
                                    drop = FALSE] # genotypes by locus
                al <- 1:dim(thesegen)[2] # numbered alleles
                for(s in 1:S){
                        if(all(!is.na(thesegen[s,]))){
                                tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
                                tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
                                tab[tabrows,L] <- rep(al, times = thesegen[s,])
                        }
                }
        }
        
        # export table
        write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

