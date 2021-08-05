#code for analysis of microsatellite loci (from file formatted for STRUCTURE)

library(pegas)
library(hierfstat)
library(xlsx)
library(sendplot)
library("mmod")
library("reshape2")
library("ggplot2")
library(gtools)
library(car)
library(gridExtra)
library(grid)
library(diveRsity)
library(assignPOP)

#roseate tern color palette
myCol=c("darkorange1","pink", "lightgrey","lightskyblue", "black")
#################################################################################################
#pop analyses: adegenet (https://popgen.nescent.org/startMicrosatellite.html)

#input: converted -9s to 0s &: genotypes: 94; markers: 12
#1: NE, 2: VI, 3: FL, 4: AZ
rost <- read.structure("ROST_structure.str")
summary(rost)

#set pops for adegenet
rost_strat<-read.table("pops.csv", header=TRUE)
strata(rost)=rost_strat 
setPop(rost)=~ï..pop
#why coded as 0s and 1s?!? 
(as.matrix(rost))[c(1:30), 1:12]

#test for HWE
#b=0 for rarefaction 
hw.test(rost, B = 10000)

#hweby population 
lapply(seppop(rost), hw.test)

##############################################################

#diveRsity

#basic statistics (allelic richness etc using either rarefaction or bootstrapping)
#fis: inbreeding coefficient
#ar: allelic richness
#genepop after removing 3 crappy fl samples 
basicStats(infile = "ROST_genepop.txt", outfile = "ROST_msats_basicstats", fis_ci = TRUE,
           ar_ci = TRUE, fis_boots = 10000, ar_boots = 10000,
           mc_reps = 10000, rarefaction = FALSE, ar_alpha = 0.05,fis_alpha = 0.05)

#anova to compare pop genetic parameters
compare=read.csv("anova_ar.csv")

#model blocked by locus
#allelic richness
m2 = aov(fis~ï..site+locus, data=compare)
summary(m2)
#not significantly different

#observed heterozygosity
m3 = aov(obs_het~ï..site+locus, data=compare)
summary(m3)
#loci are different (would expect this) but not site

#expected heterozygosity
m4 = aov(exp_het~ï..site+locus, data=compare)
summary(m4)
#loci are different (would expect this) but not site

AZ=compare %>%
  filter(ï..site=="AZ")
FL=compare %>%
  filter(ï..site=="FL")
NE=compare %>%
  filter(ï..site=="NE")
VI=compare %>%
  filter(ï..site=="VI")

#bartlett test of obv vs expected
bartlett.test(list(VI$obs_het, VI$exp_het))
bartlett.test(list(NE$obs_het, NE$exp_het))
bartlett.test(list(FL$obs_het, FL$exp_het))
bartlett.test(list(AZ$obs_het, AZ$exp_het))

#####################################################################################################################################
#differentiation 
#analysis in Adegenet from https://popgen.nescent.org/2015-12-15-microsatellite-differentiation.html 

#get confidence intervals with bootstrapping 
set.seed(20151219) # Be sure to set a seed for any random analysis!
bs_reps <- chao_bootstrap(rost, nreps = 1000)
summarise_bootstrap(bs_reps, D_Jost) # Using the D_Jost function to summarize.

#AMOVA
rost_dist=dist(rost)
rost_stra=strata(rost)
rost_amova=pegas::amova(rost_dist ~ ï..pop, data=rost_stra, nperm=10)
rost_amova
###########################################################################################

#from http://adegenet.r-forge.r-project.org/files/tutorial.pdf
#PCA file: NAs (ie missing data) replaced with mean allele frequency from each locus bc NAs not allowed 
rost_pca <- read.structure("ROST_pca.str")

#look at eigenvalues to pick best # PCs
rost.pca <- glPca(rost_pca, nf = 3)
barplot(100*rost.pca$eig/sum(rost.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#classic PCA
pca1 <- dudi.pca(rost_pca$tab, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "Eigenvalues")
s.class(pca1$li, rost_pca$pop, sub = "PCA 1-2", csub = 2)
add.scatter.eig(pca1$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "top")

#############################################################################

#DACP: https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
##https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf

#identify number of groups with k means 
#https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
grp=find.clusters(rost, max.n.clust=10)
#look at group assignment 
table(pop(rost), grp$grp)

#cross-validation to determine number of principal components to use 
set.seed(999)
#no need to remove NAs, automatically converts to mean value for you
rostx <- xvalDapc(tab(rost, NA.method = "mean"), grp$grp)
#check results 
#look for # of PCs that gives highest % of correctly predicted subsamples with lowest error
rostx[-1]
#5

#5: optimal # of PCs
#without a priori population assignment (by cluster)
dapc1 <- dapc(rost, n.pca=5, n.da=1, grp$grp)
#print contents of the object
print.dapc(dapc1)
#summary/useful info 
summary.dapc(dapc1)

scatter(dapc1,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

#with a priori population grouping

dapc2 <- dapc(rost, n.pca=5, n.da=3, pop(rost))
scatter(dapc2,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

summary.dapc(dapc2)

compoplot(dapc2, lab="", ncol=2)

#predict individual assigment to group 
predict.dapc(dapc2)

#plot of population assignment 
assignplot(dapc2)

#####################################################################################################################################

#cross-validate population assignmnet: DAPC with training population 
#training pop: 25% of each population

#subset object: 

train=rost_pca[c(2, 3, 14:17, 33:37, 77:82),]
#check samples
(as.matrix(train))[c(1:12), 1:3]
set.seed(999) #2
trainx <- xvalDapc(tab(train, NA.method = "mean"), pop(train))
#check results 
trainx[-1]
#3 PCs

new=rost_pca[c(1, 4:13, 18:32, 38:76, 83:94)]

dapc3 <- dapc(train, var.contrib = TRUE, n.pca=3, n.da=3)
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

assignpop

pop <- read.Structure( "ROST_pca.str")

#MCMC validation: 
assign.MC( pop, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.5, 0.75, 1),
           loci.sample="fst", iterations=100, model="svm", dir="Result-folder_assignpop/")
#calculate MCMC accuracy 
accuMC <- accuracy.MC(dir = "Result-folder_assignpop/")

#assignment matrix
assign.matrix( dir="Result-folder_assignpop/")

#plot accuracy 
accuracy.plot(accuMC, pop = "all")
