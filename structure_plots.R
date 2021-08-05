#code for creating plots from STRUCTURE ouput files for both SNP and microsatellite loci 

library(pophelper)
library(ggplot2)
library(gridExtra)

#STRUCTURE RUNS: SNPS: MODERN: 90% COMPLETE: THINNED (1/UCE)

#code from http://www.royfrancis.com/pophelper/articles/index.html

#input "results" files interactively as list 
slist=readQ(files=choose.files(multi=TRUE),filetype="structure")
# check class of ouput
class(slist)
# view head of first converted file
head(slist[[1]])

#check attributes
attributes(slist)

#check tabulated list
head(tabulateQ(slist))

#summarized table of structure runs
head(summariseQ(tabulateQ(slist)))

#evanno method: output likelihood etc to pick k 
em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))
em 

#add in labels 
#import population file 
#snp_pops=read.csv("snp_pops.csv")
snp_pops=read.csv("snp_pops.csv")

# length of labels equal to number of individuals?
#nrow(snp_pops)
nrow(snp_pops)

# check if labels are a character data type
#sapply(snp_pops, is.character)
sapply(snp_pops, is.character)

#use a single label set for just population 
onelabset1 <- snp_pops[,2,drop=FALSE]
head(onelabset1)

#code to plot runs separately 
p1 <- plotQ(slist[19],returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset1,subsetgrp=c("Massachusetts", "Azores", "Florida", "Virgin Islands"), grplabsize=4,linesize=0.8,pointsize=4,sortind="all")

p2 <- plotQ(slist[c(19,21)], imgoutput="join", returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset1,subsetgrp=c("Massachusetts", "Azores","Florida", "Virgin Islands"),grplabsize=4,linesize=0.8,pointsize=4)

grid.arrange(p1$plot[[1]],nrow=1)

#FULL FIGURE: all iterations 2-5 merged and aligned 
#merge q: merge together runs of the same k (11:40) to omit k = 1 and k = 5
slist=readQ(files=choose.files(multi=TRUE),filetype="structure")
slist_1 <- alignK(slist)
slist_2 <- mergeQ(slist_1)

p2 <- plotQ(slist_2, imgoutput="join", returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset1,subsetgrp=c("Massachusetts","Florida",  "Virgin Islands", "Azores"),grplabsize=6,
            linesize=0.8,pointsize=4, splab=c("K = 2","K = 3","K = 4"), splabsize=12, 
            clustercol=c("lightgrey","lightskyblue","pink", "orange"), showtitle=F, titlelab = "Roseate Tern Population Structure", titlesize=20,
            showsubtitle=F, subtitlelab = "Population structure from SNPs", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)
        
grid.arrange(p2$plot[[1]],nrow=1)

#SINGLE FIGURE:  all runs of k=2 merged and aligned: 11-20
k_2_snp=readQ(files=choose.files(multi=TRUE),filetype="structure")
k2_snp_list <- alignK(k_2_snp)
k2_snp_merge <- mergeQ(k2_snp_list)

k2_snp <- plotQ(k2_snp_merge,returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset1,subsetgrp=c("Massachusetts", "Florida", "Virgin Islands", "Azores"),grplabsize=6,
            linesize=0.8,pointsize=4, splab="K = 2", splabsize=12, 
            clustercol=c("lightgrey","lightskyblue"), showtitle=F, titlelab = "Roseate Tern Population Structure", titlesize=20,
            showsubtitle=T, subtitlelab = "Structure from SNPs", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)
grid.arrange(k2_snp$plot[[1]],nrow=1)

#################################################################################################################################

#STRUCTURE RUNS: MSATS: MODERN

#code from http://www.royfrancis.com/pophelper/articles/index.html

#input files interactively as list 
mlist=readQ(files=choose.files(multi=TRUE),filetype="structure")
# check class of ouput
class(mlist)
# view head of first converted file
head(mlist[[1]])

#check attributes
attributes(mlist)

#check tabulated list
head(tabulateQ(mlist))

#summarized table of structure runs
head(summariseQ(tabulateQ(mlist)))

#evanno method: output likelihood etc to pick k 
em <- evannoMethodStructure(summariseQ(tabulateQ(mlist)))
em 

#add in labels 
#import population file 
m_pops=read.csv("msat_pops.csv")
# length of labels equal to number of individuals?
nrow(m_pops)
# check if labels are a character data type
sapply(m_pops, is.character)
#use a single label set for just population 
onelabset2 <- m_pops[,2,drop=FALSE]
head(onelabset2)

#FULL FIGURE: all iterations 2-5 merged and aligned 
#merge q: merge together runs of the same k (11:40)
mlist=readQ(files=choose.files(multi=TRUE),filetype="structure")
mlist_1 <- alignK(mlist)
mlist_2 <- mergeQ(mlist_1)

m2 <- plotQ(mlist_2, imgoutput="join", returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset2,subsetgrp=c("Massachusetts","Florida",  "Virgin Islands", "Azores"),grplabsize=6,
            linesize=0.8,pointsize=4, splab=c("K = 2","K = 3","K = 4"), splabsize=12, 
            clustercol=c("lightskyblue","lightgrey","pink", "orange"), showtitle=F, titlelab = "Roseate Tern Population Structure", titlesize=20,
            showsubtitle=T, subtitlelab = "Population structure from microsatellites", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)

grid.arrange(m2$plot[[1]],nrow=1)

#SINGLE FIGURE:  all runs of k=2 merged and aligned, 11-20
k_2_msat=readQ(files=choose.files(multi=TRUE),filetype="structure")
k2_msat_list <- alignK(k_2_msat)
k2_msat_merge <- mergeQ(k2_msat_list)

k2_msat <- plotQ(k2_msat_merge,returnplot=T,exportplot=F,quiet=T,basesize=11,
                grplab=onelabset2,subsetgrp=c("Massachusetts","Florida",  "Virgin Islands", "Azores"),grplabsize=6,
                linesize=0.8,pointsize=4, splab="K = 2", splabsize=12, 
                clustercol=c("lightgrey","lightskyblue"), showtitle=F, titlelab = "Roseate Tern Population Structure", titlesize=20,
                showsubtitle=T, subtitlelab = "Structure from microsatellites", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)
grid.arrange(k2_msat$plot[[1]],nrow=1)


########################################################################################################

#Arrange SNP and msats into one figure 
k2_msat <- plotQ(k2_msat_merge,returnplot=T,exportplot=F,quiet=T,basesize=11,
                 grplab=onelabset2,subsetgrp=c("Massachusetts", "Azores","Florida", "Virgin Islands"),grplabsize=6,
                 linesize=0.8,pointsize=4, splab="K = 2", splabsize=12, 
                 clustercol=c("lightgrey","lightskyblue"), 
                 showsubtitle=T, subtitlelab = "a.", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)

k2_snp <- plotQ(k2_snp_merge,returnplot=T,exportplot=F,quiet=T,basesize=11,
                grplab=onelabset1,subsetgrp=c("Massachusetts", "Azores","Florida", "Virgin Islands"),grplabsize=6,
                 linesize=0.8,pointsize=4, splab="K = 2", splabsize=12, 
                 clustercol=c("lightskyblue", "lightgrey"), 
                 showsubtitle=T, subtitlelab = "b.", subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1)


grid.arrange(k2_msat$plot[[1]],k2_snp$plot[[1]], nrow=2)
