

#############################Part 1
#Instructions: set working directory and put 2 input files in there. "Invertsinput.csv","lookup.csv". Example files are in the "data" folder on github. Column headings must match.
#any taxonomic level that does not have a name associated with it must be "NA"
# all of the taxa in the invertsinput.csv must be in the lookup table
#set method to "widespread" or "numerous; set method_divide to "nodivide or "divide";

##The final file is ambiguous_APTC_results.csv in the "all" folder, which has the original taxa (InvName), the new taxa when
#ambiguous taxa are assigned at the site level (sitetaxa) and the new taxa when ambiguous taxa are assigned at the dataset level (newtaxa)
# to produce final datasets for other common methods in Meredith et al., run extracode, part 2


setwd("C:/christy/examples") # set working directory where input files are located; columns must be the same as in example)
method="numerous"  #options are "numerous" or "widespread" (see Meredith et al. 2019)
method_divide="divide" #options are "divide" or "nodivide" (see Meredith et al. 2019)

library(doBy)
library(plyr)
library(devtools)

install_github("christystarr/ambiguoustaxa")
library(ambiguoustaxa)

inverts=read.table("invertsinput.csv",sep=",",header=TRUE)
lookup=read.table("lookup.csv",sep=",",header=TRUE)

setup(lookup)
missinglevels(lookupsetup)
siteresolve(inverts,lookupready,"nodivide") #methods are "divide" or "nodivide" referring to whether the parent taxa will be divided among children at the site, or assigned to the most numerous at the site
allresolve(lookupready,invertssite_obs,invertsnew,"widespread") # methods are "numerous" or "widespread" refering to whether ambiguous taxa will be assigned to the child taxa that is most abundance,#or the child taxa found at the most sites



############################## Part 2
#after running part 1, simply select the entire code below; final files after ambiguous taxa are resolved will be put into the working directory
#creates datasets in which ambiguous taxa are resolve in different ways: APTC_S,APTC_SG,APTC_SG1,APTC_SG2,RPCK_S; puts these in working directory
#also estimates variables describing number of singletons and doubleton and number uniques and duplicates; these will be in the global environment
#################################

setwd("all")
methods.b=read.table("ambiguous_APTC_results.csv",sep=",",header=TRUE)


if(method_divide=="divide") {
  methods.b=methods.b[,c(1,3,4,5,7,2)]
}

if (method_divide != "divide"){
  methods.b=methods.b[,c(1,3,4,5,7,2)]
}

##calculate number of singletons and doubletons in original dataset
setwd("C://christy/examples")
colnames(methods.b)=c("obs","newtaxa","InvEventFK","InvName","InvCount","sitetaxa")
invertsorig=read.table("invertsinput.csv",sep=",",header=TRUE)
origsum=summaryBy(InvCount~InvName,FUN=c(sum),data=invertsorig)
origsing=subset(origsum,origsum$InvCount.sum==1)
origdoub=subset(origsum,origsum$InvCount.sum==2)
singledoub=rbind(origsing,origdoub) 
length(origsing[,1])
length(origdoub[,1])

## calculate number of uniques and duplicates in original dataset
origcount=summaryBy(InvCount~InvName+InvEventFK,FUN=c(sum),data=invertsorig)
origcount$InvCount.sum=ifelse(origcount$InvCount.sum>0,1,0)
origcount2=summaryBy(InvCount.sum~InvName,FUN=sum,data=origcount)
uniqueorig=subset(origcount2,origcount2$InvCount.sum.sum==1)
duporig=subset(origcount2,origcount2$InvCount.sum.sum==2)
uniqdup=rbind(uniqueorig,duporig)
length(uniqueorig[,1])
length(duporig[,1])


# prepare to summarize
amb_sing2=merge(methods.b,singledoub,by.x=c("newtaxa"),by.y=c("InvName"),all.x=TRUE,all.y=FALSE)
amb_sing3=merge(amb_sing2,uniqdup,by.x=c("newtaxa"),by.y=c("InvName"),all.x=TRUE,all.y=FALSE)
amb_sing3=amb_sing3[,c(2,3,5,4,6,1,7,8)]
colnames(amb_sing3)=c("obs","InvEventFK","InvCount","InvName","sitetaxa","newtaxa","sing_doub","uniq_dup")
amb_sing3$sing_doub[is.na(amb_sing3$sing_doub)] <- 0
amb_sing3$uniq_dup[is.na(amb_sing3$uniq_dup)] <- 0

if(method=="numerous"){

## create APTC_SG dataset
APTC_SG=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3,FUN=sum)
amb_sing3$newtaxa=as.character(amb_sing3$newtaxa)
amb_sing3$InvName=as.character(amb_sing3$InvName)
SG_singdoub=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=APTC_SG)
singSG=subset(SG_singdoub,SG_singdoub$InvCount.sum.sum==1)
lsingSG=length(singSG[,1])
doubSG=subset(SG_singdoub,SG_singdoub$InvCount.sum.sum>1 & SG_singdoub$InvCount.sum.sum<=2 )
ldoubSG=length(doubSG[,1])
write.table(APTC_SG,"APTC_SG.csv",sep=",",row.names=FALSE)

#create APTC_SG1 dataset
amb_sing3.SG1=amb_sing3[!(amb_sing3$sing_doub!=0 & amb_sing3$InvName!=amb_sing3$newtaxa),]
APTC_SG1=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3.SG1,FUN=sum)
test=unique(APTC_SG1$newtaxa)

APTC_SG1singdoub=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=APTC_SG1)

singSG1=subset(APTC_SG1singdoub,APTC_SG1singdoub$InvCount.sum.sum==1)
lsingSG1=length(singSG1[,1])
doubSG1=subset(APTC_SG1singdoub,APTC_SG1singdoub$InvCount.sum.sum>1 & APTC_SG1singdoub$InvCount.sum.sum<=2 )
ldoubSG1=length(doubSG1[,1])
write.table(APTC_SG1,"APTC_SG1.csv",sep=",",row.names=FALSE)

# create APTC_SG2 dataset
amb_sing3.SG2=amb_sing3[!(amb_sing3$sing_doub!=0 & amb_sing3$newtaxa==amb_sing3$sitetaxa & amb_sing3$InvName!=amb_sing3$newtaxa),]
APTC_SG2=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3.SG2,FUN=sum)
test=unique(APTC_SG2$newtaxa)

APTC_SG2singdoub=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=APTC_SG2)

singSG2=subset(APTC_SG2singdoub,APTC_SG2singdoub$InvCount.sum.sum==1)
lsingSG2=length(singSG2[,1])
doubSG2=subset(APTC_SG2singdoub,APTC_SG1singdoub$InvCount.sum.sum>1 & APTC_SG2singdoub$InvCount.sum.sum<=2 )
ldoubSG2=length(doubSG2[,1])
write.table(APTC_SG2,"APTC_SG2.csv",sep=",",row.names=FALSE)
## create APTC_S dataset

APTC_S=summaryBy(InvCount~InvEventFK + sitetaxa,FUN=sum,data=amb_sing3)
APTC_singdoub=summaryBy(InvCount~sitetaxa,FUN=sum,data=amb_sing3)
singAPTC_s=subset(APTC_singdoub,APTC_singdoub$InvCount.sum==1)
doubAPTC_s=subset(APTC_singdoub,APTC_singdoub$InvCount.sum==2)
l_singAPTC_s=length(singAPTC_s[,1])
l_doubAPTC_s=length(doubAPTC_s[,1])
write.table(APTC_S,"APTC_S.csv",sep=",",row.names=FALSE)

##create RPKC_s datast
amb_sing3$diff=ifelse(amb_sing3$sitetaxa!=amb_sing3$InvName,1,0)
RPKC_s.a=subset(amb_sing3,amb_sing3$diff==0)
RPKC_s=summaryBy(InvCount~InvEventFK + sitetaxa,FUN=sum,data=RPKC_s.a)
RPKC_s_singdoub=summaryBy(InvCount.sum~sitetaxa,FUN=sum,data=RPKC_s)
sing_RPKC_s=subset(RPKC_s_singdoub,RPKC_s_singdoub$InvCount.sum==1)
doub_RPKC_s=subset(RPKC_s_singdoub,RPKC_s_singdoub$InvCount.sum==2)
lsing_RPKCS=length(sing_RPKC_s[,1])
ldoub_RPKCS=length(doub_RPKC_s[,1])
write.table(RPKC_s,"RPKC_s.csv",sep=",",row.names=FALSE)

}

###if method = widespread instead of numerous
# create APTC_SG  dataset
if (method != "numerous") {
SGcount=summaryBy(InvCount~newtaxa+InvEventFK,FUN=c(sum),data=amb_sing3)
APTC_SG=SGcount
SGcount$InvCount.sum=ifelse(SGcount$InvCount.sum>0,1,0)
SGcount2=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=SGcount)
uniqueSG=subset(SGcount2,SGcount2$InvCount.sum.sum==1)
lSG.uniques=length(uniqueSG[,1])
dupSG=subset(SGcount2,SGcount2$InvCount.sum.sum>1 & SGcount2$InvCount.sum.sum<=2)
lSG.dups=length(dupSG[,1])

write.table(APTC_SG,"APTC_SG.csv",sep=",")
#write.table(dupSG,"origdup.csv",sep=",")


# create APTC_SG1 dataste
amb_sing3$InvName=as.character(amb_sing3$InvName)
amb_sing3$newtaxa=as.character(amb_sing3$newtaxa)
amb_sing3.SG1.uniq=amb_sing3[!(amb_sing3$uniq_dup!=0 & amb_sing3$InvName!=amb_sing3$newtaxa),]
APTC_SG1=amb_sing3.SG1.uniq
APTC_SG1.uniq=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3.SG1.uniq,FUN=sum)
#test=unique(APTC_SG1$newtaxa)
APTC_SG1.uniq$InvCount.sum=ifelse(APTC_SG1.uniq$InvCount.sum>0,1,0)
APTC_SG1.uniqdup=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=APTC_SG1.uniq)
SG1.uniques=subset(APTC_SG1.uniqdup,APTC_SG1.uniqdup$InvCount.sum.sum==1)
SG1.dups=subset(APTC_SG1.uniqdup,APTC_SG1.uniqdup$InvCount.sum.sum>1 & APTC_SG1.uniqdup$InvCount.sum.sum<=2)
lSG1.uniques=length(SG1.uniques[,1])
lSG1.dups=length(SG1.dups[,1])
write.table(APTC_SG1,"APTC_SG1.csv",sep=",")
##hmmmm

#create APTC_SG2 dataset
amb_sing3.SG2=amb_sing3[!(amb_sing3$uniq_dup!=0 & amb_sing3$newtaxa!=amb_sing3$sitetaxa),]
APTC_SG2=amb_sing3.SG2
APTC_SG2=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3.SG2,FUN=sum)
APTC_SG2.uniq=summaryBy(InvCount~InvEventFK+newtaxa,data=amb_sing3.SG2,FUN=sum)
APTC_SG2.uniq$InvCount.sum=ifelse(APTC_SG2.uniq$InvCount.sum>0,1,0)
APTC_SG2.uniqdup=summaryBy(InvCount.sum~newtaxa,FUN=sum,data=APTC_SG2.uniq)
SG2.uniques=subset(APTC_SG2.uniqdup,APTC_SG2.uniqdup$InvCount.sum.sum==1)
SG2.dups=subset(APTC_SG2.uniqdup,APTC_SG2.uniqdup$InvCount.sum.sum>1 & APTC_SG2.uniqdup$InvCount.sum.sum<=2)
lSG2.uniques=length(SG2.uniques[,1])
lSG2.dups=length(SG2.dups[,1])

write.table(APTC_SG2,"APTC_SG2.csv",sep=",")

}
######################################






