

#'A missing levels function
#'  this  part of the program finds the most common child at the InvEventFK;
#If there is a tie, a random child is selected
#' @param 
#' @keywords 
#' @export
#' @examples
#' siteresolve()

############################################
#4 this  part of the program finds the most common child at the InvEventFK;
#If there is a tie, a random child is selected
##########################################
siteresolve=function(datatable,lookupready,method_divide){
  setwd("./output")
  library(plyr)
  
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  
  inverts=datatable
  lookup=lookupready
  
  inverts$obs=1:nrow(inverts)
  inverts1=inverts
  invertsadd=inverts
  invertsnew=inverts
  lnameold="InvPHYLUM"
  lcountold="CntPHYLUM"
  
  #iterate through each level
  listcategory2=c("PHYLUM","CLASS","SUPERORDER","ORDER","SUBORDER","FAMILY","SUBFAMILY","TRIBE","GENUS_GROUP","GENUS","SPECIES_GROUP","SPECIES") 
  #listcategory2=c("PHYLUM","CLASS") 
  #listcategory2=c("PHYLUM","CLASS","SUPERORDER","ORDER","SUBORDER","FAMILY","SUBFAMILY","TRIBE","GENUS_GROUP","GENUS") 
  
  #k2="SUBORDER"
  k="PHYLUM"
  run=1
  
  #lnameold=paste0("Inv",paste0(k))
  #lcountold=paste0("Cnt",paste0(k))
  #w=which(colnames(inverts1)=="InvName")
  #w2=which(colnames(inverts1)=="InvCount")
  #colnames(inverts1)[w]=lnameold
  #colnames(inverts1)[w2]=lcountold
  
  
  
  for (k in listcategory2){
    
    #invertstest1=inverts  
    
    inverts$InvName=as.character(inverts$InvName)  
    inverts.lookup=merge(inverts,lookup, by.x="InvName",by.y="TAXON_NAME",all.x=TRUE, all.y=FALSE) # original
    inverts.lookup.b=merge(invertsnew,lookup, by.x="InvName",by.y="TAXON_NAME",all.x=TRUE, all.y=FALSE) #new after iterations
    inverts.lookup1=inverts.lookup[inverts.lookup$IDLevel==k,]#find all at level of id
    
    colnumber=which(colnames(inverts.lookup)==k) 
    numlist=na.omit(inverts.lookup1$no)
    num=mean(numlist)
    num2=num+1
    colnumberplus=colnumber+1
    colplusname=colnames( inverts.lookup1 )[colnumberplus] 
    j=colplusname  
    
    inverts.lookup2=inverts.lookup.b[inverts.lookup.b$no>num,] ## find all taxa that are higher in the hierarchy than level of interest from new dataset
    invertsmerge.test.a=subsetBy(~InvName + InvEventFK, no==min(no),data=inverts.lookup2) #coarsest ID if widespreads in lookup table
    invertsmerge.test=ddply(invertsmerge.test.a,.(InvEventFK,InvName),randomRows,1)# this should not be needed unless a spelling or other error in lookup table 
    #creating two widespread entries
    
    #mini-lookp to test for errors in lookup table
    invertsmergetest.na <- inverts.lookup[is.na(inverts.lookup$PHYLUM),]
    write.table(invertsmergetest.na,"possibleerrors.csv",sep=",")   #possible errors, no match in lookup table
    
    #############
    
    inverts.lookup2=invertsmerge.test
    
    inverts.lookup2$k=inverts.lookup2[,k] #current level
    inverts.lookup2$j=inverts.lookup2[,j] #next level  of ID
    inverts.lookup2=data.frame(inverts.lookup2)
    inverts.lookup2$j=as.character(inverts.lookup2$j)
    inverts.lookup2=cbind(as.character(inverts.lookup2$k),as.character(inverts.lookup2$InvName),as.character(inverts.lookup2$InvEventFK),as.numeric(inverts.lookup2$InvCount),as.numeric(inverts.lookup2$no),inverts.lookup2$j)
    colnames(inverts.lookup2)=c("k","InvName.y","InvEventFK","InvCount","no","j")
    inverts.lookup2=data.frame(inverts.lookup2)
    inverts.lookup2$InvCount=as.numeric(as.character(inverts.lookup2$InvCount))
    inverts.lookup2$j=as.character(inverts.lookup2$j)
    inverts.lookup2$InvName=as.character(inverts.lookup2$InvName)
    inverts.lookup2$InvEventFK=as.character(inverts.lookup2$InvEventFK)
    inverts.lookup2$no=as.numeric(as.character(inverts.lookup2$no))# added this and next two lines
    inverts.lookup2$k=as.character(inverts.lookup2$k)
    l2=length(inverts.lookup2[,1])
    #merge parents with all possible children; children have the same ID category as the parents InvName
    invertsmerge=merge(inverts,inverts.lookup2,by.x=c("InvEventFK","InvName"),by.y=c("InvEventFK","k"),all.x=TRUE,all.y=FALSE)
    
    invertsmerge$InvName=as.character(invertsmerge$InvName)
    invertsmerge$j=as.character(invertsmerge$j)
    invertsmerge$j[is.na(invertsmerge$j)] <- 0 # if no possible child for that level then give a 0
    invertsmerge$InvName.y=as.character(invertsmerge$InvName.y)
    invertsmerge$InvCount.y=as.numeric(as.character(invertsmerge$InvCount.y))
    invertsmerge$j=as.character(invertsmerge$j)
    invertsmerge$InvName.y[is.na(invertsmerge$InvName.y)] <- "none"
    invertsmerge$InvCount.y[is.na(invertsmerge$InvCount.y)] <- 0
    
    invertsmerge.a=summaryBy(InvCount.y~InvEventFK+InvName+j+obs,invertsmerge,FUN=sum) # sum all at that level of ID (could be multiple in same genus,etc)
    totalby=summaryBy(InvCount.y.sum~obs + InvName, FUN=sum,data=invertsmerge.a)
    #if divided
    
    invertsmerge=invertsmerge.a
    invertsmerge$InvCount.y=invertsmerge$InvCount.y.sum
    invertsmerge$InvName.y=invertsmerge$j
    
    #if x="nodivide"{
    invertsmerge.subset=subsetBy(~obs, subset=InvCount.y==max(InvCount.y),data=invertsmerge) #find the child at the site with the most counts
    invertsmerge.subset2=ddply(invertsmerge.subset,.(obs),randomRows,1) #if a tie then randomly choose one child
    
    invertsmerge.subset2$InvName.y=as.character(invertsmerge.subset2$InvName.y)
    invertsmerge.subset2$InvName=as.character(invertsmerge.subset2$InvName)
    prelim_subset=invertsmerge.subset2
    prelim_subset=prelim_subset[,c(4,3)]
    #}
    invertsmerge.subset2$InvName.y=ifelse(invertsmerge.subset2$InvName.y!="0",invertsmerge.subset2$InvName.y, invertsmerge.subset2$InvName)
    invertsmerge.subset2$InvCount.y=as.numeric(as.character(invertsmerge.subset2$InvCount.y))
    same=subset(invertsmerge.subset2,invertsmerge.subset2$InvName==invertsmerge.subset2$InvName.y)
    notsame=subset(invertsmerge.subset2,invertsmerge.subset2$InvName!=invertsmerge.subset2$InvName.y)
    invertsmerge.subset2=rbind(same,notsame)
    
    
    ####################This portion compares child to original file 
    
    invertsold=inverts[c(4,3)]
    invertsold=data.frame(invertsold)
    invertsmerge.subset2=merge(invertsmerge.subset2,invertsold,by.x="obs",by.y="obs",all.x=TRUE,all.y=TRUE)
    invertsmerge.subset2$obs=as.character(as.numeric(invertsmerge.subset2$obs))
    invertsmerge.subset2$InvName.y[is.na(invertsmerge.subset2$InvName.y)] <- 0
    
    
    #####################
    
    if (method_divide == "divide") {
      
      
      lname=paste0("Iterite",paste0(k))
      print(lname)
      
      invertsold=inverts[c(2,4,3)]
      invertsold=data.frame(invertsold)
      invertsmerge.subset2a=merge(invertsmerge.a,invertsold,by.x=c("obs","InvName"),by.y=c("obs","InvName"),all.x=TRUE,all.y=TRUE)
      invertsmerge.subset2aa=merge(invertsmerge.subset2a,totalby,by.x=c("obs","InvName"),by.y=c("obs","InvName"),all.x=TRUE,all.y=TRUE)
      invertsmerge.subset2aa$j=as.character(invertsmerge.subset2aa$j)
      invertsmerge.sub3=merge(invertsmerge.subset2aa,prelim_subset,by.x="obs",by.y="obs",all.x="TRUE",all.y="TRUE")
      
      invertsmerge.sub3$InvCount.y.sum.sum=ifelse(invertsmerge.sub3$InvCount.y.sum.sum==0,1,invertsmerge.sub3$InvCount.y.sum.sum)
      options(scipen=999)
      invertsmerge.sub3$percent=invertsmerge.sub3$InvCount.y.sum/invertsmerge.sub3$InvCount.y.sum.sum
      invertsmerge.sub3$new=invertsmerge.sub3$percent * invertsmerge.sub3$InvCount
      
      invertsmerge.sub3$j.y=as.character(invertsmerge.sub3$j.y)
      invertsmerge.sub3$j.x=as.character((invertsmerge.sub3$j.x))
      invertsmerge.subset2check=invertsmerge.sub3[c(1,3,2,4,5,6,8,10)]
      colnames(invertsmerge.subset2check)=c("obs","InvEventFK","InvName","j","InvCount.y.sum","InvCount.y","InvName.y","InvCount")
      
      invertsmerge.subset2=invertsmerge.subset2check
      invertsmerge.subset2$InvName.y=ifelse(invertsmerge.subset2$j=="0",invertsmerge.subset2$InvName,invertsmerge.subset2$j)
      invertsmerge.subset2$obs=as.character(as.numeric(invertsmerge.subset2$obs))
      invertsmerge.subset2$InvName.y[is.na(invertsmerge.subset2$InvName.y)] <- 0
      
      invertsmerge.subset2$InvCount=ifelse(invertsmerge.subset2$InvName==invertsmerge.subset2$InvName.y,invertsmerge.subset2$InvCount.y,invertsmerge.subset2$InvCount)
    }
    
    
    inverts=cbind(as.character(invertsmerge.subset2$InvEventFK),as.character(invertsmerge.subset2$InvName.y),as.numeric(as.character(invertsmerge.subset2$InvCount)),as.numeric(as.character(invertsmerge.subset2$obs)))
    inverts=invertsmerge.subset2[,c(2,7,8,1)]
    colnames(inverts)=c("InvEventFK","InvName","InvCount","obs")
    inverts=data.frame(inverts)
    inverts$InvCount=as.numeric(as.character(inverts$InvCount))
    inverts$obs=as.numeric(as.character(inverts$obs))
    inverts$InvName=as.character(inverts$InvName)
    inverts$InvEventFK=as.character(inverts$InvEventFK)
    inverts=unique(inverts)
    
    lname=paste0("Inv",paste0(k))
    lcount=paste0("Cnt",paste0(k))
    inverts.a=cbind(as.character(invertsmerge.subset2$InvEventFK),as.character(invertsmerge.subset2$InvName),as.character(invertsmerge.subset2$InvName.y),as.numeric(as.character(invertsmerge.subset2$InvCount)),as.numeric(as.character(invertsmerge.subset2$InvCount.y)),as.numeric(as.character(invertsmerge.subset2$obs)))
    colnames(inverts.a)=c("InvEventFK",paste(lnameold),paste(lname),paste(lcount),paste(lcountold),"obs")
    inverts.a=data.frame(inverts.a)
    inverts.a$obs=as.character(inverts.a$obs)
    #inverts.a$InvCount=as.numeric(as.character(inverts.a$InvCount))
    invertsnew=inverts[c(1,2,3)]
    
    write.table(inverts.a,(paste0(lname,".csv")),sep=",")
    
    
    lnameold=lname
    lcountold=lcount
    
    invertsnew=summaryBy(InvCount~InvEventFK+InvName,data=invertsnew,FUN=sum)
    colnames(invertsnew)=c("InvEventFK","InvName","InvCount")
    invertsnew$InvCount=as.numeric(as.character(invertsnew$InvCount))
    invertsnew$obs=0
    
    #site-level outputs
    
    
    if (method_divide != "divide") {
      #print("yes")
      inverts1=merge(inverts,inverts1,by.x=c("InvEventFK","obs"),by.y=c("InvEventFK","obs"),all.x=TRUE,all.y=TRUE)
      
      
      #write.table(inverts1,"site_iter.csv",sep=",",row.names=FALSE) 
      #colnames(inverts1)=c("InvEventFK","obs","sitetaxa","CNT13","SPECIES","CNT12","SPECIES_GROUP","CNT11","GENUS","CNT10","GENUSGROUP","CNT9","TRIBE","CNT7","SUBFAMILY","CNT7","FAMILY","CNT6","SUBORDER","CNT5","ORDER","CNT4","SUPERORDER","CNT3","CLASS","CNT2","PHYLUM","CNT1")
    }
  }
  
  #t=summaryBy(InvCount~InvEventFK+InvName+obs,data=inverts)
  if (method_divide == "divide"){
    #print("yes")
    PHYLUM=read.table("InvPHYLUM.csv",sep=",")
    
    CLASS=read.table("InvCLASS.csv",sep=",")
    
    SUPERORDER=read.table("InvSUPERORDER.csv",sep=",")
    
    ORDER=read.table("InvORDER.csv",sep=",")
    
    SUBORDER=read.table("InvSUBORDER.csv",sep=",")
    
    FAMILY=read.table("InvFAMILY.csv",sep=",")
    
    SUBFAMILY=read.table("InvSUBFAMILY.csv",sep=",")
    
    TRIBE=read.table("InvTRIBE.csv",sep=",")
    
    GENUS_GROUP=read.table("InvGENUS_GROUP.csv",sep=",")
    
    #GENUS_GROUP=read.table("InvGENUS_GROUP.csv",sep=",")
    
    GENUS=read.table("InvGENUS.csv",sep=",")
    
    SPECIES_GROUP=read.table("InvSPECIES_GROUP.csv",sep=",")
    
    SPECIES=read.table("InvSPECIES.csv",sep=",")
    
    merge1=merge(CLASS,SUPERORDER,by.x=c("obs","InvEventFK","InvCLASS","CntCLASS"),by.y=c("obs","InvEventFK","InvCLASS","CntCLASS"),all.x=TRUE,all.y=TRUE)
    
    merge2=merge(merge1,ORDER,by.x=c("obs","InvEventFK","InvSUPERORDER","CntSUPERORDER"),by.y=c("obs","InvEventFK","InvSUPERORDER","CntSUPERORDER"),all.x=TRUE,all.y=TRUE)
    
    merge3=merge(merge2,SUBORDER,by.x=c("obs","InvEventFK","InvORDER","CntORDER"),by.y=c("obs","InvEventFK","InvORDER","CntORDER"),all.x=TRUE,all.y=TRUE)
    
    merge4=merge(merge3,FAMILY,by.x=c("obs","InvEventFK","InvSUBORDER","CntSUBORDER"),by.y=c("obs","InvEventFK","InvSUBORDER","CntSUBORDER"),all.x=TRUE,all.y=TRUE)
    
    merge5=merge(merge4,SUBFAMILY,by.x=c("obs","InvEventFK","InvFAMILY","CntFAMILY"),by.y=c("obs","InvEventFK","InvFAMILY","CntFAMILY"),all.x=TRUE,all.y=TRUE)
    
    merge6=merge(merge5,TRIBE,by.x=c("obs","InvEventFK","InvSUBFAMILY","CntSUBFAMILY"),by.y=c("obs","InvEventFK","InvSUBFAMILY","CntSUBFAMILY"),all.x=TRUE,all.y=TRUE)
    
    merge7=merge(merge6,GENUS_GROUP,by.x=c("obs","InvEventFK","InvTRIBE","CntTRIBE"),by.y=c("obs","InvEventFK","InvTRIBE","CntTRIBE"),all.x=TRUE,all.y=TRUE)
    
    merge8=merge(merge7,GENUS,by.x=c("obs","InvEventFK","InvGENUS_GROUP","CntGENUS_GROUP"),by.y=c("obs","InvEventFK","InvGENUS_GROUP","CntGENUS_GROUP"),all.x=TRUE,all.y=TRUE)
    
    merge9=merge(merge8,SPECIES_GROUP,by.x=c("obs","InvEventFK","InvGENUS","CntGENUS"),by.y=c("obs","InvEventFK","InvGENUS","CntGENUS"),all.x=TRUE,all.y=TRUE)
    
    merge9.a=merge(merge9,SPECIES,by.x=c("obs","InvEventFK","InvSPECIES_GROUP","CntSPECIES_GROUP"),by.y=c("obs","InvEventFK","InvSPECIES_GROUP","CntSPECIES_GROUP"),all.x=TRUE,all.y=TRUE)
    
    merge10=merge9.a[,c(1,2,25,26,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]
    #merge4=merge(merge3,FAMILY,by.x=c("obs","InvEventFK","InvSUBORDER","CntCLASS"),by.y=c("obs","InvEventFK","InvSUBORDER","CntFAMILY"),all.x=TRUE,all.y=FALSE)
    
    #merge5=merge(merge4,SUBFAMILY,by.x=c("obs","InvEventFK","InvFAMILY","CntCLASS"),by.y=c("obs","InvEventFK","InvFAMILY","CntSUBFAMILY"),all.x=TRUE,all.y=FALSE)
    
    #inverts1=merge(inverts1,inverts.a,by.x=c("obs"),by.y="obs",all.x=TRUE,all.y=TRUE)
    #inverts1$InvCount=as.numeric(inverts1$InvCount.y)
    dir.create("./site")
    write.table(merge10,"./site/site_iter.csv",sep=",",row.names=FALSE) #most abundant child at each iteration
    #}
    
    #write.table(inverts,"invertssite_obs.csv",sep=",") #observation number with new InvName assigned if applicable
    
    #write.table(invertsnew,"invertsnew.csv",sep=",") #summaryBy of inverts, by InvEventFK and InvName
    #l1=length(inverts1[1,])
    #l1.a=l1-2
    #if (method_divide=="divide"){
    # print("yes") 
    sitelevel=merge9.a[,c(1,2,23,24,25,26)]
    #write.table(merge10,"sequence.csv",sep=",")
    colnames(sitelevel)=c("obs","InvEventFK","InvName","InvCount","sitetaxa","sitecount")
    write.table(sitelevel,"./site/sitelevel.csv",sep=",",row.names=FALSE) # final site output- most abundant child at site
    #}
  }
  
  if (method_divide !="divide"){
    colnames(inverts1)=c("InvEventFK","obs","sitetaxa","CNT13","SPECIES","CNT12","SPECIES_GROUP","CNT11","GENUS","CNT10","GENUSGROUP","CNT9","TRIBE","CNT7","SUBFAMILY","CNT7","FAMILY","CNT6","SUBORDER","CNT5","ORDER","CNT4","SUPERORDER","CNT3","CLASS","CNT2","PHYLUM","CNT1")
    
    #print("yes")
    sitelevel=inverts1[,c(2,1,25,4,3,4)]
    #sitelevel=cbind(inverts1$obs,inverts1$InvEventFK, inverts1$InvNamePHYLUM, inverts1$InvCount.y.11,inverts1$InvName)
    colnames(sitelevel)=c("obs","InvEventFK","InvName","InvCount","sitetaxa","sitecount")
    sitelevel$sitetaxa=as.character(sitelevel$sitetaxa)
    sitelevel$InvName=as.character(sitelevel$InvName)
    
    dir.create("./site")
    write.table(inverts1,"./site/site_iter.csv",sep=",",row.names=FALSE) 
    write.table(sitelevel,"./site/sitelevel.csv",sep=",",row.names=FALSE) # final site output- most abundant child at site
  }
  assign("invertssite_obs",inverts,.GlobalEnv)
  
  assign("invertsnew",invertsnew,.GlobalEnv)
  
  write.table(invertsnew,"./site/siteinvertsnew.csv",sep=",",row.names=FALSE)
  file.remove("InvSPECIES.csv","InvSPECIES_GROUP.csv","InvGENUS.csv","InvGENUS_GROUP.csv","InvTRIBE.csv","InvSUBFAMILY.csv","InvFAMILY.csv","InvSUBORDER.csv","InvCLASS.csv","InvPHYLUM.csv","InvORDER.csv","InvSUPERORDER.csv")
  
}
