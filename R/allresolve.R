

##################################
#5 This section takes the output from the fourth section (find the most abundant child at the site), and 
#any assigns remaining parents to most abundant child in study area
###############################

allresolve=function(lookupready,invertssite_obs,invertsnew,method){
  
  #method="numerous" or "widespread with "numerous" for chao1 & "widespread for chao2  
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  library(plyr)
  #lookup=lookupready
  #inverts=datatable
  #invertsnew=summarized_data
  #data=invertssite_obs
  #brings in files from part 4
  #method="widespread"
  #method2="divide"
  #inverts=read.table("invertssite_obs.csv",header=TRUE,sep=",") 
  
  #lookup=read.table("lookup.site.csv",header=TRUE,sep=",")
  #
  #invertsnew=read.table("siteinvertsnew.csv",header=TRUE,sep=",") 
  
  inverts=invertssite_obs
  lookup=lookupready
  
  inverts.sum=summaryBy(InvCount~obs+InvEventFK+InvName,FUN=sum,data=inverts)
  
  inverts=inverts.sum[,c(2,3,4,1)]
  colnames(inverts)=c("InvEventFK","InvName","InvCount","obs")
  
  #inverts$obs=1:nrow(inverts)
  
  inverts1=inverts
  invertsadd=inverts
  #invertsnew=invertsnew
  
  listcategory2=c("PHYLUM","CLASS","SUPERORDER","ORDER","SUBORDER","FAMILY",
                  "SUBFAMILY","TRIBE","GENUS_GROUP","GENUS","SPECIES_GROUP","SPECIES")
  
  run=1
  
  for (k in listcategory2){
    #extract info from parent column k and child column j
    inverts$InvName=as.character(inverts$InvName)  
    inverts.lookup=merge(inverts,lookup, by.x="InvName",by.y="TAXON_NAME",all.x=TRUE, all.y=FALSE)
    inverts.lookup.b=merge(invertsnew,lookup, by.x="InvName",by.y="TAXON_NAME",all.x=TRUE, all.y=FALSE)
    inverts.lookup1=inverts.lookup[inverts.lookup$IDLevel==k,]
    colnumber=which(colnames(inverts.lookup)==k) 
    numlist=na.omit(inverts.lookup1$no)
    num=mean(numlist)
    num2=num+1
    colnumberplus=colnumber+1
    colplusname=colnames( inverts.lookup1 )[colnumberplus] 
    j=colplusname  
    
    inverts.lookup2=inverts.lookup.b[inverts.lookup.b$no>num,] ## find all taxa that are higher in the hierarchy than level of interest from new dataset
    invertsmerge.test.a=subsetBy(~InvName + InvEventFK, no==min(no),data=inverts.lookup2) #coarsest ID if widespreads in lookup table
    #invertsmerge.test=ddply(invertsmerge.test.a,.(InvName,InvEventFK),randomRows,1)# this should not be needed unless a spelling or other error in lookup table
    
    inverts.lookup2=invertsmerge.test.a
    inverts.lookup2$k=inverts.lookup2[,k]
    inverts.lookup2$j=inverts.lookup2[,j]
    inverts.lookup2=data.frame(inverts.lookup2)
    inverts.lookup2$j=as.character(inverts.lookup2$j)
    inverts.lookup2=cbind(as.character(inverts.lookup2$k),as.character(inverts.lookup2$InvName),as.character(inverts.lookup2$InvEventFK),as.numeric(inverts.lookup2$InvCount),as.numeric(inverts.lookup2$no),inverts.lookup2$j)
    colnames(inverts.lookup2)=c("k","InvName.y","InvEventFK","InvCount","no","j")
    inverts.lookup2=data.frame(inverts.lookup2)
    inverts.lookup2$InvCount=as.numeric(as.character(inverts.lookup2$InvCount))
    inverts.lookup2$j=as.character(inverts.lookup2$j)
    l2=length(inverts.lookup2[,1])
    
    invertsmerge=merge(inverts,inverts.lookup2,by.x=c("InvName"),by.y=c("k"),all.x=TRUE,all.y=FALSE)
    
    invertsmerge$InvName=as.character(invertsmerge$InvName)
    invertsmerge$j=as.character(invertsmerge$j)
    invertsmerge$j[is.na(invertsmerge$j)] <- 0
    invertsmerge$InvCount.y=ifelse(invertsmerge$InvName==invertsmerge$j,0,invertsmerge$InvCount.y)
    invertsmerge$j=ifelse(invertsmerge$InvName==invertsmerge$j,0,invertsmerge$j)
    invertsmerge$InvName.y=as.character(invertsmerge$InvName.y)
    invertsmerge$InvCount.y=as.numeric(as.character(invertsmerge$InvCount.y))
    invertsmerge$j=as.character(invertsmerge$j)
    invertsmerge$InvName.y[is.na(invertsmerge$InvName.y)] <- "none"
    invertsmerge$InvCount.y[is.na(invertsmerge$InvCount.y)] <- 0
    
    if (method == "numerous"){
      invertsmerge.a=summaryBy(InvCount.y~InvEventFK.x+InvName+j+obs,invertsmerge,FUN=sum)
    }
    
    if (method == "widespread"){
      invertsmerge$InvCount.y=ifelse(invertsmerge$InvCount.y > 0,1,invertsmerge$InvCount.y)
      invertsmerge.a=summaryBy(InvCount.y~InvEventFK.x+InvName+j+obs,invertsmerge,FUN=sum)
    }
    
    invertsmerge=invertsmerge.a
    invertsmerge$InvCount.y=invertsmerge$InvCount.y.sum
    invertsmerge$InvName.y=invertsmerge$j
    ###################################################
    invertsmerge.subset=subsetBy(~obs+InvName, subset=InvCount.y==max(InvCount.y),data=invertsmerge) #find child with max count
    
    
    invertsmerge.subset2=ddply(invertsmerge.subset,.(obs,InvName),randomRows,1) #randomly select child if a tie
    
    invertsmerge.subset2$InvName.y=as.character(invertsmerge.subset2$InvName.y)
    invertsmerge.subset2$InvName=as.character(invertsmerge.subset2$InvName)
    invertsmerge.subset2$InvName.y=ifelse(invertsmerge.subset2$InvName.y!="0",invertsmerge.subset2$InvName.y, invertsmerge.subset2$InvName)
    invertsmerge.subset2$InvCount.y=as.numeric(as.character(invertsmerge.subset2$InvCount.y))
    same=subset(invertsmerge.subset2,invertsmerge.subset2$InvName==invertsmerge.subset2$InvName.y)
    notsame=subset(invertsmerge.subset2,invertsmerge.subset2$InvName!=invertsmerge.subset2$InvName.y)
    
    invertsmerge.subset2=rbind(same,notsame)
    invertsold=inverts[c(2,4,3)]
    invertsold=data.frame(invertsold)
    
    invertsmerge.subset2=merge(invertsmerge.subset2,invertsold,by.x=c("obs","InvName"),by.y=c("obs","InvName"),all.x=TRUE,all.y=TRUE)
    invertsmerge.subset2$obs=as.character(as.numeric(invertsmerge.subset2$obs))
    invertsmerge.subset2$InvName.y[is.na(invertsmerge.subset2$InvName.y)] <- 0
    inverts=cbind(as.character(invertsmerge.subset2$InvEventFK),as.character(invertsmerge.subset2$InvName.y),as.numeric(as.character(invertsmerge.subset2$InvCount)),as.numeric(as.character(invertsmerge.subset2$obs)))
    colnames(inverts)=c("InvEventFK","InvName","InvCount","obs")
    inverts=data.frame(inverts)
    inverts$InvCount=as.numeric(as.character(inverts$InvCount))
    inverts$obs=as.numeric(as.character(inverts$obs))
    inverts$InvName=as.character(inverts$InvName)
    
    inverts.a=cbind(as.character(invertsmerge.subset2$InvEventFK),as.character(invertsmerge.subset2$InvName),as.character(invertsmerge.subset2$InvName.y),as.numeric(as.character(invertsmerge.subset2$InvCount.y.sum)),as.numeric(as.character(invertsmerge.subset2$obs)))
    colnames(inverts.a)=c("InvEventFK","Orig","InvName","InvCount","obs")
    inverts.a=data.frame(inverts.a)
    inverts.a$InvCount=as.numeric(as.character(inverts.a$InvCount))
    #problem with not knowing the summ of before iteration
    # merge previous iteration with new interation
    lname=paste0("InvName",paste0(k))
    w=which(colnames(inverts1)=="InvName")
    
    colnames(inverts1)[w]=lname
    
    lname=paste0("InvName",paste0(k))
    name=colnames(inverts1[lname])
    
    test=merge(inverts1,inverts.a,by.x=c("obs",lname),by.y=c("obs","Orig"),all.x=TRUE,all.y=TRUE)
    #test$InvNameTRIBE=as.character(test$InvNameTRIBE)
    test$InvName=as.character(test$InvName)
    
    inverts1=test
    inverts1$InvCount=as.numeric(inverts1$InvCount.y)
    invertsnew=inverts[c(1,2,3)]
    
    #assign(paste0(lname),inverts.a)
    
    invertsnew=summaryBy(InvCount~InvEventFK+InvName,data=invertsnew,FUN=sum)
    colnames(invertsnew)=c("InvEventFK","InvName","InvCount")
    invertsnew$InvCount=as.numeric(as.character(invertsnew$InvCount))
    invertsnew$obs=0
    run=run+1
    
    lnameold=lname
    
  }
  
  
  
  ##calculate summaryfiles
  
  inverts1.a=inverts1[,c(1,49,50,15,13,17,12,20,11,23,10,26,9,29,8,32,7,35,6,38,5,41,4,44,3,47,2,51)]
  #write.table(inverts1.a,"assign_level.csv",sep=",",row.names=FALSE)
  
  dir.create("./all")
  write.table(inverts1.a,"./all/all_iter.csv",sep=",",row.names=FALSE) 
  
  
  alllevel=inverts1.a[,c(1,5,3)]
  
  colnames(alllevel)=c("obs","sitetaxa","newtaxa")
  alllevel=unique(alllevel)
  sitelevel=read.table("./site/sitelevel.csv",header=TRUE,sep=",")
  
  sitelevel$obs=as.character(sitelevel$obs)
  alllevel$obs=as.character(alllevel$obs)
  alllevel=data.frame(alllevel)
  sitelevel=data.frame(sitelevel)
  final=merge(alllevel,sitelevel,by.x=c("obs","sitetaxa"),by.y=c("obs","sitetaxa"),all.x=TRUE,all.y=TRUE)
  final.sort=final[order(final$InvName),]
  final.sort=na.omit(final.sort)
  
  write.table(final.sort,"./all/ambiguous_APTC_results.csv",sep=",",row.names=FALSE)
  
  
  
}