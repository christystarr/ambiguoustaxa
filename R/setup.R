

#' A setup function
#'
#' This section of the program fills in NAs sandwiched between non-NA's with the name of the category that is the next coarsest category
#that is not an NA; this is because program will crash if NAs are present between two non NA's
#' @param 
#' @keywords 
#' @export
#' @examples
#' setup()

######################################################################################
setup<-function(lookuptable){
  
  #write function to randomly select a row if there is a tie for the child
  randomRows = function(df,n){
    return(df[sample(nrow(df),n),])
  }
  
  names(lookup)[names(lookup)=="SUBFAM"]<-"SUBFAMILY"
  
  #change species group to be GENUS followed by group to match #TAXON_NAME, unless no SPECIES_GROUP then leave "NA"
  lookup.notna=lookup[!is.na(lookup$SPECIES_GROUP),]
  lookup.na=lookup[is.na(lookup$SPECIES_GROUP),]
  lookup.notna$SPECIES_GROUP=paste(lookup.notna$GENUS,lookup.notna$SPECIES_GROUP,sep=" ")
  lookup=rbind(lookup.na,lookup.notna)
  
  #change SPECIES to GENUS + SPECIES to match TAXON_NAME
  lookuptry.unique=unique(lookup)
  lookuptry.notna=lookuptry.unique[!is.na(lookuptry.unique$SPECIES),]
  lookuptry.na=lookuptry.unique[is.na(lookuptry.unique$SPECIES),]
  lookuptry.notna$SPECIES=paste(lookuptry.notna$GENUS,lookuptry.notna$SPECIES,sep=" ")
  lookuptry.notna$SPECIES=lookuptry.notna$TAXON_NAME
  lookup.try=rbind(lookuptry.notna,lookuptry.na)
  lookup=lookup.try
  

  ########################################
  
  
  lookup$obs3=rownames(lookup)
  listcategory=c("Phylum","Class","SuperOrder","Order","SubOrder","Family","Subfamily","Tribe","Genus_group","Genus","Species_group")
  
  newdata_old=lookup[0,]
  siteinterest_1=lookup[0,]
  
  list2=lookup$obs3
  
  for (i in list2) {
    siteinterest=subset(lookup,lookup$obs3==i)
    
    for (j in listcategory) {
      IDcategory=j
      colnumber=which(colnames(siteinterest)==toupper(IDcategory))
      colnumberadd=colnumber+1
      
      if ((is.na (siteinterest[,colnumber])) == "TRUE") {
        NonNAindex <- which(!is.na(siteinterest[,colnumberadd:13]))
        lNonNAindex=length(NonNAindex)
        
        if (lNonNAindex > 0) {
          firstNonNA <- min(NonNAindex) #+ colnumber
          colnumberadd2=colnumber+firstNonNA
          siteinterest[,toupper(IDcategory)]=siteinterest[,colnumberadd2]
        }
      }
      
    }
    
    print(siteinterest$TAXON_NAME)
    siteinterest_all=rbind(siteinterest_1,siteinterest)
    siteinterest_1=siteinterest_all
  }
  
  lookup=siteinterest_1 #new lookuptable with NA issue fixed
  #inverts=inverts
  #return(lookup)
  #assign("invertsready",inverts,.GlobalEnv)
  assign("lookupsetup",lookup,.GlobalEnv)
}

