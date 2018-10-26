

#####################package 2 missing levels
#3 
####################

#' A missing levels function
#'
#' this part of program makes a unique taxon_name for intermediate levels present in the dataset that currently
##do not have a unique taxon_name in the lookup table
#' @param 
#' @keywords 
#' @export
#' @examples
#' missinglevels()



missinglevels<-function(lookupsetup){

  lookuptry.unique=lookupsetup
  
  #phylum level
  lookup.try.phylum=lookuptry.unique
  lookup.try.phylum$TAXON_NAME=lookup.try.phylum$PHYLUM
  lookup.try.phylum=lookup.try.phylum[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.phylum[,3:13]=NA
  lookup.try.phylum=unique(lookup.try.phylum)
  lookup.try.phylum$IDLevel="Phylum"
  lookup.try.phylum$no=1
  
  ##class level
  lookup.try.class=lookuptry.unique
  lookup.try.class$TAXON_NAME=lookup.try.class$CLASS
  lookup.try.class=lookup.try.class[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.class[,4:13]=NA
  lookup.try.class=unique(lookup.try.class)
  lookup.try.class$IDLevel="Class"
  lookup.try.class$no=2
  
  #########superorder level
  lookup.try.superorder=lookuptry.unique
  lookup.try.superorder$TAXON_NAME=lookup.try.superorder$SUPERORDER
  lookup.try.superorder=lookup.try.superorder[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.superorder[,5:13]=NA
  lookup.try.superorder=unique(lookup.try.superorder)
  lookup.try.superorder$IDLevel="SuperOrder"
  lookup.try.superorder$no=3
  
  ###order level
  lookup.try.order=lookuptry.unique
  lookup.try.order$TAXON_NAME=lookup.try.order$ORDER
  lookup.try.order=lookup.try.order[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.order[,6:13]=NA
  lookup.try.order=unique(lookup.try.order)
  lookup.try.order$IDLevel="Order"
  lookup.try.order$no=4
  
  ###suborder level
  lookup.try.suborder=lookuptry.unique
  lookup.try.suborder$TAXON_NAME=lookup.try.suborder$SUBORDER
  lookup.try.suborder=lookup.try.suborder[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.suborder[,7:13]=NA
  lookup.try.suborder=unique(lookup.try.suborder)
  lookup.try.suborder$IDLevel="SubOrder"
  lookup.try.suborder$no=5
  
  ######family level
  lookup.try.family=lookuptry.unique
  lookup.try.family$TAXON_NAME=lookup.try.family$FAMILY
  lookup.try.family=lookup.try.family[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.family[,8:13]=NA
  lookup.try.family=unique(lookup.try.family)
  lookup.try.family$IDLevel="Family"
  lookup.try.family$no=6
  
  ###subfamily level
  lookup.try.subfamily=lookuptry.unique
  lookup.try.subfamily$TAXON_NAME=lookup.try.subfamily$SUBFAM
  lookup.try.subfamily=lookup.try.subfamily[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.subfamily[,9:13]=NA
  lookup.try.subfamily=unique(lookup.try.subfamily)
  lookup.try.subfamily$IDLevel="Subfamily"
  lookup.try.subfamily$no=7
  
  ########################## tribe level
  lookup.try.tribe=lookuptry.unique
  lookup.try.tribe$TAXON_NAME=lookup.try.tribe$TRIBE
  lookup.try.tribe=lookup.try.tribe[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.tribe[,10:13]=NA
  lookup.try.tribe=unique(lookup.try.tribe)
  lookup.try.tribe$IDLevel="Tribe"
  lookup.try.tribe$no=8
  
  #######################genus group level
  lookup.try.Genus_group=lookuptry.unique
  lookup.try.Genus_group$TAXON_NAME=lookup.try.Genus_group$GENUS_GROUP
  lookup.try.Genus_group=lookup.try.Genus_group[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.Genus_group[,11:13]=NA
  lookup.try.Genus_group=unique(lookup.try.Genus_group)
  lookup.try.Genus_group$IDLevel="Genus_group"
  lookup.try.Genus_group$no=9
  
  ##############genus level
  lookup.try.Genus=lookuptry.unique
  lookup.try.Genus$TAXON_NAME=lookup.try.Genus$GENUS
  lookup.try.Genus=lookup.try.Genus[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.Genus[,12:13]=NA
  lookup.try.Genus=unique(lookup.try.Genus)
  lookup.try.Genus$IDLevel="Genus"
  lookup.try.Genus$no=10
  
  ##speciesgroup
  lookup.try.Species_group=lookuptry.unique
  lookup.try.Species_group$TAXON_NAME=lookup.try.Species_group$SPECIES_GROUP
  lookup.try.Species_group=lookup.try.Species_group[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  lookup.try.Species_group[,13:13]=NA
  lookup.try.Species_group=unique(lookup.try.Species_group)
  lookup.try.Species_group$IDLevel="Species_group"
  lookup.try.Species_group$no=11
  
  ####species
  lookup.try.Species=lookuptry.unique
  lookup.try.Species$TAXON_NAME=lookup.try.Species$SPECIES
  lookup.try.Species=lookup.try.Species[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
  #lookup.try.Species[,11:11]=NA
  lookup.try.Species=unique(lookup.try.Species)
  lookup.try.Species$IDLevel="Species"
  lookup.try.Species$no=12
  
  merged=rbind(lookup.try.order,lookup.try.phylum,lookup.try.class,lookup.try.family,lookup.try.subfamily,lookup.try.tribe,lookup.try.Genus_group,lookup.try.Genus,lookup.try.Species_group,lookup.try.Species,lookup.try.superorder,lookup.try.suborder)
  merge1=unique(merged)
  merge1=merge1[!is.na(merge1$TAXON_NAME),]
  lookup=merge1 #new lookup table
  lookup$IDLevel=toupper(lookup$IDLevel)
  dir.create("./output")
  write.table(lookup,"./output/lookup.site.csv",sep=",",row.names=FALSE) 
  #write.table(lookup,"lookup.site.csv",sep=",",row.names = FALSE)
  
  #assign("inverts2",inverts,.GlobalEnv)
  
  assign("lookupready",lookup,.GlobalEnv)
  
  
}
