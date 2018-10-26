##by Christy Meredith cstarrmeredith@gmail.com

You will need to download the packages doBy and plyr to run this program

Create workingdirectory to location of input files. I would recommend creating a separate working directory for each run.
It must  have these files.
"invertinput.csv"
"lookup.csv"
"siteevent.csv"


###

The input files must  have these exact names. They should look exactly the same as in ambiguousrun_examples in terms of the same EXACT columns and column headings in the same order, or the program will not run

In lookup.csv, the lookup table,the taxonomic level should be "NA" if there is no level 

#copy folder "ambiguous program" from the ambiguousrun_examples folder into your working directory. DO NOT alter files in folder. To attempt edits, make a copy. 

###
Open program "ambiguoustaxaprogram.R"
set working directory to your working directory
Set the method = "numerous" if you want to the program to assign ambiguous taxa to the most numerous taxa
Set the method = "widespread" if you want the program to assign ambiguous taxa to the most widespread taxa (most sites)

Set the method_divide = "divide" if you want to divide among children at the site before reassigning; otherwise = "no_divide"

####
Select all the text in the "ambiguoustaxaprogram.R" and then select "Run". Program may take up to ten minutes to run. 

Created files:
lookup.site.csv; an edited lookup table which has assigned missing levels to the next level in the hierarchy
possibleerrors.csv; a table of unresolved entries due to spelling or other errors, corrections  may be needed in database

Within the site folder:
site_iter.csv; when resolving ambiguous taxa at a site, the taxa assigned at each iteration
sitelevel; for each record, the site-level taxa to which each ambiguous taxa was assigned is given
siteinvertsnew; the new dataset after the ambiguous taxa are resolved at the site level

Within the all folder:
all_iter.csv; the taxa assigned at each step in the iteration after the program has moved on to the level of the entire dataset
ambiguous_APTC_results.csv; for each initial record, including the sitelevel taxa assigned, and all level taxa assigned at the second step

##The remaining files in the "all" folder are the final files after ambiguous taxa have been resolved for the APTC_S, APTC_SG,APTC_SG1,APTC_SG2, and RPKC_s methods

# to the the list of singletons or doubletons or uniques and duplicates you have to look at the output files of the ambiguousprogram.R. WIthin R studio
this will be listed as "values" in the global environment.
For instance, for the APTC_SG method, the list of singletons will be called "singSG" and doubletons will be called "doubSG". The number of singletons 
will be "lsingSG" and the number of doubletons will be "ldoubSG"; just type in "singSG" or the variable you are interested in to get the result. Similar naming conventions are used for 
each of the ambiguous taxa methods