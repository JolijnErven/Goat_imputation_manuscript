# Load necessary libraries
library(stringr)
library(dplyr)

# This scripts returns the number of matching and discordant positions to calculate local concordance
# For Non reference concordance we printed matching and discordant non reference (het_homalt) positions
# and matching homozygous positions, for the FPR we printed the incorrect imputed heterozygotes positions

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

# Load validation dataset from CSV file from command line (output from validate_imputed_downsampled_vcfs.r)
val.GP<-data.frame(read.csv(args[1],header=F))

# Assign column names to the dataset
colnames(val.GP)<-c("chr","pos","REF","ALT","RAF","SAMPLE","TRUTH","NAME","COVERAGE","DATASET")

# Replacing "|" with "/" and normalizing genotype representation
val.GP$SAMPLE <- str_replace(str_replace(val.GP$SAMPLE,"\\Q|\\E","/"),"1/0","0/1")

# Create a new column VARIANT to represent the variant as a concatenation of REF and ALT alleles
val.GP$VARIANT<-paste0(val.GP$REF,val.GP$ALT)

# Initialize a column TYPE and classify variants as transitions or transversions
val.GP$TYPE<-"tranvs"   # Default type is transversion
val.GP[val.GP$VARIANT %in% c("GA","AG","CT","TC"),]$TYPE<-"trans"# Assign transition type

# Identify matching and discordant genotypes between SAMPLE and TRUTH
val.GP$MATCH<-"match" # Default assumption is a match
val.GP[val.GP$SAMPLE != val.GP$TRUTH,]$MATCH<-"discordant" # Mark discordant genotypes

# Script is memory intesive so run garbage collection to free memory
gc()


# Iterate through each MAF threshold
for (MAF in c(0.01, 0.05)) {
   # Assign MAF range to each row
   val.GP$MAF <- MAF
   # Filter variants within MAF threshold
   val.GP.05<-val.GP[val.GP$RAF >= MAF & val.GP$RAF <= (1-MAF),]
   
   #Free memory
   gc()
   
   # Filter for called sites (i.e., genotypes that are not missing)
   val.GP.05.sum.type<-data.frame(val.GP.05[val.GP.05$SAMPLE != "./." & val.GP.05$TRUTH != "./." ,])
   val.GP.05.sum.type$MAF <- MAF

   # Filter for matching variants
   val.GP.05.sum.type.dis<-val.GP.05.sum.type[val.GP.05.sum.type$MATCH=="discordant",]
   # Filter for discordant variants
   val.GP.05.sum.type.match<-val.GP.05.sum.type[val.GP.05.sum.type$MATCH!="discordant",]
   
   #Free memory
   gc()

   # Filter for correctly imputed homozygotes
   val.GP.05.sum.type.match.hom<-val.GP.05.sum.type.match[val.GP.05.sum.type.match$TRUTH=="1/1" | val.GP.05.sum.type.match$TRUTH == "0/0",]
   # Filter for correctly imputed non reference (non reference in Truth)
   val.GP.05.sum.type.match.het.homoalt<-val.GP.05.sum.type.match[val.GP.05.sum.type.match$TRUTH != "0/0",]
   # Filter for incorrectly imputed non reference (non reference in Truth)
   val.GP.05.sum.type.dis.het.homoalt<-val.GP.05.sum.type.dis[val.GP.05.sum.type.dis$TRUTH != "0/0",]
   # Filter for false heterozygous calls where the sample is "0/1" but the truth is homozygous ("0/0" or "1/1")
   val.GP.05.type.falseHet<-data.frame(val.GP.05[(val.GP.05$SAMPLE == "0/1") & (val.GP.05$TRUTH == "0/0" | val.GP.05$TRUTH == "1/1" ),])

   # Write processed data to CSV files with appropriate naming conventions
   write.csv(x = val.GP.05.sum.type.match.hom,file=paste0(args[1],"MAF_",MAF,"_Match_homozygous_positions.csv"), quote = F, row.names = F)
   write.csv(x = val.GP.05.sum.type.dis.het.homoalt,file=paste0(args[1],"MAF_",MAF,"_Discordant_Het_homalt_positions.csv"), quote = F, row.names = F)
   write.csv(x = val.GP.05.sum.type.match.het.homoalt,file=paste0(args[1],"MAF_",MAF,"_Match_Het_homoalt_positions.csv"), quote = F, row.names = F)
   write.csv(x = val.GP.05.type.falseHet,file=paste0(args[1],"MAF_",MAF,"_False_heterozygous_positions.csv"), quote = F, row.names = F)
   
   #Free memory
   gc()
}
