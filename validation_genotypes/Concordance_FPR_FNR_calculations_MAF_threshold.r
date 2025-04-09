# Load necessary libraries
library(stringr)
library(dplyr)

# This scripts calculates the concordance of Non reference, heterozygous, homozygous alternative, homozygous 
# false-positive and false negative heterozygous rate for each downsampled coverage, GP threshold and sample
# for MAF thresholds
# Created by Kevin Daly updated by Jolijn Erven

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

# Create dataframes
recov = data.frame()
df = data.frame()
df2 = data.frame()
df3 = data.frame()
df4 = data.frame()
df.het = data.frame()
df3.het = data.frame()
df.FN1 = data.frame()
df.FN2 = data.frame()
df5 = data.frame()
df6 = data.frame()
df7 = data.frame()
df8 = data.frame()
df.count = data.frame()
df.count.type=data.frame()
df9 = data.frame()
df10 = data.frame()
df.data = data.frame()


# Iterate through each MAF threshold
for (MAF in c(0.00,0.01, 0.05, 0.10)) {

   # Assign MAF range to each row
   val.GP$MAF <- MAF
   
   # Filter variants within MAF threshold
   val.GP.05<-val.GP[val.GP$RAF >= MAF & val.GP$RAF <= (1-MAF),]
   
   #get number of genotype imputed
   val.GP.05.recov <- val.GP.05 %>% group_by(NAME,COVERAGE,SAMPLE,DATASET) %>% summarise(n = n())
   val.GP.05.recov$MAF <- MAF
   # Append to dataframe
   recov<-rbind(recov, val.GP.05.recov)
   
   #Free memory
   gc()
   
   # Count the number of called variants per sample, coverage, and dataset for both all sites and transitions and transverions seperately (Type)
   val.GP.05.called<-val.GP.05[val.GP.05$SAMPLE != "./." ,]  %>% group_by(NAME,COVERAGE,DATASET) %>% summarise(n = n())
   val.GP.05.called$MAF<-MAF
   val.GP.05.called.type<-val.GP.05[val.GP.05$SAMPLE != "./." ,]  %>% group_by(NAME,COVERAGE,TYPE,DATASET) %>% summarise(n = n())
   val.GP.05.called.type$MAF<-MAF
   # Count the number of called variants for imputed and truth per sample, coverage, and dataset for both all sites and transitions and transverions seperately (Type)
   val.GP.05.sum<-data.frame(val.GP.05[val.GP.05$SAMPLE != "./." & val.GP.05$TRUTH != "./.",] %>% group_by(NAME,COVERAGE,MATCH,TRUTH,DATASET) %>% summarise(n = n()))
   val.GP.05.sum.type<-data.frame(val.GP.05[val.GP.05$SAMPLE != "./." & val.GP.05$TRUTH != "./." ,] %>% group_by(NAME,COVERAGE,TYPE,MATCH,TRUTH,DATASET) %>% summarise(n = n()))
   val.GP.05.sum.type$MAF <- MAF
   
   #Free memory
   gc()
   
   # Separate discordant and matching variants for both all sites and transitions and transverions seperately (Type)
   # All sites
   val.GP.05.sum.dis<-val.GP.05.sum[val.GP.05.sum$MATCH=="discordant",]
   val.GP.05.sum.match<-val.GP.05.sum[val.GP.05.sum$MATCH!="discordant",]
   val.GP.05.sum.type.dis<-val.GP.05.sum.type[val.GP.05.sum.type$MATCH=="discordant",]
   val.GP.05.sum.type.match<-val.GP.05.sum.type[val.GP.05.sum.type$MATCH!="discordant",]
   
   # Compute R-squared and concordance rates for both all sites and transitions and transverions seperately (Type)
   val.GP.05.sum.type.match$r2<-(val.GP.05.sum.type.match$n / ( val.GP.05.sum.type.dis$n + val.GP.05.sum.type.match$n))^2
   val.GP.05.sum.match$r2<-(val.GP.05.sum.match$n / ( val.GP.05.sum.dis$n + val.GP.05.sum.match$n))^2
   val.GP.05.sum.type.match$concord<-(val.GP.05.sum.type.match$n / ( val.GP.05.sum.type.dis$n + val.GP.05.sum.type.match$n))
   val.GP.05.sum.match$concord<-(val.GP.05.sum.match$n / ( val.GP.05.sum.dis$n + val.GP.05.sum.match$n))
   
   # Filter for heterozygous genotypes (0/1)
   val.GP.05.sum.type.match.het<-val.GP.05.sum.type.match[val.GP.05.sum.type.match$TRUTH=="0/1",]
   val.GP.05.sum.match.het<-val.GP.05.sum.match[val.GP.05.sum.match$TRUTH=="0/1",]
   val.GP.05.sum.match.het.transv<-val.GP.05.sum.type.match.het[val.GP.05.sum.type.match.het$TYPE!="trans",]
   val.GP.05.sum.match.het.trans<-val.GP.05.sum.type.match.het[val.GP.05.sum.type.match.het$TYPE=="trans",]
   
   #Free memory
   gc()
   
   # Filter for non-reference (0/1, 1/1)
   val.GP.05.sum.het.homoalt<-val.GP.05.sum[val.GP.05.sum$TRUTH != "0/0",] %>% group_by(NAME,COVERAGE,MATCH,DATASET)  %>% summarise(N = sum(n))
   val.GP.05.sum.het.homoalt$MAF<-MAF
   # Separate discordant and matching variants
   val.GP.05.sum.het.homoalt.dis<-val.GP.05.sum.het.homoalt[val.GP.05.sum.het.homoalt$MATCH=="discordant",]
   val.GP.05.sum.het.homoalt.match<-val.GP.05.sum.het.homoalt[val.GP.05.sum.het.homoalt$MATCH!="discordant",]
   # Compute R-squared and concordance rates
   val.GP.05.sum.het.homoalt.match$r2<-(val.GP.05.sum.het.homoalt.match$N / ( val.GP.05.sum.het.homoalt.match$N + val.GP.05.sum.het.homoalt.dis$N))^2
   val.GP.05.sum.het.homoalt.match$concord<-(val.GP.05.sum.het.homoalt.match$N / ( val.GP.05.sum.het.homoalt.match$N + val.GP.05.sum.het.homoalt.dis$N))
   # Filter for non-reference (0/1, 1/1) for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het.homoalt<-val.GP.05.sum.type[val.GP.05.sum.type$TRUTH != "0/0",] %>% group_by(NAME,COVERAGE,TYPE,DATASET,MATCH,MAF)  %>% summarise(N = sum(n))
   # Separate discordant and matching variants for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het.homoalt.dis<-val.GP.05.sum.type.het.homoalt[val.GP.05.sum.type.het.homoalt$MATCH=="discordant",]
   val.GP.05.sum.type.het.homoalt.match<-val.GP.05.sum.type.het.homoalt[val.GP.05.sum.type.het.homoalt$MATCH!="discordant",]
   # Compute R-squared and concordance rates for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het.homoalt.match$r2<-(val.GP.05.sum.type.het.homoalt.match$N / ( val.GP.05.sum.type.het.homoalt.match$N + val.GP.05.sum.type.het.homoalt.dis$N))^2
   val.GP.05.sum.type.het.homoalt.match$concord<-(val.GP.05.sum.type.het.homoalt.match$N / ( val.GP.05.sum.type.het.homoalt.match$N + val.GP.05.sum.type.het.homoalt.dis$N)
   
   #Free memory
   gc()
   
   # Filter for heterozygotes (0/1)
   val.GP.05.sum.het<-val.GP.05.sum[val.GP.05.sum$TRUTH == "0/1",] %>% group_by(NAME,COVERAGE,MATCH,DATASET)  %>% summarise(N = sum(n))
   val.GP.05.sum.het$MAF<-MAF
   # Separate discordant and matching variants
   val.GP.05.sum.het.dis<-val.GP.05.sum.het[val.GP.05.sum.het$MATCH=="discordant",]
   val.GP.05.sum.het.match<-val.GP.05.sum.het[val.GP.05.sum.het$MATCH!="discordant",]
   # Compute R-squared and concordance rates
   val.GP.05.sum.het.match$r2<-(val.GP.05.sum.het.match$N / ( val.GP.05.sum.het.match$N + val.GP.05.sum.het.dis$N))^2
   val.GP.05.sum.het.match$concord<-(val.GP.05.sum.het.match$N / ( val.GP.05.sum.het.match$N + val.GP.05.sum.het.dis$N))
   # Filter for heterozygotes (0/1) for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het<-val.GP.05.sum.type[val.GP.05.sum.type$TRUTH == "0/1",] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF)  %>% summarise(N = sum(n))
   # Separate discordant and matching variants for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het.dis<-val.GP.05.sum.type.het[val.GP.05.sum.type.het$MATCH=="discordant",]
   val.GP.05.sum.type.het.match<-val.GP.05.sum.type.het[val.GP.05.sum.type.het$MATCH!="discordant",]
   # Compute R-squared and concordance rates for transitions and transverions seperately (Type)
   val.GP.05.sum.type.het.match$r2<-(val.GP.05.sum.type.het.match$N / ( val.GP.05.sum.type.het.match$N + val.GP.05.sum.type.het.dis$N))^2
   val.GP.05.sum.type.het.match$concord<-(val.GP.05.sum.type.het.match$N / ( val.GP.05.sum.type.het.match$N + val.GP.05.sum.type.het.dis$N))
   
   #Free memory
   gc()
   
   # Filter for homozygotes (0/0, 1/1)
   val.GP.05.sum.hom<-val.GP.05.sum[val.GP.05.sum$TRUTH == "0/0" | val.GP.05.sum$TRUTH == "1/1" ,] %>% group_by(NAME,COVERAGE,MATCH,DATASET)  %>% summarise(N = sum(n))
   val.GP.05.sum.hom$MAF<-MAF
   # Separate discordant and matching variants
   val.GP.05.sum.hom.dis<-val.GP.05.sum.hom[val.GP.05.sum.hom$MATCH=="discordant",]
   val.GP.05.sum.hom.match<-val.GP.05.sum.hom[val.GP.05.sum.hom$MATCH!="discordant",]
   # Compute R-squared and concordance rates
   val.GP.05.sum.hom.match$r2<-(val.GP.05.sum.hom.match$N / ( val.GP.05.sum.hom.match$N + val.GP.05.sum.hom.dis$N))^2
   val.GP.05.sum.hom.match$concord<-(val.GP.05.sum.hom.match$N / ( val.GP.05.sum.hom.match$N + val.GP.05.sum.hom.dis$N))
   # Filter for homozygotes (0/0, 1/1) for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom<-val.GP.05.sum.type[val.GP.05.sum.type$TRUTH == "0/0" | val.GP.05.sum.type$TRUTH == "1/1",] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF)  %>% summarise(N = sum(n))
   # Separate discordant and matching variants for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom.dis<-val.GP.05.sum.type.hom[val.GP.05.sum.type.hom$MATCH=="discordant",]
   val.GP.05.sum.type.hom.match<-val.GP.05.sum.type.hom[val.GP.05.sum.type.hom$MATCH!="discordant",]
   # Compute R-squared and concordance rates for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom.match$r2<-(val.GP.05.sum.type.hom.match$N / ( val.GP.05.sum.type.hom.match$N + val.GP.05.sum.type.hom.dis$N))^2
   val.GP.05.sum.type.hom.match$concord<-(val.GP.05.sum.type.hom.match$N / ( val.GP.05.sum.type.hom.match$N + val.GP.05.sum.type.hom.dis$N))
   
   # Free memory
   gc()
    
	# Filter for homozygotes alternative (1/1)
   val.GP.05.sum.hom.alt<-val.GP.05.sum[val.GP.05.sum$TRUTH == "1/1" ,] %>% group_by(NAME,COVERAGE,MATCH,DATASET)  %>% summarise(N = sum(n))
   val.GP.05.sum.hom.alt$MAF<-MAF
   # Separate discordant and matching variants
   val.GP.05.sum.hom.alt.dis<-val.GP.05.sum.hom.alt[val.GP.05.sum.hom.alt$MATCH=="discordant",]
   val.GP.05.sum.hom.alt.match<-val.GP.05.sum.hom.alt[val.GP.05.sum.hom.alt$MATCH!="discordant",]
   # Compute R-squared and concordance rates
   val.GP.05.sum.hom.alt.match$r2<-(val.GP.05.sum.hom.alt.match$N / ( val.GP.05.sum.hom.alt.match$N + val.GP.05.sum.hom.alt.dis$N))^2
   val.GP.05.sum.hom.alt.match$concord<-(val.GP.05.sum.hom.alt.match$N / ( val.GP.05.sum.hom.alt.match$N + val.GP.05.sum.hom.alt.dis$N))
   # Filter for homozygotes alternative(1/1) for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom.alt<-val.GP.05.sum.type[val.GP.05.sum.type$TRUTH == "1/1",] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF)  %>% summarise(N = sum(n))
   # Separate discordant and matching variants for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom.alt.dis<-val.GP.05.sum.type.hom.alt[val.GP.05.sum.type.hom.alt$MATCH=="discordant",]
   val.GP.05.sum.type.hom.alt.match<-val.GP.05.sum.type.hom.alt[val.GP.05.sum.type.hom.alt$MATCH!="discordant",]
   # Compute R-squared and concordance rates for transitions and transverions seperately (Type)
   val.GP.05.sum.type.hom.alt.match$r2<-(val.GP.05.sum.type.hom.alt.match$N / ( val.GP.05.sum.type.hom.alt.match$N + val.GP.05.sum.type.hom.alt.dis$N))^2
   val.GP.05.sum.type.hom.alt.match$concord<-(val.GP.05.sum.type.hom.alt.match$N / ( val.GP.05.sum.type.hom.alt.match$N + val.GP.05.sum.type.hom.alt.dis$N))
   
   # Free memory
   gc()

   # Append results to respective data frames
   df.data<-rbind(df.data, val.GP.05.sum.type)
   df<-rbind(df,val.GP.05.sum.het.homoalt.match)
   df3<-rbind(df3,val.GP.05.sum.type.het.homoalt.match)
   df5<-rbind(df5,val.GP.05.sum.hom.match)
   df6<-rbind(df6,val.GP.05.sum.type.hom.match)
   df.count<-rbind(df.count,val.GP.05.called)
   df.count.type<-rbind(df.count.type,val.GP.05.called.type)
   df7<-rbind(df7,val.GP.05.sum.hom.alt.match)
   df8<-rbind(df8,val.GP.05.sum.type.hom.alt.match)
   df9<-rbind(df9,val.GP.05.sum.type.match)
   df10<-rbind(df10,val.GP.05.sum.match)
   df.het<-rbind(df.het,val.GP.05.sum.het.match)
   df3.het<-rbind(df3.het,val.GP.05.sum.type.het.match)


   # Summarize incorrectly imputed heterozygotes (false heterozygotes)
   val.GP.05.falseHet<-data.frame(val.GP.05[val.GP.05$SAMPLE == "0/1" & (val.GP.05$TRUTH == "0/0" | val.GP.05$TRUTH == "1/1" ),] %>% group_by(NAME,COVERAGE,MATCH,DATASET,MAF) %>% summarise(n = n()))
   # Summarize all imputed heterozygotes where truth is called
   val.GP.05.allHet<-data.frame(val.GP.05[val.GP.05$SAMPLE == "0/1" & val.GP.05$TRUTH != "./.",] %>% group_by(NAME,COVERAGE,DATASET,MAF) %>% summarise(n = n()))
   # Compute the false heterozygous rate
   val.GP.05.falseHet$falseHetRate<-val.GP.05.falseHet$n/val.GP.05.allHet$n
   # Compute the false positive rate (FPR) for false heterozygous calls
   val.GP.05.falseHet$FPR<-100*(val.GP.05.falseHet$n/(val.GP.05.sum.hom.match$N + val.GP.05.falseHet$n))
   # Summarize incorrectly imputed heterozygotes (false heterozygotes) for transitions and transverions seperately (Type)
   val.GP.05.type.falseHet<-data.frame(val.GP.05[(val.GP.05$SAMPLE == "0/1") & (val.GP.05$TRUTH == "0/0" | val.GP.05$TRUTH == "1/1" ),] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF) %>% summarise(n = n()))
   # Summarize all imputed heterozygotes where truth is called for transitions and transverions seperately (Type)
   val.GP.05.type.allHet<-data.frame(val.GP.05[val.GP.05$TRUTH == "0/1" & val.GP.05$SAMPLE != "./.",] %>% group_by(NAME,COVERAGE,DATASET,MAF,TYPE) %>% summarise(n = n()))
   # Compute the false heterozygous rate for transitions and transverions seperately (Type)
   val.GP.05.type.falseHet$falseHetRate<-val.GP.05.type.falseHet$n/val.GP.05.type.allHet$n
   
   # Try computing FPR and handle errors if division by zero occurs for transitions and transverions seperately (Type)
   mtry <- try(val.GP.05.type.falseHet$FPR<-100*(val.GP.05.type.falseHet$n/(val.GP.05.sum.type.hom.match$N + val.GP.05.type.falseHet$n)))
   if (!inherits(mtry, "try-error")) {
     val.GP.05.type.falseHet$FPR<-100*(val.GP.05.type.falseHet$n/(val.GP.05.sum.type.hom.match$N + val.GP.05.type.falseHet$n))}
   else {
    val.GP.05.type.falseHet$FPR<- "NA"
   }

   # Summarize incorrectly imputed sites that are heterozygous in TRUTH
   val.GP.05.falseHetNeg<-data.frame(val.GP.05[(val.GP.05$SAMPLE == "0/0"  | val.GP.05$SAMPLE == "1/1") & (val.GP.05$TRUTH == "0/1"  ),] %>% group_by(NAME,COVERAGE,MATCH,DATASET,MAF) %>% summarise(n = n()))
   # Summarize all heterozygotes in the TRUTH dataset which are not missing in imputed
   val.GP.05.allHet.Neg<-data.frame(val.GP.05[val.GP.05$TRUTH == "0/1" & val.GP.05$SAMPLE != "./.",] %>% group_by(NAME,COVERAGE,DATASET,MAF) %>% summarise(n = n()))
   # Compute the false negative heterozygous rate
   val.GP.05.falseHetNeg$falseHetRate<-val.GP.05.falseHetNeg$n/val.GP.05.allHet.Neg$n
   # Summarize incorrectly imputed sites that are heterozygous in TRUTH for transitions and transverions seperately (Type)
   val.GP.05.type.falseHetNeg<-data.frame(val.GP.05[(val.GP.05$SAMPLE == "0/0"  | val.GP.05$SAMPLE == "1/1") & (val.GP.05$TRUTH == "0/1"  ),] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF) %>% summarise(n = n()))
   # Summarize all heterozygotes in the TRUTH dataset which are not missing in imputed for transitions and transverions seperately (Type)
   val.GP.05.type.allHet.Neg<-data.frame(val.GP.05[val.GP.05$TRUTH == "0/1" & val.GP.05$SAMPLE != "./.",] %>% group_by(NAME,COVERAGE,TYPE,MATCH,DATASET,MAF) %>% summarise(n = n()))
   # Compute the false negative heterozygous rate for transitions and transverions seperately (Type)
   val.GP.05.type.falseHetNeg$falseHetRate<-val.GP.05.type.falseHetNeg$n/val.GP.05.type.allHet.Neg$n
   
   # Free memory
   gc()
    
   # Append results to respective data frames
   df2<-rbind(df2,val.GP.05.falseHet)
   df4<-rbind(df4,val.GP.05.type.falseHet)
   df.FN1<-rbind(df.FN1,val.GP.05.falseHetNeg)
   df.FN2<-rbind(df.FN2,val.GP.05.type.falseHetNeg)

}

# write results to respective files
write.csv(x = df,file=paste0(args[1],"_MAF_threshold_concordance.csv"), quote = F, row.names = F)
write.csv(x = df2,file=paste0(args[1],"_MAF_threshold_falseHetRate.csv"), quote = F, row.names = F)
write.csv(x = df3,file=paste0(args[1],"_MAF_threshold_concordance_type.csv"), quote = F, row.names = F)
write.csv(x = df4,file=paste0(args[1],"_MAF_threshold_falseHetRate_type.csv"), quote = F, row.names = F)
write.csv(x = df.het,file=paste0(args[1],"_MAF_threshold_het_concordance.csv"), quote = F, row.names = F)
write.csv(x = df3.het,file=paste0(args[1],"_MAF_threshold_het_concordance_type.csv"), quote = F, row.names = F)
write.csv(x = df5,file=paste0(args[1],"_MAF_threshold_concordance_hom.csv"), quote = F, row.names = F)
write.csv(x = df6,file=paste0(args[1],"_MAF_threshold_concordance_hom_type.csv"), quote = F, row.names = F)
write.csv(x = df.count,file=paste0(args[1],"_MAF_threshold_count.csv"), quote = F, row.names = F)
write.csv(x = df.count.type,file=paste0(args[1],"_MAF_threshold_count_type.csv"), quote = F, row.names = F)
write.csv(x = df7,file=paste0(args[1],"_MAF_threshold_concordance_hom_alt.csv"), quote = F, row.names = F)
write.csv(x = df8,file=paste0(args[1],"_MAF_threshold_concordance_hom_alt_type.csv"), quote = F, row.names = F)
write.csv(x = df10,file=paste0(args[1],"_MAF_threshold_concordance_hom_all.csv"), quote = F, row.names = F)
write.csv(x = df9,file=paste0(args[1],"_MAF_threshold_concordance_hom_all_type.csv"), quote = F, row.names = F)
write.csv(x = df.FN1,file=paste0(args[1],"_MAF_threshold_false_negHetRate.csv"), quote = F, row.names = F)
write.csv(x = df.FN2,file=paste0(args[1],"_MAF_threshold_false_negHetRate_type.csv"), quote = F, row.names = F)
write.csv(x = df.data,file=paste0(args[1],"_MAF_threshold_genotype_information.csv"), quote = F, row.names = F)

#recovery
write.csv(x = recov,file=paste0(args[1],"_all_sites_geno_information_MAF_threshold.csv"), quote = F, row.names = F)
