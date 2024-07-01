# Clear environment variables and load necessary libraries
rm(list=ls()) # Remove all objects in the workspace
library(tidyverse) # Load tidyverse package for data manipulation and visualization
library(vcfR) # Load vcfR package for handling VCF files
library(UpSetR) # Load UpSetR package for creating UpSet plots
library(ggplot2) # Load ggplot2 package for advanced plotting
library(dplyr) # Load dplyr package for data manipulation
library(limma) # Load limma package for linear models
library(VennDiagram) # Load VennDiagram package for creating Venn diagrams
library(RColorBrewer) # Load RColorBrewer package for color palettes
library(ggVennDiagram) # Load ggVennDiagram package for creating Venn diagrams with ggplot2
library(maftools) # Load maftools package for analyzing cancer genomics data
library(viridis) # Load viridis package for color palettes
library(readxl) # Load readxl package for reading Excel files
library(gplots) # Load gplots package for various types of plots

# Define categories of chromosome ruptures
early_rupture=c('13', '18', '21') # Early rupture categories
late_rupture=c('1','2','3','4','5','6','17','19','X') # Late rupture categories
mid_rupture=c('7','8','9','10','11','12','14','15','16','20','22') # Mid rupture categories

# Read and preprocess chromothripsis data from an Excel file
chromothripsis=read_excel('Path to Supplementary table from Cortes-Ciriano et al Nature Genetic 2020') # Replace with actual path
chromothripsis=chromothripsis[chromothripsis$type_chromothripsis!='With other complex events',] # Filter out rows with certain conditions
chromothripsis=chromothripsis[chromothripsis$chromo_label!='Linked to low confidence',] # Further filter based on chromo_label
chromothripsis=chromothripsis[chromothripsis$chromo_label!='Linked to high confidence',] # Additional filtering
chromothripsis=chromothripsis[chromothripsis$chromo_label!='Low confidence',] # Final filtering
chromothripsis$chromothripsis=ifelse(chromothripsis$chromo_label == 'High confidence','chromothripsis','no') # Assign values based on condition
chromothripsis$rupture=ifelse(chromothripsis$Chr %in% early_rupture,'EarlyRupture',ifelse(chromothripsis$Chr %in% mid_rupture, 'MidRupture','LateRupture')) # Assign rupture type
chromothripsis$status=paste(chromothripsis$rupture,chromothripsis$chromothripsis,sep='_') # Combine rupture and chromothripsis labels
colnames(chromothripsis)[2]='chromosome' # Rename column

# Process single nucleotide variant (SNV) data
snv_list='List of paths to downloaded.tsv.gz files containing SNVs in the PCAWG cohort' # Placeholder for file paths
for(k in 1:length(snv_list)){
  snv=as.data.frame(fread(snv_list[k],header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA"))) # Read SNV data
  snv=snv[snv$mutation_type=="single base substitution",] # Filter for single base substitutions
  snv=snv[snv$variation_calling_algorithm=='PCAWG Consensus SNV-MNV caller',] # Filter based on variation calling algorithm
  a=unique(snv$icgc_donor_id) # Get unique donor IDs
  a=a[a %in% chromothripsis$icgc_donor_id] # Filter donors present in chromothripsis data
  snv=snv[snv$icgc_donor_id %in% a,] # Filter SNV data based on donor ID
  snv1=merge(snv, chromothripsis,by=c('icgc_donor_id','chromosome')) # Merge SNV and chromothripsis data
  snv1$sample=paste(snv1$icgc_donor_id,snv1$status,sep='_') # Create sample identifier
  a3=unique(snv1$sample) # Get unique samples
  for (i in 1:length(a3)){ # Loop through each sample
    trim=snv1[snv1$sample==a3[i],] # Subset data for current sample
    trim=trim[,c('chromosome','chromosome_start','sample','mutated_from_allele','mutated_to_allele')] # Select relevant columns
    colnames(trim)=c('Chr','Start','sample_name','Ref','Alt') # Rename columns
    trim$Chr=paste0('chr',trim$Chr) # Prefix chromosome numbers with 'chr'
    trim=trim[duplicated(trim$Start),] # Remove duplicate entries
    if(nrow(trim)>10){ # Conditionally write trimmed data to file if row count exceeds 10
      write.table(trim,file=paste0('/fh/fast/ha_g/projects/PCAWG/High_Confidence_not_complex_events/',a3[i],'.vcf'), quote = FALSE,sep="\t",col.names = FALSE,row.names = FALSE)
    }
  }
}

# Define mutational signatures
mutsig=c('SBS1'= 'Age',
         'SBS2'= 'APOBEC',
         'SBS3'= 'HR deficiency',
         'SBS4'= 'Tobacco consumption',
         'SBS5'= 'Age',
         'SBS6'= 'Defective DNA mismatch repair',
         'SBS7a'= 'Ultraviolet light exposure',
         'SBS7b'= 'Ultraviolet light exposure',
         'SBS7c'= 'Ultraviolet light exposure',
         'SBS7d'= 'Ultraviolet light exposure',
         'SBS8'= 'HR deficiency',
         'SBS9'= 'Polymerase eta somatic hypermutation activity',
         'SBS10a'= 'Polymerase epsilon exonuclease domain mutations',
         'SBS10b'= 'Polymerase epsilon exonuclease domain mutations',
         'SBS10c'= 'Defective POLD1 proofreading',
         'SBS10d'= 'Defective POLD1 proofreading',
         'SBS11'= 'Chemotherapy treatment',
         'SBS12'= 'Unknown',
         'SBS13'= 'APOBEC',
         'SBS14'= 'Defective DNA mismatch repair',
         'SBS15'= 'Defective DNA mismatch repair',
         'SBS16'= 'Unknown',
         'SBS17a'= 'Unknown',
         'SBS17b'= 'Unknown',
         'SBS18'= 'Damage by ROS',
         'SBS19'= 'Unknown',
         'SBS20'= 'Defective DNA mismatch repair',
         'SBS21'= 'Defective DNA mismatch repair',
         'SBS22a'= 'Aristolochic acid exposure',
         'SBS23'= 'Unknown',
         'SBS24'= 'Aflatoxin exposure',
         'SBS25'= 'Chemotherapy treatment',
         'SBS26'= 'Defective DNA mismatch repair',
         'SBS27'= 'Sequencing Artifacts',
         'SBS28'= 'Unknown',
         'SBS29'= 'Tobacco consumption',
         'SBS30'= 'BER Deficiency',
         'SBS31'= 'Chemotherapy treatment',
         'SBS32'= 'Immunosuppressive drug treatment',
         'SBS33'= 'Unknown',
         'SBS34'= 'Unknown',
         'SBS35'= 'Chemotherapy treatment',
         'SBS36'= 'BER Deficiency',
         'SBS37'= 'Unknown',
         'SBS38'= 'Ultraviolet light exposure',
         'SBS39'= 'Unknown',
         'SBS40a'= 'Unknown',
         'SBS40b'= 'Unknown',
         'SBS40c'= 'Unknown',
         'SBS41'= 'Unknown',
         'SBS42'= 'Haloalkane exposure',
         'SBS43'= 'Sequencing Artifacts',
         'SBS44'= 'Defective DNA mismatch repair',
         'SBS45'= 'Sequencing Artifacts',
         'SBS46'= 'Sequencing Artifacts',
         'SBS47'= 'Sequencing Artifacts',
         'SBS48'= 'Sequencing Artifacts',
         'SBS49'= 'Sequencing Artifacts',
         'SBS50'= 'Sequencing Artifacts',
         'SBS51'= 'Sequencing Artrops',
         'SBS52'= 'Sequencing Artifacts',
         'SBS53'= 'Sequencing Artifacts',
         'SBS54'= 'Sequencing Artifacts',
         'SBS55'= 'Sequencing Artifacts',
         'SBS56'= 'Sequencing Artifacts',
         'SBS57'= 'Sequencing Artifacts',
         'SBS58'= 'Sequencing Artifacts',
         'SBS59'= 'Sequencing Artifacts',
         'SBS60'= 'Sequencing Artifacts',
         'SBS84'= 'LymphoidCancer',
         'SBS85'= 'LymphoidCancer',
         'SBS86'= 'Chemotherapy treatment',
         'SBS87'= 'Immunosuppressive drug treatment',
         'SBS88'= 'Colibactin exposure',
         'SBS89'= 'Unknown',
         'SBS90'= 'Chemotherapy treatment',
         'SBS91'= 'Unknown',
         'SBS92'= 'Tobacco consumption',
         'SBS93'= 'Unknown',
         'SBS94'= 'Unknown',
         'SBS95'= 'Sequencing Artifacts',
         "SBS96" = "Unknown",
         "SBS97" = "Unknown",
         "SBS98" = "Unknown",
         "SBS99" = "Chemotherapy treatment")

# Convert mutational signature dictionary to a data frame
mutsig=(as.data.frame(mutsig))

# Load and preprocess signature data
sig=as.data.frame(fread('Path to SigProfiler output',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA"))) # Replace with actual path
sig1=sig[,2:87] # Select relevant columns
rownames(sig1)=sig$Samples # Set row names to Samples
# Reshape the dataframe
sig2=sig1[,colSums(sig1)!=0] # Keep columns with non-zero sums
sig2=sig2/rowSums(sig2) # Normalize columns
sig2$id=rownames(sig2) # Add row names as identifiers
sig1_long <- melt(sig2) # Melt the dataframe
# Merge melted signature data with mutational signature data
sig1_long=merge(sig1_long,mutsig,by.x='variable',by.y='sbs',all.x=TRUE) # Merge on variable and sbs columns

# Filter out rows where id contains 'no'
sig1_long1=sig1_long[str_detect(sig1_long$id, 'no'),] # Exclude rows with 'no' in id

# Assign rupture timing based on id pattern
sig1_long1$rupture <- ifelse(grepl("EarlyRupture", sig1_long1$id), "Early Rupture",
                             ifelse(grepl("MidRupture", sig1_long1$id), "Mid Rupture",
                                    ifelse(grepl("LateRupture", sig1_long1$id), "Late Rupture", NA))) # Conditional assignment of rupture timing

# Filter out rows where value is not zero
sig1_long1=sig1_long1[sig1_long1$value!=0,] # Exclude rows with zero value

# Identify unique mutational signatures
aetiology=unique(sig1_long1$mutsig) # Extract unique mutational signatures

#Track signatures that show a significant difference based on the rupture timing of the chromosomes
track=0
for(i in 1:length(aetiology))
{
  test=sig1_long1[sig1_long1$mutsig==aetiology[i],]
  anova_result <- summary(aov(value ~ rupture, data = test))[[1]]$`Pr(>F)`[1]
  if(anova_result<0.05)
  {
    track=c(track,i)
  }
}
# Example plot for showcasing significant differences in mutational signatures based on rupture timing of micronuclei containing certain chromosomes
test=sig1_long1[sig1_long1$mutsig==aetiology[2],] # Filter for the second most common signature

# Perform ANOVA to compare means across rupture timings for the selected signature
anova_result=(aov(value ~ rupture, data = test)) # ANOVA analysis
TukeyHSD(anova_result) # Tukey HSD post-hoc test to identify significant differences

# Plot means across rupture timings for the selected signature
plotmeans(value ~ rupture, data = test) # Plot means

# Subset data for APOBEC mutational signature
apobec=sig1_long[sig1_long$mutsig=='APOBEC',] # Filter for APOBEC signature

# Extract rupture timing and sample information from id
apobec$rupture=str_extract(apobec$id, "(?<=_).*?(?=_|$)")   # Extract rupture timing
apobec$sample=str_extract(apobec$id, "[^_]+") # Extract sample name
apobec$status=str_extract(apobec$id, "_(.*)") # Extract status

# Filter for rows where id contains 'chromothripsis'
apobec1=apobec[str_detect(apobec$id, 'chromothripsis'),] # Include rows with 'chromothripsis' in id

# Aggregate mean value by status
averaged_data <- aggregate(value ~ status, apobec1, mean) # Calculate mean values

# Bar plot of aggregated mean values by status
barplot(averaged_data$value,names.arg = averaged_data$status) # Bar plot

# Filter for rows where rupture is not 'MidRupture'
apobec3=apobec1[apobec1$rupture!='MidRupture',] # Exclude rows with 'MidRupture'

# Filter for rows where id contains 'no'
apobec2=apobec[str_detect(apobec$id, 'no'),] # Include rows with 'no' in id

# Aggregate mean value by rupture for filtered data
averaged_data <- aggregate(value ~ rupture, apobec3, mean) # Calculate mean values

# T-test comparing values across rupture timings for filtered data
t.test(value~rupture,data = apobec3) # T-test

# ANOVA to compare means across rupture timings for all data
res_aov <- aov(value~rupture,data = apobec3) # ANOVA analysis
summary(res_aov) # Summary of ANOVA results

# Repeat ANOVA for original data excluding 'no' rows
res_aov <- aov(value~rupture,data = apobec1) # ANOVA analysis
summary(res_aov) # Summary of ANOVA results

# Boxplot and jitter plot of APOBEC signature analysis
ggplot(apobec1, aes(x=status, y=value)) +
  geom_boxplot() + # Boxplot
  scale_fill_viridis(discrete = TRUE, alpha=0.6) + # Color scale
  geom_jitter(color="black", size=0.4, alpha=0.9) + # Jitter plot
  theme(plot.title = element_text(size=11)) + # Theme settings
  ggtitle("APOBEC Signature Analysis") + # Title
  xlab("Rupture Timing") # X-axis label
