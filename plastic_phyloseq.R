
#-----MAKING A PHYLOSEQ OBJECT-----#

#load data into R for phyloseq analysis

#load sample data
sdata <- as.csv("/home/lgschaer/old/Plastic_Deg/Plastic_16S/PlasticDegSampleIDs.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

#Here I am just expanding my metadata from one column into all of the columns I need for my metadata
#You can skip this!
sdata2 <- sdata %>% 
  mutate(SampleInfo = SampleID) %>%
  separate(SampleInfo, into = c("Inocula", "Substrate", "Replicate"), fill = "right") %>%
  mutate(Replicate = ifelse(Substrate=="R1"|Substrate=="R2"|Substrate=="R3"|Substrate=="R4", Substrate, Replicate),
         Substrate = ifelse(Substrate!="TPA"&Substrate!="BPA"&Substrate!="TA", "Inocula", Substrate),
         Replicate = ifelse(is.na(Replicate), "Blank", Replicate),
         Substrate = ifelse(is.na(Substrate), "Blank", Substrate))
View(sdata2)

write.csv(sdata2, "/home/lgschaer/old/Plastic_Deg/Plastic_16S/plastic_meta.csv")
#Stop skipping here.

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Plastic_Deg/Plastic_16S/dada2_output/seqtab.rds")
colnames(sequence_table) <- NULL                                   #remove column names
sequence_table <- as.data.frame(sequence_table)                    #convert to data frame format
sequence_table[1:4,1:4]                                            #view a portion to make sure it looks how we expect
class(sequence_table)                                              #check class, should be data frame 

#Here I am fixing the sample names 
#You can skip this too
sequence_table <- rownames_to_column(sequence_table, var = "SampleID")    #change the rownames to a column with header "Sample_ID"
sequence_table[1:4,1:4]                                            #view again
class(sequence_table)  
sampleNames <- sequence_table$SampleID
sampleNames
sequence_table <- column_to_rownames(sequence_table, var = "SampleID") 

sdata3 <- sampleNames %>%
  as_tibble() %>%
  mutate(rownames = value) %>%
  separate(value, into = c("SampleID", "SampleNumber"), sep = "_S") %>%
  left_join(sdata2, by = "SampleID") %>%
  mutate(SampleID = rownames) %>%
  unite(Inocula_Substrate, Inocula, Substrate, sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "rownames")
head(sdata3)
#Stop skipping here

#make nonzero subset to remove all columns with zero taxa counts
subset_all <- as.matrix(sequence_table)                                #change to matrix format
m <- (colSums(subset_all, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- subset_all[, m]                                         #all the non-zero columns


#load taxa table
taxa_table <- readRDS("/home/lgschaer/old/Plastic_Deg/Plastic_16S/dada2_output/taxa.rds")
taxa_table <- as.matrix(taxa_table)                                #change to matrix format
taxa_table[1:5,1:5]                                                #view to make sure everything looks good


#make phyloseq object
samdata = sample_data(sdata3)                                      #define sample data
sample_names(samdata)
colnames(nonzero) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(taxa_table)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table

taxa_names(taxtab)
seqtab[1:4,1:4]
taxa_names(seqtab)

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all

#sequencing depth before rarefying
TPA <- subset_samples(phyloseq_object_all, Substrate == "TPA")
median(sample_sums(TPA))
TA <- subset_samples(phyloseq_object_all, Substrate == "TA")
median(sample_sums(TA))
BPA <- subset_samples(phyloseq_object_all, Substrate == "BPA")
median(sample_sums(BPA))


#sample counts before rarefying
sample_counts <- sample_data(phyloseq_object_all) %>%
  group_by(Substrate) %>%
  mutate(Count = 1) %>%
  summarise(SumCount = sum(Count))
sample_counts

#remove blanks
phyloseq_object_all_nob <- phyloseq_object_all %>% subset_samples(Substrate != "Blank")

#Look at metadata variables
head(sample_data(phyloseq_object_all_nob))

#group samples together
merged_phyloseq = merge_samples(phyloseq_object_all_nob, group = c("Inocula_Substrate"), fun = sum)
merged_phyloseq

#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(merged_phyloseq, sample_sums(merged_phyloseq) > 500)

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

#noblanks_notmerged_rarefied <- rarefy_samplesover1000_all

pps <- rarefy_samplesover1000_all
pps

#NOTE: PAY ATTENTION TO RAREFYING STEP BEFORE RUNNING BELOW, I DID MANY ITERATIONS CHANGING WHETHER BELOW WERE MERGED OR NOT AND RAREFIED OR NOT
#assign colors to samples
colors <- c("Orange", "Purple", "Green", "Magenta", "Black", "Blue", "Red", "White")

#filter out eukaryotes and mitochondria
justbacteria_notmerged <- phyloseq_object_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) 
justbacteria_notmerged

#filter out eukaryotes and mitochondria
justbacteria_notrare_merged <- merged_phyloseq %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) 
justbacteria_notrare_merged

#filter out eukaryotes and mitochondria
justbacteria_rare_notmerged <- noblanks_notmerged_rarefied %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(Substrate != "BPA")
justbacteria_rare_notmerged

#filter out eukaryotes and mitochondria
justbacteria <- merged_phyloseq %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) 
justbacteria

#filter out eukaryotes and mitochondria
justbacteria_rare_merged <- merged_phyloseq %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) 
justbacteria_rare_merged

#saveRDS file
#saveRDS(justbacteria, file = "/home/lgschaer/old/Plastic_Deg/Plastic_16S/plastic_phyloseq_merged_notrarified.rds")

#-----ALPHA DIVERISTY-----#

#Alpha diversity, just observed, all substrates and inocula grouped together
justbacteria_rare_notmerged %>%                                      #phyloseq object
  plot_richness(
    x = "Substrate",                                                 #compare diversity of datatype
    measures = c("Observed")) +                                      #choose diversity measures
  geom_boxplot(aes(fill = Substrate), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_point() +                                                     #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                                #set fill colors
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#Subset to only inlude Inocula samples
justsediment <- justbacteria_rare_notmerged %>%
  subset_samples(Substrate == "Inocula")
justsediment

#Alpha diversity, just observed, all substrates and inocula grouped together
justsediment %>%                                      #phyloseq object
  plot_richness(
    x = "Inocula",                                                 #compare diversity of datatype
    measures = c("Observed")) +                                      #choose diversity measures
  geom_boxplot(aes(fill = Substrate), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_point() +                                                     #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                                #set fill colors
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position



#-----BETA DIVERSITY-----#


#t-SNE plot

tsne <- tsne_phyloseq(justbacteria_rare_notmerged, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())
summary(tsne)

justbacteria

#tSNE Plot
plot_tsne_phyloseq(justbacteria_rare_notmerged, tsne, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate, shape = Inocula), size = 5) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8, 2))+
  scale_fill_manual(values = colors)+
                    #breaks = sample_types,
#                    labels = sample_labels) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

warnings()

# PERMANOVA
all_pcoa <- ordinate(
  physeq = justbacteria_rare_notmerged, 
  method = "PCoA", 
  distance = "bray"
)

#install.packages("colorspace")
library(colorspace)

sed_colors <- terrain.colors(7, alpha = 1)
sed_colors <- sequential_hcl(7)
sed_colors <- terrain_hcl(7, h = c(600, 40), c = 96, l = c(20, 90))
sed_colors

#plot
plot_ordination(
  physeq = justbacteria_rare_notmerged,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(shape = Substrate, fill = Inocula), size = 10) +                         #sets fill color to sampletype
  scale_fill_manual(values = sed_colors) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command


#--------------------------------------------------------#
# MAKING TAXA PLOTS

#Summarize abundance of each class
classabundance <- justbacteria_notrare_merged %>%
  tax_glom(taxrank = "Class") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class) 
head(classabundance)

#Select and summarize necessary variables
classOfInocula <- classabundance %>%
  select(Phylum, Class, Sample, Abundance) %>%
  separate(Sample, into = c("Inocula", "Substrate")) %>%
  filter(Substrate == "Inocula") %>%
  mutate(
    Phylum = as.character(Phylum),                                #Change Phylum and Class cols to character vectors
    Class = as.character(Class),
    Taxa = ifelse(Phylum == "Proteobacteria", Class, Phylum),     #Create "Taxa" column which shows Class for Proteobacteria, Phylum for all other phyla
    Taxa.1p = ifelse(Abundance < 0.01, "<1%", Taxa),              #Label taxa present at low abundance "< 1%" for both Taxa and Class
    Class.1p = ifelse(Abundance < 0.01, "<1%", Class)
  ) 
head(classOfInocula)


#MAKING A TAXA PLOT BY CLASS

#make color vector

colors10 <- c(
  "black",   "orchid1",      "red",       "green",       "coral1",           "cyan",
  "grey47",  "yellow",       "blue",      "darkgreen",   "palegoldenrod",    "darkblue",
  "grey77",  "orange",       "darkcyan",  "magenta",     "mediumpurple1",    "tan4",
  "white",   "yellowgreen",  "firebrick", "dodgerblue",  "purple4", "lawngreen", "lightpink"
)  


justInocula <- classOfInocula %>%
  filter(Substrate != "Blank") %>%
  group_by(Inocula, Phylum, Class.1p) %>%
  summarise(Abundance = sum(Abundance))
head(justInocula)

#just the inocula
ggplot(justInocula)+
  geom_col(mapping = aes(x = Inocula, y = Abundance, fill = Class.1p), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  #scale_x_discrete(
  #  breaks = b,
   # labels = la,
    #limits = li)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

#MAKING A TAXA PLOT OF CULTURES SHOWING FAMILY

#Summarize abundance of each family
familyabundance <- justbacteria_notrare_merged %>%
  tax_glom(taxrank = "Family") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family) 
head(familyabundance)

#Select and summarize necessary variables
familyOfEnrich <- familyabundance %>%
  select(Phylum, Class, Order, Family, Sample, Abundance) %>%
  separate(Sample, into = c("Inocula", "Substrate")) %>%
  filter(Substrate != "Inocula") %>%
  mutate(
    Phylum = as.character(Phylum),                                #Change Phylum and Class cols to character vectors
    Class = as.character(Class),
    Class.1p = ifelse(Abundance < 0.01, "<1%", Class),
    Order = as.character(Order),
    Order.1p = ifelse(Abundance < 0.01, "<1%", Order),
    Family = as.character(Family),
    Family.1p = ifelse(Abundance < 0.01, "<1%", Family)
  ) 
head(familyOfEnrich)

substrateFamily <- familyOfEnrich %>%
  filter(Substrate != "Blank" & Substrate != "BPA") %>%
  group_by(Substrate, Phylum, Family.1p) %>%
  summarise(Abundance = sum(Abundance))
head(substrateFamily)

substrateClass <- familyOfEnrich %>%
  filter(Substrate != "Blank" & Substrate != "BPA") %>%
  group_by(Substrate, Phylum, Class.1p) %>%
  summarise(Abundance = sum(Abundance))
head(substrateClass)


#just the enrichments at class level
ggplot(substrateClass)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Class.1p), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = c("black", "red", "green", "darkblue", "white")) +
  #scale_x_discrete(
  #  breaks = b,
  # labels = la,
  #limits = li)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

#just the substrates at family level
ggplot(substrateFamily)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Family.1p), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))



#save CSV of abundance of taxa throughout ww treatment process
head(all2)

csv <- all2 %>%
  select(SampleID, Phylum, Class, Abundance) %>%
  filter(Abundance > 0.01)%>%
  arrange(desc(Abundance))
head(csv)

#write_csv(csv, "/home/lgschaer/old/wastewater_analysis/taxa_abundance.csv")

csv2 <- all2 %>%
  select(SampleID, Phylum, Abundance) %>%
  filter(Abundance > 0.01) %>%
  group_by(SampleID, Phylum) %>%
  summarize(
    SumAb = sum(Abundance)
  )
head(csv2)

View(csv2)

#write_csv(csv2, "/home/lgschaer/old/wastewater_analysis/taxa_abundance2.csv")

### Other Figures


