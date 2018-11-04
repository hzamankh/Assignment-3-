#I want to evaluate phylogeny of tigers to see how they are related to each other if possible relate this to the geographic regions they live. Tigers (Panthera tigris) is the largest cat species and are classified in the genus Panthera with the lion, leopard and other big cats. Results of genetic analysis indicate that about 2.88 million years ago, the tiger and the snow leopard lineages diverged from the other Panthera species, and that both may be more closely related to each other than to the lion, leopard and jaguar (Wikipedia).Actually I searched and found these after finishing the R coding when I wanted to add comments and because I wasn't satisfied with the results I have got.

#All of the tigers are from the same species and they are only different in subspecies. That is why I want to compare them based on mitocondrial DNA, which is usually used for such analysis. Cytochrome c oxidase I is the main subunit of the cytochrome c oxidase complex and I choose this becuase found the most data available on tigers on based on these gene



#First step is installing packages needed

#install.packages("tidyverse")
library(tidyverse)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
#Downloading desired data directly from BOLD website. Panthera is the genus for the big cats.
#isntall.packages("msa")
library(msa)
#install.packages("ape")
library(ape)
#biocLite("DECIPHER")
library(DECIPHER)

Tigers <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Panthera&format=tsv")

#Saving the data to a file for further use

write_tsv(Tigers, "Tigers_BOLD_data.tsv")


#check and choose the desired variables exist from dataframe. Tried to keep geological data but it was not available for many of the samples

names(Tigers)
Tigers2 <- Tigers [c(1, 22, 24, 47, 48, 55, 70, 72)]

#Check the markers available for tigers to choose for phylogeny. COI-5P is the most common one for the Panthra tigris. There are other markers but I choose this one.

grouped_Tigers <- Tigers2 %>%
  group_by(markercode, species_name) %>%
  summarize(n=length(species_name)) %>%
  arrange(desc(n)) %>%
  print()


#Select the Panthera tigris from other big cats in the dataframe and with the markercode COI-5P

Tigers_selected <- subset(Tigers2, species_name == "Panthera tigris" & markercode == "COI-5P")


#Unfortunately many of the samples dose not have subspecies name. So, I need the delete them from the data

Tigers_final <- drop_na(Tigers_selected, subspecies_name)

#make the final dataframe consisting name of subspecies and nucleotide sequences
Tigers_final2 <- Tigers_final [c(3, 8)]


#change the dataframe to the DNAStringSet for the further analysis

Tigers_final_StringSet <- DNAStringSet(Tigers_final2$nucleotides)


#Alignment using Muscle which is algorithm with good accuracy

Tigers_Alignment <- DNAStringSet(msa(Tigers_final_StringSet, "Muscle"))

print(Tigers_Alignment, show= "complete")


#Adding subspecies name to the alignment results

names(Tigers_Alignment) <- Tigers_final2$subspecies_name

print(Tigers_Alignment, show= "complete")


#Change the alignment data to DNAbin to use for distance matrix

dnaBin_Tigers_Alignment <- as.DNAbin(Tigers_Alignment)


#Making the distance matrix. I used N model which is the raw and shows the number of nucleotide differ between each two of sequences. Also tried JC69 but only saw small differences which didn't have significant effect on the final result and interpretation

distanceMatrix_Tigers <- dist.dna(dnaBin_Tigers_Alignment, model = "N", as.matrix = TRUE, 
                                   pairwise.deletion = TRUE)


#Drawing phylogeny using UPGMA which is based on average-linkage clustering. Also tried single method but the results became more unclear

Tigers_clusters <- IdClusters(distanceMatrix_Tigers,
                                             method = "UPGMA",
                                             cutoff= 0.03,
                                             showPlot = TRUE,
                                             type = "both",
                                             verbose = TRUE)


#Save the phylogenic tree to a tiff file. Sorry because the labels are too long and dose not completely show. 

tiff("Tigers clusters.tiff", width = 2000, height = 3000, units = "px", res = 300)
IdClusters(distanceMatrix_Tigers,
           method = "UPGMA",
           cutoff= 0.03,
           showPlot = TRUE,
           type = "dendrogram",
           verbose = TRUE)
dev.off()


#Trying Neighbor Joining phylogeny but the result was not satisfactory. Also tried other packages like phangorn but could not get beautiful tree, so because of lack of the time I could not use more time for this part. I think it is good if you tried some features of making plots and phylogenies in the class

NJ_Tigers <- nj(distanceMatrix_Tigers)

plot(NJ_Tigers, 'u')


tiff("NJ Tigers.tiff", width = 2000, height = 3000, units = "px", res = 300)
plot(NJ_Tigers, 'u')
dev.off()


#About the interpretation of the results: First, I am not satisfied with the results. I thought that different subspecies should be at different branches and same subspecies at the same branches actually with some exceptions. I wish too see which subspecies are closer together. But based on what I have got, the Bengal tiger (P. t. tigris) samples spread through all branches, and one of the is outlier to the whole others. The same story is correct about Siberian tiger (P. t. altaica). The only subspecies that have similarity to each other is Indochinese tiger (P. t. corbetti). The South China tiger (P. t. amoyensis) is the most distance one but actually it is considered the most ancient one and I expected it to be the outlier. I searched and found an article (Phylogeography and Genetic Ancestry of Tigers (Panthera tigris), Shu-Jin Luo et al. 2004) which have complete different results than mine. They used 3 different markers or analysis which are different than mine. There are some reasons might cause my result: Low number of samples, bad analytic methods and algorthms and most important one bad marker for analysis




