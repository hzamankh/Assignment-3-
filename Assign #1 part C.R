
#My question: Comparison between distribution of Rodents in the Northern and Southern parts of Canada
#I want to compare geographic distribution of Rodents in different regions of Canada and also between Northern and Southern parts based on the latitude that the sample was taken. Sorry if it is sound ridiculous.

#First Adding packages needed

library(tidyverse)
library(vegan)
library(iNEXT)
library(reshape2)

#Downloading desired data directly from BOLD website

Rodentia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Rodentia&geo=Canada&format=tsv")

#Saving the data to a file for further use

write_tsv(Rodentia, "Rodentia_BOLD_data.tsv")

#There are many variables exist in this dataframe, to check them:

names(Rodentia)

#Most of these data are not required here, so choose some of them to make data smaller

Rodentia2 <- Rodentia [c(1, 8, 10, 12, 14, 16, 18, 20, 22, 24, 40, 47, 51, 55, 56, 57, 70, 71, 72)]

#This sub data still large but better than the last on. To check what we have here:

str(Rodentia2)

names(Rodentia2)

#Some columns are still not useful such as phylum_name or class_name and order_name because all we have here are RODENTS, just keep them to check if something absurd exist in the data. Some other columns like subspecies_name and elev have many NA valuable and are not useful. So, the data will be cut to a smaller version again

Rodentia3 <- Rodentia2 [c(1, 2, 5, 6, 7, 8, 9, 12, 15, 16, 17, 18, 19)]

str(Rodentia3)

names(Rodentia3)

#Now just draw some graphs of data to give us a visual summary

hist(Rodentia3$lat)

#Because the Province data where Rodents found is not numeric, it is not possible to draw a histogram of it with labels, so ggplot is better option

png("Regionplot.png", units="px", width=6000, height=4800, res=300)
ggplot(data.frame(Rodentia3$province_state), aes(x=Rodentia3$province_state)) + geom_bar()
dev.off()

#Now check the maximum and minimum of latitidue

which.max(Rodentia3$lat)
Rodentia3$lat [165]

which.min(Rodentia3$lat)
Rodentia3$lat [631]

#Checking which provinces has the max and min number of samples

which.max(table(Rodentia3$province_state))

which.min(table(Rodentia3$province_state))

#Actually another province also has frequency of 1 which this command dose not show!!!!! for now I don't want to find solution to this

#The last two line show us which states have the max and min numbers but for further use and also ease of use it is better to save these data

Province_Rodentia3 <- as.data.frame(table(Rodentia3$province_state))
names(Province_Rodentia3) [1] <- paste("Province")
names(Province_Rodentia3) [2] <- paste("Frequency")

#Now to show the numbers of samples from max and min states:

Province_Rodentia3$Frequency [8]
Province_Rodentia3$Frequency [10]

#Based on lat max and min and some geographic search, 60 degree lat is considered the criteria for deviding Canada to Northern and Southern parts

#First need to get rid of data without lat for the following analysis

Rodentia4 <- drop_na(Rodentia3, lat)

#Devide data into two groups based on lat

Northern <- filter(Rodentia4, Rodentia4$lat > 60)
Southern <- filter(Rodentia4, Rodentia4$lat < 60)
length(Northern)
    
chisq.test(Rodentia4$family_name, Rodentia4$lat)

#Filter data to only contain genus name and province they found for vegdist analysis

Rodentia5 <- Rodentia3 [c(6, 9)]

#Get rid of NA elements 

Rodentia6 <- drop_na(Rodentia5, genus_name)
Rodentia6 <- drop_na(Rodentia6, province_state)

#For using vegdist we need to group data by genus and province name in a way that we have a matrix with rows and columns are these two. The first step is grouppi 

Rodentia6.grouped <- Rodentia6 %>%
  group_by(genus_name, province_state) %>%
  summarize(n=length(genus_name)) %>%
  print()

#Next step is to move province name to the columns and filling the table which took lots of time to find the answer and solve it.

temp <- dcast(Rodentia6.grouped, genus_name ~ province_state)

#for analysing NA value should be replace with 0 (I am not sure but the tutorial I found online the file contain only numbers so I did the same)

Rodentia_final <- replace(temp, is.na(temp), 0)

#Using vegdist

Result <- vegdist(data.matrix(Rodentia_final))
hist(Result)

#I don't know much about the interpretation of the results

