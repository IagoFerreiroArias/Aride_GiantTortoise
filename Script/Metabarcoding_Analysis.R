# Scripts: Metabarcoding analysis of Aldabra giant tortoise diet
# Author: Iago Ferreiro-Arias
# Date: 5th May, 2024

library(ggplot2)
library(dplyr)
library(lme4)
library(sjPlot)
library(tidyr)
library(ggpubr)
library(bipartite)
library(tidyr)

#Import reads
metabarcoding<- read.csv("Data/Metabarcoding_Aride.csv")

#Rename tortoise ID

metabarcoding$Tortoise_ID <- ifelse(metabarcoding$Tortoise_ID =="T1", "T01",
                             ifelse(metabarcoding$Tortoise_ID =="T2", "T02",
                             ifelse(metabarcoding$Tortoise_ID =="T3", "T03",
                             ifelse(metabarcoding$Tortoise_ID =="T5", "T05",
                             ifelse(metabarcoding$Tortoise_ID =="T6", "T06",
                             ifelse(metabarcoding$Tortoise_ID =="T7", "T07",
                             ifelse(metabarcoding$Tortoise_ID =="T8", "T08",
                                    metabarcoding$Tortoise_ID)))))))

metabarcoding <- metabarcoding %>% select(!'Unassigned.__.__.__.__.__.__' )%>% 
                                   select(!'index' ) %>% select(!GPS_attached) %>%
                                   select(!Fecal_sample_number)

metabarcoding <- metabarcoding %>% group_by(Tortoise_ID) %>%
  pivot_longer(cols = 1:120, names_to = "Meta_ID", values_to = "Reads")

#Create Level_ID and Plant_ID columns
metabarcoding <- metabarcoding %>%
  mutate(Level_ID = case_when(
    grepl("\\.s__", Meta_ID) ~ "Species",
    grepl("\\.g__", Meta_ID) ~ "Genus",
    grepl("\\.f__", Meta_ID) ~ "Family",
    grepl("\\.o__", Meta_ID) ~ "Order",
    grepl("\\.c__", Meta_ID) ~ "Class",
    grepl("\\.p__", Meta_ID) ~ "Phylum",
    grepl("\\.k__", Meta_ID) ~ "Kingdom",
    TRUE ~ NA_character_
  )) %>%
  mutate(Plant_ID = case_when(
    Level_ID == "Species" ~ gsub("\\.", "_", gsub("^.+\\.s__", "", Meta_ID)),
    Level_ID == "Genus" ~ gsub("\\.__", "", gsub("^.+\\.g__", "", Meta_ID)),
    Level_ID == "Family" ~ gsub("\\.__\\.\\__", "", gsub("^.+\\.f__", "", Meta_ID)),
    Level_ID == "Order" ~ gsub("\\.__\\.\\__\\.\\__", "", gsub("^.+\\.o__", "", Meta_ID)),
    Level_ID == "Class" ~ gsub("\\.__\\.\\__\\.\\__\\.\\__", "", gsub("^.+\\.c__", "", Meta_ID)),
    Level_ID == "Phylum" ~ gsub("\\.__\\.\\__\\.\\__\\.\\__\\.\\__", "", gsub("^.+\\.p__", "", Meta_ID)),
    Level_ID == "Kingdom" ~ gsub("\\.__\\.\\__\\.\\__\\.\\__\\.\\__\\.\\__", "", gsub("^.+\\.k__", "", Meta_ID))
  ))

#Plant_ID = unknown are always Phyllum =Streptophyta & Class = unknown
metabarcoding$Level_ID <- ifelse(metabarcoding$Plant_ID=="unknown", "Phylum", metabarcoding$Level_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="unknown", "Streptophyta", metabarcoding$Plant_ID)

#Filter number of reads to mitigate ppotential error during sequencing (Reads>25)
metabarcoding <- metabarcoding %>% filter(Reads>25)

summary <-metabarcoding %>% group_by(Level_ID) %>% 
                            summarise("Total_Reads" = sum(Reads),
                                      "Unique_Reads"= length(unique(Plant_ID)))

summary$Perc_Reads <- round((summary$Total_Reads / sum(metabarcoding$Reads))*100,2)

#Categorize reads of leftovers
leftovers <- c("Aegilops", "Averrhoa_carambola", "Cannabis_sativa", "Capsicum_frutescens",
              "Coffea","Cucumis_melo","Cucumis_sativus", "Juglans_regia", "Lactuca_serriola",
              "Musa_x_paradisiaca","Prunus_armeniaca","Prunus_persica","Prunus",
              "Sesamum_indicum","Sinapis_alba","Triticum")

metabarcoding$Leftovers <- ifelse(metabarcoding$Plant_ID %in% leftovers, 1, 0)
rm(leftovers)

metabarcoding %>% group_by(Leftovers) %>% summarise("Count"=sum(Reads)) 
#1.51% of reads correspond to leftovers
#Filter them

metabarcoding <- metabarcoding %>% filter(Leftovers==0) 
# Calculate total number of reads per tortoise ID
total_reads <- metabarcoding %>% group_by(Tortoise_ID) %>%summarise(Total_Reads = sum(Reads))

#Calculate relative frequency of reads for each plant species for each tortoise
metabarcoding <- metabarcoding %>% left_join(total_reads, by = "Tortoise_ID") %>%
                                   mutate(Rel_Frequency = round(Reads / Total_Reads, 2))

#Filter only identifications at genus and species level
metabarcoding <- metabarcoding %>% filter(Level_ID=="Genus" | Level_ID=="Species")

#Insert tortoise weight
metabarcoding <- left_join(metabarcoding, tortoises[,1:2], by="Tortoise_ID")