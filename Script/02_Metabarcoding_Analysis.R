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
library(igraph)
library(effect.lndscp) #devtools::install_github("pedroj/effectiveness_pckg")

#Import reads
metabarcoding<- read.csv("Data/Metabarcoding_Aride.csv")

#Rename tortoise IDz
metabarcoding$Tortoise_ID <- ifelse(metabarcoding$Tortoise_ID =="T1", "T01",
                             ifelse(metabarcoding$Tortoise_ID =="T2", "T02",
                             ifelse(metabarcoding$Tortoise_ID =="T3", "T03",
                             ifelse(metabarcoding$Tortoise_ID =="T5", "T05",
                             ifelse(metabarcoding$Tortoise_ID =="T6", "T06",
                             ifelse(metabarcoding$Tortoise_ID =="T7", "T07",
                             ifelse(metabarcoding$Tortoise_ID =="T8", "T08",
                                    metabarcoding$Tortoise_ID)))))))

metabarcoding <- metabarcoding %>% select(!'Unassigned.__.__.__.__.__.__' )%>% 
                                   select(!'index' ) %>% select(!GPS_attached)

metabarcoding <- metabarcoding %>% group_by(Tortoise_ID, Fecal_sample_number) %>%
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

#Filter number of reads to mitigate potential error during sequencing (Reads>25)
metabarcoding <- metabarcoding %>% filter(Reads>25)

#Calculate mean of reads across replicates of fecal samples
metabarcoding <- metabarcoding %>%group_by(Tortoise_ID,Fecal_sample_number, Level_ID,Plant_ID) %>%
                    summarise(Reads=mean(Reads))

#Categorize reads of leftovers
leftovers <- c("Aegilops", "Averrhoa_carambola", "Cannabis_sativa", "Capsicum_frutescens",
              "Coffea","Cucumis_melo","Cucumis_sativus", "Juglans_regia", "Lactuca_serriola",
              "Musa_x_paradisiaca","Prunus_armeniaca","Prunus_persica","Prunus",
              "Sesamum_indicum","Sinapis_alba","Triticum", "Anacardium_occidentale",
              "Solanum", "Pistacia_vera")

metabarcoding$Leftovers <- ifelse(metabarcoding$Plant_ID %in% leftovers, 1, 0)
rm(leftovers)

metabarcoding %>% group_by(Leftovers) %>% summarise("Count"=sum(Reads)) 
#1.5% of reads correspond to leftovers
#Filter them

summary <-metabarcoding %>%filter(Leftovers==0) %>%  group_by(Level_ID) %>% 
  summarise("Total_Reads" = sum(Reads),
            "Unique_Reads"= length(unique(Plant_ID)))

summary$Perc_Reads <- round((summary$Total_Reads / sum(metabarcoding$Reads))*100,2)
writexl::write_xlsx(summary, "Results/Metabacoding_Reads_Summary.xlsx")
rm(summary)#Clean environment

### Correct metabarcoding ids ######

metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Arenaria_longipedunculata","Arenaria",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Asystasia_nemorum","Asystasia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Boerhavia_diffusa","Boerhavia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Cocculus_orbiculatus","Cocculus",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Euphorbia_boivinii","Euphorbia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Euphorbia_falcata","Euphorbia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Ficus","Ficus_reflexa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Ficus_natalensis","Ficus_reflexa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Ficus_saussureana","Ficus_reflexa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Gossypium_anomalum","Gossypium",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Lavandula_multifida","Lavandula",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Mastixia_chinensis","Mastixia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Molinia_caerulea","Molinia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Morella_faya","Morella",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Musa_acuminata","Musa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Musa_balbisiana","Musa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Musa_ornata","Musa",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Pavonia_immitis","Pavonia",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Pericallis_murrayi","Pericallis",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Pisonia","Pisonia_grandis",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Pisonia_subcordata","Pisonia_grandis",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Platanus_orientalis","Platanus",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Psidium_cattleyanum","Psidium",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Taraxacum_officinale","Taraxacum",metabarcoding$Plant_ID)
metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="Urtica_dioica","Urtica",metabarcoding$Plant_ID)
#metabarcoding$Plant_ID <- ifelse(metabarcoding$Plant_ID=="","",metabarcoding$Plant_ID)

# Remove sp that are likely misidentifications and not present on island

metabarcoding <- metabarcoding %>% filter(Plant_ID!="Morella")%>% 
                                   filter(Plant_ID!="Molinia")%>%
                                   filter(Plant_ID!="Kleinia_neriifolia")


metabarcoding <- metabarcoding %>% filter(Leftovers==0) %>% 
                group_by(Tortoise_ID, Fecal_sample_number, Level_ID,Plant_ID) %>%
                summarise("Reads"=mean(Reads))

str(metabarcoding)

samples <- metabarcoding %>% group_by(Tortoise_ID) %>%
  summarise(unique_fecal_samples = n_distinct(Fecal_sample_number))

metabarcoding <- left_join(metabarcoding, samples, by="Tortoise_ID")

metabarcoding <-metabarcoding %>% group_by(Tortoise_ID, Plant_ID) %>% summarise(Reads = mean(Reads),
                                                     Samples=unique(unique_fecal_samples),
                                                     Level_ID=unique(Level_ID))
rm(samples)

# Calculate total number of reads per tortoise ID
total_reads <- metabarcoding %>% group_by(Tortoise_ID) %>%summarise(Total_Reads = sum(Reads))

#Calculate relative frequency of reads for each plant species for each tortoise
metabarcoding <- metabarcoding %>% left_join(total_reads, by = "Tortoise_ID")

#Relative read abundance: Formula from Deagle et al 2018 Molecular Ecology
metabarcoding$RRA <- round(1/metabarcoding$Samples*(metabarcoding$Reads/metabarcoding$Total_Reads),5)

#Filter only identifications at genus and species level
metabarcoding <- metabarcoding %>% filter(Level_ID=="Genus" | Level_ID=="Species")

#import sampling effort to standarize RRA
effort <- read.csv("Data/Sampling_effort_Tortoise.csv")
metabarcoding <- left_join(metabarcoding,effort,by="Tortoise_ID")

#standarize by 100 hours of sampling
metabarcoding$RRA_ST <- (metabarcoding$RRA /metabarcoding$Effort)
metabarcoding$Reads_ST <- round((1/metabarcoding$Samples *metabarcoding$Reads) /metabarcoding$Effort) 


write.csv(metabarcoding, "Data/CLEAN/Metabarcoding_Browsing.csv",row.names = FALSE)
#rm(list=ls())


#### Browsing effectiveness #####
native <- c( "Pisonia_grandis",  "Morinda_citrifolia",  "Commelina_benghalensis", 
  "Cordia_subcordata",  "Cyperus_rotundus", "Ficus_reflexa",  "Glinus_oppositifolius", 
  "Ochrosia_oppositifolia", "Boerhavia",  "Boerhavia_diffusa", "Brachiaria_ramosa", 
  "Cyperus", "Cyperus_polystachyos", "Digitaria_ciliaris",  "Pisonia", "Acalypha_lanceolata", 
  "Portulaca_oleracea",  "Ipomoea_pes_caprae",  "Scaevola_taccada")

browsing_efectiveness<-metabarcoding %>%
  group_by(Tortoise_ID) %>%
  filter(!Plant_ID %in% c("Ficus_reflexa", "Pisonia_grandis", "Artocarpus_altilis", 
                  "Morinda_citrifolia", "Terminalia_catappa", "Euphorbia_pyrifolia")) %>%
  summarise(
    Samples = unique(Samples),
    Total_Reads = sum(Reads),
    Native_Reads = sum(Reads[Plant_ID %in% native]),
    Prop_Native_Reads = (Native_Reads /Total_Reads)*100
  )



browsing_efectiveness <- left_join(browsing_efectiveness, effort, by="Tortoise_ID")

browsing_efectiveness$Reads_Sample <- (browsing_efectiveness$Total_Reads/ browsing_efectiveness$Samples)

browsing_efectiveness <- browsing_efectiveness %>%
  add_row(Tortoise_ID = "T11", Prop_Native_Reads = 0, Reads_Sample =0)%>%
  add_row(Tortoise_ID = "TXX", Prop_Native_Reads = 0, Reads_Sample =316200)%>% # plot isolines to 100% of x axis
  add_row(Tortoise_ID = "TX1", Prop_Native_Reads = 100, Reads_Sample =316200) # plot isolines to 100% of x axis

brow_eff<-effectiveness_plot(browsing_efectiveness$Prop_Native_Reads, 
                   log10(browsing_efectiveness$Reads_Sample),
                   label = browsing_efectiveness$Tortoise_ID,  
                   myxlab = "Proportion of native plants (%)", 
                   myylab = "log10(No. reads / sample)") + 
                   scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits=c(0, 5.5), breaks=c(0, 1, 2, 3, 4, 5))

write.csv(browsing_efectiveness[,c("Tortoise_ID","Effort","Samples")],"Data/CLEAN/Sampling_effort_Tortoise.csv",row.names = FALSE)


detritivory_efectiveness<-metabarcoding %>%
  group_by(Tortoise_ID) %>%
  filter(Plant_ID %in% c("Ficus_reflexa", "Pisonia_grandis", "Artocarpus_altilis", 
                          "Morinda_citrifolia", "Terminalia_catappa", "Euphorbia_pyrifolia")) %>%
  summarise(
    Samples = unique(Samples),
    Total_Reads = sum(Reads),
    Native_Reads = sum(Reads[Plant_ID %in% native]),
    Prop_Native_Reads = (Native_Reads /Total_Reads)*100
  )

detritivory_efectiveness <- left_join(detritivory_efectiveness, effort, by="Tortoise_ID")

detritivory_efectiveness$Reads_Sample <- (detritivory_efectiveness$Total_Reads/ detritivory_efectiveness$Samples)

detritivory_efectiveness <- detritivory_efectiveness %>%
  add_row(Tortoise_ID = "TXX", Prop_Native_Reads = 0, Reads_Sample =316200)%>% # plot isolines to 100% of x axis
  add_row(Tortoise_ID = "TX1", Prop_Native_Reads = 100, Reads_Sample =316200) # plot isolines to 100% of x axis

det_eff<-effectiveness_plot(detritivory_efectiveness$Prop_Native_Reads, 
                   log10(detritivory_efectiveness$Reads_Sample),
                   label = detritivory_efectiveness$Tortoise_ID,  
                   myxlab = "Proportion of native plants (%)", 
                   myylab = "log10(No. reads / sample)") + 
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits=c(0, 5.5), breaks=c(0, 1, 2, 3, 4, 5))

rm(list = setdiff(ls(), c("brow_eff", "det_eff")))