library(dplyr)
library(tidyr)
library(igraph)
library(bipartite)
library(ggplot2)
library(ggpubr)

# Browsing matrix
metabarcoding <- read.csv("Data/CLEAN/Metabarcoding_Browsing.csv")
browsing_focal <- read.csv("Data/CLEAN/Browsing_Focal_Observations.csv")

unique(browsing_focal)
browsing_focal$Plant_ID<- gsub(" ","_", browsing_focal$Plant_Species)
unique(metabarcoding$Plant_ID)

metabarcoding <- metabarcoding %>%
  mutate(Plant_ID = ifelse(!grepl("_", Plant_ID), paste0(Plant_ID, "_sp."), Plant_ID))


metabarcoding <- metabarcoding %>%
  mutate(Plant_ID = ifelse(!grepl("_", Plant_ID), paste0(Plant_ID, "_sp."), Plant_ID))

all(browsing_focal$Plant_ID %in% metabarcoding$Plant_ID)
unique(browsing_focal[-which(browsing_focal$Plant_ID %in% metabarcoding$Plant_ID), "Plant_ID"])

browsing_matrix <- left_join(metabarcoding, browsing_focal, by=c("Tortoise_ID","Plant_ID"))

browsing_focal <- browsing_focal %>% 
                   filter(Plant_ID =="Panicum_brevifolium" |
                          Plant_ID == "Asystasia_gangetica" |
                          Plant_ID == "Kyllinga_polyphylla"|
                          Plant_ID == "Megathyrsus_maximus"|
                          Plant_ID == "Digitaria_horizontalis"|
                          Plant_ID =="Euphorbia_pyrifolia"|
                          Plant_ID =="Amaranthus_dubius")

                   
browsing_focal$RRA_ST <- 0 #Zero reads 
browsing_focal <- browsing_focal %>%
  mutate(Freq_Interaction = rowMeans(select(., RRA_ST, Freq_browsing), na.rm = TRUE)) %>%
  select(Tortoise_ID, Plant_ID, Freq_Interaction)

browsing_matrix <- browsing_matrix %>%
  mutate(Freq_Interaction = rowMeans(select(., RRA_ST, Freq_browsing), na.rm = TRUE)) %>%
  select(Tortoise_ID, Plant_ID, Freq_Interaction)

str(browsing_matrix)

browsing_matrix <- rbind(browsing_matrix, browsing_focal)
str(browsing_matrix)
browsing_matrix$Freq_Interaction <- as.numeric(browsing_matrix$Freq_Interaction)
rm(browsing_focal)

browsing_matrix <- browsing_matrix %>%
  distinct()

browsing_matrix<- browsing_matrix %>% filter(!Plant_ID %in% c("Ficus_reflexa", "Pisonia_grandis", "Artocarpus_altilis", 
                          "Morinda_citrifolia", "Terminalia_catappa", "Euphorbia_pyrifolia"))%>%
  tidyr::pivot_wider(names_from =Tortoise_ID , values_from = Freq_Interaction, values_fill = 0)


browsing_matrix<-as.data.frame(browsing_matrix)
row.names(browsing_matrix) <- browsing_matrix$Plant_ID
browsing_matrix <- browsing_matrix %>% select(-Plant_ID)
#browsing_matrix <- log10(browsing_matrix + 1) # for better visualization of links
row.names(browsing_matrix) <- sub(
  "^([A-Za-z]{4})[a-zA-Z]*_(.{3}).*", # Expresión regular para extraer 4 primeras letras del género y 3 del epíteto
  "\\1_\\2",                          # Reemplaza por las primeras 4 letras del género y las primeras 3 del epíteto
  row.names(browsing_matrix)
)
#rm(list=c("browsing_focal","metabarcoding"))

colnames(browsing_matrix)
browsing_matrix <- browsing_matrix[, c("T01", "T02", "T03", "T05", "T06", "T07", "T08", "T10", "T11", "T12")]

#Plot browsing interaction matrix
bipartite::plotweb(browsing_matrix, text.rot = 90)



detritivory_focal <- read.csv("Data/CLEAN/Detritivory_Focal_Observations.csv")
unique(detritivory_focal)
detritivory_focal$Plant_ID<- gsub(" ","_", detritivory_focal$Plant_Species)
all(detritivory_focal$Plant_ID %in% metabarcoding$Plant_ID)
unique(detritivory_focal[-which(detritivory_focal$Plant_ID %in% metabarcoding$Plant_ID), "Plant_ID"])

metabarcoding_detritivory<- metabarcoding %>% filter(Plant_ID %in% c("Ficus_reflexa", "Pisonia_grandis", "Artocarpus_altilis", 
                                                              "Morinda_citrifolia", "Terminalia_catappa", "Euphorbia_pyrifolia"))



detritivory_matrix <- left_join(metabarcoding_detritivory, detritivory_focal, by=c("Tortoise_ID","Plant_ID"))
detritivory_matrix <- detritivory_matrix %>%
  mutate(Freq_Interaction = rowMeans(select(., RRA_ST, Freq_detritivory), na.rm = TRUE)) %>%
  select(Tortoise_ID, Plant_ID, Freq_Interaction)

detritivory_focal$RRA_ST <- 0 #Zero reads 
detritivory_focal <- detritivory_focal %>% filter(Plant_ID=="Euphorbia_pyrifolia") %>%
  mutate(Freq_Interaction = rowMeans(select(., RRA_ST, Freq_detritivory), na.rm = TRUE)) %>%
  select(Tortoise_ID, Plant_ID, Freq_Interaction)

detritivory_matrix <- rbind(detritivory_matrix, detritivory_focal)
str(detritivory_matrix)
detritivory_matrix$Freq_Interaction <- as.numeric(detritivory_matrix$Freq_Interaction)
rm(detritivory_focal)

detritivory_matrix <- detritivory_matrix %>% distinct()
detritivory_matrix$Freq_Interaction <- as.numeric(detritivory_matrix$Freq_Interaction)

detritivory_matrix<- detritivory_matrix %>%
  tidyr::pivot_wider(names_from =Tortoise_ID , values_from = Freq_Interaction, values_fill = 0)

detritivory_matrix <- as.data.frame(detritivory_matrix)
row.names(detritivory_matrix) <- detritivory_matrix$Plant_ID
detritivory_matrix <- detritivory_matrix %>% select(!Plant_ID)

row.names(detritivory_matrix) <- sub(
  "^([A-Za-z]{4})[a-zA-Z]*_(.{3}).*", # Expresión regular para extraer 4 primeras letras del género y 3 del epíteto
  "\\1_\\2",                          # Reemplaza por las primeras 4 letras del género y las primeras 3 del epíteto
  row.names(detritivory_matrix)
)

#### Create frugivory and seed dispersal matrix ####
frugivory <- read.csv("Data/CLEAN/Frugivory_Focal_Observations.csv")
frugivory <- frugivory[-which(is.na(frugivory$Freq_frugivory)),]
frugivory$Plant_ID <- gsub(" ","_",frugivory$Plant_Species)

seed_dispersal <- read.csv("Data/CLEAN/Seeds_Frequency.csv")

frugivory_matrix <- left_join(seed_dispersal, frugivory, by=c("Tortoise_ID","Plant_ID"))
frugivory_matrix[which(is.na(frugivory_matrix$Freq_frugivory)),"Freq_frugivory"]<-0
frugivory_matrix <- frugivory_matrix %>%
  mutate(Freq_Interaction = rowMeans(select(., Seeds_Frequency, Freq_frugivory), na.rm = TRUE)) %>%
  select(Tortoise_ID, Plant_ID, Freq_Interaction)

frugivory_matrix<- frugivory_matrix %>%
  tidyr::pivot_wider(names_from = Tortoise_ID, values_from = Freq_Interaction, values_fill = 0)
frugivory_matrix<-as.data.frame(frugivory_matrix)
row.names(frugivory_matrix) <- frugivory_matrix$Plant_ID
frugivory_matrix <- frugivory_matrix %>% select(-Plant_ID)
#frugivory_matrix <- log10(frugivory_matrix + 1) # for better visualization of links
row.names(frugivory_matrix) <- sub(
  "^([A-Za-z]{4})[a-zA-Z]*_(.{3}).*", # Expresión regular para extraer 4 primeras letras del género y 3 del epíteto
  "\\1_\\2",                          # Reemplaza por las primeras 4 letras del género y las primeras 3 del epíteto
  row.names(frugivory_matrix)
)
rm(list=c("frugivory","seed_dispersal"))

#Plot browsing interaction matrix


native <-  c("Ficu_ref", "Mori_cit", "Piso_gra", "Euph_pyr", "Comm_ben", 
                      "Cord_sub", "Cype_rot", "Glin_opp", "Ochr_opp", "Boer_sp.", 
                      "Brac_ram", "Cype_sp.", "Cype_pol", "Digi_cil", "Digi_set", 
                      "Port_ole", "Ipom_pes", "Scae_tac", "Mega_max", "Kyll_pol")

colors_plants<-ifelse(row.names(detritivory_matrix) %in% native, "#a2d1a3", "#7F7699")
bipartite::plotweb(detritivory_matrix, method="normal", text.rot= 90, labsize=0.5,
                   col.high = "#2D2D2D", col.low =  colors_plants, 
                   bor.col.high =  "#2D2D2D", bor.col.low = colors_plants,
                   empty=FALSE)

colors_plants<-ifelse(row.names(browsing_matrix) %in% native, "#a2d1a3", "#7F7699")
bipartite::plotweb(browsing_matrix, method="normal", text.rot= 90, labsize=0.5,
                   col.high = "#2D2D2D", col.low = colors_plants,
                   bor.col.high =  "#2D2D2D", bor.col.low = colors_plants,
                   empty=FALSE)

colors_plants<-ifelse(row.names(frugivory_matrix) %in% native, "#a2d1a3", "#7F7699")
bipartite::plotweb(frugivory_matrix, method="normal", text.rot= 90, labsize=0.5,
                   col.high = "#2D2D2D",col.low = colors_plants,
                   bor.col.high =  "#2D2D2D", bor.col.low = colors_plants,
                   empty=FALSE)


#Plot node properties for plant species
plant_metrics_browsing<-as.data.frame(bipartite::specieslevel(browsing_matrix,index=c("degree", "species strength", "closeness"),
                        level="lower"))
plant_metrics_browsing$species <- row.names(plant_metrics_browsing)
plant_metrics_browsing$species.strength <- round(plant_metrics_browsing$species.strength,3)
row.names(plant_metrics_browsing)<-NULL
head(plant_metrics_browsing)

top_10_degree <- plant_metrics_browsing %>%  arrange(desc(degree)) %>%  slice(1:5)
top_10_strength <- plant_metrics_browsing %>%  arrange(desc(species.strength)) %>%  slice(1:5)
top_10_closeness <- plant_metrics_browsing %>%  arrange(desc(weighted.closeness)) %>% slice(1:5)

color_plants <- ifelse(top_10_degree$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_degree$species
plot_degree <- ggplot(top_10_degree, aes(y = reorder(species, desc(degree)), x = degree, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  labs(x = "Degree", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )
print(plot_degree)

color_plants <- ifelse(top_10_strength$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_strength$species
plot_strength <- ggplot(top_10_strength, aes(y = reorder(species, desc(species.strength)), x = species.strength, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 8, by = 2), limits = c(0, 8)) +
  labs(x = "Strength", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

color_plants <- ifelse(top_10_closeness$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_closeness$species
plot_closeness  <- ggplot(top_10_closeness, aes(y = reorder(species, desc(weighted.closeness)), x = weighted.closeness, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(x = "Closeness", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )
browsing_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=3)
print(browsing_ranking)
ggsave(plot=browsing_ranking, filename="Manuscript/Figures_Paper/Base/Browsing_Metrics_Plants.pdf",
       width=9.27,height=6.65, units="cm")
#Plot node properties for frugivory
plant_metrics_frugivory<-as.data.frame(bipartite::specieslevel(frugivory_matrix,index=c("degree", "species strength", "closeness"),
                                                              level="lower"))
plant_metrics_frugivory$species <- row.names(plant_metrics_frugivory)
plant_metrics_frugivory$species.strength <- round(plant_metrics_frugivory$species.strength,3)
row.names(plant_metrics_frugivory)<-NULL
head(plant_metrics_frugivory)

top_10_degree <- plant_metrics_frugivory %>%  arrange(desc(degree)) %>%  slice(1:5)
top_10_strength <- plant_metrics_frugivory %>%  arrange(desc(species.strength)) %>%  slice(1:5)
top_10_closeness <- plant_metrics_frugivory %>%  arrange(desc(weighted.closeness)) %>% slice(1:5)

color_plants <- ifelse(top_10_degree$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_degree$species
plot_degree <- ggplot(top_10_degree, aes(y = reorder(species, desc(degree)), x = degree, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  labs(x = "Degree", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )
print(plot_degree)

color_plants <- ifelse(top_10_strength$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_strength$species
plot_strength <- ggplot(top_10_strength, aes(y = reorder(species, desc(species.strength)), x = species.strength, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 8, by = 2), limits = c(0, 8)) +
  labs(x = "Strength", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

color_plants <- ifelse(top_10_closeness$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_closeness$species
plot_closeness  <- ggplot(top_10_closeness, aes(y = reorder(species, desc(weighted.closeness)), x = weighted.closeness, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(x = "Closeness", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

frugivory_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=3)
ggsave(plot=frugivory_ranking, filename="Manuscript/Figures_Paper/Base/Frugivory_Metrics_Plants.pdf",
       width=9.27,height=6.65, units="cm")


#Plot node properties for detritivory
plant_metrics_detritivory<-as.data.frame(bipartite::specieslevel(detritivory_matrix,index=c("degree", "species strength", "closeness"),
                                                               level="lower"))
plant_metrics_detritivory$species <- row.names(plant_metrics_detritivory)
plant_metrics_detritivory$species.strength <- round(plant_metrics_detritivory$species.strength,3)
row.names(plant_metrics_detritivory)<-NULL
head(plant_metrics_detritivory)

top_10_degree <- plant_metrics_detritivory %>%  arrange(desc(degree)) %>%  slice(1:5)
top_10_strength <- plant_metrics_detritivory %>%  arrange(desc(species.strength)) %>%  slice(1:5)
top_10_closeness <- plant_metrics_detritivory %>%  arrange(desc(weighted.closeness)) %>% slice(1:5)

color_plants <- ifelse(top_10_degree$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_degree$species
plot_degree <- ggplot(top_10_degree, aes(y = reorder(species, desc(degree)), x = degree, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  labs(x = "Degree", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )
print(plot_degree)

color_plants <- ifelse(top_10_strength$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_strength$species
plot_strength <- ggplot(top_10_strength, aes(y = reorder(species, desc(species.strength)), x = species.strength, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 8, by = 2), limits = c(0, 8)) +
  labs(x = "Strength", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

color_plants <- ifelse(top_10_closeness$species %in% native, "#a2d1a3", "#7F7699")
names(color_plants) <- top_10_closeness$species
plot_closeness  <- ggplot(top_10_closeness, aes(y = reorder(species, desc(weighted.closeness)), x = weighted.closeness, fill = species)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(x = "Closeness", y = NULL) + 
  theme_bw() +
  scale_fill_manual(values = color_plants) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),                          
    axis.title.x = element_text(size = 14),                          
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

detritivory_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=3)
ggsave(plot=detritivory_ranking, filename="Manuscript/Figures_Paper/Base/Detritivory_Metrics_Plants.pdf",
       width=9.27,height=6.65, units="cm")

ranking_plants <- ggarrange(detritivory_ranking,browsing_ranking,frugivory_ranking,nrow=3)
print(ranking_plants)
ggsave(plot=ranking_plants, filename = "Figures/Ranking_Plants.pdf",
       width = 15, height = 29, units="cm")


rm(list=c("top_10_degree","top_10_strength","top_10_closeness", "ranking_plants",
          "plot_degree","plot_strength","plot_closeness", "plant_metrics_browsing",
          "plant_metrics_frugivory"))


### Individual rol3 of tortoises
tortoise_frugivory<-as.data.frame(bipartite::specieslevel(frugivory_matrix,
                                                          index=c("normalised degree", 
                                                                  "species strength", 
                                                                  "closeness"),
                                                               level="higher"))
tortoise_frugivory <- tortoise_frugivory %>% rename("degree" = normalised.degree)
tortoise_frugivory$species <- row.names(tortoise_frugivory)
tortoise_frugivory$species.strength <- round(tortoise_frugivory$species.strength,3)
row.names(tortoise_frugivory)<-NULL

plot_degree <- ggplot(tortoise_frugivory, aes(x = reorder(species, degree), y = degree)) +
  geom_bar(stat = "identity", fill = "#f2c772") + coord_flip() +
  scale_y_continuous(limits = c(0, 20)) +
  labs(x = "Tortoise ID", y = "Degree") + theme_bw()

plot_strength <- ggplot(tortoise_frugivory, aes(x = reorder(species, species.strength), y = species.strength)) +
  geom_bar(stat = "identity", fill = "#f2c772") +coord_flip() + scale_y_continuous(limits = c(0, 15)) +
  labs(x = "Tortoise ID", y = "Strength") + theme_bw()

plot_closeness <- ggplot(tortoise_frugivory, aes(x = reorder(species, weighted.closeness), y = weighted.closeness)) +
  geom_bar(stat = "identity", fill = "#f2c772") + coord_flip() +scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Tortoise ID", y = "Closeness") + theme_bw() 


frugivory_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=1)

#Tortoise browsing ranking
tortoise_browsing<-as.data.frame(bipartite::specieslevel(browsing_matrix,
                                                          index=c("normalised degree", 
                                                                  "species strength", 
                                                                  "closeness"),
                                                          level="higher"))
tortoise_browsing <- tortoise_browsing %>% rename("degree" = normalised.degree)
tortoise_browsing$species <- row.names(tortoise_browsing)
tortoise_browsing$species.strength <- round(tortoise_browsing$species.strength,3)
row.names(tortoise_browsing)<-NULL

plot_degree <- ggplot(tortoise_browsing, aes(x = reorder(species, degree), y = degree)) +
  geom_bar(stat = "identity", fill = "#a5b893") + coord_flip() +
  scale_y_continuous(breaks = seq(0, 80, by = 10), limits = c(0, 25)) +
  labs(x = "Tortoise ID", y = "Degree") + theme_bw()

plot_strength <- ggplot(tortoise_browsing, aes(x = reorder(species, species.strength), y = species.strength)) +
  geom_bar(stat = "identity", fill = "#a5b893") +coord_flip() + scale_y_continuous(limits = c(0, 15)) +
  labs(x = "Tortoise ID", y = "Strength") + theme_bw()

plot_closeness <- ggplot(tortoise_browsing, aes(x = reorder(species, weighted.closeness), y = weighted.closeness)) +
  geom_bar(stat = "identity", fill = "#a5b893") + coord_flip() +scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Tortoise ID", y = "Closeness") + theme_bw() 

browsing_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=1)

#Tortoise detritivory ranking
tortoise_detritivory<-as.data.frame(bipartite::specieslevel(detritivory_matrix,
                                                         index=c("normalised degree", 
                                                                 "species strength", 
                                                                 "closeness"),
                                                         level="higher"))
tortoise_detritivory <- tortoise_detritivory %>% rename("degree" = normalised.degree)
tortoise_detritivory$species <- row.names(tortoise_detritivory)
tortoise_detritivory$species.strength <- round(tortoise_detritivory$species.strength,3)
row.names(tortoise_detritivory)<-NULL

plot_degree <- ggplot(tortoise_detritivory, aes(x = reorder(species, degree), y = degree)) +
  geom_bar(stat = "identity", fill = "#bea777") + coord_flip() +
  scale_y_continuous(breaks = seq(0, 80, by = 10), limits = c(0, 25)) +
  labs(x = "Tortoise ID", y = "Degree") + theme_bw()

plot_strength <- ggplot(tortoise_detritivory, aes(x = reorder(species, species.strength), y = species.strength)) +
  geom_bar(stat = "identity", fill = "#bea777") +coord_flip() + scale_y_continuous(limits = c(0, 15)) +
  labs(x = "Tortoise ID", y = "Strength") + theme_bw()

plot_closeness <- ggplot(tortoise_detritivory, aes(x = reorder(species, weighted.closeness), y = weighted.closeness)) +
  geom_bar(stat = "identity", fill = "#bea777") + coord_flip() +scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Tortoise ID", y = "Closeness") + theme_bw() 

detritivory_ranking<-ggarrange(plot_degree, plot_strength, plot_closeness, ncol=1)

tortoise_ranking<-ggarrange(detritivory_ranking, browsing_ranking, frugivory_ranking, ncol=3)
ggsave(plot=tortoise_ranking, filename = "Manuscript/Figures_Paper/Base/Ranking_Tortoises.pdf",
       width = 21.7, height = 29, units="cm")



tortoise_detritivory$Interaction <- "Detritivory"
tortoise_frugivory$Interaction<- "Frugivory and Seed Dispersal"
tortoise_browsing$Interaction <- "Browsing and Grazing"

tortoise_ranking <- rbind(tortoise_detritivory, tortoise_browsing,  tortoise_frugivory)
tortoise_ranking$Interaction <- factor(tortoise_ranking$Interaction, levels=c("Detritivory", "Browsing and Grazing", "Frugivory and Seed Dispersal"))
str(tortoise_ranking)

tortoise_degree<-ggplot(tortoise_ranking, aes(x = Interaction, y = degree, group = species, color = species)) +
  geom_line(size=1.5, alpha=0.5) + geom_point(size=3, alpha=0.5) +  
  labs(x = "Interaction", y = "Degree", color="Tortoise ID") +theme_bw() +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                                  "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                                  "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                                  "T12"="#66B752"))+
  scale_x_discrete(labels = c("Detritivory" = "Detritivory",
                              "Browsing and Grazing" = "Browsing\nand Grazing",
                              "Frugivory and Seed Dispersal" = "Frugivory and \nSeed Dispersal"))

tortoise_strength<-ggplot(tortoise_ranking, aes(x = Interaction, y = species.strength, group = species, color = species)) +
  geom_line(size=1.5, alpha=0.5) + geom_point(size=3, alpha=0.5) +
  labs(x = "Interaction", y = "Strength", color="Tortoise ID") +theme_bw()+
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))+ scale_y_continuous(limits = c(0,15))+
scale_x_discrete(labels = c("Detritivory" = "Detritivory",
                            "Browsing and Grazing" = "Browsing\nand Grazing",
                            "Frugivory and Seed Dispersal" = "Frugivory and \nSeed Dispersal"))

tortoise_closeness<-ggplot(tortoise_ranking, aes(x = Interaction, y = closeness, group = species, color = species)) +
  geom_line(size=1.5, alpha=0.5) + geom_point(size=3, alpha=0.5) +
  labs(x = "Interaction", y = "Closeness", color="Tortoise ID") +
  theme_bw()  + scale_y_continuous(limits=c(0,1))+
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) +
  scale_x_discrete(labels = c("Detritivory" = "Detritivory",
      "Browsing and Grazing" = "Browsing\nand Grazing",
      "Frugivory and Seed Dispersal" = "Frugivory and \nSeed Dispersal"))

ggarrange(tortoise_degree,tortoise_strength, tortoise_closeness, ncol=3,
          common.legend = TRUE)


tortoise_browsing <- tortoise_browsing %>% rename("browsing_degree"=degree,
                                                  "browsing_strength"=species.strength,
                                                  "browsing_closeness"=weighted.closeness)

tortoise_frugivory<- tortoise_frugivory %>% rename("frugivory_degree"=degree,
                                                 "frugivory_strength"=species.strength,
                                                "frugivory_closeness"=weighted.closeness)

tortoise_detritivory<- tortoise_detritivory %>% rename("detritivory_degree"=degree,
                                                   "detritivory_strength"=species.strength,
                                                   "detritivory_closeness"=weighted.closeness)

tortoises_properties <- left_join(tortoise_browsing, tortoise_frugivory, by="species")
tortoises_properties <- left_join(tortoises_properties, tortoise_detritivory, by="species")
cor(tortoises_properties$browsing_degree,tortoises_properties$frugivory_degree)
cor(tortoises_properties$browsing_strength,tortoises_properties$frugivory_strength)
cor(tortoises_properties$browsing_closeness,tortoises_properties$frugivory_closeness)

weight <- read.csv("Data/Tortoise_Traits.csv", sep=";")
weight <- weight %>% rename("species"=Tortoise_ID)
tortoises_properties<- left_join(tortoises_properties, weight, by="species")

library(smplot2)
degree_cor1<-ggplot(data=tortoises_properties, aes(y=frugivory_degree, x=browsing_degree, alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) +theme_bw() + sm_statCorr(corr_method = "spearman")+
  labs(x="Degree Browsing and Grazing", y= "Degree Frugivory and SD") +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                               "T12"="#66B752")) 

strenght_cor1<-ggplot(data=tortoises_properties, aes(y=frugivory_strength, x=browsing_strength,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight),color=species)) +sm_statCorr(corr_method = "spearman") + theme_bw() +
  labs(x="Strength Browsing and Grazing", y= "Strength Frugivory and SD",alpha = NULL)+
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) 
  
closeness_cor1<-ggplot(data=tortoises_properties, aes(y=frugivory_closeness, x=round(browsing_closeness,4),alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(x="Closeness Browsing and Grazing", y= "Closeness Frugivory and SD")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))


degree_cor2<-ggplot(data=tortoises_properties, aes(y=detritivory_degree, x=browsing_degree, alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) +theme_bw() + sm_statCorr(corr_method = "spearman")+
  labs(x="Degree Browsing and Grazing", y= "Degree Detritivory") +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) +
  theme(legend.position = "none") 

strenght_cor2<-ggplot(data=tortoises_properties, aes(y=detritivory_strength, x=browsing_strength,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) +sm_statCorr(corr_method = "spearman") + theme_bw() +
  labs(x="Strength Browsing and Grazing", y= "Strength Detritivory")+
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) +
  theme(legend.position = "none") 

closeness_cor2<-ggplot(data=tortoises_properties, aes(y=detritivory_closeness, x=browsing_closeness,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(x="Closeness Browsing and Grazing", y= "Closeness Detritivory")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))


degree_cor3<-ggplot(data=tortoises_properties, aes(y=detritivory_degree, x=frugivory_degree, alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) +theme_bw() + sm_statCorr(corr_method = "spearman")+
  labs(x="Degree Frugivory and SD", y= "Degree Detritivory") +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) +
  theme(legend.position = "none") 

strenght_cor3<-ggplot(data=tortoises_properties, aes(y=detritivory_strength, x=frugivory_strength,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) +sm_statCorr(corr_method = "spearman") + theme_bw() +
  labs(x="Strength Frugivory and SD", y= "Strength Detritivory")+
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) +
  theme(legend.position = "none") 

closeness_cor3<-ggplot(data=tortoises_properties, aes(y=detritivory_closeness, x=frugivory_closeness,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=species)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(x="Closeness Frugivory and SD", y= "Closeness Detritivory")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))


correlations <- ggarrange(degree_cor2,strenght_cor2, closeness_cor2,
                          degree_cor3, strenght_cor3, closeness_cor3,
                          degree_cor1, strenght_cor1, closeness_cor1, ncol=3, nrow=3,
                          common.legend = TRUE)

ggsave(plot=correlations, filename = "Figures/Topological_Rol_Tortoises.pdf",  
       width = 29, height = 21, units="cm")

rm(list=c("ranking_tortoises", "correlations","closeness_cor","degree_cor",
          "strenght_cor", "weight",  "plot_closeness",
          "plot_degree","plot_strength","tortoise_browsing","tortoise_frugivory",
          "browsing_ranking","frugivory_ranking", "tortoise_ranking",
          "tortoise_degree","tortoise_strength","tortoise_closeness",
          "tortoise_detritivory")) #clean environment

rm(list=setdiff(ls(), c("browsing_matrix", "frugivory_matrix", "detritivory_matrix")))

######Structural role analysis #######

##### Browsing and grazing  #######
set.seed(1233)

mod_browsing<-bipartite::computeModules(browsing_matrix)
plotModuleWeb(mod_browsing, labsize=0.5) #5 modules
?plotModuleWeb
random_browsing<-bipartite::nullmodel(browsing_matrix, N=1000, method="vaznull")
modularities <- numeric(1000)
num_modules <- numeric(1000)

for (i in 1:1000) {
  modules <- computeModules(random_browsing[[i]])
  modularities[i] <- modules@likelihood
  rm(modules)
}

#observed modularity
mod_browsing@likelihood #0.34
#expected modularity
mean(modularities) #0.44
# deviate from random?: z-score < 2
(mod_browsing@likelihood - mean(modularities))/sd(modularities) #-8.38

cz_browsing<-czvalues(mod_browsing, weighted=TRUE, level="higher")
null.mod.list <- sapply(random_browsing, computeModules)
null.cz <- lapply(null.mod.list, czvalues)
null.cs <- sapply(null.cz, function(x) x$c) # c-values across all individuals in nulls
cthreshold_browsing<-quantile(null.cs, 0.95)
null.zs <- sapply(null.cz, function(x) x$z) # z-values across all individuals in nulls
zthreshold_browsing<-quantile(null.zs, 0.95,na.rm=TRUE)

#create data.frame from list
cz_browsing <- data.frame(
  Tortoise_ID = names(cz_browsing$c),   
  c_browsing = cz_browsing$c,           
  z_browsing = cz_browsing$z,
  Interaction= "Browsing"
)

rm(list=c("mod_browsing","null.cs","null.cz","null.mod.list","null.zs", 
          "random_browsing", "modularities","num_modules", "i"))
gc()


##### FRUGIVORY AND SEED DISPERSAL #######
mod_frugivory<-bipartite::computeModules(frugivory_matrix)
plotModuleWeb(mod_frugivory, labsize=0.5)
random_frugivory<-bipartite::nullmodel(frugivory_matrix, N=1000, method="vaznull")
modularities <- numeric(1000)
num_modules <- numeric(1000)


for (i in 1:1000) {
  modules <- computeModules(random_frugivory[[i]])
  modularities[i] <- modules@likelihood
  rm(modules)
}

#observed
mod_frugivory@likelihood # 0.44
#expected
mean(modularities) # 0.16
# deviate from random?: z-score < 2
(mod_frugivory@likelihood - mean(modularities))/sd(modularities) # 3.66


cz_frugivory<-czvalues(mod_frugivory, weighted=TRUE, level="higher")
null.mod.list <- sapply(random_frugivory, computeModules)
null.cz <- lapply(null.mod.list, czvalues)
null.cs <- sapply(null.cz, function(x) x$c) # c-values across all individuals in nulls
cthreshold_frugivory<-quantile(null.cs, 0.95)
null.zs <- sapply(null.cz, function(x) x$z) # z-values across all individuals in nulls
zthreshold_frugivory<-quantile(null.zs, 0.95,na.rm=TRUE)

#create data.frame from list
cz_frugivory <- data.frame(
  Tortoise_ID = names(cz_frugivory$c),   
  c_frugivory = cz_frugivory$c,           
  z_frugivory = cz_frugivory$z,
  Interaction= "Frugivory"
)

gc()
##### Detritivory #######
mod_detritivory<-bipartite::computeModules(detritivory_matrix)
plotModuleWeb(mod_detritivory, labsize=0.5)
random_detritivory<-bipartite::nullmodel(detritivory_matrix, N=1000, method="vaznull")
modularities <- numeric(1000)
num_modules <- numeric(1000)


for (i in 1:1000) {
  modules <- computeModules(random_detritivory[[i]])
  modularities[i] <- modules@likelihood
  rm(modules)
}

#observed
mod_detritivory@likelihood #  0.3
#expected
mean(modularities) # 0.13

# deviate from random?: z-score < 2
(mod_detritivory@likelihood - mean(modularities))/sd(modularities) # 8.88

cz_detritivory<-czvalues(mod_detritivory, weighted=TRUE, level="higher")
null.mod.list <- sapply(random_detritivory, computeModules)
null.cz <- lapply(null.mod.list, czvalues)
null.cs <- sapply(null.cz, function(x) x$c) # c-values across all individuals in nulls
cthreshold_detritivory<-quantile(null.cs, 0.95)
null.zs <- sapply(null.cz, function(x) x$z) # z-values across all individuals in nulls
zthreshold_detritivory<-quantile(null.zs, 0.95,na.rm=TRUE)

#create data.frame from list
cz_detritivory <- data.frame(
  Tortoise_ID = names(cz_detritivory$c),   
  c_detritivory = cz_detritivory$c,           
  z_detritivory = cz_detritivory$z,
  Interaction= "Detritivory"
)

cz_detritivory$z_detritivory<-ifelse(is.na(cz_detritivory$z_detritivory),0,cz_detritivory$z_detritivory)
cz_browsing$z_browsing<-ifelse(is.na(cz_browsing$z_browsing),0,cz_browsing$z_browsing)
cz_frugivory$z_frugivory<-ifelse(is.na(cz_frugivory$z_frugivory),0,cz_frugivory$z_frugivory)


weight <- read.csv("Data/Tortoise_Traits.csv", sep=";")
cz_detritivory<- left_join(cz_detritivory, weight, by="Tortoise_ID")

zc_detritivory <- ggplot(cz_detritivory, aes(y=z_detritivory, x=c_detritivory,)) + 
  geom_point(aes(size=log10(Weight), alpha=0.5, color=Tortoise_ID)) +theme_bw() + labs(x="Among-module connectivity (c)",
                                                                    y="Within-module connectivity (z)")+
  geom_vline(xintercept = cthreshold_detritivory, linetype="dashed", linewidth=0.5,colour="#393939")+
  geom_hline(yintercept = zthreshold_detritivory, linetype="dashed",linewidth=0.5, colour="#393939")+
  #geom_text(aes(label = Tortoise_ID), size = 3, hjust = -0.1, vjust = 1.5) +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))  

cz_browsing<- left_join(cz_browsing, weight, by="Tortoise_ID")
zc_browsing <- ggplot(cz_browsing, aes(y=z_browsing, x=c_browsing)) + 
  geom_point(aes(size=log10(Weight), alpha=0.5, color=Tortoise_ID)) +theme_bw() + labs(x="Among-module connectivity (c)",
                                                        y="Within-module connectivity (z)")+
  geom_vline(xintercept = cthreshold_browsing, linetype="dashed", linewidth=0.5,colour="#393939")+
  geom_hline(yintercept = zthreshold_browsing, linetype="dashed",linewidth=0.5, colour="#393939")+
  #geom_text(aes(label = Tortoise_ID), size = 3, hjust = -0.1, vjust = 1.5) +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) 

cz_frugivory<- left_join(cz_frugivory, weight, by="Tortoise_ID")
zc_frugivory <- ggplot(cz_frugivory , aes(y=z_frugivory , x=c_frugivory)) + 
  geom_point(aes(size=log10(Weight), alpha=0.5, color=Tortoise_ID)) +theme_bw() + labs(x="Among-module connectivity (c)",
                                                        y="Within-module connectivity (z)")+
  geom_vline(xintercept = cthreshold_frugivory , linetype="dashed", linewidth=0.5,colour="#393939")+
  geom_hline(yintercept = zthreshold_frugivory , linetype="dashed",linewidth=0.5, colour="#393939")+
  #geom_text(aes(label = Tortoise_ID), size = 3, hjust = -0.1, vjust = 1.5) +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752")) 
  

zc_plot <- ggarrange(zc_detritivory, zc_browsing,zc_frugivory, ncol=3)
zc_plot

cz_browsing<-cz_browsing %>% rename("c_browsing" =c_values, "z_browsing"=z_values)
cz_frugivory<-cz_frugivory %>% rename("c_frugivory"=c_values, "z_frugivory"=z_values)
cz_detritivory<-cz_detritivory %>% rename("c_detritivory"=c_values, "z_detritivory"=z_values)


cz_values <- left_join(weight, cz_browsing[,c("c_browsing","z_browsing","Tortoise_ID")], by="Tortoise_ID")
cz_values<- left_join(cz_values, cz_frugivory[,c("c_frugivory","z_frugivory","Tortoise_ID")], by="Tortoise_ID")
cz_values<- left_join(cz_values, cz_detritivory[,c("c_detritivory","z_detritivory","Tortoise_ID")], by="Tortoise_ID")

###

z1<-ggplot(data=cz_values, aes(y=z_detritivory, x=z_browsing,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="z Detritivory", x= "z Browsing and Grazing")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

z2<-ggplot(data=cz_values, aes(y=z_detritivory, x=z_frugivory,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="z Detritivory", x= "z Frugivory and SD")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

z3<-ggplot(data=cz_values, aes(y=z_browsing, x=z_frugivory,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="z Browsing and Grazing", x= "z Frugivory and SD")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

c1<-ggplot(data=cz_values, aes(y=c_detritivory, x=c_browsing,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="c Detritivory", x= "c Browsing and Grazing")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

c2<-ggplot(data=cz_values, aes(y=c_detritivory, x=c_frugivory,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="c Detritivory", x= "c Frugivory and SD")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

c3<-ggplot(data=cz_values, aes(y=c_browsing, x=c_frugivory,alpha=0.5)) + 
  geom_point(aes(size=log10(Weight), color=Tortoise_ID)) + sm_statCorr(corr_method = "spearman") +theme_bw() +
  labs(y="c Browsing and Grazing", x= "c Frugivory and SD")+
  theme(legend.position = "none") +
  scale_color_manual(values=c("T01"="#F2766D","T02"="#F1280D","T03"="#73AAC8",
                              "T05"="#A2D6C0","T06"="#FACE1B","T07"="#584897",
                              "T08"="#E89E58","T10"="#F86713","T11"="#86439A",
                              "T12"="#66B752"))

zc_plots <- ggarrange(zc_detritivory, zc_browsing,zc_frugivory,
                      z1,z2,z3,c1,c2,c3, ncol=3, nrow=3)


ggsave(plot=zc_plots, filename = "Figures/Structural_Rol_Tortoises.pdf",  
       width = 29, height = 21, units="cm")
gc()
