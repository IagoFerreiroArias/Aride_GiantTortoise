{
    "id": "F93021E6",
    "path": null,
    "project_path": null,
    "type": "r_dataframe",
    "hash": "0",
    "contents": "",
    "dirty": false,
    "created": 1736970195978.0,
    "source_on_save": false,
    "relative_order": 2,
    "properties": {
        "expression": "tracks",
        "caption": "tracks",
        "totalObservations": "3738",
        "displayedObservations": "3738",
        "variables": "8",
        "cacheKey": "EC24FC5C",
        "object": "tracks",
        "environment": "",
        "contentUrl": "grid_resource/gridviewer.html?env=&obj=tracks&cache_key=EC24FC5C&max_display_columns=50",
        "preview": "0",
        "source_window_id": "",
        "Source": "Source"
    },
    "folds": "",
    "lastKnownWriteTime": 2322206432585609839,
    "encoding": "",
    "collab_server": "",
    "source_window": "",
    "last_content_update": 1736970195978,
    "read_only": false,
    "read_only_alternatives": []
}                                                                             ts = "hour"))

#Create unique id to remove repeated rows and overestimate sampling effort
raw$id<- paste0(raw$Tortoise_ID,"_",raw$DateTime_Start)
effort <- raw[-which(duplicated(raw$id)),]
effort <- effort %>% group_by(Tortoise_ID) %>% summarise(Effort=round(sum(Effort),2))
write.csv(effort,"Data/CLEAN/Sampling_effort_Tortoise.csv",row.names = FALSE)

sum(effort$Effort)
mean(effort$Effort)
sd(effort$Effort)

# Estimate frequency of interactions: browsing and frugivory
str(raw)
levels(as.factor(raw$Plant_Material))
levels(as.factor(raw$Plant_Type))

raw$Plant_Species <- ifelse(raw$Plant_Species=="Kyllinga polyphylla ", "Kyllinga polyphylla", raw$Plant_Species)
raw$Plant_Species <- ifelse(raw$Plant_Species=="Ficus sp.", "Ficus reflexa", raw$Plant_Species)


str(raw)

behaviours <- raw %>%
  group_by(Tortoise_ID) %>%
  summarise(
    Frugivory = sum(Plant_Material == "Fruits", na.rm = TRUE),
    Browsing_Grazing = sum(Plant_Material == "Leaves" & Notes == "Fresh leaves", na.rm = TRUE),
    Detritivory = sum(Plant_Material == "Leaves" & Notes == "Fallen leaves", na.rm = TRUE)
  ) %>%
  mutate(
    Total = rowSums(select(., Frugivory, Browsing_Grazing, Detritivory), na.rm = TRUE)
  ) %>%
  mutate(Freq_Frugivory = round((Frugivory/Total)*100,2),
         Freq_BrowsingGrazing = round((Browsing_Grazing/Total)*100,2),
         Freq_Detritivory = round((Detritivory/Total)*100,2))

ggtern(data = behaviours, aes(x = Freq_Frugivory, y = Freq_BrowsingGrazing, z = Freq_Detritivory)) +
  geom_point(aes(color = Tortoise_ID), size = 4) +  # Añadir puntos con color por Tortoise_ID
  geom_text(aes(label = Tortoise_ID), size = 4, hjust = 0.5, vjust = -0.5) +  # Añadir etiquetas con Tortoise_ID
  labs(x = "Frugivory (%)", y = "Browsing and Grazing (%)", z = "Detritivory (%)") +
  theme_bw()  

rm(behaviours)

#Standarize browsing interactions by sampling effort
browsing_focal <- raw %>% filter(Plant_Material=="Leaves") %>% 
                      filter(Notes!="Fallen leaves") %>% #do not take into account detritivory
                      group_by(Tortoise_ID, Plant_Species) %>%
                      summarise(N_Interactions=n()) %>% 
                      filter(Plant_Species!="None")

browsing_focal <- left_join(effort,browsing_focal, by="Tortoise_ID")

#Standarize to 100 hours of sampling
browsing_focal$Freq_browsing <- (browsing_focal$N_Interactions/browsing_focal$Effort)* 100
write.csv(browsing_focal, "Data/CLEAN/Browsing_Focal_Observations.csv", row.names = FALSE)

#Standarize frugivory interactions by sampling effort
frugivory_focal <- raw %>% filter(Plant_Material=="Fruits") %>% 
  group_by(Tortoise_ID, Plant_Species) %>%
  summarise(N_Interactions=n()) %>% 
  filter(Plant_Species!="None")

frugivory_focal <- left_join(effort,frugivory_focal, by="Tortoise_ID")
frugivory_focal$Freq_frugivory <- (frugivory_focal$N_Interactions/frugivory_focal$Effort)*100
write.csv(frugivory_focal, "Data/CLEAN/Frugivory_Focal_Observations.csv", row.names = FALSE)

#Standarize frugivory interactions by sampling effort
detritivory_focal <- raw %>% filter(Plant_Material=="Leaves")%>% 
  filter(Notes=="Fallen leaves") %>%
  group_by(Tortoise_ID, Plant_Species) %>%
  summarise(N_Interactions=n()) %>% 
  filter(Plant_Species!="None")

detritivory_focal <- left_join(effort,detritivory_focal, by="Tortoise_ID")
detritivory_focal$Freq_detritivory <- (detritivory_focal$N_Interactions/detritivory_focal$Effort)*100
write.csv(detritivory_focal, "Data/CLEAN/Detritivory_Focal_Observations.csv", row.names = FALSE)

unique(detritivory_focal$Plant_Species)
