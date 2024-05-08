# Scripts: Seed dispersal of Aldabra giant tortoise on Aride
# Author: Iago Ferreiro-Arias
# Date: 5th May, 2024

library(ggplot2)
library(dplyr)
library(lme4)
library(sjPlot)
library(tidyr)
library(ggpubr)

#Import viridis#Import data
tortoises <- read.csv("Data/Tortoise_Traits.csv", sep=";")
seeds <- read.csv("Data/Seeds_Consumed.csv", sep=",")

#Calculate total number of seeds per scat
seeds$N_seeds <- rowSums(seeds[,4:12])
summary(seeds$N_seeds) #descriptive stats of number of seeds per scat

#Calculate total seeds per individual

tortoise_seeds <- seeds %>% group_by(Tortoise_ID) %>% summarise("Total_seeds"=sum(N_seeds), 
                                                                "Scat_Size"=mean(Scat_WetWeight))
tortoise_seeds <- left_join(tortoise_seeds, tortoises, by="Tortoise_ID")
summary(tortoise_seeds$Total_seeds)

#Relationship between number of seeds consumed and weight of the tortoise
model <- glm(Total_seeds~Weight, family = "poisson", data=tortoise_seeds)
summary(model)
g1 <- sjPlot::plot_model(model, terms= "Weight", type="pred", color="#8aba97") + theme_bw() +
        ylab("Number of seeds consumed") + xlab("Tortoise weigth (kg)") +
  theme(plot.title = element_blank())
performance::r2(model)
rm(model)

#Relationship between number of seeds consumed and weight of the scat
model2 <- glm(Total_seeds~Scat_Size, family = "poisson", data=tortoise_seeds)
summary(model2)
g2<-sjPlot::plot_model(model2, terms= "Scat_Size", type="pred", color="#8cdbce") + theme_bw() +
  ylab("Number of seeds consumed") + xlab("Wet weight of feacal sample (gr)") +
  theme(plot.title = element_blank())
performance::r2(model2)
rm(model2)

#Relationship between scat weight and weight of the tortoise 
hist(log10(tortoise_seeds$Scat_Size))
model3 <- lm(log10(Scat_Size) ~ log10(Weight) ,data=tortoise_seeds)
summary(model3)
g3<-sjPlot::plot_model(model3, terms= "Weight", type="pred", color="#697dd1") + theme_bw() +
  ylab("Wet weight of faecal sample (gr)") + xlab("Tortoise weight (kg)") +
  theme(plot.title = element_blank())
performance::r2(model3)
rm(model3)

# Relationship between individual feeding areas and tortoise weight
model4 <- lm(log10(Feeding_Area/1000000) ~ log10(Weight) ,data=tortoises)
summary(model4)
g4<-sjPlot::plot_model(model4, terms= "Weight", type="pred", color="#6b8185") + theme_bw() +
  ylab("Feeding area (km2)") + xlab("Tortoise weight (kg)") +
  theme(plot.title = element_blank())
performance::r2(model4)
rm(model4)

# Create summary tamble of seeds consumed by each tortoise and each plant species
summary_seeds <- seeds %>% group_by(Tortoise_ID) %>%
  pivot_longer(cols = 4:12, names_to = "Species", values_to = "N_Seeds_Species")

#Order levels of tortoises
ggplot(summary_seeds, aes(y=Tortoise_ID, x=N_Seeds_Species, fill=Species)) + geom_col() +
  theme_bw() + ylab("Tortoise ID") + xlab("Number of seeds consumed") + labs(fill="Plant species")+
  scale_fill_manual(values=c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", 
                             "#fee08b", "#abdda4","#e6f598", "#66c2a5", "#3288bd"))

#Species richness of seeds in function of feeding areas

seeds_rich <- summary_seeds %>% group_by(Tortoise_ID, Species)%>% 
              summarise("N_Seeds_Species"=sum(N_Seeds_Species))

seeds_rich <- seeds_rich  %>% group_by(Tortoise_ID) %>%
  summarise("Seed_Richness" = sum(N_Seeds_Species != 0, na.rm = TRUE))

seeds_rich <- left_join(seeds_rich, tortoises, by="Tortoise_ID")

# Relationship between seed richness on tortoise diet and tortoise weight 
model5 <- glm(Seed_Richness~Weight, family = "poisson", data=seeds_rich)
summary(model5)
g5 <- sjPlot::plot_model(model5, terms= "Weight", type="pred", color="#8aba97") + theme_bw() +
  ylab("Seed richness on tortoise diet") + xlab("Tortoise weigth (kg)") +
  theme(plot.title = element_blank())
performance::r2(model5)
rm(model5)

# Relationship between seed richness on tortoise diet and size of feeding areas
model6 <- glm(Seed_Richness~Feeding_Area, family = "poisson", data=seeds_rich)
summary(model6)
g6 <- sjPlot::plot_model(model6, terms= "Weight", type="pred", color="#8aba97") + theme_bw() +
  ylab("Seed richness on tortoise diet") + xlab("Feeding area (m2)") +
  theme(plot.title = element_blank())
g6
performance::r2(model6)
rm(model6)

ggarrange(g1,g2,g3,g4,g5,g6,ncol=3, nrow=2)
rm(list=c("g1","g2","g3","g4","g5","g6", "seeds", "seeds_rich", "summary_seeds", 
          "tortoise_seeds"))
