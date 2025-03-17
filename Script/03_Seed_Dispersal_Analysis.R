# Scripts: Seed dispersal of Aldabra giant tortoise on Aride
# Author: Iago Ferreiro-Arias
# Date: 5th May, 2024

library(ggplot2)
library(dplyr)
library(lme4)
library(sjPlot)
library(tidyr)
library(ggpubr)
library(effect.lndscp) #devtools::install_github("pedroj/effectiveness_pckg")

#Import data
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


#Seed dispersal effectiveness of different individuas
native <- c("Ficus_reflexa", "Morinda_citrifolia", "Cordia_subcordata", "Megathyrsus_maximus")
introduced <- c("Terminalia_catappa", "Phillantus_amarus", "Syzygium_cumini", "Ficus_benghalensis", "Ananas_comosus")

str(seeds)
sde<- seeds %>%
  group_by(Tortoise_ID) %>%
  summarise(Total_Native = sum(across(all_of(native)), na.rm = TRUE),
            Total_Seeds = sum(N_seeds, na.rm = TRUE),
            Prop_Native = (Total_Native / Total_Seeds)*100)

effort <- read.csv("Data/CLEAN/Sampling_effort_Tortoise.csv")
sde<- left_join(sde,effort,by="Tortoise_ID")
sde$Seeds_Rate <- (sde$Total_Seeds / 3) # 3 scats per tortoise

sde <- sde %>% 
  add_row(Tortoise_ID = "TXX", Prop_Native = 0, Seeds_Rate =316200)%>% # plot isolines to 100% of x axis
  add_row(Tortoise_ID = "TX1", Prop_Native = 100, Seeds_Rate =316200) # plot isolines to 100% of x axis
# plot isolines to 100% of x axis
#log10(x) = 5.5 => x= 316200

seed_efe<-effectiveness_plot(sde$Prop_Native, log10(sde$Seeds_Rate),
                   label = sde$Tortoise_ID,  
                   myxlab = "Proportion of native plants (%)", 
                   myylab = "log10(No. seeds / sample)") +
                   scale_x_continuous(limits = c(0, 100)) +
                   scale_y_continuous(limits=c(0, 5.5),
                                      breaks=c(0, 1, 2, 3, 4, 5))

sde_merged <-ggpubr::ggarrange(det_eff,brow_eff, seed_efe, nrow=3)
ggsave(plot=sde_merged, filename="Effectiveness.pdf", width=14, height = 21, units="cm")

#Relationship between number of seeds consumed and weight of the tortoise
model <- glm(Total_seeds~Weight, family = "poisson", data=tortoise_seeds)
summary(model)
g1 <- sjPlot::plot_model(model, terms= "Weight", type="eff", color="#8aba97", show.values = TRUE) + theme_bw() +
        ylab("Number of seeds consumed") + xlab("Tortoise weigth (kg)") +
  theme(plot.title = element_blank())
?plot_model()
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
hist(seeds_rich$Seed_Richness)
model5 <- glm(Seed_Richness~log10(Weight), family = "poisson", data=seeds_rich)
summary(model5)
g5 <- sjPlot::plot_model(model5, terms= "Weight", type="pred", color="#8aba97") + theme_bw() +
  ylab("Seed richness on tortoise diet") + xlab("Tortoise weigth (kg)") +
  theme(plot.title = element_blank())
performance::r2(model5)
rm(model5)

# Relationship between seed richness on tortoise diet and size of feeding areas
model6 <- glm(Seed_Richness~Feeding_Area, family = "poisson", data=seeds_rich)
summary(model6)
g6 <- sjPlot::plot_model(model6, terms= "Feeding_Area", type="pred", color="#8aba97") + theme_bw() +
  ylab("Seed richness on tortoise diet") + xlab("Feeding area (m2)") +
  theme(plot.title = element_blank())
g6
performance::r2(model6)
rm(model6)

ggarrange(g1,g2,g3,g4,g5,g6,ncol=3, nrow=2)

### Rearrange data.frame for adjacency matrix 

str(seeds)

seeds <- seeds %>%
  dplyr::select(-Scat_WetWeight, -N_seeds) %>%  
  pivot_longer(cols = -c(Tortoise_ID, Sample_ID),  
               names_to = "Plant_ID",  
               values_to = "Seeds")  %>%
     group_by(Tortoise_ID, Plant_ID) %>% summarise(Seeds=mean(Seeds))

effort <- read.csv("Data/CLEAN/Sampling_effort_Tortoise.csv")
seeds<- left_join(seeds,effort,by="Tortoise_ID")
seeds$Seeds_Frequency <- (seeds$Seeds/seeds$Effort)*100
write.csv(seeds, "Data/CLEAN/Seeds_Frequency.csv", row.names = FALSE)

rm(list=c("g1","g2","g3","g4","g5","g6", "seeds", "seeds_rich", "summary_seeds", 
          "tortoise_seeds"))
rm(list = ls())


