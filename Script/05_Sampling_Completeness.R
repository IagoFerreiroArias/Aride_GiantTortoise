library(vegan)
library(ggplot2)


effort <- read.csv("Data/Sampling_effort_Tortoise.csv",sep=",")
str(effort)
str(t(browsing_matrix))
str(t(frugivory_matrix))

browsing_accum <- specaccum(t(browsing_matrix), method = "random")
browsing_accum <- data.frame(Tortoises = browsing_accum$sites,
                           Interactions = browsing_accum$richness,
                           Lower = browsing_accum$richness - browsing_accum$sd,
                           Upper = browsing_accum$richness + browsing_accum$sd,
                           "Interaction"="Browsing and Grazing")
browsing_accum$Hours <- cumsum(effort$Effort)

dim(t(frugivory_matrix))
frugivory_accum <- specaccum(t(frugivory_matrix), method = "random")
frugivory_accum <- data.frame(Tortoises = frugivory_accum$sites,
                            Interactions = frugivory_accum$richness,
                             Lower = frugivory_accum$richness - frugivory_accum$sd,
                             Upper = frugivory_accum$richness + frugivory_accum$sd,
                             "Interaction"="Frugivory and Seed Dispersal")
frugivory_accum$Hours <- cumsum(effort$Effort)

detritivory_accum <- specaccum(t(detritivory_matrix), method = "random")
detritivory_accum <- data.frame(Tortoises = detritivory_accum$sites,
                              Interactions = detritivory_accum$richness,
                              Lower = detritivory_accum$richness - detritivory_accum$sd,
                              Upper = detritivory_accum$richness + detritivory_accum$sd,
                              "Interaction"="Detritivory")
detritivory_accum$Hours <- cumsum(effort$Effort)

accum_data <- rbind(detritivory_accum, browsing_accum,frugivory_accum)

g1<-ggplot(accum_data, aes(x = Tortoises, y = Interactions, color = Interaction)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Interaction), alpha = 0.4) +
  geom_line(linewidth = 1) + theme_bw() +
  labs(x = "Number of individuals sampled", y = "Number of distinct pairwise interactions",
       fill="Interaction type", colour="Interaction type") +
  scale_color_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772")) +
  scale_fill_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772"))+
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(1, 10)) +  
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80)) +
  theme(legend.position = "top", legend.direction = "horizontal")

ggplot(accum_data, aes(x = Hours, y = Interactions, color = Interaction, group = Interaction)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Interaction), alpha = 0.2) +
  geom_line(size = 1) + theme_bw() +
  labs(x = "Number of hours following individual tortoises",
       y = "Number of distinct pairwise interactions",
       fill="Interaction type", colour="Interaction type") +
  scale_color_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772")) +
  scale_fill_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772"))+
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80)) +
  theme(legend.position = "top", legend.direction = "horizontal")

g2<-ggplot(accum_data, aes(x = Tortoises, y = Interactions, color = Interaction)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Interaction), alpha = 0.4) +
  geom_line(linewidth = 1) + theme_bw() +
  labs(x = "Number of tortoise individuals sampled", y = "Number of distinct pairwise interactions",
       fill="Interaction type", colour="Interaction type") +
  scale_color_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772")) +
  scale_fill_manual(values = c("Detritivory"="#cca983","Browsing and Grazing" = "#c0c7b7", "Frugivory and Seed Dispersal" = "#f2c772"))+
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(1, 10),  
    name = "Sampling effort (hours)", 
    labels = c("28.45", "60.53", "66.16", "96.58", "107.85", "120.93", "156.06", "175.56", "204.04", "215.82"))+
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80)) +
  theme(legend.position = "top", legend.direction = "horizontal")

ggpubr::ggarrange(g1,g2,ncol=2, common.legend = TRUE)

### Venn diagram depicting interaction detection for different methods
int_det <- read.csv("Data/CLEAN/Interaction_Detection_Methods.csv", sep=";")
str(int_det)

library(ggvenn)

det_collected <- int_det$Plant_Species[int_det$Collected == 1]
det_metabarcoding <- int_det$Plant_Species[int_det$Metabarcoding == 1]
det_focal_observation <- int_det$Plant_Species[int_det$Focal_Observation == 1]
det_faeces_examination <- int_det$Plant_Species[int_det$Faeces_Examination == 1]

detected_methods <- list(
  Collected = det_collected,
  Metabarcoding = det_metabarcoding,
  Focal_Observation = det_focal_observation,
  Faeces_Examination = det_faeces_examination
)

# Generar el diagrama de Venn
ggvenn(detected_methods, fill_color = c("#36313d","#73648c","#e39f4b", "#994854"), show_percentage = FALSE)

