# Required packages
library(ggplot2)
library(dplyr)

##############################################################
# Power (RCRXO vs CRXO)
##############################################################
simulation_results_RCRXO <- read.csv("simulation_results_RCRXO.csv")
simulation_results_CRXO <- read.csv("simulation_results_CRXO.csv")

filtered_data_RCRXO_4 <- simulation_results_RCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 4) %>%
  mutate(Design = "RCRXO-4")

filtered_data_RCRXO_5 <- simulation_results_RCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 5) %>%
  mutate(Design = "RCRXO-5")

filtered_data_RCRXO_6 <- simulation_results_RCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 6) %>%
  mutate(Design = "RCRXO-6")

filtered_data_CRXO <- simulation_results_CRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  mutate(Design = "CRXO")

ggplot() +
  geom_line(data = filtered_data_RCRXO_4, aes(x= ICC, y = Alpha, group = Design, color = "RCRXO-4", linetype = "RCRXO-4")) +
  geom_line(data = filtered_data_RCRXO_5, aes(x= ICC, y = Alpha, group = Design, color = "RCRXO-5", linetype = "RCRXO-5")) +
  geom_line(data = filtered_data_RCRXO_6, aes(x= ICC, y = Alpha, group = Design, color = "RCRXO-6", linetype = "RCRXO-6")) +
  geom_line(data = filtered_data_CRXO, aes(x= ICC, y = Alpha, group = Design, color = "CRXO-4", linetype = "CRXO-4")) +
  facet_wrap(~ CAC, scales = "fixed", labeller = label_bquote(cols = "CAC = "* .(CAC))) +
  labs(title = "Statistical Power Across CAC Levels",
       x = "ICC",
       y = "Power",
       color = "Design",
       linetype = "Design") +  # Legend for design type
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.title.x = element_text(margin = margin(t = 10), size = 12)) + 
  scale_color_manual(values = c("RCRXO-4" = "#1f77b4", "RCRXO-5" = "#6A0DAD",
                                "RCRXO-6" = "#800000", "CRXO-4" = "black")) +
  scale_linetype_manual(values = c("RCRXO-4" = "solid", "RCRXO-5" = "solid",
                                   "RCRXO-6" = "solid", "CRXO-4" = "dotted"))




##############################################################
# Power (MCRXO vs CRXO)
##############################################################

simulation_results_MCRXO <- read.csv("simulation_results_MCRXO.csv")
simulation_results_CRXO <- read.csv("simulation_results_CRXO.csv")

filtered_data_MCRXO_4 <- simulation_results_MCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 4) %>%
  mutate(Design = "MCRXO-4")

filtered_data_MCRXO_6 <- simulation_results_MCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 6) %>%
  mutate(Design = "MCRXO-6")

filtered_data_MCRXO_8 <- simulation_results_MCRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  filter(Time == 8) %>%
  mutate(Design = "MCRXO-8")

filtered_data_CRXO <- simulation_results_CRXO %>%
  filter(TreatmentEffect == 0.2) %>%
  mutate(Design = "CRXO")

ggplot() +
  geom_line(data = filtered_data_MCRXO_4, aes(x= ICC, y = Alpha, group = Design, color = "MCRXO-4", linetype = "MCRXO-4")) +
  geom_line(data = filtered_data_MCRXO_6, aes(x= ICC, y = Alpha, group = Design, color = "MCRXO-6", linetype = "MCRXO-6")) +
  geom_line(data = filtered_data_MCRXO_8, aes(x= ICC, y = Alpha, group = Design, color = "MCRXO-8", linetype = "MCRXO-8")) +
  geom_line(data = filtered_data_CRXO, aes(x= ICC, y = Alpha, group = Design, color = "CRXO-4", linetype = "CRXO-4")) +
  facet_wrap(~ CAC, scales = "fixed", labeller = label_bquote(cols = "CAC = "* .(CAC))) +
  labs(title = "Statistical Power Across CAC Levels",
       x = "ICC",
       y = "Power",
       color = "Design",
       linetype = "Design") +  # Legend for design type
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.title.x = element_text(margin = margin(t = 10), size = 12)) + 
  scale_color_manual(values = c("MCRXO-4" = "#1f77b4", "MCRXO-6" = "#6A0DAD",
                                "MCRXO-8" = "#800000", "CRXO-4" = "black")) +
  scale_linetype_manual(values = c("MCRXO-4" = "solid", "MCRXO-6" = "solid",
                                   "MCRXO-8" = "solid", "CRXO-4" = "dotted"))
