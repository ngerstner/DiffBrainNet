##################################################
## Project: DexStim Mouse Brain
## Date: 30.11.2021
## Author: Nathalie
##################################################
# Numbers/Percentages of DE and hub genes in vCA1
# including also number of unique ones

library(ggplot2)

# correct numbers? --> DE percentage different than before
DE_genes <- 465
hub_genes <- 260

unique_DE <- 25
unique_hub <- 24

# data frame for plotting
df <- data.frame("type" = c("DE", "DE", "hub", "hub"),
                 "unique" = c(TRUE, FALSE, TRUE, FALSE),
                 "number" = c(unique_DE, DE_genes - unique_DE,
                              unique_hub, hub_genes - unique_hub),
                 "text" = c(paste0(round(unique_DE/DE_genes*100,1),"%"), "",
                            paste0(round(unique_hub/hub_genes*100,1), "%"), ""))

# barplot
ggplot(df, aes(x = type, y = number, fill = type)) +
  geom_bar(position = "stack", stat="identity", aes(alpha = unique)) +
  geom_text(aes(y = number, label = text), vjust = 1.5) +
  scale_alpha_manual("",
                     breaks = c(FALSE, TRUE),
                     values = c(1.0,0.5),
                     labels = c("Shared genes", "Unique genes")) +
  xlab("") +
  ylab("# DE and hub genes in vCA1") +
  scale_fill_manual("",
                    breaks = c("DE", "hub"),
                    labels = c("DE genes", "Hub genes"),
                    values = c("orange", "darkred")) +
  scale_x_discrete(breaks = c("DE", "hub"),
                   labels = c("DE genes", "Hub genes")) +
  theme_light() +
  theme(text = element_text(size= 14)) +
  guides(fill = "none")

ggsave(filename = "08_FigureIII_A.pdf", width = 8, height = 6)
