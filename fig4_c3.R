setwd("./nextalign_4_derep_mpf/")
################################################################################################################################

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(phylotate)
library(tidytree)
library(ggnewscale)
library(stringr)
library(dplyr)
library(RColorBrewer)


tre <- read.nexus("4_joint_mugration_region/annotated_tree.nexus")
meta <- read.csv("../metadata_dataset210303_clade.tsv", sep = "\t", header = T)
head(meta)

t <- as_tibble(tre)

u <- c("Uruguay" = "#495057", "EpiCoV/GISAID" = "#e9ecef")


library(RColorBrewer)
coul <- brewer.pal(11, "Spectral") 
coul <- colorRampPalette(coul)(16)


head(meta)
colnames(meta)[1] <- c("label")
meta <- merge(t, meta, by = "label")
meta <- as.data.frame(meta[order(meta$date, decreasing = T),])


p <- ggtree(tre, size =0.1, color = "black", layout = "fan", mrsd = "2021-02-23") %<+% meta +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  geom_rootpoint(position = "identity") +
  labs(colour="Clado") +
  guides(colour = guide_legend(override.aes = list(size = 10))) 
p


clade <- select(meta, label, clade)
clade <- clade[!duplicated(clade$label),]
row.names(clade) <- clade$label
clade$label <- NULL

count <- select(meta, label, Region)
count <- count[!duplicated(count$label),]
row.names(count) <- count$label
count$label <- NULL
head(count)

count <- count %>% 
  mutate(dataset = ifelse(grepl("Uruguay", Region), "Uruguay", "EpiCoV/GISAID"))

count$Region <- NULL

p1 <- gheatmap(p, count, offset=0.01, width=.1, colnames = F, color = FALSE) + 
  scale_fill_manual(values = u, name = "Dataset")

p1


p2 <- p1 + new_scale_fill()
t <- gheatmap(p2, clade, offset=0.1, width=.1, colnames = F, color = F) +
  scale_fill_manual(values = coul, name = "Clado")
t


png("../Figuras/Figura4.png", res = 600, height = 20, width = 20, units = 'cm')
t
dev.off()

