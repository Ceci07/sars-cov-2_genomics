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
library(tidyr)


setwd("./nextalign_4_derep_mpf/")
tre <- read.nexus("7_marginal_mugration_region/annotated_tree.nexus")
meta <- read.csv("../metadata_dataset210303_clade.tsv", sep = "\t", header = T)
head(meta)

t <- as_tibble(tre)
u <- c("Uruguay" = "#495057", "EpiCoV/GISAID" = "#e9ecef")


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


meta <- meta %>%
  mutate(detection = case_when(
    endsWith(clade, "19A") ~ "Detectado",
    endsWith(clade, "19B") ~ "Detectado",
    endsWith(clade, "20B") ~ "Detectado",
    endsWith(clade, "20D") ~ "Detectado"))

meta$detection <- meta$detection %>% replace_na('Sin detectar')
meta <- meta[!duplicated(meta$label),]

det <- select(meta, detection, label)
row.names(det) <- det$label
det$label <- NULL
head(det)

v <- c("#6c757d", "white")

p3 <- t + new_scale_fill()
tt <- gheatmap(p3, det, offset=0.2, width=.1, colnames = F, color = F) +
  scale_fill_manual(values = v, name = "Detección en marzo 2020")
tt


##############################################################################################

cls <- c("#ff0a47", "#0039a3", "#bfde26", "#00b2ea", "#60c888", "#8734d5", "#e50056", "#ff7b30", "#354f52",
         "#f94144", "#f3722c", "#f8961e", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1",
         "#ea88a3", "#dd365b", "#a32b45", "#f2f1ff", "#bfbdd7", "#f4ecbf", "#fed92e", "#c8a221", "#e2547b", "#ad3b64")


tre <- read.beast("9_joint_mugration_region/annotated_tree.nexus")
t <- as.data.frame(as_tibble(tre))
t[grepl("Uruguay", t$Region),]


d <- ggtree(tre, aes(color = Region), size =1, mrsd = "2021-02-23") +
  scale_colour_manual(values = cls) +
  theme_tree2() + 
  geom_rootpoint(position = "identity") +
  xlim_tree(xlim = 2021) 
d  + geom_text(aes(label=node), hjust=-.3)


f <- d %<+% meta +
  geom_tippoint(size = 2) +
  geom_rootpoint(position = "identity") +
  guides(colour = guide_legend("Region", override.aes = list(size = 2, shape=19))) +
  scale_x_ggtree(labels = c("2020-01", "2020-07", "2021-01", "2021-07"), breaks = c(2020.0, 2020.5, 2021.0, 2021.5)) +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) 
f 

g2 <-f + new_scale_fill()

gg <- gheatmap(g2, clade, offset=0.15, width=0.05, colnames = FALSE, color = FALSE) +
  scale_y_continuous(expand=c(0, 0.4)) +
  scale_fill_manual(values = coul) +
  labs(fill = "Clado") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
  scale_x_ggtree(labels = c("2020-01", "2020-07", "2021-01"), breaks = c(2020.0, 2020.5, 2021.0)) 
gg


meta <- meta[!duplicated(meta$label),]
meta2 <- meta

rownames(meta2) <- meta[,1]
clade2 <- rep('Other', dim(meta2)[1])

library(data.table)
target <- meta$clade %like% '19A|19B|20B|20D'
clade2[target] <- 'target'
meta2$clade2 <- clade2


target2 <- meta$clade %like% '20A'
clade2[target2] <- 'target2'
meta2$clade2 <- clade2

target3 <- meta$clade %like% '20D'
clade2[target3] <- 'target3'
meta2$clade2 <- clade2

target4 <- meta$clade %like% '19A|19B|20A'
clade2[target4] <- 'target4'
meta2$clade2 <- clade2

target5 <- meta$clade %like% '20B'
clade2[target5] <- 'target5'
meta2$clade2 <- clade2

target6 <- grepl("Gamma", meta$clade)
clade2[target6] <- 'target6'
meta2$clade2 <- clade2

grp <- list()
uph <- unique(meta2$clade2)
oth <- c()
nam <- c()
r   <- 1

for (u in 1:length(uph)) {
  
  w <- which(meta2$clade2 == uph[u])
  
  if (length(w) < 1) {
    
    oth <- c(oth, w)
    
  } else {
    
    s <- meta2[w,1]
    grp[[r]] <- s
    nam <- c(nam, uph[u])
    r <- r + 1
  }
}

names(grp) <- uph

other.clade <- MRCA(tre,  grp$Other) 

head(meta)


e <- collapse(d, node =3444) %<+% meta +
  geom_tippoint(size = 2) +
  geom_rootpoint(position = "identity") +
  geom_point2(aes(subset=(node==3444)), shape=23, size=6, fill="white", color = "black") +
  guides(colour = guide_legend("Region", override.aes = list(size = 2, shape=19))) +
  scale_x_ggtree(labels = c("2020-01", "2020-07", "2021-01", "2021-07"), breaks = c(2020.0, 2020.5, 2021.0, 2021.5)) +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) 
e

meta <- meta %>% 
  mutate(level = case_when(str_detect(country, "Uruguay" ) ~ "Uruguay",
                           TRUE ~ "EpiCoV/GISAID"))

meta <- meta[!duplicated(meta$label),]

dataset <- select(meta, label, level)
row.names(dataset) <- dataset$label
dataset$label <- NULL
head(dataset)

u <- c("Uruguay" = "#495057", "EpiCoV/GISAID" = "#e9ecef")

head(dataset)

g1 <- gheatmap(e, dataset, offset=0.1, width=0.05, colnames= F, color = F) +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 18), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = u) + 
  labs(fill = "Muestras")
g1 

head(dataset)
head(clade)

g2 <- g1 + new_scale_fill()

gg <- gheatmap(g2, clade, offset=0.15, width=0.05, colnames = F, color = F) +
  scale_y_continuous(expand=c(0, 1.4)) +
  scale_fill_manual(values = c("#e76f51", "#2a9d8f", "#e9c46a", "#3a0ca3", "#0081a7")) +
  labs(fill = "Clado") +
  theme(axis.title.x = element_text(size = 12), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))  
gg 




png("../Figuras/Figura5.png", res = 600, height = 40, width = 40, units = 'cm')
gg
dev.off()
