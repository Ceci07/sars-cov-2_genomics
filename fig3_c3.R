library(stringr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggsignif)
library(ggpubr)

setwd("./")

df <- read.table("metadata_uy_tree.tsv", sep = "\t", header = T) 
head(df)
unique(df$clade)
df$status <- gsub("Este estudio", "Institut Pasteur de Montevideo", df$status)
df$status <- gsub("EpiCov/GISAID", "EpiCoV/GISAID", df$status)

list <- df %>%
  group_by(accession, status, lineage, clade, aaSubstitutions) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)


library(RColorBrewer)
coul <- brewer.pal(2, "Set3") 
coul <- colorRampPalette(coul)(2)

coul <- c("Institut Pasteur de Montevideo" = "#00916e",
          "EpiCoV/GISAID" = "#ffcf00")
 

l <- ggplot(list, aes(axis1 = ID, axis2= Linaje, axis3=Clado, y = n)) +
  geom_alluvium(aes(fill = status), aes.bind=F, width = 1/4) +
  geom_stratum(width = 1/4, fill = "white", color = "black") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size =5) +
  scale_x_discrete(limits = c("Muestra", "Linaje PANGO", "Clado Nextstrain"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = coul) +
  guides(fill=guide_legend(title="")) +
  labs(y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20)) + 
  ggtitle("") 
l


tab <- as.data.frame(table(df$aaSubstitutions, df$date))
tab <- tab %>% filter(Freq > 0)
tab[order(tab$Var2, decreasing = F),]
unique(tab$Var1)
head(tab)
colnames(tab) <- c("aaSubstitutions", "date", "Freq")

clade <- select(df, aaSubstitutions, clade)

head(clade)
head(tab)
tab <- merge(tab, clade, by = "aaSubstitutions")


cl <- c("ORF8:L84S" = "#f94144", 
        "N:S197L,N:S310N,ORF1a:F3071Y,ORF3a:G196V,ORF8:L84S,S:I197V" ="#f3722c",
        "N:S197L,N:S310N,ORF1a:F3071Y,ORF8:L84S" ="#f8961e", 
        "ORF1a:P1096S,ORF1b:P314L,S:D614G" = "black",
        "ORF1a:R550S,ORF8:L84S" ="#f9c74f",
        "ORF1a:T2457K,ORF8:L84S" ="#90be6d",
        "ORF1b:P1427L,ORF1b:Y1464C,ORF8:L84S" = "#43aa8b",
        "ORF1b:P314L,S:D614G" = "#4d908e",
        "N:D22G,ORF1a:S944L,ORF1a:L3606F,ORF1b:D54Y,ORF1b:V1493L,ORF3a:G251V,ORF9b:I19V" ="#577590",
        "ORF3a:G11E,ORF8:L84S" ="#277da1", 
        "N:R203K,N:G204R,N:I292T,ORF1b:P314L,ORF6:I33T,S:D614G" = "#f4f1de",
        "N:R203K,N:G204R,ORF1a:G3278S,ORF1a:L3606F,ORF1b:P314L,S:D614G" = "#e07a5f",
        "N:R203K,N:G204R,ORF1a:T1246I,ORF1a:G3278S,ORF1a:L3606F,ORF1b:P314L,S:D614G" = "#4f772d",
        "N:R203K,N:G204R,ORF1a:T1246I,ORF1a:G3278S,ORF1a:L3606F,ORF1b:P314L,S:D614G", "#3d405b",
        "N:T265I,ORF1a:L3606F,ORF3a:G251V" ="#81b29a")

coul <- brewer.pal(12, "Set3") 
coul <- colorRampPalette(coul)(14)


p <- ggplot(tab, aes(fill=aaSubstitutions, y=Freq, x=date)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = coul) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  ylab("Abundancia realtiva")  + xlab("Fecha") +
  guides(fill=guide_legend(ncol=2, title="Sustituciones AA"))
p  


fig <- ggarrange(l, p, nrow = 2, labels = c("A", "B"), font.label=list(color="black",size=25))
fig

dir.create("Figuras")
png("Figuras/Figura3.png", res = 600, height = 50, width = 55, units = 'cm')
fig
dev.off()
