
# normalised enrichment score vs pathway plotted

library(tidyverse)


# Plot cluster 2, ONS76-CBO (non-malignant granule neurons) FGSEA
fgseaResTidy.clus2.GOBP <- 
  readRDS("outputs50K/ONS76-CBO/cluster2.ONS76Cbo.multifgsea.presto.GOBiolProcessGenes.rds")
o76cbo.clus2.ggplot <- ggplot(fgseaResTidy.clus2.GOBP %>% filter(padj < 0.05) %>% head(n=20), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "") + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 30)) # can change bar colours using ggplot

ggsave(file = "./outputs50K/ONS76-CBO/figures/NESplotCLust2FGSEAons76Cbo.png",
       units = c("in"),
       width = 25,
       height = 8,
       dpi = 300,
       o76cbo.clus2.ggplot)

# Plot cluster 3, ONS76-CBO (malignant granule precursors) FGSEA
fgseaResTidy.clus3.GOBP <- 
  readRDS("outputs50K/ONS76-CBO/cluster3.ONS76Cbo.multifgsea.presto.GOBiolProcessGenes.rds")

# only plot the top 20 pathways
o76cbo.clus3.ggplot <- ggplot(fgseaResTidy.clus3.GOBP %>% filter(padj < 0.05) %>% head(n=20), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "") + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 30)) # can change bar colours using ggplot

ggsave(file = "./outputs50K/ONS76-CBO/figures/NESplotClust3FGSEAons76Cbo.png",
       units = c("in"),
       width = 25,
       height = 20,
       dpi = 300,
       o76cbo.clus3.ggplot)

# Plot cluster 10, ONS76-CBO (malignant granule precursors) FGSEA
fgseaResTidy.clus10.GOBP <- 
  readRDS("outputs50K/ONS76-CBO/cluster10.ONS76Cbo.multifgsea.presto.GOBiolProcessGenes.rds")

# only plot the top 20 pathways
o76cbo.clus10.ggplot <- ggplot(fgseaResTidy.clus10.GOBP %>% filter(padj < 0.05) %>% head(n=20), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "") + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 30)) # can change bar colours using ggplot

ggsave(file = "./outputs50K/ONS76-CBO/figures/NESplotClust10FGSEAons76Cbo.png",
       units = c("in"),
       width = 25,
       height = 20,
       dpi = 300,
       o76cbo.clus10.ggplot)


# Plot cluster 6, CBO (non-malignant granule neurons) FGSEA
fgseaResTidy.clus6.cbo.GOBP <- 
  readRDS("outputs50K/CBO/cluster6.GranuleNeurons.Cbo.multifgsea.presto.GOBiolProcessGenes.rds")

# only plot the top 20 pathways
cbo.clus6.ggplot <- ggplot(fgseaResTidy.clus6.cbo.GOBP %>% filter(padj < 0.05) %>% head(n=20), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "") + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 30)) # can change bar colours using ggplot

ggsave(file = "./outputs50K/CBO/figures/NESplotClust6FGSEAcbo.png",
       units = c("in"),
       width = 25,
       height = 20,
       dpi = 300,
       cbo.clus6.ggplot)

###############################################################################################

# ONS76 cells in G1: coculture (SCT integrated). Spheroid culture not plotted as no pathway enrichment, monolayer
# also not plotted as only 1 pathway, 'GOBP_DNA_BIOSYNTHETIC_PROCESS' had a padj < 0.05

fgseaResTidy.cocult.o76.comb.GOBP <- 
  readRDS("outputs50K/INTEGRATE/ONS76onlyCocultureG1.pathwayEnrichment.G1integrate.Mono.Spher.Cocult.multifgsea.presto.GOBP.rds")

# only plot the top 20 pathways
o76.g1.cocult.ggplot <- ggplot(fgseaResTidy.cocult.o76.comb.GOBP %>% filter(padj < 0.05) %>% head(n=20), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "") + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 30)) # can change bar colours using ggplot

ggsave(file = "./outputs50K/INTEGRATE/figures_ONS76/NESplotONS76onlyCoculture.Integrated.mono.spher.cocult.FGSEA.png",
       units = c("in"),
       width = 25,
       height = 20,
       dpi = 300,
       o76.g1.cocult.ggplot)



