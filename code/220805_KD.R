library(readxl)
library(tidyverse)
library(janitor)
library(tidylog)
library(vegan)
library(rstatix)
library(ape)
library(RColorBrewer)
library(decontam)
library(phyloseq)
library(microbiome)
library(patchwork)
library(corncob)


set.seed(1851)
permutations = 9999

# Load in feature counts
feature.count <- read_excel("data/results.xlsx", sheet = "Sample_Summary", skip = 13) 
colnames(feature.count) <-  c("sample_id", "feature.count")
  rename(featureid = "#OTU ID")

# Load in dada2 table
df <- read_excel("data/results.xlsx", sheet = "Annotated_table") %>% 
  rename(featureid = "#OTU ID")

# load metadata
metadata <- read_excel("data/results.xlsx", sheet = "Sheet1") %>% 
  select(sample, Brand, Batch, Sample_ID) %>% 
  clean_names() %>% 
  mutate(sample_type = case_when(brand == "lysis buffer" ~ "NEG",
                                 TRUE ~ "POS"))

# pull taxonony and table from df
taxonomy <- df %>% 
  select(featureid, Taxon)

taxonomy_long <- taxonomy %>% 
  separate(col = Taxon, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ") %>% 
  pivot_longer(-featureid,
               names_to = "level",
               values_to = "taxon") %>% 
  group_by(featureid) %>% 
  fill(taxon)

table <- df %>% 
  select(featureid, metadata$sample_id) %>% 
  column_to_rownames(var = "featureid")





## DECONTAM
OTU = otu_table(table, taxa_are_rows = TRUE)

metdata_ps = sample_data(metadata %>% 
                           column_to_rownames(var = "sample_id"),
                         errorIfNULL = T)

taxonomy_wide <- taxonomy_long %>% 
  pivot_wider(names_from = "level",
              values_from = "taxon") %>% 
  column_to_rownames(var = "featureid")
tax = tax_table(taxonomy_wide %>% as.matrix())

ps <- phyloseq(OTU, tax, metdata_ps)

sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "NEG"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
2460+48
keep <- contamdf.prev %>% 
  filter(contaminant == "FALSE") %>% 
  row.names() %>% 
  as.array() 

filt.table <- table %>% filter(rownames(table) %in% keep) %>% as.data.frame() %>% 
  rownames_to_column(var = "featureid") %>% 
  as.data.frame()

metadata %>%  # Get N for each group
  group_by(brand) %>% 
  tally()

table_long <- filt.table %>% 
  pivot_longer(-featureid,
               names_to = "sample_id",
               values_to = "count")

table_long %>% # get rarefaction depth, seems too low to rarefy
  group_by(sample_id) %>% 
  summarise(sum = sum(count)) %>% 
  arrange(sum) %>% 
  inner_join(., metadata, by= "sample_id")



# ALPHA STATS
brands <- c("Alcon", "Anvision", "HyVisc")

alpha.stats <- table_long %>% 
  group_by(sample_id) %>% 
  summarize(richness = richness(count, index = c("observed")) %>% as.double(),
            choa1 = richness(count, index = c("chao1")) %>% as.double(),
            simpson = vegan::diversity(count, 
                                index = "simpson"),
            shannon = vegan::diversity(count, 
                                index = "shannon",
                                base = exp(1)),
            feature.count = sum(count)) %>% 
  inner_join(., metadata, by = "sample_id") 

alpha.long <- alpha.stats %>%
  select(richness, choa1, shannon, simpson, feature.count, brand) %>% 
  pivot_longer(-brand,
               names_to = "alpha.stat",
               values_to = "value")

chao1 <- alpha.long %>% 
  filter(alpha.stat == "choa1") %>% 
  filter(brand %in% brands) 

sink("stats/chao1.txt")  

print("Shapiro Wilkes test for normality")
shapiro_test(chao1$value) # p = 0.0209

print("Kruskal Wallis Test")
print("chao1 ~ brand")
chao1 %>% 
  kruskal_test(value~brand)

sink()


chao1.plot <- alpha.long  %>% 
  filter(alpha.stat == "choa1") %>% 
  filter(brand %in% brands) %>% 
  ggplot(aes(x = brand, y = value)) +
  geom_boxplot(width = 0.5) +
  geom_point(size = 2,
             aes(color = brand)) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(y = "Chao1 Index",
       caption = "p = 0.697") +
  theme(axis.title.y = element_text(face = "bold", color = "black", size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = 'none',
        plot.caption = element_text(face = "bold", color = "black", size = 12)) 

# ggsave("plots/chao1.png",
#        width = 4,
#        height = 4,
#        units = c("in"),
#        bg = "white")


## Beta diversity

metadata_filt <- metadata %>% 
    filter(brand %in% brands)
table_filt <- table %>% 
  select(contains(metadata_filt$sample_id))

table_t <- table_filt %>% t() # Transpose table

table.transformed <- table_t^1/4  #quarter root transformation of table

dist_bc <- vegdist(table.transformed, method= "bray")  # Calculate bray-curtis distances
dist_j <- vegdist(table.transformed, method= "jaccard")  # Calculate jaccard  distances

pcoa_bc <- pcoa(dist_bc, correction = "cailliez")   # Calculate PCoA with BC distances
pcoa_j <- pcoa(dist_j, correction = "cailliez")   # Calculate PCoA with J distances

p_var_bc <- pcoa_bc$values$Eigenvalues/pcoa_bc$trace*100 # Variance
p_var_j <- pcoa_j$values$Eigenvalues/pcoa_j$trace*100 # Variance


pcoa_vectors_bc <- pcoa_bc$vectors %>% as_tibble(rownames = "sampleid") %>% 
  select(sampleid, Axis.1,Axis.2)
colnames(pcoa_vectors_bc) <- c("sample_id", "PCo1", "PCo2") # Rename values

pcoa_vectors_j <- pcoa_j$vectors %>% as_tibble(rownames = "sampleid") %>% 
  select(sampleid, Axis.1,Axis.2)
colnames(pcoa_vectors_j) <- c("sample_id", "PCo1", "PCo2") # Rename values


variance_rep_bc <- round(p_var_bc[1:2],2) ## Pull variance represented for axis lables
variance_rep_j <- round(p_var_j[1:2],2) ## Pull variance represented for axis lables

  
#Join with metadata for plotting
pcoa.metadata_bc <- inner_join(pcoa_vectors_bc, metadata, by = "sample_id")
pcoa.metadata_j <- inner_join(pcoa_vectors_j, metadata, by = "sample_id")

samples_filt <- metadata_filt %>% 
  filter(sample_id %in% attributes(dist_bc)$Labels)


sink("stats/adonis_output.txt")
print("Bray-Curtis Distances")
adonis.brand_bc <- adonis(dist_bc ~ brand, 
                             permutations = permutations, 
                             data = samples_filt)
adonis.brand_bc$aov.tab
print("Jaccard Distances")
adonis.brand_j <- adonis(dist_j ~ brand, 
                          permutations = permutations, 
                          data = samples_filt)
adonis.brand_j$aov.tab
sink()
# no sig differences in either distance

bc_pcoa <- pcoa.metadata_bc %>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(shape = brand,
                 color = brand),
             alpha = 0.9,
             size = 4) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(caption = paste0("F: ", round(adonis.brand_bc$aov.tab$F.Model, 2), 
                        ", p = ", round(adonis.brand_bc$aov.tab$`Pr(>F)`, 2)), 
       x = paste0("PCo1 - ", variance_rep_bc[1], "%"),
       y = paste0("PCo2 - ", variance_rep_bc[2], "%")) +
  theme(axis.title = element_text(face = "bold", color = "black", size = 14),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        plot.caption = element_text(face = "bold", color = "black", size = 12)) 


j_pcoa <- pcoa.metadata_j%>% 
  ggplot(aes(x = PCo1, y = PCo2)) +
  geom_point(aes(shape = brand,
                 color = brand),
             alpha = 0.9,
             size = 4) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(caption = paste0("F: ", round(adonis.brand_j$aov.tab$F.Model, 2), 
                        ", p = ", round(adonis.brand_j$aov.tab$`Pr(>F)`, 2)), 
       x = paste0("PCo1 - ", variance_rep_j[1], "%"),
       y = paste0("PCo2 - ", variance_rep_j[2], "%")) +
  theme(axis.title = element_text(face = "bold", color = "black", size = 14),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        plot.caption = element_text(face = "bold", color = "black", size = 12)) 

bc_pcoa + theme(legend.position = "none") + j_pcoa + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold",
                                size = 14))
  
ggsave("plots/bc_j_pcoa.png",
       height = 4,
       width = 8,
       units = c("in"),
       bg = "white")

# Taxonomy
fill_palette <- colorRampPalette(brewer.pal(11, "Spectral")) # expand spectral palette

table_long
taxonomy_long

rel_abund_table <- table_long %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  left_join(., taxonomy_long, by = "featureid") %>% 
  mutate(taxon = substring(taxon, 4)) %>% 
  left_join(., metadata, by = "sample_id")

rel_abund_table %>% 
  filter(sample_id == "KD001_A09")

plot_taxa_bar_plot <- function(taxa_level = "Phylum", cutoff = 0.05){
  data_to_plot <- rel_abund_table %>% 
    filter(brand %in% brands) %>% 
    filter(level == taxa_level) %>% 
    group_by(sample_id, taxon) %>% 
    summarize(rel_abund = sum(rel_abund)) %>% 
    left_join(., metadata, by = "sample_id")
  
  taxon_pool <- data_to_plot %>% 
    group_by(taxon) %>% 
    summarise(pool = max(rel_abund) < cutoff, .groups = 'drop') 
  
  data_to_plot <- inner_join(data_to_plot, taxon_pool, by = "taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) %>% 
    group_by(sample_id, brand, taxon) %>% 
    summarise(rel_abund = sum(rel_abund)) %>% 
    arrange(desc(rel_abund))
  # return(data_to_plot)
  

  taxa_list <- unique(data_to_plot$taxon)
  taxa_list <- taxa_list[!taxa_list %in% c("Other")]

  data_to_plot$taxon <- factor(data_to_plot$taxon,
                               levels = c(taxa_list, "Other"))
  # return(data_to_plot)

  fill_length <- data_to_plot %>%
    ungroup() %>%
    distinct(taxon) %>%
    pull(taxon) %>%
    length()

  data_to_plot %>% ggplot(aes(x = sample_id, y = rel_abund)) +
    geom_col(aes(fill = taxon)) +

    scale_fill_manual(values = c(fill_palette(fill_length-1), "#666666")) +
    scale_y_continuous(expand = c(0,0),
                       labels = scales::percent) +

    facet_grid(~brand,
               scales = "free",
               space = "free",
               switch = "both") +
    labs(y = "Relative Abundance") +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(face = "bold",
                                 color = "black",
                                 size = 12),

      strip.background = element_rect(fill = NA,
                                      color = NA),
      strip.placement = "outside",
      strip.text = element_text(face = "bold",
                                color = "black",
                                size = 14),

      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold",
                                  color = "black",
                                  size = 14),
      axis.text.y = element_text(face = "bold",
                                 color = "black",
                                 size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

}

f <- plot_taxa_bar_plot("Family", 0.05)
p <- plot_taxa_bar_plot("Phylum", 0.05)

(p/f ) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold",
                                size = 14),
        legend.justification = "left")

# ggsave("plots/p_f_taxa_barplot.png",
#        height = 8,
#        width = 10,
#        units = c("in"),
#        bg = "white")


kw_test <- rel_abund_table %>% 
  filter(brand %in% brands) %>% 
  group_by(featureid) %>% 
  kruskal_test(rel_abund ~ brand) 
  
kw_test %>% filter(p < 0.05) %>% 
  filter(featureid == "65d43491988bfe557da4d86a5ba25dae")

### CORE MICROBIOME
# https://microbiome.github.io/tutorials/Core.html
tmp <- ps.rel 
taxa_names(tmp)[1:5]
  
ps.rel <- microbiome::transform(ps, "compositional") 
otu_to_taxonomy(data = ps.rel, OTU = rownames(otu_table(ps.rel)) ,level = "Genus")


prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

core <- plot_core(ps.rel,
               plot.type = "heatmap", 
               colours = rev(brewer.pal(5, "Spectral")),
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(ps.rel, sort = TRUE)[100]) 

keep_taxonomy <- core$data %>%
  distinct(Taxa)

taxa_label <- taxonomy_long %>% 
  filter(level == "Genus") %>% 
  select(featureid, taxon) %>% 
  mutate(taxon = substring(taxon, 4)) %>% 
  rename(Taxa = "featureid") %>% 
  filter(Taxa %in% keep_taxonomy$Taxa)

core$data <- core$data %>% 
  inner_join(., taxa_label, by = "Taxa") %>% 
  mutate(Taxa = taxon) %>% 
  select(Taxa, DetectionThreshold, Prevalence) %>% 
  arrange(DetectionThreshold, Prevalence)

core$data %>% 
  ggplot(aes(x = DetectionThreshold, y = reorder(Taxa, Prevalence), 
             fill = Prevalence)) +
  geom_tile() +
  scale_fill_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                       guide = guide_colorbar(ticks = F), 
                       labels = scales::percent, 
                       limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(face = "bold.italic", 
                                   color = "black", size = 10),
        axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        legend.title = element_text(face = "bold", hjust = 0.5, color = "black")) +
  labs(x = "Detection Threshold\nRelative Abundance (%)")

  # ggsave("plots/core_microbiome.png",
  #        width = 7,
  #        height = 7,
  #        units = c("in"),
  #        bg = "white")

core$data %>% 
  as_tibble() %>% 
  filter(as.double(DetectionThreshold) > 0.05) %>% 
  filter(Prevalence > 0.5)
  



## SPECIFIC TAXA
  
pull_asv <- function(ASV) {
  table_long %>% 
  inner_join(., taxonomy_long) %>%
  filter(level == "Genus") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund= count/sum(count) ) %>% 
  filter(featureid == ASV) %>% 
  group_by(sample_id, taxon) %>% 
  summarize(rel_abund = sum(rel_abund)) %>% 
  inner_join(., metadata, by = "sample_id") %>% 
  filter(brand %in% brands)
}


plot_asv <-  function(table, taxa, p_val) {
  table %>% 
  ggplot(aes(x = brand, y = rel_abund)) +
  geom_boxplot(width = 0.5) +
  geom_point(size = 2,
             aes(color = brand)) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = paste(taxa, "\nRelative Abundance"),
       caption = paste0("p = ", p_val)) +
  theme(axis.title.y = element_text(face = "bold", color = "black", size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = 'none',
        plot.caption = element_text(face = "bold", color = "black", size = 12)) 
}

staph_table <- pull_asv("65d43491988bfe557da4d86a5ba25dae")
staph_plot <- plot_asv(staph_table, "Staphylococcus", 0.116)

# ggsave("plots/staph_asv.png",
#        width = 4,
#        height = 4,
#        units = c("in"),
#        bg = "white")


sink("stats/specific_taxa/staphylococcus.txt")
print("summary stats")
staph_table %>%
  group_by(brand) %>% 
  summarize(mean = mean(rel_abund) *100,
            sd = sd(rel_abund) * 100)

print("Shaprio-Wilkes")
shapiro_test(staph_table$rel_abund)
print("ANOVA: rel_abund ~ brand")
aov(staph_table$rel_abund ~ staph_table$brand) %>% 
  summary()
sink()

cuti_table <- pull_asv("5a7b179b1b45f0fe2282f260bf073f60")
cuti_plot <- plot_asv(cuti_table,"Cutibacterium", 0.783)
# ggsave("plots/cutibacterium_asv.png",
#        width = 4,
#        height = 4,
#        units = c("in"),
#        bg = "white")

sink("stats/specific_taxa/cutibacterium.txt")
print("summary stats")
cuti_table %>%
  group_by(brand) %>% 
  summarize(mean = mean(rel_abund) *100,
            sd = sd(rel_abund) * 100)

print("Shaprio-Wilkes")
shapiro_test(cuti_table$rel_abund)
print("ANOVA: rel_abund ~ brand")
aov(cuti_table$rel_abund ~ cuti_table$brand) %>% 
  summary()
print("KW: rel_abund ~ brand")
kruskal.test(cuti_table$rel_abund ~ cuti_table$brand)
sink()

flavo_table <- pull_asv("68262884f1021c896f8e1bf7348d773c")
flavo_plot <- plot_asv(flavo_table,"Flavobacterium", 0.285)

sink("stats/specific_taxa/flavobacterium.txt")
print("summary stats")
flavo_table %>%
  group_by(brand) %>% 
  summarize(mean = mean(rel_abund) *100,
            sd = sd(rel_abund) * 100)

print("Shaprio-Wilkes")
shapiro_test(flavo_table$rel_abund)
print("ANOVA: rel_abund ~ brand")
aov(flavo_table$rel_abund ~ flavo_table$brand) %>% 
  summary()
print("KW: rel_abund ~ brand")
kruskal.test(flavo_table$rel_abund ~ flavo_table$brand)
sink()


lacto_table <- pull_asv("95fdd816723ca482a5caba10bea171c8")
lacto_plot <- plot_asv(flavo_table,"Lactobacillus", 0.478)

sink("stats/specific_taxa/lactobacillus.txt")
print("summary stats")
lacto_table %>%
  group_by(brand) %>% 
  summarize(mean = mean(rel_abund) *100,
            sd = sd(rel_abund) * 100)

print("Shaprio-Wilkes")
shapiro_test(lacto_table$rel_abund)
print("ANOVA: rel_abund ~ brand")
aov(lacto_table$rel_abund ~ lacto_table$brand) %>% 
  summary()
print("KW: rel_abund ~ brand")
kruskal.test(lacto_table$rel_abund ~ lacto_table$brand)
sink()


staph_plot + cuti_plot + flavo_plot + lacto_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold",
                                size = 14))
ggsave("plots/dominant_bacteria.png",
       width = 7,
       height = 7,
       units = c("in"),
       bg = "white")



## FEATURE COUNTS
feature.count.plot <- feature.count %>% 
  inner_join(., metadata, by = "sample_id") %>% 
  filter(brand %in% brands) %>% 
  ggplot(aes(x = brand, y = feature.count)) +
  geom_boxplot(width = 0.5) +
  geom_point(size = 2,
             aes(color = brand)) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_log10() +
  labs(y = "Log10 Feature Count",
       caption = "p = 0.179") +
  theme(axis.title.y = element_text(face = "bold", color = "black", size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = 'none',
        plot.caption = element_text(face = "bold", color = "black", size = 12)) 

ggsave("plots/feature.count.png",
       width = 4,
       height = 4,
       units = c("in"),
       bg = "white")


feature.count %>% 
  inner_join(., metadata, by = "sample_id") %>% 
  filter(brand %in% brands) %>% 
  shapiro_test(feature.count) # p < 0.001

sink(file = "stats/feature.count.kw.txt")

print("Shapiro-Wilkes test for normality")

feature.count %>% 
  inner_join(., metadata, by = "sample_id") %>% 
  filter(brand %in% brands) %>% 
  shapiro_test(feature.count) # p < 0.001

print("Kruskal-Wallis Test")

print("feature.count ~ brand")

feature.count %>% 
  inner_join(., metadata, by = "sample_id") %>% 
  kruskal_test(feature.count ~ brand) 

sink()


feature.count.plot + chao1.plot +plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold",
                                size = 14))
ggsave("plots/feature.count_chao1.png",
       width = 7,
       height = 4,
       units = c("in"),
       bg = "white")

