# BUSCO Synteny

Using [chromsyn](https://github.com/slimsuite/chromsyn/blob/main/Walkthrough.md), which investigates synteny using BUSCO-level variation between species. First we just need to run busco on all the genomes. **Currently using hymenoptera_db12**.

**Ultimate outputs:**

## 1. Number of syntenic blocks >100kb and containing >2 BUSCO genes are in each comparison

![synteny panel C](/figures/synteny_panelC.png)

___

# Code

First, ensure proper input files for chromsyn are made.

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <target>"
    exit 1
fi

#module load miniconda
#source activate chromsyn

# Prep busco db 
BUSCO_DB=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/busco_synteny/busco_downloads
LINEAGE=hymenoptera_odb12
#busco --download ${LINEAGE} --download_path ${BUSCO_DB}

target=$1 

WD=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/
GENOME="${WD}/${target}.chr.fa"
cd ${WD}/busco_synteny

# Generate inputs: seq input
python /90daydata/coffea_pangenome/fly_alignments/chr_assemblies/telociraptor/code/telociraptor.py seqin=${GENOME} basefile=${target} i=-1 tweak=F telonull=T

# busco 
busco -f -o run_${target} -i ${GENOME} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
cp -v run_${target}/run_hymenoptera_odb12/full_table.tsv ${target}.busco5.tsv
```

Assess BUSCO conservation across species:

* This will plot only BUSCO regions that are > 100kb

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/busco_synteny
cd ${WD}

> busco.fofn > sequences.fofn

for i in $(cat Samples.list); do 
    echo -e "${i} busco_synteny/${i}.busco5.tsv" >> busco.fofn
    echo -e "${i} busco_synteny/${i}.chr.telomeres.tdt" >> sequences.fofn
done 

# focus=Quadrastichus_erythrinae  minregion=-1 for all busco
Rscript ~/symlinks/software/chromsyn/chromsyn.R \
    orient=none \
    plotdir=chromsyn_n8_format basefile=chromsyn_n8_format \
    minregion=100000 \
    ygap=3 ypad=0.1 labelsize=0.75 opacity=0.8 \
    pdfheight=5 pdfwidth=6
```

## Plot

Plot alignments in R:

- Plot BUSCO completeness per genome by counting BUSCO (C+S) hits.
- Load chromsyn region data and compute pairwise BUSCO block counts between genomes (Genome × HitGenome).
- Plot log‑scaled heatmap of BUSCO block sharing.

```R
setwd('/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/')
library(tidyverse)
library(openxlsx)
library(scales)
library(RColorBrewer)

# colors
ord <- read_tsv('Samples.list',col_names = F)
ord <- ord %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", X1) )
spcols <- brewer.pal(8,'Set2')
show_col(spcols)
spcol <- spcols[c(3,5,1,8,6,2,7,4)]
show_col(spcol)
spdf <- data.frame(Species = ord$Species, col=spcol)  

f <- read.xlsx('chromsyn_n8_format_db12/chromsyn_n8_format_db12.chromsyn.xlsx',sheet='BUSCO')
f <- f %>% mutate( Genome = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Genome) ) 
f$Genome <- factor(f$Genome,levels=rev(spdf$Species))
counts <- f %>% count(Genome) 
counts
cp <- counts %>% 
  ggplot(aes(y = Genome, x = n, fill = Genome)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n,
                x = n * 1.02),   #push labels slightly to the right
            size = 1.5, hjust = 0) +
  scale_fill_manual(values = spdf$col, breaks = spdf$Species) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  xlab("BUSCO Genes (C+S)") +
  ylab("")
cp

ggsave('figures/20260105_BUSCO_hits.pdf',cp,height=1.75,width=2)

# parasitoid-specific busco? this won't really tell you much since these are highly conserved anyways 
parasitoids <- c("Q_erythrinae","N_vitripennis", "B_kinseyi")
non_parasitoids <- c("D_longicaudata", "A_mellifera", "V_crabro", "S_invicta", "N_pinetum")

busco_matrix <- f %>% select(Genome, BuscoID) %>% 
  mutate(value = 1) %>% distinct() %>% 
  tidyr::pivot_wider(names_from = Genome, values_from = value, values_fill = 0) %>% 
  dplyr::select( BuscoID, 
                 dplyr::all_of(parasitoids), 
                 dplyr::all_of(non_parasitoids) )

pdf('figures/20260102_UpSetAll.pdf',height=5,width=11)
UpSetR::upset(
  as.data.frame(busco_matrix),
  sets = c(parasitoids, non_parasitoids),
  nintersects = NA,
  text.scale=0.25,
  order.by='freq'
)
dev.off()

# Simple upset
busco_parasitoids <- f %>% filter(Genome %in% parasitoids) %>% distinct(BuscoID) %>% mutate(ID = 'Parasitoids')
busco_nonparasitoids <- f %>% filter(Genome %in% non_parasitoids) %>% distinct(BuscoID) %>% mutate(ID = 'Non_Parasitoids')
busco_simple_matrix <- rbind(busco_parasitoids,busco_nonparasitoids) %>% 
  select(ID, BuscoID) %>% 
  mutate(value = 1) %>% distinct() %>% 
  tidyr::pivot_wider(names_from = ID, values_from = value, values_fill = 0)

pdf('figures/20260102_UpSetHighLevel.pdf',height=2,width=4.5)
UpSetR::upset(
  as.data.frame(busco_simple_matrix),
  order.by='freq',
  text.scale=0.75
)
dev.off()
# n=3 parasitoid only, but highly skeptical 

# which genes
unique_parasitoid_buscos <- setdiff(busco_parasitoids$BuscoID, busco_nonparasitoids$BuscoID)
unique_parasitoid_buscos

# read in busco files to see genes
b1 <- read_tsv('busco_synteny/Apis_mellifera.busco5.tsv',skip = 2) %>% dplyr::rename('BuscoID' = `# Busco id`)
b1 %>% filter(BuscoID %in% unique_parasitoid_buscos)

# how about differences betwen parasitoids and non parasitoids?
counts$group <- ifelse(counts$Genome %in% parasitoids, "parasitoid", "non_parasitoid")
aggregate(n ~ group, counts, function(x) c(mean = mean(x), range = range(x)))
t.test(n ~ group, data = counts)

# Welch Two Sample t-test
# 
# data:  n by group
# t = 4.5899, df = 2.1177, p-value = 0.03975
# alternative hypothesis: true difference in means between group non_parasitoid and group parasitoid is not equal to 0
# 95 percent confidence interval:
#   47.63187 812.23479
# sample estimates:
#   mean in group non_parasitoid     mean in group parasitoid 
# 5794.600                     5364.667 


#Heatmap
h <- read.xlsx('chromsyn_n8_format_db12/chromsyn_n8_format_db12.chromsyn.xlsx',sheet='Regions')
h <- h %>% mutate( Genome = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Genome),
                   HitGenome = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", HitGenome) )
h2 <- h %>% group_by(Genome,HitGenome) %>% summarize(len = sum(Length), .groups='drop')
h2 <- h %>% group_by(Genome,HitGenome) %>% summarize(len = sum(BuscoID))
h2$Genome <- factor(h2$Genome, levels = rev(spdf$Species)) 
h2$HitGenome <- factor(h2$HitGenome, levels = rev(spdf$Species))

# only keep lower triangle
h2_lower <- h2 %>%
  filter(as.integer(Genome) <= as.integer(HitGenome))

# label
h2_lower <- h2_lower %>%
  dplyr::mutate(
    label = ifelse(len >= 1000,
                   paste0(round(len/1000, 1), "k"),
                   as.character(len))
  )

# Genome      HitGenome  len label
# 1  D_longicaudata   Q_erythrinae  227   227
# 2        V_crabro D_longicaudata  239   239
# 3       N_pinetum D_longicaudata  241   241
# 4  D_longicaudata  N_vitripennis  246   246
# 5       S_invicta D_longicaudata  277   277
# 6     A_mellifera D_longicaudata  345   345
# 7  D_longicaudata      B_kinseyi  357   357
# 8       N_pinetum       V_crabro  371   371
# 9       N_pinetum   Q_erythrinae  386   386
# 10      N_pinetum  N_vitripennis  394   394
# 11       V_crabro  N_vitripennis  429   429
# 12       V_crabro   Q_erythrinae  468   468
# 13      N_pinetum      S_invicta  488   488
# 14      S_invicta  N_vitripennis  492   492
# 15      S_invicta   Q_erythrinae  535   535
# 16      N_pinetum      B_kinseyi  550   550
# 17       V_crabro      B_kinseyi  560   560
# 18      N_pinetum    A_mellifera  568   568
# 19    A_mellifera   Q_erythrinae  582   582
# 20    A_mellifera  N_vitripennis  590   590
# 21      B_kinseyi  N_vitripennis  686   686
# 22      S_invicta      B_kinseyi  743   743
# 23    A_mellifera      B_kinseyi  812   812
# 24      B_kinseyi   Q_erythrinae  834   834
# 25      S_invicta       V_crabro  881   881
# 26       V_crabro    A_mellifera 1086  1.1k
# 27      S_invicta    A_mellifera 1261  1.3k
# 28  N_vitripennis   Q_erythrinae 2482  2.5k


# Make a heatmap
hm <- ggplot(h2_lower, aes(x = Genome, y = HitGenome, fill = len)) +
  geom_tile(color = "white") +
  scale_fill_continuous(low='yellow',high='red', trans='log10') +
  geom_text(aes(label = label), color = "black", size = 2) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = "",
    fill = "BUSCO Blocks\n >100 Kb",
  )
hm

ggsave('figures/20260105_HeatmapChromsynBuscoHits100kb.pdf',hm,height=3,width=4.5)

```
