# Repeat-level variation

Repeats are characterized with [earlgrey v6.3.5](https://github.com/TobyBaril/EarlGrey/), using the hymenopteran dfam database. 

Ultimate outputs: 

1. High‑level repeat landscape across species:

![highlevel panel A](figures/highlevel_panelA.png)

2. Excluding B. kinseyi:

![highlevel panel B](figures/highlevel_panelB.png)

3. Genome size ~ repeat variation:

![divergence panel AB](figures/divergence_panelAB.png)

4. Divergence among TE DNA and LTR subclasses:

![divergence panel C](figures/divergence_panelC.png)

4. Model results: divergence across species LMM:

![divergence panel D](figures/divergence_panelD.png)



Configuration:

```
Ensure you have downloaded hymenoptera dfam database here:
/project/coffea_pangenome/Software/Merondun/.conda/envs/eg6/share/RepeatMasker/Libraries/famdb

library_path=$(which RepeatMasker | sed 's|/bin/RepeatMasker|/share/RepeatMasker/Libraries/famdb|g')
cd $(which RepeatMasker | sed 's|/bin/RepeatMasker|/share/RepeatMasker/|g')

# Then, add those partitions to repeatmasker:
perl ./configure -libdir $(which RepeatMasker | sed 's|/bin/RepeatMasker|/share/RepeatMasker/Libraries/|g') \
    		-trf_prgm $(which RepeatMasker | sed 's|/bin/RepeatMasker|/bin/trf|g') \
   			-rmblast_dir $(which RepeatMasker | sed 's|/bin/RepeatMasker|/bin|g') \
    		-hmmer_dir $(which RepeatMasker | sed 's|/bin/RepeatMasker|/bin|g') \
    		-abblast_dir $(which RepeatMasker | sed 's|/bin/RepeatMasker|/bin|g') \
    		-crossmatch_dir $(which RepeatMasker | sed 's|/bin/RepeatMasker|/bin|g') \
    		-default_search_engine rmblast
    		
```

Run earlgrey:

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=112Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

#module load miniconda
#source activate earlgrey
#source activate eg6

SAMPLE=$1
WD=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies

mkdir -p ${WD}/repeats
cd ${WD}/repeats

FASTA="${WD}/${SAMPLE}.chr.fa" 

echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ~~~~\e[0m"
# Run earlgrey with eudicotyledons repeatmasker search time, generating soft-masked genome and run helitrons. 
earlGrey -r hymenoptera -d yes -e yes -t ${t} -g ${FASTA} -s ${SAMPLE} -o ${WD}/repeats
```

Copy the useful outputs:

```bash
for i in $(cat ../Samples.list); do echo $i; cp ${i}_EarlGrey/${i}_summaryFiles/${i}.highLevelCount.kable ${i}_EarlGrey/${i}_summaryFiles/${i}.*gff results/; done
```

And also the softmasked genome output, in case we want to run cactus later:

```bash
for i in $(cat ../Samples.list); do echo $i; cp ${i}_EarlGrey/${i}_summaryFiles/${i}.softmask* ../softmasked/; done
```

Summarize output gffs:

* Loads repeat GFFs from earlgrey standardizes TE classifications.

* Joins chromosome sizes to compute genome‑wide repeat proportions.

* Merge overlaps and assign each base to a single TE class using a hierarchy.

* Calculates TE coverage, copy number, family counts, and non‑repeat fraction per species.

* Exports two files: a TE summary table and a coordinate‑level file with Kimura divergence.

```R
setwd('/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/repeats/results')
library(tidyverse)
library(ggpubr)
library(GenomicRanges)
library(scales)

#### File for assigning TE hierarchical ordering etc ####
te_order <-c("LTR","LINE","SINE","DNA","Penelope","Rolling Circle",
             "Other (Simple Repeat, Microsatellite, RNA)","Unclassified","Non-Repeat")

##### Import gff #####
files <- list.files('.', pattern = "gff", full.names = TRUE)

# import and parse names
gff <- map_dfr(files, ~ {
  tib <- read_tsv(.x, col_types = cols(),col_names = F) %>% 
    mutate(filename = basename(.x)) %>% 
    separate(filename, into = c("Species","del"),sep='\\.', remove = FALSE) %>% 
    #mutate(Species = paste(Genus, Species)) %>% 
    select(-filename,-del,-X2,-X6,-X8)
})
gff
names(gff) <- c('chr','type','start','end','Strand','ID','Species')

# add chr lengths
chrs <- read_tsv('Chr_lengths.txt',col_names = F)
names(chrs) <- c('chr','ChrSize','Species')
chrs <- chrs %>% 
  group_by(Species) %>% 
  summarize(GSize = sum(ChrSize)) %>% 
  ungroup %>% 
  mutate(phylo_order = as.integer(factor(Species, levels = unique(Species))))

# classify + NAME
gff2 <- gff %>%
  mutate(
    NAME = str_match(ID, "(?<=ID=)[^;]+")[,1],
    Classification = case_when(
      str_detect(type, "^LTR") ~ "LTR",
      str_detect(type, regex("(^|/)Penelope", ignore_case = TRUE)) ~ "Penelope",
      str_detect(type, "^LINE") ~ "LINE",
      str_detect(type, "^SINE") ~ "SINE",
      str_detect(type, "^DNA") ~ "DNA",
      str_detect(type, "Helitron|Rolling") ~ "Rolling Circle",
      type %in% c("Satellite","Simple_repeat","Low_complexity","rRNA","tRNA","snRNA","srpRNA","Microsatellite","Tandem","tandem") ~ "Other (Simple Repeat, Microsatellite, RNA)",
      type %in% c("Unknown","Unclassified","noCat") ~ "Unclassified",
      TRUE ~ "Unclassified"
    ),
    len = end - start + 1
  ) %>% left_join(chrs)

### Extract divergence ###
divdat <- gff2 %>% 
  filter(!grepl('Non|Other|Unclass',Classification)) %>% 
  mutate(div = gsub('.*kimura80=','',ID)) %>% 
  select(chr,start,end,Strand,Classification,div,Species,GSize,phylo_order) %>%
  arrange(phylo_order)

### Extract prop
# priority map for overlapping regions: assign an overlapping genomic region only ONE repeat class so that our proportion doesn't go above 100% 
prio <- c("LTR","LINE","Penelope","SINE","DNA","Rolling Circle",
          "Other (Simple Repeat, Microsatellite, RNA)","Unclassified")
prio_map <- setNames(seq_along(prio), prio)

# minimal input cols: Species, Chromosome, start, end, Classification
g <- gff2 %>%
  transmute(seqname = paste(Species, chr, sep = "|"),
            start, end, Classification)

gr_all <- GRanges(seqnames = g$seqname,
                  ranges   = IRanges(g$start, g$end),
                  Classification = g$Classification,
                  pr = unname(prio_map[g$Classification]))

# shrink: reduce per class (huge speedup cmpard to genomic ranges on full )
gr_by_class <- split(gr_all, mcols(gr_all)$Classification)
gr_reduced <- do.call(c, lapply(names(gr_by_class), function(cl) {
  if (length(gr_by_class[[cl]]) == 0) return(GRanges())
  rr <- reduce(gr_by_class[[cl]], ignore.strand = TRUE)
  mcols(rr)$Classification <- cl
  mcols(rr)$pr <- prio_map[[cl]]
  rr
}))

# disjoin once across both chrs 
dj <- disjoin(gr_reduced, with.revmap = TRUE)

# fast class assignment by min priority (vectorized)
rev <- mcols(dj)$revmap
prv <- mcols(gr_reduced)$pr
min_pr <- vapply(rev, function(ix) min(prv[as.integer(ix)]), integer(1))
cls <- names(prio_map)[min_pr]

# now just count 
cov_by_class <- tibble(
  seqname = as.character(seqnames(dj)),
  Classification = cls,
  Coverage = as.numeric(width(dj))
) %>%
  tidyr::separate(seqname, into = c("Species","chr"), sep = "\\|") %>%
  group_by(Species, Classification) %>%
  summarise(Coverage = sum(Coverage), .groups = "drop")

# add copies/family counts and chromosome sizes (ord)
counts <- gff2 %>%
  mutate(NAME = stringr::str_match(ID, "(?<=Name=)[^;]+")[,1]) %>%
  group_by(Species, Classification) %>%
  summarise(Copies = n(), Family_Count = n_distinct(NAME), .groups = "drop")

# estimate proportion by and grab non-repeat
hl_rep <- cov_by_class %>%
  left_join(counts, by = c("Species","Classification")) %>%
  left_join(chrs) %>%
  mutate(Proportion = 100 * Coverage / GSize)

# non-repeat rows
nonrep <- hl_rep %>%
  group_by(Species, GSize, phylo_order) %>%
  summarise(Covered = sum(Coverage, na.rm = TRUE), .groups = "drop") %>%
  transmute(Classification = "Non-Repeat",
            Coverage = pmax(GSize - Covered, 0),
            Copies = NA_real_,
            Proportion = 100 * Coverage / GSize,
            GSize, Family_Count = NA_real_,
            Species, phylo_order)

# final combined tibble
hl <- bind_rows(hl_rep, nonrep) %>%
  select(Classification, Coverage, Copies, Proportion, GSize, Family_Count,
         Species, phylo_order) %>% 
  arrange(phylo_order) %>% 
  mutate(Classification = gsub(',','',Classification))

hl %>% group_by(Classification,Species) %>% summarize(Copies = sum(Copies), Family_Count = sum(Family_Count))

# Save summary file !
write.csv(hl,file='Summary_GFF_20251229.csv',quote=F,row.names=F)

# Save coordinate file! 
simp <- gff2 %>% 
  mutate(KimuraDiv = as.numeric(ifelse(grepl('Non|Other|Unclass',Classification),NA,gsub('.*kimura80=','',ID)))) %>% 
  select(chr,start,end,Strand,Species,phylo_order,GSize,Classification,KimuraDiv) %>% 
  mutate(Classification = gsub(',','',Classification),) %>% arrange(phylo_order) 

# Save divergence file! 
write.csv(simp,file='Repeats_Coordinates_Divergence_20251229.csv',quote=F,row.names=F)
```

## Plot

This script will summarize and plot the files above:

* Plot high‑level repeat landscape across species (Proportion, Coverage, Family_Count), including full dataset and a version excluding B. kinseyi.

* Compute total repeat coverage per genome (excluding Non‑Repeat) and run Spearman correlations between repeat bp and genome size.

* Fit PGLS models (λ estimated by ML) for log(GSize) ~ log(Repeats), both with full tree and with B. kinseyi removed; extract fitted values and λ profiles.

* Plot genome size vs repeat content with LM trendline, species‑specific colors, and annotated correlation statistics.

* Generate divergence distributions (Kimura) for LTR and DNA elements using ridge densities + boxplots across species.

* Fit mixed‑effects model for log‑transformed divergence (Species fixed, TE class random), compute EMMs, and test whether Q. erythrinae differs from others.

```R
setwd('/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/repeats/results')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(vegan)
library(broom)
library(scales)

##### High Level #####
high_level <- read_csv('Summary_GFF_20251229.csv') %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) )

ord <- read_tsv('../../Samples.list',col_names = F) %>% dplyr::rename(Species=X1) %>%  mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) )
high_level$Species <- factor(high_level$Species,levels=rev(ord$Species))
high_level$Classification <- factor(high_level$Classification,levels=c('Non-Repeat','Unclassified','Other (Simple Repeat Microsatellite RNA)','DNA','Penelope','Rolling Circle','LTR','LINE'))
cols <- high_level %>% dplyr::select(Classification) %>% distinct %>% mutate(Color = brewer.pal(8,'Paired'))

# Plot landscape 
all <- high_level %>% 
  mutate(Coverage = Coverage / 1e6) %>% 
  pivot_longer(c(Proportion,Coverage,Family_Count)) %>%
  filter( !(name == 'Family_Counts' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat Microsatellite RNA)")) ) %>% 
  ggplot(aes(y=Species,x=value,fill=Classification))+
  geom_bar(stat='identity',position=position_stack())+
  theme_bw(base_size=8)+
  facet_grid(.~name,scales='free',space='free_y')+
  scale_fill_manual(values=cols$Color,breaks=cols$Classification)+
  theme(strip.text.y = element_text(angle = 0))+
  ylab('')+xlab('Distinct Classifications')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
all

ggsave('../../figures/20251229_RepeatsHighLevelSummary.pdf',
       all,dpi=300,height=4,width=8)


# Excluding B. kinseyi 
excl <- high_level %>% 
  mutate(Coverage = Coverage / 1e6) %>% 
  filter(Species != 'B_kinseyi') %>% 
  pivot_longer(c(Proportion,Coverage,Family_Count)) %>%
  filter( !(name == 'Family_Counts' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat Microsatellite RNA)")) ) %>% 
  ggplot(aes(y=Species,x=value,fill=Classification))+
  geom_bar(stat='identity',position=position_stack())+
  theme_bw(base_size=8)+
  facet_grid(.~name,scales='free',space='free_y')+
  scale_fill_manual(values=cols$Color,breaks=cols$Classification)+
  theme(strip.text.y = element_text(angle = 0))+
  ylab('')+xlab('Distinct Classifications')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
excl

ggsave('../../figures/20251229_RepeatsHighLevelSummary_No_Belonocnema.pdf',
       excl,dpi=300,height=4,width=8)

#### PGLS ####
# Significance genome size ~ LTRs 
rep_df <- high_level %>% filter(!grepl('Non-Repeat',Classification)) %>% group_by(Species,GSize,phylo_order) %>% summarize(Repeats = sum(Coverage))
rep_df

# Spearman correlation
cor_res <- cor.test(
  ~ Repeats + GSize,
  data = rep_df,
  method = "spearman"
) %>% tidy()
cor_res
# A tibble: 1 × 5
# estimate statistic p.value method                          alternative
# <dbl>     <dbl>   <dbl> <chr>                           <chr>      
#   1    0.929      6.00 0.00223 Spearman's rank correlation rho two.sided  

# pgls
library(caper)
pgls_in <- rep_df %>% dplyr::select(Species,Repeats,GSize) %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) ) %>% as.data.frame
nwk <- read.tree('../../cactus/cactus_rooted_tree_full.nwk')
comp <- comparative.data(phy = nwk, data = pgls_in, names.col = Species, vcv = TRUE, na.omit = FALSE )

pgls_model <- pgls( log(GSize) ~ log(Repeats) , data = comp, lambda = "ML" ) 
summary(pgls_model)

pgls_line <- data.frame(
  Repeats = pgls_in$Repeats,
  fitted = fitted(pgls_model)
)

profile <- pgls.profile(pgls_model) 
plot(profile)

# Call:
#   pgls(formula = log(GSize) ~ log(Repeats), data = comp, lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.34596 -0.20404 -0.07719  0.25265  0.35830 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.000
# lower bound : 0.000, p = 1    
# upper bound : 1.000, p = 0.077322
# 95.0% CI   : (NA, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  10.444096   1.494202  6.9897 0.0004269 ***
#   log(Repeats)  0.496757   0.081409  6.1020 0.0008827 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2938 on 6 degrees of freedom
# Multiple R-squared: 0.8612,	Adjusted R-squared: 0.8381 
# F-statistic: 37.23 on 1 and 6 DF,  p-value: 0.0008827 

# Format annotation text
label_text <- sprintf(
  "Spearman rho = %.3f\np = %.2e",
  cor_res$estimate,
  cor_res$p.value
)

cp <- pgls_in %>% mutate(across(c(Repeats, GSize), ~ .x / 1e9))

# Axis limits for positioning
x_min <- min(cp$GSize, na.rm = TRUE)
y_max <- max(cp$Repeats, na.rm = TRUE)

# Plot with lm line and annotation
spcols <- brewer.pal(8,'Set2')
show_col(spcols)
spcol <- spcols[c(3,5,1,8,6,2,7,4)]
show_col(spcol)
spdf <- data.frame(Species=ord$Species,col=spcol) %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) ) 
library(ggrepel)
rp <- ggplot(cp, aes(x = GSize, y = Repeats, col=Species, group=NA,label=Species)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_point(size=2) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values=spdf$col,breaks=spdf$Species)+
  annotate("text", x = x_min, y = y_max, label = label_text,
           hjust = 0, vjust = 1, size = 3.5) +
  ylab('Repeat Coverage')+
  theme_bw(base_size=8)+
  theme(legend.position='none')
rp
ggsave('../../figures/20251229_RepeatsHighLevelSummary-Repeats-GenomeSize.pdf',
       rp,dpi=300,height=4,width=4)

#### PGLS wihtout B kinseyi ####
# Significance genome size ~ LTRs 
rep_df <- high_level %>% filter(!grepl('Non-Repeat',Classification) & Species != 'B_kinseyi') %>% group_by(Species,GSize,phylo_order) %>% summarize(Repeats = sum(Coverage))
rep_df

# Spearman correlation
cor_res <- cor.test(
  ~ Repeats + GSize,
  data = rep_df,
  method = "spearman"
) %>% tidy()
cor_res
# # A tibble: 1 × 5
# estimate statistic p.value method                          alternative
# <dbl>     <dbl>   <dbl> <chr>                           <chr>      
#   1    0.893      6.00  0.0123 Spearman's rank correlation rho two.sided  
  
# pgls
pgls_in <- rep_df %>% dplyr::select(Species,Repeats,GSize) %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) ) %>% as.data.frame
nwk <- read.tree('../../cactus/cactus_rooted_tree_full.nwk')
nwk2 <- drop.tip(nwk,'B_kinseyi')
plot(nwk2)
comp <- comparative.data(phy = nwk2, data = pgls_in, names.col = Species, vcv = TRUE, na.omit = FALSE )

pgls_model <- pgls( log(GSize) ~ log(Repeats) , data = comp, lambda = "ML" ) 
summary(pgls_model)

# Call:
#   pgls(formula = log(GSize) ~ log(Repeats), data = comp, lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32712 -0.06551  0.02450  0.10344  0.22675 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.000
# lower bound : 0.000, p = 1    
# upper bound : 1.000, p = 0.15747
# 95.0% CI   : (NA, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  13.585950   1.536706  8.8410 0.0003076 ***
#   log(Repeats)  0.318951   0.085522  3.7294 0.0135788 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2028 on 5 degrees of freedom
# Multiple R-squared: 0.7356,	Adjusted R-squared: 0.6827 
# F-statistic: 13.91 on 1 and 5 DF,  p-value: 0.01358 

# Format annotation text
label_text <- sprintf(
  "Spearman rho = %.3f\np = %.2e",
  cor_res$estimate,
  cor_res$p.value
)

cp <- pgls_in %>% mutate(across(c(Repeats, GSize), ~ .x / 1e9))

# Axis limits for positioning
x_min <- min(cp$GSize, na.rm = TRUE)
y_max <- max(cp$Repeats, na.rm = TRUE)

# Plot with lm line and annotation
spcols <- brewer.pal(8,'Set2')
show_col(spcols)
spcol <- spcols[c(3,5,8,6,2,7,4)]
show_col(spcol)
spdf <- data.frame(Species=ord$Species,col=spcol) %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) ) 
library(ggrepel)
rp <- ggplot(cp, aes(x = GSize, y = Repeats, col=Species, group=NA,label=Species)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_point(size=2) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values=spdf$col,breaks=spdf$Species)+
  annotate("text", x = x_min, y = y_max, label = label_text,
           hjust = 0, vjust = 1, size = 3.5) +
  ylab('Repeat Coverage')+
  theme_bw(base_size=8)+
  theme(legend.position='none')
rp
ggsave('../../figures/20251229_RepeatsHighLevelSummary-Repeats-GenomeSize-NoBkinseyi.pdf',
       rp,dpi=300,height=2.75,width=3.5)


###### Divergence Summaries ######
library(ggridges)
t <- read_csv('Repeats_Coordinates_Divergence_20251229.csv')
t <- t %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", Species) )
t$Species <- factor(t$Species,levels=rev(ord$Species))

hists <- t %>% 
  filter(grepl('LTR|DNA',Classification)) %>% 
  pivot_longer(Classification) %>% 
  ggplot(aes(x = KimuraDiv, y = Species, fill = Species)) +
  # half‑width density curves, nudged slightly up
  geom_density_ridges(
    aes(x = KimuraDiv, y = Species),
    rel_min_height = 0.01,
    scale = 0.4,
    alpha = 0.6,
    color = NA,
    position = position_nudge(y = 0.1)
  )+ 
  # thin boxplots, nudged slightly down
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
               position = position_nudge(y = -0.1))+
  scale_fill_manual(values=spdf$col,breaks=spdf$Species)+
  theme_bw(base_size = 8) +
  coord_cartesian(xlim = c(0,0.9))+
  facet_grid(.~value,scales='free')+
  theme(legend.position = "none") +
  ylab('')
hists

ggsave('../../figures/20251229_RepeatsDivergence.pdf',
       hists,dpi=300,height=3,width=3.5)

# is Q ery significantly lower than others?
library(lme4) 
library(emmeans)

# log xform
lmm_df <- t %>% 
  filter(grepl('LTR|LINE|DNA',Classification)) %>% 
  mutate(KD_log = log1p(KimuraDiv)) 

# Mixed model: Species as fixed, TE family as random (if available) 
m1 <- lmer(KD_log ~ Species + (1 | Classification), data = lmm_df)
summary(m1)

emm <- emmeans(m1, ~ Species)
emm

contrast(emm, method = "trt.vs.ctrl", ref = "Q_erythrinae", adjust = "BH")

plot(emm)
emm_df <- as.data.frame(emm) 

emms <- ggplot(emm_df, aes(y = Species, x = emmean,col=Species)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.2) + 
  theme_bw(base_size=8) + 
  scale_color_manual(values=spdf$col,breaks=spdf$Species)+
  xlab("Estimated marginal mean Kimura divergence (log scale)")+
  theme(legend.position='none') + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
emms
ggsave('../../figures/20251229_RepeatsDivergence-EMMs.pdf',
       emms,dpi=300,height=2,width=2)
```



# BUSCO Synteny

Using [chromsyn](https://github.com/slimsuite/chromsyn/blob/main/Walkthrough.md), which investigates synteny using BUSCO-level variation between species. First we just need to run busco on all the genomes.

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
LINEAGE=hymenoptera_odb10
#busco --download ${LINEAGE} --download_path ${BUSCO_DB}

target=$1 

WD=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/
GENOME="${WD}/${target}.chr.fa"
cd ${WD}/busco_synteny

# Generate inputs: seq input
python /90daydata/coffea_pangenome/fly_alignments/chr_assemblies/telociraptor/code/telociraptor.py seqin=${GENOME} basefile=${target} i=-1 tweak=F telonull=T

# busco 
busco -f -o run_${target} -i ${GENOME} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
cp -v run_${target}/run_hymenoptera_odb10/full_table.tsv ${target}.busco5.tsv
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

> busco.fofn > gaps.fofn > sequences.fofn > tidk.fofn

for i in $(grep -v 'Belonocnema' Samples.list); do 
    echo -e "${i} busco_synteny/${i}.busco5.tsv" >> busco.fofn
    echo -e "${i} busco_synteny/${i}.chr.telomeres.tdt" >> sequences.fofn
done 

# focus=Quadrastichus_erythrinae  minregion=-1 for all busco
Rscript ~/symlinks/software/chromsyn/chromsyn.R \
    orient=none \
    plotdir=chromsyn_n8_format basefile=chromsyn_n8_format \
    minregion=100000 \
    ygap=3 ypad=0.1 labelsize=0.75 opacity=0.8 \
    pdfscale=0.3
```

Plot alignments in R:

- Plot BUSCO completeness per genome by counting BUSCO (C+S) hits.
- Load chromsyn region data and compute pairwise BUSCO block counts between genomes (Genome × HitGenome).
- Plot log‑scaled heatmap of BUSCO block sharing.

```R
setwd('/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/chromsyn_n8_format/')
library(tidyverse)
library(openxlsx)
library(scales)
library(RColorBrewer)

# colors
ord <- read_tsv('../Samples.list',col_names = F)
ord <- ord %>% mutate( Species = gsub("^([A-Za-z])[A-Za-z]*_(.*)$", "\\1_\\2", X1) )
spcols <- brewer.pal(8,'Set2')
show_col(spcols)
spcol <- spcols[c(3,5,1,8,6,2,7,4)]
show_col(spcol)
spdf <- data.frame(Species = ord$Species, col=spcol)  

f <- read.xlsx('chromsyn_n8_format.chromsyn.xlsx',sheet='BUSCO')
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

ggsave('../figures/20260102_BUSCO_hits.pdf',cp,height=1.75,width=2)

# parasitoid-specific busco? this won't really tell you much since these are highly conserved anyways 
parasitoids <- c("Q_erythrinae","N_vitripennis", "B_kinseyi", "D_longicaudata")
non_parasitoids <- c("A_mellifera", "V_crabro", "S_invicta", "N_pinetum")

busco_matrix <- f %>% select(Genome, BuscoID) %>% 
  mutate(value = 1) %>% distinct() %>% 
  tidyr::pivot_wider(names_from = Genome, values_from = value, values_fill = 0) %>% 
  dplyr::select( BuscoID, 
                 dplyr::all_of(parasitoids), 
                 dplyr::all_of(non_parasitoids) )

pdf('../figures/20260102_UpSetAll.pdf',height=5,width=11)
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

pdf('../figures/20260102_UpSetHighLevel.pdf',height=2,width=4.5)
UpSetR::upset(
  as.data.frame(busco_simple_matrix),
  order.by='freq',
  text.scale=0.75
)
dev.off()
# n=10 parasitoid only, but highly skeptical 

# which genes
unique_parasitoid_buscos <- setdiff(busco_parasitoids$BuscoID, busco_nonparasitoids$BuscoID)
unique_parasitoid_buscos

# read in busco files to see genes
b1 <- read_tsv('../busco_synteny/Apis_mellifera.busco5.tsv',skip = 2) %>% dplyr::rename('BuscoID' = `# Busco id`)
b1 %>% filter(BuscoID %in% unique_parasitoid_buscos)

#Heatmap
h <- read.xlsx('../chromsyn_n8_format.chromsyn.xlsx',sheet='Regions')
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

ggsave('../figures/20260102_HeatmapChromsynBuscoHits100kb.pdf',hm,height=3,width=4.5)

```

