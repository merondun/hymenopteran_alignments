# Gene-based Alignments

Instead of examining BUSCO-level variation, this will examine whole-genic variation between genomes. Examine syntenic gene blocks where there are at least 10 genes in each block after filtering for reciprocal best hits. 

**The ultimate outputs will be:**

## 1. Synteny blocks among genomes using all genes: reciprocal best hits

![genic panel A](/figures/genic_panelA.png)

## 2. Dot plots of those reciprocal best hits:

![genic panel B and C](/figures/genic_panelBC.png)

___

## Prep

First, rename the `gff` chromsomes so that they are in a sensible 'chr' format. This will also exclude any non-chr scaffolds. 

```bash
for f in *.gff; do
  awk 'BEGIN{OFS="\t"}
       NR==FNR {map[$1]=$2; next}
       /^#/ {print; next}
       {
         if ($1 in map) {
           $1 = map[$1]
           print
         }
       }' chr_map.tsv "$f" > "${f%.gff}.renamed.gff"
done
```

Using large map file:

```bash
NC_087225.1	chr1
NC_087226.1	chr2
NC_087227.1	chr3
NC_087228.1	chr4
...
```

Note that for Nasonia sp., for some reason they put 'Curated Genomic psuedogene' as field 2-3 in a subset of ~30 genes which will cause everything to fail because the gff fields are shifted:

```
awk '$8=="-"' Nasonia_vitripennis.gff
chr1    Curated Genomic pseudogene      1781380 1782914 .       -       .       ID=gene-Or4PSE;Dbxref=GeneID:100463191,NASONIABASE:NV22749;Name=Or4PSE;description=odorant      receptor        4       pseudogene;gbkey=Gene;gene=Or4PSE;gene_biotype=pseudogene;gene_synonym=NvOr4PSE;pseudo=true
chr1    Curated Genomic pseudogene      2730126 2731601 .       -       .       ID=gene-Or235PSE;Dbxref=GeneID:100463271;Name=Or235PSE;description=odorant      receptor        235     pseudogene;gbkey=Gene;gene=Or235PSE;gene_biotype=pseudogene;gene_synonym=NvOr235PSE;pseudo=true
chr1    Curated Genomic pseudogene      2731955 2733487 .       -       .       ID=gene-Or234PSE;Dbxref=GeneID:100463270;Name=Or234PSE;description=odorant      receptor        234     pseudogene;gbkey=Gene;gene=Or234PSE;gene_biotype=pseudogene;gene_synonym=NvOr234PSE;pseudo=true
chr1    Curated Genomic pseudogene      2740382 2742013 .       -       .       ID=gene-Or231PSE;Dbxref=GeneID:100463269,NASONIABASE:NV22738;Name=Or231PSE;description=odorant  receptor        231     pseudogene;gbkey=Gene;gene=Or231PSE;gene_biotype=pseudogene;gene_synonym=NvOr231PSE;pseudo=true
```

Implement this insanity to fix, for some reason I really struggled to fix it with sed. 

```bash
awk 'BEGIN{OFS="\t"}
     /^#/ {print; next}
     {
       if ($2=="Curated" && $3=="Genomic") {
         $2="Curated_Genomic"
         #shift pseudogene into column 3
         $3=$4
         $4=$5; $5=$6; $6=$7; $7=$8; $8=$9
         # rebuild attributes from field 10 onward
         attr=""
         for(i=10;i<=NF;i++){attr=attr $i (i<NF?OFS:"")}
         $9=attr
         NF=9
       }
       print
     }' Nasonia_vitripennis.renamed.gff > fixed.gff

```

Moving on...

Extract cds and essential files for single longest transcript per gene:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/90daydata/coffea_pangenome/fly_alignments/gffs
mkdir -p ${WD}/only_longest_transcript_per_gene
cd ${WD}

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample>"
    exit 1
fi
# module load miniconda
# source activate isoseq_ann

SAMPLE=$1

TARGET=/90daydata/coffea_pangenome/fly_alignments/chr_assemblies/${SAMPLE}.chr.fa

agat_sp_keep_longest_isoform.pl --gff renamed/${SAMPLE}.renamed.gff -o only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf
gffread only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf -o only_longest_transcript_per_gene/${SAMPLE}.gff3 --keep-genes -O -g ${TARGET} -w only_longest_transcript_per_gene/${SAMPLE}.fa
TransDecoder.LongOrfs -t only_longest_transcript_per_gene/${SAMPLE}.fa --output_dir only_longest_transcript_per_gene/
TransDecoder.Predict -t only_longest_transcript_per_gene/${SAMPLE}.fa --output_dir only_longest_transcript_per_gene/ --single_best_only

```

Submit serial via sample:

```bash
cat Samples.list | xargs -I {} sbatch -J ann_{} 01_Extra_Single_Transcript.sh {} 
```

Compress and prepare files for jcvi alignment:

```bash
for i in $(cat ../Samples.list); do 
	gzip ${i}*cds -c > ../input_jcvi/${i}.cds.fa.gz; 
	gzip ${i}.gff3 -c > ../input_jcvi/${i}.gff3.gz; 	
done
```

Ensure the fasta headers from the genes match the 4th col of the .bed:

```bash
sed -i 's/.p1.*//g' *.cds;

for f in *.gff3; do
  awk 'BEGIN{OFS="\t"} {$3 = ($3=="RNA" ? "mRNA" : $3); print}' "$f" > "${f%.gff3}.fixRNA.gff3"
done
```

Convert to bed with jcvi:

```bash
# source activate jcvi
for i in $(cat  Samples.list); do 
	echo "Working on ${i}"
	python -m jcvi.formats.gff bed --type=mRNA --key=ID ${i}.fixRNA.gff3 -o ${i}.bed
	python -m jcvi.formats.fasta format ${i}.cds.fa.gz ${i}.cds
done 
```

## Alignments

Synteny with last. Retain only RBH, and only plot blocks with >15 genes (minspan)

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

reference=$1  #reference
target=$2 #target
echo "Working on $target aligned to $reference"

WD=/90daydata/coffea_pangenome/fly_alignments/gffs/input_jcvi/
cd ${WD}

rm ${reference}*last* ${reference}*lifted* ${reference}*pdf

# Align and extract. Note the formatting must be EXACT with file names, cscore = RBH according to jcvi documentation
python -m jcvi.compara.catalog ortholog ${reference} ${target} --cscore=.99 --no_strip_names 

# for dotplots
magick -density 300 ${reference}.${target}.pdf ../output_jcvi/${reference}_${target}.png

# for filtering anchors 
python -m jcvi.compara.synteny screen --minsize=10 --simple ${reference}.${target}.anchors ${reference}.${target}.anchors.minsize10

```

## Plot 

Create a layout file:

```
#       y,      xstart, xend,   rotation,       color,  label,  va,     bed
        0.89,   0.2,    0.95,   0,      ,       Quadrastichus_erythrinae,       top,    Quadrastichus_erythrinae.bed
        0.78,   0.2,    0.95,   0,      ,       Nasonia_vitripennis,    top,    Nasonia_vitripennis.bed
        0.67,   0.2,    0.95,   0,      ,       Belonocnema_kinseyi,    top,    Belonocnema_kinseyi.bed
        0.56,   0.2,    0.95,   0,      ,       Diachasmimorpha_longicaudata,   top,    Diachasmimorpha_longicaudata.bed
        0.44,   0.2,    0.95,   0,      ,       Apis_mellifera, top,    Apis_mellifera.bed
        0.33,   0.2,    0.95,   0,      ,       Vespa_crabro,   top,    Vespa_crabro.bed
        0.22,   0.2,    0.95,   0,      ,       Solenopsis_invicta,     top,    Solenopsis_invicta.bed
        0.11,   0.2,    0.95,   0,      ,       Neodiprion_pinetum,     bottom, Neodiprion_pinetum.bed
#       edges
e,      0,      1,      Quadrastichus_erythrinae.Nasonia_vitripennis.anchors.simple
e,      1,      2,      Nasonia_vitripennis.Belonocnema_kinseyi.anchors.simple
e,      2,      3,      Belonocnema_kinseyi.Diachasmimorpha_longicaudata.anchors.simple
e,      3,      4,      Diachasmimorpha_longicaudata.Apis_mellifera.anchors.simple
e,      4,      5,      Apis_mellifera.Vespa_crabro.anchors.simple
e,      5,      6,      Vespa_crabro.Solenopsis_invicta.anchors.simple
e,      6,      7,      Solenopsis_invicta.Neodiprion_pinetum.anchors.simple
```

And the chrs.txt file. I took this ordering from the chromsyn ordering. 

```
chr1,chr2,chr3,chr4,chr5
chr4,chr2,chr1,chr3,chr5
chr1,chr10,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9
chr10,chr13,chr18,chr19,chr5,chr6,chr1,chr15,chr2,chr17,chr14,chr4,chr11,chr12,chr16,chr20,chr7,chr8,chr3,chr9
LG14,LG8,LG10,LG5,LG9,LG12,LG15,LG1,LG16,LG2,LG13,LG4,LG7,LG11,LG6,LG3
chr2,chr19,chr10,chr21,chr25,chr18,chr24,chr1,chr11,chr12,chr13,chr14,chr16,chr17,chr20,chr22,chr23,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr15
chr1,chr12,chr13,chr14,chr2,chr6,chr7,chr9,chr3,chr10,chr11,chr5,chr4,chr15,chr16,chr8
chr6,chr4,chr5,chr1,chr3,chr7,chr2
```

And create karyoplot:

```bash
python -m jcvi.graphics.karyotype \
	--format png --font Arial --seed 1 \
	-o ../output_jcvi/20260102_gene_alignments_min10.pdf \
	chrs.txt chr_layout.txt --basepair
```

For re-filtering:

```bash
cat Pairs.list  | xargs -L1 ./02_Alignments.sh; cd input_jcvi; ./03_Plot.sh
```

Dotplots:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

reference=$1  #reference
target=$2 #target
echo "Working on $target aligned to $reference"

WD=/90daydata/coffea_pangenome/fly_alignments/gffs/input_jcvi/
cd ${WD}

python -m jcvi.graphics.dotplot --figsize 4x4 --minfont 8 --font Arial --style white ${reference}.${target}.anchors.minize10
```
