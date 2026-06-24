# Target: To compare genetic variations between the wild-type strain (SRR21747981) and the mutant-type strain (SRR21747980) of Saccharomyces cerevisiae using Whole Genome Sequencing (WGS)
# Data source article: Lu et al., 2023, Omics Sequencing of Saccharomyces cerevisiae Strain with Improved Capacity for Ethanol Production
# Phenotype: The wild-type strain produced 13.72% ethanol, whereas the mutant-type strain showed a higher ethanol production of 16.13%
# Note: The original article used BWA, GATK, and ANNOVAR, while this re-analysis used BWA, FreeBayes, and SnpEff. This analysis was conducted for learning purposes only and is not intended for publication.

# Setup

conda create -n wgs python=3.9 -y
conda activate wgs

conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install -y sra-tools fastqc multiqc fastp bwa samtools freebayes

conda install -c conda-forge openjdk=21
conda install -c bioconda snpeff

prefetch --version
fastq-dump --version | head -1
fastqc --version
fastp --version
bwa 2>&1 | head -3
samtools --version | head -1
freebayes --version 2>&1 | head -1
snpEff -version

# Data Preparation

## Data Preparation for wild-type
mkdir -p ~/wgs_yeast/01_wild_type/01_DataPreparation
cd wgs_yeast/01_wild_type/01_DataPreparation

### Download Data
prefetch SRR21747981
fastq-dump --split-files SRR21747981/SRR21747981.sra

### Quality Control (fastqc) and Trimming (fastp)
fastqc SRR21747981_1.fastq -o .
fastqc SRR21747981_2.fastq -o .

explorer.exe SRR21747981_1_fastqc.html
explorer.exe SRR21747981_2_fastqc.html

fastp -i SRR21747981_1.fastq -I SRR21747981_2.fastq -o wild_clean_1.fastq -O wild_clean_2.fastq\
 --detect_adapter_for_pe --html fastp_wild_report.html --json fastp_wild_report.json
explorer.exe fastp_wild_report.html

fastqc wild_clean_1.fastq -o .
fastqc wild_clean_2.fastq -o .

explorer.exe wild_clean_1_fastqc.html
explorer.exe wild_clean_2_fastqc.html

## Data Preparation for mutant-type
mkdir -p ~/wgs_yeast/02_mutant_type/01_DataPreparation
cd wgs_yeast/02_mutant_type/01_DataPreparation

### Download Data
prefetch SRR21747980
fastq-dump --split-files SRR21747980/SRR21747980.sra

### Quality Control (fastqc) and Trimming (fastp)
fastqc SRR21747980_1.fastq -o .
fastqc SRR21747980_2.fastq -o .

explorer.exe SRR21747980_1_fastqc.html
explorer.exe SRR21747980_2_fastqc.html

fastp -i SRR21747980_1.fastq -I SRR21747980_2.fastq -o mutant_clean_1.fastq -O mutant_clean_2.fastq\
 --detect_adapter_for_pe --html fastp_mutant_report.html --json fastp_mutant_report.json
explorer.exe fastp_mutant_report.html

fastqc mutant_clean_1.fastq -o .
fastqc mutant_clean_2.fastq -o .

explorer.exe mutant_clean_1_fastqc.html
explorer.exe mutant_clean_2_fastqc.html

# Mapping

## Mapping for wild-type
cd ~/wgs_yeast/01_wild_type/
mkdir -p 02_Mapping/
cd 02_Mapping/
cp ../01_DataPreparation/wild_clean_1.fastq .
cp ../01_DataPreparation/wild_clean_2.fastq .

### Download Reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip GCF_000146045.2_R64_genomic.fna.gz
mv GCF_000146045.2_R64_genomic.fna reference.fna
bwa index reference.fna

### Mapping with BWA-MEM
bwa mem reference.fna wild_clean_1.fastq wild_clean_2.fastq > wild_aligment.sam

## Mapping for mutant-type
cd ~/wgs_yeast/02_mutant_type/
mkdir -p 02_Mapping/
cd 02_Mapping/
cp ../01_DataPreparation/mutant_clean_1.fastq .
cp ../01_DataPreparation/mutant_clean_2.fastq .

### Get Reference from wild-type mapping directory
cp ../../01_wild_type/02_Mapping/reference* .

### Mapping with BWA-MEM
bwa mem reference.fna mutant_clean_1.fastq mutant_clean_2.fastq > mutant_alignment.sam

# Convert, Sort, and Index

## Convert, Sort, and Index for wild-type
cd ~/wgs_yeast/01_wild_type/
mkdir -p 03_Convert-Sort/
cd 03_Convert-Sort/
cp ../02_Mapping/wild_aligment.sam .

samtools view -Sb wild_aligment.sam > wild_alignment.bam
samtools sort wild_aligment.bam -o wild_aligment.sorted.bam
samtools index wild_aligment.sorted.bam

## Convert, Sort, and Index for mutant-type
cd ~/wgs_yeast/02_mutant_type/
mkdir -p 03_Convert-Sort/
cd 03_Convert-Sort/
cp ../02_Mapping/mutant_alignment.sam .

samtools view -Sb mutant_alignment.sam > mutant_alignment.bam
samtools sort mutant_alignment.bam -o mutant_alignment.sorted.bam
samtools index mutant_alignment.sorted.bam

# Variant Calling

## Variant Calling for wild-type
cd ~/wgs_yeast/01_wild_type/
mkdir -p 04_VariantCalling
cd 04_VariantCalling
cp ../02_Mapping/reference.fna .
cp ../03_Convert-Sort/wild_aligment.sorted.bam* .

### Mapping with freebayes
freebayes -f reference.fna --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10\
 wild_aligment.sorted.bam > wild_variant_raw.vcf
### Filtering and normalization with bcftools
bcftools filter -e 'QUAL < 30 || INFO/DP < 10 || AF < 0.1' wild_variant_raw.vcf -o wild_variant_filtered.vcf
bcftools norm -f reference.fna -m- wild_variant_filtered.vcf -o wild_variant_norm.vcf

grep "#CHROM" wild_variant_norm.vcf && grep -v "^#" wild_variant_norm.vcf | head -5

## Variant Calling for mutant-type
cd ~/wgs_yeast/02_mutant_type/
mkdir -p 04_VariantCalling
cd 04_VariantCalling
cp ../02_Mapping/reference.fna .
cp ../03_Convert-Sort/mutant_alignment.sorted.bam* .

### Mapping with freebayes
freebayes -f reference.fna --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10\
 mutant_alignment.sorted.bam > mutant_variant_raw.vcf
### Filtering and normalization with bcftools
bcftools filter -e 'QUAL < 30 || INFO/DP < 10 || AF < 0.1' mutant_variant_raw.vcf -o mutant_variant_filtered.vcf
bcftools norm -f reference.fna -m- mutant_variant_filtered.vcf -o mutant_variant_norm.vcf

grep "#CHROM" mutant_variant_norm.vcf && grep -v "^#" mutant_variant_norm.vcf | head -5

# Create intersection of variants between wild-type and mutant-type
cd ~/wgs_yeast/
mkdir -p 03_var_compare_wild_mutant
cd 03_Comparison
cp ../01_wild_type/04_VariantCalling/wild_variant_norm.vcf .
cp ../02_mutant_type/04_VariantCalling/mutant_variant_norm.vcf .
mv wild_variant_norm.vcf wild.vcf
mv mutant_variant_norm.vcf mutant.vcf

## Compress and index VCF files
bgzip wild.vcf
bgzip mutant.vcf
tabix -p vcf wild.vcf.gz
tabix -p vcf mutant.vcf.gz

## Create intersection using bcftools isec
bcftools isec -p wild_mutant_isec wild.vcf.gz mutant.vcf.gz

# Annotation for mutant-specific variants
cd ~/wgs_yeast/
mkdir -p 04_annot_mutant_specific_var
cd 04_annot_mutant_specific_var
cp ../03_var_compare_wild_mutant/wild_mutant_isec/0001.vcf .
mv 0001.vcf mutant_specific_var.vcf

## Check the chromosome names in the VCF file and the SnpEff database
bcftools query -f '%CHROM\n' mutant_specific_var.vcf | sort -u
snpEff databases | grep -i "Saccharomyces_cerevisiae" | grep -i "R64"
snpEff dump -bed R64-1-1.99 | cut -f1 | sort -u

## Because the chromosome names in the VCF file and SnpEff database are different, we need to rename VCF file before annotation
### Check the description of the chromosome name in the reference file that is used by the VCF file
cp ../02_mutant_type_var/04_VariantCalling/reference.fna .
grep "^>" reference.fna | head -n 17

### Create a chromosome name mapping file (chr_map.txt) and rename the chromosome names in the VCF file
nano chr_map.txt
NC_001133.9 I
NC_001134.8 II
NC_001135.5 III
NC_001136.10 IV
NC_001137.3 V
NC_001138.5 VI
NC_001139.9 VII
NC_001140.6 VIII
NC_001141.2 IX
NC_001142.9 X
NC_001143.9 XI
NC_001144.5 XII
NC_001145.3 XIII
NC_001146.8 XIV
NC_001147.6 XV
NC_001148.4 XVI
NC_001224.1 Mito
bcftools annotate --rename-chrs chr_map.txt -O v -o mutant_specific_var_renamed.vcf mutant_specific_var.vcf

## Annotate the mutant-specific variants using SnpEff
snpEff ann -v -stats snpEff_summary.html -csvStats snpEff_summary.csv R64-1-1.99 mutant_specific_var_renamed.vcf > mutant_specific_var_renamed_annotated.vcf
explorer.exe snpEff_summary.html

## Extract annotation summary for all variants using SnpSift
snpSift extractFields -s "," mutant_specific_var_renamed_annotated.vcf CHROM POS REF ALT\
 ANN[*].GENEID "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > annotation_summary_all.tsv

## Next: Clean the annotation summary using python script 'Script_annotation_summary_high_moderate.py' to show only 'High' and 'Moderate' Impact
## Next: Compare the cleaned annotation summary with the list of genes that are known to be involved in ethanol production in S. cerevisiae from the original article

## snpEff summary visualizations
cd ~/wgs_yeast/

### Download and install ngi_visualizations
git clone https://github.com/ewels/ngi_visualizations.git
conda install -y setuptools matplotlib numpy
cd ngi_visualizations
python setup.py install

### Run the ngi_visualizations script
sed -i 's/axes.spines.itervalues()/axes.spines.values()/g' ngi_visualizations/snpEff/snpEff_plots.py
python ngi_visualizations/snpEff/snpEff_plots.py ~/wgs/mission04/snpEff_summary.csv

explorer.exe .