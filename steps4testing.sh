1. Create virtual conda env.
# Export the isoquant conda env
# conda env export -n pacbio > environment-pacbio.yaml

2. Install:
# can use the yaml file
    Python=3.8.20
    samtools=1.10 or 1. 9
    openssl=1.1.1w     
    minimap2=2.26 (ref genome) or 2.22 (conda env)
    hisat2=2.2.1
    pandas-1.3.5
    numpy-1.21.6
    
3. Download and check the bam file.
wget https://downloads.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/FullLengthReads/flnc.bam

4. Split the entire bam file into multiple bam files. Each one small bam file contains reads from one sample.
    Split into 6 bam files representing 6 samples.
subreads/
|__ COLO8299T
|   |__ m54026_190120_000756.subreads.bam 
|   |__ m54026_190120_000756.subreads.bam.pbi 
|   |__ m54119_190202_095143.subreads.bam 
|   |__ m54119_190202_095143.subreads.bam.pbi 
|   |__ m54119_190203_061153.subreads.bam 
|   |__ m54119_190203_061153.subreads.bam.pbi 
|__ COLO829BL
|   |__ m54019_190120_021709.subreads.bam 
|   |__ m54019_190120_021709.subreads.bam.pbi 
|   |__ m54119_190131_171128.subreads.bam 
|   |__ m54119_190131_171128.subreads.bam.pbi
|   |__ m54119_190201_133141.subreads.bam 
|   |__ m54119_190201_133141.subreads.bam.pbi 
|__ md5sum.txt

5. Convert each bam to fasta individually, and extract the transcript ids and corresponding raw long read ids, output into a txt file.

# Search and extract the reads identified in a certain sample, and generate the subset sam file.
export sampleID=m54119_190202_095143; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
# Convert the sam file to fasta file.
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

export sampleID=m54026_190120_000756; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

export sampleID=m54119_190203_061153; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

export sampleID=m54019_190120_021709; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

export sampleID=m54119_190131_171128; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

export sampleID=m54119_190201_133141; samtools view transcripts.bam -h |perl -anE 'if
 (/^\@/){print}else{my @aa=split(",",$F[12]);@aa=grep{/$ENV{sampleID}/}@aa;
$F[12]=join(",",@aa);
$F[12]="im:Z:".$F[12] unless ($F[12]=~/^im:Z/ and @aa>0);
print join("\t",@F),"\n" if @aa>0;}
' >$sampleID.sam 
bash ./PreparePacbioClusteredBAM2TAGET.sh $sampleID.sam $sampleID

6. Prepare the reference genome
# Downloaded hg38.fa and gtf files from GENECODE or Ensembl

7. Generate CPM files from fasta. 
#These codes prepare CPM files for TAGET; please update the config file for each sample.
ls *fasta | perl -ane 'chomp;$b=$_;$b=~s/\.fasta/_for_TAGET.tpm.txt/;print $b,"\n";`python Py4LongReadRNASeq/DetermineExp4ClusteredReads.py -f $_ -o $b`'

8. Prepare the configure file and add the CPM path into the configure file for each sample.

9. Run TAGET to refine mapping and quantification using the configure file prepared in the last step.
ls *TransAnnot.Config-* | perl -ane '
chomp;
my $t=$_;
`python ./RefinedLongReadMappingAndQuantification.py -c $t`;
'

10. Merge outputs, i.e. identified isoform for each sample.
# Script to prepare configure file
ls PATH_TO_FOLDER_OF_RESULTS | perl -F'\/' -MFile::Basename -ane 'chomp;my $f=basename($_);$d=dirname($_);$f=~s/\.stat//;my $s=$f;$s=~s/\.fa\.annot//;print $s,"\t",$_,"\t","$d/$f".".bed","\t","$d/$f".".db.pickle","\n";' > merge.config

# Merge
python ./TransAnnotMerge.py -c merge.config -o outputdir -m FLC
# The output includes the gene.exp, transcript.exp and others.

# Note: these results were generated with the original gtf file downloaded from GENECODE or Ensembl.
# We can rerun the Isoquant to generate a new gft file containing the known and novel transctipts, as a new gtf input to TAGET to get the gene.exp and transcript.exp again. This time, we don't have time to run this analysis and display the new results.

9. Optional: Isoquant
# Create the Isoquant conda env
# https://github.com/ablab/IsoQuant
conda create -c conda-forge -c bioconda -n isoquant python=3.8 isoquant -y
# IsoQuant version 3.7.1
# test
conda activate isoquant 
isoquant.py --test
# TEST PASSED CORRECTLY 
# Export the isoquant conda env
conda env export -n isoquant > environment-isoquant.yaml

# Download the fastq file
wget https://downloads.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/FullLengthReads/flnc.fastq
# Run Isoquant
# Don't need to create the output directory in advance, the isoquant function will create it and output the results.
conda activate isoquant
isoquant.py --reference /home/lead/notebooks/database4TAGET/hg38.fa \
--genedb /home/lead/notebooks/database4TAGET/hg38.ensembl.gtf \
--fastq /home/lead/notebooks/flnc.fastq \
--data_type pacbio_ccs -o output_isoquant
# Results are in: /home/lead/notebooks/output_isoquant/ and /home/lead/notebooks/output_isoquant/OUT
# The new gft file is /home/lead/notebooks/output_isoquant/OUT/OUT.transcript_models.gtf



10. Compare with the results, including gene.exp and transcript.exp, produced by the PacBio pipeline provided by the PacBio company.
# To save time, we performed the analysis in the St. Jude HPC this time.
#####################Pacbio updated pipeline, including IsoQuant#############
module load parallel/20240222
export USER=zcheng
export HG38=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq/RefGenomes/Human_hg38_Gencode_v39
export REFS=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq

mkdir -p /scratch_space/$USER/KinnexBulkIsoSeqAnalysis
cd /scratch_space/$USER/KinnexBulkIsoSeqAnalysis

singularity pull docker://quay.io/biocontainers/pbskera:1.4.0--hdfd78af_0
singularity pull docker://quay.io/biocontainers/lima:2.13.0--h9ee0642_0
singularity pull docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0
singularity pull docker://quay.io/pacbio/pbmm2:1.17.0_build1
singularity pull docker://quay.io/biocontainers/pbfusion:0.5.1--hdfd78af_0

# Run ISOSEQ REFINE
# Don't need to run this step, because we have downloaded the flnc.bam in the beginning.
# wget https://downloads.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/FullLengthReads/flnc.bam
# If we don't have the refined file, flnc.bam, the user needs to run this code.
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq refine --require-polya -j 16 \
	02_COLO829Pacbio.{1}.fl.{2}.bam \
	$REFS/bulkRNA/02_lima-primers/IsoSeq_v2_primers_12.fasta \
	03_COLO829Pacbio.{1}.flnc.bam

### x.sh
# Use previously generated or downloaded flnc BAM;
https://downloads.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/
https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/Kinnex_IsoSeq_Pipelines/KinnexBulkIsoSeq.lsf
# extract reads for each sample;
# Search and extract the reads identified in a certain sample, and generate the subset sam file.
export sampleID=$1; samtools view flnc.bam -h |perl -anE 'if
 (/^\@/){print}else{print if /$ENV{sampleID}/}
' |samtools view -Sb - -o $sampleID.bam;
samtools index $sampleID.bam
# put the above into x.sh to run it for each sample in the cluster;

bsub_Grace_Next -n 1 -m 20 "bash x.sh m54019_190120_021709"
bsub_Grace_Next -n 1 -m 20 "bash x.sh m54026_190120_000756"
bsub_Grace_Next -n 1 -m 20 "sleep 40s;bash x.sh m54119_190202_095143"
bsub_Grace_Next -n 1 -m 20 "sleep 100s;bash x.sh m54119_190203_061153"
bsub_Grace_Next -n 1 -m 20 "sleep 200s;bash x.sh m54119_190131_171128"
bsub_Grace_Next -n 1 -m 20 "sleep 300s;bash x.sh m54119_190201_133141"
 
### y.sh
export USER=zcheng
export HG38=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq/RefGenomes/Human_hg38_Gencode_v39
export REFS=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq
realpath $1.bam > 03_COLO829Pacbio.$1.flnc.fofn
# Run PBMM2
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pbmm2:1.17.0_build1 \
	pbmm2 align -j 16 --preset ISOSEQ --sort \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	03_COLO829Pacbio.$1.flnc.fofn \
	04_COLO829Pacbio.$1.align.bam
#put the above codes into y.sh to run it for each sample;

#Run all samples using the y.sh
ls m*bam | perl -pe 's/\.bam//' | perl -ane 'chomp;`bsub_Grace_Next -n 1 -m 100 "bash y.sh $_"`'
	
realpath 04_COLO829Pacbio.*.align.bam > 04_COLO829Pacbio.align.fofn
perl -pe 's/.*kinnex.//g' 04_COLO829Pacbio.align.fofn | perl -pe 's/.align.bam//g' > 04_COLO829Pacbio.align.labels

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant_generateYAML.py -b 04_COLO829Pacbio.align.fofn -l 04_COLO829Pacbio.align.labels -e 05_COLO829Pacbio.isoquant -o 05_isoquant4COLO829.yaml
rm 05_COLO829Pacbio.isoquant -rf 
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0 \
	isoquant.py -t 64 -d pacbio --yaml 05_isoquant4COLO829.yaml \
	-r $HG38/human_GRCh38_no_alt_analysis_set.fasta \
	-g $HG38/gencode.v39.annotation.sorted.gtf.db --complete_genedb \
	-o 05_COLO829Pacbio.isoquant \
	--sqanti_output

ln -s 05_COLO829Pacbio.isoquant/05_COLO829Pacbio.isoquant/05_COLO829Pacbio.isoquant.* .

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant2pigeon.py \
	--gtf 05_COLO829Pacbio.isoquant.transcript_models.gtf \
	--tsv 05_COLO829Pacbio.isoquant.transcript_model_grouped_counts.tsv \
	--output 06_COLO829Pacbio.pigeon.transcript_model_grouped_counts.csv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon classify -j 64 -o 06_COLO829Pacbio.pigeon \
	05_COLO829Pacbio.isoquant.transcript_models.gtf \
	$HG38/gencode.v39.annotation.sorted.gtf \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	--flnc 06_COLO829Pacbio.pigeon.transcript_model_grouped_counts.csv \
	--cage-peak $HG38/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
	--poly-a $HG38/polyA.list.txt \
	--coverage $HG38/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon report -j 64 06_COLO829Pacbio.pigeon_classification.txt 06_COLO829Pacbio.pigeon_classification.report.txt

cat 05_COLO829Pacbio.isoquant.extended_annotation.gtf | grep transcript | perl -pe 's/.*gene_id.{2}(ENSG\d+\.\d+).; transcript_id .(ENST\d+\.\d+).;.*; .*gene_name "(\S+)".*transcript_name "(\S+)".*/$1\t$2\t$3\t$4/' |sort -u |grep gene_id -v >ensembl_gene2tx.txt

#Clean the file headers;
perl -i.bak -pe 's/\/[^\t]+\/04_COLO829Pacbio.//g' 05_COLO829Pacbio.isoquant.transcript_model_grouped_counts.tsv
perl -i.bak -pe 's/\/[^\t]+\/04_COLO829Pacbio.//g' 05_COLO829Pacbio.isoquant.gene_grouped_counts.tsv
ln -s 05_COLO829Pacbio.isoquant.transcript_model_grouped_counts.tsv 07_isoquant.isoforms4COLO829.matrix
ln -s 05_COLO829Pacbio.isoquant.gene_grouped_counts.tsv 07_isoquant.genes4COLO829.matrix
#Sample matrix, including sample grp and sample ids;
#COLO8299T	m54026_190120_000756 
#COLO8299T	m54119_190202_095143  
#COLO8299T	m54119_190203_061153
#COLO829BL	m54019_190120_021709
#COLO829BL	m54119_190131_171128
#COLO829BL	m54119_190201_133141

11. DEG analysis
11.1 For output from PacBio pipeline
11.1.1 Analyse the isoform count matrix.
# Modify the isoform matrix (transcript.exp) from Pacbio pipeline to map the transcript id to the gene name
vlookup 07_isoquant.isoforms4COLO829.matrix 1 ensembl_gene2tx.txt 2 3,4 y |perl -ane '$F[0]=$F[-1] unless /NaN/;print join("\t",@F[0..($#F-2)]),"\n";' >07_isoquant.isoforms4COLO829.matrix.new

# Run DESeq2
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.isoforms4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.isoforms

# Select the top 20 isoforms.
head -n21 07_deseq2.isoforms/07_isoquant.isoforms4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.isoforms.new.txt

11.1.2 The same strategy applied to the gene.exp matrix
# Modify the gene matrix (gene.exp) from Pacbio pipeline to map the gene id to the gene name
vlookup 07_isoquant.genes4COLO829.matrix 1 ensembl_gene2tx.txt 1 3 y |perl -ane '$F[0]=$F[-1] unless /NaN/;print join("\t",@F[0..($#F-1)]),"\n";'|SortFileByCols.sh - '-k1,1 -u ' 1  >07_isoquant.genes4COLO829.matrix.new 

# Run DESeq2
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.genes4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.genes

# Select the top 20 genes.
head -n21 07_deseq2.genes/07_isoquant.genes4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.genes.new.txt


11.2 Analyze output from our modified TAGET pipeline using the same strategy used for the output from the standard PacBio pipeline.
11.2.1 Analyse the isoform count matrix.
cat ~/working_scripts/Py4LongReadRNASeq/outputdir/transcript.exp |delete_column 1 |perl -pe 's/Transcript/#feature_id/'|SortFileByCols.sh - '-k1,1 -u ' 1 >07_TAGET.isoforms4COLO829.matrix.new

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_TAGET.isoforms4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.isoforms_TAGET

head -n21 07_deseq2.isoforms_TAGET/07_TAGET.isoforms4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.isoforms.TAGET.new.txt

11.2.2 Analyse the gene count matrix.
cat ~/working_scripts/Py4LongReadRNASeq/outputdir/gene.exp |perl -pe 's/Gene/#feature_id/'|SortFileByCols.sh - '-k1,1 -u ' 1 >07_TAGET.genes4COLO829.matrix.new
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_TAGET.genes4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.genes_TAGET

head -n21 07_deseq2.genes_TAGET/07_TAGET.genes4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.genes.TAGET.new.txt

# Now we obtained the DE results of isoforms and genes, from our modified TAGET and the standard PacBio pipelines, respectively.

