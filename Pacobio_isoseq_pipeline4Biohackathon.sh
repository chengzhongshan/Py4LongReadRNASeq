
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
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq refine --require-polya -j 16 \
	02_COLO829Pacbio.{1}.fl.{2}.bam \
	$REFS/bulkRNA/02_lima-primers/IsoSeq_v2_primers_12.fasta \
	03_COLO829Pacbio.{1}.flnc.bam
#Use previously generated flnc BAM;
https://downloads.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/
https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/Kinnex_IsoSeq_Pipelines/KinnexBulkIsoSeq.lsf
#extract reads for each sample;
# Search and extract the reads identified in a certain sample, and generate the subset sam file.
export sampleID=$1; samtools view flnc.bam -h |perl -anE 'if
 (/^\@/){print}else{print if /$ENV{sampleID}/}
' |samtools view -Sb - -o $sampleID.bam;
samtools index $sampleID.bam
#put the above into x.sh to run it for each sample in the cluster;
bsub_Grace_Next -n 1 -m 20 "bash x.sh m54019_190120_021709"
bsub_Grace_Next -n 1 -m 20 "bash x.sh m54026_190120_000756"
bsub_Grace_Next -n 1 -m 20 "sleep 40s;bash x.sh m54119_190202_095143"
bsub_Grace_Next -n 1 -m 20 "sleep 100s;bash x.sh m54119_190203_061153"
bsub_Grace_Next -n 1 -m 20 "sleep 200s;bash x.sh m54119_190131_171128"
bsub_Grace_Next -n 1 -m 20 "sleep 300s;bash x.sh m54119_190201_133141"

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

#put these into y.sh to run it for each sample;
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
vlookup 07_isoquant.isoforms4COLO829.matrix 1 ensembl_gene2tx.txt 2 3,4 y |perl -ane '$F[0]=$F[-1] unless /NaN/;print join("\t",@F[0..($#F-2)]),"\n";' >07_isoquant.isoforms4COLO829.matrix.new

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.isoforms4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.isoforms

head -n21 07_deseq2.isoforms/07_isoquant.isoforms4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.isoforms.new.txt

for i in `cat 07_deseq2.top20.isoforms.txt`
	do 
	grep $i 06_COLO829Pacbio.pigeon_classification.txt >> 07_deseq2.top20.isoforms.classify.txt
	done
vlookup 07_isoquant.genes4COLO829.matrix 1 ensembl_gene2tx.txt 1 3 y |perl -ane '$F[0]=$F[-1] unless /NaN/;print join("\t",@F[0..($#F-1)]),"\n";'|SortFileByCols.sh - '-k1,1 -u ' 1  >07_isoquant.genes4COLO829.matrix.new 

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.genes4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.genes

head -n21 07_deseq2.genes/07_isoquant.genes4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.genes.new.txt

for i in `cat 07_deseq2.top20.genes.txt`
	do 
	grep $i 05_COLO829Pacbio.isoquant.transcript_models.gtf | grep -P "\tgene\t" >> 07_deseq2.top20.genes.classify.txt
	done

#remove previous results of heart tissue;
#rm -rf  *heart*
