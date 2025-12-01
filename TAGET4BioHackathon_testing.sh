cd /research/groups/cab/projects/automapper/common/zhongshan/test/NGS_Pipelines/TAGET/BioHackathon

cd /home/zcheng/working_scripts/Py4LongReadRNASeq
#python ./DetermineExp4ClusteredReads.py
ls *fasta | perl -ane 'chomp;$b=$_;$b=~s/\.fasta/_for_TAGET.tpm.txt/;print $b,"\n";`python ./DetermineExp4ClusteredReads.py -f $_ -o $b`'
#ls ../*fasta | perl -pe 's/.fasta//;s/\.\.\///'|grep test -v |grep m54119_190202_095143 -v | perl -ane 'chomp;my $t=$_;open F,"TransAnnot.Config-m54119_190202_095143";open O,">TransAnnot.Config-$t" or die;while (<F>){my $l=$_;$l=~s/m54119_190202_095143/$t/g;print O $l;};close F;close O;' 

ls *fasta | perl -pe 's/.fasta//;'| perl -ane 'chomp;my $t=$_;open F," TransAnnot.Config";open O,">TransAnnot.Config-$t" or die;while (<F>){my $l=$_;$l=~s/SAMPLE/$t/g;print O $l;};close F;close O;' 


ls *TransAnnot.Config-* | perl -ane '
chomp;
my $t=$_;
`
bsub_Grace_Next -n 1 -m 300 "
source activate /research_jude/rgs01_jude/groups/kennegrp/projects/kennegrp_cab/common/Pacbio_Isoseq/Shaohua_HiFi_Isoseq_New/TAGET
module load minimap2/2.26 samtools/1.17 
module load hisat/2.2.1
python3 ./RefinedLongReadMappingAndQuantification.py -c $t
"
`;

'


#Merge output from each sample.
# Script to prepare configure file
ls m*/*stat | perl -F'\/' -MFile::Basename -ane 'chomp;my $f=basename($_);$d=dirname($_);$f=~s/\.stat//;my $s=$f;$s=~s/\.annot//;print $s,"\t",$_,"\t","$d/$f".".bed","\t","$d/$f".".db.pickle","\n";' >merge.config

# Merge
#it is necessary to create the output dir;
mkdir outputdir;
python ./TransAnnotMerge.py -c merge.config -o outputdir -m TPM
mkdir outputdir4flc
python ./TransAnnotMerge.py -c merge.config -o outputdir4flc -m FLC

cat outputdir/transcript.exp  |grep '(Gene|SBDS)' -P | Check_Header_DataRow

##########DEG2 analysis;
cd /scratch_space/zcheng/KinnexBulkIsoSeqAnalysis
export USER=zcheng
export HG38=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq/RefGenomes/Human_hg38_Gencode_v39
export REFS=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq

cat ~/working_scripts/Py4LongReadRNASeq/outputdir4flc/transcript.exp |delete_column 1 |perl -pe 's/Transcript/#feature_id/'|SortFileByCols.sh - '-k1,1 -u ' 1 >07_TAGET.isoforms4COLO829.matrix.new

cat ~/working_scripts/Py4LongReadRNASeq/outputdir4flc/gene.exp |perl -pe 's/Gene/#feature_id/'|SortFileByCols.sh - '-k1,1 -u ' 1 >07_TAGET.genes4COLO829.matrix.new
 
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_TAGET.isoforms4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.isoforms_TAGET

head -n21 07_deseq2.isoforms_TAGET/07_TAGET.isoforms4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.isoforms.TAGET.new.txt


singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_TAGET.genes4COLO829.matrix.new \
	--method DESeq2 \
	--samples_file 07_samples.matrix4COLO829.txt \
	--output 07_deseq2.genes_TAGET

head -n21 07_deseq2.genes_TAGET/07_TAGET.genes4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.genes.TAGET.new.txt


tar -zcvf DEG_Rests.tgz 07_deseq2.*


##compare TAGET and Pacbio results;

07_deseq2.genes_TAGET/07_TAGET.genes4COLO829.matrix.new.COLO8299T_vs_COLO829BL.DESeq2.DE_results



#So it is feasible to add transcript names and genesymbols for each long-read;
#and then merge these regions by transcript names;
cat m54119_190131_171128/*.annot.cluster.transcript| perl -ane 'if ($.==1){print;}else{my @aa=split(",",$F[-1]);foreach my $a(@aa){print join("\t",@F[0..$#F-1],$a),"\n";}}' | vlookupMany2Many - -1 m54119_190131_171128/*.annot.bed 4 1,2,3 y | perl -ane 'print join("\t",@F[5..7],"exon",$F[1],$F[2].".".$F[4]),"\n";' >m54119_190131_171128for_matlab.bed
#cat SBDS_tgt.bed | bed2gff_used_by_matlab.pl -
cat m54119_190131_171128for_matlab.bed| bed2gff_used_by_matlab.pl - >m54119_190131_171128.gff
cat /research/groups/cab/projects/automapper/common/zhongshan/test/NGS_Pipelines/TAGET/hg38.ensembl.gtf >>m54119_190131_171128.gff

#Note: CLC genomic workbench can load GTF file with genomic reference;
#which would be even better in terms of its pretty gene tracks;

perl -i -pe 's/^/chr/' m54119_190131_171128.gff

#In matlab;
%!perl -i -pe 's/^/chr/' SBDS_tgt.gff
x=GTFAnnotationMy('m54119_190131_171128.gff');

%For SBDS;
chr='chr7';
Start=66977703;
Stop=67003493;

NumOfTgtPlots=0;
add_path2MATLAB_PATH
[ha,gene_exon_info]=GTF_gene_track4OtherPlots(x,chr,Start,Stop,NumOfTgtPlots,[1],0.01,1,'SBDS',1,500,0.1);
trsnames=flipud(unique(gene_exon_info.GeneName));
trsnames=flipud(trsnames);
%do not need to remove the extra prefix in the trsnames;
%trsnames=regexprep(trsnames,'^SBDS.','');
numgrps=grp2idx(contains(trsnames,'ENST'));
[ha,gene_exon_info]=GTF_gene_track4OtherPlots(x,chr,Start,Stop,NumOfTgtPlots,[1],0.01,1,'SBDS',1,500,0.1,[],trsnames,1,numgrps,[0.9 0.9 0.9;0.8 1 1]);


%For CALU
%Chromosome 7: 128,739,292-128,773,400
chr='chr7';
Start=128739292;
Stop=128773400;
genesymbol='CALU';
NumOfTgtPlots=0;
add_path2MATLAB_PATH
[ha,gene_exon_info]=GTF_gene_track4OtherPlots(x,chr,Start,Stop,NumOfTgtPlots,[1],0.01,1,genesymbol,1,500,0.1);
trsnames=flipud(unique(gene_exon_info.GeneName));
%trsnames=flipud(trsnames);
numgrps=grp2idx(contains(trsnames,'ENST'));
[ha,gene_exon_info]=GTF_gene_track4OtherPlots(x,chr,Start,Stop,NumOfTgtPlots,[1],0.01,1,genesymbol,1,500,0.1,[],trsnames,1,numgrps,[0.9 0.9 0.9;0.8 1 1]);
