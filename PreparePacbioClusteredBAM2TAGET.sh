#!/bin/bash
if [ "$#" -lt 2 ]; then
 echo ""
 echo "Usage: PreparePacbioClusteredBAM2TAGET.sh
(1) bam_fullpath (pacbio isoseq cluster2 generated 
    trancript sam/bam fullpath)
(2) output_prefix (output prefix for two files, 
    including one fasta file and the other one
    flnc_count.txt for the software TAGET!)
(3) isbam (default value is 1, otherwise provide value 0 for
    processing the input as sam!)

"
 echo "";
 exit 0;
fi
export bam=$1;
export output=$2;

if [ "$#" -lt 3 ]; then
    export SamPar='';
fi
if [ "$#" -eq 3 && "$3" -eq 0 ]; then
    export SamPar='-S';
fi

samtools view $SamPar $bam |perl -ane '
$F[12]=~s/im:Z://;
my @aa=split(",",$F[12]);
print "$F[0],",scalar @aa,",$F[12]\n"
' >$output.flnc_count.txt
 samtools view $SamPar $bam -h |samtools fasta -F 0x900 - >$output.tmp

cat $output.tmp | perl -anE '
chomp;
$l=$_;
my @hds=split("\\||:",$l) if $l=~/^>/;
my $id=$hds[0];
   $id=~s/^\>//;
#print join("\n",@hds),"\n" and exit;
if ($.==1){
our %FLC;
open F,"$ENV{output}.flnc_count.txt" or die;
while(my $k=<F>){
#next if $k=~/^id/;
chomp $k;
my ($ID,$cnt)=split(",",$k);
$FLC{$ID}=$cnt;
#print $ID,"\t",$FLC{$ID},"\n";
}
close F;
print join(" ",@hds," full_length_coverage=".$FLC{$id}),"\n";
}else{
if ($l=~/^>/){
print STDERR "the $id is not defiend in the flnc_count.txt file:\n",
"$ENV{sample_dir}/$ENV{sample_prefix}.flnc_count.txt\n" and exit unless defined $FLC{$id};
print join(" ",@hds," full_length_coverage=".$FLC{$id}),"\n";
}else{
print $l,"\n";
}
}
' >$output.fasta
rm $output.tmp

