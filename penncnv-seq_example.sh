# penncnv-seq_example.sh
# C: Aug 14, 2015
# M: Dec 23, 2016
# A: Leandro Lima <lelimaufc@gmail.com>

#penncnv_dir=~/Dropbox/programming/penncnv-seq (@ data cruncher)


if [ "$#" -ne 5 ]; then
    echo -e "** Illegal number of parameters.**\n"
    echo -e "Usage:"
    echo -e "\t./penncnv-seq_example.sh [penncnv_dir] [penncnv_ref_dir] [genome_version] [population] [reference.fasta] [bam_file]\n"
    echo -e "Accepted populations: ALL, AFR, AMR, EAS, EUR and SAS"
    echo -e "Accepted genome versions: hg19 and hg38"
    echo
    echo -e "Example:"
    echo -e "./penncnv-seq_example.sh ~/PennCNV ~/reference hg19 AFR chr21.fasta chr21.bam"
    echo 
    exit 1
fi

# PennCNV directory
PennCNV_dir=$1

# PennCNV reference directory
ref_dir=$2

# Genome version
gv=$3

# Population
pop=$4

# reference path here (e.g. references/hg19/chr1.fa)
ref_fasta=$5

# Define your BAM file here
bam=$6

# Define your sample name
sample=`echo $bam | awk -F'/' '{print $NF}' | perl -pe 's/.bam//g'`

# Prepare files
echo " \
    perl $PennCNV_dir/convert_map2signal.pl \
        $bam \
        $ref_fasta \
    --outfile $PWD/$sample;
    rm $PWD/$sample.read1 $PWD/$sample.read2;" > map2signal.$sample.sh
# qsub -l h_vmem=1G -cwd -V map2signal.$sample.sh
bash map2signal.$sample.sh


# Define 
baf_tmp=$sample.baflrr.tmp.txt
pfb_bed=$ref_dir/$gv""_$pop.sites.2015_08.chrom$chrom.bed
baf_minus_pfb=$sample.chr$chrom.tmp_baf_minus_pfb.txt
read3=$sample.read3
pfb_markers_with_baflrr=$sample.chr$chrom.pfb_markers_with_baflrr.txt
echo "cat $read3 | tail -n +2 | awk '{region=split(\$1, array, \"-\"); print \"chr\"array[1]\"\\t\"array[2]\"\\t\"\$4\"\\t\"\$7}' | perl -pe 's/:/\t/g; s/chrchr/chr/g' > $baf_tmp" > penncnv1.$sample.chr$chrom.sh
echo "bedtools subtract  -a $baf_tmp -b $pfb_bed > $baf_minus_pfb" >> penncnv1.$sample.chr$chrom.sh
echo "bedtools intersect -a $pfb_bed -b $baf_tmp -wb | cut -f1-3,7,8 > $pfb_markers_with_baflrr" >> penncnv1.$sample.chr$chrom.sh
echo -e "Name\tB Allele Freq\tLRR\tChr\tPosition" > $sample.chr$chrom.baflrr
echo "cat $baf_minus_pfb $pfb_markers_with_baflrr | sort -k2n | awk -v chrom=$chrom '{print \$1\":\"\$2\"-\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"chrom\"\\t\"\$2}' >> $sample.chr$chrom.baflrr" >> penncnv1.$sample.chr$chrom.sh
echo "rm  $baf_minus_pfb $pfb_markers_with_baflrr $baf_tmp" >> penncnv1.$sample.chr$chrom.sh
#qsub -V -cwd -hold_jid map2signal.$sample.chr$chrom.sh penncnv1.$sample.chr$chrom.sh 
bash penncnv1.$sample.chr$chrom.sh 

# Run PennCNV
echo " \
    perl $PennCNV_dir/detect_cnv.pl \
        $sample.chr$chrom.baflrr \
        -wgs \
        -hmmfile wgs.hmm \
        -pfb $ref_dir/$gv""_$pop.sites.2015_08.pfb \
        -output $PWD/$sample.rawcnv" > detect_cnv.$sample.sh
# qsub -l h_vmem=1G -cwd -V detect_cnv.$sample.sh
bash detect_cnv.$sample.sh

