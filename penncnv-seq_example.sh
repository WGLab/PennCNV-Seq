# penncnv-seq_example.sh
# C: Aug 14, 2015
# M: Dec  9, 2016
# A: Leandro Lima <lelimaufc@gmail.com>

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters."
    echo -e "You have to set it with the path to the reference (e.g. references/hg19.fa)"

    exit 1
fi



# place your reference path here (e.g. references/hg19.fa)
ref=$1

# Define your BAM file here
bam=$2

# Define your sample name
sample=`echo $bam | awk -F'/' '{print $NF}' | perl -pe 's/.bam//g'`


echo " \
    perl convert_map2signal.pl \
        $bam \
        $ref \
    --outfile $PWD/$sample;
    rm $PWD/$sample.read1 $PWD/$sample.read2;" > map2signal.$sample.sh
# qsub -l h_vmem=1G -cwd -V map2signal.$sample.sh
bash qsub -l h_vmem=1G -cwd -V map2signal.$sample.sh

# Define 
baf_tmp=$sample.chr$chrom.baflrr.tmp.txt
pfb_txt=$HOME/references/penncnv/$gv""_$pop.sites.2015_08.chrom$chrom.pfb.txt
baf_minus_pfb=$sample.chr$chrom.tmp_baf_minus_pfb.txt
read3=$sample.chr$chrom.exome.read3
pfb_markers_with_baflrr=$sample.chr$chrom.pfb_markers_with_baflrr.txt
echo "cat $read3 | tail -n +2 | awk '{region=split(\$1, array, \"-\"); print \"chr\"array[1]\"\\t\"array[2]\"\\t\"\$4\"\\t\"\$7}' | perl -pe 's/:/\t/g; s/chrchr/chr/g' > $baf_tmp" >> penncnv1.$sample.chr$chrom.sh
echo "bedtools subtract  -a $baf_tmp -b $pfb_txt > $baf_minus_pfb" >> penncnv1.$sample.chr$chrom.sh
echo "bedtools intersect -a $pfb_txt -b $baf_tmp -wb | cut -f1-3,8,9 > $pfb_markers_with_baflrr" >> penncnv1.$sample.chr$chrom.sh
echo -e "Name\tB Allele Freq\tLRR\tChr\tPosition" > $sample.chr$chrom.baflrr
echo "cat $baf_minus_pfb $pfb_markers_with_baflrr | sort -k2n | awk -v chrom=$chrom '{print \$1\":\"\$2\"-\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"chrom\"\\t\"\$2}' >> $sample.chr$chrom.baflrr" >> penncnv1.$sample.chr$chrom.sh
echo "rm  $baf_minus_pfb $pfb_markers_with_baflrr $baf_tmp" >> penncnv1.$sample.chr$chrom.sh
qsub -V -cwd -hold_jid map2signal.$sample.chr$chrom.sh penncnv1.$sample.chr$chrom.sh 



