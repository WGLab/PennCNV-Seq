# svgen_simulations_penncnvseq.sh
# C: June 10, 2016 
# M: June  2, 2017
# A: Leandro Lima <lelimaufc@gmail.com>

SVGen_dir=/home/llima/tools/SVGen/
ref_dir=/home/llima/references/
email=lelimaufc@gmail.com

gv=hg38 # genome version
# populations = AFR (African), AMR (Ad Mixed American), EAS (East Asian), EUR (European) and SAS (South Asian) 
pop_mother=EUR
pop_father=AFR

echo "pop_mother = $pop_mother" >  parents_population.txt
echo "pop_father = $pop_father" >> parents_population.txt

# working_dir=/home/llima # TO DO: change this in the code, both for reference and analysis
# cd /home/llima/projects/genome_simulations/SVGen_tests/
# cd $gv

# cd /home/llima/projects/genome_simulations/SVGen_tests/for_PennCNV_Seq/sample_1

ln -s ~/tools/SVGen/reference

# 1kb, 1.5kb, 2kb, 2.5kb, 3kb, 3.5kb, 4kb, 5kb, 6kb, 8kb, 10kb, 20kb, 30kb, 40kb, 50kb, 75kb, 100kb, 150kb, 200kb, 500kb, 1mb, 5mb
echo -e "1000\n1500\n2000\n2500\n3000\n3500\n4000\n5000\n6000\n8000\n10000\n20000
30000\n40000\n50000\n75000\n100000\n150000\n200000\n500000\n1000000\n5000000" > SV_lengths.txt

# exome
echo -e "10000\n20000\n30000\n40000\n50000\n75000\n100000\n150000\n200000\n500000\n1000000\n5000000" > SV_lengths.txt

for sample in {1..20}; do
    echo
    mkdir sample_$sample
    cd sample_$sample
    # Generate SVs' BED
    for chrom in {1..22} X; do
        echo chrom $chrom
        python -c "import sys, random; numbers = open(sys.argv[1]).read().split(); print '\n'.join(random.sample(numbers, 5))" ../SV_lengths.txt > temp_lens_del.txt
        python -c "import sys, random; numbers = open(sys.argv[1]).read().split(); print '\n'.join(random.sample(numbers, 5))" ../SV_lengths.txt > temp_lens_dup.txt
        python $SVGen_dir/simulate_SV_BED.py \
            --dup_lens temp_lens_dup.txt \
            --del_lens temp_lens_del.txt \
            --chroms $chrom \
            --distance 1000000 \
            --chrom_lens ../reference/chrom_lengths_$gv.txt \
            --gaps ../reference/gaps_$gv.txt \
            -o mother_SV_$gv""_chr$chrom.bed
    done
    echo chrom $chrom
    python -c "import sys, random; numbers = open(sys.argv[1]).read().split(); print '\n'.join(random.sample(numbers, 4))" ../SV_lengths.txt > temp_lens_del.txt
    python -c "import sys, random; numbers = open(sys.argv[1]).read().split(); print '\n'.join(random.sample(numbers, 4))" ../SV_lengths.txt > temp_lens_dup.txt
    python $SVGen_dir/simulate_SV_BED.py \
        --dup_lens temp_lens_dup.txt \
        --del_lens temp_lens_del.txt \
        --chroms $chrom \
        --distance 100000 \
        --chrom_lens ../reference/chrom_lengths_$gv.txt \
        --gaps ../reference/gaps_$gv.txt \
        -o mother_SV_$gv""_chr$chrom.bed
    rm temp_lens_del.txt temp_lens_dup.txt

    cat mother_SV_*bed | sort -k1V -k2n > SVs_mother_$gv.bed
    rm mother_SV_*bed

    # Simulating CNVs with 4 copies
    for chrom in {1..22}; do grep ^$chrom SVs_mother_$gv.bed | awk '$3-$2 < 90000 && $4 == "dup"' | head -n2; done > SVs_father_$gv.bed
    # Simulating CNVs with 0 copies
    for chrom in {1..22}; do grep ^$chrom SVs_mother_$gv.bed | awk '$3-$2 < 90000 && $4 == "del"' | head -n2; done >> SVs_father_$gv.bed
    # Simulating loss-of-heterozigosity (LOH)
    for chrom in {1..22}; do grep ^$chrom SVs_mother_$gv.bed | grep dup | tail -n2 | awk '$3-$2+1 >= 5000000 {print $1"\t"$2"\t"$3"\tdel"}'; done >> SVs_father_$gv.bed
    for chrom in {1..22}; do grep ^$chrom SVs_mother_$gv.bed | grep del | tail -n2 | awk '$3-$2+1 >= 5000000 {print $1"\t"$2"\t"$3"\tdup"}'; done >> SVs_father_$gv.bed

    sort -k1,1V -k2,2n SVs_father_$gv.bed > bla; mv bla SVs_father_$gv.bed

    # Merging parents BEDs
    grep -E 'del|dup' SVs_*.bed | sort | uniq | cut -f2 -d: | sort -k1V -k2n | bedtools merge -i - -c 4,4 -o count,collapse > parents_CNVs_collapsed.bed
    awk '$5=="del,del"                     {print $1"\t"$2"\t"$3"\tcn=0\t0\t*\t"$2"\t"$3"\t128,0,0"}' parents_CNVs_collapsed.bed > parents_zerocopies.bed
    awk '$5=="del" && ($1=="X" || $1=="Y") {print $1"\t"$2"\t"$3"\tcn=1\t0\t*\t"$2"\t"$3"\t255,0,0"}' parents_CNVs_collapsed.bed >> parents_zerocopies.bed
    awk '$5=="del" &&  $1!="X" && $1!="Y"  {print $1"\t"$2"\t"$3"\tcn=1\t0\t*\t"$2"\t"$3"\t255,0,0"}' parents_CNVs_collapsed.bed > parents_onecopy.bed
    # Although chroms. X and Y have only one copy, the dups. in these chroms. were categorized as 3 copies, just to simplify
    awk '$5=="dup"                         {print $1"\t"$2"\t"$3"\tcn=3\t0\t*\t"$2"\t"$3"\t0,0,255"}' parents_CNVs_collapsed.bed > parents_threecopies.bed
    awk '$5=="dup,dup"                     {print $1"\t"$2"\t"$3"\tcn=4\t0\t*\t"$2"\t"$3"\t0,0,128"}' parents_CNVs_collapsed.bed > parents_fourcopies.bed
    awk '$5=="del,dup" || $5=="dup,del" {print $1"\t"$2"\t"$3"\tLOH\t0\t*\t"$2"\t"$3"\t128,128,128"}' parents_CNVs_collapsed.bed > parents_LOH.bed
    # Putting CNVs together again
    cat parents_*cop*bed parents_LOH.bed | sort -k1,1V -k2n > s$sample""_all_parents_CNVs.bed
    cat parents_zerocopies.bed  parents_onecopy.bed    | sort -k1,1V -k2n > s$sample""_all_parents_dels.bed
    cat parents_threecopies.bed parents_fourcopies.bed | sort -k1,1V -k2n > s$sample""_all_parents_dups.bed
    cd ..
done

# cat SVs_father_$gv.bed SVs_mother_$gv.bed > parents_SV_$gv.bed

# Inserting SNVs
for chrom in {1..22} X Y; do
    echo "python $SVGen_dir/insert_SNVs_select_indels.py \
            --fasta_input $HOME/references/$gv/chr$chrom.fa \
            --fasta_output mother_$gv""_chr$chrom.fa \
            --freq_file $ref_dir/annovar/$gv""_$pop_mother"".sites.2015_08.chrom$chrom.txt \
            --chrom $chrom \
            --vcf_output simulated_SNVs_mother_chr$chrom.vcf" | qsub -V -cwd -l h_vmem=2G -N ins_SNV_mot_$gv""_chr$chrom
    echo "python $SVGen_dir/insert_SNVs_select_indels.py \
            --fasta_input $HOME/references/$gv/chr$chrom.fa \
            --fasta_output father_$gv""_chr$chrom.fa \
            --freq_file $ref_dir/annovar/$gv""_$pop_father"".sites.2015_08.chrom$chrom.txt \
            --chrom $chrom \
            --vcf_output simulated_SNVs_father_chr$chrom.vcf" | qsub -V -cwd -l h_vmem=2G -N ins_SNV_fat_$gv""_chr$chrom
done

rm *mother*Y*

# Inserting SVs
for chrom in {1..22} X Y; do
    echo "python $SVGen_dir/insert_SVs_and_indels.py \
            -i mother_$gv""_chr$chrom.fa \
            -o mother_$gv""_chr$chrom""_SVs.fa \
            --chrom_lens reference/chrom_lengths_$gv.txt \
            --chrom $chrom \
            --bed SVs_mother_$gv.bed \
            -v" | qsub -cwd -V -N mother.simSV_""$chrom
    echo "python $SVGen_dir/insert_SVs_and_indels.py \
            -i father_$gv""_chr$chrom.fa \
            -o father_$gv""_chr$chrom""_SVs.fa \
            --chrom_lens reference/chrom_lengths_$gv.txt \
            --chrom $chrom \
            --bed SVs_father_$gv.bed \
            -v" | qsub -cwd -V -N father.simSV_""$chrom
done

# simulating short reads
# paired-end
mkdir -p short_read/paired_end
cov=20
p=1 # threads
rm sample_*/short_read/paired_end/*
for sample in {1..10}; do
    cd sample_$sample/short_read/paired_end
    for chrom in {1..22}; do
        echo "python $SVGen_dir/create_reads.py -pe --read_label mother -i ../../mother_"$gv"_chr"$chrom"_SVs.fa -o mother_"$gv"_chr"$chrom"_SVs.fq --cov $(($cov/2)) --read_len 100 --snp_rate 0.01 --del_rate 0.0001 --ins_rate 0.0001 --fast_sample 100 -v" > s$sample.mo.pe.$gv.$chrom.sh
        echo "bwa mem -t $p /home/llima/references/"$gv"/chr"$chrom".fa mother_"$gv"_chr"$chrom"_SVs1.fq mother_"$gv"_chr"$chrom"_SVs2.fq | samtools view - -Sb > mother_"$gv"_chr"$chrom"_SVs.bam" >> s$sample.mo.pe.$gv.$chrom.sh
        echo "rm mother_"$gv"_chr"$chrom"_SVs1.fq mother_"$gv"_chr"$chrom"_SVs2.fq" >> s$sample.mo.pe.$gv.$chrom.sh
        echo "samtools sort mother_"$gv"_chr"$chrom"_SVs.bam -@ $p -T mother_"$gv"_chr"$chrom"_SVs.temp -O bam -o mother_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.mo.pe.$gv.$chrom.sh
        echo "samtools index mother_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.mo.pe.$gv.$chrom.sh
        qsub -l h_vmem=8G -pe smp $p -cwd -V s$sample.mo.pe.$gv.$chrom.sh
        echo "python $SVGen_dir/create_reads.py -pe --read_label father -i ../../father_"$gv"_chr"$chrom"_SVs.fa -o father_"$gv"_chr"$chrom"_SVs.fq --cov $(($cov/2)) --read_len 100 --snp_rate 0.01 --del_rate 0.0001 --ins_rate 0.0001 --fast_sample 100 -v" > s$sample.fa.pe.$gv.$chrom.sh
        echo "bwa mem -t $p /home/llima/references/"$gv"/chr"$chrom".fa father_"$gv"_chr"$chrom"_SVs1.fq father_"$gv"_chr"$chrom"_SVs2.fq | samtools view - -Sb > father_"$gv"_chr"$chrom"_SVs.bam" >> s$sample.fa.pe.$gv.$chrom.sh
        echo "rm father_"$gv"_chr"$chrom"_SVs1.fq father_"$gv"_chr"$chrom"_SVs2.fq" >> s$sample.fa.pe.$gv.$chrom.sh
        echo "samtools sort father_"$gv"_chr"$chrom"_SVs.bam -@ $p -T father_"$gv"_chr"$chrom"_SVs.temp -O bam -o father_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.fa.pe.$gv.$chrom.sh
        echo "samtools index father_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.fa.pe.$gv.$chrom.sh
        qsub -l h_vmem=8G -pe smp $p -cwd -V s$sample.fa.pe.$gv.$chrom.sh
        echo "samtools merge -@ $p child_"$gv"_chr"$chrom"_SVs.sort.bam father_"$gv"_chr"$chrom"_SVs.sort.bam mother_"$gv"_chr"$chrom"_SVs.sort.bam" > s$sample.ch.pe.$gv.$chrom.sh
        echo "samtools index child_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.ch.pe.$gv.$chrom.sh
        echo "rm father_"$gv"_chr"$chrom"_SVs.sort.bam* mother_"$gv"_chr"$chrom"_SVs.sort.bam*" >> s$sample.ch.pe.$gv.$chrom.sh
        qsub -l h_vmem=8G -pe smp $p -cwd -V -hold_jid s$sample.fa.pe.$gv.$chrom.sh -hold_jid s$sample.mo.pe.$gv.$chrom.sh s$sample.ch.pe.$gv.$chrom.sh
    done
    chrom=X
    echo "python $SVGen_dir/create_reads.py -pe --read_label mother -i ../../mother_"$gv"_chr"$chrom"_SVs.fa -o mother_"$gv"_chr"$chrom"_SVs.fq --cov $cov --read_len 100 --snp_rate 0.01 --del_rate 0.0001 --ins_rate 0.0001 --fast_sample 100 -v" > s$sample.mo.pe.$gv.$chrom.sh
    echo "bwa mem -t $p /home/llima/references/"$gv"/chr"$chrom".fa mother_"$gv"_chr"$chrom"_SVs1.fq mother_"$gv"_chr"$chrom"_SVs2.fq | samtools view - -Sb > mother_"$gv"_chr"$chrom"_SVs.bam" >> s$sample.mo.pe.$gv.$chrom.sh
    echo "rm mother_"$gv"_chr"$chrom"_SVs1.fq mother_"$gv"_chr"$chrom"_SVs2.fq" >> s$sample.mo.pe.$gv.$chrom.sh
    echo "samtools sort mother_"$gv"_chr"$chrom"_SVs.bam -@ $p -T mother_"$gv"_chr"$chrom"_SVs.temp -O bam -o child_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.mo.pe.$gv.$chrom.sh
    echo "rm mother_"$gv"_chr"$chrom"_SVs.bam" >> s$sample.mo.pe.$gv.$chrom.sh
    echo "samtools index child_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.mo.pe.$gv.$chrom.sh
    qsub -l h_vmem=8G -pe smp $p -cwd -V s$sample.mo.pe.$gv.$chrom.sh
    chrom=Y
    echo "python $SVGen_dir/create_reads.py -pe --read_label father -i ../../father_"$gv"_chr"$chrom"_SVs.fa -o father_"$gv"_chr"$chrom"_SVs.fq --cov $cov --read_len 100 --snp_rate 0.01 --del_rate 0.0001 --ins_rate 0.0001 --fast_sample 100 -v" > s$sample.fa.pe.$gv.$chrom.sh
    echo "bwa mem -t $p /home/llima/references/"$gv"/chr"$chrom".fa father_"$gv"_chr"$chrom"_SVs1.fq father_"$gv"_chr"$chrom"_SVs2.fq | samtools view - -Sb > father_"$gv"_chr"$chrom"_SVs.bam" >> s$sample.fa.pe.$gv.$chrom.sh
    echo "rm father_"$gv"_chr"$chrom"_SVs1.fq father_"$gv"_chr"$chrom"_SVs2.fq" >> s$sample.fa.pe.$gv.$chrom.sh
    echo "samtools sort father_"$gv"_chr"$chrom"_SVs.bam -@ $p -T father_"$gv"_chr"$chrom"_SVs.temp -O bam -o child_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.fa.pe.$gv.$chrom.sh
    echo "samtools index child_"$gv"_chr"$chrom"_SVs.sort.bam" >> s$sample.fa.pe.$gv.$chrom.sh
    qsub -l h_vmem=4G -pe smp $p -cwd -V s$sample.fa.pe.$gv.$chrom.sh
    cd ../../../
done
