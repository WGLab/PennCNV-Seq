#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use File::Basename;

our ($verbose, $help, $man);
our ($bamfile, $reffile);
our ($outfile, $smooth, $mapq_threshold, $freqfile, $step, $region);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'smooth=f'=>\$smooth, 'mapq_threshold=i'=>\$mapq_threshold,
	'freqfile=s'=>\$freqfile, 'step=s'=>\$step, 'region=s'=>\$region) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($bamfile, $reffile) = @ARGV;
$outfile ||= $bamfile;
$smooth ||= 0.5;
$mapq_threshold ||= 20;
$region ||= '';

my %valistep;
my $meancov;

if ($step) {
	my @step = split (/,/, $step);
	for my $nextstep (@step) {
		$nextstep =~ m/^\d+$/ or pod2usage ("Erron in argument: the -step argument should be comma-delimited numbers");
		$valistep{$nextstep}++;
	}
}

if ($valistep{1} or not $step) {
	readVariantInfo ($bamfile, $reffile, "$outfile.read1");
}
#my ($coverage, $region) = convertBam2Read ($bamfile, $reffile, "$outfile.read1");		#generate READ file
#convertRead2Signal ("$outfile.read1", "$outfile.read2", $freqfile, $coverage, $region);		#generate LRR/BAF file

if ($valistep{2} or not $step) {
	addPairedDistance ($bamfile, "$outfile.read1", "$outfile.read2");
}

if ($valistep{3} or not $step) {
	addLRRBAF ("$outfile.read2", "$outfile.read3", $meancov);
}




sub readVariantInfo {
	my ($bamfile, $reffile, $readoutfile) = @_;
	my $command;
	
	if ($region) {
		(-f "$bamfile.bai" ||-f "$bamfile.crai") or die "Error: for -r argument, the BAM/CRAM file must be indexed first\n";
		$command = "bcftools mpileup -Ou -f $reffile -r $region $bamfile  | bcftools call -c - |";
	} else {
		$command = "bcftools mpileup -Ou -f $reffile $bamfile $region | bcftools call -c - |";
	}
	#$command = 'temp.vcf';				#for debugging purposes
	
	print STDERR "NOTICE: Running samtools <$command>\n";
	
	my ($countsite, $countcov, $countsnp, $sumcov) = (0, 0, 0, 0);		#count of sites examined, sum of coverage for all bases in a region
	my ($prechr, $prepos, $precov) = (-1);
	
	print STDERR "NOTICE: Starting variant calling and processing variant information from BAM/CRAM file (input=$bamfile output=$readoutfile)\n";
	open (VAR, $command) or die "Error: cannot read from command output: <$command>\n";
	open (OUT, ">$readoutfile") or die "Error: cannot write to output file $readoutfile: $!\n";
	print OUT "Name\tCoverage\tBAC\tBAF\tLength\n";
	while (<VAR>) {
		m/^#/ and next;
		
		my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split (/\t/, $_);
		my ($cov, $bac);			#depth coverage of current position, alternative allele count at current position

		$info =~ m/^INDEL/ and next;		#do not consider indels for the moment (for samtools mpileup)
		if ($info =~ m/DP=(\d+)/) {
			$cov = $1;
		} else {
			die "Error in VCF record (DP not found): <$_>\n";
		}
		$countsite++;				#count number of sites examined
		$countcov += $cov;

		if ($prechr ne $chr) {			#encountered a new chromosome
			$prechr = $chr;
			$precov = 0;
			$prepos = $pos;
			$sumcov = 0;
		}
		
		if ($alt ne '.') {			#SNP marker
			$countsnp++;
			if ($info =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/) {
				$bac = $3+$4;
			} else {
				die "Error in VCF record (DP4 not found): <$_>\n";
			}
			if ($prepos < $pos) {		#if last record is less than the current position
				my $mean_cov = sprintf ("%.2f", $sumcov / ($pos-$prepos));	#average coverage from the last record to the previous position
				print OUT $chr, ":", $prepos, "-", $pos-1, "\t", $mean_cov, "\t", 0, "\t", 2, "\t", $pos-$prepos, "\n";	#chr:start-end coverage BAC length
			}
			print OUT $chr, ":", $pos, "-", $pos, "\t", $cov, "\t", $bac, "\t", sprintf ("%.2f", $bac/$cov), "\t", 1, "\n";	#chr:start-end coverage BAC length
			$prepos = $pos+1;		#move the prepos point to the next position to handle writing two SNPs in a row
			$precov = 0;
			$sumcov = 0;
		} else {				#segment marker
			if ($precov==0 and $cov==0) {	#two consecutive marker with DP=0
				next;
			} elsif ($precov == 0 and $cov or abs ($cov - $precov)/$precov > $smooth) {		#coverage fractional chagne is more than the smooth threshold
				if ($prepos < $pos) {
					my $mean_cov = sprintf ("%.2f", $sumcov / ($pos-$prepos));
					print OUT $chr, ":", $prepos, "-", $pos-1, "\t", $mean_cov, "\t", 0, "\t", 2, "\t", $pos-$prepos, "\n";	#chr:start-end coverage BAC length
				}
				#print STDERR "Found smooth checkpoint prepos=$prepos pos=$pos sumcov=$sumcov cov=$cov precov=$precov\n";
				$prepos = $pos;
				$precov = $cov;
				$sumcov = $cov;
				#print STDERR "Assigned smooth checkpoint prepos=$prepos pos=$pos sumcov=$sumcov cov=$cov precov=$precov\n";
			} else {
				$sumcov += $cov;
				#print STDERR "pos=$pos sumcov=$sumcov precov=$precov cov=$cov\n";
			}
		}
		#$countsite > 20 and last;
	}
	$meancov = $countcov/$countsite;
	print STDERR "NOTICE: Finished examining $countsite sites (including $countsnp SNPs) with mean coverage $meancov\n";
}









=head1
sub convertBam2Read {
	my ($bamfile, $reffile, $read1file) = @_;
	print STDERR "NOTICE: Starting step 1: generating Coverage, B_Allele_Count and Length information\n";
	
	my ($countline) = (0);		#number of lines in mpileup file
	my %region;			#key=chr value=1 (store chromosomes that have been read in the BAM file)
	
	my (@field);
	my ($prechr, $prepos, $precov) = (-1);
	
	my ($command);
	my ($sumcov) = (0);		#sum of coverage for a sliding window

	$command = "samtools mpileup -f $reffile $bamfile |";

	my ($bpcount, $bpcovcount, $bplncovcount, $bpcovcount2);
	print STDERR "NOTICE: Collecting and analyzing output from: <$command>\n";
	open (FH, $command) or die "Error: cannot open system command: <$command>\n";
	open (OUT, ">$read1file") or die "Error: cannot write to temporary file $read1file: $!\n";
	while (<FH>) {
		#2       10074   A       0               
		#2       10211   T       8       ..,,,,.,        ;?26:=6:
		#2       10212   A       9       ....,.,.+2AT,   6;=:?.23;
		#2       10213   C       9       ..,,.,+1a,^:.^8.        =<64=57>;
		#2       10214   C       11      ...,,,.,,+1a..  1><73/>:72<
		#2       10215   C       13      .$...,,,.,,..^-,        24>A1.2/4/=;.
		#2       10216   G       8       ..,.,..,        44207429
		#2       10217   A       18      ...,,,,.,.,..,,^:.^$.^!,        3;94E<;=725<99>=<>
		#2       10218   A       15      ...,,,.,..,..,^!,       <>=@856>>;:>??5
		#2       10219   C       17      ...,,,,,,.,,.,,,^!,     0<1?84::::;1:;180
		#2       10220   C       21      ...,,,,.,,..,,.,,,,^!.^!,       ;=3?:;:>0?=;:15;0:7<1
		#2       10221   C       21      ..,,,,,a,..,..,,,,.,^!, =A9<<</2;=;<81>6:7;:7
		#2       10222   T       20      ...,,,.,,..,..,,,.,,    1;;<46?43?>=>6@6:@54
		#2       10223   A       25      ...,,,,.,.,..,,..,,,,.,,,       4;;6C<<:A2=<;:>13A6=8=:@7
		chomp;
		$countline++;				#number of lines in the pileup file
		#@field = split (/\t/, $_);              #chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities
		my ($chr, $pos, $ref, $cov, $read, $qual) = split (/\t/, $_);              #chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities

		$region{$chr} ||= 1;		#register this region/chromosome for examination in freq file later on
		
		if ($prechr ne $chr) {	#encountered a new chromosome
			$precov = $cov;
			$prepos = $pos;
			$sumcov = 0;
		}
		$cov or next;			#this position has no read so coverage is zero
		
		#if ($preregion ne $field[0]) {		#encountered new chromosome
		#	$precv = $field[3];
		#	$prepos = $field[1];
		#	$sumcv = 0;
		#}
		#($curcv, $curseq) = @field[3,4];
		#$curcv or next;				#this position has no read
		
		#defined $curseq or die "Error: cannot find defined seq in <$_>\n";
		#$curseq =~ s/^\@//;			#MAQ-generated alignment typically has a heading @, which should be eliminated
		
		if ($curcv) {
			my $same = $curseq =~ tr/,/,/ + $curseq =~ tr/././;
			if ($curcv != $same) {
				#$curseq =~ s/[\+\-]1.//g;
				#$curseq =~ s/[\+\-]2..//g;				#most 2 insertion/deletion are found in MAQ alignment for short reads
				if ($curseq =~ m/[\+\-](\d+)/) {
					while ($curseq =~ m/[\+\-](\d+)/) {
						my $delete = $1;
						$curseq =~ s/[\+\-](\d+).{$delete}//;
					}
				}
				
				my $count_a = $curseq =~ tr/Aa/Aa/;
				my $count_g = $curseq =~ tr/Gg/Gg/;
				my $count_c = $curseq =~ tr/Cc/Cc/;
				my $count_t = $curseq =~ tr/Tt/Tt/;
				my $count_max = $count_a;
				$count_g and $count_max < $count_g and $count_max = $count_g;
				$count_c and $count_max < $count_c and $count_max = $count_c;
				$count_t and $count_max < $count_t and $count_max = $count_t;
				$bacount = $count_max;
			} else {
				$bacount = 0;
			}
			
			if ($bacount > $curcv) {					#for example, when ,+12aataaataaata.,+12aataaataaata,,,,C.,,,,,,,, is present
				$bacount = $curcv-$same;
			}
			
			$baf = 0;
			my $deviation;
			if ($bacount) {
				#if ($bacount/$curcv>0.04 and $bacount/$curcv<0.5) {
				#	$deviation = kc::bitest ($curcv, $bacount, 0.04);
				#	$deviation < 0.05 and $baf = $bacount / $curcv;
				#} elsif ($bacount/$curcv<0.96 and $bacount/$curcv>0.5) {
				#	$deviation = kc::bitest ($curcv, $bacount, 0.96);
				#	$deviation < 0.05 and $baf = $bacount / $curcv;
				#}
				
				if ($bacount/$curcv > 0.2 and $bacount/$curcv < 0.8 and $bacount > 3) {
					$baf = $bacount/$curcv;
				}
			}
		}
		
		if ($curcv and $baf > 0) {			#BIG CHANGE EHRE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if ($prepos < $field[1]) {
				my $mean_cv = $sumcv/($field[1]-$prepos);
				$mean_cv =~ tr/././ and $mean_cv = sprintf("%.2f", $mean_cv);
				print OUT $field[0], ":", $prepos, "-", $field[1]-1, "\t", $mean_cv, "\t", 0, "\t", $field[1]-$prepos, "\n";
			}
			print OUT $field[0], ":", $field[1], "-", $field[1], "\t", $curcv, "\t", $bacount, "\t", 1, "\n";
			$prepos = $field[1]+1;
			$sumcv = 0;
		} else {
			if ($precv == 0 and $curcv or $smooth >= 1 and $curcv - $precv > $smooth || $precv-$curcv > $smooth or $smooth < 1 and $precv and ($curcv - $precv)/$precv > $smooth || ($precv-$curcv)/$precv > $smooth) {
				if ($prepos < $field[1]) {
					my $mean_cv = $sumcv/($field[1]-$prepos);
					$mean_cv =~ tr/././ and $mean_cv = sprintf("%.2f", $mean_cv);
					print OUT $field[0], ":", $prepos, "-", $field[1]-1, "\t", $mean_cv, "\t", 0, "\t", $field[1]-$prepos, "\n";
				}
				$prepos = $field[1];
				$sumcv = $curcv;
			} else {
				$sumcv += $curcv;		#continue adding the cv to the sum
			}
		}

		$precv = $curcv;
		
		$preregion = $field[0];
		$countcv += $curcv;
		$countline % 1_000_000 == 1 and print STDERR "NOTICE: Processing line $countline (chr=$field[0] base=$field[1])\n";
		$countline == $maxline and last;
		
		if ($curcv < 100 and $field[2] ne 'N') {
			$bpcovcount+=$curcv;
			$bplncovcount+=log($curcv || 1e-9);
			$bpcount++;
			$bpcovcount2+=$curcv*$curcv;
		}
		
	}
	
	if (not $countline) {
		die "Error: cannot execute system command (program in path? memory issue?): <$command>\n";
	}
	if ($prepos < $field[1]) {
		my $mean_cv = $sumcv/($field[1]-$prepos);
		$mean_cv =~ tr/././ and $mean_cv = sprintf("%.2f", $mean_cv);
		print OUT $field[0], ":", $prepos, "-", $field[1]-1, "\t", $mean_cv, "\t", 0, "\t", $field[1]-$prepos, "\n";
	}
	
	if ($bpcount<=10) {
		1;
	} else {
	
		my $var2 = ($bpcovcount2 - $bpcovcount*$bpcovcount/$bpcount) / ($bpcount-1);
		my $k2 = ($bpcovcount/ $bpcount)*($bpcovcount/ $bpcount) / $var2;
		my $theta2 = $var2 / ($bpcovcount/ $bpcount);
		
		close (OUT);
		close (FH);
		
		open (LOG, ">$outfile.log") or die "Error: cannot write to log file $outfile.log: $!\n";
		print LOG "NOTICE: Processing total of $countline positions ($bpcount non-N bases) in ", scalar (keys %region), " chromosomes/regions of genome, with overall base count of $countcv, representing ", $countcv/$bpcount, " coverage\n";
		print LOG "NOTICE: Estimated Gamma distribution parameter: mean=", sprintf ("%.3f", $k2*$theta2), " variance=", sprintf ("%.3f", $var2), " k=", sprintf ("%.3f", $k2), " theta=", sprintf ("%.3f", $theta2), "\n";
		close (LOG);
	
		print STDERR "NOTICE: Processing total of $countline positions ($bpcount non-N bases) in ", scalar (keys %region), " chromosomes/regions of genome, with overall base count of $countcv, representing ", sprintf ("%.2f", $countcv/$bpcount), "X coverage\n";
		print STDERR "NOTICE: Estimated Gamma distribution parameter: mean=", sprintf ("%.3f", $k2*$theta2), " variance=", sprintf ("%.3f", $var2), " k=", sprintf ("%.3f", $k2), " theta=", sprintf ("%.3f", $theta2), "\n";
		print STDERR "NOTICE: Read coverage values have been written to $outfile.read file\n";
	}
	
	$coverage = sprintf("%.3g", $countcv/$countline);
	return ($coverage, \%region);
}
=cut


sub addLRRBAF {
	my ($readinfile, $readoutfile, $meancov) = @_;
	$meancov ||= 30;
	print STDERR "NOTICE: Start adding LRR/BAF information to the read file (input=$readinfile output=$readoutfile)\n";
	open (IN, $readinfile) or die "Error: cannot read from inputfile $readinfile: $!\n";
	open (OUT, ">$readoutfile") or die "Error: cannot write to outputfile $readoutfile: $!\n";
	$_ = <IN>;
	chomp;
	m/^Name\tCoverage/ or die "Error: invalid header line in inputfile $readinfile: <$_>\n";
	print OUT $_, "\tLRR\tChr\tPosition\n";
	
	while (<IN>) {
		chomp;
		my ($region, $coverage, $bac, $baf, $length, $ped) = split (/\t/, $_);
		$coverage > 0 or $coverage = 1;	#deal with coverage=0 situation
		$region =~ m/(chr)?(\w+):(\d+)-(\d+)/ or die "Error: invalid region specifier <$region>\n";
		print OUT $_, "\t", log($coverage/$meancov), "\t", "$2\t$3\n";
	}
}
	

sub addPairedDistance {
	my ($bamfile, $readinfile, $readoutfile) = @_;
	print STDERR "NOTICE: Start adding paired end distance information to the read file (input=$readinfile output=$readoutfile)\n";
	open (IN, $readinfile) or die "Error: cannot read from inputfile $readinfile: $!\n";
	open (OUT, ">$readoutfile") or die "Error: cannot write to outputfile $readoutfile: $!\n";
	open (SAM, "bcftools view $bamfile |") or die "Error: cannot read from SAM input\n";
	$_ = <IN>;
	chomp;
	m/^Name\tCoverage/ or die "Error: invalid header line in inputfile $readinfile: <$_>\n";
	print OUT $_, "\t", "PairedEndDistance\n";
	
	my ($prername, $prepos, $pretlen);
	while (<IN>) {
		chomp;
		my ($region, $coverage, $bac, $baf, $length) = split (/\t/, $_);
		$region =~ m/^(chr)?(\w+):(\d+)-(\d+)$/ or die "Error: invalid region specifier in inputfile: <$region>\n";
		my ($chr, $start, $end) = ($2, $3, $4);
		my @dist;		#an array storing query template length (paired end distance) information
		
		if (defined $prername) {
			if ($prername eq $chr and $prepos >= $start and $prepos <= $end) {
				push @dist, $pretlen;
			}
			($prername, $prepos, $pretlen) = (undef, undef, undef);
		}
		
		while (<SAM>) {
			my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split (/\t/, $_);
			$tlen < 0 and $tlen = -$tlen;		#sometimes it is negative since the first read is located after the paired read
			$tlen or next;				#sometimes tlen is zero
			$mapq > $mapq_threshold or next;			#only consider highly reliable alignment
			
			$rname =~ s/^chr//;
			if ($rname ne $chr) {
				next;
			}
			if ($rnext ne '=') {
				next;
			}
			if ($pos < $start) {
				next;
			}
			if ($pos >= $start and $pos <= $end) {
				push @dist, $tlen;
			}
			if ($pos > $end) {
				($prername, $prepos, $pretlen) = ($rname, $pos, $tlen);		#save predist information
				last;		#end this SAM reading cycle
			}
		}
		if (@dist) {
			print OUT join ("\t", $region, $coverage, $bac, $baf, $length, join (",", @dist)), "\n";
		} else {
			print OUT join ("\t", $region, $coverage, $bac, $baf, $length, "."), "\n";
		}
	}
}
		
	

sub convertRead2Signal {
	my ($read1file, $read2file, $freqfile, $coverage, $region) = @_;
	print STDERR "NOTICE: Starting step 2: adding PFB information to read file\n";
	
	my (@record, %snploc, %snppfb);
	
	my ($countline, $preregion, $prepos) = (0);
	
	if (defined $freqfile) {
		open (FREQ, $freqfile) or die "Error: cannot read freqfile $freqfile: $!\n";
		print STDERR "NOTICE: Reading freqfile $freqfile for population frequency of B allele ...";
		while (<FREQ>) {
			@record = split (/\t/, $_);
			$region->{$record[0]} or next;
			push @{$snploc{$record[0]}}, $record[1];
			push @{$snppfb{$record[0]}}, sprintf ("%.3f", $record[2] eq $record[10] ? ($record[9]/($record[9]+$record[11])) : ($record[11]/($record[9]+$record[11])));
			$countline++;
			$prepos ||= $record[1];
			$preregion ||= $record[0];
			$preregion eq $record[0] or $prepos = 0;
			$record[1] >= $prepos or die "Error: the freq file is not sorted ($record[1] occur after $prepos in region $record[0])\n";
			$prepos = $record[1];
			$preregion = $record[0];
		}
		close (FREQ);
		print STDERR " Done with $countline frequency records\n";
	} else {
		print STDERR "NOTICE: No external allele frequency information will be used in signal extraction (you can use --freqfile to supply the information)\n";
	}
	
	my (%index);
	my ($count1, $count2, $count3) = qw/0 0 0/;
	open (READ, $read1file) or die "Error: cannot read from sigfile $read1file\n";
	open (SIG, ">$read2file") or die "Error: cannot write to lrr/baf file $read2file.read2: $!\n";
	print SIG "Name\tCoverage\tB Allele Count\tPFB\tLength\n";
	while (<READ>) {		
		@record = split (/\t/, $_);
		$record[0] =~ m/(\w+):(\d+)-(\d+)/ or die "Error: invalid record found in signalfile $read1file: <$_>\n";
		my ($chr, $start, $end, $curpfb) = ($1, $2, $3);

		
		if ($start == $end) {
			my $curindex = $index{$chr} || 0;
			my $curloc = $snploc{$chr};
			my $flag_add;
			while (defined $curloc->[$curindex] and $curloc->[$curindex] < $start) {
				$curindex++;
				$flag_add++;
			}
			if (defined $curloc->[$curindex] and $curloc->[$curindex] == $start) {
				$curpfb = $snppfb{$chr}->[$curindex];
				$count1++;
			} else {
				if ($record[2]) {
					$curpfb = 0;			
				} else {
					$curpfb = 2;			#this is a sinle-bp CN marker
				}
				$flag_add and $curindex--;
				$count2++;
			}
			$index{$chr} = $curindex;
			
			#if ($freqout) {
			#	my $newbaf = $record[2]/$record[1];
			#	length ($newbaf)>4 and $newbaf = sprintf("%.2f", $newbaf);
			#	print SIG join ("\t", @record[0..1]), "\t", $newbaf, "\t", $curpfb, "\t", $end-$start+1, "\n";
			#} else {
				print SIG join ("\t", @record[0..2]), "\t", $curpfb, "\t", $end-$start+1, "\n";
			#}
		} else {
			$curpfb = 2;
			$count3++;
			print SIG join ("\t", @record[0..2]), "\t", $curpfb, "\t", $end-$start+1, "\n";
		}
		
		
	}
	close (READ);
	close (SIG);
	print STDERR "NOTICE: Finished writting $count1 SNP markers with population frequency information, $count2 SNP markers without population freq info, and $count3 segment markers\n";
	
	open (LOG, ">>$outfile.log") or die "Error: cannot write to log file $outfile.log: $!\n";
	print LOG "NOTICE: Finished writting $count1 SNP markers with population frequency information, $count2 SNP markers without population freq info, and $count3 segment markers\n";
	close (LOG);
}


=head1 SYNOPSIS

 convert_map2signal.pl [arguments] <BAM/CRAMfile> <FASTAfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --outfile <string>		prefix of output files (default=<BAMfile>)
 	    --smooth <int>		smooth level (default=0.5)
 	    --freqfile <file>		a file annotating allele frequency in reference population
 	    --mapq_threshold <int>	mapping quality threshold to use (default=20)
 	    --step <string>		comma-delimited step to use
 	    --region <string>		BAM region to extract signal

 Function: convert genome mapping file for short read sequences to signal 
 intensity values, which can be used for CNV calling by algorithms such as PennCNV

 Example: convert_map2signal.pl chr22.bam chr22.fasta

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=back

=head1 DESCRIPTION

This program is used to convert genome mapping file for short read sequences 
into signal intensity scores, such that these values can be subsequently used in 
follow-up copy number variation analysis. The genome mapping file is typically 
also refered to as alignment file, and they could be in various formats. This 
program works on the so-called pileup format directly, but users can also feed 
other popular format (such as MAQ alignment format), or BAM alignment format, 
and the program can call external programs to decode these map files into pile-
up format and analyze them. The decoding procedure is performed on-the-fly, so 
it does not consume much disk space.

