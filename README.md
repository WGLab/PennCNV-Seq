# PennCNV-Seq
Adaptation of the original PennCNV algorithm for whole-genome sequencing data.

## PennCNV-Seq installation

### Install PennCNV (at least version v1.0.4)

   Releases can be found [here](https://github.com/WGLab/PennCNV/releases).
   Installation instructions at [http://penncnv.openbioinformatics.org/en/latest/user-guide/install/](http://penncnv.openbioinformatics.org/en/latest/user-guide/install/)
   

### Install BEDtools

   Instructions at [http://bedtools.readthedocs.io/en/latest/content/installation.html](http://bedtools.readthedocs.io/en/latest/content/installation.html)
   
### Download PennCNV-Seq scripts

	git clone git@github.com:WGLab/PennCNV-Seq.git
	cd PennCNV-Seq
	
	# Define your download options (run the command below to see the parameters)
	./download_and_format_database.sh
	

## Run PennCNV-Seq example

	cd PennCNV-Seq
	./penncnv-seq_example.sh [penncnv_dir] [penncnv_ref_dir] [genome_version] [population] [reference.fasta] [bam_file]
	

### Visualize PennCNV results in genome browsers (IGV or UCSC website)

	python penncnv2bed.py results.rawcnv > results.bed
	
