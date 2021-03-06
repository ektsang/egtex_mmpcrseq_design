## To keep track of order in which to run scripts

Setup
=====

Prerequisite software (version used)
------------------------------------

* Bedtools (v2.25.0)
* Vcftools (UNKNOWN)


Variables used throughout the pipeline
--------------------------------------

To make these accessible to all scripts, add the following paths to `.bashrc`.  
Adjust the paths for your particular setup.

### GTEx directories/files
```
export GTEXv6=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12
export GTEXv7=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7
export GTEX_WGSv7=${GTEXv7}/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.PIR.vcf.gz # used the version without PIR when selecting sites
export GTEX_WGSv7vep=${GTEXv7}/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Indiv_GATK_HaplotypeCaller.siteOnly.vep.vcf.gz
export GTEX_EXOMEv7=${GTEXv7}/genotypes/WES/variant_calls/GTEx_Analysis_2016-01-15_v7_ExomeSeq_603Indiv_GATK_HaplotypeCaller.vcf.gz
export GTEX_OMNIv6=${GTEXv6}/genotypes/OMNI_arrays/variant_calls/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz
export GTEX_RNAv7=${GTEXv7}/rna-seq
export GTEX_SAMPLESv7=${GTEXv7}/sample_annotations/GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt
export GTEX_SUBJECTSv7=${GTEXv7}/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt
```

### Other files and directories
```
export EGTEXDIR=/srv/scratch/etsang/projects/egtex_mmPCR_2015
export GENCODE19=/mnt/lab_data/montgomery/shared/annotations/gencode.v19.annotation.gtf
export BLASTDIR=/mnt/lab_data/montgomery/etsang/BLAST/gencode19
export GENOMEDIR=/mnt/lab_data/montgomery/etsang/genome
export KGDIR=/mnt/lab_data/montgomery/shared/1KG
export BIAS=/mnt/lab_data/montgomery/etsang/ASE/gem_biased/EUR01_50bp_result_stats_05bias.txt
export MAPPABILITY=/mnt/lab_data/montgomery/shared/mappability/wgEncodeCrgMapabilityAlign75mer.bedGraph
```
`EGTEXDIR` is the directory under which all the results will go  
`GENCODE19` is from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz  
`BLASTDIR` is populated with files generated by running `primer_design/make_blast_db.sh` (see below in Primer design section)  
`GENOMEDIR` must contain hg19.fa, pipeline writes here so make a link to hg19.fa in desired directory where masked genome will be written  
`KGDIR` is read only and must contain the 1000 Genomes phase 3v5a VCFs  
`BIAS` is from ftp://jungle.unige.ch/Allelic_map_bias/single_end_transcriptome_based_GEM/  
`MAPPABILITY` is from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/, turned into bedgraph by UCSC's bigWigToBedGraph tool  

Directory structure
-------------------
```
$EGTEXDIR
	figures
	geneSelection
	primerDesign
	sampleInfo
		genotypes
	selectedSites
```

Directory setup
---------------

Make the following sub directories to hold results.

```
mkdir ${EGTEXDIR}/figures
mkdir ${EGTEXDIR}/geneSelection
mkdir ${EGTEXDIR}/sampleInfo
mkdir ${EGTEXDIR}/selectedSites
mkdir ${EGTEXDIR}/primerDesign
```

Necessary input files
---------------------

* master.platemap.txt (move it or link to it in `${EGTEXDIR}/sampleInfo`)


Site selection
==============

Sample information and genotypes
--------------------------------

### Extract subject/sample information
```
bash experimental_layout/process_platemap.sh
bash experimental_layout/subset_sample_annotations.sh
bash experimental_layout/get_individual_info_86inds.sh
```
These scripts generate several files that are used downstream.

### Subset GTEx vcf files to the individuals in eGTEX
```
bash experimental_layout/retrieve_genotypes.sh &> experimental_layout/retrieve_genotypes.log
```
Takes a few hours.

### Process the vcf(s) into bed files with the relevant info
```
bash experimental_layout/retrieve_sites.sh &> experimental_layout/retrieve_sites.log
```
Relies on `experimental_layout/retrieve_sites.py`.
Can take several hours (depending on input file size).

Gene Selection
--------------
The following steps rely on a bunch of files that aren't provided here.
Instead, the final output files are given alongside the code to generate them.

### Get list of cancer genes based on several panels (run on laptop)
```
Rscript gene_selection/select_cancer_genes.R
```
See `gene_selection/readme.cancer.txt`.

### Get set of eQTL genes
```
Rscript gene_selection/select_eqtl_genes.R &> gene_selection/select_eqtl_genes.log
```
Among other files, produces `selected.genes.v2.txt`, which is provided under `gene_selection`. See note below.

### Combine gene lists from different sources
```
Rscript gene_selection/combine_gene_lists.R &> gene_selection/combine_gene_lists.log
```
This generates `selected.genes.txt`, which is provided under `gene_selection`. See note below.

**Note**: For subsequent code to run, link to both selected.genes files from `${EGTEXDIR}/geneSelection`.

Site selection within chosen genes
----------------------------------

### Process gtf file into information for the selected genes
```
bash experimental_layout/retrieve_genes.sh &> experimental_layout/retrieve_genes.log
```
Relies on `experimental_layout/retrieve_genes.py`.
Takes several minutes but should be < 30.

### Process transcript information for the selected genes
```
Rscript experimental_layout/process_transcripts.R
Rscript experimental_layout/process_transcripts.R .v2
```
Data can be reloaded from `experimental_layout/process_transcripts.[v2].RData` once run once.
Takes several minutes but should be < 1h.

### Intersect the created gene and variant files and select variants
```
bash experimental_layout/get_sites_in_exons.sh &> experimental_layout/get_sites_in_exons.log
```
Relies on `experimental_layout/get_sites_in_exons.py`.  
Afterwards did some manual curation of the sites.
For some of this manual curation, used `experimental_layout/print_close_sites.py`.

Primer Design
-------------

### Generate annotation files for primer design with ordered transcripts and exons
```
python primer_design/produce_ordered_annotation.py primer_design/UCSC_TableBrowser_gencode19 \
       ${EGTEXDIR}/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt \
       primer_design/UCSC_TableBrowser_gencode19
```

### Make blast database and link it to the primer_design directory
```
bash primer_design/make_blast_db.sh
ln -s ${BLASTDIR}/hg19.gencode19.blast.db* primer_design/
```

### Process and mask genome
```
bash primer_design/mask_genome.sh &> primer_design/mask_genome.log
```
Relies on `primer_design/split_sequences.py`.
Takes several hours to run.

### Design primer pools
#### First pass, run in groups of 100
```
firstPassFile=${EGTEXDIR}/primerDesign/selected.variants.firstpass.master.txt
tail -n +2 ${EGTEXDIR}/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt | \
     awk 'BEGIN{FS=OFS="\t"}{print $1,$5,$7,$7,$4}' | sort -R > $firstPassFile
bash primer_design/primer_design_pass.sh firstpass 11 0 $firstPassFile
```

#### Combine for second pass and runa second pass with all remaining primers
```
bash primer_design/cleanup_pools_prepare_next_pass.sh firstpass secondpass
mkdir ${EGTEXDIR}/primerDesign/secondpass
python primer_design/run_yamRTPCR_cDNA.py \
       --poolSize 11 \
       --numberPools 96 \
       --geneAnnotation UCSC_TableBrowser_gencode19_selected_exons UCSC_TableBrowser_gencode19_selected_transcripts \
       --existingPools ${EGTEXDIR}/primerDesign/firstpass.pools.txt \
       --directory ${EGTEXDIR}/primerDesign/secondpass \
       --inputSiteFile ${EGTEXDIR}/primerDesign/selected.variants.secondpass.master.txt \
       &> ${EGTEXDIR}/primerDesign/secondpass/log.txt
```

#### Get additional variants to make additional pools
Make new set of variants (see `gene_selection/select_eqtl_genes.R`). 
Then ran variants of the pipeline above.  
Pruned site file is `${EGTEXDIR}/selectedSites/genes.v2.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt`.

#### Make new annotation file including the new sites
```
python primer_design/produce_ordered_annotation.py primer_design/UCSC_TableBrowser_gencode19 \
       ${EGTEXDIR}/selectedSites/genes.v2.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt \
       primer_design/UCSC_TableBrowser_gencode19_V2
```

#### Run third pass using these new variants
```
mkdir ${EGTEXDIR}/primerDesign/thirdpass
thirdPassFile=${EGTEXDIR}/primerDesign/selected.variants.thirdpass.master.txt
tail -n +2 ${EGTEXDIR}/selectedSites/genes.v2.intersect.WGS79.OMNI86.selected.variants.ordered.firstDesign.txt | \
     awk 'BEGIN{FS=OFS="\t"}{print $1,$5,$7,$7,$4}' | sort -R > $thirdPassFile
python primer_design/run_yamRTPCR_cDNA.py \
       --poolSize 11 \
       --numberPools 5 \
       --counterStart 104 \
       --geneAnnotation UCSC_TableBrowser_gencode19_V2_selected_exons UCSC_TableBrowser_gencode19_V2_selected_transcripts \
       --directory ${EGTEXDIR}/primerDesign/thirdpass \
       --inputSiteFile $thirdPassFile &> ${EGTEXDIR}/primerDesign/thirdpass/log.txt
```

#### Combine completed pools
```
cat ${EGTEXDIR}/primerDesign/secondpass/pools.txt ${EGTEXDIR}/primerDesign/thirdpass/pools.txt \
    > ${EGTEXDIR}/primerDesign/almostFinalPools.txt
```

Identified site where the primers mapped to a different gene:
```
cat ${EGTEXDIR}/primerDesign/almost.final.pools.txt | tr '_' ' ' | awk '$2!=$8'
$ Pool5 STAT4 4  ACTGAGAACAAAGCATTGTAATGT	59.83	GCATTTCTCCTCCCACTTGA	60.39	VGLL1 ENST00000370634.3:563-967	404	UCSC TableBrowser gencode19 selected exons
```
Manually removed that entry and replaced it.
```
mkdir ${EGTEXDIR}/primerDesign/fourthpass
python primer_design/run_yamRTPCR_cDNA.py \
       --poolSize 11 \
       --numberPools 96 \
       --geneAnnotation UCSC_TableBrowser_gencode19_V2_selected_exons UCSC_TableBrowser_gencode19_V2_selected_transcripts \
       --existingPools ${EGTEXDIR}/primerDesign/almost.final.pools.txt \
       --directory ${EGTEXDIR}/primerDesign/fourthpass \
       --inputSiteFile ${EGTEXDIR}/primerDesign/thirdpass/remaining.sites.txt \
       &> ${EGTEXDIR}/primerDesign/fourthpass/log.txt
```

#### Prepare designed primers to be ordered
Link to completed pools and add adapters
```
ln -s ${EGTEXDIR}/primerDesign/fourthpass/pools.txt ${EGTEXDIR}/primerDesign/final.pools.txt
python primer_design/add_adaptor.py ${EGTEXDIR}/primerDesign/final.pools.txt | \
       sort -V > ${EGTEXDIR}/primerDesign/final.pools.adaptors.txt
```

Make spreadsheet (this outputs a format similar to what IDT requires and also makes a site info file for all the sites with primers designed).
```
Rscript primer_design/format_pools_into_order.R
```

#### Extract amplicon positions
```
python primer_design/get_amplicon_coordinates.py \
       ${EGTEXDIR}/selectedSites/genes.intersect.WGS79.OMNI86.selected.variants.primerpools.txt
       ${EGTEXDIR}/primerDesign/final.pools.txt primer_design/UCSC_TableBrowser_gencode19.bed
```
