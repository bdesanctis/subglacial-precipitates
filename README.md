Note: very much still in progress

# The ancient subglacial precipitate microbiome

Supplementary materials, results files and code for subglacial precipitate metagenomics paper. 

Many intermediate and output data files are hosted at (add link) and code scripts are on this github in the code folder. 
Note that there is also raw data for a new sequencing round, new_subglacial_raw_data_250730, and associated metadata in the same folder, which has not been processed yet. 

Raw data is on the erda under fq/raw. In the folders excluding new_subglacial_raw_data_250730, we have 1.12 billion raw reads across the 30 non-control libraries covering 26 unique samples, and 61 million raw reads across the 5 controls. From the raw files, I merged lanes, then ran fastp and sga for QC (see code/1-qc.sh). I then removed reads which mapped against some common contaminants, and merged libraries from the same rock in cases where it made sense. The resulting "clean" fastqs are in fq/mapped_to_contams_merged. 

Functional annotation and MAG creation results are also in the erda, in the functional_annotation folder. Methods are described in the paper and code is in this github code folder (see code/7-idba-kegg.sh and paper methods). 

The "clean" fastqs were then mapped against gtdb v226.0, which includes 715230 bacteria and and 17245 archaea genomes, concatenated into 15 fastas of up to 64GB each, totalling 488GB of fastas. This database and its bowtie2 indices and accessory files can be found here https://sid.erda.dk/cgi-sid/ls.py?share_id=iFxbh56MuJ&current_dir=gtdb226&flags=f . There is further database information here https://gtdb.ecogenomic.org/stats/r226 . Accessory files were made with [gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump). Bams for each sample were merged and read-sorted and run through ngsLCA (see code/). These output bams and lca files are in ... [fill in] 


# to do

put up all the combined tsvs (includes taxa with 50 total reads or more) and also individual tsvs and subs files (down to 5 reads per sample)
put up the new kronas - to host publicly soon 
put up per-phyla abundance matrix, and then again for just the subglacial phyla renormalized
put up the nmds results, both loadings and per-sample nmds points
metadata when it's all done

