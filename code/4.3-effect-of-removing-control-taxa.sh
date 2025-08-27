#!/bin/bash
# Directories for the two file sets
DIR1="bamdam-gtdb"
DIR2="bamdam-gtdb-with-control-taxa"

# Print header
printf "Sample\tLCAlines\tLCAlineswithcontroltax\tTSVlines\tTSVlineswithcontroltaxa\n"

# Loop over the LCA files in the first folder (assumed to end in _filtered.lca)
for file in ${DIR1}/*_filtered.lca; do
    # Extract sample name by stripping the suffix _filtered.lca
    sample=$(basename "$file" _filtered.lca)

    # Count lines for LCA file in DIR1
    lca1=$(wc -l < "$file")

    # TSV file in DIR1 (assumed to be sample.tsv)
    tsv1_file="${DIR1}/${sample}.tsv"
    tsv1=$( [ -f "$tsv1_file" ] && wc -l < "$tsv1_file" || echo "NA" )

    # In DIR2 the filenames usually have _includecontroltax, except for S14 and S28.
    if [ -f "${DIR2}/${sample}_filtered_includecontroltax.lca" ]; then
         lca2_file="${DIR2}/${sample}_filtered_includecontroltax.lca"
         tsv2_file="${DIR2}/${sample}_includecontroltax.tsv"
    else
         lca2_file="${DIR2}/${sample}_filtered.lca"
         tsv2_file="${DIR2}/${sample}.tsv"
    fi

    lca2=$( [ -f "$lca2_file" ] && wc -l < "$lca2_file" || echo "NA" )
    tsv2=$( [ -f "$tsv2_file" ] && wc -l < "$tsv2_file" || echo "NA" )

    # Print out the table row
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$lca1" "$lca2" "$tsv1" "$tsv2"
done
