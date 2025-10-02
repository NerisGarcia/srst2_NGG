#=========================================================================
# SRST2-dev - new database
#=========================================================================
# Using the AMRFinderPlus database as reference:
# /home/ngarcia/miniconda3/envs/genomics_env/share/amrfinderplus/data/2024-07-22.1

# Move to the working directory for the new database
cd /home/ngarcia/RESEARCH/Projects/2024_AMR_Enteroccoco/Goals/srst2_NGG/NGG_db/

# Remove all files in the current directory to start fresh
rm *

# Copy the latest AMR_CDS file from AMRFinderPlus data to the working directory
cp /home/ngarcia/miniconda3/envs/genomics_env/share/amrfinderplus/data/latest/AMR_CDS /home/ngarcia/RESEARCH/Projects/2024_AMR_Enteroccoco/Goals/srst2_NGG/NGG_db

# List files to confirm copy
ls

# Add custom variants: first, modify their names in the FASTA headers
# Example header format explained for reference

# Activate the genomics conda environment
conda activate genomics_env

# Loop through each resistance gene to process its alleles
for gene in "tet(M)" "tet(L)" "tet(S)" "erm(B)" "sat4" "aph(3')-IIIa" "ant(6)-Ia"; do
    
    # Set input and output FASTA filenames for each gene
    input_fasta="/home/ngarcia/RESEARCH/Projects/2024_AMR_Enteroccoco/Goals/Tet_resistance/data/00_tet_variant_catalog/02_final_tet_alleles/${gene}_alleles_nuc.fasta"
    output_fasta="${gene}_alleles_nuc.NGG_db.fasta"
    edited_fasta="${gene}_alleles_nuc.NGG_db.editednames.fasta"
    
    # Edit FASTA headers to insert gene name before AMR_ref_gene_catalog
    sed "s/|AMR_ref_gene_catalog/|${gene}|AMR_ref_gene_catalog/" "$input_fasta" > "$output_fasta"

    # Transform headers to AMRFinderPlus format using a custom Python script
    python ../NGG_code/transform_fasta_headers.py "$output_fasta" -o "$edited_fasta" --amrfplus-fam-tab /home/ngarcia/miniconda3/envs/genomics_env/share/amrfinderplus/data/latest/fam.tab

    # For aph(3')-IIIa, replace special characters in the header for compatibility
    if [ "$gene" == "aph(3')-IIIa" ]; then
        cp "$edited_fasta" "$edited_fasta".temp
        seqkit replace -p "aph\(3'\)-IIIa" -r "aph3_IIIa" "$edited_fasta".temp -o "$edited_fasta"
        mv "$edited_fasta"  "aph3_IIIa_alleles_nuc.NGG_db.editednames.fasta"
        rm "$edited_fasta".temp
    fi

    # For ant(6)-Ia, replace special characters in the header for compatibility
    if [ "$gene" == "ant(6)-Ia" ]; then
        cp "$edited_fasta" "$edited_fasta".temp
        seqkit replace -p "ant\(6\)-Ia" -r "ant6_Ia" "$edited_fasta".temp -o "$edited_fasta"
        mv "$edited_fasta"  "ant6_Ia_alleles_nuc.NGG_db.editednames.fasta"
        rm "$edited_fasta".temp
    fi

    # For tet and erm genes, remove parentheses from headers for compatibility
    if [ "$gene" == "tet(L)" ] || [ "$gene" == "tet(M)" ] || [ "$gene" == "tet(S)" ] || [ "$gene" == "erm(B)" ]; then
        cp "$edited_fasta" "$edited_fasta".temp
        seqkit replace -p "\(|\)" -r "" "$edited_fasta".temp -o "$edited_fasta"
        gene1=$(echo "$gene" | sed 's/[()]//g')
        mv "$edited_fasta"  "${gene1}_alleles_nuc.NGG_db.editednames.fasta"
        rm "$edited_fasta".temp
    fi

    # Remove intermediate output FASTA
    rm "$output_fasta"
    # Optionally, check headers
    #grep ">" "$edited_fasta"
done

# List files to confirm processing
ls

# Count number of sequences in AMR_CDS
grep -c ">" AMR_CDS

# Concatenate AMR_CDS and all processed allele FASTAs into a single file
cat AMR_CDS *_alleles_nuc.NGG_db.editednames.fasta > AMR_CDS_plus_tetLMS.fasta

# Check headers in the combined FASTA
grep ">" AMR_CDS_plus_tetLMS.fasta 

# Convert combined FASTA to SRST2 format using a custom Python script
python ../NGG_code/amrfinder_to_srst2.py AMR_CDS_plus_tetLMS.fasta -o AMR_CDS_plus_tetLMS.def.fasta

# Check headers in the SRST2-formatted FASTA
grep ">" AMR_CDS_plus_tetLMS.def.fasta

# Remove gaps ("-") from all sequences using seqkit
seqkit seq --remove-gaps "AMR_CDS_plus_tetLMS.def.fasta" > "AMR_CDS_plus_tetLMS.def.nogaps.fasta"

# Clean headers: remove parentheses, slashes, apostrophes, and rare symbols, but keep dots and colons
seqkit replace -p '[()/\\'\''|,;[\]{}]' -r '' "AMR_CDS_plus_tetLMS.def.nogaps.fasta" -o "AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta"

# Check headers in the cleaned FASTA
grep ">" "AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta"

# Count number of sequences in the cleaned FASTA
grep -c ">" "AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta"

# Search for specific gene (blaCAU) in the headers, showing context
grep ">" "AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta" | grep -B 10 -A 10 "blaCAU"

# Deactivate genomics environment
conda deactivate
#-----------

# Activate SRST2 development environment
conda activate srst2-dev
#conda install -c conda-forge biopython

# Cluster similar sequences using CD-HIT-EST at 90% identity
cd-hit-est -i AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta -o AMR_CDS_plus.cdhit90 -d 0 > AMR_CDS_plus.cdhit90.stdout

# Show CD-HIT output
cat AMR_CDS_plus.cdhit90.stdout

# Convert CD-HIT cluster file to CSV using a custom Python script
python ../database_clustering/cdhit_to_csv.py --cluster_file AMR_CDS_plus.cdhit90.clstr --infasta AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta --outfile AMR_CDS_plus.clustered.csv

# Generate clustered gene database FASTA from CSV using a custom Python script
python ../database_clustering/csv_to_gene_db.py -t AMR_CDS_plus.clustered.csv -o AMR_CDS_plus.clustered.fasta -f AMR_CDS_plus_tetLMS.def.nogaps.nosym.fasta -c 4
# The output FASTA (seqs_clustered.fasta) is ready for SRST2 (--gene_db).

# Build Bowtie2 index for the clustered gene database
bowtie2-build AMR_CDS_plus.clustered.fasta AMR_CDS_plus.clustered.fasta
