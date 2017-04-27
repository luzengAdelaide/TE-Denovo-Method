#!/bin/bash

# Extract ortholog genes between species, two by two

# Extract exon information from each genome
# make index for each fa
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load BEDTools/2.25.0-GCC-5.3.0-binutils-2.25
module load BLAST/2.2.26-Linux_x86_64 

# Index fasta sequences
for i in *.fa
do
samtools faidx $i
done

# Extract transcript information from gtf files
awk '{if($3 == "exon") print }' Gallus_gallus.Galgal4.83.gtf > gallus_transcript.gff
awk '{if($3 == "exon") print }' Monodelphis_domestica.BROADO5.83.gtf > opossum_transcript.gff
awk '{if($3 == "exon") print }' Ornithorhynchus_anatinus.OANA5.83.gtf > platypus_transcript.gff
awk '{if($3 == "exon") print }' Mus_musculus.GRCm38.83.gtf > mouse_transcript.gff
awk '{if($3 == "exon") print }' hg19v37.refgene.gtf > human_transcript.gff

# Change format catered to gid
for i in *_transcript.gff
do
	perl format_change.pl $i > $i.use
done

# Extract sequences
bedtools getfasta -s -fi Gallus_gallus.Galgal4.dna.toplevel.fa -bed gallus_transcript.gff.use -name -fo gallus_exon.fasta
bedtools getfasta -s -fi Monodelphis_domestica.BROADO5.dna.toplevel.fa -bed opossum_transcript.gff.use -name -fo opossum_exon.fasta
bedtools getfasta -s -fi Mus_musculus.GRCm38.dna.toplevel.fa -bed mouse_transcript.gff.use -name -fo mouse_exon.fasta
bedtools getfasta -s -fi Ornithorhynchus_anatinus.OANA5.dna.toplevel.fa -bed platypus_transcript.gff.use -name -fo platypus_exon.fasta
bedtools getfasta -s -fi hg19v37.fa -bed human_transcript.gff.use -name -fo human_exon.fasta


# Use blastn to run them
# Human and Chicken
module load blast+/2.2.30 
blastn -query human_exon.fasta.use -subject gallus_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_hg_gal.out
blastn -query gallus_exon.fasta.use -subject gallus_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_gal_hg.out

# Human and bearded dragon
blastn -query human_exon.fasta.use -subject Pogona_vitticeps.male.cds -outfmt 6 -word_size 9 -out blastn_hg_bg.out
blastn -query Pogona_vitticeps.male.cds -subject human_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_bg_hg.out

# Human and mouse
blastn -query human_exon.fasta.use -subject mouse_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_hg_mm.out
blastn -query mouse_exon.fasta.use -subject human_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_mm_hg.out

# Human and platypus
blastn -query human_exon.fasta.use -subject platypus_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_hg_oana.out
blastn -query platypus_exon.fasta.use -subject human_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_oana_hg.out

# Human and Opossum
blastn -query human_exon.fasta.use -subject opossum_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_hg_mdo.out
blastn -query opossum_exon.fasta.use -subject human_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_mdo_hg.out

# Human and Echidna
blastn -query human_exon.fasta.use -subject tachy_cds.fa.use -outfmt 6 -word_size 9 -out blastn_hg_tachy.out
blastn -query tachy_cds.fa.use -subject human_exon.fasta.use -outfmt 6 -word_size 9 -out blastn_tachy_hg.out

################################################

# Way to get orthologs between paired blastn output
cat blastn_hg_gal.out blastn_gal_hg.out | ./recip > uniq_hg_gal
cat blastn_hg_mdo.out blastn_mdo_hg.out | ./recip > uniq_hg_mdo
cat blastn_hg_oana.out blastn_oana_hg.out | ./recip > uniq_hg_oana
cat blastn_hg_tachy.out blastn_tachy_hg.out | ./recip > uniq_hg_tachy
cat blastn_hg_bg.out blastn_bg_hg.out | ./recip > uniq_hg_bg

# Then make the human id in the first column
for i in uniq_hg_*
do
	echo $i
	perl change_position.sh $i
done


# Extract ortholog gene names
perl merge_same_id.pl uniq_hg_bg.use uniq_hg_gal.use > ortho_hg_bg_gal
perl merge_same_id.pl ortho_hg_bg_gal uniq_hg_mdo.use > ortho_hg_bg_gal_mdo
perl merge_same_id.pl ortho_hg_bg_gal_mdo uniq_hg_oana.use > ortho_hg_bg_gal_mdo_oana
perl merge_same_id.pl  ortho_hg_bg_gal_mdo_oana uniq_hg_tachy.use > ortho_hg_bg_gal_mdo_oana_tachy


# Extract name for each species
sed 's/\s/\t/g' ortho_hg_bg_gal_mdo_oana_tachy > ortho_hg_bg_gal_mdo_oana_tachy.use
cut -f 1 ortho_hg_bg_gal_mdo_oana_tachy.use > ortho_human
cut -f 2 ortho_hg_bg_gal_mdo_oana_tachy.use > ortho_tachy
cut -f 6 ortho_hg_bg_gal_mdo_oana_tachy.use > ortho_oana
cut -f 10 ortho_hg_bg_gal_mdo_oana_tachy.use > ortho_mdo
cut -f 14 ortho_hg_bg_gal_mdo_oana_tachy.use > ortho_gal
cut -f 18 ortho_hg_bg_gal_mdo_oana_tachy.use> ortho_bg


# Extract gene expression from ortholog genes 
# For other species, change name first, as I change name to difference species when running blast
sed 's/gal_//g' ortho_gal > ortho_gal.use
sed 's/human_//g' ortho_human > ortho_hg.use
sed 's/mdo_//g' ortho_mdo > ortho_mdo.use
sed 's/oana_//g' ortho_oana > ortho_oana.use

# Move into directory
mv ortho_human hg_RSEM/
mv ortho_bg bg_RSEM/
mv ortho_gal gal_RSEM/
mv ortho_oana oana_RSEM/
mv ortho_mdo mdo_RSEM/

# Pogona does not need to change name
 ./extract_gene_expression.sh


 