#!/bin/bash

# This pipeline is from Atma Ivancevic 'L1-dynamics'
# Look at RT in all L2 sequences
# From both repbase library and de novo library

# Go the directory with L2 sequences
cd /data/rc003/lu/consenuses_library/RT_analysis

# Translate sequences into amino acids
for i in *.fa
do
    usearch -fastx_findorfs $i -ntout $i.out -aaout $i.translated -mincodons 700
done

###### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database
for i in *.translated
do
    echo $i
    hmmscan --domtblout $i.out /mnt/lu/pfam/Pfam-A.hmm $i > $i.log
    grep '>' "$i" | sed 's/>//g' > "$i"_ORF2_queries.txt
    sed -i '/^#/ d' "$i".out
    awk '{print $1 "\t" $4}' "$i".out > "$i"_ORF2_filtered.tmp
    awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' "$i"_ORF2_filtered.tmp > "$i"_ORF2_sep.tmp
    cat "$i"_ORF2_sep.tmp | cut -f 1,2| awk '{if ( $1=="RVT_1" || $1=="RVT_3" ) print $2 }' | sort -u > "$i"_confirmedORF2.txt
    cat 2kb_bg.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "2kb_bg.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 2kb_bg.L2.fa.translated_confirmedORF2.fasta
    cat 2kb_mdo.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "2kb_mdo.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 2kb_mdo.L2.fa.translated_confirmedORF2.fasta
    cat 2kb_oana.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "2kb_oana.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 2kb_oana.L2.fa.translated_confirmedORF2.fasta
    cat 2kb_tua.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "2kb_tua.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > 2kb_tua.L2.fa.translated_confirmedORF2.fasta
    cat repbase_seq.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "repbase_seq.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > repbase_seq.L2.fa.translated_confirmedORF2.fasta
    cat RMK.L2.fa.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "RMK.L2.fa.translated_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > RMK.L2.fa.translated_confirmedORF2.fasta
    comm -2 -3 <(sort "$i"_ORF2_queries.txt) <(sort "$i"_confirmedORF2.txt) > "$i"_unconfirmedORF2.txt
    awk '$2 != prev {if (NR != 1) print prev; prev=$2; delete a};!($1 in a){a[$1]++; printf "%s,", $1}; END {print prev}' "$i"_ORF2_filtered.tmp > "$i"_ORF2_merged.tmp
    cat "$i"_ORF2_merged.tmp \
        | awk -F "," '{$NF = ""; print $0}' \
        | sed 's/ /,/g' \
        | sed 's/\,$//g' > 1sthalf.tmp
    cat "$i"_ORF2_merged.tmp \
        | awk -F "," '{print $NF}' > 2ndhalf.tmp
    paste 2ndhalf.tmp 1sthalf.tmp > "$i"_ORF2_domains.txt
    wc -l "$i"_ORF2_queries.txt \
        | awk '{print "# ORF2 candidates: " $1}' \
        > "$i"_ORF2_summary.txt
    wc -l "$i"_confirmedORF2.txt \
        | awk '{print "\n# candidates that contain RVT: " $1}' \
        >> "$i"_ORF2_summary.txt
    grep -c '>' "$i"_confirmedORF2.fasta \
        | awk '{print "\n# extracted FASTA seqs (should be same): " $1}' \
        >> "$i"_ORF2_summary.txt
    awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$i"_confirmedORF2.txt "$i"_ORF2_domains.txt \
        | awk '{print $2}' \
        | tr "," "\n" \
        | grep -v "^\s*$" \
        | sort \
        | uniq -c \
        | sort -bnr \
        | awk '{print $1 " " $2}' \
        | sed '1s/^/\nORF2 domains:\n/' \
        >> "$i"_ORF2_summary.txt
    wc -l "$i"_unconfirmedORF2.txt \
        | awk '{print "\n# candidates that do not contain RVT: " $1}' \
        >> "$i"_ORF2_summary.txt
    awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$i"_unconfirmedORF2.txt "$i"_ORF2_domains.txt \
        | awk '{print $2}' \
        | tr "," "\n" \
        | grep -v "^\s*$" \
        | sort \
        | uniq -c \
        | sort -bnr \
        | awk '{print $1 " " $2}' \
        | sed '1s/^/non-ORF2 domains:\n/' \
        >> "$i"_ORF2_summary.txt

done


##### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database
hmmscan --domtblout "$1"_confirmedORF2.fasta.out /data01/Protein_dbs/Pfam-A.hmm \
"$1"_confirmedORF2.fasta \
> "$1"_confirmedORF2.fasta.log

# pull out a BED file containing the location of the RT domains
# check that the RT domains validates as RVT_1 or RVT_3
# and that it is full-size (e.g. expect it to be around 255 amino acids,
# so set a restriction of >200 amino acids
# note: $20 and $21 are the envelope coordinates on the query protein seq
# so even if the alignment doesn't stretch as far, or scores low in some regions
# this 'envelope' spans the entire suspected RT domain
cat "$1"_confirmedORF2.fasta.out \
| awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' \
> "$1"_RT_domain.bed

# retrieve the fasta using this BED file
# check that you are not off by 1!! (e.g. if HMMer is 1-based and BEDTools is 0-based)
# no need to use strand since they are all facing the same way
fastaFromBed -fi "$1"_confirmedORF2.fasta -bed "$1"_RT_domain.bed -fo "$1"_RT_domain.fasta

# move to a easy-to-find location
mv "$1"_RT_domain.fasta /data01/Protein_dbs/RT_domains/
cd RT_domains

#!/bin/bash

for i in *_confirmedORF2.fasta
do
    hmmscan --domtblout "$i".out /mnt/lu/Pfam/Pfam-A.hmm \
       "$i" \
       > "$i".log
    cat "$i".out \
       | awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' \
       > "$i"_RT_domain.bed
    fastaFromBed -fi "$i" -bed "$i"_RT_domain.bed -fo "$i"_RT_domain.fasta
    mv "$i"_RT_domain.fasta RT_domains/
done

# Change the name of each fasta sequences
sed 's/>/&bg_/' 2kb_bg.L2.fa.translated_confirmedORF2.fasta_RT_domain.fasta > bg.L2_RT_domain.fasta
sed 's/>/&oana_/' 2kb_oana.L2.fa.translated_confirmedORF2.fasta_RT_domain.fasta > oana.L2_RT_domain.fasta
sed 's/>/&tua_/' 2kb_tua.L2.fa.translated_confirmedORF2.fasta_RT_domain.fasta > tua.L2_RT_domain.fasta
sed 's/>/&RMK_/' RMK.L2.fa.translated_confirmedORF2.fasta_RT_domain.fasta > RMK.L2_RT_domain.fasta
sed 's/>/&repbase_/' repbase_seq.L2.fa.translated_confirmedORF2.fasta_RT_domain.fasta > repbase.L2_RT_domain.fasta

# Make a directory to put original fasta sequences in
mkdir original_data
mv *confirmedORF2.fasta_RT_domain.fasta original_data

# Combine these file together
cat *_RT_domain.fasta > combined_L2_RT_domain.fasta

# Run muscle against fasta sequences
muscle -in combined_L2_RT_domain.fasta -out combined_L2_RT_domain.afa -maxiters 2




# Deal with a test tree that contain 4-5kb CR1, 4-5kb Repbase L2, 2-3kb denovo library L2
# Use README_TREE.sh in /data/rc003/lu/consenuses_library/L2_analysis/tree_building (phoenix)
# This is analysis after that
##################VM##################

for i in *.translated
do
    echo $i
    hmmscan --domtblout $i.out /mnt/lu/Pfam/Pfam-A.hmm $i > $i.log
done

grep '>' combined_unique.ORF2.candidate.translated  | sed 's/>//g' > combined_unique.ORF2.queries.txt
sed -i '/^#/ d' combined_unique.ORF2.candidate.translated.out 
awk '{print $1 "\t" $4}' combined_unique.ORF2.candidate.translated.out > combined_unique_ORF2_filtered.tmp
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' combined_unique_ORF2_filtered.tmp > combined_unique_ORF2_sep.tmp
cat *_ORF2_sep.tmp | cut -f 1,2| awk '{if ( $1=="RVT_1" || $1=="RVT_3" ) print $2 }' | sort -u > combined_unique_confirmedORF2.txt
cat combined_unique.ORF2.candidate.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "combined_unique_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > combined_unique_confirmedORF2.fasta
comm -2 -3 <(sort *ORF2.queries.txt) <(sort *_confirmedORF2.txt) > combined_unique_unconfirmedORF2.txt
awk '$2 != prev {if (NR != 1) print prev; prev=$2; delete a}; !($1 in a){a[$1]++; printf "%s,", $1}; END {print prev}' *_ORF2_filtered.tmp > combined_unique_ORF2_merged.tmp
cat *_ORF2_merged.tmp | awk -F "," '{$NF = ""; print $0}' | sed 's/ /,/g' | sed 's/\,$//g' > 1sthalf.tmp
cat *_ORF2_merged.tmp | awk -F "," '{print $NF}' > 2ndhalf.tmp
paste 2ndhalf.tmp 1sthalf.tmp > combined_unique_ORF2_domains.txt

###### Write out a summary log file ######
# print the original number of ORF2 candidate seqs
wc -l *_ORF2_queries.txt \
| awk '{print "# ORF2 candidates: " $1}' \
> *_ORF2_summary.txt

# print how many confirmed ORF2 there are
wc -l *_confirmedORF2.txt \
| awk '{print "\n# candidates that contain RVT: " $1}' \
>> *_ORF2_summary.txt

# check that the FASTA file of confirmed ORF2 is the same
grep -c '>' *_confirmedORF2.fasta \
| awk '{print "\n# extracted FASTA seqs (should be same): " $1}' \
>> *_ORF2_summary.txt

# print the domains in the confirmed ORF2 seqs
awk 'NR==FNR{c[$1]++;next};c[$1]>0' *_confirmedORF2.txt *_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/\nORF2 domains:\n/' \
>> *_ORF2_summary.txt

# print how many ORF2 do not contain RVT
wc -l *_unconfirmedORF2.txt \
| awk '{print "\n# candidates that do not contain RVT: " $1}' \
>> *_ORF2_summary.txt

# print the domains in the unconfirmed ORF2 candidates
awk 'NR==FNR{c[$1]++;next};c[$1]>0' *_unconfirmedORF2.txt *_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/non-ORF2 domains:\n/' \
>> *_ORF2_summary.txt
mv \*_ORF2_summary.txt combined_unique_ORF2_summary.txt

# Confirmed RT domain
hmmscan --domtblout combined_unique_confirmedORF2.fasta.out /mnt/lu/Pfam/Pfam-A.hmm combined_unique_confirmedORF2.fasta > combined_unique_confirmedORF2.fasta.log

cat ombined_unique_confirmedORF2.fasta.out | awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' > ombined_unique_RT_domain.bed
fastaFromBed -fi combined_unique_confirmedORF2.fasta -bed combined_unique_RT_domain.bed -fo combined_unique_RT_domain.fasta

muscle -in combined_unique_RT_domain.fasta -out combined_unique_RT_domain.afa -maxiters 2

###################LEEUWENHOEK###################
/scratch/luAnalysis/consensus_library/L2_real_tree
fasttree combined_unique_RT_domain.afa > RT_domain.tree



############## L2 RT Domain, combine with our L2 and Repbase L2 (Anolis, crocodile and turtle)
cd /data/rc003/lu/consenuses_library/L2_analysis/RT_analysis
perl extract_specific_fasta.pl L2_Repbase.fa.use > test
mv test L2_Repbase.fa.use

usearch -fastx_findorfs L2_RT_analysis.fa -ntout L2_RT_analysis.out -aaout L2_RT_analysis.translated -orfstyle 7 -mincodons 500

# Transfer to VM
cd /mnt/lu/consensus_library/RT_analysis






