#!/bin/bash

# This pipeline is from Atma Ivancevic 'L1-dynamics'
# Use emboss getorf to find and extract open reading frames
# HMMer confirmation for amino acid ORF2 candidates

# Load emboss
module load EMBOSS/6.5.7-GCC-5.3.0-binutils-2.25

# $1 = species name

# Translate the candidate sequences into amino acid sequences
usearch -fastx_findorfs 2kb_"$1".L2.fa -ntout "$1".ORF2.candidate.fasta -aaout "$1".ORF2.candidate.translated -mincodons 800
# -fastx_findorfs means input is the fasta sequences
# use default ofrstyle(5), 7 is the most loose number

###### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database
hmmscan --domtblout "$1".L2.candidate.translated.out pfam/Pfam-A.hmm "$1".L2.candidate.translated > "$1".L2.candidate.translated.log

###### Make a list of all the query seq names (for reference) ######
# list all the headers from the orf candidates that were queried against Pfam
grep '>' "$1".ORF2.candidate.translated | sed 's/>//g' > "$1"_ORF2_queries.txt

###### Confirm ORF2 candidates by recognising RT domains ######
# remove comment lines
sed -i '/^#/ d' "$1".ORF2.candidate.translated.out

# extract the fields of interest
awk '{print $1 "\t" $4}' "$1".ORF2.candidate.translated.out > "$1"_ORF2_filtered.tmp

# separate each ORF candidate sequence by a blank line
# ie. group all the domains found in each candidate seq
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' "$1"_ORF2_filtered.tmp > "$1"_ORF2_sep.tmp

# grab the headers from output that contain "RVT_*" (ORF2 reverse transcriptase)
cat "$1"_ORF2_sep.tmp | cut -f 1,2| awk '{if ( $1=="RVT_1" || $1=="RVT_3" ) print $2 }' | sort -u > "$1"_confirmedORF2.txt

##### Extract FASTA sequences from confirmed ORF2 ######
# extract the amino acid seq of the confirmed ORFs from the original query FASTA file
cat "$1"_extended_ORF2_candidates_translation.fasta | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < #"'$1'_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > "$1"_confirmedORF2.fasta
#cat combined_unique.ORF2.candidate.translated | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -F"\t" 'BEGIN{while((getline k < "combined_unique_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > combined_unique_confirmedORF2.fasta

###### Keep track of the candidates without confirmed RT domains ######
# list headers that did not contain confirmed seqs
# (either had no recognisable domains, or none that contained RVT)
comm -2 -3 <(sort "$1"_ORF2_queries.txt) <(sort "$1"_confirmedORF2.txt) > "$1"_unconfirmedORF2.txt

###### Keep track of the domains found in each ORF2 candidate sequence ######
# (both confirmed and unconfirmed)
# merge domains according to the header name
# so that each line has a unique header
awk '$2 != prev {if (NR != 1) print prev; prev=$2; delete a};
!($1 in a){a[$1]++; printf "%s,", $1};
END {print prev}' "$1"_ORF2_filtered.tmp > "$1"_ORF2_merged.tmp

# rearrange to put into nicer format
# e.g. extract all domains and seperate them by commas
cat "$1"_ORF2_merged.tmp \
| awk -F "," '{$NF = ""; print $0}' \
| sed 's/ /,/g' \
| sed 's/\,$//g' > 1sthalf.tmp

# then extract the corresponding seq name
cat "$1"_ORF2_merged.tmp \
| awk -F "," '{print $NF}' > 2ndhalf.tmp

# paste together (tab-delimited) with col1=header, col2=domains
paste 2ndhalf.tmp 1sthalf.tmp > "$1"_ORF2_domains.txt

###### Write out a summary log file ######
# print the original number of ORF2 candidate seqs
wc -l "$1"_ORF2_queries.txt \
| awk '{print "# ORF2 candidates: " $1}' \
> "$1"_ORF2_summary.txt

# print how many confirmed ORF2 there are
wc -l "$1"_confirmedORF2.txt \
| awk '{print "\n# candidates that contain RVT: " $1}' \
>> "$1"_ORF2_summary.txt

# check that the FASTA file of confirmed ORF2 is the same
grep -c '>' "$1"_confirmedORF2.fasta \
| awk '{print "\n# extracted FASTA seqs (should be same): " $1}' \
>> "$1"_ORF2_summary.txt

# print the domains in the confirmed ORF2 seqs
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_confirmedORF2.txt "$1"_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/\nORF2 domains:\n/' \
>> "$1"_ORF2_summary.txt

# print how many ORF2 do not contain RVT
wc -l "$1"_unconfirmedORF2.txt \
| awk '{print "\n# candidates that do not contain RVT: " $1}' \
>> "$1"_ORF2_summary.txt

# print the domains in the unconfirmed ORF2 candidates
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_unconfirmedORF2.txt "$1"_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/non-ORF2 domains:\n/' \
>> "$1"_ORF2_summary.txt

###### Remove unnecessary files ######
# i.e. all files with extension .tmp
rm *.tmp

