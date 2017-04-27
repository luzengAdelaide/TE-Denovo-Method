#!/bin/bash

# This script used for calculate the repeat content in target genome from CENSOR output,
# including specific repeat types

$1 = "genome name" (tuatara, human, bearded dragon, anolis, echidna, opossum, platypus, chicken)

# Go to directory contains censor output 
cd /data/rc003/lu/$i/censor_combine_new

# Make a directory to store all the censor output files
mkdir maps
mv *.map maps
mkdir logs
mv *.log logs
mkdir censor_process
mv *.fa.* censor_process
mkdir output_analysis
mkdir slurm_output
mv slurm* slurm_output

# Merge the map files together, and move to a empty to do further analysis
cd maps
cat *.map > '$1'_combine.map
mv '$1'_combine.map '$1'/output_analysis

# CENSOR will ignore header character after space, so when running CENSOR I replaced the space to '*' 

############################## Tuatara ###########################
sed 's/*/\t/g' tuatara_combine.map > tuatara_combine.map.use
###################################################################
# Extract repeat type from censor output
awk '{if($5=="LTR") print $0}' tuatara_combine.map.use > LTR.repeat.map 
awk '{if($5=="DNA") print $0}' tuatara_combine.map.use > DNA.repeat.map 
awk '{if($5=="CR1") print $0}' tuatara_combine.map.use> CR1.repeat.map 
awk '{if($5=="Copia") print $0}' tuatara_combine.map.use > Copia.repeat.map 
awk '{if($5=="hAT") print $0}' tuatara_combine.map.use > hAT.repeat.map 
awk '{if($5=="Harbinger") print $0}' tuatara_combine.map.use > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' tuatara_combine.map.use > Helitron.repeat.map 
grep 'Mariner' tuatara_combine.map.use > Mariner.repeat.map 
awk '{if($5=="SINE" || $5=="SINE2/tRNA" || $5=="SINE1/7SL") print $0}' tuatara_combine.map.use > SINE.repeat.map 
awk '{if($5=="Gypsy") print $0}' tuatara_combine.map.use > Gypsy.repeat.map 
awk '{if($5=="Polinton") print $0}' tuatara_combine.map.use > Polinton.repeat.map 
awk '{if($5=="RTEX") print $0}' tuatara_combine.map.use > RTEX.repeat.map 
awk '{if($5=="L2") print $0}' tuatara_combine.map.use > L2.repeat.map 
awk '{if($5=="L1") print $0}' tuatara_combine.map.use > L1.repeat.map 
awk '{if($5=="DIRS") print $0}' tuatara_combine.map.use > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' tuatara_combine.map.use > SAT.repeat.map
awk '{if($5=="BEL") print $0}' tuatara_combine.map.use > BEL.repeat.map
awk '{if($5=="Penelope") print $0}' tuatara_combine.map.use > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' tuatara_combine.map.use > Jockey.repeat.map
grep ')n' tuatara_combine.map.use > simple.repeat.map
awk '{if($5=="Daphne") print $0}' tuatara_combine.map.use > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' tuatara_combine.map.use > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' tuatara_combine.map.use > RTE.repeat.map
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' tuatara_combine.map.use   > ERV.repeat.map
awk '{if($5=="SAT") print $0}' tuatara_combine.map.use > SAT.repeat.map
awk '{if($5=="Transposable") print $0}' tuatara_combine.map.use > transposable.repeat.map
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt # Repeat type name
    wc -l $i > count.txt # Copy number for each repeat type
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt # total length for each repeat type
    awk '{len+=$3-$2+1}{print (len/4271823035)*100}' $i|tail -1 > test3.txt #total proportion for each repeat type in genome
    paste test.txt count.txt test2.txt test3.txt >> tuatara_repeat_constitution.txt 
done



###################### chicken #####################
sed 's/*/\t/g' chicken_combine.map > chicken_combine.map.use
#####################################################
awk '{if($5=="CR1") print $0}' chicken_combine.map.use  > CR1.repeat.map 
awk '{if($5=="Gypsy") print $0}' chicken_combine.map.use  > Gypsy.repeat.map 
awk '{if($5=="hAT") print $0}' chicken_combine.map.use  > hAT.repeat.map 
awk '{if($5=="L1") print $0}' chicken_combine.map.use  > L1.repeat.map 
awk '{if($5=="DNA") print $0}' chicken_combine.map.use  > DNA.repeat.map 
grep ')n' chicken_combine.map.use  > simple.repeat.map
grep 'Mariner' chicken_combine.map.use  > Mariner.repeat.map 
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' chicken_combine.map.use  > ERV.repeat.map
awk '{if($5=="Copia") print $0}' chicken_combine.map.use  > Copia.repeat.map 
awk '{if($5=="EnSpm/CACTA") print $0}' chicken_combine.map.use  > EnSpm.repeat.map 
awk '{if($5=="BEL") print $0}' chicken_combine.map.use  > BEL.repeat.map
awk '{if($5=="Harbinger") print $0}' chicken_combine.map.use  > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' chicken_combine.map.use  > Helitron.repeat.map 
awk '{if($5=="Tx1") print $0}' chicken_combine.map.use  > Tx1.repeat.map 
awk '{if($5=="Polinton") print $0}' chicken_combine.map.use  > Polinton.repeat.map 
awk '{if($5=="LTR") print $0}' chicken_combine.map.use  > LTR.repeat.map 
awk '{if($5=="MuDR") print $0}' chicken_combine.map.use  > MuDR.repeat.map 
awk '{if($5=="L2") print $0}' chicken_combine.map.use  > L2.repeat.map 
awk '{if($5=="DIRS") print $0}' chicken_combine.map.use  > DIRS.repeat.map
awk '{if($5=="Transposable") print $0}' chicken_combine.map.use  > Transposable.repeat.map
grep 'Charlie' chicken_combine.map.use > Charlie.repeat.map
awk '{if($5=="SINE" || $5=="SINE2/tRNA") print $0}' chicken_combine.map.use  > SINE.repeat.map 
awk '{if($5=="RTEX") print $0}' chicken_combine.map.use  > RTEX.repeat.map 
awk '{if($5=="DIRS") print $0}' chicken_combine.map.use  > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' chicken_combine.map.use  > SAT.repeat.map
awk '{if($5=="Penelope") print $0}' chicken_combine.map.use  > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' chicken_combine.map.use  > Jockey.repeat.map
awk '{if($5=="Daphne") print $0}' chicken_combine.map.use  > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' chicken_combine.map.use  > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' chicken_combine.map.use  > RTE.repeat.map

awk '{if($5=="Sola") print $0}' chicken_combine.map.use  > Sola.repeat.map
awk '{if($5=="Kolobok") print $0}' chicken_combine.map.use  > Kolobok.repeat.map
awk '{if($5=="Daphne") print $0}' chicken_combine.map.use  > Daphne.repeat.map
awk '{if($5=="piggyBac") print $0}' chicken_combine.map.use  > PiggyBac.repeat.map
awk '{if($5=="I") print $0}' chicken_combine.map.use  > I.repeat.map
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt  # Repeat type name
    wc -l $i > count.txt
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt
    awk '{len+=$3-$2+1}{print (len/1046932099)*100}' $i|tail -1 > test3.txt
    paste count.txt test.txt test2.txt test3.txt >> chicken_repeat_constitution2.txt
done


###################### Echidna #####################
sed 's/*/\t/g' tachy_combine.map > tachy_combine.map.use
############################################################
awk '{if($5=="CR1") print $0}' tachy_combine.map.use  > CR1.repeat.map 
awk '{if($5=="Gypsy") print $0}' tachy_combine.map.use  > Gypsy.repeat.map 
awk '{if($5=="hAT") print $0}' tachy_combine.map.use  > hAT.repeat.map 
awk '{if($5=="L1") print $0}' tachy_combine.map.use  > L1.repeat.map 
awk '{if($5=="DNA") print $0}' tachy_combine.map.use  > DNA.repeat.map 
grep ')n' tachy_combine.map.use  > simple.repeat.map
grep 'Mariner' tachy_combine.map.use  > Mariner.repeat.map 
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' tachy_combine.map.use  > ERV.repeat.map
awk '{if($5=="Copia") print $0}' tachy_combine.map.use  > Copia.repeat.map 
awk '{if($5=="EnSpm/CACTA") print $0}' tachy_combine.map.use  > EnSpm.repeat.map 
awk '{if($5=="BEL") print $0}' tachy_combine.map.use  > BEL.repeat.map
awk '{if($5=="Harbinger") print $0}' tachy_combine.map.use  > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' tachy_combine.map.use  > Helitron.repeat.map 
awk '{if($5=="Tx1") print $0}' tachy_combine.map.use  > Tx1.repeat.map 
awk '{if($5=="Polinton") print $0}' tachy_combine.map.use  > Polinton.repeat.map 
awk '{if($5=="LTR") print $0}' tachy_combine.map.use  > LTR.repeat.map 
awk '{if($5=="MuDR") print $0}' tachy_combine.map.use  > MuDR.repeat.map 
awk '{if($5=="L2") print $0}' tachy_combine.map.use  > L2.repeat.map 
awk '{if($5=="DIRS") print $0}' tachy_combine.map.use  > DIRS.repeat.map
awk '{if($5=="Transposable") print $0}' tachy_combine.map.use  > Transposable.repeat.map
grep 'Charlie' tachy_combine.map.use > Charlie.repeat.map
awk '{if($5=="SINE" || $5=="SINE2/tRNA") print $0}' tachy_combine.map.use  > SINE.repeat.map 
awk '{if($5=="RTEX") print $0}' tachy_combine.map.use  > RTEX.repeat.map 
awk '{if($5=="DIRS") print $0}' tachy_combine.map.use  > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' tachy_combine.map.use  > SAT.repeat.map
awk '{if($5=="Penelope") print $0}' tachy_combine.map.use  > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' tachy_combine.map.use  > Jockey.repeat.map
awk '{if($5=="Daphne") print $0}' tachy_combine.map.use  > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' tachy_combine.map.use  > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' tachy_combine.map.use  > RTE.repeat.map

awk '{if($5=="Sola") print $0}' tachy_combine.map.use  > Sola.repeat.map
awk '{if($5=="Kolobok") print $0}' tachy_combine.map.use  > Kolobok.repeat.map
awk '{if($5=="Daphne") print $0}' tachy_combine.map.use  > Daphne.repeat.map
awk '{if($5=="piggyBac") print $0}' tachy_combine.map.use  > PiggyBac.repeat.map
awk '{if($5=="I") print $0}' tachy_combine.map.use  > I.repeat.map
############################################################
grep 'L2' CR1.repeat.map > CR1-L2.repeat.map
sed '/L2/d' CR1.repeat.map > CR1-noL2.repeat.map
grep 'MON1' SINE.repeat.map > SINE-MON1.repeat.map
sed '/MON1/d' SINE.repeat.map > SINE-noMON1.repeat.map
sed '/hAT/d' Charlie.repeat.map > test
mv test Charlie.repeat.map 
sed '/Mariner/d' DNA.repeat.map > test
mv test DNA.repeat.map 
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt
    wc -l $i > count.txt
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt
    awk '{len+=$3-$2+1}{print (len/2020007912)*100}' $i|tail -1 > test3.txt
    paste test.txt count.txt test2.txt test3.txt >> echidna_repeat_constitution2.txt
done



###################### Platypus #####################
sed 's/*/\t/g' oana_combine.map > oana_combine.map.use
############################################################
awk '{if($5=="CR1") print $0}' oana_combine.map.use  > CR1.repeat.map 
awk '{if($5=="Gypsy") print $0}' oana_combine.map.use  > Gypsy.repeat.map 
awk '{if($5=="hAT") print $0}' oana_combine.map.use  > hAT.repeat.map 
awk '{if($5=="L1") print $0}' oana_combine.map.use  > L1.repeat.map 
awk '{if($5=="DNA") print $0}' oana_combine.map.use  > DNA.repeat.map 
grep ')n' oana_combine.map.use  > simple.repeat.map
grep 'Mariner' oana_combine.map.use  > Mariner.repeat.map 
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' oana_combine.map.use  > ERV.repeat.map
awk '{if($5=="Copia") print $0}' oana_combine.map.use  > Copia.repeat.map 
awk '{if($5=="EnSpm/CACTA") print $0}' oana_combine.map.use  > EnSpm.repeat.map 
awk '{if($5=="BEL") print $0}' oana_combine.map.use  > BEL.repeat.map
awk '{if($5=="Harbinger") print $0}' oana_combine.map.use  > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' oana_combine.map.use  > Helitron.repeat.map 
awk '{if($5=="Tx1") print $0}' oana_combine.map.use  > Tx1.repeat.map 
awk '{if($5=="Polinton") print $0}' oana_combine.map.use  > Polinton.repeat.map 
awk '{if($5=="LTR") print $0}' oana_combine.map.use  > LTR.repeat.map 
awk '{if($5=="MuDR") print $0}' oana_combine.map.use  > MuDR.repeat.map 
awk '{if($5=="L2") print $0}' oana_combine.map.use  > L2.repeat.map 
awk '{if($5=="DIRS") print $0}' oana_combine.map.use  > DIRS.repeat.map
awk '{if($5=="Transposable") print $0}' oana_combine.map.use  > Transposable.repeat.map
grep 'Charlie' oana_combine.map.use > Charlie.repeat.map
awk '{if($5=="SINE" || $5=="SINE2/tRNA") print $0}' oana_combine.map.use  > SINE.repeat.map 
awk '{if($5=="RTEX") print $0}' oana_combine.map.use  > RTEX.repeat.map 
awk '{if($5=="DIRS") print $0}' oana_combine.map.use  > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' oana_combine.map.use  > SAT.repeat.map
awk '{if($5=="Penelope") print $0}' oana_combine.map.use  > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' oana_combine.map.use  > Jockey.repeat.map
awk '{if($5=="Daphne") print $0}' oana_combine.map.use  > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' oana_combine.map.use  > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' oana_combine.map.use  > RTE.repeat.map

awk '{if($5=="Sola") print $0}' oana_combine.map.use  > Sola.repeat.map
awk '{if($5=="Kolobok") print $0}' oana_combine.map.use  > Kolobok.repeat.map
awk '{if($5=="Daphne") print $0}' oana_combine.map.use  > Daphne.repeat.map
awk '{if($5=="piggyBac") print $0}' oana_combine.map.use  > PiggyBac.repeat.map
awk '{if($5=="I") print $0}' oana_combine.map.use  > I.repeat.map
############################################################
grep 'L2' CR1.repeat.map > CR1-L2.repeat.map
sed '/L2/d' CR1.repeat.map > CR1-noL2.repeat.map
grep 'MON1' SINE.repeat.map > SINE-MON1.repeat.map
sed '/MON1/d' SINE.repeat.map > SINE-noMON1.repeat.map
sed '/hAT/d' Charlie.repeat.map > test
mv test Charlie.repeat.map 
sed '/Mariner/d' DNA.repeat.map > test
mv test DNA.repeat.map 
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt
    wc -l $i > count.txt
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt
    awk '{len+=$3-$2+1}{print (len/2073148626)*100}' $i|tail -1 > test3.txt
    paste test.txt count.txt test2.txt test3.txt >> platypus_repeat_constitution2.txt
done




###################### human #####################
sed 's/*/\t/g' human_combine.map > human_combine.map.use
############################################################
awk '{if($5=="CR1") print $0}' human_combine.map.use  > CR1.repeat.map 
awk '{if($5=="Gypsy") print $0}' human_combine.map.use  > Gypsy.repeat.map 
awk '{if($5=="hAT") print $0}' human_combine.map.use  > hAT.repeat.map 
awk '{if($5=="L1") print $0}' human_combine.map.use  > L1.repeat.map 
awk '{if($5=="DNA") print $0}' human_combine.map.use  > DNA.repeat.map 
grep ')n' human_combine.map.use  > simple.repeat.map
grep 'Mariner' human_combine.map.use  > Mariner.repeat.map 
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' human_combine.map.use  > ERV.repeat.map
awk '{if($5=="Copia") print $0}' human_combine.map.use  > Copia.repeat.map 
awk '{if($5=="EnSpm/CACTA") print $0}' human_combine.map.use  > EnSpm.repeat.map 
awk '{if($5=="BEL") print $0}' human_combine.map.use  > BEL.repeat.map
awk '{if($5=="Harbinger") print $0}' human_combine.map.use  > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' human_combine.map.use  > Helitron.repeat.map 
awk '{if($5=="Tx1") print $0}' human_combine.map.use  > Tx1.repeat.map 
awk '{if($5=="Polinton") print $0}' human_combine.map.use  > Polinton.repeat.map 
awk '{if($5=="LTR") print $0}' human_combine.map.use  > LTR.repeat.map 
awk '{if($5=="MuDR") print $0}' human_combine.map.use  > MuDR.repeat.map 
awk '{if($5=="L2") print $0}' human_combine.map.use  > L2.repeat.map 
awk '{if($5=="DIRS") print $0}' human_combine.map.use  > DIRS.repeat.map
awk '{if($5=="Transposable") print $0}' human_combine.map.use  > Transposable.repeat.map
grep 'Charlie' human_combine.map.use > Charlie.repeat.map
awk '{if($5=="SINE" || $5=="SINE2/tRNA" || $5=="SINE1/7SL") print $0}' human_combine.map.use  > SINE.repeat.map 
awk '{if($5=="RTEX") print $0}' human_combine.map.use  > RTEX.repeat.map 
awk '{if($5=="DIRS") print $0}' human_combine.map.use  > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' human_combine.map.use  > SAT.repeat.map
awk '{if($5=="Penelope") print $0}' human_combine.map.use  > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' human_combine.map.use  > Jockey.repeat.map
awk '{if($5=="Daphne") print $0}' human_combine.map.use  > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' human_combine.map.use  > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' human_combine.map.use  > RTE.repeat.map

awk '{if($5=="Sola") print $0}' human_combine.map.use  > Sola.repeat.map
awk '{if($5=="Kolobok") print $0}' human_combine.map.use  > Kolobok.repeat.map
awk '{if($5=="Daphne") print $0}' human_combine.map.use  > Daphne.repeat.map
awk '{if($5=="piggyBac") print $0}' human_combine.map.use  > PiggyBac.repeat.map
awk '{if($5=="I") print $0}' human_combine.map.use  > I.repeat.map
awk '{if($5=="Endogenous") print $0}' human_combine.map.use  > Endogenous.repeat.map
############################################################
sed '/hAT/d' Charlie.repeat.map > test
mv test Charlie.repeat.map 
sed '/Mariner/d' DNA.repeat.map > test
mv test DNA.repeat.map 
grep 'ALU\|Alu' SINE.repeat.map > Alu.repeat.map
grep 'MIR' SINE.repeat.map > MIR.repeat.map
sed '/ALU\|Alu/d' SINE.repeat.map > test
mv test SINE.repeat.map 
sed '/MIR/d' SINE.repeat.map > test
mv test SINE.repeat.map
grep 'SVA' SINE.repeat.map > SVA.repeat.map
sed '/SVA/d' SINE.repeat.map > test
mv test SINE.repeat.map
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt
    wc -l $i > count.txt
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt
    awk '{len+=$3-$2+1}{print (len/3095677412)*100}' $i|tail -1 > test3.txt
    paste test.txt count.txt test2.txt test3.txt >> human_repeat_constitution2.txt
done



###################### Opossum #####################
sed 's/*/\t/g' mdo_combine.map > mdo_combine.map.use
############################################################
awk '{if($5=="CR1") print $0}' mdo_combine.map.use  > CR1.repeat.map 
awk '{if($5=="Gypsy") print $0}' mdo_combine.map.use  > Gypsy.repeat.map 
awk '{if($5=="hAT") print $0}' mdo_combine.map.use  > hAT.repeat.map 
awk '{if($5=="L1") print $0}' mdo_combine.map.use  > L1.repeat.map 
awk '{if($5=="DNA") print $0}' mdo_combine.map.use  > DNA.repeat.map 
grep ')n' mdo_combine.map.use  > simple.repeat.map
grep 'Mariner' mdo_combine.map.use  > Mariner.repeat.map 
awk '{if($5=="ERV1" || $5=="ERV2" || $5=="ERV3") print $0}' mdo_combine.map.use  > ERV.repeat.map
awk '{if($5=="Copia") print $0}' mdo_combine.map.use  > Copia.repeat.map 
awk '{if($5=="EnSpm/CACTA") print $0}' mdo_combine.map.use  > EnSpm.repeat.map 
awk '{if($5=="BEL") print $0}' mdo_combine.map.use  > BEL.repeat.map
awk '{if($5=="Harbinger") print $0}' mdo_combine.map.use  > Harbinger.repeat.map 
awk '{if($5=="Helitron") print $0}' mdo_combine.map.use  > Helitron.repeat.map 
awk '{if($5=="Tx1") print $0}' mdo_combine.map.use  > Tx1.repeat.map 
awk '{if($5=="Polinton") print $0}' mdo_combine.map.use  > Polinton.repeat.map 
awk '{if($5=="LTR") print $0}' mdo_combine.map.use  > LTR.repeat.map 
awk '{if($5=="MuDR") print $0}' mdo_combine.map.use  > MuDR.repeat.map 
awk '{if($5=="L2") print $0}' mdo_combine.map.use  > L2.repeat.map 
awk '{if($5=="DIRS") print $0}' mdo_combine.map.use  > DIRS.repeat.map
awk '{if($5=="Transposable") print $0}' mdo_combine.map.use  > Transposable.repeat.map
grep 'Charlie' mdo_combine.map.use > Charlie.repeat.map
awk '{if($5=="SINE" || $5=="SINE2/tRNA" || $5=="SINE1/7SL") print $0}' mdo_combine.map.use  > SINE.repeat.map 
awk '{if($5=="RTEX") print $0}' mdo_combine.map.use  > RTEX.repeat.map 
awk '{if($5=="DIRS") print $0}' mdo_combine.map.use  > DIRS.repeat.map
awk '{if($5=="SAT2" || $5=="SAT") print $0}' mdo_combine.map.use  > SAT.repeat.map
awk '{if($5=="Penelope") print $0}' mdo_combine.map.use  > Penelope.repeat.map
awk '{if($5=="Jockey") print $0}' mdo_combine.map.use  > Jockey.repeat.map
awk '{if($5=="Daphne") print $0}' mdo_combine.map.use  > Daphne.repeat.map
awk '{if($5=="Non-LTR") print $0}' mdo_combine.map.use  > Non-LTR.repeat.map
awk '{if($5=="RTE") print $0}' mdo_combine.map.use  > RTE.repeat.map
awk '{if($5=="Sola") print $0}' mdo_combine.map.use  > Sola.repeat.map
awk '{if($5=="Kolobok") print $0}' mdo_combine.map.use  > Kolobok.repeat.map
awk '{if($5=="Daphne") print $0}' mdo_combine.map.use  > Daphne.repeat.map
awk '{if($5=="piggyBac") print $0}' mdo_combine.map.use  > PiggyBac.repeat.map
awk '{if($5=="I") print $0}' mdo_combine.map.use  > I.repeat.map
awk '{if($5=="Endogenous") print $0}' mdo_combine.map.use  > Endogenous.repeat.map
############################################################
sed '/Charlie/d' hAT.repeat.map > test
mv test hAT.repeat.map 
sed '/Mariner/d' DNA.repeat.map > test
mv test DNA.repeat.map 
grep 'THER' SINE.repeat.map >THER.repeat.map
sed '/THER/d' SINE.repeat.map > test
mv test SINE.repeat.map
grep 'MIR' SINE.repeat.map > MIR.repeat.map
sed '/MIR/d' SINE.repeat.map > test
mv test SINE.repeat.map
grep 'SINE-1' SINE.repeat.map > SINE1.repeat.map
sed '/SINE-1/d' SINE.repeat.map > test
mv test SINE.repeat.map
##########################################################
# Calculate repeat content in each genome with respect to specific repeat types
for i in *.map;
do
    echo $i > test.txt
    wc -l $i > count.txt
    awk '{len+=$3-$2+1}{print len}' $i|tail -1 > test2.txt
    awk '{len+=$3-$2+1}{print (len/3605631728)*100}' $i|tail -1 > test3.txt
    paste test.txt count.txt test2.txt test3.txt >> opossum_repeat_constitution2.txt
done


