# Build super consensus for platypus-like tuatara L2 and non-platypus like tuatara L2

muscle -in ht.nonplat.fa.use -out ht.nonplat.fa
piler -cons ht.nonplat.afa -out ht.nonplat.cons.fa -label ht.nonplat.gb.cons

muscle -in ht.plat.fa.use -out ht.plat.afa
piler -cons ht.plat.afa -out ht.plat.cons.fa -label ht.plat.gb.cons

# Use RepeatMasker to calculate divergence rate for platypus-like tuatara L2 and non-platypus like tuatara L2
RepeatMasker -pa 16 -a -nolow -lib ht.plat.cons.fa ht.plat.fa.use
RepeatMasker -pa 16 -a -nolow -lib ht.nonplat.cons.fa ht.nonplat.fa.use
/usr/local/src/RepeatMasker_latest/util/calcDivergenceFromAlign.pl -s ht.plat.divsum ht.plat.fa.use.align 
/usr/local/src/RepeatMasker_latest/util/calcDivergenceFromAlign.pl -s ht.nonplat.divsum ht.nonplat.fa.use.align 

# Combine two divergence output
paste ht.plat.divsum ht.nonplat.divsum | column -s $'\t' -t > tmp
awk '{print $1 "\t" $2 "\t" $4}' tmp > plat_nonplat_divsum