vcf="..." # input VCF file with phased GT, LRR, and BAF
pfx="..." # output prefix
thr="..." # number of threads to use
crt="..." # tab delimited file with call rate information (first column sample ID, second column call rate)
sex="..." # tab delimited file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="..." # directory where output files will be generated
mkdir -p $dir
