#Get just CEPH samples
awk -F "\t" '{if ($4 == "CEU") print $1}' igsr_samples.tsv > CEPH_samples.txt

#Extract CEPH samples from each chromosome of 1000 genomes.
for i in `seq 1 22; echo X; echo Y`; do echo "vcfkeepsamples ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $(cat european_samples.txt) > EUR.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"; done;

#Use plink to generate a triangular matrix of SNP LD scores, chromosome by chromosome.
for i in `seq 1 22`; do echo "bsub -J ${i} plink2 --vcf CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --r gz triangle yes-really --out CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes" ; done;

#Generate list of duplicate rsIDs
for i in `seq 22`; do echo "grep -v '#' CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf | cut -f3 | uniq -cd | awk '{ print $2 }' > chr${i}_dup_rs.txt; echo '.' >> chr${i}_dup_rs.txt" | bsub -J ${i}duplist -o ${i}duplist.out -T 8; done;

#Exclude duplicate rsIDs in the VCF files
for i in `seq 22`; do echo "vcftools --vcf CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --exclude chr${i}_dup_rs.txt --recode --stdout > CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nodup.vcf" | bsub -J ${i}nodup -o ${i}nodup.out -T 8; done;


#Generate plink files from VCF files without duplicates
for i in `seq 22`; do echo "vcftools --vcf CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nodup.vcf --plink --out > CEPH.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nodup" | bsub -J ${i}plink -o ${i}plink.out -T 16; done;


#Generate R values given a list of SNPs
plink2 --vcf <(vcftools --vcf CEPH.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --snps rsIDs.txt --recode --stdout) --r square --out testLD

