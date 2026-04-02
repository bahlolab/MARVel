module load bcftools/1.23 
cd /vast/scratch/users/wang.lo/1000G/PCA
###################################################
# one sample test
vcf=/vast/projects/bahlo_mtDNA/Results_1000G/HG00096/mitoHPC/mutect2.mutect2.03.merge.vcf
sample=HG00096
bcftools view $vcf | grep -v "^#" | wc -l
bcftools view $vcf -v snps | grep -v "^#" | wc -l
bcftools view $vcf -i 'FORMAT/AF>=0.95' | grep -v "^#" | wc -l

bcftools view $vcf \
    -v snps -m2 -M2 \
| bcftools view -i 'FORMAT/AF>=0.95' \
-Oz -o ${sample}.clean.vcf.gz

bcftools index -t ${sample}.clean.vcf.gz
bcftools view ${sample}.clean.vcf.gz | grep -v "^#" | wc -l
bcftools view ${sample}.clean.vcf.gz | head
bcftools view ${sample}.clean.vcf.gz | grep 14569


###################################################
# STEP 1 — Create a list of all 1000G sample VCFs
ls /vast/projects/bahlo_mtDNA/Results_1000G/*/mitoHPC/mutect2.mutect2.03.merge.vcf \
> 1000G.mtDNA.vcf.list

# STEP 2 — Clean each sample VCF (important)
mkdir clean_vcfs

while read vcf; do
    sample=$(basename "$(dirname "$(dirname "$vcf")")")
    echo $sample
    
    bcftools view $vcf \
        -v snps -m2 -M2 \
    | bcftools view -i 'FORMAT/AF>=0.95' \
    -Oz -o clean_vcfs/${sample}.clean.vcf.gz

    bcftools index -t clean_vcfs/${sample}.clean.vcf.gz
    
done < 1000G.mtDNA.vcf.list

# STEP 3 — Merge all samples into one multi-sample VCF
ls clean_vcfs/*.vcf.gz > merge.list

bcftools merge \
    -m none \
    -Oz \
    -o 1000G.mtDNA.homoplasmy.merged.vcf.gz \
    -l merge.list

bcftools index -t 1000G.mtDNA.homoplasmy.merged.vcf.gz

bcftools merge \
  --missing-to-ref \
  -m none \
  -Oz -o 1000G.mtDNA.homoplasmy.merged.miss2ref.vcf.gz \
  -l merge.list

bcftools index -t 1000G.mtDNA.homoplasmy.merged.miss2ref.vcf.gz

# STEP 4 — Remove D-loop 
printf "chrM\t1\t576\nchrM\t16024\t16569\n" > dloop.hg38.bed

bcftools view \
    -T ^dloop.hg38.bed \
    1000G.mtDNA.homoplasmy.merged.miss2ref.vcf.gz \
    -Oz -o 1000G.mtDNA.noDloop.vcf.gz

bcftools index -t 1000G.mtDNA.noDloop.vcf.gz

# remove related samples
bcftools query -l 1000G.mtDNA.noDloop.vcf.gz | sort > vcf.samples.final.txt
ped=/vast/projects/bahlo_mtDNA/1000G/20130606_g1k.ped
awk '$3=="0" && $4=="0" {print $2}' $ped | sort > ped.founders.base.txt
sed 's/\.final$//' vcf.samples.final.txt | sort > vcf.samples.base.txt

comm -12 ped.founders.base.txt vcf.samples.base.txt > keep.base.txt

awk '{print $0".final"}' keep.base.txt | sort > keep.final.txt

bcftools view -S keep.final.txt \
  1000G.mtDNA.noDloop.vcf.gz \
  -Oz -o 1000G.mtDNA.unrelated.vcf.gz

bcftools index -t 1000G.mtDNA.unrelated.vcf.gz
# 2589

# STEP 5 — Build stable PCA panel (common sites only)
bcftools +fill-tags \
    1000G.mtDNA.unrelated.vcf.gz \
    -Oz -o 1000G.mtDNA.tags.vcf.gz \
    -- -t AF,F_MISSING
    
bcftools index -t 1000G.mtDNA.tags.vcf.gz
bcftools view 1000G.mtDNA.tags.vcf.gz | grep -v "^#" | wc -l

bcftools query -f '%POS\t%INFO/AF\t%INFO/F_MISSING\n' \
  1000G.mtDNA.tags.vcf.gz | head
  
bcftools view \
    -i 'INFO/AF>=0.01 && INFO/AF<=0.99 && F_MISSING<=0.01' \
    1000G.mtDNA.tags.vcf.gz \
    -Oz -o 1000G.mtPCA.panel.vcf.gz

bcftools index -t 1000G.mtPCA.panel.vcf.gz
bcftools view 1000G.mtPCA.panel.vcf.gz | grep -v "^#" | wc -l
# 317
bcftools query -f '%POS\t%INFO/AF\t%INFO/F_MISSING\n' \
  1000G.mtPCA.panel.vcf.gz | head
  
bcftools query -f '%INFO/AF\n' 1000G.mtPCA.panel.vcf.gz | head

# STEP 6 — Export genotype matrix
bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' \
    1000G.mtPCA.panel.vcf.gz \
    > 1000G.mtPCA.GT.tsv

{
  echo -ne "Variant"
  bcftools query -l 1000G.mtPCA.panel.vcf.gz | \
    sed 's/\.final$//' | \
    awk '{printf "\t%s", $0}'
  echo

  bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' \
    1000G.mtPCA.panel.vcf.gz
} > 1000G.mtPCA.GT.withHeader.tsv

bcftools query -l 1000G.mtPCA.panel.vcf.gz | \
  sed 's/\.final$//' > mtPCA.sample_ids.txt
  
# STEP 7 — Compute PCA loadings (R)

