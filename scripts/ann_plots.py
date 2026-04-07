#!/opt/miniconda3/bin/python3


from cyvcf2 import VCF
import pandas as pd
import matplotlib.pyplot as plt

vcf_file = "/home/biouser/learning/variant_calling/HCC/results/annotation/HCC_filtered_indels_extracted_ann.vcf.gz"  # bgzipped recommended

vcf = VCF(vcf_file)

records = []

for variant in vcf:
    chrom = variant.CHROM
    pos = variant.POS
    
    ann = variant.INFO.get('ANN')
    
    if ann is None:
        continue
    
    # ANN can have multiple transcripts
    annotations = ann.split(',')
    
    for entry in annotations:
        fields = entry.split('|')
        
        if len(fields) < 4:
            continue
        
        effect = fields[1]
        impact = fields[2]
        gene = fields[3]
        
        records.append({
            "CHROM": chrom,
            "POS": pos,
            "GENE": gene,
            "EFFECT": effect,
            "IMPACT": impact
        })

df = pd.DataFrame(records)

# -------------------------
# 🔥 Filter important variants
# -------------------------
df = df[df["IMPACT"].isin(["HIGH", "MODERATE"])]

# -------------------------
# 📊 Gene frequency
# -------------------------
top_genes = df["GENE"].value_counts().head(10)

top_genes.plot(kind="bar")
plt.title("Top Mutated Genes")
plt.xlabel("Gene")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("/home/biouser/learning/variant_calling/HCC/results/annotation/plots/gene_frequency.png")
plt.clf()

# -------------------------
# 📊 Impact distribution
# -------------------------
df["IMPACT"].value_counts().plot(kind="bar")
plt.title("Impact Distribution")
plt.xlabel("Impact")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("/home/biouser/learning/variant_calling/HCC/results/annotation/plots/impact_distribution.png")
plt.clf()

# -------------------------
# 📊 Variant type distribution
# -------------------------
df["EFFECT"].value_counts().head(10).plot(kind="bar")
plt.title("Variant Type Distribution")
plt.xlabel("Effect")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("/home/biouser/learning/variant_calling/HCC/results/annotation/plots/variant_types.png")
plt.clf()

# -------------------------
# 💾 Save filtered data
# -------------------------
df.to_csv("/home/biouser/learning/variant_calling/HCC/results/annotation/plots/significant_indels.csv", index=False)

print("Analysis complete. Files saved.")