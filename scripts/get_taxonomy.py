import time
import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa

Entrez.email = "kamikiritsu@google.com"  # 请改成你的真实邮箱

input_file = "/Users/kamikiritsu/Desktop/GH7 origin/GH7_metadata.xlsx"
output_file = "/Users/kamikiritsu/Desktop/GH7 origin/GH7_metadata_enriched_v2.xlsx"

print("Initializing local NCBI taxonomy database...")
ncbi = NCBITaxa()

df = pd.read_excel(input_file)
accessions = df["Accession"].tolist()

def get_taxonomy(acc):
    time.sleep(0.5)
    print(f"Processing {acc}")
    try:
        handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        taxid = None
        for rec in records:
            if "GBSeq_taxonomy" in rec:
                taxid = rec["GBSeq_taxonomy"]
                break
        if not taxid:
            print(f"  -> No taxid found for {acc}")
            return {"Phylum": None, "Class": None, "Order": None, "Family": None}
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        res = {"Phylum": None, "Class": None, "Order": None, "Family": None}
        for tid in lineage:
            rank = ranks[tid]
            name = names[tid]
            if rank == "phylum":
                res["Phylum"] = name
            elif rank == "class":
                res["Class"] = name
            elif rank == "order":
                res["Order"] = name
            elif rank == "family":
                res["Family"] = name
        return res
    except Exception as e:
        print(f"  -> Error on {acc}: {e}")
        return {"Phylum": None, "Class": None, "Order": None, "Family": None}

results = [get_taxonomy(acc) for acc in accessions]
taxa_df = pd.DataFrame(results)
df_enriched = pd.concat([df, taxa_df], axis=1)
df_enriched.to_excel(output_file, index=False)
print(f"\nDone! Enriched table saved to: {output_file}")
