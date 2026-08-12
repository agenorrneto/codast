import pandas as pd
import argparse
import subprocess
import requests, sys
import gget
from pathlib import Path
from tqdm import tqdm

"""
InterProSCan result
0 - Protein accession (e.g. P51587)
1 - Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
2 - Sequence length (e.g. 3418)
3 - Analysis (e.g. Pfam / PRINTS / Gene3D)
4 - Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
5 - Signature description (e.g. BRCA2 repeat profile)
6 - Start location
7 - Stop location
8 - Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
9 - Status - is the status of the match (T: true)
10 - Date - is the date of the run
11 - InterPro annotations - accession (e.g. IPR002093)
12 - InterPro annotations - description (e.g. BRCA2 repeat)
13 - GO annotations with their source(s), e.g. GO:0005515(InterPro)|GO:0006302(PANTHER)|GO:0007195(InterPro,PANTHER). This is an optional column; only displayed if the --goterms option is switched on
14 - Pathways annotations, e.g. REACT_71. This is an optional column; only displayed if the --pathways option is switched on
"""

SERVER = "https://www.ebi.ac.uk/proteins/api"
EXT_BASIS = "/coordinates/location/"
MANE = pd.read_csv("/home/agenor/phd/projeto/MANE.GRCh38.v1.4.summary.txt", sep ="\t", header=0)

def filter_domain(protein_list, database):
    """
    
    """
    ip_res = pd.read_csv(protein_list, sep="\t", header=None)
    ip_res = ip_res[ip_res[3] == f"{database}"]
    ip_domain = ip_res[ip_res[5].str.contains("domain", case=False, na=False)]
    print("filtered")

    return ip_domain

def gcoord_translator(ip_result, out_directory):

    """
    
    """
    results = []

    for i in tqdm(range(len(ip_result)), desc="Processing rows"):
        try:
            ensembl_id = ip_result.iloc[i][0].replace("9606.", "")
            domain_annottation = str(ip_result.iloc[i][12])
            domain_range = [int(ip_result.iloc[i][6]), int(ip_result.iloc[i][7])]

            #get info
            id_info = gget.info(f"{ensembl_id}")
            uniprot_id = id_info["uniprot_id"].to_string(index=False)

            ext = EXT_BASIS +f"{uniprot_id}:{domain_range[0]}-{domain_range[1]}?format=json"
            r = requests.get(SERVER+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decoded = r.json()

            #Calculate gene lenght and gene_domain lenght / gene lenght ratio
            MANE["NCBI_GeneID"] = MANE.NCBI_GeneID.str.replace("GeneID:", "")

            gene_start = int(MANE[MANE["NCBI_GeneID"] == id_info.loc[ensembl_id, "ncbi_gene_id"]].chr_start.values[0])
            gene_end = int(MANE[MANE["NCBI_GeneID"] == id_info.loc[ensembl_id, "ncbi_gene_id"]].chr_end.values[0])

            gene_range = gene_start - gene_end
            if gene_range < 0:
                gene_range = abs(gene_range)
                gene_range = gene_range + 1
            else:
                gene_range = gene_range + 1

            domain_range_nt = int(decoded['locations'][0]["geneStart"]) - int(decoded['locations'][0]["geneEnd"])
            if domain_range_nt < 0:
                domain_range_nt = abs(domain_range_nt)
                domain_range_nt = domain_range_nt + 1
            else:
                domain_range_nt = domain_range_nt + 1

            dom_gene_ratio = domain_range_nt / gene_range



            results.append({
                        "ensemblTrans_id": ensembl_id,
                        "ensembleGene_id": decoded['locations'][0]["ensemblGeneId"],
                        "uniprot_id": uniprot_id,
                        "gene_name": id_info.loc[ensembl_id,"primary_gene_name"],
                        "ncbi_gene_id": id_info.loc[ensembl_id, "ncbi_gene_id"],
                        "genomic_coordinates": str(decoded['locations'][0]["chromosome"]) + ":" + str(decoded['locations'][0]["geneStart"]) + "-" + str(decoded['locations'][0]["geneEnd"]),
                        "gene_lenght": gene_range,
                        "gene_star": gene_start,
                        "gene_end": gene_end,
                        "domain_start": domain_range[0],
                        "domain_end": domain_range[1],
                        "pfam_domain": ip_result.iloc[i][4],
                        "domain_range_nt": domain_range_nt,
                        "dom_gene_ratio": dom_gene_ratio,
                        "domain_annotations": domain_annottation
                    })
            print("it worked!")
        except Exception as e:
            print(f"Error processing row {i} (ensemble_id={ensembl_id}): {e}")
            continue
    
    final_df = pd.DataFrame(results)

    return final_df

        
    "Q9NZT1:81-142?format=json"





def init_argparser():
    '''
    Initializes an argument parser
    '''

    parser = argparse.ArgumentParser(
        usage='%(prog)s',
        description='Get genomics coordinates for annottaded domains', 
    )

    parser.add_argument('-ip', '--interpro_result', nargs=1)
    parser.add_argument('-db', '--database', nargs=1)
    parser.add_argument('-out_dir', '--output_directory', nargs=1)

    return parser

def main() -> None:
    
    """
    Main function for [program's name].
    """

    parser = init_argparser()
    args = parser.parse_args()

    print("initializing")


    filtered = filter_domain(f"{args.interpro_result[0]}", 
                             f"{args.database[0]}")
    final_result = gcoord_translator(filtered, 
                              f"{args.output_directory[0]}")
    final_result.to_csv(Path(f"{args.output_directory[0]}","genomic_coordinates_batch_2.csv"))


if __name__ == '__main__':
    main()