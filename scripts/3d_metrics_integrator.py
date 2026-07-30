import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import gget
import requests
import os
from Bio.PDB.DSSP import make_dssp_dict
import time

RES_MAX_ACC = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
               'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
               'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
               'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
               'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}

SS_MAP = {"H": "H",
"B": "B",
 "E":"B",
 "G":"H",
 "I":"H",
 "P":"H",
 "T":"C",
 "S": "B",
 "-": "C"}


def download_af_model(uniprot_id, output_dir):
    # Variables for API
    file_name = f"AF-{uniprot_id}-F1-model_v6.cif"
    url = f"https://alphafold.ebi.ac.uk/files/{file_name}"
    output_path = os.path.join(output_dir, file_name)

    # Check if the file exists
    if os.path.exists(output_path):
        print("File already exists")
    else:

        print(f"Attempting to download: {file_name}...")

        # Request
        response = requests.get(url, stream=True)
        
        # Status Code 200 = OK
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Success! Saved to: {output_path}")
        else:
            print(f"Error {response.status_code}: Could not find model for ID {uniprot_id}")


def run_dssp(cif_file, output_dir):
    # Conf
    URL_PDB = "https://pdb-redo.eu/dssp/do"
    input_file_path = os.path.join(output_dir, cif_file) 
    output_file = f"{cif_file}.dssp"
    output_path = os.path.join(output_dir, output_file)                  

    if os.path.exists(output_path):
        print("DSSP file already exists")
    else:
        # Prepare the payload
        form_data = {
            'format': 'dssp'
        }

        # Open the file and send the POST request
        print(f"Sending {cif_file} to DSSP server...")

        with open(input_file_path, 'rb') as f:
            # 'files' matches the file upload (-F "data=<filename")
            # The key 'data' corresponds to the field name the server expects
            files_payload = {'data': f}

            try:
                response = requests.post(URL_PDB, data=form_data, files=files_payload)

                # Status Code 200 = OK
                if response.status_code == 200:
                    # Write the result to the output file
                    with open(output_path, 'wb') as out_f:
                        out_f.write(response.content)
                    print(f"Success! Output saved to: {output_path}")
                else:
                    print(f"Error {response.status_code}: {response.text}")
            except Exception as e:
                print(f"An error occurred: {e}")


def process_dssp(file_path):
    # make_dssp_dict parses the file into a dictionary
    # Key = (chain, res_id), Value = (aa, sec_struct, acc)
    parsed_data = []
    if os.path.exists(file_path):
        dssp_dict, keys = make_dssp_dict(file_path)

        for key in keys:
            chain_id = key[0]
            res_id   = key[1][1] # BioPython residue ID tuple

            # dssp_dict returns a tuple.
            # Index 0 = Amino Acid
            # Index 1 = Secondary Structure
            # Index 2 = ACC (Accessibility)
            data = dssp_dict[key]

            parsed_data.append({
                'chain': chain_id,
                'residue': res_id,
                'AA': data[0],
                'SS': data[1],
                'ACC': data[2]
            })
    else:
        parsed_data.append({
                'chain': np.nan,
                'residue': "",
                'AA': "",
                'SS': "",
                'ACC': ""
            })
        

    return parsed_data


def integrator_3d(data, raw):
    """
    Docstring for integrator_3d
    """
    class_3d_list = []
    df = pd.read_csv(data)

    for uid in df["uniprot_id"]:
        # Downloading AlphaFold for each protein
        download_af_model(uid, output_dir=raw)

        # Checking the existence of the AlphaFold model and running DSSP
        af_file = f'AF-{uid}-F1-model_v6.cif'
        file_path = os.path.join(raw, af_file)
        if os.path.exists(file_path):
            print(f"AlphaFold model for {uid} already exists")
            try:
                run_dssp(af_file, output_dir=raw)
            except:
                print(f"Error running DSSP for {uid}")
        else:
            print(f"AlphaFold model for {uid} not found")

    # Try to reduce the loops

    for i in range(len(df)):
        id = df.iloc[i]["uniprot_id"]
        aa_start = int(df.iloc[i]["aa_start"])
        aa_end = int(df.iloc[i]["aa_end"])

        info_3d = process_dssp(f"{raw}/AF-{id}-F1-model_v6.cif.dssp")

        try:
            domain_3d_data = pd.DataFrame(info_3d[aa_start - 1:aa_end])
            domain_3d_data = domain_3d_data[['residue', 'AA', 'SS', 'ACC']]

            domain_3d_data["max_acc"] = domain_3d_data["AA"].map(RES_MAX_ACC)
            domain_3d_data["rsa"] = domain_3d_data["ACC"] / domain_3d_data["max_acc"]
            domain_3d_data["ss_3class"] = domain_3d_data["SS"].map(SS_MAP)
            domain_3d_data["burial_status"] = np.where(domain_3d_data["rsa"] > 0.25, "E", "B")

            proportions_rsa = domain_3d_data.burial_status.value_counts(normalize=True)
            proportions_ss = domain_3d_data.ss_3class.value_counts(normalize=True)

            if proportions_rsa.max() >= 0.7:
                dominant_burial_status = proportions_rsa.idxmax()
            else:
                dominant_burial_status = "MEB" # Mixed Exposed/Buried
            
            if proportions_ss.max() >= 0.7:
                dominant_ss = proportions_ss.idxmax()
            else:
                dominant_ss = "MSS" # Mixed Secondary Structure

            class_3d_list.append({
                'dominant_burial_status': dominant_burial_status,
                'dominant_ss': dominant_ss
            })

        except IndexError:
            print(f'IndexError: {aa_start} to {aa_end} is out of range of {id} model (range: {len(info_3d)})')

            class_3d_list.append({
                'dominant_burial_status': "out of range",
                'dominant_ss': "out of range"
            })

    
    class_3d_df = pd.DataFrame(class_3d_list)

    class_3d_df = class_3d_df.reset_index()
    print(class_3d_df)
    final_df = pd.concat([df, class_3d_df], axis=1)

    return final_df



def init_argparser():
    '''
    Initializes an argument parser
    '''

    parser = argparse.ArgumentParser(
        usage='%(prog)s',
        description='Integrating 3D data for protein domains', 
    )

    parser.add_argument('-gc', '--genomic_coordinates')
    parser.add_argument('-out_dir', '--output_directory')
    parser.add_argument('-raw_dir', '--rawdata_directory')

    return parser

def main() -> None:
    
    """
    Main function for 3d_metrics_integrator.
    """

    parser = init_argparser()
    args = parser.parse_args()

    print("Initializing")

    threed_df = integrator_3d(f"{args.genomic_coordinates}", 
                            raw=f"{args.rawdata_directory}")
    threed_df.to_csv(Path(f"{args.output_directory}","domain_3d_metrics.csv"))


if __name__ == '__main__':
    main()