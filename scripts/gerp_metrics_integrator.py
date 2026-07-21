import numpy as np
import argparse
from pathlib import Path
import pyBigWig
import pandas as pd

BW_GERP = pyBigWig.open("/home/agenor/phd/projeto/gerp_conservation_scores.homo_sapiens.GRCh38.bw")
HG_CHR = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y"
]

# g_coords = pd.read_csv("/home/agenor/phd/projeto/output/results_merged.csv")

def integrator(gc):
    '''
    Integration of bigWig file with GERP++ cpnservation scores and protein domain genomic coordinates.
    It also generates some summary statistics.
    '''

    g_coords = pd.read_csv(gc)
    
    chr_coords = g_coords['genomic_coordinates'].str.split(':', expand=True)[0].astype(int, errors="ignore")
    pos_coords = g_coords['genomic_coordinates'].str.split(':', expand=True)[1].str.split('-', expand=True).astype(int)
    pos_coords.rename({0:"start", 1:"end"}, axis=1, inplace=True)

    index_g_coord = pos_coords.merge(chr_coords, left_index=True, right_index=True)
    index_g_coord.rename({0:"chr"}, axis=1, inplace=True)

    index_g_coord["start_sort"] = np.minimum(index_g_coord["start"], index_g_coord["end"])
    index_g_coord["end_sort"] = np.maximum(index_g_coord["start"], index_g_coord["end"])

    gerp_means = []
    for r in index_g_coord.itertuples():
        if r.chr in HG_CHR:
            gerp_means.append(BW_GERP.values(f"{r.chr}", r.start_sort, r.end_sort, numpy=True).mean())
        else:
            gerp_means.append(np.nan)

    g_coords["gerp_mean"] = gerp_means

    return g_coords

    

def init_argparser():
    '''
    Initializes an argument parser
    '''

    parser = argparse.ArgumentParser(
        usage='%(prog)s',
        description='Integrates GERP++ conservation scores with protein domain genomic coordinates', 
    )

    parser.add_argument('-gc', '--genomic_coordinates', nargs=1)
    parser.add_argument('-out_dir', '--output_directory', nargs=1)

    return parser


def main() -> None:

    parser = init_argparser()
    args = parser.parse_args()

    print("Initializing")

    final_result = integrator(f"{args.genomic_coordinates[0]}")
    final_result.to_csv(Path(f"{args.output_directory[0]}","gerp_integration.csv"))


if __name__ == '__main__':
    main()