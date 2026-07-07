import numpy as np
import argparse
from pathlib import Path
import pyBigWig
import pandas as pd

BW_GERP = pyBigWig.open("/home/agenor/phd/projeto/gerp_conservation_scores.homo_sapiens.GRCh38.bw")

# g_coords = pd.read_csv("/home/agenor/phd/projeto/output/results_merged.csv")

def integrator(gc, output):
    '''
    Integration of bigWig file with GERP++ cpnservation scores and protein domain genomic coordinates.
    It also generates some summary statistics.
    '''

    g_coords = pd.read_csv(gc)
    
    chr_coords = g_coords['genomic_coordinates'].str.split(':', expand=True)[0].astype(int)
    pos_coords = g_coords['genomic_coordinates'].str.split(':', expand=True)[1].str.split('-', expand=True).astype(int)


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



    final_result = integrator(gc, 
                              f"{args.output_directory[0]}")
    final_result.to_csv(Path(f"{args.output_directory[0]}","gerp_integration.csv"))


if __name__ == '__main__':
    main()