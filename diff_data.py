import pandas as pd
import click as ck

@ck.command()
def main():
    """
    Script for preparing and analyzing differential expression data.
    """
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    control_columns = ['CRC-10-ENAS-P8',]
    treatment_columns = ['CRC-10-ENA-P8']

if __name__ == "__main__":
    main()