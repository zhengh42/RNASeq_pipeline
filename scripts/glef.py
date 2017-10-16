#! /bin/env python
import argparse
import pandas as pd
import numpy as np
import sys

def main(args):
    gtable = pd.read_table(args.ginput).set_index('Name')
    ttable = pd.read_table(args.tinput).set_index('Name')
    tgmap = pd.read_table(args.tgmap, names=['t', 'g']).set_index('t')
    gene_lengths = {}
    j = 0

    # Map over all gene groups (a gene and its associated transcripts)
    for g, txps in tgmap.groupby('g').groups.iteritems():
        if j % 500 == 1:
            print("Processed {} genes".format(j))
        j += 1
        # The set of transcripts present in our salmon index
        tset = []
        for t in txps:
            if t in ttable.index:
                tset.append(t)
        # If at least one of the transcripts was present
        if len(tset) > 0:
            # The denominator is the sum of all TPMs
            totlen = ttable.loc[tset,'TPM'].sum()
            # Turn the relative TPMs into a proper partition of unity
            if totlen > 0:
                tpm_fracs = ttable.loc[tset, 'TPM'].values / ttable.loc[tset,'TPM'].sum()
            else:
                tpm_fracs = np.ones(len(tset)) / float(len(tset))
            # Compute the gene's effective length as the abundance-weight
            # sum of the transcript lengths
            elen = 0.0
            for i,t in enumerate(tset):
                elen += tpm_fracs[i] * ttable.loc[t, 'EffectiveLength']
            gene_lengths[g] = elen

    # Give the table an effective length field
    gtable['EffectiveLength'] = gtable.apply(lambda r : gene_lengths[r.name] if r.name in gene_lengths else 1.0, axis=1)
    # Write it to the output file
    gtable.to_csv(args.output, sep='\t', index_label='Name')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute gene-level effective lengths")
    parser.add_argument('--ginput', type=str, help='gene level input table')
    parser.add_argument('--tinput', type=str, help='transcript level input table')
    parser.add_argument('--tgmap', type=str, help='transcript -> gene mapping')
    parser.add_argument('--output', type=str, help='output table with extra column')

    args = parser.parse_args()
    main(args)

