#!/usr/bin/env python

import argparse

import scipy.io
import scipy.sparse

import anndata as ad
import pegasusio as pgio
import demuxEM

def main():
    parser = argparse.ArgumentParser(description='Run demuxEM')
    parser.add_argument("in_h5ad")
    parser.add_argument("out_mtx")
    parser.add_argument("--min_signal", type=float, default=2)
    parser.add_argument("--n_threads", type=int, default=1)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.in_h5ad)

    adata = pgio.MultimodalData(adata)
    adata2 = adata.copy()

    demuxEM.estimate_background_probs(adata2)
    demuxEM.demultiplex(adata, adata2,
                        min_signal=args.min_signal,
                        n_threads=args.n_threads)

    assignments_mat = get_assignments(adata)
    scipy.io.mmwrite(args.out_mtx, assignments_mat)

def get_assignments(adata):
    ret = scipy.sparse.dok_array(adata.shape, dtype=int)
    feature2idx = dict(zip(adata.var.index, range(adata.shape[1])))

    for i, (_, row) in enumerate(adata.obs.iterrows()):
        if row['demux_type'] == 'unknown' or row['assignment'] == '':
            continue

        for feature in row['assignment'].split(','):
            ret[i, feature2idx[feature]] = 1

    return ret.T

if __name__ == '__main__':
    main()
