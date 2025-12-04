#!/usr/bin/env python

import argparse
import crispat
import pandas as pd
import anndata as ad

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run crispat poisson model')
    parser.add_argument("in_h5ad")
    parser.add_argument("out_prefix")
    # crispat will try to use all CPUs on the node (more than LSF
    # allocated) unless we manually specify this
    parser.add_argument("--cpus", type=int, required=True, default=0)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.in_h5ad)
    # avoid crashing due to subsample size bigger than anndata. The
    # default value is 15000
    subsample_size = min(15000, round(adata.shape[0] / 2))
    del adata

    if args.cpus > 1:
        nproc = args.cpus
    else:
        nproc = None

    crispat.ga_poisson(
        args.in_h5ad,
        args.out_prefix,
        subsample_size=subsample_size,
        parallelize=args.cpus > 1,
        n_processes = nproc
    )
