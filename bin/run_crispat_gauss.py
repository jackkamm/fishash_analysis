#!/usr/bin/env python

import argparse
import crispat
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run crispat gauss model')
    parser.add_argument("in_h5ad")
    parser.add_argument("out_prefix")
    args = parser.parse_args()

    crispat.ga_gauss(args.in_h5ad,
                     args.out_prefix,
                     nonzero=True)
