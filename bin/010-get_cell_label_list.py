#!/usr/bin/env python


import argparse
import h5py
import numpy as np

# read_elem has moved between anndata versions
try:
    from anndata.io import read_elem  # anndata >= 0.11
except ImportError:
    try:
        from anndata.experimental import read_elem  # anndata 0.8-0.10
    except ImportError:
        from anndata._io.specs import read_elem  # older


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Get a list of all celltypes in an anndataframe.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '--cell_label',
        action='store',
        dest='cell_label',
        default='cluster',
        help='Anndata cell type label name in obs slot. (default: %(default)s)'
    )

    options = parser.parse_args()

    # Read only obs (X, layers, obsm etc. are never loaded into memory)
    with h5py.File(options.h5, "r") as f:
        obs = read_elem(f["obs"])

    if options.cell_label not in obs.columns:
        raise KeyError(
            f"{options.cell_label} not in obs. Available columns: {list(obs.columns)}"
        )

    np_array = np.sort(obs[options.cell_label].unique().astype(str))
    np.savetxt('cell_labels.csv', np_array, fmt='%s')


if __name__ == '__main__':
    main()
