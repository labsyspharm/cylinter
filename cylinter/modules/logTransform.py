import numpy as np

from ..utils import input_check, read_markers, reorganize_dfcolumns


def logTransform(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    abx_channels_mod = data[abx_channels].copy()
    abx_channels_mod = np.log10(abx_channels_mod + 0.00000000001)
    data.loc[:, abx_channels] = abx_channels_mod

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
