import numpy as np

from ..utils import input_check, read_markers, reorganize_dfcolumns


def logTransform(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    abx_channels_mod = data[abx_channels].copy()
    abx_channels_mod = np.log10(abx_channels_mod + 0.001)
    data.loc[:, abx_channels] = abx_channels_mod
    
    # clip cells with zero-valued signal intensities to the Nth percentile of the
    # distribution (not considering the zero-valued signals themselves).     
    # percentile = 5
    # percentiles = (
    #     data.loc[:, abx_channels][data.loc[:, abx_channels] > 0.0].quantile(q=percentile / 100)
    # )
    # data.loc[:, abx_channels] = (
    #     data.loc[:, abx_channels].clip(lower=percentiles, upper=None, axis=1)
    # )

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
