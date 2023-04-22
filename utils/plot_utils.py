"""Utils for plotting"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

import constants


def plotPrmRgbOnCt(
    ctArray: np.ndarray, prmAllArray: np.ndarray, sliceNum: int, path: str
):
    """Plot single slice of RGB color array overlaid on HRCT.

    Args:
        ctArray (np.array): HRCRT image to overlay PRGM RGB image onto
        prmAllArray (np.array): PRM image sorted into bins denoting PRM voxel classification
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        path (str): path to save final image to
    """

    # create custom discrete colormap and bounds
    cmap = colors.ListedColormap(["black", "green", "yellow", "red", "purple"])
    bounds = [0, 1, 2, 3, 4, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # take single slice of PRM binned image and mask out un-binned regions
    prmAllArraySlice = prmAllArray[:, sliceNum, :]
    prmAllArraySlice = np.ma.masked_where(prmAllArraySlice < 1, prmAllArraySlice)

    # plot image overlaid on corresponding HRCT slice
    plt.figure()
    plt.imshow(ctArray[:, sliceNum, :], cmap="gray")
    plt.imshow(prmAllArraySlice, cmap=cmap, norm=norm, alpha=0.7)
    plt.axis("off")
    plt.savefig(path, transparent=True, bbox_inches="tight", pad_inches=-0.05, dpi=600)
    plt.clf()
    plt.close()
