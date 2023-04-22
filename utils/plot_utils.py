"""Utils for plotting"""
import matplotlib.pyplot as plt
import numpy as np

import constants


def plotPrmRgbOnCt(ctArray: np.ndarray, prmAllArray: np.ndarray, path: str):
    """Plot single slice of RGB color array overlaid on HRCT.

    Args:
        ctArray (np.array): HRCRT image to overlay PRGM RGB image onto
        prmAllArray (np.array): PRM image sorted into bins denoting PRM voxel classification
        path (str): path to save final image to
    """

    plt.figure()
    plt.imshow(prmRgbArray[:, constants.prmProcessing.PLOT_SLICENUM, :, :], cmap="gray")
    plt.axis("off")
    plt.savefig(path, transparent=True, bbox_inches="tight", pad_inches=-0.05, dpi=600)
    plt.clf()
    plt.close()


def plotPrmRgbImage2(expArray, prmAllArray, path):
    cmap = colors.ListedColormap(["black", "green", "yellow", "red", "purple"])
    bounds = [0, 1, 2, 3, 4, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    prmAllArraySlice = prmAllArray[:, constants.prmProcessing.PLOT_SLICENUM, :]
    prmAllArraySlice = np.ma.masked_where(prmAllArraySlice < 1, prmAllArraySlice)
    plt.figure()
    plt.imshow(expArray[:, constants.prmProcessing.PLOT_SLICENUM, :], cmap="gray")
    plt.imshow(prmAllArraySlice, cmap=cmap, norm=norm, alpha=0.7)
    plt.axis("off")
    plt.savefig(path, transparent=True, bbox_inches="tight", pad_inches=-0.05, dpi=600)
    plt.clf()
    plt.close()
