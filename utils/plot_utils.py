"""Utils for plotting"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

import constants


def findPlotSliceNum(maskArray: np.ndarray):
    """Find optimal slice number to use for demonstrative plots of thoracic cavity.

    Args:
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity

    Returns:
        sliceNum (int): slice number along anterior-posterior dimension

    Finds slice along anterior-posterior dimension with largest thoracic cavity volume.
    """

    # create array with mask size of each slice
    numSlices = maskArray.shape[1]
    maskSizeBySlice = np.zeros(numSlices)
    for i in range(numSlices):
        maskSizeBySlice[i] = (maskArray[:, i, :] > 0).sum()

    sliceNum = np.argmax(maskSizeBySlice)

    return sliceNum


def plotPrmRgbOnCt(
    ctArray: np.ndarray, prmAllArray: np.ndarray, sliceNum: int, path: str
):
    """Plot single slice of RGB color array overlaid on HRCT.

    Args:
        ctArray (np.array): HRCT image to overlay PRGM RGB image onto
        prmAllArray (np.array): PRM image sorted into bins denoting PRM voxel classification
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        path (str): path to save final image to
    """

    # create custom discrete colormap and bounds
    cmap = colors.ListedColormap(["green", "yellow", "red", "purple"])
    bounds = [1, 2, 3, 4, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # take single slice of PRM binned image and mask out un-binned regions
    prmAllArraySlice = prmAllArray[:, sliceNum, :]
    prmAllArraySlice = np.ma.masked_where(prmAllArraySlice < 1, prmAllArraySlice)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
    plt.imshow(ctArray[:, sliceNum, :], cmap="gray")
    plt.imshow(prmAllArraySlice, cmap=cmap, norm=norm)
    plt.axis("off")

    # add colorbar
    imgRatio = prmAllArraySlice.shape[0] / prmAllArraySlice.shape[1]
    cbar = plt.colorbar(ticks=[1.5, 2.5, 3.5, 4.5], fraction=0.047 * imgRatio, pad=0)
    cbar.ax.set_yticklabels(["norm", "fSAD", "emph", "empt\nemph"])

    # save figure and close it
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=600)
    plt.clf()
    plt.close()


def plotTopoOnCt(
    ctArray: np.ndarray,
    maskArray: np.ndarray,
    topoArray: np.ndarray,
    sliceNum: int,
    path: str,
):
    """Plot single slice of local topology map overlaid on HRCT.

    Args:
        ctArray (np.array): HRCT image to overlay PRGM RGB image onto
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity
        topoArray (np.arry): 3D map of local topology
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        path (str): path to save final image to
    """

    # take single slice of topoArray and mask out regions outside of maskArray
    topoArraySlice = topoArray[:, sliceNum, :]
    maskArraySlice = maskArray[:, sliceNum, :]
    topoArraySlice = np.ma.masked_where(maskArraySlice <= 0, topoArraySlice)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
    plt.imshow(ctArray[:, sliceNum, :], cmap="gray")
    plt.imshow(topoArraySlice, cmap="plasma")
    plt.axis("off")

    # add colorbar
    imgRatio = topoArraySlice.shape[0] / topoArraySlice.shape[1]
    cbar = plt.colorbar(fraction=0.047 * imgRatio, pad=0)

    # save figure and close it
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=600)
    plt.clf()
    plt.close()
