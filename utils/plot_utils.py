"""Utils for plotting"""
import pdb

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from PIL import Image

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


def cropImg(img: np.ndarray, mask: np.ndarray, pad: int):
    """Crop image based on mask boundaries.

    Args:
        img (np.array): nxm image
        mask (np.array): nxm mask specifying ROI in image
        pad (int): number of voxels to leave outside of the mask

    Returns:
        imgCrop (np.array): cropped image
        maskCrop (np.array): cropped mask
    """

    # get max and min of ROI along each dimension
    xMin = min(np.argwhere(mask > 0)[:, 1]) - pad
    xMax = max(np.argwhere(mask > 0)[:, 1]) + pad
    yMin = min(np.argwhere(mask > 0)[:, 0]) - pad
    yMax = max(np.argwhere(mask > 0)[:, 0]) + pad

    # create PIL Image objects from img and mask arrays
    imgI = Image.fromarray(img.astype(float))
    maskI = Image.fromarray(mask.astype(float))
    imgCropI = imgI.crop((xMin, yMin, xMax + 1, yMax + 1))
    maskCropI = maskI.crop((xMin, yMin, xMax + 1, yMax + 1))

    # convert back to np.arrays
    imgCrop = np.asarray(imgCropI)
    maskCrop = np.asarray(maskCropI)

    return imgCrop, maskCrop


def plotPrmRgb(
    maskArray: np.ndarray,
    prmAllArray: np.ndarray,
    sliceNum: int,
    path: str,
):
    """Plot single slice of RGB color array.

    Args:
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity
        prmAllArray (np.array): PRM image sorted into bins denoting PRM voxel classification
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        path (str): path to save final image to
    """

    # create custom discrete colormap and bounds
    cmap = colors.ListedColormap(["green", "yellow", "red", "purple"])
    bounds = [1, 2, 3, 4, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # take single slice of HRCT, mask, and PRM map
    maskArraySlice = maskArray[:, sliceNum, :]
    prmAllArraySlice = prmAllArray[:, sliceNum, :]

    # crop HRCT and PRM images
    prmAllArraySlice, _ = cropImg(
        prmAllArraySlice, maskArraySlice, constants.proc.PLOT_PAD
    )

    # mask out un-binned regions
    prmAllArraySlice = np.ma.masked_where(prmAllArraySlice < 1, prmAllArraySlice)

    # orient slice
    prmAllArraySlice = np.rot90(prmAllArraySlice, k=1)
    prmAllArraySlice = np.flip(prmAllArraySlice, axis=1)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
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


def plotTopo(
    maskArray: np.ndarray,
    topoArray: np.ndarray,
    sliceNum: int,
    mapUnit: str,
    path: str,
):
    """Plot single slice of local topology map.

    Args:
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity
        topoArray (np.arry): 3D map of local topology
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        mapUnit (str): string specifying units of topology metric
        path (str): path to save final image to
    """

    # take single slice of HRCT, mask, and topology map
    topoArraySlice = topoArray[:, sliceNum, :]
    maskArraySlice = maskArray[:, sliceNum, :]

    # crop HRCT, mask, and topology map
    topoArraySlice, maskArraySliceK = cropImg(
        topoArraySlice, maskArraySlice, constants.proc.PLOT_PAD
    )

    # mask out regions outside of maskArray
    topoArraySlice = np.ma.masked_where(maskArraySliceK <= 0, topoArraySlice)

    # orient slice
    topoArraySlice = np.rot90(topoArraySlice, k=1)
    topoArraySlice = np.flip(topoArraySlice, axis=1)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
    plt.imshow(topoArraySlice, cmap="plasma")
    plt.axis("off")

    # add colorbar
    imgRatio = topoArraySlice.shape[0] / topoArraySlice.shape[1]
    cbar = plt.colorbar(fraction=0.047 * imgRatio, pad=0)
    cbar.set_label(label=mapUnit, size=12)

    # save figure and close it
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=600)
    plt.clf()
    plt.close()


def plotPrmRgbOnCt(
    ctArray: np.ndarray,
    maskArray: np.ndarray,
    prmAllArray: np.ndarray,
    sliceNum: int,
    path: str,
):
    """Plot single slice of RGB color array overlaid on HRCT.

    Args:
        ctArray (np.array): HRCT image to overlay PRGM RGB image onto
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity
        prmAllArray (np.array): PRM image sorted into bins denoting PRM voxel classification
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        path (str): path to save final image to

    NOTE: this function is currently unused
    """

    # create custom discrete colormap and bounds
    cmap = colors.ListedColormap(["green", "yellow", "red", "purple"])
    bounds = [1, 2, 3, 4, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # take single slice of HRCT, mask, and PRM map
    ctArraySlice = ctArray[:, sliceNum, :]
    maskArraySlice = maskArray[:, sliceNum, :]
    prmAllArraySlice = prmAllArray[:, sliceNum, :]

    # crop HRCT and PRM images
    prmAllArraySlice, _ = cropImg(
        prmAllArraySlice, maskArraySlice, constants.proc.PLOT_PAD
    )
    ctArraySlice, _ = cropImg(ctArraySlice, maskArraySlice, constants.proc.PLOT_PAD)

    # mask out un-binned regions
    prmAllArraySlice = np.ma.masked_where(prmAllArraySlice < 1, prmAllArraySlice)

    # orient slices
    prmAllArraySlice = np.rot90(prmAllArraySlice, k=1)
    prmAllArraySlice = np.flip(prmAllArraySlice, axis=1)
    ctArraySlice = np.rot90(ctArraySlice, k=1)
    ctArraySlice = np.flip(ctArraySlice, axis=1)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
    plt.imshow(ctArraySlice, cmap="gray")
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
    mapUnit: str,
    path: str,
):
    """Plot single slice of local topology map overlaid on HRCT.

    Args:
        ctArray (np.array): HRCT image to overlay PRGM RGB image onto
        maskArray (np.array): segmentation mask where positive integers denote thoracic cavity
        topoArray (np.arry): 3D map of local topology
        sliceNum (int): index of slice along the anterior to posterior dimension to plot
        mapUnit (str): string specifying units of topology metric
        path (str): path to save final image to

    NOTE: this function is currently unused
    """

    # take single slice of HRCT, mask, and topology map
    topoArraySlice = topoArray[:, sliceNum, :]
    maskArraySlice = maskArray[:, sliceNum, :]
    ctArraySlice = ctArray[:, sliceNum, :]

    # crop HRCT, mask, and topology map
    topoArraySlice, maskArraySliceK = cropImg(
        topoArraySlice, maskArraySlice, constants.proc.PLOT_PAD
    )
    ctArraySlice, _ = cropImg(ctArraySlice, maskArraySlice, constants.proc.PLOT_PAD)

    # mask out regions outside of maskArray
    topoArraySlice = np.ma.masked_where(maskArraySliceK <= 0, topoArraySlice)

    # orient slices
    topoArraySlice = np.rot90(topoArraySlice, k=1)
    topoArraySlice = np.flip(topoArraySlice, axis=1)
    ctArraySlice = np.rot90(ctArraySlice, k=1)
    ctArraySlice = np.flip(ctArraySlice, axis=1)

    # plot image overlaid on corresponding HRCT slice
    plt.subplots()
    plt.imshow(ctArraySlice, cmap="gray")
    plt.imshow(topoArraySlice, cmap="plasma")
    plt.axis("off")

    # add colorbar
    imgRatio = topoArraySlice.shape[0] / topoArraySlice.shape[1]
    cbar = plt.colorbar(fraction=0.047 * imgRatio, pad=0)
    cbar.set_label(label=mapUnit, size=12)

    # save figure and close it
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=600)
    plt.clf()
    plt.close()
