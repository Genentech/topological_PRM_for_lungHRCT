"""Utils for image processing."""
from typing import Dict

import numpy as np
import scipy
from quantimpy import minkowski as mk


def medFilt(image: np.ndarray, kernelSize: int):
    """Apply moving 2D median filter.

    Args:
        image (np.array): input image
        kernelSize (int): size of moving window

    Returns:
        imageFilt (np.array): filtered image
    """
    imageFilt = np.zeros(image.shape)
    for s in range(0, image.shape[0]):
        imageFilt[s, :, :] = scipy.signal.medfilt2d(
            image[s, :, :], kernel_size=kernelSize
        )

    return imageFilt


def bin2rgb(binImage: np.ndarray, colMap: Dict[int, np.ndarray]):
    """Convert image of bin numbers to rgb.

    Args:
        binImage (np.array): input image where each element denotes a bin number
        colMap (Dict[int, np.ndarray]): dictionary mapping bin number to RGB values

    Returns:
        rgbImage (np.array): RGB image of shape (x, y, z, 3)
    """
    rgbImage = np.zeros((binImage.shape[0], binImage.shape[1], binImage.shape[2], 3))
    for binNum in colMap.keys():
        rgbImage[binImage == binNum] = colMap[binNum]

    return rgbImage


def calcMkFns(binaryImage: np.ndarray, pixDims: np.ndarray):
    """Calculate Minkowski functionals for a binary image.

    Args:
        binaryImage (np.array): image of zeros and ones denoting PRM regions
        pixDims (np.array): pixel dimensions of binaryImage

    Returns:
        mkFnsArray (np.array): output array from quantimpy minkowski.functionals
    """

    # convert input binary array to boolean array
    boolImage = np.array(binaryImage, dtype=bool)

    # temporarily scale pix dims up (quantimpy requires them to be >=1)
    pixDimsScale = 10
    pixDimsRescale = pixDims * pixDimsScale

    # calculate global Mk fns using quantimpy and rescale to match original pixel dim units
    mkFnsArray = np.array(mk.functionals(boolImage, pixDimsRescale))
    mkFnsArray[0] /= pixDimsScale**3  # volume
    mkFnsArray[1] /= pixDimsScale**2  # surface area
    mkFnsArray[2] /= pixDimsScale  # curvature

    return mkFnsArray


def createMkFnDict(mkFnsArray: np.ndarray, label: str):
    """Create a dictionary containing Minkowski functionals.

    Args:
        mkFnsArray (np.array): output array from quantimpy minkowski.functionals
        label (str): string to append to dictionary entries denoting image type
                        (e.g. emph, fSAD, norm, emptemph)

    Returns:
        mkFnsDict (dict):
    """

    # store in mkFnsArray in dictionary with label label
    mkFnsDict = {}
    mkFnsDict["vol_" + label] = mkFnsArray[0]
    mkFnsDict["surf_area_" + label] = mkFnsArray[1]
    mkFnsDict["curv_" + label] = mkFnsArray[2]
    mkFnsDict["euler_" + label] = mkFnsArray[3]

    return mkFnsDict
