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


def calcGlobalMkFns(binaryImage: np.ndarray, pixDims: np.ndarray, imgType: str):
    """Calculate global Minkowski functionals.

    Args:
        binaryImage (np.array): image of zeros and ones denoting PRM regions
        pixDims (np.array): pixel dimensions of binaryImage
        imgType (str): string to append to dictionary entries denoting image type
                        (e.g. emph, fSAD, norm, emptemph)

    Returns:
        globalMkFnsDict (dict): dictionary containing value of global Minkowski functionals
                            (volume, surface area, curvature, and the Euler characteristic)
    """

    # convert input binary array to boolean array
    boolImage = np.array(binaryImage, dtype=bool)

    # calculate global Minkowski functionals using quantimpy and store in dictionary
    globalMkFnsArray = mk.functionals(boolImage, pixDims)
    globalMkFnsDict = {}
    globalMkFnsDict["vol_global_" + imgType] = globalMkFnsArray[0]
    globalMkFnsDict["surf_area_global_" + imgType] = globalMkFnsArray[1]
    globalMkFnsDict["curv_global" + imgType] = globalMkFnsArray[2]
    globalMkFnsDict["euler_global_" + imgType] = globalMkFnsArray[3]

    return globalMkFnsDict
