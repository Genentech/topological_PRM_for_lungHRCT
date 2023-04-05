"""Utils for image processing."""
from typing import Dict

import numpy as np
import scipy
from quantimpy import minkowski as mk


def medFilt(image: np.array, kernelSize: int):
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


def bin2rgb(binImage: np.array, colMap: Dict[int, np.ndarray]):
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


def calcMinkowskiFns(binaryImage: np.array, pixDims: np.array):
    """Calculate global Minkowski functionals.

    Args:
        binaryImage (np.array): image of zeros and ones denoting PRM regions
        pixDims (np.array): pixel dimensions of binaryImage

    Returns:
        globalMkFns (dict): dictionary containing value of global Minkowski functionals
                            (volume, surface area, curvature, and the Euler characteristic)
    """

    # convert input binary array to boolean array
    boolImage = np.array(binaryImage, dtype=bool)

    # calculate global Minkowski functionals using quantimpy and store in dictionary
    globalMkFnsArray = mk.functionals(boolImage, pixDims)
    globalMkFnsArray["vol"] = globalMkFnsArray[0]
