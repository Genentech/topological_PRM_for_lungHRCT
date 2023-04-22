"""Utils for image processing."""
from typing import Dict

import numpy as np
import scipy.signal
from quantimpy import minkowski as mk


def normalizeCt(image: np.ndarray, mask: np.ndarray):
    """Normalize CT image.

    Args:
        image (np.array): CT image to be normalized
        mask (np.array): segmentation mask of CT image, must have same shape as image

    Returns:
        imageNorm (np.array): normalized CT image

    Sets positive voxels to zero and normalizes CT based on min and max voxels in mask.
    """

    # set positive voxels to zero and normalize image
    image[image > 0] = 0
    minVox = min(image[mask >= 1].flatten())
    maxVox = max(image[mask >= 1].flatten())
    imageNorm = 1000 / (maxVox - minVox) * (image - minVox) - 1000

    return imageNorm


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

    NOTE: this function is currently unused in the pipeline
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
