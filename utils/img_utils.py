"""Utils for image processing."""
from typing import Dict

import numpy as np
import scipy


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
