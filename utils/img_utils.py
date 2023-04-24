"""Utils for image processing."""
import math
from typing import Dict

import numpy as np
import scipy.signal
from quantimpy import minkowski as mk

import constants


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


def genLowResGrid(highResImgShape: np.ndarray, windowRadius: int, gridRes: int):
    """Create low resolution grid for output image after passing moving window
    over high resolution image every jth voxel (j=gridRes).

    Args:
        highResImgShape (np.array): nxnxn shape of high res image that moving window will be passed over
        windowRadius (int): half the length of one side of nxnxn moving window
        gridRes (int): interval (in voxels) between nxnxn moving windows

    Returns:
        lowResGrid (np.ndarray): low resolution array of zeros for storing output after passing moving window
                                    over high resolution image every jth voxel
    """
    x = math.ceil((highResImgShape[0] - windowRadius * 2) / gridRes) + 1
    y = math.ceil((highResImgShape[1] - windowRadius * 2) / gridRes) + 1
    z = math.ceil((highResImgShape[2] - windowRadius * 2) / gridRes) + 1
    lowResGrid = np.zeros(x, y, z)

    return lowResGrid


def genLocalTopoMaps(binaryImage: np.ndarray, pixDims: np.ndarray):
    """Generate 3D maps of local topology for a binary image.

    Args:
        binaryImage (np.array): image of zeros and ones denoting PRM regions
        pixDims (np.array): pixel dimensions of binaryImage

    Returns:
        volMap (np.array): 3D map of local volume
        areaMap (np.array): 3D map of local surface area
        curvMap (np.array): 3D map of local curvature
        eulerMap (np.array): 3D map of local Euler-Poincare characteristic

    Note: output maps are lower resolution than input binaryImage.
    """

    # initialize low resolution arrays in which to store local topology calculations
    volMap = genLowResGrid(
        binaryImage.shape,
        constants.topoMapping.WIND_RADIUS,
        constants.topoMapping.GRID_RES,
    )
    areaMap = genLowResGrid(
        binaryImage.shape,
        constants.topoMapping.WIND_RADIUS,
        constants.topoMapping.GRID_RES,
    )
    curvMap = genLowResGrid(
        binaryImage.shape,
        constants.topoMapping.WIND_RADIUS,
        constants.topoMapping.GRID_RES,
    )
    eulerMap = genLowResGrid(
        binaryImage.shape,
        constants.topoMapping.WIND_RADIUS,
        constants.topoMapping.GRID_RES,
    )

    # pass moving window over binaryImage every jth voxel and calculate local topology for each window

    # get indices of the center of each window along each dimension in high resolution input binary image
    iIdxHighRes = range(
        constants.topoMapping.WIND_RADIUS,
        binaryImage.shape[0] - constants.topoMapping.WIND_RADIUS + 1,
        constants.topoMapping.GRID_RES,
    )
    jIdxHighRes = range(
        constants.topoMapping.WIND_RADIUS,
        binaryImage.shape[1] - constants.topoMapping.WIND_RADIUS + 1,
        constants.topoMapping.GRID_RES,
    )
    kIdxHighRes = range(
        constants.topoMapping.WIND_RADIUS,
        binaryImage.shape[2] - constants.topoMapping.WIND_RADIUS + 1,
        constants.topoMapping.GRID_RES,
    )

    # initialize indices for storing local topology metrics in low resolution output maps
    iIdxLowRes = 0
    jIdxLowRes = 0
    kIdxLowRes = 0

    for i in iIdxHighRes:
        for j in jIdxHighRes:
            for k in kIdxHighRes:
                # create local binary from 3D window of input binary image
                localBinaryImage = binaryImage[
                    (i - constants.topoMapping.WIND_RADIUS) : (
                        i + constants.topoMapping.WIND_RADIUS
                    ),
                    (j - constants.topoMapping.WIND_RADIUS) : (
                        j + constants.topoMapping.WIND_RADIUS
                    ),
                    (k - constants.topoMapping.WIND_RADIUS) : (
                        k + constants.topoMapping.WIND_RADIUS
                    ),
                ]

                # calculate Minkowski functionals for local binary image
                mkFnsArray = calcMkFns(localBinaryImage, pixDims)

                # store Minkowski functionals
                volMap[iIdxLowRes, jIdxLowRes, kIdxLowRes] = mkFnsArray[0]
                areaMap[iIdxLowRes, jIdxLowRes, kIdxLowRes] = mkFnsArray[1]
                curvMap[iIdxLowRes, jIdxLowRes, kIdxLowRes] = mkFnsArray[2]
                eulerMap[iIdxLowRes, jIdxLowRes, kIdxLowRes] = mkFnsArray[3]
