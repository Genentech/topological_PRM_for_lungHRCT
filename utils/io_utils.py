"""Import and export util functions."""
import logging
from typing import Dict

import nibabel as nib
import numpy as np
import pandas as pd


def readFiles(path: str):
    """Read in image as np.array.

    Args:
        path (str): path of file

    Returns:
        outArray (np.array): image contents of file
        pixDims (np.array): pixel dimensions (in mm)

    If pixel dimensions available, extract them.
    Supported file formats: .nii
    """

    # initialize pixDims array
    pixDims = np.array([])

    # extract file name from path
    fName = path.split("/")[-1]

    # check if file is .nii and read accordingly
    if ".nii" in fName:
        # load nifti
        niiImg = nib.load(path)

        # convert to RAS orientation
        niiImg = nib.as_closest_canonical(niiImg)

        # extract pixel dimensions in mm
        pixDims = niiImg.header["pixdim"][1:4]

        # convert to numpy array
        outArray = np.array(niiImg.dataobj)
    else:
        logging.warning("Registered HRCT file format is unsupported, must be .nii")

    return outArray, pixDims


def saveAsNii(inArray: np.ndarray, path: str, pixDims=None):
    """Save np.array as nifti to path.

    Args:
        inArray (np.array): array to be saved as nifti
        path (str): path to save nifti to
        pixDims (np.array): optional argument speficying image pixel dimensions

    If pixel dimensions provided, add those to nifti header.
    """
    # if pixel dimensions provided, add them to nifti header
    header = nib.Nifti1Header()
    if pixDims is not None:
        header["pixdim"][1:4] = pixDims

    # create nifti image and save it
    outNiiImg = nib.Nifti1Image(inArray, np.eye(4), header)
    nib.save(outNiiImg, path)


def saveStatsCsv(statsDict: Dict[str, float], path: str):
    """Save statistics dictionary as csv.

    Args:
        statsDict (dict): dictionary containing stats names and values, first column must be "sid"
        path (str): path to save csv to
    """
    # convert dict to pandas DataFrame, set index to sid, and save as csv
    statsDf = pd.DataFrame.from_dict([statsDict])
    statsDf.set_index("sid", inplace=True)
    statsDf.to_csv(path)


def createMkFnDict(mkFnsArray: np.ndarray, label: str):
    """Create a dictionary containing Minkowski functionals.

    Args:
        mkFnsArray (np.array): output array from quantimpy minkowski.functionals
        label (str): string to append to dictionary entries denoting image type
                        (e.g. emph, fSAD, norm, emptemph)

    Returns:
        mkFnsDict (dict): dictionary of Minkowski functionals
    """

    # store in mkFnsArray in dictionary with label label
    mkFnsDict = {}
    mkFnsDict["vol_" + label] = mkFnsArray[0]
    mkFnsDict["surf_area_" + label] = mkFnsArray[1]
    mkFnsDict["curv_" + label] = mkFnsArray[2]
    mkFnsDict["euler_" + label] = mkFnsArray[3]

    return mkFnsDict
