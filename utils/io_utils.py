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
        pixDims (np.array): pixel dimensions (in mm), if available

    If pixel dimensions available, extract them.
    Supported file formats: .nii
    """

    # initialize pixDims array
    pixDims = np.array([])

    # extract file name from path
    fName = path.split("/")[-1]

    # check if file is .nii and read accordingly
    if ".nii" in fName:
        niiImg = nib.load(path)
        pixDims = niiImg.header["pixdim"][1:4]
        outArray = np.array(niiImg.dataobj)
    else:
        logging.warning("Registered HRCT file format is unsupported, must be .nii")

    return outArray, pixDims


def saveAsNii(inArray: np.array, path: str, pixDims=None):
    """Save np.array as nifti to path.

    Args:
        inArray (np.array): array to be saved as nifti
        path (str): path to save nifti to
        pixDims (np.array): optional argument speficying image pixel dimensions

    If pixel dimensions provided, add those to nifti header
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
