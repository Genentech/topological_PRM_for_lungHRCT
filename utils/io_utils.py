"""Import and export util functions."""
import logging

import nibabel as nib
import numpy as np
import scipy


def readFiles(path: str):
    """Read file from path as np.array.

    Args:
        path (str): path of file

    Returns:
        outArray (np.array): image contents of file

    Supported file formats: .nii
    """

    # extract file name from path
    fName = path.split("/")[-1]

    # check if file is .nii and read accordingly
    if ".nii" in fName:
        outArray = np.array(nib.load(path).dataobj)
    else:
        logging.warning("Registered HRCT file format is unsupported, must be .nii")

    return outArray
