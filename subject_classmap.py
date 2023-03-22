"""Class for generating HRCT PRM and topoligcal maps of the lungs."""
import nibabel as nib
import numpy as np

import constants
from utils import img_utils, io_utils


class Subject(object):
    """Class for generating HRCT PRM and topological maps.

    Attributes:
        config (ConfigParser):
        expArray (np.array): expiratory HRCT in HU
        inspRegArray (np.array): inspiratory hrct in HU registered to expiratory HRCT
        maskArray (np.array): segmentation of thoracic cavity
        expArrayFilt (np.array): expiratory image with median filter applied
        inspRegArrayFilt (np.array): inspiratory image with median filter applied
    """

    def __init__(self, config):
        """Initialize Subject object."""
        self.config = config
        self.expArray = np.array([])
        self.inspRegArray = np.array([])
        self.maskArray = np.array([])
        self.expArrayFilt = np.array([])
        self.inspRegArrayFilt = np.array([])

    def readCtFiles(self):
        """Read in files and convert to np.array.

        Read in files containing expiratory and inspiratory HRCTs and mask (.nii).
        """
        self.expArray = io_utils.readFiles(self.config["inFiles"]["exp"])
        self.inspRegArray = io_utils.readFiles(self.config["inFiles"]["insp_reg"])
        self.maskArray = io_utils.readFiles(self.config["inFiles"]["mask"])

    def dimOutsideVoxels(self):
        """Dim voxels outside of thoracic cavity."""
        self.expArray[self.maskArray == 0] = constants.prePrmProcessing.DIM_OUTSIDE_VAL
        self.inspRegArray[
            self.maskArray == 0
        ] = constants.prePrmProcessing.DIM_OUTSIDE_VAL

    def orientImages(self):
        """Put HRCTs into proper orientation.

        Rotate 180 degrees and reflect images
        """
        self.expArray = np.rot90(np.swapaxes(self.expArray, 0, 2), 2)[:, ::-1, :]
        self.inspRegArray = np.rot90(np.swapaxes(self.inspRegArray, 0, 2), 2)[
            :, ::-1, :
        ]
        self.maskArray = np.rot90(np.swapaxes(self.maskArray, 0, 2), 2)[:, ::-1, :]

    def applyMedFilts(self):
        """Apply moving 2D median filter to images."""

        self.expArrayFilt = img_utils.medFilt(
            self.expArray, constants.prePrmProcessing.MEDFILT_KERNEL_SIZE
        )
        self.inspRegArrayFilt = img_utils.medFilt(
            self.inspRegArray, constants.prePrmProcessing.MEDFILT_KERNEL_SIZE
        )
