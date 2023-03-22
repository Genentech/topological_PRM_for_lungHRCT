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
        self.normArray (np.array): image denoting normal regions from PRM
        self.emphArray (np.array): image denoting emphysema regions from PRM
        self.fSadArray (np.array): image denoting regions of non-emphysematous air trapping from PRM
    """

    def __init__(self, config):
        """Initialize Subject object."""
        self.config = config
        self.expArray = np.array([])
        self.inspRegArray = np.array([])
        self.maskArray = np.array([])
        self.expArrayFilt = np.array([])
        self.inspRegArrayFilt = np.array([])
        self.normArray = np.array([])
        self.emphArray = np.array([])
        self.fSadArray = np.array([])

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

    def excludeVoxels(self):
        """Exclude voxels above and below certain thresholds.

        Exclude voxels from mask that fall above upperThresh to
        minimize the contribution of blood vessels and airways.
        """

        # exclude voxels in filtered expiratory image
        self.maskArray[
            self.expArrayFilt < constants.prePrmProcessing.EXCLUDE_VOX_LOWERTHRESH
        ] = 0
        self.maskArray[
            self.expArrayFilt > constants.prePrmProcessing.EXCLUDE_VOX_UPPERTHRESH
        ] = 0

        # exclude voxels in filtered inspiratory image
        self.maskArray[
            self.inspRegArrayFilt < constants.prePrmProcessing.EXCLUDE_VOX_LOWERTHRESH
        ] = 0
        self.maskArray[
            self.inspRegArrayFilt > constants.prePrmProcessing.EXCLUDE_VOX_UPPERTHRESH
        ] = 0

    def classifyVoxelsPrm(self):
        """Create maps of emph, fSAD, and norm voxeles.

        Classify voxels into emph, fSAD, or normal based on
        expiratory and inspiratory images and HU thresholds.
        """

        # get indices of normal, emph, and fSAD regions within mask
        normIdx = np.argwhere(
            (self.expArrayFilt > constants.prmThresholds.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.prmThresholds.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        emphIdx = np.argwhere(
            (self.expArrayFilt < constants.prmThresholds.EXP_THRESH)
            & (self.inspRegArrayFilt < constants.prmThresholds.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        fSadIdx = np.argwhere(
            (self.expArrayFilt < constants.prmThresholds.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.prmThresholds.INSP_THRESH)
            & (self.maskArray >= 1)
        )

        # fill arrays with zeros
        self.normArray = np.zeros(self.expArrayFilt.shape)
        self.emphArray = np.zeros(self.expArrayFilt.shape)
        self.fSadArray = np.zeros(self.expArrayFilt.shape)

        # set regions of interest =1
        self.normArray[normIdx[:, 0], normIdx[:, 1], normIdx[:, 2]] = 1
        self.emphArray[emphIdx[:, 0], emphIdx[:, 1], emphIdx[:, 2]] = 1
        self.fSadArray[fSadIdx[:, 0], fSadIdx[:, 1], fSadIdx[:, 2]] = 1
