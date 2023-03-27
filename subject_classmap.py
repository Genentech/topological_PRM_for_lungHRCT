"""Class for generating HRCT PRM and topoligcal maps of the lungs."""
from os.path import join

import nibabel as nib
import numpy as np

import constants
from utils import img_utils, io_utils, plot_utils


class Subject(object):
    """Class for generating HRCT PRM and topological maps.

    Attributes:
        config (ConfigParser): configuration file
        subjID (str): subject ID
        outDir (str): main directory to save output files to
        expArray (np.array): expiratory HRCT in HU
        inspRegArray (np.array): inspiratory hrct in HU registered to expiratory HRCT
        maskArray (np.array): segmentation of thoracic cavity
        pixDims (np.array): pixel dimenstions in mm
        expArrayFilt (np.array): expiratory image with median filter applied
        inspRegArrayFilt (np.array): inspiratory image with median filter applied
        self.normArray (np.array): image denoting normal regions from PRM
        self.fSadArray (np.array): image denoting regions of non-emphysematous air trapping from PRM
        self.emphArray (np.array): image denoting emphysema regions from PRM
        self.emptEmphArray (np.array): image denoting regions of "emptying emphysema" from PRM
        self.prmAllArray (np.array): image denoting all four voxel classifications from PRM
    """

    def __init__(self, config):
        """Initialize Subject object."""
        self.config = config
        self.subjID = config["subjInfo"]["subjID"]
        self.outDir = config["io"]["outDir"]
        self.expArray = np.array([])
        self.inspRegArray = np.array([])
        self.maskArray = np.array([])
        self.pixDims = np.array([])
        self.expArrayFilt = np.array([])
        self.inspRegArrayFilt = np.array([])
        self.normArray = np.array([])
        self.fSadArray = np.array([])
        self.emphArray = np.array([])
        self.emptEmphArray = np.array([])
        self.prmAllArray = np.array([])

    def readCtFiles(self):
        """Read in files and convert to np.array.

        Read in files containing expiratory and inspiratory HRCTs and mask (.nii).
        Generate array denoting thoracic cavity outline.
        """
        self.expArray, self.pixDims = io_utils.readFiles(self.config["io"]["inFileExp"])
        self.inspRegArray, _ = io_utils.readFiles(self.config["io"]["inFileInspReg"])
        self.maskArray, _ = io_utils.readFiles(self.config["io"]["inFileMask"])

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
        """Create maps of norm, fSAD, emph, and emptying emph voxels.

        Classify voxels into norm, fSAD, emph, or emptying emph based on
        expiratory and inspiratory images and HU thresholds.
        """

        # get indices of normal, fSAD, emph, and emptying emph regions within mask
        normIdx = np.argwhere(
            (self.expArrayFilt > constants.prmProcessing.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.prmProcessing.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        fSadIdx = np.argwhere(
            (self.expArrayFilt <= constants.prmProcessing.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.prmProcessing.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        emphIdx = np.argwhere(
            (self.expArrayFilt <= constants.prmProcessing.EXP_THRESH)
            & (self.inspRegArrayFilt <= constants.prmProcessing.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        emptEmphIdx = np.argwhere(
            (self.expArrayFilt > constants.prmProcessing.EXP_THRESH)
            & (self.inspRegArrayFilt < constants.prmProcessing.INSP_THRESH)
            & (self.maskArray >= 1)
        )

        # fill arrays with zeros
        self.normArray = np.zeros(self.expArrayFilt.shape)
        self.fSadArray = np.zeros(self.expArrayFilt.shape)
        self.emphArray = np.zeros(self.expArrayFilt.shape)
        self.emptEmphArray = np.zeros(self.expArrayFilt.shape)
        self.prmAllArray = np.zeros(self.expArrayFilt.shape)

        # for individual classification arrays, set regions of interest =1
        self.normArray[normIdx[:, 0], normIdx[:, 1], normIdx[:, 2]] = 1
        self.fSadArray[fSadIdx[:, 0], fSadIdx[:, 1], fSadIdx[:, 2]] = 1
        self.emphArray[emphIdx[:, 0], emphIdx[:, 1], emphIdx[:, 2]] = 1
        self.emptEmphArray[emptEmphIdx[:, 0], emptEmphIdx[:, 1], emptEmphIdx[:, 2]] = 1

        # for combined classification, set each region to different number (using IMBIO standard)
        self.prmAllArray[
            normIdx[:, 0], normIdx[:, 1], normIdx[:, 2]
        ] = constants.prmProcessing.CLASSIFICATION_NUM_NORM
        self.prmAllArray[
            fSadIdx[:, 0], fSadIdx[:, 1], fSadIdx[:, 2]
        ] = constants.prmProcessing.CLASSIFICATION_NUM_FSAD
        self.prmAllArray[
            emphIdx[:, 0], emphIdx[:, 1], emphIdx[:, 2]
        ] = constants.prmProcessing.CLASSIFICATION_NUM_EMPH
        self.prmAllArray[
            emptEmphIdx[:, 0], emptEmphIdx[:, 1], emptEmphIdx[:, 2]
        ] = constants.prmProcessing.CLASSIFICATION_NUM_EMPTEMPH

    def savePrmNiis(self):
        """Save PRM maps as niftis."""

        # create path + file names
        normArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_NORM + self.subjID,
        )
        fSadArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_FSAD + self.subjID,
        )
        emphArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_EMPH + self.subjID,
        )
        emptEmphArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_EMPTEMPH + self.subjID,
        )
        prmAllArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_ALL + self.subjID,
        )

        # save arrays as niftis
        io_utils.saveAsNii(self.normArray, normArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.fSadArray, fSadArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.emphArray, emphArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.emptEmphArray, emptEmphArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.prmAllArray, prmAllArrayOutPath, self.pixDims)

    def genPrmColor(self):
        """Generate RGB images of PRM maps."""

        # convert binned PRM array to RGB color array
        self.prmAllArrayColor = img_utils.bin2rgb(
            self.prmAllArray, constants.prmProcessing.PRM_BIN2RGB
        )

        # plot representative slice of RGB color array
        prmAllArrayColorOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_ALL + "color_" + self.subjID + ".png",
        )
        plot_utils.plotPrmRgbImage(self.prmAllArrayColor, prmAllArrayColorOutPath)
