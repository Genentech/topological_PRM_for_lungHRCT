"""Class for generating HRCT PRM and topoligcal maps of the lungs."""
import logging
import os
import pdb
from os.path import exists, join

import nibabel as nib
import numpy as np
from quantimpy import minkowski as mk

import constants
from utils import img_utils, io_utils, plot_utils


class Subject(object):
    """Class for generating HRCT PRM and topological maps.

    Attributes:
        config (ConfigParser): configuration file
        subjID (str): subject ID
        outDir (str): main directory to save output files to
        expArray (np.array): expiratory HRCT in HU
        expArrayPlotting (np.array): expiratory HRCT in HU to use for plotting
        inspRegArray (np.array): inspiratory hrct in HU registered to expiratory HRCT
        maskArray (np.array): segmentation of thoracic cavity
        pixDims (np.array): pixel dimenstions in mm
        plotSliceNum (int): slice number along anterior-posterior dimnension to use for plotting
        expArrayFilt (np.array): expiratory image with median filter applied
        inspRegArrayFilt (np.array): inspiratory image with median filter applied
        prmAllArray (np.array): image denoting all four voxel classifications from PRM
        prmStats (dict): key PRM metrics/statistics
        binArrayDict (dict): binary image array of each bin category in image
        topoMapsHiResDict (dict): normalized 3D topological maps for each bin category in image
        topologyStatsGlobal (dict): global topology (Minkowski functionals) metrics for all prm maps
        topologyStatsLocal (dict): mean of local topology (Minkowski functionals) metrics for all prm maps
    """

    def __init__(self, config):
        """Initialize Subject object."""
        self.config = config
        self.subjID = config["subjInfo"]["subjID"]
        self.outDir = config["io"]["outDir"]
        self.expArray = np.array([])
        self.expArrayPlotting = np.array([])
        self.inspRegArray = np.array([])
        self.maskArray = np.array([])
        self.pixDims = np.array([])
        self.plotSliceNum = int
        self.expArrayFilt = np.array([])
        self.inspRegArrayFilt = np.array([])
        self.prmAllArray = np.array([])
        self.prmStats = {}
        self.binArrayDict = {}
        self.topoMapsHiResDict = {}
        self.topologyStatsGlobal = {}
        self.topologyStatsLocal = {}

    def readCtFiles(self):
        """Read in files and convert to np.array.

        Read in files containing expiratory and inspiratory HRCTs and mask (.nii).
        Save a copy of expiratory HRCT for plotting.
        """
        self.expArray, self.pixDims = io_utils.readFiles(self.config["io"]["inFileExp"])
        self.inspRegArray, _ = io_utils.readFiles(self.config["io"]["inFileInspReg"])
        self.maskArray, _ = io_utils.readFiles(self.config["io"]["inFileMask"])

        # ensure that HRCT arrays are in usable data type
        self.expArray = self.expArray.astype(float)
        self.inspRegArray = self.inspRegArray.astype(float)

        # make separate copy of expiratory HRCT to use for plotting
        self.expArrayPlotting = np.copy(self.expArray)

    def dimOutsideVoxels(self):
        """Dim voxels outside of thoracic cavity.

        NOTE: currently unused.
        """
        self.expArray[self.maskArray == 0] = constants.preProc.DIM_OUTSIDE_VAL
        self.inspRegArray[self.maskArray == 0] = constants.preProc.DIM_OUTSIDE_VAL

    def orientImages(self):
        """Put HRCTs into proper orientation.

        Rotate 180 degrees and reflect images
        """
        self.expArray = np.rot90(np.swapaxes(self.expArray, 0, 2), 2)[:, ::-1, :]
        self.expArrayPlotting = np.rot90(np.swapaxes(self.expArrayPlotting, 0, 2), 2)[
            :, ::-1, :
        ]
        self.inspRegArray = np.rot90(np.swapaxes(self.inspRegArray, 0, 2), 2)[
            :, ::-1, :
        ]
        self.maskArray = np.rot90(np.swapaxes(self.maskArray, 0, 2), 2)[:, ::-1, :]

    def applyMedFilts(self):
        """Apply moving 2D median filter to images."""

        self.expArrayFilt = img_utils.medFilt(
            self.expArray, constants.preProc.MEDFILT_KERNEL_SIZE
        )
        self.inspRegArrayFilt = img_utils.medFilt(
            self.inspRegArray, constants.preProc.MEDFILT_KERNEL_SIZE
        )

    def excludeVoxels(self):
        """Exclude voxels above and below certain thresholds.

        Exclude voxels from mask that fall above UPPERTHRESH and below LOWERTHRESH to
        minimize the contribution of blood vessels and airways.
        """

        # exclude voxels in filtered expiratory image
        self.maskArray[
            self.expArrayFilt < constants.preProc.EXCLUDE_VOX_LOWERTHRESH
        ] = 0
        self.maskArray[
            self.expArrayFilt > constants.preProc.EXCLUDE_VOX_UPPERTHRESH
        ] = 0

        # exclude voxels in filtered inspiratory image
        self.maskArray[
            self.inspRegArrayFilt < constants.preProc.EXCLUDE_VOX_LOWERTHRESH
        ] = 0
        self.maskArray[
            self.inspRegArrayFilt > constants.preProc.EXCLUDE_VOX_UPPERTHRESH
        ] = 0

    def classifyVoxelsPrm(self):
        """Create maps of norm, fSAD, emph, and emptying emph voxels.

        Classify voxels into norm, fSAD, emph, or emptying emph based on
        expiratory and inspiratory images and HU thresholds.
        """

        # get indices of normal, fSAD, emph, and emptying emph regions within mask
        normIdx = np.argwhere(
            (self.expArrayFilt > constants.proc.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.proc.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        fSadIdx = np.argwhere(
            (self.expArrayFilt < constants.proc.EXP_THRESH)
            & (self.inspRegArrayFilt > constants.proc.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        emphIdx = np.argwhere(
            (self.expArrayFilt < constants.proc.EXP_THRESH)
            & (self.inspRegArrayFilt < constants.proc.INSP_THRESH)
            & (self.maskArray >= 1)
        )
        emptEmphIdx = np.argwhere(
            (self.expArrayFilt > constants.proc.EXP_THRESH)
            & (self.inspRegArrayFilt <= constants.proc.INSP_THRESH)
            & (self.maskArray >= 1)
        )

        # for combined classification, set each region to different number (using IMBIO standard)
        self.prmAllArray = np.zeros(self.expArrayFilt.shape)
        self.prmAllArray[
            normIdx[:, 0], normIdx[:, 1], normIdx[:, 2]
        ] = constants.proc.PRM_NUM_NORM
        self.prmAllArray[
            fSadIdx[:, 0], fSadIdx[:, 1], fSadIdx[:, 2]
        ] = constants.proc.PRM_NUM_FSAD
        self.prmAllArray[
            emphIdx[:, 0], emphIdx[:, 1], emphIdx[:, 2]
        ] = constants.proc.PRM_NUM_EMPH
        self.prmAllArray[
            emptEmphIdx[:, 0], emptEmphIdx[:, 1], emptEmphIdx[:, 2]
        ] = constants.proc.PRM_NUM_EMPTEMPH

    def readPrmFile(self):
        """Read in PRM file and extract indvidual PRM maps."""
        self.prmAllArray, self.pixDims = io_utils.readFiles(
            self.config["io"]["inFilePrm"]
        )
        self.prmAllArray = np.rot90(self.prmAllArray, axes=(0, 2))

    def genMaskFromPrm(self):
        """Generate binary mask from PRM map."""

        self.maskArray = (self.prmAllArray > 0).astype(int)

    def genDictOfImageArrays(self):
        """Separate binned PRM image into a list of arrays.

        Each array in the list is a binary image for one of the bins.
        """

        # loop over bin numbers and extract individual binary image arrays
        for binNum in constants.proc.BINS:
            binArray = np.ascontiguousarray((self.prmAllArray == binNum).astype(int))
            self.binArrayDict[binNum] = binArray

    def calcPrmStats(self):
        """Calculate percentage of voxels in each PRM classification."""

        # get number of voxels in thoracic cavity
        numMaskVoxels = (self.maskArray > 0).sum()

        # define subject ID into stats dictionary
        self.prmStats["sid"] = self.subjID

        # loop over each bin number and get percentage of voxels in each bin
        for binNum in constants.proc.BINS:
            self.prmStats[constants.proc.BIN_DICT[binNum] + "_prct"] = (
                100 * np.count_nonzero(self.prmAllArray == binNum) / numMaskVoxels
            )

        # create out path and save stats as csv
        prmStatsOutPath = join(
            self.outDir, constants.outFileNames.PRM_STATS + self.subjID + ".csv"
        )
        io_utils.saveStatsCsv(self.prmStats, prmStatsOutPath)

    def savePrmNiis(self):
        """Save PRM maps as niftis."""

        # if it doesn't already exists, create directory for PRM niftis
        if not exists(join(self.outDir, constants.outFileNames.PRM_DIR)):
            os.mkdir(join(self.outDir, constants.outFileNames.PRM_DIR))

        # create path + file name
        prmAllArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_ALL + self.subjID,
        )

        # save array as nifti
        io_utils.saveAsNii(self.prmAllArray, prmAllArrayOutPath, self.pixDims)

    def plotPrmColor(self):
        """Generate RGB images of PRM maps."""

        # plot representative slice of RGB color array
        self.plotSliceNum = plot_utils.findPlotSliceNum(
            self.maskArray
        )  # get optimal slice number to plot
        prmAllArrayColorOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_ALL + "color_" + self.subjID + ".png",
        )
        plot_utils.plotPrmRgb(
            self.maskArray,
            self.prmAllArray,
            self.plotSliceNum,
            prmAllArrayColorOutPath,
        )

    def calcTopologyGlobal(self):
        """Calculate global Minkowski functionals for all PRM maps.

        Calculates global volume, surface area, curvature, and the Euler characteristic
        for classification maps of all PRM regions.
        """

        # loop over each bin number
        for binNum in constants.proc.BINS:
            # retrieve bin image array
            binArray = self.binArrayDict[binNum]

            # get array containing global mk fns for each prm region
            binGlobal = img_utils.calcMkFnsNorm(binArray, self.maskArray, self.pixDims)

            # create dictionary from array
            binGlobalDict = io_utils.createMkFnDict(
                binGlobal, "global_" + constants.proc.BIN_DICT[binNum]
            )

            # merge with main global topology stats dictionary
            self.topologyStatsGlobal.update(binGlobalDict)

    def genLocalTopoMaps(self):
        """Generate 3D maps of local topology features for all PRM maps.

        Calculates maps of local volume, surface area, curvature, and the Euler characteristic
        for classification maps of all PRM regions.
        """

        # loop over each bin number
        for binNum in constants.proc.BINS:
            # generate low resolution 3D local topolgy maps
            binTopoMaps = img_utils.genLowResTopoMaps(
                self.binArrayDict[binNum], self.maskArray, self.pixDims
            )
            logging.info(
                constants.proc.BIN_DICT[binNum]
                + " low resolution local topology mapping complete."
            )

            # generate high resolution 3D maps using interpolation
            self.topoMapsHiResDict[binNum] = img_utils.resizeTopoMaps(
                self.binArrayDict[binNum].shape, self.maskArray, binTopoMaps
            )
            logging.info(
                constants.proc.BIN_DICT[binNum]
                + " high resolution local topology mapping interpolation complete."
            )

    def calcMeanLocalTopoStats(self):
        """Calculate whole lung mean of each topology metric in local topology maps."""

        # loop over each bin number
        for binNum in constants.proc.BINS:
            # calculate whole lung mean of each topology map
            binMeanLocal = img_utils.meanLocalTopo(
                self.topoMapsHiResDict[binNum], self.maskArray
            )
            # create individual dictionary
            binLocalDict = io_utils.createMkFnDict(
                binMeanLocal, "local_" + constants.proc.BIN_DICT[binNum]
            )
            # merge with main local topology stats dictionary
            self.topologyStatsLocal.update(binLocalDict)

    def saveTopoNiis(self):
        """Save topology maps as niftis."""

        # if it doesn't already exists, create directory for topology niftis
        if not exists(join(self.outDir, constants.outFileNames.TOPO_DIR)):
            os.mkdir(join(self.outDir, constants.outFileNames.TOPO_DIR))

        # loop over each bin number
        for binNum in constants.proc.BINS:
            # loop over topology arrays and save them as niftis
            for i in range(self.topoMapsHiResDict[binNum].shape[0]):
                outPath = join(
                    self.outDir,
                    constants.outFileNames.TOPO_DIR,
                    "prm_"
                    + constants.proc.BIN_DICT[binNum]
                    + "_"
                    + constants.outFileNames.TOPO[i]
                    + self.subjID,
                )
                io_utils.saveAsNii(
                    self.topoMapsHiResDict[binNum][i, :, :, :], outPath, self.pixDims
                )

    def plotTopoColor(self):
        """Plot slice of select local topology maps."""

        # specify units of topology map
        mapUnit = "m$^{-1}$"

        # plot surface area density for each bin number
        for binNum in constants.proc.BINS:
            # norm surface area density
            binAreaOutPath = join(
                self.outDir,
                "prm_"
                + constants.proc.BIN_DICT[binNum]
                + "_"
                + constants.outFileNames.TOPO[1]
                + "color_"
                + self.subjID
                + ".png",
            )
            plot_utils.plotTopo(
                self.maskArray,
                self.topoMapsHiResDict[binNum][1, :, :, :],
                self.plotSliceNum,
                mapUnit,
                binAreaOutPath,
            )

    def saveTopologyStats(self):
        """Save combined global and local topology metrics in CSV."""

        # create combined dictionary with global and local topology metrics
        topologyStats = {}
        topologyStats["sid"] = self.subjID
        topologyStats.update(self.topologyStatsGlobal)
        topologyStats.update(self.topologyStatsLocal)

        # save combined dictionary as CSV
        topologyStatsOutPath = join(
            self.outDir, constants.outFileNames.TOPO_STATS + self.subjID + ".csv"
        )
        io_utils.saveStatsCsv(topologyStats, topologyStatsOutPath)
