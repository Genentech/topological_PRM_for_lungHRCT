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
        normArray (np.array): image denoting normal regions from PRM
        fSadArray (np.array): image denoting regions of non-emphysematous air trapping from PRM
        emphArray (np.array): image denoting emphysema regions from PRM
        emptEmphArray (np.array): image denoting regions of "emptying emphysema" from PRM
        prmAllArray (np.array): image denoting all four voxel classifications from PRM
        prmStats (dict): key PRM metrics/statistics
        topologyGlobal (dict): global topology (Minkowski functionals) metrics for all prm maps
        topologyStatsLocal (dict): mean of local topology (Minkowski functionals) metrics for all prm maps
        self.normTopoMapsHiRes (np.array): normalized 3D topological maps of normal PRM voxels
        self.fSadTopoMapsHiRes (np.array): normalized 3D topological maps of fSAD PRM voxels
        self.emphTopoMapsHiRes (np.array): normalized 3D topological maps of emphysema PRM voxels
        self.emptEmphTopoMapsHiRes (np.array): normalized 3D topological maps of emptying emphysema PRM voxels
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
        self.normArray = np.array([])
        self.fSadArray = np.array([])
        self.emphArray = np.array([])
        self.emptEmphArray = np.array([])
        self.prmAllArray = np.array([])
        self.prmStats = {}
        self.topologyStatsGlobal = {}
        self.topologyStatsLocal = {}
        self.normTopoMapsHiRes = np.array([])
        self.fSadTopoMapsHiRes = np.array([])
        self.emphTopoMapsHiRes = np.array([])
        self.emptEmphTopoMapsHiRes = np.array([])

    def readCtFiles(self):
        """Read in files and convert to np.array.

        Read in files containing expiratory and inspiratory HRCTs and mask (.nii).
        Save a copy of expiratory HRCT for plotting.
        """
        self.expArray, self.pixDims = io_utils.readFiles(self.config["io"]["inFileExp"])
        self.inspRegArray, _ = io_utils.readFiles(self.config["io"]["inFileInspReg"])
        self.maskArray, _ = io_utils.readFiles(self.config["io"]["inFileMask"])

        # make separate copy of expiratory HRCT to use for plotting
        self.expArrayPlotting = np.copy(self.expArray)

    def dimOutsideVoxels(self):
        """Dim voxels outside of thoracic cavity."""
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
        ] = constants.proc.CLASSIFICATION_NUM_NORM
        self.prmAllArray[
            fSadIdx[:, 0], fSadIdx[:, 1], fSadIdx[:, 2]
        ] = constants.proc.CLASSIFICATION_NUM_FSAD
        self.prmAllArray[
            emphIdx[:, 0], emphIdx[:, 1], emphIdx[:, 2]
        ] = constants.proc.CLASSIFICATION_NUM_EMPH
        self.prmAllArray[
            emptEmphIdx[:, 0], emptEmphIdx[:, 1], emptEmphIdx[:, 2]
        ] = constants.proc.CLASSIFICATION_NUM_EMPTEMPH

    def calcPrmStats(self):
        """Calculate key PRM statistics."""

        # calculate percentage of voxels in each PRM classification
        numMaskVoxels = (self.maskArray > 0).sum()
        self.prmStats["sid"] = self.subjID
        self.prmStats["PRM_norm_prct"] = 100 * np.sum(self.normArray) / numMaskVoxels
        self.prmStats["PRM_fSAD_prct"] = 100 * np.sum(self.fSadArray) / numMaskVoxels
        self.prmStats["PRM_emph_prct"] = 100 * np.sum(self.emphArray) / numMaskVoxels
        self.prmStats["PRM_emptemph_prct"] = (
            100 * np.sum(self.emptEmphArray) / numMaskVoxels
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

        # create path + file names
        normArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_NORM + self.subjID,
        )
        fSadArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_FSAD + self.subjID,
        )
        emphArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_EMPH + self.subjID,
        )
        emptEmphArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_EMPTEMPH + self.subjID,
        )
        prmAllArrayOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_DIR,
            constants.outFileNames.PRM_ALL + self.subjID,
        )

        # save arrays as niftis
        io_utils.saveAsNii(self.normArray, normArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.fSadArray, fSadArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.emphArray, emphArrayOutPath, self.pixDims)
        io_utils.saveAsNii(self.emptEmphArray, emptEmphArrayOutPath, self.pixDims)
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
        plot_utils.plotPrmRgbOnCt(
            self.expArrayPlotting,
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

        # get individual arrays containing global mk fns for each prm region
        normGlobal = img_utils.calcMkFnsNorm(
            self.normArray, self.maskArray, self.pixDims
        )
        fSadGlobal = img_utils.calcMkFnsNorm(
            self.fSadArray, self.maskArray, self.pixDims
        )
        emphGlobal = img_utils.calcMkFnsNorm(
            self.emphArray, self.maskArray, self.pixDims
        )
        emptEmphGlobal = img_utils.calcMkFnsNorm(
            self.emptEmphArray, self.maskArray, self.pixDims
        )

        # create individual dictionaries from arrays
        normGlobalDict = io_utils.createMkFnDict(normGlobal, "global_norm")
        fSadGlobalDict = io_utils.createMkFnDict(fSadGlobal, "global_fSAD")
        emphGlobalDict = io_utils.createMkFnDict(emphGlobal, "global_emph")
        emptEmphGlobalDict = io_utils.createMkFnDict(emptEmphGlobal, "global_emptemph")

        # merge indivudal dictionaries
        self.topologyStatsGlobal.update(normGlobalDict)
        self.topologyStatsGlobal.update(fSadGlobalDict)
        self.topologyStatsGlobal.update(emphGlobalDict)
        self.topologyStatsGlobal.update(emptEmphGlobalDict)

    def genLocalTopoMaps(self):
        """Generate 3D maps of local topology features for all PRM maps.

        Calculates maps of local volume, surface area, curvature, and the Euler characteristic
        for classification maps of all PRM regions.
        """

        # generate low resolution 3D local topolgy maps
        normTopoMaps = img_utils.genLowResTopoMaps(
            self.normArray, self.maskArray, self.pixDims
        )
        logging.info("Normal low resolution local topology mapping complete.")
        fSadTopoMaps = img_utils.genLowResTopoMaps(
            self.fSadArray, self.maskArray, self.pixDims
        )
        logging.info("fSAD low resolution local topology mapping complete.")
        emphTopoMaps = img_utils.genLowResTopoMaps(
            self.emphArray, self.maskArray, self.pixDims
        )
        logging.info("Emphysema low resolution local topology mapping complete.")
        emptEmphTopoMaps = img_utils.genLowResTopoMaps(
            self.emptEmphArray, self.maskArray, self.pixDims
        )
        logging.info(
            "Emptying emphysema low resolution local topology mapping complete."
        )

        # generate high resolution 3D maps using interpolation
        self.normTopoMapsHiRes = img_utils.resizeTopoMaps(
            self.normArray.shape, self.maskArray, normTopoMaps
        )
        self.fSadTopoMapsHiRes = img_utils.resizeTopoMaps(
            self.fSadArray.shape, self.maskArray, fSadTopoMaps
        )
        self.emphTopoMapsHiRes = img_utils.resizeTopoMaps(
            self.emphArray.shape, self.maskArray, emphTopoMaps
        )
        self.emptEmphTopoMapsHiRes = img_utils.resizeTopoMaps(
            self.emptEmphArray.shape, self.maskArray, emptEmphTopoMaps
        )

    def calcMeanLocalTopoStats(self):
        """Calculate whole lung mean of each topology metric in local topology maps."""

        # calculate whole lung mean of each topology map
        normMeanLocal = img_utils.meanLocalTopo(self.normTopoMapsHiRes, self.maskArray)
        fSadMeanLocal = img_utils.meanLocalTopo(self.fSadTopoMapsHiRes, self.maskArray)
        emphMeanLocal = img_utils.meanLocalTopo(self.emphTopoMapsHiRes, self.maskArray)
        emptEmphMeanLocal = img_utils.meanLocalTopo(
            self.emptEmphTopoMapsHiRes, self.maskArray
        )

        # create individual dictionaries
        normLocalDict = io_utils.createMkFnDict(normMeanLocal, "local_norm")
        fSadLocalDict = io_utils.createMkFnDict(fSadMeanLocal, "local_fSAD")
        emphLocalDict = io_utils.createMkFnDict(emphMeanLocal, "local_emph")
        emptEmphLocalDict = io_utils.createMkFnDict(emptEmphMeanLocal, "local_emptemph")

        # merge indivudal dictionaries
        self.topologyStatsLocal.update(normLocalDict)
        self.topologyStatsLocal.update(fSadLocalDict)
        self.topologyStatsLocal.update(emphLocalDict)
        self.topologyStatsLocal.update(emptEmphLocalDict)

    def saveTopoNiis(self):
        """Save topology maps as niftis."""

        # if it doesn't already exists, create directory for topology niftis
        if not exists(join(self.outDir, constants.outFileNames.TOPO_DIR)):
            os.mkdir(join(self.outDir, constants.outFileNames.TOPO_DIR))

        # loop over normal topology arrays and save them as niftis
        for i in range(self.normTopoMapsHiRes.shape[0]):
            outPath = join(
                self.outDir,
                constants.outFileNames.TOPO_DIR,
                constants.outFileNames.PRM_NORM
                + constants.outFileNames.TOPO[i]
                + self.subjID,
            )
            io_utils.saveAsNii(
                self.normTopoMapsHiRes[i, :, :, :], outPath, self.pixDims
            )

        # loop over fSAD topology arrays and save them as niftis
        for i in range(self.fSadTopoMapsHiRes.shape[0]):
            outPath = join(
                self.outDir,
                constants.outFileNames.TOPO_DIR,
                constants.outFileNames.PRM_FSAD
                + constants.outFileNames.TOPO[i]
                + self.subjID,
            )
            io_utils.saveAsNii(
                self.fSadTopoMapsHiRes[i, :, :, :], outPath, self.pixDims
            )

        # loop over emphysema topology arrays and save them as niftis
        for i in range(self.emphTopoMapsHiRes.shape[0]):
            outPath = join(
                self.outDir,
                constants.outFileNames.TOPO_DIR,
                constants.outFileNames.PRM_EMPH
                + constants.outFileNames.TOPO[i]
                + self.subjID,
            )
            io_utils.saveAsNii(
                self.emphTopoMapsHiRes[i, :, :, :], outPath, self.pixDims
            )

        # loop over emptying emphysema topology arrays and save them as niftis
        for i in range(self.emptEmphTopoMapsHiRes.shape[0]):
            outPath = join(
                self.outDir,
                constants.outFileNames.TOPO_DIR,
                constants.outFileNames.PRM_EMPTEMPH
                + constants.outFileNames.TOPO[i]
                + self.subjID,
            )
            io_utils.saveAsNii(
                self.emptEmphTopoMapsHiRes[i, :, :, :], outPath, self.pixDims
            )

    def plotTopoColor(self):
        """Plot slice of select local topology maps."""

        # specify topology map type
        mapType = "surface area density"

        # norm surface area density
        normAreaOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_NORM
            + constants.outFileNames.TOPO[1]
            + "color_"
            + self.subjID
            + ".png",
        )
        plot_utils.plotTopoOnCt(
            self.expArrayPlotting,
            self.maskArray,
            self.normTopoMapsHiRes[1, :, :, :],
            self.plotSliceNum,
            mapType,
            normAreaOutPath,
        )

        # fSAD surface area density
        fSadAreaOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_FSAD
            + constants.outFileNames.TOPO[1]
            + "color_"
            + self.subjID
            + ".png",
        )
        plot_utils.plotTopoOnCt(
            self.expArrayPlotting,
            self.maskArray,
            self.fSadTopoMapsHiRes[1, :, :, :],
            self.plotSliceNum,
            mapType,
            fSadAreaOutPath,
        )

        # emph surface area density
        emphAreaOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_EMPH
            + constants.outFileNames.TOPO[1]
            + "color_"
            + self.subjID
            + ".png",
        )
        plot_utils.plotTopoOnCt(
            self.expArrayPlotting,
            self.maskArray,
            self.emphTopoMapsHiRes[1, :, :, :],
            self.plotSliceNum,
            mapType,
            emphAreaOutPath,
        )

        # emptying emphysema surface area density
        emptEmphAreaOutPath = join(
            self.outDir,
            constants.outFileNames.PRM_EMPTEMPH
            + constants.outFileNames.TOPO[1]
            + "color_"
            + self.subjID
            + ".png",
        )
        plot_utils.plotTopoOnCt(
            self.expArrayPlotting,
            self.maskArray,
            self.emptEmphTopoMapsHiRes[1, :, :, :],
            self.plotSliceNum,
            mapType,
            emptEmphAreaOutPath,
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
