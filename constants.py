"""Define important constants used throughout the pipeline."""


class prePrmProcessing(object):
    """Constants for preprocessing before PRM, after registration.

    Attributes:
        DIM_OUTSIDE_VAL (float): HU value to set voxels outside of mask to
        MEDFILT_KERNEL_SIZE (int): size of moving window for median filtering of images
        EXCLUDE_VOX_LOWERTHRESH (float): HU value below which to exclude voxels from mask
        EXCLUDE_VOX_UPPERTHRESH (float): HU value above which to exclude voxels from mask
    """

    DIM_OUTSIDE_VAL = -2000
    MEDFILT_KERNEL_SIZE = 3
    EXCLUDE_VOX_LOWERTHRESH = -1000
    EXCLUDE_VOX_UPPERTHRESH = -500


class prmThresholds(object):
    """Thresholds for expiratory and inspiratory images that define PRM regions.

    Attributes:
        EXP_THRESH (float): HU value marking experitaroy threshold for PRM classification
        INSP_THRESH (float): HU value marking inspiratory threshold for PRM classification
    """

    EXP_THRESH = -856
    INSP_THRESH = -950


class outFileNames(object):
    """Establish standard file names for output.

    Attributes:
        PRM_NORM (str): file name for map of norm regions from PRM
        PRM_EMPH (str): file name for map of emph regions from PRM
        PRM_FSAD (str): file name for map of fSAD regions from PRM
    """

    PRM_NORM = "prm_norm"
    PRM_EMPH = "prm_emph"
    PRM_FSAD = "prm_fsad"
