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


class prmProcessing(object):
    """Constants for creating PRM maps.

    Attributes:
        EXP_THRESH (float): HU value marking experitaroy threshold for PRM classification
        INSP_THRESH (float): HU value marking inspiratory threshold for PRM classification
        CLASSIFICATION_NUM_NORM (int): number assigned to PRM norm voxels in combined PRM map
        CLASSIFICATION_NUM_FSAD (int): number assigned to PRM fSAD voxels in combined PRM map
        CLASSIFICATION_NUM_EMPH (int): number assigned to PRM emph voxels in combined PRM map
        CLASSIFICATION_NUM_EMPTEMPH (int): number assigned to PRM emptying emph voxels in combined PRM map
        PRM_BIN2RGB (dict): dictionary mapping PRM bin number classification to RGB value
        PLOT_SLICENUM (int): slice number along the anterior-posterior dimension for plotting PRM color image
    """

    EXP_THRESH = -856
    INSP_THRESH = -950
    CLASSIFICATION_NUM_NORM = 1
    CLASSIFICATION_NUM_FSAD = 2
    CLASSIFICATION_NUM_EMPH = 3
    CLASSIFICATION_NUM_EMPTEMPH = 4
    PRM_BIN2RGB = {1: [0.4, 0.8, 0], 2: [1, 1, 0], 3: [0.8, 0, 0], 4: [0.6, 0.2, 1]}
    PLOT_SLICENUM = 200


class topoMapping(object):
    """Constants for generating PRM topological maps.

    Attributes:
        WIND_RADIUS (int): half the length of one side of nxnxn moving window for calculating local topology, rounded up
        GRID_RES (int): interval (in voxels) between nxnxn moving windows for calculating local topology
    """

    WIND_RADIUS = 11
    GRID_RES = 5


class outFileNames(object):
    """Establish standard file names for output.

    Attributes:
        PRM_NORM (str): file name for map of norm regions from PRM
        PRM_FSAD (str): file name for map of fSAD regions from PRM
        PRM_EMPH (str): file name for map of emph regions from PRM
        PRM_EMPTEMPH (str): file name for map of emptying emph regions from PRM
        PRM_ALL (str): file name for map of all regions from PRM
        PRM_STATS (str): file name for csv containing prm stats
        TOPO_STATS (str): file name for csv containing global and local topology metrics

    Intended that subject ID will be appended to file name.
    """

    PRM_NORM = "prm_norm_"
    PRM_FSAD = "prm_fsad_"
    PRM_EMPH = "prm_emph_"
    PRM_EMPTEMPH = "prm_emptemph_"
    PRM_ALL = "prm_all_"
    PRM_STATS = "prm_stats_"
    TOPO_STATS = "topology_stats_"
