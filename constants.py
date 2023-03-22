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
