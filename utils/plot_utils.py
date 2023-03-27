"""Utils for plotting"""
import matplotlib.pyplot as plt
import numpy as np

import constants


def plotPrmRgbImage(prmRgbArray: np.array, path: str):
    """Plot single slice of RGB color array.

    Args:
        prmRgbArray (np.array): RGB image denoting different PRM regions of lung
        path (str): path to save RGB image to
    """

    plt.figure()
    plt.imshow(prmRgbArray[:, constants.prmProcessing.PLOT_SLICENUM, :, :], cmap="gray")
    plt.axis("off")
    plt.savefig(path, transparent=True, bbox_inches="tight", pad_inches=-0.05, dpi=600)
    plt.clf()
    plt.close()
