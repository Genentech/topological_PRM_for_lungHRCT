"""Run ANTsPy registration on CTs.

Args:
    regfiles (str): path to csv with columns ["sid","fixedCt","movingCt","fixedMask","movingMask"].
                    if multiple rows of file names, registrations are run as a batch process.
                        sid (str): subject ID
                        fixedCt (str): path to fixed CT (.nii)
                        movingCt (str): path to moving CT (.nii)
                        fixedMask (str): path to fixed mask (.nii)
                        movingMask (str): path to movign mask (.nii)
    outdir (str): path to directory to save registered files to

Outputs (files are saved to outDir/sid/):
    -<fixed CT file name>.<moving CT file name>.warped.nii.gz: warped moving CT registered to fixed CT
    -<fixed CT file name>.<moving CT file name>.jac.nii.gz: jacobian determinant image

Usage:
    python /scripts/register_antspy.py --regfiles <path_to_csv> --outdir <path_to_outdir>
"""
import argparse
import concurrent.futures
import time
from os import mkdir
from os.path import exists, join

import ants
import nibabel as nib
import pandas as pd

# set command line flags/args
parser = argparse.ArgumentParser(description="Run batch ANTsPy registration")
parser.add_argument(
    "--regfiles",
    type=str,
    metavar="",
    required=True,
    help="path of csv containing reg file paths",
)
parser.add_argument(
    "--outdir",
    type=str,
    metavar="",
    required=True,
    help="directory to save registered images to",
)
args = parser.parse_args()


def registerCOPDGene(fixedCtFile, movingCtFile, fixedMaskFile, movingMaskFile, outDir):
    """Register movingCtFile to fixedCtFile.

    Args:
        fixedCtFile (str): file path of fixed CT
        movingCtFile (str): file path of moving CT to be registered to fixed CT
        fixedMaskFile (str): file path of fixed CT mask
        movingMaskFile (str): file path of moving CT mask
        outDir (str): path of directory to save output files to
    """

    # read in files
    fixed = ants.image_read(fixedCtFile)
    moving = ants.image_read(movingCtFile)
    fixedMask = ants.image_read(fixedMaskFile)
    movingMask = ants.image_read(movingMaskFile)

    # create copy of original moving HRCT to apply transform to
    movingAntsOG = moving.clone()

    # set masked regions of segmentations =1 (auto converts to np array)
    fixedMask[fixedMask > 1] = 1
    movingMask[movingMask > 1] = 1

    # normalize HRCTs and apply masks
    fixed = ants.iMath(fixed, "Normalize")
    moving = ants.iMath(moving, "Normalize")
    fixed[fixedMask == 0] = 0
    moving[movingMask == 0] = 0

    # perform registration and apply transform to unmasked moving image
    tx = ants.registration(fixed=fixed, moving=moving, type_of_transform="ElasticSyN")
    movingWarped = ants.apply_transforms(
        fixed=fixed, moving=movingAntsOG, transformlist=tx["fwdtransforms"]
    )

    # generate jacobian determinant image
    jac = ants.create_jacobian_determinant_image(fixed, tx["fwdtransforms"][0], 1)

    # if directory for sid doesn't exist in outDir, make one
    sid = fixedCtFile.split("/")[-1].split("_")[0]
    subjDir = join(outDir, sid)
    if not exists(subjDir):
        mkdir(subjDir)

    # convert to Nifti1Images
    movingWarpedNii = ants.to_nibabel(movingWarped)
    jacNii = ants.to_nibabel(jac)

    # save warped moving image and jacobian determinant image as nii
    movingWarpedFile = (
        fixedCtFile.split("/")[-1].split(".")[0]
        + "."
        + movingCtFile.split("/")[-1].split(".")[0]
        + ".warped.nii.gz"
    )
    jacFile = (
        fixedCtFile.split("/")[-1].split(".")[0]
        + "."
        + movingCtFile.split("/")[-1].split(".")[0]
        + ".jac.nii.gz"
    )
    nib.save(movingWarpedNii, join(subjDir, movingWarpedFile))
    nib.save(jacNii, join(subjDir, jacFile))

    movingCtFileOnly = movingCtFile.split("/")[-1]
    print(f"file {movingCtFileOnly} registered")


def main():

    # read command line args
    regFilesDf = pd.read_csv(args.regfiles)
    outDir = args.outdir
    outDirList = [outDir for idx in range(len(regFilesDf))]

    # read in csv of registration files
    fixedCtFiles = regFilesDf["fixedCt"].tolist()
    movingCtFiles = regFilesDf["movingCt"].tolist()
    fixedMaskFiles = regFilesDf["fixedMask"].tolist()
    movingMaskFiles = regFilesDf["movingMask"].tolist()

    t1 = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(
            registerCOPDGene,
            fixedCtFiles,
            movingCtFiles,
            fixedMaskFiles,
            movingMaskFiles,
            outDirList,
        )
    t2 = time.perf_counter()
    elapsedTime = t2 - t1

    print(f"time for all registrations: {elapsedTime} seconds")


if __name__ == "__main__":
    main()
