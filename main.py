"""Run HRCT PRM + topoligcal mapping pipeline."""
import argparse
import glob
import logging
import time
from configparser import ConfigParser
from os.path import join

from subject_classmap import Subject

# set logging level to that logging.info statements are printed
logging.basicConfig(level=logging.INFO)

# set command line flags/args
parser = argparse.ArgumentParser(description="Run HRCT voxel-wise lung image analysis")
parser.add_argument(
    "--config",
    type=str,
    metavar="",
    required=True,
    help="single subject: config file path. batch processing: directory of config files",
)
parser.add_argument(
    "--batch",
    action="store_true",
    required=False,
    help="indicate to run a batch process",
)
parser.add_argument(
    "--glbl",
    action="store_true",
    required=False,
    help="indicate to calc only prm and global topology",
)
args = parser.parse_args()


def processSubject(config, args):
    """Process a single subject from config file.

    Args:
        config: configuration file containing input/ouput file locations
        args: ArgumentParser object containing command line flags
    """

    # record start time for processing subject
    t1 = time.perf_counter()

    subject = Subject(config)

    logging.info("*****Processing subject %s*****" % subject.subjID)

    if config.has_option("io", "inFilePrm"):
        # if PRM file specified in config, read in PRM map

        logging.info("Reading in PRM map")
        subject.readPrmFile()

        if config.has_option("io", "inFileMask"):
            # if mask file specified in config, use it
            subject.readMaskFile()
        else:
            # if nomask file in config, generate mask from binned regions in prm map
            subject.genMaskFromPrm()

    elif config.has_option("io", "inFileExp"):
        # if HRCT file specified in config, generate PRM maps from HRCT

        # generate PRM maps
        logging.info("Generating PRM maps")
        subject.readCtFiles()
        subject.applyMedFilts()
        subject.excludeVoxels()
        subject.classifyVoxelsPrm()

    # save PRM maps, calculate PRM stats, and plot representative slice of PRM map
    subject.savePrmNiis()
    subject.calcPrmStats()
    subject.plotPrmColor()

    # calculate global topology metrics
    logging.info("Calculating global topology metrics")
    subject.genDictOfImageArrays()
    subject.calcTopologyGlobal()
    subject.saveTopologyStats()

    if not args.glbl:
        # if 'glbl' flag not specified, process local topology

        # generate local PRM topology maps
        logging.info("Generating PRM topology maps")
        subject.genLocalTopoMaps()
        subject.saveTopoNiis()
        subject.calcMeanLocalTopoStats()
        subject.plotTopoColor()

    # save topology stats
    subject.saveTopologyStats()

    logging.info("Program complete")

    # record end time for processing subject
    t2 = time.perf_counter()
    elapsedTime = (t2 - t1) / 60

    logging.info(f"Program runtime: {elapsedTime} mins")


def main():
    """Run PRM and topological mapping HRCT analysis.

    If no batch processing indicated, run analysis with single config file.
    If batch processing indicated, run analysis on each config file in specified directory.
    """

    if not args.batch:
        # read in config file
        config = ConfigParser()
        config.read(args.config)

        # process subject
        processSubject(config, args)

    else:
        # get list of config files in specified directory
        configList = glob.glob(join(args.config, "*.ini"))

        # loop over config files
        for configDir in configList:
            # read in config file
            config = ConfigParser()
            config.read(configDir)

            # process subject
            processSubject(config, args)


if __name__ == "__main__":
    main()
