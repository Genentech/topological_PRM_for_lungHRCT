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
        args: ArgumentParse object containing command line flags
    """

    # record start time for processing subject
    t1 = time.perf_counter()

    subject = Subject(config)

    logging.info("*****Processing subject %s*****" % subject.subjID)

    if config.has_option("io", "inFileExp"):
        # if HRCT file location is specified in config, generate PRM maps from HRCT

        # generate PRM maps
        logging.info("Generating PRM maps")
        subject.readCtFiles()
        subject.orientImages()
        subject.applyMedFilts()
        subject.excludeVoxels()
        subject.classifyVoxelsPrm()
        subject.savePrmNiis()
    elif config.has_option("io", "inFilePrm"):
        # if PRM file location is specified in config, read in PRM map

        logging.info("Reading in PRM maps")
        # subject.readPrmFile()
        # subject.genMaskFromPrm()

    # calculate PRM stats and plot representative slice of PRM map
    subject.calcPrmStats()
    subject.plotPrmColor()

    if not args.glbl:
        # if 'glbl' flag not specified, process local topology

        # generate local PRM topology maps
        logging.info("Generating PRM topology maps")
        subject.genLocalTopoMaps()
        subject.saveTopoNiis()
        subject.calcMeanLocalTopoStats()
        subject.plotTopoColor()

    # calculate global topology metrics and save topology stats
    logging.info("Calculating global topology metrics")
    subject.calcTopologyGlobal()
    subject.saveTopologyStats()

    logging.info("Program complete")

    # record end time for processing subject
    t2 = time.perf_counter()
    elapsedTime = (t2 - t1) / 60

    logging.info(f"Program runtime: {elapsedTime} mins")


# def genPrmTopoMaps(config):
#     """Run full pipeline from registered HRCT files.

#     Calculate PRM maps, global topology metrics, and local topology maps.
#     """

#     # record start time for processing subject
#     t1 = time.perf_counter()

#     subject = Subject(config)

#     # generate PRM maps
#     logging.info("Generating PRM maps")
#     subject.readCtFiles()
#     subject.orientImages()
#     subject.applyMedFilts()
#     subject.excludeVoxels()
#     subject.classifyVoxelsPrm()
#     subject.calcPrmStats()
#     subject.savePrmNiis()
#     subject.plotPrmColor()

#     # calculate global topology metrics
#     logging.info("Calculating global topology metrics")
#     subject.calcTopologyGlobal()

#     # generate PRM topology maps
#     logging.info("Generating PRM topology maps")
#     subject.genLocalTopoMaps()
#     subject.saveTopoNiis()
#     subject.calcMeanLocalTopoStats()
#     subject.plotTopoColor()
#     subject.saveTopologyStats()

#     logging.info("Program complete")

#     # record end time for processing subject
#     t2 = time.perf_counter()
#     elapsedTime = (t2 - t1) / 60

#     logging.info(f"Program runtime: {elapsedTime} mins")


# def genPrmGlobalTopo(config):
#     """Run just PRM and global topology pipeline from registered HRCT files.

#     Calculate PRM maps and global topology metrics.
#     """
#     # record start time for processing subject
#     t1 = time.perf_counter()

#     subject = Subject(config)

#     # generate PRM maps
#     logging.info("Generating PRM maps")
#     subject.readCtFiles()
#     subject.orientImages()
#     subject.applyMedFilts()
#     subject.excludeVoxels()
#     subject.classifyVoxelsPrm()
#     subject.calcPrmStats()
#     subject.savePrmNiis()
#     subject.plotPrmColor()

#     # calculate global topology metrics
#     logging.info("Calculating global topology metrics")
#     subject.calcTopologyGlobal()
#     subject.saveTopologyStats()

#     logging.info("Program complete")

#     # record end time for processing subject
#     t2 = time.perf_counter()
#     elapsedTime = (t2 - t1) / 60

#     logging.info(f"Program runtime: {elapsedTime} mins")


def main():
    """Run PRM and topological mapping HRCT analysis.

    If no batch processing indicated, run analysis with single config file.
    If batch processing indicated, run analysis on each config file in specified directory.
    """

    if not args.batch:
        # read in config file
        config = ConfigParser()
        config.read(args.config)

        # run pipeline from registered HRCT files
        logging.info("*****Processing subject %s*****" % config["subjInfo"]["subjID"])
        if not args.glbl:
            # run full pipeline from registered HRCT files
            logging.info(
                "Running full pipeline: PRM maps, global topolgoy metrics, local topology maps"
            )
            genPrmTopoMaps(config)
        else:
            # run just PRM mapping and global topology metrics
            logging.info("Running just PRM maps, global topology metrics")
            genPrmGlobalTopo(config)
    else:
        # get list of config files in specified directory
        configList = glob.glob(join(args.config, "*.ini"))
        for configDir in configList:
            # read in config file
            config = ConfigParser()
            config.read(configDir)

            # run pipeline from registered HRCT files
            logging.info(
                "*****Processing subject %s*****" % config["subjInfo"]["subjID"]
            )
            if not args.glbl:
                # run full pipeline from registered HRCT files
                logging.info(
                    "Running full pipeline: PRM maps, global topolgoy metrics, local topology maps"
                )
                genPrmTopoMaps(config)
            else:
                # run just PRM mapping and global topology metrics
                logging.info("Running just PRM maps, global topology metrics")
                genPrmGlobalTopo(config)


if __name__ == "__main__":
    main()
