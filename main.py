"""Run HRCT PRM + topoligcal mapping pipeline."""
import argparse
import glob
import logging
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
args = parser.parse_args()


def genPrmTopoMaps(config):
    """Run pipeline from registered HRCT files."""

    subject = Subject(config)

    # generate PRM maps
    logging.info("Generating PRM maps")
    subject.readCtFiles()
    subject.normalizeAllCt()
    subject.dimOutsideVoxels()
    subject.orientImages()
    subject.applyMedFilts()
    subject.excludeVoxels()
    subject.classifyVoxelsPrm()
    subject.calcPrmStats()
    subject.savePrmNiis()
    subject.genPrmColor()

    # calculate global topology metrics
    logging.info("Calculating global topology metrics")
    subject.calcTopologyGlobal()
    subject.saveTopologyStats()

    # generate PRM topology maps
    logging.info("Generating PRM topology maps")
    subject.genLocalTopoMaps()
    subject.saveTopoNiis()

    logging.info("Program complete")


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
        genPrmTopoMaps(config)
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
            genPrmTopoMaps(config)


if __name__ == "__main__":
    main()
