"""Run HRCT PRM + topoligcal mapping pipeline."""
import argparse
import logging
from configparser import ConfigParser

from subject_classmap import Subject

# set logging level to that logging.info statements are printed
logging.basicConfig(level=logging.INFO)

# set command line flags/args
parser = argparse.ArgumentParser(description="Run HRCT voxel-wise lung image analysis")
parser.add_argument(
    "--config", type=str, metavar="", required=True, help="configuration file path"
)
args = parser.parse_args()


def prmMapping(config):
    subject = Subject(config)
    subject.readCtFiles()
    subject.dimOutsideVoxels()
    subject.orientImages()
    subject.applyMedFilts()
    subject.excludeVoxels()
    logging.info("PRM complete")


def topoMapping(config):
    subject = Subject(config)
    logging.info("Topological mapping complete")


def main():
    # read in config file
    config = ConfigParser()
    config.read(args.config)

    # generate PRM maps
    logging.info("Generating PRM maps")
    prmMapping(config)

    # generate topological maps
    logging.info("Generating topological maps")
    topoMapping(config)


if __name__ == "__main__":
    main()
