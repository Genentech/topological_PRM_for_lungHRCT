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


def prmTopoMapping(config):
    """Run pipeline from registered HRCT files."""

    subject = Subject(config)

    # generate PRM maps
    logging.info("Generating PRM maps")
    subject.readCtFiles()
    subject.dimOutsideVoxels()
    subject.orientImages()
    subject.applyMedFilts()
    subject.excludeVoxels()
    subject.classifyVoxelsPrm()
    subject.calcPrmStats()
    subject.savePrmNiis()
    subject.genPrmColor()
    logging.info("PRM complete")

    # generate global and local topology metrics and maps
    logging.info("Generating topological maps")
    subject.calcTopologyGlobal()
    subject.saveTopologyStats()
    logging.info("Topological mapping complete")


def main():
    """Run PRM and topological mapping HRCT analysis."""

    # read in config file
    config = ConfigParser()
    config.read(args.config)

    # run pipeline from registered HRCT files
    prmTopoMapping(config)


if __name__ == "__main__":
    main()
