"""Run HRCT PRM + topoligcal mapping pipeline."""
import logging
from configparser import ConfigParser

from absl import app

from subject_classmap import Subject

# read config file
config = ConfigParser()
config.read("config/config.txt")
import argparse

def prmMapping(config):
    subject = Subject(config)
    subject.readCtFiles()
    subject.dimOutsideVoxels()
    subject.orientImages()
    subject.applyMedFilts()

    logging.info("PRM complete")


def topoMapping(config):
    subject = Subject(config)
    logging.info("Topological mapping complete")


def main(argv):
    logging.info("Generating PRM maps")
    prmMapping(config)

    logging.info("Generating topological maps")
    topoMapping(config)


if __name__ == "__main__":
    app.run(main)
