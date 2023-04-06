"""Write config file."""
from configparser import ConfigParser

config = ConfigParser()
config.optionxform = str  # make keys case sensitive
config["subjInfo"] = {"subjID": "000001"}
config["io"] = {
    "inFileExp": "path/to/expiratory_image.nii.gz",
    "inFileInspReg": "path/to/inspiratory_registered_image.nii.gz",
    "inFileMask": "path/to/mask.nii.gz",
    "outDir": "path/to/main/out_directory/",
}

with open("config/config_demo.ini", "w") as f:
    config.write(f)
