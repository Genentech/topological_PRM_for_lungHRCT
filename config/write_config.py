"""Write config file."""
from configparser import ConfigParser

config = ConfigParser()
config["inFiles"] = {
    "exp": "path/to/expiratory_image.nii.gz",
    "insp_reg": "path/to/inspiratory_registered_image.nii.gz",
    "mask": "path/to/mask.nii.gz",
}

with open("config/config_demo.ini", "w") as f:
    config.write(f)
