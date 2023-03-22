"""Write config file."""
from configparser import ConfigParser

config = ConfigParser()
config["inFiles"] = {
    "exp": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDST2OP/data/sandbox/test_CODPGene_registered/10005Q/10005Q_EXP_STD_NJC_COPD.10005Q_INSP_B31f_330_NJC_COPD2.warped.nii.gz",
    "insp_reg": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDST2OP/data/sandbox/test_CODPGene_registered/10005Q/10005Q_EXP_STD_NJC_COPD.10005Q_INSP_B31f_330_NJC_COPD2.warped.nii.gz",
    "mask": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDST2OP/data/sandbox/test_CODPGene_registered/10005Q/10005Q_EXP_STD_NJC_COPD.njh.lung_segmentation.nii.gz",
}

with open("config.ini", "w") as f:
    config.write(f)
