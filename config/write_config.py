"""Write config file."""
from configparser import ConfigParser

config = ConfigParser()
config.optionxform = str  # make keys case sensitive
config["subjInfo"] = {"subjID": "10161E_phase3"}
config["io"] = {
    "inFileExp": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDGene/data/sandbox/test_CODPGene_registered/10161E_phase3/10161E_EXP_STD_NJC_COPD.10161E_EXP_B31f_351_NJC_COPD3.warped.nii.gz",
    "inFileInspReg": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDGene/data/sandbox/test_CODPGene_registered/10161E_phase3/10161E_EXP_STD_NJC_COPD.10161E_INSP_B31f_351_NJC_LD_COPD3.warped.nii.gz",
    "inFileMask": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDGene/data/sandbox/test_CODPGene_registered/10161E_phase3/10161E_EXP_STD_NJC_COPD.njh.lung_segmentation.nii.gz",
    "outDir": "/Users/bechtea2/Documents/projects/COPD_longitudinal/COPDGene/data/sandbox/test_CODPGene_registered/10161E_phase3",
}

with open("config/config4.ini", "w") as f:
    config.write(f)
