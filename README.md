# topological_PRM_for_lungHRCT

## Description
A pipeline for analyzing inspiratory and expiratory HRCT lung images through quantitative metrics and 3D maps. Currently, this pipeline provides 3D parametric response maps (PRM), global topology metrics, and 3D maps of local topology metrics for regions characterized by normal lung structure, emphysema, functional small airways disease (fSAD), and emptying emphysema. Pipeline can also process global and local topology for an existing PRM map. Analysis methods are adapted from methodology described in Hoff et al. (2017)<sup>1</sup>.

NOTE: this pipeline calculates topology metrics from Minkowski functionals computed with the QuantImPy Python package<sup>2</sup>, which may differ from Minkowski measures in Hoff et al. (2017)<sup>1</sup>.


## Table of contents
1. [Installation](#1-installation)
2. [Usage](#2-usage)
3. [Outputs](#3-outputs)
4. [Authors and acknowledgement](#4-authors-and-acknowledgment)

## 1. Installation
The following instructions are intended for Mac and Linux systems. Prior to installation, download or clone this git repository. This pipeline is written for use with Python 3.9.1.

### 1.1 Conda installation
We recommend creating a virtual environment to insall the necessary Python libraries required to run this pipeline. A conda distribution is required to create a virtual environment. If conda is already installed on your system, skip this section.

#### 1.1.1 Conda installation: Intel Mac and Linux systems
The Anaconda or Miniconda distirbution can be installed. To install the Anaconda distribution, follow the instructions in this [link](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html). To install the Miniconda distribution, follow the instructions in this [link](https://docs.anaconda.com/free/anaconda/install/)

#### 1.1.2 Conda installation: Apple Silicon Mac
Install Miniforge using step 2 of this [link](https://caffeinedev.medium.com/how-to-install-tensorflow-on-m1-mac-8e9b91d93706)

### 1.2 Install Python packages in a virtual environment

#### 1.2.1 Create a virual environment
Create a conda virtual environment using the following command in terminal, where `<envName>` can be any environment name of your choosing.
```bash
conda create --name <envName> python=3.9.1
```
Confirm that the environment is available:
```bash
conda env list
```
To activate the environment:
```bash
conda activate <envName>
```

#### 1.2.2 Install Python packages
Before installing packages, activate the virtual environment. The list of required packages can be found in `setup/requirements.txt`. To install the pacakges, navigate to the main program directory in terminal and execute the following command.
```bash
pip install -r setup/requirements.txt
```

## 2. Usage
The pipeline currently takes the following required inputs:
- Registered inspiratory and expiratory HRCTs in hounsfield units (HU) (supported file formats: .nii)
- Segmentation mask with positive integers denoting regions of lung parenchyma (supported file formats: .nii)

### 2.1 Create subject configuration file(s)
Configuration (config) files (.ini) specify subject ID, inspiratory HRCT file path, expiratory file path, input PRM map (if available), and the path to save output files to. If the path to an existing PRM map is specified in the config file, the pipeline will skip generation of PRM map from HRCTs.
<br /> 
<br />`config/config_demo.ini` shows config file structure and fields. If an existing PRM map is provided, the required fields are: subjID, inFilePrm, and outDir.
 <br />
<br />`scripts/write_config.py` can be used to create a config file for a single subject and can be adapted to create config files for a batch of subjects.
<br />
<br />If an existing PRM map is provided, the required fields are: `subjID`, `inFilePrm`, and `outDir`.
<br />
<br />If no existing PRM map is provided, the required fields are: `subjID`, `inFileExp`, `inFileInspReg`, `inFileMask`, and `outDir`.

### 2.2 Process a single subject
To process a single subject run, active your virtual environment, navigate to the main program directory, and run the following command.
```bash
python main.py --config <path-to-subject-config-file>
```

### 2.3 Process a batch of subjects
First create a config file for each subejct and place them all in one directory. To process the batch of subjects, activate your virtual environment, navigate to the main program directory in terminal, and run the following command in terminal.
```bash
python main.py --batch --config <path-to-config-file-directory>
```

### 2.4 Compute only PRM maps and global topology metrics
Processing only PRM maps and global topology metrics significantly cuts down computation time. This is useful if local topology maps are not needed. To do this, add the following flag in the command line when processing a single subject or a batch: `--glbl`
```bash
python main.py --glbl --config <path-to-subject-config-file> --glbl
```
or
```bash
python main.py --glbl --batch --config <path-to-config-file-directory> --glbl
```

## 3. Outputs
- PRM
    -  Separate 3D PRM maps of normal lung structure, emphysema, fSAD, and emptying emphysema (.nii)
    - Combined 3D PRM map of normal lung structure, emphysema, fSAD, and emptying emphysema (.nii). PRM classifications are assigned the following values -> normal: 1, fSAD: 2, emphysema: 3, emptying emphysema: 4
    - Colorcoded PRM image of a representative slice along the anterior-posterior dimension (.png)
    - Percentage of lung parenchyma voxels in each PRM classification (.csv)
- Topology
    - Separate 3D maps of local topology density metrics for voxels characterized by normal lung structure, emphysema, fSAD, and emptying emphysema (.nii)
    - Images of surface area density for a reprsentative slice along the anterior-posterior for all PRM maps (.png)
    - Global and mean local topology metrics and whole-lung mean local topology metrics for each PRM classification (.csv)
    - Note: output topology density metrics have the following units (pipeline assumes voxel dimensions in HRCT NIfTI header are in mm)
        - Fractional volume (volume density): unitless
        - Surface area density: m<sup>-1</sup>
        - Integral mean curvature density: m<sup>-2</sup>
        - Euler-Poincare characteristic density: unitless


## 4. Authors and acknowledgment
Sources:
>1. Hoff, B.A., Pompe, E., Galbán, S. et al. CT-Based Local   Distribution Metric Improves Characterization of COPD. Sci Rep 7, 2999 (2017). https://doi.org/10.1038/s41598-017-02871-1

>2. Arnout M.P. Boelens, and Hamdi A. Tchelepi, QuantImPy: Minkowski functionals and functions with Python. SoftwareX, Volume 16 100823, ISSN 2352-7110 (2021). https://doi.org/10.1016/j.softx.2021.100823

Related publications:
>1. Geng, R., Grimbaldeston, M.A., Coimbra, A., Bell, L.C. et al. PRM and Topological Imaging Features From HRCT Can Identify PRISm Patients Who Progress to COPD. Annual ATS meeting, Washington D.C., May 19-24, 2023. Oral Presentation. 

Developers: Aryil Bechtel, Ruiqi Geng<br /> 
Credits: Laura Bell <br /> 
Correspondence: Aryil Bechtel (bechtel.aryil@gene.com), Laura Bell (bell.laura@gene.com)
