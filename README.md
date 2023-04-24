# HRCT_lung_image_analysis_pipeline

## Description
A pipeline for analyzing inspiratory and expiratory HRCT lung images through quantitative metrics and 3D maps. Currently, this pipeline provides 3D parametric response maps (PRM) and global topology metrics for regions characterized by normal lung structure, emphysema, functional small airways disease (fSAD), and emptying emphysema.

## Table of contents
1. [Installation](#1-installation)
2. [Usage](#2-usage)
3. [Outputs](#3-outputs)
4. [Authors and acknowledgement](#4-authors-and-acknowledgment)
5. [Project status](#5-project-status)

## 1. Installation
The following instructions are intended for Mac and Linux systems. Prior to installation, download or clone this git repository. This pipeline is written for use with Python 3.9.1.

### 1.1 Conda installation
We recommend creating a virtual environment to insall the necessary Python libraries required to run this pipeline. A conda distribution is required to create a virtual environment. If conda is already installed on your system, skip this section.

##### 1.1.1 Conda installation: Intel Mac and Linux systems
The Anaconda or Miniconda distirbution can be installed. To install the Anaconda distribution, follow the instructions in this [link](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html). To install the Miniconda distribution, follow the instructions in this [link](https://docs.anaconda.com/free/anaconda/install/)

##### 1.1.2 Conda installation: Apple Silicon Mac
Install Miniforge using step 2 of this [link](https://caffeinedev.medium.com/how-to-install-tensorflow-on-m1-mac-8e9b91d93706)

### 1.2 Install Python packages in a virtual environment

##### 1.2.1 Create a virual environment
Create a conda virtual environment using the following command in terminal, where `<envName>` can be any environment name of your choosing.
```bash
conda create --name <envName> python=3.9.1
```
Confirm that the environment is available.
```bash
conda env list
```
To activate the environment.
```bash
conda activate <envName>
```

##### 1.2.2 Install Python packages
Before installing packages, activate the virtual environment. The list of required packages can be found in `setup/requirements.txt`. To install the pacakges, navigate to the main program directory in terminal and execute the following command.
```bash
pip install -r setup/requirements.txt
```

## 2. Usage
The pipeline currently takes the following required inputs:
- Registered inspiratory and expiratory HRCTs in hounsfield units (HU) (supported file formats: .nii)
- Segmentation mask with positive integers denoting regions of lung parenchyma (supported file formats: .nii)

### 2.1 Create subject configuration file(s)
Configuration (config) files (.ini) specify subject ID, inspiratory HRCT file path, expiratory file path, and the path to save out files to. 
<br /> 
<br />`config/config_demo.ini` shows the required config file structure and fields.
 <br />
<br />`scripts/write_config.py` can be used to create a config file for a single subject and can be adapted to create config files for a batch of subjects.

### 2.2 Process a single subject
To process a single subject run, active your virtual environment, navigate to the main program directory, and run the following command.
```bash
python main.py --config <path-to-subject-config-file>
```

### 2.3 Process a batch of subjects
First create a config file for each subejct and place them all in one directory. To process the batch of subjects, activate your virtual environment, navigate to the main program directory in terminal, and run the following command in terminal.
```bash
python main.py --config <path-to-config-file-directory>
```

## 3. Outputs
- Separate 3D PRM maps of normal lung structure, emphysema, fSAD, and emptying emphysema (.nii)
- Combined 3D PRM map of normal lung structure, emphysema, fSAD, and emptying emphysema (.nii). PRM classifications are assigned the following values -> normal: 1, fSAD: 2, emphysema: 3, emptying emphysema: 4
- Colorcoded PRM image of a representative slice along the anterior-posterior dimension (.png)
- Percentage of lung parenchyma voxels in each PRM classification (.csv)
- Global topology metrics for each PRM classification: volume, surface area, mean curvature length, Euler-Poincare characteristic (.csv)

## 4. Authors and acknowledgment
Developers: Aryil Bechtel, Riqui Geng<br /> 
Credits: Laura Bell<br /> 
Correspondence: Aryil Bechtel (bechtel.aryil@gene.com), Laura Bell (bell.laura@gene.com)

## 5. Project status
This pipeline is currently regularly updated to improve existing functions and include new analysis features.
