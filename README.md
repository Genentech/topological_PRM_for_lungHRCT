# HRCT_lung_image_analysis_pipeline

## Description
A pipeline for analyzing inspiratory and expiratory HRCT lung images through quantitative metrics and 3D maps. Currently, this pipeline provides 3D parametric response maps (PRM) and global topology metrics for regions characterized by normal lung structure, emphysema, functional small airways disease (fSAD), and emptying emphysema.

## Table of contents
1. [Installation](#installation)

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
Confirm that the environment is available.
```bash
conda env list
```
To activate the environment.
```bash
conda activate <envName>
```

#### 1.2.2 Install Python packages
Before installing packages, activate the virtual environment. The list of required packages can be found in `setup/requirements.txt`. To install the pacakges, navigate to the main program directory and execute the following command.
```bash
pip install -r setup/requirements.txt
```

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
