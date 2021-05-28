# ctools
This repository hosts part of the refactoring of the ctools pipeline for CTA-GRB-WG:
https://github.com/thomasgas/ctools_pipe.

Specifically the part dedicated to GRB's ignificance computation. 

## Environment
To create a virtual environment with all required dependencies:

```bash
conda env create --name <envname> --file=environment.yaml
```
Note that you should already have anaconda installed: https://www.anaconda.com/
## Configuration file

Under cfg you can find a sample configuration file. Description of each parameter is commented within. This file will serve as input when running the code.

## Models

XML models for GRB afterglows can be obtained with ctools_pipe (https://github.com/thomasgas/ctools_pipe). The position of the source has to be set at (ra,dec)=(0, 0.75), in order to be consistent with the 'region selection' step. 
To produce xml files you are supposed to have, at your disposal fits files of the events. These are needed to extract trigger time for single events. 


## Visibility tables

Visibility tables must be given as an input. These tables are generated with runCatVisibility.py: https://github.com/cta-rta/ctools-sci-grb. Pay attention to the tables details, to make sure they fit with your needs (Moon, threshold altitude)

## Compiling the code

After adjusting the configuration file to your needs, you can run the code as follows:

```bash
python GRB_significance.py -f cfg/config.yaml
```

