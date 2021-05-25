# ctools
Significance of GRB afterglows

## Configuration file

Under cfg you can find a sample configuration file. Description of each parameter is commented within. This file will serve as input when running the code.
## Models
XML models for GRB afterglows are obtained with Thomas pipeline (https://github.com/cta-rta/ctools-sci-grb) where the position of the source has been set in (ra,dec)=(0, 0.75) in order to be consistent with the 'region selection' step. 


## Compiling the code

After adjusting the configuration file to your needs, you can run the code as follows:

```bash
python TOW_evth_works.py -f cfg/config.yaml
```
N.B. Some functions are imported from the file TOW_functions.py. Make sure they are in the same directory when compiling.
<HR>

