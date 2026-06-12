# read-rheol 

Brice Saint-Michel, Telespazio Belgium S.R.L for ESA.

This programme aims to provide a unified (as best as possible) framework to read rheological data coming from different rheometer manufacturers (Anton Paar, Malvern / Netzsch, TA Instruments), i.e. to take the `.csv` or `.txt` exports of their respective softwares (the dreaded Anton Paar Rheocompass, rSpace, Trios) and spit out :

- a `data` table  containing all the relevant variables as a function of time (caveat : this means that you have to include them in the export template, e.g. for Kinexus rheometers), barring maybe LAOS (Large Amplitude Oscillatory Shear) data that might become its own table.
- a `units` dictionary listing all the units for the table data
- a `info` dictionary listing all the metadata of the rheometer acquisition, i.e. everything that is not the table in the `.csv` or `.txt` file. This file does not exist for Malvern rheometers.

Most of the heavy lifting is done by the [`rheol_functions.py`](./rheol_functions.py), notably by the `read_trios` and `format_trios` (`read_antonpaar`, `format_antonpaar`, ...) functions. You can straight up import these functions and call them.

## Notes :

- go to [`read_rheol.ipynb`](./read_rheol.ipynb) to read your files and plot the results.

- for now, the code only works for TA Instruments TRIOS files. If you wish to use the code for other brands of rheometer, please use a previous version of this code (latest `cleanup` commit from 2023). This version of the code is different and ignores most of the file metadata.

- There are a few demonstration files to check if the code is working on this package.

- Anton Paar Rheocompass natively formats their files with UTF-16 LE, which causes a ton of compatibility issues. 

- Don't forget to specify the correct separator when reading your input file with `read_trios`, since Trios plays some tricks on the column separators (`.txt` files should go for `\t`, `.csv` files should go for `,`)
