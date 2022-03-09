# How to use the functions
These functions are available to convert data from Riess et al. (2016, 2019, 2021a), Anand et al. (2021) and the Pantheon dataset to the expected format.

## Riess_to_data.py
To convert the files (.csv/.fits) from R16, R19 and R21 to the expected format, run from the working directory the command:
`python3 other_functions/Riess_to_data.py R16_Cepheids.csv R16_SNe.csv R19_LMC.csv R21_MW.csv`

Example for my files:
`python3 other_functions/Riess_to_data.py data/Riess/Cepheids_R16.fits data/Riess/SN_R16_Table5.csv data/Riess/Cepheids_LMC_R19.fits data/Riess/Cepheids_MW_R21.csv
