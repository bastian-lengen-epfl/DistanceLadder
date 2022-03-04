# How to use the functions
These functions are available to convert data from Riess et al. (2016, 2019, 2021), Anand et al. (2021) and the Pantheon dataset to the expected format.

## Riess_to_Data.py
To convert the files (.csv/.fits) from R16, R19 and R21 to the expected format, run from the working directory the command:
`python3 Functions/Usefull/Riess_to_Data.py R16_Cepheids.csv R16_SNe.csv R19_LMC.csv R21_MW.csv`

Example for my files:
`python3 Functions/Usefull/Riess_to_Data.py Data/Riess/Cepheids_R16.fits Data/Riess/SN_R16_Table5.csv Data/Riess/Cepheids_LMC_R19.fits Data/Riess/Cepheids_MW_R21.csv`


## Anand_to_Data.py
To convert the files (.csv/.fits)  from Anand et al. (2021) to the expected format, run from the working directory the command:
`python3 Functions/Usefull/Anand_to_Data.py SN_and_TRGB.csv`

Example for my files:
`python3 Functions/Usefull/Anand_to_Data.py Data/Anand/SN_TRGB.csv`


## Pantheon_to_Data.py
To convert the files (.txt) from the Pantheon dataset to the expected format, run from the working directory the command:
`python3 Functions/Usefull/Anand_to_Data.py Pantheon.txt`

Example for my files:
`python3 Functions/Usefull/Anand_to_Data.py Data/Pantheon/Pantheon.txt`
