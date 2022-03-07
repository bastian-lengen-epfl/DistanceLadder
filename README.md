# Distance ladder

This folder contains all the scripts to process multivariate regression on the Cepheid- and/or TRGB-based distanced 
ladder. The models used are initially based on Riess (2021b) but modified according to our needs.
## Fit setup and data pre-processing
Two thing have to be done before fitting the distance ladder:
* Modify the `fit_parameters.py` in order to choose what model you want to fit.
* Include your data in the `/data` folder. According to your needs the following `.csv` have to be in this folder:
  * `/data/Cepheids.csv` containing the data from Cepheids in SN-host galaxies.
  * `/data/Cepheids_MW.csv` containing the data from the MW-Cepheids.
  * `/data/Cepheids_anchors.csv` containing the data from the Cepheids in anchor galaxies.
  * `/data/TRGB.csv` containing the data for TRGB in SN-host galaxies.
  * `/data/TRGB_anchors.csv` containing the data from the TRGB in anchor galaxies.
  * `/data/SNe_Cepheids.csv` containing the data from SNe in Cepheid-host galaxies.
  * `/data/SNe_TRGB.csv` containing the data from SNe in TRGB-host galaxies.
  * `/data/SNe_Hubble.csv` containing the data from high-z SNe for the redshift-magnitude diagram.

### Cepheids data 
The `.csv` have to be in a specific form in order to run the code. The columns name are not revelant but their order is.

For the `/data/Cepheids.csv`the columns have to be in that order:
1. Galaxy host
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you will not use the K-corrections

For the `/data/Cepheids_MW.csv`:
1. Galaxy host (Here `'MW'`)
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift (Here `z=0`)
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you will not use the K-corrections
8. Parallax
9. Uncertainty of the parallax

For the `/data/Cepheids_anchors.csv`:
1. Galaxy host
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you will not use the K-corrections
8. Geometric distance modulus of the anchor galaxy
9. Uncertainty of the geometric distance modulus


### TRGB data 
For the `/data/TRGB.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. Color term V-I. It can be a column full of `NaN` if the model choose don't use the color term.

For the `/data/TRGB_anchors.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. Color term V-I. It can be a column full of `NaN` if the model choose don't use the color term.
7. Geometric distance modulus of the anchor galaxy
8. Uncertainty of the geometric distance modulus

### SNe data
For the `/data/SNe_Cepheids.csv` and `/data/SNe_TRGB.csv`.
1) Galaxy host
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude

For the `/data/SNe_Hubble.csv`:
1) SN name
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude
4) Redshift
5) Uncertainty of the redshift

### Note
Note that utility functions that convert data from Riess et al. (2016, 2019, 2021a), Anand et al. (2021) and the 
Pantheon dataset in the expected format are available in the `/other_functions` folder.

## Fit
To run the simulation with the data from the `/data` folder and the parameters from the `fit_parameters.py`,
just run the command `python3 run.py` from the working directory. 
