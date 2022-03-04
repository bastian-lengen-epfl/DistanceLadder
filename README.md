# Distance ladder

This folder contains all the script to process multivariate regression on the Cepheid- and/or TRGB-based distanced ladder.

Results Blablabla

## Fit setup and data pre-processing
Two thing have to be done before fitting the distance ladder:
* Modify the `Fit_parameters.py` in order to choose what you want to include in your fit you want.
* Include your data in the `/Data` folder. According to your needs the following `.csv` have to be in this folder:
  * `/Data/Cepheids.csv` containing the data from Cepheids in SN-host galaxies.
  * `/Data/Cepheids_MW.csv` containing the data from the MW-Cepheids.
  * `/Data/Cepheids_anchors.csv` containing the data from the Cepheids in anchor galaxies.
  * `/Data/TRGB.csv` containing the data for TRGB in SN-host galaxies.
  * `/Data/TRGB_anchors.csv` containing the data from the TRGB in anchor galaxies.
  * `/Data/SNe_Cepheids.csv` containing the data from SNe in Cepheid-host galaxies.
  * `/Data/SNe_TRGB.csv` containing the data from SNe in TRGB-host galaxies.
  * `/Data/SNe_Hubble.csv` containing the data from high-z SNe for the redshift-magnitude diagram.

### Cepheids Data 
The `.csv` have to be in a specific form in order to run the code. The columns name are not revelant but their order is.

For the `/Data/Cepheids.csv`the columns have to be in that order:
1. Galaxy host
2. log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all cepheids
6. Redshift
7. Color term V-I, only used for K-corrections, can be a column full of NaN if not necessary

For the `/Data/Cepheids_MW.csv`:
1. Galaxy host (Here gal='MW')
2. log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all cepheids
6. Redshift (Here z=0)
7. Color term V-I, only used for K-corrections, can be a column full of NaN if not necessary
8. Parallax
9. Uncertainty of the parallax

For the `/Data/Cepheids_anchors.csv`:
1. Galaxy host
2. log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all cepheids
6. Redshift
7. Color term V-I, only used for K-corrections, can be a column full of NaN if not necessary
8. Geometric distance modulus of the anchor galaxy
9. Uncertainty of the geometric distance modulus


### TRGB Data 
For the `/Data/TRGB.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. A color term V-I

For the `/Data/TRGB_anchors.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. A color term V-I
7. Geometric distance modulus of the anchor galaxy
8. Uncertainty of the geometric distance modulus

### SNe Data
For the `/Data/SNe_Cepheids.csv` and `/Data/SNe_TRGB.csv`.
1) Galaxy host
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude

For the `/Data/SNe_Hubble.csv`:
1) SN name
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude
4) Redshift
5) Uncertainty of the redshift

### Note
Note that utility functions that convert data from Riess et al. (2016, 2019, 2021), Anand et al. (2021) and the Pantheon dataset in the expected format are available in the `/Functions/Usefull` folder.

## Fit
1. Idea 1
2. Idea 2

''' Some txt '''

3. Idea 3
4. Idea 4

''' Or some other text '''

3. Other Idea 3
4. Other Idea 4 