# Distance ladder

This folder contains all the scripts to process multivariate regression on the Cepheid- and/or TRGB-based distanced 
ladder. The models used are initially based on Riess (2021b) but modified according to our needs.
## Fit setup and data pre-processing
Two things have to be done before fitting the distance ladder:
* Modify the `fit_parameters.py` in order to choose what model you want to fit.
* Include your data in the `/data_static` folder. According to your needs the following `.csv` have to be in this folder:
  * `/data_static/Cepheids.csv` containing the data from Cepheids in SN-host galaxies.
  * `/data_static/Cepheids_MW.csv` containing the data from the MW-Cepheids.
  * `/data_static/Cepheids_anchors.csv` containing the data from the Cepheids in anchor galaxies.
  * `/data_static/TRGB.csv` containing the data for TRGB in SN-host galaxies.
  * `/data_static/TRGB_anchors.csv` containing the data from the TRGB in anchor galaxies.
  * `/data_static/SNe_Cepheids.csv` containing the data from SNe in Cepheid-host galaxies.
  * `/data_static/SNe_TRGB.csv` containing the data from SNe in TRGB-host galaxies.
  * `/data_static/SNe_Hubble.csv` containing the data from high-z SNe for the redshift-magnitude diagram.

### Cepheids data 
The `.csv` have to be in a specific form in order to run the code. The columns name are not revelant but their order is.

For the `/data_static/Cepheids.csv`the columns have to be in this order:
1. Galaxy host
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you don't use the K-corrections

For the `/data_static/Cepheids_MW.csv`:
1. Galaxy host. Here the name stands for the dataset used (`'MW1'` for dataset1, `'MW2'` for 2, ...).
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift (Here `z=0`)
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you don't use the K-corrections
8. Parallax
9. Uncertainty of the parallax

For the `/data_static/Cepheids_anchors.csv`:
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
For the `/data_static/TRGB.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. Color term V-I. It can be a column full of `NaN` if the model choose don't use the color term.

For the `/data_static/TRGB_anchors.csv`:
1. Galaxy host 
2. Observed I-band magnitude of the TRGB
3. Uncertainty of the observed I-band magnitude
4. Absorbtion in the I-band
5. Redshift of the host galaxy
6. Color term V-I. It can be a column full of `NaN` if the model choose don't use the color term.
7. Geometric distance modulus of the anchor galaxy
8. Uncertainty of the geometric distance modulus

### SNe data
For the `/data_static/SNe_Cepheids.csv` and `/data_static/SNe_TRGB.csv`.
1) Galaxy host
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude

For the `/data_static/SNe_Hubble.csv`:
1) SN name
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude
4) Redshift
5) Uncertainty of the redshift

### Note
If you want to generate these `.csv` from the pre-loaded data from `/data_static/Riess/`, `/data_static/Pantheon/`,
`/data_static/Anand/`, `/data_static/H1PStars/`, a script is available in the `/other_functions/` folder. 

Run, from the `./` directory, the command:
`python3 other_functions/setup.py RH [--dir 'path_to_data_static']` . 

Or the command: `python3 setup.py RH `

You can choose
if you want to use the MW Cepheids from R21 (R), from H1PStars (H) or both (RH). 

## Fit
To run the simulation with the data from the `/data_static/` folder and the parameters from the `fit_parameters.py`,
just run, from the `./` directory, the command:

`python3 code/run.py fit_name [--dir 'path_to_work_dir]'`

The results will be saved in your working directory (by default `./work_dir/`).
