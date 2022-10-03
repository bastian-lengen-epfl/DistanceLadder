# Distance ladder
This folder contains all the scripts to process multivariate regression on the Cepheid- and/or TRGB-based distanced 
ladder. The models used are initially based on Riess (2021b) but modified according to our needs.

## Requirement
The following packages are required to run this code:
* `Python`
* `Numpy`
* `Pandas`
* `Scipy`
* `Matplotlib`
* (`Astropy`) for the utility functions

---
## Fit setup and data pre-processing
Two things have to be done before fitting the distance ladder:
* Modify the `fit_parameters.py` in order to choose what model you want to fit.
* Include your data in a `/data_tmp` folder. According to your needs the following `.csv` have to be in this folder:
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

For the `/data/Cepheids.csv`the columns have to be in this order:
1. Galaxy host
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you don't use the K-corrections

For the `/data/Cepheids_MW.csv`:
1. Galaxy host. Here the name stands for the dataset used (`'MW1'` for dataset1, `'MW2'` for dataset2, ...).
2. Log10 of the pulsation period
3. Wesenheit magnitude
4. Uncertainty of the Wesenheit magnitude
5. Metallicity [Fe/H] or [O/H] (it has to be consistent between all Cepheids)
6. Redshift (Here `z=0`)
7. Color term V-I, only used for K-corrections. It can be a column full of `NaN` if you don't use the K-corrections
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
For the `/data/SNe_Cepheids.csv` and `/data_tmp/SNe_TRGB.csv`.
1) Galaxy host
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude

For the `/data/SNe_Hubble.csv`:
1) SN name
2) Apparent B peak magnitude
3) Uncertainty of the apparent B peak magnitude
4) Redshift
5) Uncertainty of the redshift


### Generating the different `.csv`
It is possible to generate these `.csv` from the preloaded data in the `/data_static` folder:

| Folder                   | Author                        | Bibcode               | Output                                                     |
|--------------------------|-------------------------------|-----------------------|------------------------------------------------------------|
| `/data_static/Riess/`    | Riess et al. (2016)           | [2016ApJ...826...56R] | Cepheids.csv<br/>Cepheids_anchors.csv<br/>SNe_Cepheids.csv |
| `/data_static/Riess/`    | Riess et al. (2019)           | [2019ApJ...876...85R] | Cepheids_anchors.csv (add LMC Cepheids)                    |
| `/data_static/Riess/`    | Riess et al. (2021a)          | [2021ApJ...908L...6R] | Cepheids_MW.csv                                            |
| `/data_static/Pantheon/` | Scolnic et al. (2018)         | [2018ApJ...859..101S] | SNe_Hubble.csv                                             |
| `/data_static/Anand/`    | Anand et al. (2021)           | [2021AJ....162...80A] | Cepheids_MW.csv                                            |
| `/data_static/H1PStars/` | Cruz Reyes & Anderson  (2022) | [TBD]                 | Cepheids_MW.csv (add OR erase those from R21)              |



To use these data, the script `/other_functions/setup.py` is available. This will generate the `.csv` files in the 
expected format from these datasets. To use this feature, simply run, from the `./` directory, the command:

`python3 other_functions/setup.py RH [--dir 'path_to_data_static'] [--dirdata 'path_to_data']` . 

You can specify what dataset you want to use for the MW Cepheids. R for the R21 dataset, H for the H1Pstars dataset
or RH for both dataset.

---
## Fit

### Fit parameters
In order to obtain the desired fit, it is best to modify the various parameters in the `/code/fit_parameters.py` file. 
This file contains all the different parameters that can be modified to obtain the desired fit.

In particular, it is possible to choose whether you want to make a Cepheid-based distance ladder, TRGB-based distance 
ladder or common distance ladder fit. Many other parameters are determined in this file, such as the consideration of 
the Redshift Leavitt Bias, an outlier rejection algorithm, the fit of the redshift-magnitude diagram, etc... The easiest
way is to look at the file, each parameter is explained in it.

### Single fit
To run the simulation with the data from the `/data/` folder and the parameters from the `/code/fit_parameters.py`,
just run, from the `./` directory, the command:

`python3 code/run.py fit_name [--dir 'path_to_work_dir]'`

The results will be saved in your working directory (by default `./work_dir/`).

### Multiple fit
It is also possible to run multiple simulations. That is, to run several simulations by scanning different values of a 
specific parameter. The main utility of this feature is to be able to scan multiple values for a second PLR break for 
long period Cepheids. This has also been used in the past to measure the impact of different parameters on the 
K-correction Cepheids and TRGB. 

As before, you must have the data in the requested format in the `/data/` folder, and modify the 
`/code/fit_parameters.py` parameters to perform the desired simulations. 

Then, just run, from the `./` directory, the command:

`python3 code/multiple_run.py fit_name [--dir 'path_to_work_dir]'`

The results will be saved in your working directory (by default `./work_dir/`).