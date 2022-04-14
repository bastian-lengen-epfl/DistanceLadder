### Cepheids ###
include_Cepheids = True
N_galaxies_Cep = 19         # Number of SN-host galaxies
N_anchors_Cep = 2           # Number of anchors
include_MW = True           # Include the MW-Cepheids and their parallax
multiple_zp = True          # Fits for different zp for different datasets
fixed_zp = False            # Force the parallax zero point to a given value and don't fit it (only  works with 1 zp)
zp = -0.014                 # zp in as. ONLY IF fixed_zp = True
sig_zp = 0.005		        # zp uncertainty in as. ONLY IF fixed_zp = True
fixed_Zw = False            # Force the metallicity slope Zw to a given value and don't fit for it
Zw = 0                      # Slope of the metallicity effect on the PLR in dex. ONLY IF fixed_Zw = True
sig_Zw = 0		            # Zw uncertainty in dex, ONLY IF fixed_Zw = True
PLR_break = True            # Include a break in the PLR at P = break_P
break_P = 10                # Period of the PLR-break in days. It also defines the pivot period.
PLR_break2 = True           # Include a second break in the PLR at P = break_P2 for SN-host galaxies
break_P2 = 35               # Period of the second PLR-break in days
added_scatter = 0           # Add dispersion to the Cepheids (See Moertsell et al. 2021, [2021arXiv210511461M])


### TRGB ###
include_TRGB = False
N_galaxies_TRGB = 12        # Number of SN-host galaxies
N_anchors_TRGB = 1          # Number of anchors
use_color = True            # Consider the color term (V-I) for the TRGB (See Anand et al. 2021, [2021AJ....162...80A])
mid_VI = 1.32               # Pivot color for the color term, usually the value from N4258


### Cepheids + TRGB ###
# only if include_Ceph & include_TRGB = True
different_mu = True         # Allow a different distance for mu_Cep and mu_TRGB


### SNe ###
fit_aB = True
aB = 0.715840               # Value for the Hubble diagram's intercept ONLY IF fit_aB = False
sig_aB = 0.001631           # Uncertainty on the intercept
z_min = 0.023               # Min redshift to consider when fitting for aB
z_max = 0.150               # Max redshift to consider when fitting for aB


### Outlier ###
outlier_rejection = True    # Include kappa-clipping outlier rejection
kappa = 2.7                 # Value for the kappa-clipping process


### Relativistic corrections ###
RLB_correction = True       # Correct for the RLB (See Anderson 2019, [2019A&A...631A.165A])
Kcorr_Cep = True            # Correct for the K-corrections (See Anderson 2022, [2022A&A...658A.148A])
EBV_Cep = 0.0               # E(B-V) term for the K-correction, within [0.0, 0.4].
Kcorr_TRGB = True           # Correct for the K-corrections (See Anderson 2022, [2022A&A...658A.148A])
EBV_TRGB = 0.0              # E(B-V) term for the K-correction, within [0.0, 0.05].
Teff_TRGB = 4200            # Teff term for the K-correction, within [3500, 6000].
logg_TRGB = 0.5             # logg term for the K-correction, within [0, 0.5].
FeH_TRGB = -1.75            # [Fe/H] term for the K-correction, within [-2., -1.5].


### Plot ###
show_plots = True          # If you want to display the plot

### Physics constants
c = 299792.458              # km/s
q0 = -0.55
j0 = 1
R = 0.386

##################################################################
###################### Multiple runs #############################
##################################################################
# These parameters are only used for the multiple_run.py

### Relativistic correction
multiple_Cep = False
EBV_Cep_multi = [0.0, 0.004, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5]    # E(B-V) term for the K-correction, within [0.0, 0.4].
multiple_TRGB = True
EBV_TRGB_multi = [0.0, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5]         # E(B-V) term for the K-correction, within [0.0, 0.05].
Teff_TRGB_multi = [3500, 4200, 6000]        # Teff term for the K-correction, within [3800, 6000].



