### Cepheids
include_Cepheids = True
N_galaxies_Cep = 19         # Number of SN-host galaxies
N_anchors_Cep = 2           # Number of anchors
include_MW = True           # Consider the MW-Cepheids and their parallax
fixed_zp = False            # Force the parallax zero point to a given value and don't fit it
zp = -0.014                 # zp in as. ONLY IF fixed_zp = True
sig_zp = 0.005		        # zp uncertainty in as. ONLY IF fixed_zp = True
fixed_Zw = False            # Force the metallicity slope Zw to a given value and don't fit for it
Zw = 0                      # Slope of the metallicity effect on the PLR in dex. ONLY IF fixed_Zw = True
sig_Zw = 0		            # Zw uncertainty in dex, ONLY IF fixed_Zw = True
PLR_break = False           # Allow for a break in the PLR at P=10d
break_P = 10                # Period of the PLR-break in days. Also define the pivot period.
added_scatter = 0           # Add dispersion to the Cepheids (See Moertsell et al. 2021, [2021arXiv210511461M])
RLB_correction = False      # Correct for the RLB (See Anderson 2019, [2019A&A...631A.165A])
Kcorr_Cep = False           # Correct for the K-corrections (See Anderson 2022, [2022A&A...658A.148A])


### TRGB
include_TRGB = False
N_galaxies_TRGB = 12        # Number of SN-host galaxies
N_anchors_TRGB = 1          # Number of anchors
Kcorr_TRGB = False          # Correct for the K-corrections (See Anderson 2022, [2022A&A...658A.148A])
use_color = True            # Consider the color term (V-I) for the TRGB (See Anand et al. 2021, [2021AJ....162...80A])
mid_VI = 1.32               # Pivot color for the color term, usually the value from N4258

### Cepheids + TRGB (only if include_Ceph & include_TRGB = True)
different_mu = True        # Allow a different distance for mu_Cep and mu_TRGB

### SNe
fit_aB = True
aB = 0.715840               # Value for the Hubble diagram's intercept ONLY IF fit_aB = False
sig_aB = 0.001631           # Uncertainty on the intercept
z_min = 0.023               # Min redshift to consider when fitting for a_B
z_max = 0.150               # Max redshift to consider when fitting for a_B

### Outlier
outlier_rejection = True    # Include kappa-clipping outlier rejection
kappa = 2.7                 # Value for the kappa-clipping process

### Plot
show_plots = True           # If you want to display the plot

### Physics constants
c = 299792.458              # km/s
q0 = -0.55
j0 = 1
