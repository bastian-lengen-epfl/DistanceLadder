### Cepheids
include_Cepheids = True
N_galaxies_Cep = 19         # Number of SN-host galaxies
N_anchors_Cep = 2           # Number of anchors
include_MW = True           # Considers the MW-Cepheids and their GAIA parallax
fixed_zp = False            # False -> Fit for the GAIA zp
zp = -0.014                 # GAIA zp in μas, ONLY IF fixed_zp = True
sig_zp = 0.005		        # zp uncertainty, ONLY IF fixed_zp = True
fixed_Zw = False            # False -> Fit for the metallicity effect on the PLR
Zw = 0                      # Slope of the metallicity effect on the PLR ONLY IF fixed_Zw = True
sig_Zw = 0		            # Zw uncertainty, ONLY IF fixed_Zw = True
PLR_break = False           # Allows for a break in the PLR at P=10d
break_P = 10                # Period of the PLR-break, also define the P0 for the term MW*(logP-logP0)
added_scatter = 0           # Add dispersion to the Cepheids (0.0682 in Moertsell et al. 2021)
RLB_correction = False      # Correct for the RLB (Anderson 2019)
Kcorr_Cep = False           # Correct for the K-corrections (Anderson 2021)


### TRGB
include_TRGB = False
N_galaxies_TRGB = 12        # Number of SN-host galaxies
N_anchors_TRGB = 1          # Number of anchors
Kcorr_TRGB = False          # Correct for the K-corrections (Anderson 2021)
use_color = True            # Consider or not the color (V-I) of the TRGB (like Anand 2021)
mid_VI = 1.32               # zero point for the color term (1.32 corresponds to the value for N4258 from Anand (2021)

### Cepheids + TRGB (only if include_Ceph & include_TRGB = True)
different_mu = True        # Allows a different distance for mu_Cep and mu_TRGB

### SNe
fit_aB = True
aB = 0.715840               # Value for the Hubble diagram's intercept ONLY IF fit_aB = False
sig_aB = 0.001631           # 0.71273 +/- 0.00176 Riess 2016, OR can be fitted
z_min = 0.023               # Min redshift to consider when fitting for a_B
z_max = 0.150               # Max redshift to consider when fitting for a_B

### Outlier
outlier_rejection = True    # 2.7-clipping process
kappa = 2.7                 # Value for the kappa-clipping process

### Plot
show_plots = True           # If you want to display the plot

### Physics constants
c = 299792.458              # km/s
q0 = -0.55
j0 = 1