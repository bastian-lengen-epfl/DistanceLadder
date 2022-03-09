'''
This module contains the functions that convert the static_data to the expected format.
Riess et al. (2016)   [2016ApJ...826...56R]     -> Cepheids.csv, Cepheids_anchors.csv and SNe_Cepheids.csv
Riess et al. (2019)   [2019ApJ...876...85R]     -> Cepheids_anchors.csv (add LMC Cepheids to the dataset)
Riess et al. (2021a)  [2021ApJ...908L...6R]     -> Cepheids_MW.csv
Anand et al. (2021)   [2021AJ....162...80A]     -> TRGB.csv, TRGB_anchors and SNe_TRGB.csv
Scolnic et al. (2018) [2018ApJ...859..101S]     -> SNe_Hubble.csv
H1PStars MW Cepheids  [...]                     -> Cepheids_MW.csv (add OR erase those from R21)
'''
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

def Riess_to_data(Cepheids_R16, SNe_R16, Cepheids_LMC_R19, Cepheids_MW_R21,  data_static_dir='../data_static/'):
    '''
    This script converts the data from R16, R19, R21 to the desired format.
    R16 CDS Table 4     -> Cepheids.csv and Cepheids_anchors.csv
    R16 Table 5         -> SNe_Cepheids.csv
    R19 LMC             -> Cepheids_anchors.csv
    R21 MW              -> Cepheids_MW.csv
    The SNe_Cepheids.csv is completed with the redshift from NED.

    :type   Cepheids_MW_R21: str
    :param  Cepheids_R16: Path to the .csv file from R16 (Cepheids).
    :type   SNe_R16: str
    :param  SNe_R16: Path to the .csv file from R16 (SNe).
    :type   Cepheids_LMC_R19: str
    :param  Cepheids_LMC_R19: Path to the .csv file from R19.
    :type   Cepheids_MW: str
    :param  Cepheids_MW: Path to the .csv file from R21.
    :type   data_static_dir: str
    :param  data_static_dir: Path to the data_static directory
    '''
    # Values used
    R = 0.386               # Wesenheit
    Z_sun = 8.824           # Z_sun allows to convert 12+log(O/H) from R16 to [O/H]
    mu_N4258 = 29.397
    sig_mu_N4258 = 0.032
    mu_LMC = 18.477
    sig_mu_LMC = 0.0263
    aB_R16 = 0.71273
    sig_aB_R16 = 0.00176

    # Redshift list from NED
    z_dict = {'MW':0, 'LMC':0.92, 'N4258':1.49, 'M101':0.80, 'N1015':8.77, 'N1309':7.12, 'N1365':5.45, \
              'N1448':3.90, 'N2442':4.89, 'N3021':5.14, 'N3370':4.27, 'N3447':3.56, 'N3972':2.84, \
              'N3982': 3.70, 'N4038':5.48, 'N4424':1.46, 'N4536':6.03, 'N4639':3.40, 'N5584':5.46, \
              'N5917':6.35, 'N7250': 3.89, 'U9391': 6.38, 'M31': -1.00}

    # Start with the R16 Cepheids
    if Cepheids_R16[-4:]=='.csv':
        tmp = pd.read_csv(Cepheids_R16, sep=',')
    elif Cepheids_R16[-5:]=='.fits':
        with fits.open(Cepheids_R16) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {Cepheids_R16} must be either a .csv or .fits')
        return 0
    Cepheids = pd.DataFrame()
    Cepheids['Gal'] = tmp['Gal']
    Cepheids['logP'] = np.log10(tmp['Per'])
    Cepheids['mW'] = tmp['F160W'] - R * tmp['V-I']
    Cepheids['sig_mW'] = tmp['sigTot']
    Cepheids['Fe/H'] = tmp['[O/H]']-Z_sun
    Cepheids['Gal'][Cepheids['Gal'] == 'M101 '] = 'M101' #Rename 'M101 ' to 'M101'
    Cepheids['Gal'][Cepheids['Gal'] == 'M31  '] = 'M31'  # Rename 'M31  ' to 'M31'
    for i in Cepheids.index:
        Cepheids.loc[i, 'z'] = z_dict[Cepheids.loc[i, 'Gal']]*1e-3
    Cepheids['V-I'] = tmp['V-I']
    Cepheids = Cepheids[~(Cepheids['Gal'] == 'M31')].reset_index(drop=True)  # Drop M31
    Cepheids_anchors = Cepheids[Cepheids['Gal'] == 'N4258'].reset_index(drop=True)
    Cepheids_anchors['mu'] = mu_N4258
    Cepheids_anchors['sig_mu'] = sig_mu_N4258
    Cepheids = Cepheids[~((Cepheids['Gal'] == 'N4258'))].reset_index(drop=True) # Drop N4258
    Cepheids = Cepheids[~((Cepheids['Gal'] == 'M101') & (Cepheids['logP'] > np.log10(35)))].reset_index(drop=True) # Drop M101 Cepheids with a period > 35d

    # Add the R19 LMC Cepheids to the anchors
    if Cepheids_LMC_R19[-4:]=='.csv':
        tmp = pd.read_csv(Cepheids_LMC_R19, sep=',')
    elif Cepheids_LMC_R19[-5:]=='.fits':
        with fits.open(Cepheids_LMC_R19) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {Cepheids_LMC_R19} must be either a .csv or .fits')
        return 0
    to_add = pd.DataFrame()
    to_add['Gal'] = ['LMC'] * len(tmp)
    to_add['logP'] = tmp['logP']
    to_add['mW'] = tmp['mWH']
    to_add['sig_mW'] = tmp['e_mWH']
    to_add['Fe/H'] = -0.30
    to_add['z'] = z_dict['LMC']*1e-3
    to_add['V-I'] = tmp['F555Wmag']-tmp['F814Wmag']
    to_add['mu'] = mu_LMC
    to_add['sig_mu'] = sig_mu_LMC
    Cepheids_anchors = pd.concat([Cepheids_anchors, to_add], ignore_index=True, sort=False)

    # R21 MW Cepheids
    if Cepheids_MW_R21[-4:] == '.csv':
        tmp = pd.read_csv(Cepheids_MW_R21, sep=',')
    elif Cepheids_MW_R21[-5:] == '.fits':
        with fits.open(Cepheids_MW_R21) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {Cepheids_MW_R21} must be either a .csv or .fits')
        return 0
    tmp = tmp[tmp['pi_EDR3'] != '⋯'].reset_index(drop=True) # Remove the Cepheids with no parallax
    Cepheids_MW = pd.DataFrame()
    Cepheids_MW['Gal'] = ['MW'] * len(tmp)
    Cepheids_MW['logP'] = tmp['logP']
    Cepheids_MW['mW'] = tmp['m_W']
    Cepheids_MW['sig_mW'] = tmp['sig_m_W']
    Cepheids_MW['Fe/H'] = tmp['[Fe/H]']
    Cepheids_MW['z'] = z_dict['MW']
    Cepheids_MW['V-I'] = tmp['V']-tmp['I']
    Cepheids_MW['pi'] = list(map(float,tmp['pi_EDR3'])) # Also convert str -> float
    Cepheids_MW['sig_pi'] = list(map(float, tmp['sig_pi_EDR3'])) # Idem

    # R16 SNe
    if SNe_R16[-4:]=='.csv':
        tmp = pd.read_csv(SNe_R16, sep=',')
    elif SNe_R16[-5:]=='.fits':
        with fits.open(SNe_R16) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {SNe_R16} must be either a .csv or .fits')
        return 0
    SNe_Cepheids = pd.DataFrame()
    SNe_Cepheids['Gal'] = tmp['Gal']
    SNe_Cepheids['mB'] = tmp['mB+5aB'] - 5 * aB_R16
    SNe_Cepheids['sig_mB'] = np.sqrt(tmp['sig'] ** 2 - (5 * sig_aB_R16) ** 2)


    # Save everything
    Cepheids.to_csv(data_static_dir+'Cepheids.csv',index=False)
    Cepheids_MW.to_csv(data_static_dir+'Cepheids_MW.csv',index=False)
    Cepheids_anchors.to_csv(data_static_dir+'Cepheids_anchors.csv', index=False)
    SNe_Cepheids.to_csv(data_static_dir + 'SNe_Cepheids.csv', index=False)
    return

def H1PStars_to_data(Cepheids_MW, erase = False,  data_static_dir='../data_static/'):
    '''
    This function converts the data from Mauricio to the desired format. It can either add the MW Cepheids to the
    existing `/data_static/Cepheids_MW.csv` or erase it.

    :type   Cepheids_MW: str
    :param  Cepheids_MW: Path to the .csv file from H1PStars.
    :type   erase: bool
    :param  erase: Choose if you want to erase the Cepheids from R21 from the Cepheids_MW.csv, otherwise complete the .csv.
    :type   data_static_dir: str
    :param  data_static_dir: Path to the data_static directory
    '''
    # Load data
    if Cepheids_MW[-4:] == '.csv':
        tmp = pd.read_csv(Cepheids_MW, sep=',')
    elif Cepheids_MW[-5:] == '.fits':
        with fits.open(Cepheids_MW) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {Cepheids_MW} must be either a .csv or .fits')
        return 0
    Cepheids_MW = pd.DataFrame()
    if erase == True:
        Cepheids_MW['Gal'] = ['MW'] * len(tmp)
    else:
        Cepheids_MW['Gal'] = ['MW2'] * len(tmp) # the 2 allows us to differenciate between the two dataset
    Cepheids_MW['logP'] = np.log10(tmp['Period'])
    Cepheids_MW['mW'] = tmp['mh']
    Cepheids_MW['sig_mW'] = tmp['mh_error']
    Cepheids_MW['Fe/H'] = tmp['FE/H']
    Cepheids_MW['z'] = 0
    Cepheids_MW['V-I'] = np.nan  # Don't use K-corrections !!!
    Cepheids_MW['pi'] = tmp['Parallax']
    Cepheids_MW['sig_pi'] = tmp['Parallax_error']

    # ADD to the previously existing
    if erase == False:
        tmp = pd.read_csv(data_static_dir+'Cepheids_MW.csv', sep=',')
        Cepheids_MW = pd.concat([tmp, Cepheids_MW]).reset_index(drop=True)
    # Save it
    Cepheids_MW.to_csv(data_static_dir+'Cepheids_MW.csv', index=False)

    return

def Anand_to_data(TRGB_SNe_Anand, data_static_dir='../data_static/'):
    '''
    This script converts the data from Anand et al. (2021) to the desired format.
    Table 2         -> TRGB.csv, TRGB_anchors.csv and SNe_TRGB.csv
    The SNe_TRGBcsv is completed it with the redshift from NED.

    :type   TRGB_SNe_Anand: str
    :param  TRGB_SNe_Anand: Path to the .csv file from Anand et al. (2021).
    :type   data_static_dir: str
    :param  data_static_dir: Path to the data_static directory
    '''
    # Values used
    mu_N4258 = 29.397
    sig_mu_N4258 = 0.032

    # Start with the TRGB DataFrame
    TRGB = pd.DataFrame()
    SN = pd.DataFrame()

    # z from NED
    z_dict = {'M66': 0.00241, 'M96': 0.00296, 'M101': 0.00080, 'N1316': 0.00601, 'N1365': 0.00546, 'N1404': 0.00649,
              'N1448': 0.00390, \
              'N4038': 0.00542, 'N4424': 0.00146, 'N4526': 0.00206, 'N4536': 0.00603, 'N5643': 0.00400,
              'N4258': 0.00149}

    if (TRGB_SNe_Anand[-4:]=='.csv'):
        tmp = pd.read_csv(TRGB_SNe_Anand, sep=',')
    elif (TRGB_SNe_Anand[-5:]=='.fits'):
        with fits.open(TRGB_SNe_Anand) as hdul:
            tmp = Table(hdul[1].data)
    else:
        print(f'ERROR: The file {TRGB_SNe_Anand} must be either a .csv or .fits')
        return 0
    TRGB['Gal'] = tmp['Gal']
    TRGB['m'] = tmp['m']
    TRGB['sig_m'] = tmp['sig_m']
    TRGB['A'] = tmp['A']
    TRGB['V-I'] = tmp['V-I']
    TRGB = TRGB.drop_duplicates(subset=['Gal']).reset_index(drop=True)
    for i in TRGB.index:
        TRGB.loc[i, 'z'] = z_dict[TRGB.loc[i, 'Gal']]
    TRGB_anchors = TRGB[TRGB['Gal'] == 'N4258'].reset_index(drop=True)
    TRGB_anchors['mu'] = mu_N4258
    TRGB_anchors['sig_mu'] = sig_mu_N4258
    TRGB = TRGB[~(TRGB['Gal'] == 'N4258')].reset_index(drop=True)

    # And the SNe
    SN['Gal'] = tmp['Gal']
    SN['mB'] = tmp['m_b']
    SN['sig_mB'] = tmp['sig_m_b']
    SN = SN.drop(SN.index[len(SN) - 1])


    # Save everything
    TRGB.to_csv(data_static_dir+'TRGB.csv',index=False)
    TRGB_anchors.to_csv(data_static_dir+'TRGB_anchors.csv', index=False)
    SN.to_csv(data_static_dir + 'SNe_TRGB.csv', index=False)
    return

def Pantheon_to_data(SNe_pantheon, data_static_dir='../data_static/'):
    '''

    :type   SNe_pantheon: str
    :param  SNe_pantheon: Path to the .txt file from Scolnic et al. (2021).
    :type   data_static_dir: str
    :param  data_static_dir: Path to the data_static directory
    '''
    # Create the DataFrame
    SN = pd.DataFrame()
    tmp = pd.read_csv(SNe_pantheon, sep=' ')
    SN['name'] = tmp['#name']
    SN['mB'] = tmp['mb']
    SN['sig_mB'] = tmp['dmb']
    SN['z'] = tmp['zcmb']
    SN['sig_z'] = tmp['dz']

    # Save everything
    SN.to_csv(data_static_dir + 'SNe_Hubble.csv', index=False)
    return