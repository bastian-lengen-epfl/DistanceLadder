"""
This script converts the data from Riess et al. (2016,2019,2021) to the desired format.
R16 CDS Table 4   -> `/data/Cepheids.csv` AND `/data/Cepheids_anchors.csv`
R16 Table 5     -> /data/SNe_Cepheids.csv. + complete it with the redshift from Table 1 Anderson (2019)
R19 LMC         -> added to `/data/Cepheids_anchors.csv`
R21 MW          -> `/data/Cepheids_MW.csv`
"""
import os
import argparse as ap
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

def main(Cepheids_R16, SNe_R16, Cepheids_LMC_R19, Cepheids_MW_R21, work_dir='./'):
    data_dir = work_dir + 'data/'
    if not os.path.exists(data_dir):
        print(f'I will create the {data_dir} directory for you.')
        os.mkdir(data_dir)

    # Values used
    R = 0.386               # Wesenheit
    Z_sun = 8.824           # Z_sun allows to convert 12+log(O/H) from R16 to [O/H]
    mu_N4258 = 29.397
    sig_mu_N4258 = 0.032
    mu_LMC = 18.477
    sig_mu_LMC = 0.0263
    aB_R16 = 0.71273
    sig_aB_R16 = 0.00176

    # Redshift list Anderson 2019
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
    print('Start writting the .csv files...')
    Cepheids.to_csv(data_dir+'Cepheids.csv',index=False)
    Cepheids_MW.to_csv(data_dir+'Cepheids_MW.csv',index=False)
    Cepheids_anchors.to_csv(data_dir+'Cepheids_anchors.csv', index=False)
    SNe_Cepheids.to_csv(data_dir + 'SNe_Cepheids.csv', index=False)
    print('Sucess')




if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Reformat the Riess et al. (2016) data to the desired format.",
                               formatter_class=ap.RawTextHelpFormatter)
    help_Cepheids_R16 = 'File (.csv/.fits) containing the Cepheids data (CDS Table 4 R16)'
    help_SNe_R16 = 'File (.csv/.fits) containing the SNe data (Table 5 R16)'
    help_Cepheids_LMC_R19 = 'File (.csv/.fits) containing the LMC Cepheids data (CDS Table 2 R19)'
    help_Cepheids_MW_R21 = 'File (.csv/.fits) containing the MW Cepheids data (Table 1 R21)'
    help_work_dir = "Name of the work directory"
    parser.add_argument(dest='Cepheids_R16', type=str,
                        metavar='Cepheids_R16', action='store',
                        help=help_Cepheids_R16)
    parser.add_argument(dest='SNe_R16', type=str,
                        metavar='SNe_R16', action='store',
                        help=help_SNe_R16)
    parser.add_argument(dest='Cepheids_LMC_R19', type=str,
                        metavar='Cepheids_LMC_R19', action='store',
                        help=help_Cepheids_LMC_R19)
    parser.add_argument(dest='Cepheids_MW_R21', type=str,
                        metavar='Cepheids_MW_R21', action='store',
                        help=help_Cepheids_MW_R21)
    parser.add_argument('--dir', dest='work_dir', type=str,
                        metavar='', action='store', default='./',
                        help=help_work_dir)
    args = parser.parse_args()
    main(args.Cepheids_R16, args.SNe_R16, args.Cepheids_LMC_R19, args.Cepheids_MW_R21, work_dir=args.work_dir)
