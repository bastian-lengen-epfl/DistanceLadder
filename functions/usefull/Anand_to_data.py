"""
This script converts the data from Anand et al. (2021) to the desired format.
Table 2   -> `/data/TRGB.csv` AND `/data/TRGB_anchors.csv` AND `/data/SNe_TRGB.csv`
The SNe DataFrame is completed with the redshift from NED.
"""
import os
import argparse as ap
import pandas as pd
from astropy.io import fits
from astropy.table import Table

def main(TRGB_SNe_Anand, work_dir='./'):
    data_dir = work_dir + 'data/'
    if not os.path.exists(data_dir):
        print(f'I will create the {data_dir} directory for you.')
        os.mkdir(data_dir)

    # Values used
    mu_N4258 = 29.397
    sig_mu_N4258 = 0.032

    # Start with the TRGB DataFrame
    TRGB = pd.DataFrame()
    SN = pd.DataFrame()

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
    print('Start writting the .csv files...')
    TRGB.to_csv(data_dir+'TRGB.csv',index=False)
    TRGB_anchors.to_csv(data_dir+'TRGB_anchors.csv', index=False)
    SN.to_csv(data_dir + 'SNe_TRGB.csv', index=False)
    print('Sucess')




if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Reformat the Riess et al. (2016) data to the desired format.",
                               formatter_class=ap.RawTextHelpFormatter)
    help_TRGB_SNe_Anand = 'CSV file containing the TRGB and SNe data (Table 2 Anand 2021)'
    help_work_dir = "Name of the work directory"
    parser.add_argument(dest='TRGB_SNe_Anand', type=str,
                        metavar='TRGB_SNe_Anand', action='store',
                        help=help_TRGB_SNe_Anand)
    parser.add_argument('--dir', dest='work_dir', type=str,
                        metavar='', action='store', default='./',
                        help=help_work_dir)
    args = parser.parse_args()
    main(args.TRGB_SNe_Anand, work_dir=args.work_dir)
