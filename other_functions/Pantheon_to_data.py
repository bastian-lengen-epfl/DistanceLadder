"""
This script converts the data from the Pantheon to the desired format.
Pantheon   -> `/data/SNe_Hubble.csv`
"""
import os
import argparse as ap
import pandas as pd


def main(SNe_pantheon, work_dir='./'):
    data_dir = work_dir + 'data/'
    if not os.path.exists(data_dir):
        print(f'I will create the {data_dir} directory for you.')
        os.mkdir(data_dir)

    SN = pd.DataFrame()
    tmp = pd.read_csv(SNe_pantheon, sep=' ')
    SN['name'] = tmp['#name']
    SN['mB'] = tmp['mb']
    SN['sig_mB'] = tmp['dmb']
    SN['z'] = tmp['zcmb']
    SN['sig_z'] = tmp['dz']


    # Save everything
    print('Start writting the .csv files...')
    SN.to_csv(data_dir + 'SNe_Hubble.csv', index=False)
    print('Sucess')




if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Reformat the Riess et al. (2016) data to the desired format.",
                               formatter_class=ap.RawTextHelpFormatter)
    help_SNe_Hubble = 'CSV file containing the SNe data (Pantheon)'
    help_work_dir = "Name of the work directory"
    parser.add_argument(dest='SNe_Hubble', type=str,
                        metavar='SNe_Hubble', action='store',
                        help=help_SNe_Hubble)
    parser.add_argument('--dir', dest='work_dir', type=str,
                        metavar='', action='store', default='./',
                        help=help_work_dir)
    args = parser.parse_args()
    main(args.SNe_Hubble, work_dir=args.work_dir)
