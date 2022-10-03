'''
This script the data from the ./data_static/*/*.csv, convert them to the expected format, and save the resulting
.csv in the ./data folder.
'''
import os
import argparse as ap
from converting_functions import *


def setup(MW_dataset, data_static_dir, data_dir):
    ### Process Riess, Anand and Pantheon
    print('Converting Riess data...')
    Riess_to_data('Riess/Cepheids_R16.fits', 'Riess/SN_R16_Table5.csv', 'Riess/Cepheids_LMC_R19.fits',\
                  'Riess/Cepheids_MW_R21.csv', data_static_dir=data_static_dir, data_dir=data_dir)
    print('Converting Anand data...')
    Anand_to_data('Anand/SN_TRGB.csv', data_static_dir=data_static_dir, data_dir=data_dir)
    print('Converting Pantheon data...')
    Pantheon_to_data('Pantheon/Pantheon.txt', data_static_dir=data_static_dir, data_dir=data_dir)

    ### Makes sure user choosed R/H/RH for the dataset setup
    while MW_dataset not in ['R', 'H', 'RH']:
        print(
            'ERROR: You did not choose what dataset to use for the MW Cepheids (Riess \'R\', H1PStars \'H\', or both \'RH\').')
        MW_dataset = input('Please, enter what setup you want (R/H/RH):')

    ### Choose what to do for the .csv
    if MW_dataset=='R':
        pass
    elif MW_dataset=='H':
        print('Converting H1PStars data (erase R21 MW Cepheids)...')
        H1PStars_to_data('H1PStars/Cepheids_MW.csv', erase=True, data_static_dir=data_static_dir,
                         data_dir=data_dir)
    elif MW_dataset=='RH':
        print('Converting H1PStars data (add to R21 MW Cepheids)...')
        H1PStars_to_data('H1PStars/Cepheids_MW.csv', erase=False, data_static_dir=data_static_dir,
                         data_dir=data_dir)
    else:
        print('ERROR: Could not find which setup you wanted (R/H/RH)')
    print('Done.')
    return

if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Convert the data from the data_static/ folder to the desired format,\
                                            and then save the new .csv in the ./data.",
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument(type=str, dest='MW_dataset', metavar='MW_dataset',
                        nargs='?', const=' ',
                        help='Choose the dataset for the MW Cepheids: Riess \'R\', H1PStars \'H\', or both Riess \'RH\'')
    parser.add_argument('--dir', type=str, dest='data_static_dir', metavar='data_static_dir',
                        default='./data_static/',
                        help='Name of the data_static directory')
    parser.add_argument('--dirdata', type=str, dest='data_dir', metavar='data_dir',
                        default='./data/',
                        help='Name of the data directory')

    args = parser.parse_args()
    setup(args.MW_dataset, args.data_static_dir, args.data_dir)

