'''
This module contains functions that display various plots.
'''
import os
import numpy as np
import pandas as pd
import fit_parameters as fp
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import ceil

def change_plot_parameters():
    '''
    A function that can be called in order to change the rcParams from matplotlib.
    '''
    plt.style.use('classic')
    mpl.rcParams['axes.titlesize'] = 24
    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['lines.markersize'] = 10
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['figure.figsize'] = (8, 5)
    return

def plot_individual_PL(DF_dict, q_dict, DF_dict_outliers, work_dir='./'):
    '''
    Display the PLR for each Cepheid-host galaxies. The relation in the plot is a 2D linear relation.
    Therefore, the Wesenheit magnitude is corrected for metallicity, and for the zp offset for the
    MW Cepheids.

    Parameters
    ----------
    DF_dict : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have been fitted.
    q_dict : dict of numpy array
        A dictionary that contains the fit parameters and their uncertainty.
    DF_dict_outlier : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have not been fitted.
    work_dir : str
        working directory, by default ./
    '''

    ### Check and create the figure directory
    fig_dir = work_dir + 'figure/'
    if not os.path.exists(fig_dir):
        print(f'I will create the {fig_dir} directory for you.')
        os.mkdir(fig_dir)

    # Usefull function
    def plot_ij(is_anchor):
        '''
        Plot the individual PL relation for each galaxy

        Parameters
        ----------
        is_anchor : bool
            True for anchor galaxy.
        '''

        x_lim = [0.5, 2.1]
        # Plot errorbars
        marker = 'v' if is_anchor==True else 'o'
        ax[i][j].errorbar(logP, L, marker=marker, ms=6, mfc="r", mec="k", ls="", c="k", lw=0.6)
        if fp.outlier_rejection == True:
            ax[i][j].plot(logP_out, L_out, marker=marker, ms=6, mfc="lime", mec="k", ls="", c="k", lw=0.6)
        # Plot model
        if galaxy[:2] == 'MW':
            if fp.PLR_break == False:
                bs, bl = b, b
            else:
                bs = q_dict['bs'][0]
                bl = q_dict['bl'][0]
            ax[i][j].plot([x_lim[0], np.log10(fp.break_P), x_lim[1]], mW +
                          [bs * (x_lim[0]-np.log10(fp.break_P)), 0, bl * (x_lim[1]-np.log10(fp.break_P))], zorder=10)
        else:
            if fp.PLR_break == False:
                bs, bl = b, b
            else:
                bs = q_dict['bs'][0]
                bl = q_dict['bl'][0]
            ax[i][j].plot([x_lim[0], np.log10(fp.break_P), x_lim[1]], mW + mu +
                          [bs*(x_lim[0]-np.log10(fp.break_P)), 0, bl*(x_lim[1]-np.log10(fp.break_P))], zorder=10)
        if (is_anchor == False) and (fp.PLR_break2 == True):
            bL = q_dict['bL'][0]
            ax[i][j].plot([np.log10(fp.break_P2), x_lim[1]], mW + mu + bl * np.log10(fp.break_P2/fp.break_P) +
                          [0, bL * (x_lim[1]-np.log10(fp.break_P2))], color='midnightblue', ls='--', zorder=11)
        # Set scale and parameters
        ax[i][j].set_title(galaxy, fontsize=9, fontweight="bold")
        ax[i][j].set_xlabel('log(P)', fontsize=8)
        ax[i][j].xaxis.set_label_coords(.5, -.22)
        if ((fp.include_MW == True) and (galaxy[:2] == 'MW')):
            ax[i][j].set_ylabel('$M_W-Z_W$[Fe/H]', fontsize=8)
        else:
            ax[i][j].set_ylabel('$m_W-Z_W$[Fe/H]', fontsize=8)
        ax[i][j].tick_params(axis='x', labelsize=8)
        y_lim = ax[i][j].get_ylim()
        ax[i][j].plot(np.log10([fp.break_P, fp.break_P]), [y_lim[0], y_lim[1]], c='k', ls='--', lw=0.5) # break line
        if (fp.PLR_break2 == True) and (is_anchor == False):
            ax[i][j].plot(np.log10([fp.break_P2, fp.break_P2]), [y_lim[0], y_lim[1]], c='k', ls='--',
                          lw=0.5)  # break line
        ax[i][j].set_xlim(x_lim)
        ax[i][j].set_ylim(y_lim)
        ax[i][j].invert_yaxis()
        ax[i][j].tick_params(axis='y', labelsize=8)
        return
    # increment ij
    def incr_ij(i, j):
        '''
        Increment i and j in order to match the subplots order 
        '''
        j += 1
        if j == 5:
            i += 1
            j = 0
        return i,j

    ### Create the Nx5 grid of subplot
    if fp.include_MW == True:
        N_MW_dataset = len(DF_dict['Cepheids_MW']['Gal'].drop_duplicates().reset_index(drop=True))
    else:
        N_MW_dataset = 0
    N_lines_Cep = ceil(fp.N_galaxies_Cep/5)
    N_lines_anc = ceil((fp.N_anchors_Cep + N_MW_dataset)/5)
    N = N_lines_Cep + N_lines_anc
    fig, ax = plt.subplots(nrows=N, ncols=5)
    fig.tight_layout()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.6)
    fig.set_figheight(7)
    fig.set_figwidth(12)

    ### Load the fit parameters
    mW = q_dict['MW'][0]
    if fp.PLR_break == True:
        bs = q_dict['bs'][0]
        bl = q_dict['bl'][0]
    else:
        b = q_dict['b'][0]
    if fp.PLR_break2 == True:
        bL = q_dict['bL']
    if fp.fixed_Zw == False:
        Zw = q_dict['Zw'][0]
    else:
        Zw = fp.Zw
    if fp.include_MW == True:
        if fp.fixed_zp == False:
            if fp.multiple_zp == True:
                zp_list = []
                for i in range(N_MW_dataset):
                    zp_list.append(q_dict[f'zp{i+1}'][0])
            else:
                zp = q_dict['zp'][0]
        else:
            zp = fp.zp

    ### Individual plot for each galaxy:
    i, j = 0, 0
    # Load the Cepheids DF
    Cepheids = DF_dict['Cepheids']
    if bool(DF_dict_outliers) == True:
        Cepheids_out = DF_dict_outliers['Cepheids']
    else:
        Cepheids_out = pd.DataFrame(columns=Cepheids.columns) #Empty DF if no outliers
    # Iterate over each galaxy for the Cepheids DF
    galaxies_Cep = Cepheids['Gal'].drop_duplicates().reset_index(drop=True)
    for galaxy in galaxies_Cep:
        Filtered = Cepheids[Cepheids['Gal']==galaxy]
        Filtered_out = Cepheids_out[Cepheids_out['Gal']==galaxy]
        logP = Filtered['logP']
        logP_out = Filtered_out['logP']
        L = Filtered['mW'] - Zw * Filtered['M/H'] # Correct for metallicity
        L_out = Filtered_out['mW'] - Zw * Filtered_out['M/H'] # Correct for metallicity
        err = Filtered['sig_mW']
        err_out = Filtered_out['sig_mW']
        mu = q_dict[f'mu_{galaxy}'][0]
        plot_ij(is_anchor=False)
        i,j = incr_ij(i,j)
        
    # Remove empty subplots
    while j!=0:
        ax[i][j].remove()
        i,j = incr_ij(i,j)
    # Load the Cepheids_anchors DF
    Cepheids_anchors = DF_dict['Cepheids_anchors']
    if bool(DF_dict_outliers) == True:
        Cepheids_anchors_out = DF_dict_outliers['Cepheids_anchors']
    else:
        Cepheids_anchors_out = pd.DataFrame(columns=Cepheids_anchors.columns)  # Empty DF if no outliers
    # Iterate over each galaxy for the Cepheids_anchors DF
    anchors = Cepheids_anchors['Gal'].drop_duplicates().reset_index(drop=True)
    for galaxy in anchors:
        Filtered = Cepheids_anchors[Cepheids_anchors['Gal'] == galaxy]
        Filtered_out = Cepheids_anchors_out[Cepheids_anchors_out['Gal'] == galaxy]
        logP = Filtered['logP']
        logP_out = Filtered_out['logP']
        L = Filtered['mW'] - Zw * Filtered['M/H'] # Correct for metallicity
        L_out = Filtered_out['mW'] - Zw * Filtered_out['M/H'] # Correct for metallicity
        err = Filtered['sig_mW']
        err_out = Filtered_out['sig_mW']
        mu = Cepheids_anchors[Cepheids_anchors['Gal'] == galaxy].iloc[0]['mu'] + q_dict[f'Dmu_{galaxy}'][0]
        plot_ij(is_anchor=True)
        i,j = incr_ij(i,j)
    # MW Cepheids
    if fp.include_MW == True:
        # Load the Cepheids_MW DF
        Cepheids_MW = DF_dict['Cepheids_MW']
        list_MW = Cepheids_MW['Gal'].drop_duplicates().reset_index(drop=True)
        if bool(DF_dict_outliers) == True:
            Cepheids_MW_out = DF_dict_outliers['Cepheids_MW']
        else:
            Cepheids_MW_out = pd.DataFrame(columns=Cepheids_MW.columns)  # Empty DF if no outliers
        for galaxy in list_MW:
            Filtered = Cepheids_MW[Cepheids_MW['Gal'] == galaxy]
            Filtered_out = Cepheids_MW_out[Cepheids_MW_out['Gal'] == galaxy]
            logP = Filtered['logP']
            logP_out = Filtered_out['logP']
            if fp.multiple_zp == True:
                zp_index = list_MW[list_MW == galaxy].index[0]
                zp_tmp = zp_list[zp_index]
            else:
                zp_tmp = zp
            L = Filtered['mW'] - Zw * Filtered['M/H'] \
                - 10 + 5 * np.log10(Filtered['pi']) \
                + 5 * zp_tmp / np.log(10) / Filtered['pi'] # Correct for metallicity and zp and distance
            L_out = Filtered_out['mW'] - Zw * Filtered_out['M/H'] \
                    - 10 + 5 * np.log10(Filtered_out['pi']) \
                    + 5 * zp_tmp / np.log(10) / Filtered_out['pi'] # Correct for metallicity and zp and distance
            err = Filtered['sig_mW']
            err_out = Filtered_out['sig_mW']
            plot_ij(is_anchor=True)
            i,j = incr_ij(i,j)
    # Remove empty subplots
    while j!=0:
        ax[i][j].remove()
        i,j = incr_ij(i,j)
    # Save
    print('Saving the PL_individual plot...')
    plt.savefig(fig_dir+'PL_individual.png',bbox_inches='tight', dpi=200)
    return

def plot_global_PL(DF_dict, q_dict, DF_dict_outliers, work_dir='./'):
    '''
    Display the global PLR (absolute magnitude) for each all the Cepheids. The relation in the plot is a
    2D linear relation. Therefore, the Wesenheit magnitude is corrected for metallicity, distance, and also for the
    zp offset for the MW Cepheids.

    Parameters
    ----------
    DF_dict : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have been fitted.
    q_dict : dict of numpy array
        A dictionary that contains the fit parameters and their uncertainty.
    DF_dict_outlier : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have not been fitted.
    work_dir : str
        working directory, by default ./
    '''
    ### Check and create the figure directory
    fig_dir = work_dir + 'figure/'
    if not os.path.exists(fig_dir):
        print(f'I will create the {fig_dir} directory for you.')
        os.mkdir(fig_dir)

    ### Load the fit parameters
    mW = q_dict['MW'][0]
    if fp.PLR_break == True:
        bs = q_dict['bs'][0]
        bl = q_dict['bl'][0]
    else:
        b = q_dict['b'][0]
    if fp.PLR_break2 == True:
        bL = q_dict['bL'][0]
    if fp.fixed_Zw == False:
        Zw = q_dict['Zw'][0]
    else:
        Zw = fp.Zw
    if fp.include_MW == True:
        if fp.fixed_zp == False:
            if fp.multiple_zp == True:
                zp_list = []
                for i in range(len(DF_dict['Cepheids_MW']['Gal'].drop_duplicates().reset_index(drop=True))):
                    zp_list.append(q_dict[f'zp{i + 1}'][0])
            else:
                zp = q_dict['zp'][0]
        else:
            zp = fp.zp

    ### need to shift from mW to absolute magnitude MW for each galaxy
    [logP, logP_out, M, M_out, err, err_out, anchor_flag, anchor_flag_out] = 8*[[]]
    # Load the Cepheids DF
    Cepheids = DF_dict['Cepheids']
    if bool(DF_dict_outliers) == True:
        Cepheids_out = DF_dict_outliers['Cepheids']
    else:
        Cepheids_out = pd.DataFrame(columns=Cepheids.columns)  # Empty DF if no outliers
    # Iterate over each galaxy for the Cepheids DF
    galaxies_Cep = Cepheids['Gal'].drop_duplicates().reset_index(drop=True)
    for galaxy in galaxies_Cep:
        [mu, dmu] = q_dict[f'mu_{galaxy}']
        Filtered = Cepheids[Cepheids['Gal'] == galaxy]
        Filtered_out = Cepheids_out[Cepheids_out['Gal'] == galaxy]
        logP = np.append(logP, Filtered['logP'])
        logP_out = np.append(logP_out, Filtered_out['logP'])
        M = np.append(M, Filtered['mW'] - Zw * Filtered['M/H'] - mu)
        M_out = np.append(M_out, Filtered_out['mW'] - Zw * Filtered_out['M/H'] - mu)
        err = np.append(err,np.sqrt(Filtered['sig_mW']**2 + dmu**2))
        err_out = np.append(err_out,np.sqrt(Filtered_out['sig_mW']+dmu**2))
        anchor_flag = np.append(anchor_flag, np.zeros(len(Filtered)))
        anchor_flag_out = np.append(anchor_flag_out, np.zeros(len(Filtered_out)))
    # Load the Cepheids_anchors DF
    Cepheids_anchors = DF_dict['Cepheids_anchors']
    if bool(DF_dict_outliers) == True:
        Cepheids_anchors_out = DF_dict_outliers['Cepheids_anchors']
    else:
        Cepheids_anchors_out = pd.DataFrame(columns=Cepheids_anchors.columns)  # Empty DF if no outliers
    # Iterate over each galaxy for the Cepheids DF
    anchors = Cepheids_anchors['Gal'].drop_duplicates().reset_index(drop=True)
    for galaxy in anchors:
        Filtered = Cepheids_anchors[Cepheids_anchors['Gal'] == galaxy]
        Filtered_out = Cepheids_anchors_out[Cepheids_anchors_out['Gal'] == galaxy]
        logP = np.append(logP, Filtered['logP'])
        logP_out = np.append(logP_out, Filtered_out['logP'])
        M = np.append(M, Filtered['mW'] - Zw * Filtered['M/H'] - Filtered['mu'])
        M_out = np.append(M_out, Filtered_out['mW'] - Zw * Filtered_out['M/H'] - Filtered_out['mu'])
        err = np.append(err, np.sqrt(Filtered['sig_mW'] ** 2 + Filtered['sig_mu'] ** 2))
        err_out = np.append(err_out, np.sqrt(Filtered_out['sig_mW'] + Filtered_out['sig_mu'] ** 2))
        anchor_flag = np.append(anchor_flag, np.ones(len(Filtered)))
        anchor_flag_out = np.append(anchor_flag_out, np.ones(len(Filtered_out)))
    # Load the Cepheids_MW DF
    if fp.include_MW == True:
        Cepheids_MW = DF_dict['Cepheids_MW']
        list_MW = Cepheids_MW['Gal'].drop_duplicates().reset_index(drop=True)
        if bool(DF_dict_outliers) == True:
            Cepheids_MW_out = DF_dict_outliers['Cepheids_MW']
        else:
            Cepheids_MW_out = pd.DataFrame(columns=Cepheids_MW.columns)  # Empty DF if no outliers
        for galaxy in list_MW:
            Filtered = Cepheids_MW[Cepheids_MW['Gal'] == galaxy]
            Filtered_out = Cepheids_MW_out[Cepheids_MW_out['Gal'] == galaxy]
            logP = np.append(logP, Filtered['logP'])
            logP_out = np.append(logP_out, Filtered_out['logP'])
            if fp.multiple_zp == True:
                zp_index = list_MW[list_MW == galaxy].index[0]
                zp_tmp = zp_list[zp_index]
            else:
                zp_tmp = zp
            M = np.append(M, Filtered['mW'] - Zw * Filtered['M/H'] \
                         - 10 + 5 * np.log10(Filtered['pi']) \
                         + 5 * zp_tmp / np.log(10) / Filtered['pi'])
            M_out = np.append(M_out, Filtered_out['mW'] - Zw * Filtered_out['M/H'] \
                                 - 10 + 5 * np.log10(Filtered_out['pi']) \
                                 + 5 * zp_tmp / np.log(10) / Filtered_out['pi'])
            err = np.append(err, np.sqrt(Filtered['sig_mW'] ** 2 + Filtered['sig_pi'] ** 2))
            err_out = np.append(err_out, np.sqrt(Filtered_out['sig_mW'] + Filtered_out['sig_pi'] ** 2))
            anchor_flag = np.append(anchor_flag, np.ones(len(Filtered)))
            anchor_flag_out = np.append(anchor_flag_out, np.ones(len(Filtered_out)))

    ### Plot
    #  Create the figure
    fig, ax = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [7, 3]})
    fig.set_figheight(7)
    fig.set_figwidth(12)

    # Top panel
    ax[0].set_title('Global PL relation for the absolute magnitude', fontsize=16)
    ax[0].plot(logP[anchor_flag==0], M[anchor_flag==0], marker='o', ms=6, mfc="r", mec="k", ls="", c="k", lw=3)
    ax[0].plot(logP[anchor_flag==1], M[anchor_flag==1], marker='v', ms=6, mfc="r", mec="k", ls="", c="k", lw=3)
    ax[0].plot(logP_out[anchor_flag_out==0], M_out[anchor_flag_out==0], marker='o', ms=6, mfc="lime", mec="k",
                   ls="", c="k", lw=3)
    ax[0].plot(logP_out[anchor_flag_out==1], M_out[anchor_flag_out==1], marker='v', ms=6, mfc="lime", mec="k",
                   ls="", c="k", lw=3)
    ax[0].legend(['SN-hosts', 'Anchors', 'SN-hosts + outlier', 'Anchors + outlier'])
    xmin, xmax = ax[0].get_xlim()
    ymin, ymax = ax[0].get_ylim()
    if fp.PLR_break == False:
        bs, bl = b, b
    ax[0].plot([xmin, np.log10(fp.break_P), xmax], mW +
               [bs * (xmin-np.log10(fp.break_P)), 0, bl * (xmax-np.log10(fp.break_P))], c='tab:blue', ls='-', lw=3)
    if fp.PLR_break2 == True:
        ax[0].plot([np.log10(fp.break_P2), xmax], mW + bl*np.log10(fp.break_P2/fp.break_P)+
                   [0, bL*(xmax-np.log10(fp.break_P2))], c='midnightblue', ls='--', lw=3)
        ax[0].plot(2 * [np.log10(fp.break_P2)], [ymin, ymax], c='k', ls='--', lw=1.5)
    ax[0].plot(2*[np.log10(fp.break_P)], [ymin, ymax], c='k', ls='--', lw=1.5)
    ax[0].set_xlim([xmin, xmax])
    ax[0].set_ylim([ymin, ymax])
    ax[0].invert_yaxis()
    ax[0].tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    ax[0].set_ylabel('$M_W-Z_W$[Fe/H] [mag]', fontsize=14)

    # Bottom panel
    error = np.zeros(len(logP))
    for i in range(0, np.size(logP)):
        if logP[i] < np.log10(fp.break_P):
            error[i] = M[i] - mW - bs * (logP[i]-np.log10(fp.break_P))
        else:
            if fp.PLR_break2 == True:
                error[i] = M[i] - mW - bl * (np.log10(fp.break_P2/fp.break_P)) - bL*(logP[i] - np.log10(fp.break_P2))
            else:
                error[i] = M[i] - mW - bl * (logP[i] - np.log10(fp.break_P))
    error_out = np.zeros(len(logP_out))
    for i in range(0, len(logP_out)):
        if logP_out[i] < 0:
            error_out[i] = M_out[i] - mW - bs * (logP_out[i]-np.log10(fp.break_P))
        else:
            if fp.PLR_break2 == True:
                error_out[i] = M_out[i] - mW - bl * (np.log10(fp.break_P2 / fp.break_P)) - bL * (
                               logP_out[i] - np.log10(fp.break_P2))
            else:
                error_out[i] = M_out[i] - mW - bl * (logP_out[i] - np.log10(fp.break_P))

    ax[1].plot(logP[anchor_flag==0], error[anchor_flag==0], marker='o', ms=6, mfc="none", mec="firebrick", ls="", c="k", lw=0)
    ax[1].plot(logP[anchor_flag==1], error[anchor_flag==1], marker='v', ms=6, mfc="none", mec="firebrick", ls="", c="k", lw=0)
    ax[1].plot(logP_out[anchor_flag_out==0], error_out[anchor_flag_out==0], marker='o', ms=6, mfc="none", mec="lime", ls="", c="k", lw=0)
    ax[1].plot(logP_out[anchor_flag_out==1], error_out[anchor_flag_out==1], marker='v', ms=6, mfc="none", mec="lime", ls="", c="k", lw=0)
    ymin, ymax = ax[1].get_ylim()
    ax[1].legend(['SN-hosts', 'Anchors', 'SN-hosts + outlier', 'Anchors + outlier'])
    ax[1].plot([xmin, xmax], [0, 0], c='k', ls='--', lw=1.5)
    ax[1].plot([0, 0], [ymin, ymax], c='k', ls='--', lw=1.5)
    ax[1].set_xlim([xmin, xmax])
    ax[1].set_ylim([ymin, ymax])
    ax[1].invert_yaxis()
    ax2 = ax[1].twiny()  # ax1 and ax2 share y-axis
    ax2.plot(logP, error, '.', markersize=0)
    ax[1].set_xlabel('logP', fontsize=14)
    ax[1].set_ylabel('$\Delta$M$_W$ [mag]', fontsize=14)

    # Save
    print('Saving the PL_global plot...')
    plt.savefig(fig_dir + 'PL_global.png', bbox_inches='tight', dpi=200)
    return

def plot_SNe(DF_dict, q_dict, DF_dict_outliers, work_dir='./'):
    '''
    Display the global redshift-magnitude plot for the SNe_Hubble DataFrame.

    Parameters
    ----------
    DF_dict : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have been fitted.
    q_dict : dict of numpy array
        A dictionary that contains the fit parameters and their uncertainty.
    DF_dict_outlier : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have not been fitted.
    work_dir : str
        working directory, by default ./
    '''
    ### Useful functions
    def logczexpansion(z):
        return np.log10(
            fp.c * z * (1 + 1 / 2 * (1 - fp.q0) * z - 1 / 6 * (1 - fp.q0 - 3 * fp.q0 ** 2 + fp.j0) * z ** 2))

    ### Check and create the figure directory
    fig_dir = work_dir + 'figure/'
    if not os.path.exists(fig_dir):
        print(f'I will create the {fig_dir} directory for you.')
        os.mkdir(fig_dir)

    ### Load the fit parameters
    if fp.fit_aB==False:
        aB = fp.aB
    else:
        aB = (q_dict['5logH0'][0] - q_dict['MB'][0] - 25)/5

    ### Load the SNe DF
    SNe = DF_dict['SNe_Hubble']
    if bool(DF_dict_outliers) == True:
        SNe_out = DF_dict_outliers['SNe_Hubble']
    else:
        SNe_out = pd.DataFrame(columns=SNe.columns)  # Empty DF if no outliers

    filter = ((SNe['z']<fp.z_max)&(SNe['z']>fp.z_min))
    SNe_out = pd.concat([SNe_out, SNe[~filter]], ignore_index=True)
    SNe = SNe[filter].reset_index(drop=True)

    ### Plot
    # Create the figure
    fig, ax = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [3, 1]})
    fig.set_figheight(7)
    fig.set_figwidth(12)
    fig.subplots_adjust(wspace=0, hspace=0)

    # Top panel
    ax[0].plot(logczexpansion(SNe_out['z']), 0.2 * SNe_out['mB'], marker='.', ms=12, mfc="lime", mec="k", ls="", c="k", lw=3)
    ax[0].plot(logczexpansion(SNe['z']), 0.2 * SNe['mB'], marker='.', ms=12, mfc="r", mec="k", ls="", c="k", lw=3)
    tmp = ax[0].get_xlim()
    ax[0].plot(tmp, np.array(tmp) - aB, 'k', lw=2)  # slope 1 by mean
    ax[0].set_xlim(tmp)
    tmp = ax[0].get_ylim()
    ax[0].plot(logczexpansion(np.array([fp.z_min, fp.z_min])), tmp, c='k', ls='--', lw=1)
    ax[0].text(logczexpansion(fp.z_min) + 0.05, tmp[1] - 0.15, f'z={fp.z_min}', size=16)
    ax[0].plot(logczexpansion(np.array([fp.z_max, fp.z_max])), tmp, c='k', ls='--', lw=1)
    ax[0].text(logczexpansion(fp.z_max) + 0.05, tmp[1] - 0.15, f'z={fp.z_max}', size=16)
    ax[0].set_ylim(tmp)
    ax[0].tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    ax[0].set_ylabel('0.2m$_B$ [mag]', size=16)
    ax[0].tick_params(axis='both', which='major', labelsize=12)

    # Bottom panel
    ax[1].plot(logczexpansion(SNe_out['z']), 0.2 * SNe_out['mB'] - (logczexpansion(SNe_out['z']) - aB), \
               marker='D', ms=5, mfc="none", mec="limegreen", ls="", c="k", lw=0)
    ax[1].plot(logczexpansion(SNe['z']), 0.2 * SNe['mB'] - (logczexpansion(SNe['z']) - aB), \
               marker='D', ms=5, mfc="none", mec="r", ls="", c="k", lw=0)
    tmp = [ax[0].get_xlim()[0], ax[0].get_xlim()[1]]
    ax[1].plot(tmp, [0, 0], c='k', ls='--', lw=3)
    ax[1].set_xlim(tmp)
    tmp = ax[1].get_ylim()
    ax[1].plot(logczexpansion(np.array([fp.z_min, fp.z_min])), tmp, c='k', ls='--', lw=1)
    ax[1].plot(logczexpansion(np.array([fp.z_max, fp.z_max])), tmp, c='k', ls='--', lw=1)
    ax[1].set_ylim(tmp)
    ax[1].tick_params(axis='both', which='major', labelsize=12)
    ax2 = ax[1].twiny()  # ax1 and ax2 share y-axis
    ax2.plot(logczexpansion(SNe_out['z']), 0.2 * SNe_out['mB'] - (logczexpansion(SNe_out['z']) - aB), '.', markersize=0)
    ax[1].set_xlabel('log{cz[1+0.5(1-q$_0$)z-(1/6)(1-q$_0$-3q$_0^2$+1)z$^2$]}', size=16)
    ax[1].set_ylabel('$\Delta$0.2m$_B$ [mag]', size=16)

    # Save
    print('Saving the resdhift-magnitude plot...')
    plt.savefig(fig_dir + 'redshift_magnitude.png', bbox_inches='tight', dpi=200)
    return