"""
This module contains functions that display various plots.
"""
import os
import numpy as np
import pandas as pd
import fit_parameters as fp
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import ceil

def change_plot_parameters():
    '''
    A function that can be called in order to change the rcParams from matplotlib. This function is called
    everytime a plot function is called in order to obtains an homogeneity with the plots produced.
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
    Display the PLR for each Cepheids-host galaxies. The relation in the plot is a 2D linear relation.
    To do so, the Wesenheit magnitude is corrected for metallicity, and also for the zp offset for the
    MW Cepheids.

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that was fitted.
    :type   q_dict: dictionary
    :param  q_dict: Dictionary that contains the parameters of the fits.
    :type   DF_dict_outliers: dictionary of pandas DataFrame
    :param  DF_dict_outliers: Dictionary that contains the DataFrame of the outliers. By default an empty dict.
    :type   work_dir: string
    :param  work_dir: working directory, by default ./
    '''
    ### Check and create the figure directory
    fig_dir = work_dir + 'figure/'
    if not os.path.exists(fig_dir):
        print(f'I will create the {fig_dir} directory for you.')
        os.mkdir(fig_dir)

    # Usefull function
    def plot_ij():
        '''
        Plot the individual PL relation for each galaxy
        '''
        x_lim = [-0.5, 1.1]
        # Plot errorbars
        ax[i][j].errorbar(logP-1, L, yerr=err, marker='.', ms=10, mfc="r", mec="k", ls="", c="k", lw=0.6)
        if fp.outlier_rejection == True:
            ax[i][j].plot(logP_out-1, L_out, yerr=err_out, marker='.', ms=10, mfc="lime", mec="k", ls="", c="k", lw=0.6)
        # Plot model
        if galaxy == 'MW':
            if fp.PLR_break == False:
                bs, bh = b, b
            ax[i][j].plot([x_lim[0], 0, x_lim[1]], [mW + bs * x_lim[0], mW, mW + bh * x_lim[1]], zorder=10)
        else:
            if fp.PLR_break == False:
                bs, bh = b, b
            ax[i][j].plot([x_lim[0], 0, x_lim[1]], [mW + bs * x_lim[0], mW, mW + bh * x_lim[1]] + mu, zorder=10)
        # Set scale and parameters
        ax[i][j].set_title(galaxy, fontsize=9, fontweight="bold")
        ax[i][j].set_xlabel('log(P)-1', fontsize=8)
        ax[i][j].xaxis.set_label_coords(.5, -.22)
        if ((fp.include_MW == True) and (galaxy=='MW')):
            ax[i][j].set_ylabel('$M_W$ corrected', fontsize=8)
        else:
            ax[i][j].set_ylabel('$m_W$ corrected', fontsize=8)
        ax[i][j].tick_params(axis='x', labelsize=8)
        y_lim = ax[i][j].get_ylim()
        ax[i][j].plot(np.log10([fp.break_P, fp.break_P])-1, [y_lim[0], y_lim[1]], c='k', ls='--', lw=0.5) # break line
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
    N_lines_Cep = ceil(fp.N_galaxies_Cep/5)
    N_lines_anc = ceil((fp.N_anchors_Cep + (fp.include_MW==True))/5)
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
    if fp.fixed_Zw == False:
        Zw = q_dict['Zw'][0]
    else:
        Zw = fp.Zw
    if fp.include_MW == True:
        if fp.fixed_zp == False:
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
        L = Filtered['mW'] - Zw * Filtered['M/H']
        L_out = Filtered_out['mW'] - Zw * Filtered_out['M/H']
        err = Filtered['sig_mW']
        err_out = Filtered_out['sig_mW']
        mu = q_dict[f'mu_{galaxy}'][0]
        plot_ij()
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
        L = Filtered['mW'] - Zw * Filtered['M/H']
        L_out = Filtered_out['mW'] - Zw * Filtered_out['M/H']
        err = Filtered['sig_mW']
        err_out = Filtered_out['sig_mW']
        mu = Cepheids_anchors[Cepheids_anchors['Gal'] == galaxy].iloc[0]['mu'] + q_dict[f'Dmu_{galaxy}'][0]
        plot_ij()
        i,j = incr_ij(i,j)
    # MW Cepheids
    if fp.include_MW == True:
        # Load the Cepheids_MW DF
        Cepheids_MW = DF_dict['Cepheids_MW']
        if bool(DF_dict_outliers) == True:
            Cepheids_MW_out = DF_dict_outliers['Cepheids_MW']
        else:
            Cepheids_MW_out = pd.DataFrame(columns=Cepheids_MW.columns)  # Empty DF if no outliers
        # Only 1 galaxy -> no iteration
        galaxy = 'MW'
        logP = Cepheids_MW['logP']
        logP_out = Cepheids_MW_out['logP']
        L = Cepheids_MW['mW'] - Zw * Cepheids_MW['M/H'] \
            - 10 + 5 * np.log10(Cepheids_MW['pi']) \
            + 5 * zp / np.log(10) / Cepheids_MW['pi']
        L_out = Cepheids_MW_out['mW'] - Zw * Cepheids_MW_out['M/H'] \
                - 10 + 5 * np.log10(Cepheids_MW_out['pi']) \
                + 5 * zp / np.log(10) / Cepheids_MW_out['pi']
        err = Cepheids_MW['sig_mW']
        err_out = Cepheids_MW_out['sig_mW']
        plot_ij()
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
    2D linear relation. To do so, the Wesenheit magnitude is corrected for metallicity, and also for the
    zp offset for the MW Cepheids.

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that was fitted.
    :type   q_dict: dictionary
    :param  q_dict: Dictionary that contains the parameters of the fits.
    :type   DF_dict_outliers: dictionary of pandas DataFrame
    :param  DF_dict_outliers: Dictionary that contains the DataFrame of the outliers. By default an empty dict.
    :type   work_dir: string
    :param  work_dir: working directory, by default ./
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
    if fp.fixed_Zw == False:
        Zw = q_dict['Zw'][0]
    else:
        Zw = fp.Zw
    if fp.include_MW == True:
        if fp.fixed_zp == False:
            zp = q_dict['zp'][0]
        else:
            zp = fp.zp

    ### need to shift from mW to absolute magnitude MW for each galaxy
    [logP, logP_out, M, M_out, err, err_out] = 6*[[]]
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
    # Load the Cepheids_MW DF
    Cepheids_MW = DF_dict['Cepheids_MW']
    if bool(DF_dict_outliers) == True:
        Cepheids_MW_out = DF_dict_outliers['Cepheids_MW']
    else:
        Cepheids_MW_out = pd.DataFrame(columns=Cepheids_MW.columns)  # Empty DF if no outliers
    # No need to iterate
    logP = np.append(logP, Cepheids_MW['logP'])
    logP_out = np.append(logP_out, Cepheids_MW_out['logP'])
    M = np.append(M, Cepheids_MW['mW'] - Zw * Cepheids_MW['M/H'] \
                     - 10 + 5 * np.log10(Cepheids_MW['pi']) \
                     + 5 * zp / np.log(10) / Cepheids_MW['pi'])
    M_out = np.append(M_out, Cepheids_MW_out['mW'] - Zw * Cepheids_MW_out['M/H'] \
                             - 10 + 5 * np.log10(Cepheids_MW_out['pi']) \
                             + 5 * zp / np.log(10) / Cepheids_MW_out['pi'])
    err = np.append(err, np.sqrt(Cepheids_MW['sig_mW'] ** 2 + Cepheids_MW['sig_pi'] ** 2))
    err_out = np.append(err_out, np.sqrt(Cepheids_MW_out['sig_mW'] + Cepheids_MW_out['sig_pi'] ** 2))

    ### Plot
    #  Create the figure
    fig, ax = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [7, 3]})
    fig.set_figheight(7)
    fig.set_figwidth(12)

    # Top panel
    ax[0].set_title('Global PL relation for the absolute magnitude', fontsize=16)
    ax[0].plot(logP-1, M, marker='.', ms=12, mfc="r", mec="k", ls="", c="k", lw=3)
    ax[0].plot(logP_out-1, M_out, marker='.', ms=12, mfc="lime", mec="k", ls="", c="k", lw=3)
    xmin, xmax = ax[0].get_xlim()
    ymin, ymax = ax[0].get_ylim()
    if fp.PLR_break == False:
        bs, bl = b, b
    ax[0].plot([xmin, 0, xmax], [mW + bs * xmin, mW, mW + bl * xmax], c='tab:blue', ls='-', lw=3)
    ax[0].plot([0, 0], [ymin, ymax], c='k', ls='--', lw=1.5)
    ax[0].set_xlim([xmin, xmax])
    ax[0].set_ylim([ymin, ymax])
    ax[0].invert_yaxis()
    ax[0].tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    ax[0].set_ylabel('M$_W$ corrected [mag]', fontsize=14)

    # Bottom panel
    error = np.zeros(len(logP))
    for i in range(0, np.size(logP)):
        if logP[i] < 0:
            error[i] = M[i] - mW - bs * (logP[i]-1)
        else:
            error[i] = M[i] - mW - bl * (logP[i]-1)
    error_out = np.zeros(len(logP_out))
    for i in range(0, len(logP_out)):
        if logP_out[i] < 0:
            error_out[i] = M_out[i] - mW - bs * (logP_out[i]-1)
        else:
            error_out[i] = M_out[i] - mW - bl * (logP_out[i]-1)

    ax[1].plot(logP-1, error, marker='D', ms=5, mfc="none", mec="firebrick", ls="", c="k", lw=0)
    ax[1].plot(logP_out-1, error_out, marker='D', ms=5, mfc="none", mec="lime", ls="", c="k", lw=0)
    ymin, ymax = ax[1].get_ylim()
    ax[1].plot([xmin, xmax], [0, 0], c='k', ls='--', lw=1.5)
    ax[1].plot([0, 0], [ymin, ymax], c='k', ls='--', lw=1.5)
    ax[1].set_xlim([xmin, xmax])
    ax[1].set_ylim([ymin, ymax])
    ax[1].invert_yaxis()
    ax2 = ax[1].twiny()  # ax1 and ax2 share y-axis
    ax2.plot(logP-1, error, '.', markersize=0)
    ax[1].set_xlabel('logP-1', fontsize=14)
    ax[1].set_ylabel('$\Delta$M$_W$ [mag]', fontsize=14)

    # Save
    print('Saving the PL_global plot...')
    plt.savefig(fig_dir + 'PL_global.png', bbox_inches='tight', dpi=200)
    return

def plot_SNe(DF_dict, q_dict, DF_dict_outliers)