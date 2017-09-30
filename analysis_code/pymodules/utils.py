"""Python utilities for iPython notebook analysis.

Juhye Lee."""

import os
import subprocess
import glob
import math
import time
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from colour import Color
import functools
import dms_tools2

import matplotlib
import matplotlib.lines as mlines
matplotlib.use("Pdf")
import pylab as plt
from IPython.display import Image, display

colorcycle = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

orange = '#e69f00'
skyblue = '#56b4e9'
bluishgreen = '#009e73'
yellow = '#f0e442'
blue = '#0072b2'
vermillion = '#d55e00'
reddishpurple = '#cc79a7'
darkvermillion = '#873b00'
lightvermillion = '#ff8524'

    
def MapDiffselColorToSite(diffselfile, scriptfile, script_type='pymol', 
                          map_type='abs_diffsel', colors = ['#ffff99', '#990000'], sitecolorfile=None,
                          script_preamble=False, restrict_to_chain=False, prefsfile=None):
    '''Produces a colormapping based on differential selection and writes a script for `pymol` or `chimera`
    to color a structure by this colormapping.
    
    Uses the data in *diffselfile*, which can be either a mutdiffsel or a sitediffsel file, 
    depending on the specified *map_type*. 
    
    Writes a python script to *scriptfile* for the molecular visualization program specified 
    in *script_type* (either `pymol` or `chimera`).
    
    *colors* is a list of two colors defined by hex or r,g,b codes.
    
    How the mapping from differential selection to color is determined is specified by the 
    following choices for *map_type*:
    
        * `abs_diffsel`: plot total absolute differential selection at each site. 
        Total absolute differential selection will be converted to color by interpolating between 
        the two *colors*, which will be used to show the minimum and maximum values, respectively.
        Requries *diffselfile* to be a '*sitediffsel.txt' file.
        
        * `positive_diffsel`: plot the total positive differential selection at each site. 
        Absolute positive differential selection will be converted to color by interpolating between 
        the two *colors*, which will be used to show the minimum and maximum values, respectively.
        Requries *diffselfile* to be a '*sitediffsel.txt' file.
        
        * `negative_diffsel`: plot the total negative differential selection at each site. 
        Absolute negative differential selection will be converted to color by interpolating between 
        the two *colors*, which will be used to show the (absolute) minimum and maximum values, respectively.
        Requries *diffselfile* to be a '*sitediffsel.txt' file.
        
        * `max_pos_mutdiffsel`: plot the maximum positive mutdiffsel at each site (the value for the most 
        strongly enriched mutation at each site).
        Requries *diffselfile* to be a '*mutdiffsel.txt' file.
        
        * `positive_diffsel_by_Neff`: plot the total positive differential selection at each site,
        divided by Neff, the effective number of amino-acid mutations (exp(h) where h is the Shannon entropy of
        the amino-acid preference distribution at the site). This normalizes the total amount of
        differential selection at a site by the mutational tolerance at that site. If using this option, a
        path to a valid preferences file *prefsfile* must be provided as an additional keyword argument.
        
    
    If you want to restrict the recoloring to a single chain, specify the name of the chain to color to 
    *restrict_to_chain*. This currently only works for pymol scripts, not chimera scripts.
        
    Optionally, provide a path for a *sitecolorfile* to which the (site, hexcolor, rgbcolor) mappings are
    written.
    
    Optionally, if a *script_preamble* is provided, it is written to the top of the script before the 
    commands for colormapping. This may be a useful place to add other pymol or chimera commands to load 
    the pdb file, orient the view, etc.
    '''

    # read in the data from the diffselfile:
    df = pd.read_csv(diffselfile)
    df = df.dropna()
    column_names = list(df)
    if column_names == ['site', 'wildtype', 'mutation', 'mutdiffsel']:
        filetype = 'mutdiffsel'
    elif column_names == ['site', 'abs_diffsel', 'positive_diffsel', 'negative_diffsel']:
        filetype = 'sitediffsel'
    else:
        raise ValueError('diffsel file does not have appropriate set of column identifiers')
    
    # establish the color spectrum in hex and rgb.
    n_subdivisions = 500 # the color spectrum will be divided into this many discrete colors
    color1 = Color(colors[0])
    color2 = Color(colors[1])
    hex_spectrum = [c.hex for c in color1.range_to(color2, n_subdivisions)]
    hex_spectrum_dict = dict([(i, hex_spectrum[i]) for i in range(len(hex_spectrum))]) 
    rgb_spectrum = [c.rgb for c in color1.range_to(color2, n_subdivisions)]
    rgb_spectrum_dict = dict([(i, rgb_spectrum[i]) for i in range(len(rgb_spectrum))])
    
    # generate the site ==> colorindex mapping based on *map_type*:
    if map_type == 'abs_diffsel':
        assert filetype == 'sitediffsel'
        min_diff = df.min()['abs_diffsel']  
        max_diff = df.max()['abs_diffsel']  # the min and max will be mapped to color1 and color2, respectively
        range_diff = max_diff - min_diff
        df['colorindex'] =  (df.abs_diffsel - min_diff)/range_diff*(n_subdivisions-1)
        
    elif map_type == 'negative_diffsel':
        assert filetype == 'sitediffsel'
        abs_neg_col = abs(df['negative_diffsel'])
        df = df.assign(abs_neg=abs_neg_col)
        min_diff = df.min()['abs_neg']  
        max_diff = df.max()['abs_neg']  # the min and max will be mapped to color1 and color2, respectively
        range_diff = max_diff - min_diff
        df['colorindex'] =  (df.abs_neg - min_diff)/range_diff*(n_subdivisions-1)
            
    elif map_type == 'positive_diffsel':
        assert filetype == 'sitediffsel'
        min_diff = df.min()['positive_diffsel']  
        max_diff = df.max()['positive_diffsel']  # the min and max will be mapped to color1 and color2, respectively
        range_diff = max_diff - min_diff
        df['colorindex'] =  (df.positive_diffsel - min_diff)/range_diff*(n_subdivisions-1)
    
    elif map_type == 'positive_diffsel_by_Neff':
        assert filetype == 'sitediffsel'
        # get entropies, make Neff dict, and normalize positive_diffsel values by Neff:
        (sites, wts, pi_means, pi_95credint, h) = dms_tools.file_io.ReadPreferences(prefsfile)
        neff_dict = {} # keyed by int(site) so that it can be mapped to the diffsel dataframe. 
        # this will be problematic with non-integer site strings... 
        # maybe best fix will be to cast the site column of the dataframe as string,
        # and see if the map function will work with a dictionary of strings as sites.
        for site, entropy in h.iteritems():
            neff_dict[int(site)] = 2.**entropy
        df['neff'] = df['site'].map(neff_dict)
        df['positive_diffsel_by_Neff'] = df['positive_diffsel']/df['neff']
        
        min_diff = df.min()['positive_diffsel_by_Neff']  
        max_diff = df.max()['positive_diffsel_by_Neff']  # the min and max will be mapped to color1 and color2, respectively
        range_diff = max_diff - min_diff
        df['colorindex'] =  (df.positive_diffsel_by_Neff - min_diff)/range_diff*(n_subdivisions-1)
        
        # since there are no preferences for site 1, the colorindex is currently NaN when using Neff.
        # fill to 0.
        df['colorindex'] = df['colorindex'].fillna(0)
    
    elif map_type == 'max_pos_mutdiffsel':
        assert filetype == 'mutdiffsel'
        # sort mutdiffsel file by diffsel in case it wasn't already:
        df.sort_values('mutdiffsel', ascending=False, inplace=True)
        # a new dataframe to only store the max mutdiffsel from each site
        newdf = pd.DataFrame(columns=['site', 'wt', 'mut', 'diffsel']) 
        for row in df.itertuples():
            if not any(newdf['site'].astype(int) == int(row[1])):
                newdf = newdf.append(pd.DataFrame([row[1:]], columns = ['site', 'wt', 'mut', 'diffsel']))
        newdf['site'] = newdf['site'].apply(int)
        newdf['diffsel'] = newdf['diffsel'].apply(float)
        # replace negative values with zero at sites where there was no positive selection:
        newdf['diffsel'] = (newdf['diffsel']).clip_lower(0)
        min_diff = newdf.min()['diffsel']  
        max_diff = newdf.max()['diffsel']  # the min and max will be mapped to color1 and color2, respectively
        range_diff = max_diff - min_diff
        newdf['colorindex'] =  (newdf.diffsel - min_diff)/range_diff*(n_subdivisions-1)
        df = newdf        

    else:
        raise ValueError("%s is not a recognized map_type." % map_type)
        
    # add a column for colors for each site    
    df['colorindex'] = df['colorindex'].astype(int) # round to nearest index
    df['hex'] = df['colorindex'].map(hex_spectrum_dict)
    df['rgb'] = df['colorindex'].map(rgb_spectrum_dict)        
    site_color_mapping = pd.concat([df['site'], df['hex'], df['rgb']], axis=1)
    
    if sitecolorfile:
        site_color_mapping.to_csv(sitecolorfile, index=False)
    
    # write out the script to *scriptfile*:
    f = open(scriptfile, 'w')
    
    if script_preamble:
        for line in script_preamble:
            f.write(line)
    
    if script_type == 'chimera':
        f.write("from chimera import runCommand\n")
        for i in range(len(df.index)):
            f.write("runCommand(color %s :%s)" % (site_color_mapping.iloc[i]['hex'], site_color_mapping.iloc[i]['site']) )
    elif script_type == 'pymol':
        for i in range(len(df.index)):
            rgblist = [min(1, c) for c in site_color_mapping.iloc[i]['rgb']]
            f.write("cmd.set_color(\'color%s\', \'%s\')\n" % (site_color_mapping.iloc[i]['site'], rgblist))
            if restrict_to_chain:
                f.write("cmd.color(\'color%s\', \'chain %s and resi %s\')\n" % (site_color_mapping.iloc[i]['site'], restrict_to_chain, site_color_mapping.iloc[i]['site']))
            else:
                f.write("cmd.color(\'color%s\', \'resi %s\')\n" % (site_color_mapping.iloc[i]['site'], site_color_mapping.iloc[i]['site']))
    else:
        raise ValueError("script_type must be chimera or pymol.")
    f.close()
    

def AvgAndStDevDiffsel(files, outfile, method='median', include_stdev=True):
    'method should be "median" or "mean"'
    diffsels = [pd.read_csv(f) for f in files]
    diffsel = functools.reduce(lambda left,right: pd.merge(left,right,on=['site', 'wildtype', 'mutation']), diffsels)
    diffsel_columns = ['diffsel_{0}'.format(i) for i in range(len(diffsels))]
    diffsel.columns = ['site', 'wildtype', 'mutation'] + diffsel_columns
    
    if method == 'mean':
        diffsel['mutdiffsel'] = diffsel[diffsel_columns].mean(axis=1, skipna=False)
    elif method == 'median':
        diffsel['mutdiffsel'] = diffsel[diffsel_columns].median(axis=1, skipna=False)
    else:
        raise ValueError('not a proper averaging method')
    
    if include_stdev:
        diffsel['stdev'] = diffsel[diffsel_columns].std(axis=1, skipna=False)
    
    diffsel = diffsel.sort_values('mutdiffsel', ascending=False)
    diffsel = diffsel.drop(diffsel_columns, axis=1)
    diffsel.to_csv(outfile, index=False)
    
    
def SumCodonToAA(codondict):
    aas = dms_tools2.AAS
    aadict = dict([(aa, 0) for aa in aas])
    for codon in dms_tools2.CODONS:
        aa = dms_tools2.CODON_TO_AA[codon]
        if aa in aas:
            aadict[aa] += codondict[codon]
    return aadict
    
##### NEW VERSION, compatible with dms_tools2 ######
def calculate_phi(args):
    '''Compute phi given the qRT-PCR escape fraction `gamma` and
    codon counts files `mockcounts` and `selectedcounts`.
    *args* is a dictionary with the following keys: `gamma`, `mockcounts`, `selectedcounts`, and `outfile`.
    Returns a pandas df and writes to outfile in the same format as mutdiffsel files.
    Code is heavily borrowed from the dms_diffselection script, since this is a very similar calculation.
    
    New option 3/21/17:  if a 'errorcounts' arg is provided, perform the analagous error-correction procedure that is used in differential selection.
    New option 3/22/17: if the 'gammacorrect' arg is passed as True, gamma is subtracted from phi so that this reflects fraction escape above the average mutation (testing)
    '''
    
    gamma = float(args['gamma'])
    chartype = 'codon'
    countcharacters = dms_tools2.CODONS
    characters = dms_tools2.AAS
    
    mockdf = pd.read_csv(args['mockcounts'])
    seldf = pd.read_csv(args['selectedcounts'])
    
    sites = mockdf['site'].tolist()
    assert len(sites) == len(seldf.site)
    
    wts = mockdf['wildtype'].tolist()
    assert wts == seldf['wildtype'].tolist()
    
    if args['errorcontrolcounts']:
        errdf = pd.read_csv(args['errorcontrolcounts'])
        assert len(sites) == len(errdf.site)
    
    if not args['pseudocount']:
        pseudocount = 10
        print ('using default pseudocount of 10')
    else:
        pseudocount = args['pseudocount']
    
    datacolumns = ['site', 'wildtype', 'mutation', 'mutdiffsel']
    data = dict([(column, []) for column in datacolumns])
    
    for r in range(len(sites)):
        wt = mockdf['wildtype'][r]
        mockcounts = mockdf.iloc[r, 2:]
        selcounts = seldf.iloc[r, 2:]
        
        # error control counts
        if args['errorcontrolcounts']:
            errcounts = errdf.iloc[r, 2:]
            n_err_sum = float(sum([errcounts[x] for x in countcharacters]))
            assert n_err_sum, "No counts for error control at site {0}".format(r+1)
            n_mock_sum = float(sum([mockcounts[x] for x in countcharacters]))
            assert n_mock_sum, "No counts for mock at site {0}".format(r+1)
            n_selected_sum = float(sum([selcounts[x] for x in countcharacters]))
            assert n_selected_sum, "No counts for selected at site {0}".format(r+1)
            for x in countcharacters:
                epsilon = errcounts[x] / n_err_sum
                if x == wt:
                    assert epsilon > 0, "No counts of wildtype in error control at site {0}".format(r)
                    mockcounts[x] /= epsilon
                    selcounts[x] /= epsilon
                else:
                    mockcounts[x] = max(0, n_mock_sum * (mockcounts[x] / n_mock_sum - epsilon))
                    selcounts[x] = max(0, n_selected_sum * (selcounts[x] / n_selected_sum - epsilon))
        
        # translate to AA
        mock = SumCodonToAA(mockcounts)
        selected = SumCodonToAA(selcounts)
        wt = dms_tools2.CODON_TO_AA[wt]
        
        # sum counts at site
        ntotalmock = float(sum([mock[c] for c in characters]))
        ntotalselected = float(sum([selected[c] for c in characters]))
        
        # scale pseudocounts
        if ntotalmock > ntotalselected:
            pseudocount_selected = pseudocount
            pseudocount_mock = ntotalmock / ntotalselected * pseudocount
        else:
            pseudocount_mock = pseudocount
            pseudocount_selected = ntotalselected / ntotalmock * pseudocount
        
        for x in characters:
            nxmock = float(mock[x])
            nxselected = float(selected[x])
            
            rxphi = gamma*((nxselected + pseudocount_selected) / (ntotalselected + pseudocount_selected)) / ((nxmock + pseudocount_mock) / (ntotalmock + pseudocount_mock))
            
            if args['gammacorrect']:
                rxphi -= gamma

            data['site'].append(r)
            data['wildtype'].append(wt)
            data['mutation'].append(x)
            data['mutdiffsel'].append(rxphi)
    
    phi = pd.DataFrame(data, columns=datacolumns)
    phi.to_csv(args['outfile'], index=False)
    return phi

def PlotRankedMutations(labels_files, outfileprefix, include_stdev=False, rank_lims=False, y_lims=False, title=False, colorcycle=colorcycle, alpha=0.6, make_legend=True, figsize=(5,4), ylabel='differential selection', convert_to_enrichment=False):
    '''labels_files is a list of tuples of (label, diffsel_file).
    Those diffsel_files must have a stdev column if include_stdev=True.
    
    To keep files friendly to existing dms_tools programs like dms_merge and dms_logoplot, phi files use the header `mutdiffsel` for phi (but are named appropriately)
    
    the prefix in outfile will be saved in the plots directory with .pdf added.
    rank_lims sets the x-axis limit (mutation ranks)'''
    
    fig = plt.figure(figsize=figsize)
    
    for i, (difflabel, difffile) in enumerate(labels_files):
        df = pd.read_csv(difffile).dropna()
        if not convert_to_enrichment:
            plt.plot(df['mutdiffsel'], marker='.', label = difflabel, linewidth=0.6, mew=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...
        else:
            plt.plot(2**df['mutdiffsel'], marker='.', label = difflabel, linewidth=0.6, mew=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...

    plt.xlabel('mutation rank')
    plt.ylabel(ylabel) # default is differential selection, but can change to plot rank-ordered phis.
    
    spineOffset = {'left': 4, 'bottom': 4}    
    [spine.set_position(('outward',spineOffset[loc])) if loc in ['left','bottom'] else spine.set_color('none') for loc, spine in plt.gca().spines.items() ] 
    plt.gca().tick_params(axis='x', direction='out')
    plt.gca().tick_params(axis='y', direction='out')
    plt.gca().get_xaxis().tick_bottom()
    plt.gca().get_yaxis().tick_left()
    
    if rank_lims:
        plt.xlim(rank_lims)
    if y_lims:
        plt.ylim(y_lims)
    if title:
        plt.title(title)
    if make_legend:
        legend = plt.legend(fontsize=10.5)
        
    if include_stdev:
        for i, (difflabel, difffile) in enumerate(labels_files):
            df = pd.read_csv(difffile).dropna()
            plt.gca().errorbar(range(0, len(df.index)),
                               df['mutdiffsel'],
                               yerr=df['stdev'], 
                               marker=None, color=colorcycle[i], 
                               alpha=alpha, capsize=0, elinewidth=0.9)
        
    outfile = '{0}.pdf'.format(outfileprefix)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    
def make_single_site_mutdiffselfile(infile, site, outfile):
    '''create a mutdiffsel.csv file containing just one site from the infile mutdiffsel'''
    df = pd.read_csv(infile)
    df = df[df["site"] == site]
    df.to_csv(outfile, index=False, na_rep='NaN')