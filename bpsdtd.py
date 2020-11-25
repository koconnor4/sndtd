__author__ = 'rodney'
# 2014.09.09 S.Rodney
# plotting BPS DTD curves

import numpy as np
import sys
import os
from matplotlib import pyplot as pl
import matplotlib.ticker

from pytools import plotsetup
from pytools import colorpalette as cp

from astropy.io import ascii
from astropy import table

thisfile = sys.argv[0]
if 'ipython' in thisfile : thisfile = __file__
_thisdir = os.path.dirname( os.path.abspath( thisfile ) )

DDmid = ascii.read( os.path.join( _thisdir, 'bps_dtd_dd.mid.dat' ),
                         format='commented_header', header_start=-1,
                         data_start=0 )

SDmid = ascii.read( os.path.join( _thisdir, 'bps_dtd_sd.mid.dat' ),
                         format='commented_header', header_start=-1,
                         data_start=0 )

DDedges = ascii.read( os.path.join( _thisdir, 'bps_dtd_dd.edges.dat' ),
                         format='commented_header', header_start=-1,
                         data_start=0 )

SDedges = ascii.read( os.path.join( _thisdir, 'bps_dtd_sd.edges.dat' ),
                         format='commented_header', header_start=-1,
                         data_start=0 )


grouplist = [ 'Ruiter', 'Bours/Toonen','Wang/Han',  'Claeys', 'Mennekens',  'Yungelson'  ]
colorlist = [  cp.red,    cp.purple,   cp.darkblue,   cp.gold,  cp.green,    cp.lightblue ]
reflist =   [ 'Ruiter, Belcynski \& Fryer 2009', 'Bours, Toonen \& Nelemans 2013','Wang \& Han 2012', 'Claeys et al. 2014','Mennekens et al. 2010','Yungelson 2010' ]

def plot_for_chicago( modelclass='DD', style='mid'):
    """ make a figure showing BPS DTDs from 6 DD models, plus t-1
    :return:
    """
    if modelclass=='DD' and style=='mid' : model = DDmid
    elif modelclass=='DD' and style=='edges' : model = DDedges
    elif modelclass=='SD' and style=='mid' : model = SDmid
    elif modelclass=='SD' and style=='edges' : model = SDedges

    plotsetup.presfig( figsize=[15,8] )
    pl.clf()

    tMyr = model['t']
    tGyr = tMyr / 1000.
    # tGyr[1:] -= (tGyr[1:] - tGyr[:-1])/2.
    for group, color, ref in zip( grouplist, colorlist, reflist ) :
        dtd = np.where( model[group]>0, model[group], 1e-22 )
        # dtd = model[group]
        pl.plot( tGyr, dtd, marker=' ',ls='-', color=color, label=ref )

    tGyr = np.arange( 0.04, 13.7, 0.1 )
    tminus1top = 1.4e-13 * tGyr**-1
    tminus1bot = 0.7e-13 * tGyr**-1
    pl.fill_between( tGyr, tminus1top, tminus1bot, lw=0,
                     color='k', alpha=0.3, label='_nolegend_'  )

    pl.legend( loc='upper right', frameon=False, ncol=3, fontsize=18,
               bbox_to_anchor=[0.99,0.97] )
    ax = pl.gca()
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_ylim( 1e-17, 8e-11 )
    ax.set_xlim( 0.01, 15. )
    ax.set_xticks( [0.1, 1.0, 10.0] )
    ax.set_xticklabels( ['0.1', '1', '10'] )
    ax.set_xlabel('Delay time (Gyr)')
    ax.set_ylabel(r'SN Ia Rate ( yr$^{-1}$ M$_{\odot}^{-1}$)')

    ax.text( 0.05,3e-12, r'DTD$\propto t^{-1}$', color='0.5',
             ha='right',va='top',fontweight='heavy', fontsize=32,
             rotation=-15)


    fig = pl.gcf()
    fig.subplots_adjust( left=0.15, bottom=0.15, right=0.95, top=0.95 )
    pl.draw()








