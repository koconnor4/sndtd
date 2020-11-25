__author__ = 'rodney'


from matplotlib import (pyplot as pl,ticker)
import numpy as np
from pytools import (colorpalette as cp,plotsetup)

def sfr( z, z0=1.243, A=-0.997, B=0.241, C=0.180 ):
    """  Double power-law fit a la Behroozi et al 2013"""
    return( C / ( 10**(A*(z-z0)) + 10**(B*(z-z0)) )   )


def prepfig():
    from matplotlib import rcParams
    plotsetup.presfig( figsize=[4,4])
    rcParams['xtick.labelsize'] = 20
    rcParams['ytick.labelsize'] = 20
    fig = pl.gcf()
    fig.subplots_adjust(left=0.15,bottom=0.18, right=0.95,top=0.95)
    pl.clf()

def mk_sfr_fig():
    prepfig()

    z = np.arange( 0, 3, 0.01 )
    sfrlinear = 0.8* np.arange( 0, 3, 0.01 )
    ax = pl.gca()
    ax.plot( z, sfrlinear, lw=3.5, ls='-', color=cp.darkbluegray )
    ax.set_xlabel('Redshift', fontsize=24 )
    ax.set_ylabel('Star Formation Rate', fontsize=24 )


    ax.set_ylim( 0, 2 )
    ax.set_xlim( 0, 2 )
    # ax.set_xticklabels( [] )
    ax.set_yticklabels( [] )
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 1.0 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1.0 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    pl.draw()

def mk_dtd_fig():
    prepfig()

    t = np.arange( 0.04, 10, 0.01 )
    dtdgauss = 1.4 * np.exp( -( (t-3)**2 / 0.5**2 )) # (1./(2*np.pi*0.5) ) *
    dtdt1 = 2. / (25.*t)
    ax = pl.gca()
    ax.plot( t, dtdgauss, lw=3.5, ls='-', color=cp.darkgreen )
    ax.plot( t, dtdt1, lw=3.5, ls='--', color=cp.darkred )

    ax.set_xlabel('Delay Time', fontsize=20 )
    ax.set_ylabel('Number of SN', fontsize=20 )


    ax.set_ylim( 0, 2 )
    ax.set_xlim( 0, 5 )
    # ax.set_xticklabels( [] )
    ax.set_yticklabels( [] )
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 2.0 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 1.0 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1.0 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    pl.draw()

def mk_snr_fig():
    prepfig()

    from pytools import cosmo
    t = np.arange( 0.04, 10, 0.01 )

    z = cosmo.zfromt( t )
    sfrlinear = z + 0.1

    dtdgauss = 0.7 * np.exp( -( (t-2.5)**2 / 1.5**2 )) # (1./(2*np.pi*0.5) ) *
    dtdt1 = 1. / (25.*t)

    snrgauss = np.convolve( sfrlinear, dtdgauss, mode='full' )[:len(t)]
    snrt1 = np.convolve( sfrlinear, dtdt1, mode='full' )[:len(t)]

    ax = pl.gca()
    ax.plot( z, 1.8*snrgauss/snrgauss.max(), lw=3.5, ls='-', color=cp.darkgreen )
    ax.plot( z, 2.4 * snrt1/snrt1.max(), lw=3.5, ls='--', color=cp.darkred )

    ax.set_xlabel('Redshift', fontsize=20 )
    ax.set_ylabel('SN Rate', fontsize=20 )


    ax.set_ylim( 0, 2 )
    ax.set_xlim( 0, 4 )
    # ax.set_xticklabels( [] )
    ax.set_yticklabels( [] )
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 1.0 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1.0 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    pl.draw()



