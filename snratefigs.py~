from matplotlib import ticker
from matplotlib import pyplot as pl
from hstsnpipe.tools.figs import colors, plotsetup
from cosmo import agez, zfromt
import candelsRates, candelsDTD
import numpy as np

def mkCSFRfig():

    import sys
    import os

    thisfile = sys.argv[0]
    if 'ipython' in thisfile : thisfile = __file__
    _thisdir = os.path.dirname( os.path.abspath( thisfile ) )
    _csfrdatfile = os.path.join( _thisdir, 'csfrB13.dat' )
    
    # read in the data (already includes systematic error estimates)
    z,log10csfr, errplus, errminus =  np.loadtxt(_csfrdatfile,unpack=True)

    # the Cole et al 2001 double power-law model
    sfrModel = lambda z, z0,a,b,c : c / ( 10**(a*(z-z0)) + 10**(b*(z-z0)))

    # Best fit model values
    (z0,a,b,c) = (1.243,-0.997,0.241,0.180)

    sfrobs = 10**(log10csfr)
    dsfrobs = errplus * sfrobs / (2.5 * np.log10(np.e) )
    zmod = np.arange(0.0001, 11, 0.01 )
    sfrmod = sfrModel( zmod, z0, a, b, c )
    # logsfrmod = np.log10( sfrmod )

    # plot the points and overlay the best-fit model and uncertainties
    plotsetup.usetex()
    plotsetup.presfig()
    pl.clf()
    ax = pl.gca()

    ax.errorbar( z, sfrobs/0.11, dsfrobs/0.11, ls=' ', marker='o', color=colors.lightblue, mew=0.8,mec='k',mfc=colors.lightblue,capsize=0)
    # pl.plot( z, sfrobs, ls=' ', marker='o', color='0.4', mew=0.8,mec='k',mfc='0.5')

    ax.plot( zmod, sfrmod/0.11, color=colors.blue, ls='-', lw=2.6)

    ax = pl.gca()
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"Cosmic Star Formation Rate (Normalized)")
    ax.set_xlim( 0, 11)
    ax.set_ylim( 0, 1.98)

    fig = pl.gcf()
    fig.subplots_adjust( left=0.08, bottom=0.12, top=0.88, right=0.95 )
    axtop = ax.twiny()
    ageticks = np.array( [13,1.0,3.0,0.5] )
    zageticks = zfromt( ageticks )
    axtop.set_xticks( zageticks )
    axtop.set_xticklabels( ageticks )
    axtop.set_xlabel('Age of Universe [Gyr]')
    axtop.set_xlim( ax.get_xlim() )


    ax.xaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )

    pl.draw()



def mkRatesFigGND( zrange=[0.01,3.5], grey=False,
                 ageaxis=True, showCSFR=False, lines=False, **kwargs ):
    """ 
    Making a figure for a public talk showing the observed rates vs z 
    before candels + clash
    """
    from ratetable import ALLHST, LOWZ, RATELISTGROUND, RATELISTHST
    import dtdconvolve

    plotsetup.presfig(figsize=[8,5])
    pl.clf()

    ax = pl.gca()

    if grey : c = '0.7'
    else : c = colors.maroon
    for R in RATELISTGROUND : 
        R.mfc=c
        R.mec=c
        R.color='0.5'
        R.ms=8
        R.marker='o'
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100)
   
    if showCSFR: 
        # the Cole et al 2001 double power-law model
        sfrModel = lambda z, z0,a,b,c : c / ( 10**(a*(z-z0)) + 10**(b*(z-z0)))

        # Best fit model values
        (z0,a,b,c) = (1.243,-0.997,0.241,0.180)

        zmod = np.arange(0.0001, 11, 0.01 )
        sfrmod = sfrModel( zmod, z0, a, b, c ) / 0.07
        pl.plot( zmod, sfrmod, color=colors.blue, ls='-', marker=' ')
    elif not grey : 
        ax.text(0.95,0.9,'Ground-based SN Surveys',fontsize='x-large',ha='right',va='top', transform=ax.transAxes)

    if lines : 
        z = np.arange( zrange[0], zrange[1], 0.1 )
        for data,c,ls in zip(['GROUND','HST'],[colors.blue,colors.green],['--','-.']) : 
            if data == 'GROUND' : 
                eta, fp = 1e-4, 0.75
            elif data == 'HST' : 
                eta, fp = 3e-4, 0.15

            snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=False, 
                                      t0=0.5, tmin=0.04, f0=fp, eta=eta, **kwargs )
            line = ax.plot( z, snr, color=c, ls=ls, zorder=2000 )

    fig = pl.gcf()
    fig.subplots_adjust( left=0.08, bottom=0.12, top=0.88, right=0.95 )

    ax.set_xlabel('Redshift')
    ax.set_ylabel('Cosmic Thermonuclear SN Rate')
    ax.set_xlim([-0.05, zrange[-1]+0.05])
    ax.set_ylim([0, 1.98])
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_ticklabels([])

    if ageaxis : 
        axtop = ax.twiny()
        axtop.set_xlim( ax.get_xlim() )
        ageticks = np.array( [13,8,5,3] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')

    pl.draw()
    return( ax )

def mkRatesFigHST( zrange=[0.01,3.5], 
                   redpoints=False, meanpoints=False, lines=False, 
                   ffsn=False,
                 **kwargs ):
    """ 
    Making a figure for public talks showing
    observed rates and t-1 power law models with different values for the 
    fraction of SNIa that are prompt (t<500 Myr).
    """
    from cosmo import agez, zfromt
    import dtdconvolve
    import candelsRates, candelsDTD
    from ratetable import ALLHST, HIZHST, HIZFF, LOWZ, RATELISTGROUND, RATELISTHST

    ax = mkRatesFigGND( grey=True, showCSFR=lines )

    # ALLHST.plot(zorder=100)
    if redpoints : c=colors.darkred
    else : c='0.7'
    for R in RATELISTHST : 
        R.mfc=c
        R.mec=c
        R.color=c
        R.ms=8
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100)

    if meanpoints and not ffsn: 
        HIZHST.plot(zorder=100)
    elif ffsn : 
        HIZFF.plot(zorder=100)


    if lines : 
        z = np.arange( zrange[0], zrange[1], 0.1 )
        for data,c,ls in zip(['GROUND','HST'],[colors.maroon,colors.green],['--','-.']) : 
            if data == 'GROUND' : 
                eta, fp = 1.84e-4, 0.498
            elif data == 'HST' : 
                eta, fp = 3.47e-4, 0.046

            snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=False, 
                                      t0=0.5, tmin=0.04, f0=fp, eta=eta, **kwargs )
            line = ax.plot( z, snr, color=c, ls=ls, zorder=2000 )


    fig = pl.gcf()
    fig.subplots_adjust( left=0.08, bottom=0.12, top=0.88, right=0.95 )

    ax.set_xlim([-0.05, zrange[-1]+0.05])
    ax.set_ylim([0, 1.98])

    ax.yaxis.set_ticklabels([])
    pl.draw()
    return( ax )


def mkSNRpdfs():

    mkRatesFigGND(grey=False, showCSFR=False)
    pl.savefig('snr1.pdf')

    mkRatesFigGND(grey=False, showCSFR=True)
    pl.savefig('snr2.pdf')

    mkRatesFigHST(redpoints=True,  meanpoints=False, lines=False, ffsn=False)
    pl.savefig('snr3.pdf')

    mkRatesFigHST(redpoints=False, meanpoints=False, lines=False, ffsn=False)
    pl.savefig('snr4.pdf')

    mkRatesFigHST(redpoints=False, meanpoints=True, lines=False, ffsn=False)
    pl.savefig('snr5.pdf')

    mkRatesFigHST(redpoints=False, meanpoints=True, lines=True, ffsn=False)
    pl.savefig('snr6.pdf')

    mkRatesFigHST(redpoints=False, meanpoints=False, lines=True, ffsn=True)
    pl.savefig('snr7.pdf')
