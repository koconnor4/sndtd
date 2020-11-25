# 2013.10.25 S.Rodney
# Fitting the Behroozi:2013 CSFR(z) compilation 
#  with a double-power law model 
import numpy as np
import sys
import os
from matplotlib import pyplot as pl
import matplotlib.ticker


from astropy.io import ascii
from astropy import table

import astropy.units as u
#from astropy.cosmology import FlatLambdaCDM
#from scipy.interpolate import interp1d

#LCDM = FlatLambdaCDM( 70, 0.3 )
#zsamples = np.logspace( -3, 0.7, num=50, base=10 )
#tsamples = LCDM.age( zsamples )
#getzfromt = interp1d( tsamples, zsamples, bounds_error=False, fill_value=tsamples[-1], assume_sorted=True )

#def zfromt( t ):
#    return( getzfromt( t ) )

from pytools import cosmo
zfromt = lambda t : cosmo.zfromt( t, 70, 0.3, 'Gyr' )


thisfile = sys.argv[0]
if 'ipython' in thisfile : thisfile = __file__
_thisdir = os.path.dirname( os.path.abspath( thisfile ) )
_B13datfile = os.path.join( _thisdir, 'sfrB13_dietSalpeter.fw.dat' )
_HB06datfile = os.path.join( _thisdir, 'sfrHB06_dietSalpeter.fw.dat' )
_HB10datfile = os.path.join( _thisdir, 'sfrHB10_dietSalpeter.fw.dat' )
_SFRcurvesfile = os.path.join( _thisdir, 'sfrCurves_dietSalpeter.dat' )


HB06 = ascii.read( _HB06datfile , format='fixed_width')  # z, dz, sfr, err
HB10 = ascii.read( _HB10datfile , format='fixed_width')  # z, dz, sfr, err
B13 = ascii.read( _B13datfile   , format='fixed_width')  # z, sfr, err
sfrcurves = ascii.read( _SFRcurvesfile , format='commented_header', header_start=-1, data_start=0)


"""
NOTE on IMF conversions :
to convert SFR from a Salpeter (1955) IMF to a Chabrier (2003) IMF :

   Mcha = Msal / 1.8 = 0.56 * Msal

  ==> log10( SFRcha ) = log10( SFRsal ) - 0.25

to convert from Salpeter to diet-Salpeter (Bell & de Jong 2001,2003) :

   Mdie = Msal / 1.43 = 0.7 * Msal

 ==> log10( SFRdie ) = log10( SFRsal ) - 0.15


to convert from Salpeter to Kroupa (following Reddy & Steidel 2009):

   Mkro = Msal / 1.7 = 0.59 * Msal

 ==> log10( SFRkro ) = log10( SFRsal ) - 0.16


"""

from pytools import plotsetup


def sfrplot_for_chicago( slide=1 ):
    plotsetup.presfig(figsize=[12,6])
    pl.clf()

    """
Best fit model parameters for fitting the B13 compilation
with a double-power-law SFR(z) model :
       z0    a     b     c
high  1.293 -0.980  0.219  0.269
mid   1.234 -0.979  0.236  0.176
low   1.176 -0.980  0.254  0.116
sharp 1.121 -1.054  0.287  0.157
flat  1.464 -0.867  0.203  0.218
"""

    if slide==1 : alpha=0.5
    else : alpha=0.15
    pl.errorbar( HB06['z']+np.random.uniform(-0.05,0.05,len(HB06['z'])), HB06['sfr'], HB06['err'], marker='o', color='k', alpha=alpha, ls=' ', label='HB06' )
    pl.errorbar( HB10['z']+np.random.uniform(-0.05,0.05,len(HB10['z'])), HB10['sfr'], [HB10['errsfrneg'],HB10['errsfrpos']] , marker='o', color='r', alpha=alpha, ls=' ', label='HB10' )

    if slide == 2 :
        pl.errorbar( B13['z']+np.random.uniform(-0.05,0.05,len(B13['z'])),  B13['sfr'], [B13['errneg'],B13['errpos']], marker='o', color='b', alpha=0.5, ls=' ', label='B13' )

    for zrange,ls, alpha in zip( [[0,8.2], [8.2,10.2]], ['-','--'], [0.3,0.15] ):
        z = np.arange( zrange[0], zrange[1], 0.01 )
        HB06fit = double_powerlaw( z, z0=0.840, A=-1.311, B=0.085, C=0.143 ) *1.82/1.43
        B13fit = double_powerlaw( z, z0=1.243, A=-0.979, B=0.236, C=0.176 ) *1.82/1.43
        B13high = double_powerlaw( z, z0=1.293, A=-0.980, B=0.219, C=0.269 ) *1.82/1.43
        B13low = double_powerlaw( z, z0=1.176, A=-0.980, B=0.254, C=0.116 ) *1.82/1.43

        if slide==1 :
            pl.plot( z, HB06fit, ls=ls, marker=' ', color='k', label='Hopkins+Beacom 2006' )
        else :
            pl.plot( z, HB06fit, ls=ls, marker=' ', color='0.5', label='Hopkins+Beacom 2006' )
            pl.plot( z, B13fit, ls=ls, marker=' ', color='b', label='Behroozi+ 2013' )
            pl.fill_between( z, B13low, B13high, color='c', alpha=alpha )

    ax = pl.gca()
    ax.set_yscale('log', basey=10 )
    ax.set_xlim( -0.3, 9.95 )
    ax.set_ylim( 0.001, 0.5 )
    ax.set_xlabel('Redshift')
    ax.set_ylabel(r'$\dot{\rho}_{*}$ \ [\  M$_{\odot}$ yr$^{-1}$  Mpc$^{-3}$\  ]')
    ax.set_yticks( [0.001, 0.01, 0.1 ] )
    ax.set_yticklabels( ['0.001', '0.01', '0.1' ] )

    # ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))

    axtop = ax.twiny()
    axtop.set_xlim( ax.get_xlim() )
    ageticks = np.array( [13,6,3, 1,0.5] )

    zageticks = zfromt( ageticks )
    axtop.set_xticks( zageticks )
    axtop.set_xticklabels( ageticks )
    axtop.set_xlabel('Age of Universe (Gyr)')

    fig = pl.gcf()
    fig.subplots_adjust( left=0.12, bottom=0.12, right=0.95, top=0.88 )
    pl.draw()

def plot_more():

    if logplot :
        ax = pl.gca()
        ax.set_yscale('log', base=10)

    pl.draw()


def plot_csfr_models( ):
    """ plot the best-fit double-power law models from Behroozi 2013,
    rescaled to a diet-Salpeter IMF
    :return:
    """
    z = np.arange( 0.0, 15, 0.01 )
    L08 = cole_model( z, a=0.0157, b=0.118, c=3.23, d=4.66, )
    H06 = double_powerlaw( z, z0=0.840, A=-1.311, B=0.085, C=0.143 ) *1.65/1.43
    B13 = double_powerlaw( z, z0=1.243, A=-0.997, B=0.241, C=0.180 ) *1.82/1.43

    pl.plot( z, H06, ls='-', marker=' ', color='k', label='Hopkins+Beacom 2006' )
    # pl.semilogy( z, L08, ls='-', marker=' ', color='b', label='Li+ 2008' )
    pl.plot( z, B13, ls='-', marker=' ', color='b', label='Behroozi+ 2013' )


def cole_model(  z, a=0.0157, b=0.118, c=3.23, d=4.66, ):
    """  CSFR model using the Cole et al 2001 parameterization"""
    return( (a + b* z) / ( 1 + (z/c)**d) )


def double_powerlaw( z, z0=1.243, A=-0.997, B=0.241, C=0.180 ):
    """  Double power-law fit a la Behroozi et al 2013"""
    return( C / ( 10**(A*(z-z0)) + 10**(B*(z-z0)) )   )


def sfrB13( z, dust='mid' ):
    """  The Behroozi+ 2013 SFR(z) curve, for a diet-Salpeter IMF
    :param z:
    :return:
    """
    if dust=='mid' :
        # Peter Behroozi's fit to the B13 data, before adjustment
        # i.e. assuming a Chabrier IMF
        # return( double_powerlaw(z, A= -0.997, B = 0.241, C = 0.180, z0 = 1.243 ) )

        # Steve's fit to the B13 data, including adjustment to a diet-Salpeter IMF
        return( double_powerlaw( z, z0=1.243, A=-0.979, B=0.236, C=0.176 ) *1.82/1.43 )
    elif dust=='high':
        return( double_powerlaw( z, z0=1.293, A=-0.980, B=0.219, C=0.269 ) *1.82/1.43 )
    elif dust=='low' :
        return( double_powerlaw( z, z0=1.176, A=-0.980, B=0.254, C=0.116 ) *1.82/1.43 )



def sfrHB06( z ) :
    """  The Behroozi+ 2013 fit to the HB06 SFR(z) data, converted to a diet-Salpeter IMF """
    return( double_powerlaw( z, z0=0.840, A=-1.311, B=0.085, C=0.143 ) *1.82/1.43 )


def mkfigCSFR( ):
    """ 
    make a paper-ready figure showing the CSFR data, best-fit model 
    and uncertainties
    """
    from matplotlib import pyplot as pl
    fitB13csfr( showplot=True, ageaxis=True )

def fitB13csfr( showplot=False, ageaxis=True ):
    """ 
    fit curves to the Behroozi et al 2013b Cosmic Star Formation Rate data
    showplot : make a paper-ready figure showing the best-fit model with uncertainties
    """
    from matplotlib import pyplot as pl
    from scipy import optimize as scopt
    import plotsetup
    from cosmo import zfromt

    # read in the data (already includes systematic error estimates)
    z,log10csfr, errplus, errminus =  np.loadtxt(_csfrdatfile,unpack=True)

    # in practice, errplus==errminus, so this is irrelevant
    err = (errplus + errminus)/2.

    # the Cole et al 2001 double power-law model
    sfrModel = lambda z, z0,a,b,c : c / ( 10**(a*(z-z0)) + 10**(b*(z-z0)))

    # B13 best-fit parameters for the double power law model 
    (z00,a0,b0,c0) = (1.243,-0.997,0.241,0.180)

    def chi2func( param , systerr=0 ):
        """ 
        systerr = -1 : reduce all obs to lower systematic limit
        systerr =  0 : no change in rates observations
        systerr = +1 : reduce all obs to upper systematic limit
        """
        z0,a,b,c = param
        return( ( ( log10csfr+(systerr*err) - np.log10(sfrModel(z,z0,a,b,c)) )**2 / err**2 ).sum() )

    if False : 
        # showing the chi2 variance vs each of the 4 parameters
        z0range = np.arange(1.0,2.5,0.01)
        chi2array = np.array( [ chi2func( [z0new,a0,b0,c0] ) for z0new in z0range ] )
        pl.subplot( 221 )
        pl.plot(z0range, chi2array, 'b-' )

        a0range = np.arange( 0.2*a0,2*a0,0.05*a0 )
        chi2array = np.array( [ chi2func( [z00,a0new,b0,c0] ) for a0new in a0range ] )
        pl.subplot( 222 )
        pl.plot(a0range, chi2array, 'g-' )

        b0range = np.arange( 0.2*b0,2*b0,0.05*b0 )
        chi2array = np.array( [ chi2func( [z00,a0,b0new,c0] ) for b0new in b0range ] )
        pl.subplot( 223 )
        pl.plot(b0range, chi2array, 'r-' )

        c0range = np.arange( 0.2*c0,2*c0,0.05*c0 )
        chi2array = np.array( [ chi2func( [z00,a0,b0,c0new] ) for c0new in c0range ] )
        pl.subplot( 224 )
        pl.plot(c0range, chi2array, 'm-' ) 

        return(0) 

    parambounds = [(0.5,2.5),(-2.0,0.5),(0.0,0.5),(0.0,0.5)]
    parambounds = [(0.2,3.5),(-2.5,1.5),(-0.5,1.0),(-0.5,1.0)]
    output = scopt.fmin_l_bfgs_b( chi2func, [z00,a0,b0,c0], args=(0,),
                                  fprime=None, approx_grad=1, 
                                  bounds=parambounds )
    fitparBest, chi2min, outdict = output
    z0best, abest, bbest, cbest = fitparBest

    # Fit the data points after shifting to the systematic extrema :

    # HIGH DUST : 
    # Passing args=(1,) sets systerr=+1 in the chi2func() calls, meaning 
    # that all the observed points are shifted up to the positive tip of their 
    # systematic error bars (i.e. a stronger dust correction)
    outputPlus = scopt.fmin_l_bfgs_b( chi2func, [z00+0.1,a0+0.1,b0+0.1,c0+0.1], args=(1,),
                                      fprime=None, approx_grad=1, 
                                      bounds=parambounds )
    fitparBestP, chi2minP, outdictP = outputPlus
    z0bestP, abestP, bbestP, cbestP = fitparBestP


    # LOW DUST : 
    # Passing args=(-1,) sets systerr=-1 in the chi2func() calls
    outputMinus = scopt.fmin_l_bfgs_b( chi2func, [z00-0.1,a0-0.1,b0-0.1,c0-0.1], args=(-1,),
                                      fprime=None, approx_grad=1, 
                                      bounds=parambounds )
    fitparBestM, chi2min, outdict = outputMinus
    z0bestM, abestM, bbestM, cbestM = fitparBestM

    # SHARP : 
    # Passing an array to set systerr in the chi2func() calls such that 
    # all points with z above the turnover point are shifted down and those 
    # below are shifted up
    systerrSharp = np.where( z>2, -1, 0 )
    outputSharp = scopt.fmin_l_bfgs_b( chi2func, [z00-0.1,a0-0.1,b0-0.1,c0-0.1], args=(systerrSharp,),
                                      fprime=None, approx_grad=1, 
                                      bounds=parambounds )
    fitparBestS, chi2min, outdict = outputSharp
    z0bestS, abestS, bbestS, cbestS = fitparBestS

    # FLAT : 
    # Passing an array to set systerr in the chi2func() calls such that 
    # all points with z above the turnover point are shifted up and those 
    # below are shifted down
    systerrFlat = np.where( z>2, +1, 0 )
    outputFlat = scopt.fmin_l_bfgs_b( chi2func, [z00-0.1,a0-0.1,b0-0.1,c0-0.1], args=(systerrFlat,),
                                      fprime=None, approx_grad=1, 
                                      bounds=parambounds )
    fitparBestF, chi2min, outdict = outputFlat
    z0bestF, abestF, bbestF, cbestF = fitparBestF

    print("Best fit model parameters:")
    print( "       z0    a     b     c " )
    print( "high  %5.3f %6.3f %6.3f %6.3f"%(tuple(fitparBestP) ))
    print( "mid   %5.3f %6.3f %6.3f %6.3f"%(tuple(fitparBest) ))
    print( "low   %5.3f %6.3f %6.3f %6.3f"%(tuple(fitparBestM) ))
    print( "sharp %5.3f %6.3f %6.3f %6.3f"%(tuple(fitparBestS) ))
    print( "flat  %5.3f %6.3f %6.3f %6.3f"%(tuple(fitparBestF) ))

    if showplot : 
        # plot the points and overlay the best-fit model and uncertainties
        plotsetup.usetex()
        plotsetup.halfpaperfig(1)
        pl.errorbar( z, log10csfr, [errplus,errminus], ls=' ', marker='o', color='0.4', mew=0.8,mec='k',mfc='0.5',capsize=0)

        zmod = np.arange(0.0001, 11, 0.01 )
        sfrPlus = np.log10( sfrModel( zmod, z0bestP, abestP, bbestP, cbestP ))
        sfrMinus = np.log10( sfrModel( zmod, z0bestM, abestM, bbestM, cbestM ))
        pl.fill_between( zmod, sfrMinus, sfrPlus, color='c', alpha=0.4)

        sfrFlat = np.log10( sfrModel( zmod, z0bestF, abestF, bbestF, cbestF ))
        sfrSharp = np.log10( sfrModel( zmod, z0bestS, abestS, bbestS, cbestS ))

        pl.plot( zmod, np.log10( sfrModel( zmod, z0best, abest, bbest, cbest )), color='k', ls='-', lw=2.6)
        pl.plot( zmod, sfrSharp, color='darkorange', ls='-.', lw=2.2)
        pl.plot( zmod, sfrFlat, color='forestgreen', ls='--', lw=2.2)

        ax = pl.gca()
        ax.set_xlabel("Redshift")
        ax.set_ylabel(r"log$_{10}$( $\dot{\rho}_{*}$ ) [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]")
        ax.set_xlim( 0, 11)

        if ageaxis :              
            fig = pl.gcf()
            fig.subplots_adjust( top=0.88 )
            axtop = ax.twiny()
            ageticks = np.array( [13,1.0,3.0,0.5] )
            zageticks = zfromt( ageticks )
            axtop.set_xticks( zageticks )
            axtop.set_xticklabels( ageticks )
            axtop.set_xlabel('Age of Universe (Gyr)')
            axtop.set_xlim( ax.get_xlim() )
        ax.text(0.95,0.95,'Data compilation from \n Behroozi et al. 2013',
                ha='right',va='top',fontsize=10, transform=ax.transAxes )
        pl.draw()

    if False : 
        # unbounded minimization fails to converge
        output = scopt.fmin_bfgs( chi2func, [z00,a0,b0,c0], fprime=None, args=(0,), 
                                  full_output=1, disp=1 )
        fitparBest, chi2min, gopt, invhess, func_calls, grad_calls, warnflag = output
        # take the square roots of the diagonals of the inverse hessian
        # matrix to get the standard err of the best fit parameter values
        z0err,aerr,berr,cerr = np.sqrt( np.diagonal( invhess ) )

        #pl.plot( z, np.log10(sfrMod(z,z00+z0err,a0+aerr,b0+berr,c0+cerr)), 'r:')
        #pl.plot( z, np.log10(sfrMod(z,z00-z0err,a0-aerr,b0-berr,c0-cerr)), 'r:')


