# 2013.10.25 S.Rodney
# Fitting the CANDELS rate measurements with DTD models
import numpy as np
from scipy.stats.distributions import poisson
from hstsnpipe.tools import snana
from hstsnpipe.tools.figs.colors import *
color1='m'
color2=darkgreen

import cosmo
_zFirst = 20
_H0=70
_Om=0.3
_Ode=0.7

#  * age of univ. at z of First Stars  [Gyr]
_tFS = cosmo.agez(_zFirst, H0=_H0, Om=_Om, Ode=_Ode, unit='Gyr') 
#  * age of universe today [Gyr]
_tNow = cosmo.agez( 0, H0=_H0, Om=_Om, Ode=_Ode, unit='Gyr') 


def fitSNRsystematicsSampling():
    """ explore systematic uncertainties by sampling across 
    the space of all possible combinations of data, dust and csfr
    (takes about 5 minutes on lucan)
    ( see bottom of this module for a shortcut )
    ( End result, as calculated on 2013.10.25 : 

    ==> CANDELS : fPrompt = 0.50   +0.07 -0.09 (syst) 
    ==> GROUND  : fPrompt = 0.50   +0.07 -0.09 (syst) 
    ==> HST :    fPrompt = 0.04   +0.13 -0.04 (syst) 
    ==> ALL :    fPrompt = 0.44   +0.08 -0.13 (syst)

    """
    dtdParam = {}
    for data in ['c+c','ground','all'] : 
        etalist, detalist = [], []
        fplist, dfplist = [], []
        for dust in ['mid','high','low'] : 
            for csfr in ['mid','flat','sharp']:
                eta, deta, fp, dfp = fitSNR( data, dust, csfr, Ngrid=0 )
                etalist.append(eta)
                detalist.append(deta)
                fplist.append(fp)
                dfplist.append(dfp)
        dtdParam[data] = {'eta':np.array(etalist), 
                          'deta':np.array(detalist),
                          'fP':np.array(fplist),
                          'dfP':np.array(dfplist) }
    return( dtdParam )


def fitSNR( data='ALL', dust='mid', csfr='mid', t0=0.5, tmin=0.04, 
            Ngrid=30, minchi2=False, skipN2M=True, showplots=False, 
            sneakpeek=False, Ngauss=100000):
    """ fit the set of all non-redundant rates estimates 
    with a DTD * SFR model that uses a t-1 power law after 
    500 Myr, and a constant factor from tmin=40 Myr to t0=500 Myr.
    The model efficiency eta is a free parameter
    and the fraction of prompt SNe is another free parameter.
    (eta is largely constrained by the z<1 data and fPrompt by the 
    z>1 data)

    INPUTS : 
     data : choose the data set to fit 
        'R13' : CANDELS data only
        'ALL' : use all available SNIa rates measurements
        'HST' : use only the HST surveys (GOODS,CSS,CLASH,CANDELS)
        'GROUND' : use only the ground-based surveys
        'LOWZ' : use only the ground-based surveys (excluding G11)
    
     dust : set the dust assumption, and use it to make a 
        systematic error correction for both the SNIa rates
        and the CSFR(z) curve

     t0,tmin : set the prompt-component parameters of the DTD model 
         tmin is the time in Gyr to produce the first SNIa [40 Myr]
         t0 is the time in Gyr at which the t-1 DTD shape begins [500 Myr]

      Ngrid : if Ngrid > 0 then we also construct a grid of likelihood
    values with Ngrid steps along each of the eta and fP parameter
    axes
    """
    # Get all the rate points and redshift points from all the rate
    # observations and make a set of giant arrays with all the
    # individual (single-redshift) rate measures spread out
    import time
    from math import factorial
    from numdifftools import Hessian
    from ratetable import RATELISTGROUND,RATELISTALL,RATELISTLOWZ, \
        GROUND,ALL,LOWZ,ALLHST, D08,B12,G13,R13
    RATELISTHST = [D08,B12,G13,R13]
    start = time.time()

    if Ngrid : 
        etagrid = np.linspace( 0.9e-4, 2.0e-4, Ngrid )
        fPgrid = np.linspace( 0.35, 0.8, Ngrid )

    if data.lower() in ['r13','r14'] : 
        ratelist = [ R13 ]
        composite = R13
        if Ngrid : 
            etagrid = np.linspace( 0.51e-4, 4.7e-4, Ngrid )
            fPgrid = np.linspace( 0.01, 0.8, Ngrid )
    elif data.lower() in ['c&c','c+c'] : 
        ratelist = [ G13, R13 ]
        composite = ALLHST
        if Ngrid : 
            etagrid = np.linspace( 0.51e-4, 4.7e-4, Ngrid )
            fPgrid = np.linspace( 0.01, 0.8, Ngrid )
    elif data.lower() == 'all' : 
        ratelist = RATELISTALL
        composite = ALL
    elif data.lower() == 'hst' : 
        ratelist = RATELISTHST
        composite = ALLHST
        if Ngrid : 
            etagrid = np.linspace( 0.51e-4, 4.7e-4, Ngrid )
            fPgrid = np.linspace( 0.01, 0.8, Ngrid )
    elif data.lower() == 'ground' : 
        ratelist = RATELISTGROUND
        composite = GROUND
    elif data.lower() == 'lowz' : 
        ratelist = RATELISTLOWZ
        composite = LOWZ
    robs,robsP,robsM = [] ,[],[]
    robsSysP,robsSysM = [],[]
    robsStatP,robsStatM = [],[]
    zobs,zobsP,zobsM = [] ,[],[]
    nobs = []
    for R in ratelist : 
        robs += R.rate.tolist() 
        nobs += R.ndet.tolist()
        robsP += R.errplus.tolist() 
        robsM += R.errminus.tolist() 
        robsSysP += R.errsysplus.tolist() 
        robsSysM += R.errsysminus.tolist() 
        robsStatP += R.errstatplus.tolist() 
        robsStatM += R.errstatminus.tolist() 
        zobs += R.z.tolist()
        zobsP += R.zerrplus.tolist()
        zobsM += R.zerrminus.tolist()
    isort = np.argsort( zobs ) 
    robs = np.array( robs )[isort] 
    robsP = np.array( robsP )[isort] 
    robsM = np.array( robsM )[isort] 
    robsSysP = np.array( robsSysP )[isort]
    robsSysM = np.array( robsSysM )[isort]
    robsStatP = np.array( robsStatP )[isort]
    robsStatM = np.array( robsStatM )[isort]
    zobs = np.array( zobs )[isort]
    zobsP = np.array( zobsP )[isort]
    zobsM = np.array( zobsM )[isort]
    
    # do imports and define values that don't change
    # outside of the minimizing function for efficiency
    from dtdconvolve import snrate
    from matplotlib import pyplot as pl
    from scipy import optimize as scopt
    sig2 = np.max([robsStatP,robsStatM],axis=0)**2

    # initial guesses for the DTD model parameters
    # Note that we use scaling factors so that the values
    # handled by the scipy minimizer are close to unity.
    if minchi2:
        eta0, etascale = 1., 2e-4
        fprompt0, fpscale = 1., 0.5
    else : 
        eta0, etascale = 1.,  1e-4, 
        fprompt0, fpscale = 0.5, 1.

    Ndof = len(robs) - 3

    # Set the CSFR(z) curve and the observed points based on the 
    # systematic error test case defined by the dust and csfr parameters
    SFRmodel='B13'+csfr  # B13mid, high, low, flat, sharp
    if dust=='high' : robsfit = robs+robsSysP
    elif dust=='mid' : robsfit = robs
    elif dust=='low' : robsfit = robs-robsSysM
    label='%s data, %s dust, %s CSFR(z)'%(data, dust, csfr)

    def neglikefunc( fitpar ) : 
        """ negative likelihood function to minimize. 
        fitpar[0] = eta : the efficiency factor for Ia production
        fitpar[1] = fprompt : the prompt Ia fraction
        """
        eta = fitpar[0] * etascale
        fprompt = fitpar[1] * fpscale

        # define the SNR model using a given prompt fraction
        # evaluated at each z bin center
        rmod = snrate( zobs, DTDmodel='t-1+p', SFRmodel=SFRmodel, eta=eta, f0=fprompt, t0=t0, tmin=tmin, tmax=13.5)
        # TODO : 
        # integrate the SNR model over the redshift bin to match 
        # the observed bin
        # Technically, the integration should be in time, not z

        # TODO : test how much difference we get when you 
        # just evaluate the SNR model at the bin center instead
        # of integrating

        # compute the total likelihood value for this parameter set
        negloglike = 0.
        for i in range(len(rmod)) :
            # Defining the asymmetric errors using Poisson distributions

            # The model defines lambda, using the inferred rate-to-Ndet scaling 
            # factor from the rate observations
            lam = rmod[i] * ( nobs[i] / robsfit[i] )

            # The observed number of SN detections is k : 
            k = int( round( nobs[i] ) )

            # so here is the poisson probability for this k and lambda :
            # like = (lam**k) * np.exp( -lam )  / factorial( k ) 
            ## and here is the log of that likelihood for this k and lambda : 
            # loglike = np.log( like )
            ## or, recasting for computational efficiency:  
            loglike = k * np.log( lam ) - lam - np.sum( [ np.log(ki) for ki in range(1,k+1) ] )

            negloglike -= loglike
        return( negloglike )

    if sneakpeek:
        from matplotlib import pyplot as pl
        pl.figure(1)
        pl.clf()
        eta = 2.0
        neglikefp= []
        fplist = np.arange( 0.01,0.99,0.05)
        for fp in fplist :
            neglikefp.append( neglikefunc( [eta,fp] ) )
        pl.plot( fplist, neglikefp, 'r-', marker=' ')
        pl.xlabel('fprompt')
        pl.ylabel('neg log like')
        ax1 = pl.gca()
        ax1.set_ylim(0,300)

        neglikeEta= []
        etalist = np.arange( 0.01,5.00,0.1)
        fp = 0.5
        for eta in etalist :
            neglikeEta.append( neglikefunc( [eta,fp] ) )
        pl.figure(2)
        pl.clf()
        pl.plot( etalist, neglikeEta, 'g-', marker=' ')
        pl.xlabel('eta')
        pl.ylabel('neg log like')
        ax2 = pl.gca()
        ax2.set_ylim(0,300)
        pl.draw()

        return( fplist, etalist, neglikefp, neglikeEta )

    # maximize the likelihood (i.e. minimize the negative log likelihood)
    output = scopt.fmin_l_bfgs_b( neglikefunc, [eta0,fprompt0], 
                                  fprime=None, args=(), approx_grad=True,
                                  bounds=[ (0,None), (0,0.99) ] )
    fitparBest, likemax, infodict = output
    hessian = Hessian( fun=neglikefunc )
    invhess = np.linalg.inv( hessian( fitparBest ) )

    # estimate the uncertainties from the sqrt of the diagonal elements 
    # of the inverse hessian matrix :
    etaErr, fpromptErr = np.sqrt( np.diagonal( invhess ) )
    etaBest = fitparBest[0] 
    fpromptBest = fitparBest[1]

    # rescale eta and fprompt to their natural units
    etaBest *= etascale
    etaErr *= etascale
    fpromptBest *= fpscale
    fpromptErr *= fpscale

    if Ngrid : 
        # compute the likelihood over a 2-d grid of parameter values
        # and define 68% and 95% confidence contours

        likegrid = []
        N2Mgrid = []
        for e in etagrid : 
            for fP in fPgrid : 
                likegrid.append( np.exp(-neglikefunc( [e/etascale,fP/fpscale] ) ) )
                if skipN2M : N2M = 0
                else : N2M, dN2M = getN2M( e, fP, fittoSFRobs=SFRmodel )
                N2Mgrid.append( N2M )

        likegrid = np.reshape( likegrid, (Ngrid,Ngrid) )
        N2Mgrid = np.reshape( N2Mgrid, (Ngrid,Ngrid) )

    # report   eta, chi2min, stderr
    print( label )
    print( "  likemax = %.3e "%(likemax ) )
    print( "  eta = %.5e  +-%.5e"%(etaBest, etaErr ) )
    print( "  fprompt = %.3f  +-%.3f\n"%(fpromptBest, fpromptErr ) )

    if showplots : 
        from matplotlib import pyplot as pl
        rmodA = snrate( zobs, DTDmodel='t-1+p', eta=etaBest, f0=fpromptBest, t0=0.5, tmin=0.01, tmax=13.5)
        #rmodB = snrate( zobs, DTDmodel='t-1+p', eta=etaBest+etaErr, f0=fpromptBest, t0=0.5, tmin=0.01, tmax=13.5)
        #rmodC = snrate( zobs, DTDmodel='t-1+p', eta=etaBest-etaErr, f0=fpromptBest, t0=0.5, tmin=0.01, tmax=13.5)
        #rmod1 = snrate( zobs, DTDmodel='t-1+p', eta=etaBest, f0=fpromptBest+fpromptErr, t0=0.5, tmin=0.01, tmax=13.5)
        #rmod2 = snrate( zobs, DTDmodel='t-1+p', eta=etaBest, f0=fpromptBest-fpromptErr, t0=0.5, tmin=0.01, tmax=13.5)

        pl.clf()
        ax = pl.gca()
        #ax.errorbar( zobs, robsfit, [robsStatM,robsStatP],  marker='o', capsize=0, ls=' ', elinewidth=0.5, color='0.4' )
        #ax.plot( zobs, rmod2, 'm:' )
        #ax.plot( zobs, rmodB, 'g--' )
        ax.plot( zobs, rmodA, 'r-' )
        #ax.plot( zobs, rmodC, 'g--' )
        #ax.plot( zobs, rmod1, 'm:' )
        for ipt in range(len(robsfit)):
            if nobs[ipt]>Ngauss : 
                color='r'
                marker='s'
            else : 
                color = '0.3'
                marker= 'o'
            
            ax.errorbar( [zobs[ipt]], [robsfit[ipt]], 
                         yerr=[[robsStatM[ipt]],[robsStatP[ipt]]],
                         ls=' ', color=color, marker=marker )
                             
        # composite.plot()
        pl.draw()

    end = time.time()
    print( "total execution time = %.1f seconds"%(end-start) )
    if Ngrid : return( likegrid, N2Mgrid, etagrid, fPgrid )
    else : return( etaBest, etaErr, fpromptBest, fpromptErr )

def printN2Mvalues( ):
    """ 
    1. Evaluate N2M for a set of hardcoded eta and fP values, to
    define the best-fit NIa/M* values from each of the three data sets
    of interest : Ground, C+C, All.

    2. Sum up the CANDELS+CLASH, and ALL data sets to get an estimate of 
    the total number of SNIa explosions over a hubble time, and divide
    by the total sum of all stars formed to get a DTD-independent 
    measure of NIa / M*
    """

    print( "computing the NIa/M* values with the assumed t-1 DTD form:")


    # CLASH+CANDELS :
    n2mH      = 1e3* getN2M( (2.25+0.00)*1e-4, 0.21+0.00, fittoSFRobs='B13')[0]
    n2mHpStat = 1e3* getN2M( (2.25+1.36)*1e-4, 0.21+0.34, fittoSFRobs='B13')[0] - n2mH
    n2mHmStat = 1e3* getN2M( (2.25-1.18)*1e-4, 0.21-0.21, fittoSFRobs='B13')[0] - n2mH
    n2mHpSyst = 1e3* getN2M( (2.25+0.72)*1e-4, 0.21+0.49, fittoSFRobs='B13')[0] - n2mH
    n2mHmSyst = 1e3* getN2M( (2.25-0.15)*1e-4, 0.21-0.12, fittoSFRobs='B13')[0] - n2mH
    print("CAN+CL:   NIa/M* = %.2f %+.2f %+.2f (stat) %+.2f %+.2f (syst)"%(n2mH, n2mHpStat, n2mHmStat, n2mHpSyst, n2mHmSyst ) )

    # GROUND :
    n2mG      = 1e3* getN2M( (1.38+0.00)*1e-4, 0.53+0.00, fittoSFRobs='B13')[0]
    n2mGpStat = 1e3* getN2M( (1.38+0.24)*1e-4, 0.53+0.09, fittoSFRobs='B13')[0] - n2mG
    n2mGmStat = 1e3* getN2M( (1.38-0.23)*1e-4, 0.53-0.10, fittoSFRobs='B13')[0] - n2mG
    n2mGpSyst = 1e3* getN2M( (1.38+0.25)*1e-4, 0.53+0.10, fittoSFRobs='B13')[0] - n2mG
    n2mGmSyst = 1e3* getN2M( (1.38-0.59)*1e-4, 0.53-0.26, fittoSFRobs='B13')[0] - n2mG
    print("GROUND:   NIa/M* = %.2f %+.2f %+.2f (stat) %+.2f %+.2f (syst)"%(n2mG, n2mGpStat, n2mGmStat, n2mGpSyst, n2mGmSyst ) )

    # ALL :
    n2mA      = 1e3* getN2M( (1.60+0.00)*1e-4, 0.53+0.00, fittoSFRobs='B13')[0]
    n2mApStat = 1e3* getN2M( (1.60+0.24)*1e-4, 0.53+0.09, fittoSFRobs='B13')[0] - n2mA
    n2mAmStat = 1e3* getN2M( (1.60-0.23)*1e-4, 0.53-0.10, fittoSFRobs='B13')[0] - n2mA
    n2mApSyst = 1e3* getN2M( (1.60+0.25)*1e-4, 0.53+0.10, fittoSFRobs='B13')[0] - n2mA
    n2mAmSyst = 1e3* getN2M( (1.60-0.59)*1e-4, 0.53-0.26, fittoSFRobs='B13')[0] - n2mA
    print("ALL:      NIa/M* = %.2f %+.2f %+.2f (stat) %+.2f %+.2f (syst)"%(n2mA, n2mApStat, n2mAmStat, n2mApSyst, n2mAmSyst ) )


    from hstsnpipe.tools.rates.ratetable import ALLCC, ALLHST, ALL
    from cosmo import agez
    from scipy.integrate import quad
    import sfrate

    # Integrate the CSFR to get Mtot in unites of Msun yr-1
    def SFR( t ) : return( sfrate.behroozi2013( t, abscissa='age', fitto='B13') )
    def SFRhigh( t ) : return( sfrate.behroozi2013( t, abscissa='age', fitto='B13high') )
    def SFRlow( t ) : return( sfrate.behroozi2013( t, abscissa='age', fitto='B13low') )
    Mtot = quad( SFR, 0, _tNow )[0]
    MtotHigh = quad( SFRhigh, 0, _tNow )[0]
    MtotLow = quad( SFRlow, 0, _tNow )[0]
    dMpsys = MtotHigh - Mtot
    dMmsys = Mtot - MtotLow

    print( "\n\ncomputing the NIa/M* values without assuming a DTD form:")
    for R, label in zip([ALLCC,ALLHST,ALL],['CAN+CL: ','ALL HST: ','ALL:     ']) :
        z = R.z
        dzp = R.zerr[0]
        dzm = R.zerr[1]
        dt = agez( z - dzm ) - agez( z + dzp )  # Gyr
        NIa = (R.rate * dt ).sum()
        dNIapstat = ((R.rate + R.errstatplus ) * dt ).sum()
        dNIamstat = ((R.rate - R.errstatminus) * dt ).sum()
        dNIapsys =  ((R.rate + R.errsysplus  ) * dt ).sum()
        dNIamsys =  ((R.rate - R.errsysminus ) * dt ).sum()

        # The SNIa rate carries a facotor 1e-4, and we want to report 
        # the NIa/M* values in units of 1e-3 Msun^-1, so divide by 10
        n2m = NIa / Mtot * 0.1
        n2mpStat = dNIapstat / Mtot  * 0.1
        n2mmStat = -dNIamstat / Mtot  * 0.1

        n2mpSys = n2m * np.sqrt(  (dNIapsys/NIa)**2  + (dMmsys/Mtot)**2 )
        n2mmSys = -n2m * np.sqrt(  (dNIamsys/NIa)**2  + (dMpsys/Mtot)**2 )

        print("%s  NIa/M* = %.2f %+.2f %+.2f (stat) %+.2f %+.2f (syst)"%(label, n2m, n2mpStat, n2mmStat, n2mpSys, n2mmSys ) )
        


def getN2M( eta=2e-4, fP=0.444, dt=0.05, t0=0.5, tmin=0.04, fittoSFRobs='B13', 
            showplots=True, verbose=False, debug=False) :
    """ 
    Integrate SNR(t) dt to get the total number of SNe formed (Ntot) over 
    all of cosmic history, and integrate SFR(t) dt to get the total mass 
    of all stars formed. Then divide to get the global SNIa production 
    efficiency :   NIa / M*

    INPUTS:
    dt : time step for seeding the interpolator [Gyr]
    t0 : age of universe in Gyr where the t-1 power law SNR start [0.5 Gyr]
    tmin : delay time needed for the first SNIa explosion [0.04 Gyr = 40 Myr]
    fittoSFRobs : version of the B13 compilation that sets the parameters for
                  the SFR(t) curve  [B13, B13high, B13low, B13flat, B13peaky]

    RETURNS: 
    N2M = NIa / Mtot : the number of type Ia SNe per stellar mass 
    dN2M : the error due to numerical integration

    NOTE :  This is MUCH faster but slightly less precise than getN2Mslow, 
      b/c here we use a numpy discrete convolution to sample the SNR(t) 
      function at a finite set of points, then define a linear 
      interpolator that can evaulate SNR(t) at any time.
    """
    import exceptions
    import sfrate, dtdconvolve
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    if debug: import pdb; pdb.set_trace()
    import time
    start = time.time()

    def SFR( t ) : 
        return( sfrate.behroozi2013( t, abscissa='age', fitto=fittoSFRobs) )

    # Create a time array, in units of Gyr, sampling the age of the 
    # universe from the epoch of first star formation to the present
    tcosmic = np.arange( _tFS+dt/2., _tNow, dt )

    # Fill arrays with discrete values for the Star formation history 
    # and the DTD at each time in tcosmic.
    # Note that we evaluate the SFR at times in the tcosmic array, 
    # relative to t=0 at the big bang. We evaluate the DTD at 
    # times tcosmic-tFS, relative to t=0 at the epoch of first stars.
    SFRdisc = sfrate.behroozi2013( tcosmic, abscissa='age', fitto=fittoSFRobs) 
    DTDdisc = dtdconvolve.dtd_t1_fprompt( tcosmic-_tFS, eta=eta, f0=fP, t0=t0, tmin=_tFS, tmax=_tNow )

    # Do a discrete convolution of the SFR and DTD
    SNRdisc = dt * np.convolve( SFRdisc, DTDdisc, mode='full')[:len(tcosmic)]

    # Set up an interpolator function to evaluate the SNR at any 
    # time t, relative to the big bang
    SNR = interp1d( tcosmic, SNRdisc, kind='linear', bounds_error=False, fill_value=0 )

    # use the QUADPACK routine for numerical integration
    # of the SFR(t) and SNR(t) functions to get the N/M ratio
    Mtot,dMtot = quad( SFR, 0, _tNow )[0:2]
    NIatot,dNIatot = quad( SNR, 0, _tNow )[0:2]
    N2M = NIatot/Mtot
    dN2M = abs( N2M * np.sqrt( (dNIatot/NIatot)**2 + (dMtot/Mtot)**2 ) )
    end = time.time()

    if verbose : 
        print(" USING A DISCRETE CONVOLUTION AND INTERPOLATOR ")
        print(" NIatot : %.4e +- %.4e"%(NIatot,dNIatot) )
        print(" Mtot : %.4e +- %.4e"%(Mtot,dMtot) )
        print(" N/M : %.4e +- %.4e"%(N2M, dN2M) )
        print(" %.3e sec \n\n "%(end - start ) )
    return( N2M, dN2M )

def getN2Mslow( eta=2e-6, fP=0.444, dt=0.05, t0=0.5, tmin=0.04, fittoSFRobs='B13', 
                showplots=True, verbose=False, debug=False) :
    """
    There are two ways to evaluate the SNR(t) convolution integral 
    This is the slower, but more precise method, in which we 
    numerically integrate the convolution integral using QUADPACK
    """
    import sfrate, dtdconvolve
    from scipy.integrate import quad
    if debug: import pdb; pdb.set_trace()
    import time
    start = time.time()

    # Define the integrands for the cosmic star formation rate, 
    # the DTD convolution integral, and the SNR at time t. 
    def SFR( t ) : 
        return( sfrate.behroozi2013( t, abscissa='age', fitto=fittoSFRobs) )

    def convInt( tprime, t=13.6 ) : 
        return( SFR(t-tprime) * dtdconvolve.dtd_t1_fprompt( tprime, eta=eta, f0=fP, t0=t0, tmin=_tFS, tmax=_tNow ) )

    def SNR( t ) : 
        return( quad( convInt, 0, t, args=(t,) )[0] )

    # use the QUADPACK routine for numerical integration
    Mtot,dMtot = quad( SFR, 0, _tNow )[0:2]
    NIatot,dNIatot = quad( SNR, 0, _tNow )[0:2]
    endQuad = time.time()

    # use the QUADPACK routine for numerical integration
    # of the SFR(t) and SNR(t) functions to get the N/M ratio
    Mtot,dMtot = quad( SFR, 0, _tNow )[0:2]
    NIatot,dNIatot = quad( SNR, 0, _tNow )[0:2]
    N2M = NIatot/Mtot
    dN2M = abs( N2M * np.sqrt( (dNIatot/NIatot)**2 + (dMtot/Mtot)**2 ) )
    end = time.time()

    if verbose : 
        print(" USING QUAD NUMERICAL INTEGRATION FOR CONVOLUTION INTEGRAL")
        print(" NIatot : %.4e +- %.4e"%(NIatot,dNIatot) )
        print(" Mtot : %.4e +- %.4e"%(Mtot,dMtot) )
        print(" N/M : %.4e +- %.4e"%(N2M, dN2M) )
        print(" %.3e sec \n\n "%(end - start ) )
    return( N2M, dN2M )

def plotSNRparamContoursMultiDust( likegridDict, N2MgridDict, etagridDict, fPgridDict ): 
    for label,ls,lw,c in zip( ['mid','high','low'],
                           ['solid','dashed','dashdot'],
                           [1.5, 0.8, 0.8 ],
                           ['k','darkorange','darkcyan'] ):
        plotSNRparamContours( likegridDict[label], N2MgridDict[label], etagridDict[label], fPgridDict[label], 
                              linestyles=ls, lineweights=lw, colors=c )

def plotSNRparamContours( likegrid, etagrid, fPgrid, N2Mgrid=None, newfig=False, 
                          showlabels=False, **kwarg ): 
    """ 
    plot likelihood contours in fP vs eta space

    INPUT :
    likegrid : an array of likelihood values 
    etagrid : an array of eta values (y axis)
    fPgrid  : an array of fP values (x axis)

    OPTIONAL :
    N2Mgrid : an array of N2M values (for bg coloring)
    kwarg : keyword arguments passed to contour()
    newfig : start with a new Figure instance
    showlabels : put labels on the contour lines
    """
    from matplotlib import pyplot as pl
    from matplotlib import cm
    if newfig: 
        import plotsetup
        fig = plotsetup.halfpaperfig(100)
        fig.subplots_adjust( left=0.14, right=0.91, bottom=0.12, top=0.96 )

    plotdefaults = {'colors':'k', 'linestyles':'-', 'linewidths':1.5, 'zorder':500 }
    plotargs = dict( plotdefaults.items() + kwarg.items() )

    # normalize the likelihood array so that it sums to 1 
    pgrid = likegrid/likegrid.sum()
    # convert the prob array into a "cumulative sum above" array for probability contours
    pcumgrid = scumsum( pgrid ) 
    # plot it 
    ax = pl.gca()
    ax.set_xlabel(r'f$_{P}$ : prompt SN Ia fraction')
    ax.set_ylabel(r'$\eta$ : DTD normalization factor [10$^{-4}$ yr$^{-1}$ M$_{\odot}^{-1}$]')
    CS = ax.contour( fPgrid, etagrid*1e4, pcumgrid, levels=[0.95,0.68,0], 
                     **plotargs ) #colors=cm.gist_rainbow_r )

    if showlabels :
        # Define a class that forces representation of float to look a certain way
        # This remove trailing zero so '1.0' becomes '1'
        class nf(float):
             def __repr__(self):
                 return '%.0f' % self.__float__()

        # Recast levels to new class
        CS.levels = [ nf(val*100) for val in CS.levels ]

        # Label levels with specially formatted floats
        if pl.rcParams["text.usetex"]:
            fmt = r'%r \%%'
        else:
            fmt = '%r %%'
        pl.clabel(CS, CS.levels, inline=True, fontsize=10, fmt=fmt)

    # pl.clabel(CS, inline=1, fontsize=10)
    if N2Mgrid != None : 
        imax = ax.imshow( N2Mgrid*1e3, interpolation='nearest', # aspect='equal', 
                          extent=[fPgrid[0],fPgrid[-1],etagrid[0]*1e4,etagrid[-1]*1e4],
                          vmin=0.1, vmax=2.1, zorder=-500, cmap= cm.BrBG ) #cm.Greys_r) #cm.gist_rainbow_r) 
        cbar = pl.colorbar( imax, ticks=np.arange(0.25,2.1,0.25) )
        cbar.set_label('N$_{Ia}$/M$_*$ [10$^{-3}$ M$_{\odot}^{-1}$]', rotation=-90, labelpad=18)
        cbar.ax.set_yticklabels( ['0.25','0.50','0.75','1.00','1.25','1.50','1.75','2.0+'] )


def mkDTDcontoursAll( middata=None, Ngrid=30):
    """ make a figure showing the confidence regions in eta vs fP space, 
    using the all-data constraints
    """
    if middata==None : middata = fitSNR( data='ALL', dust='mid', Ngrid=Ngrid, skipN2M=False)
    chi2mid, likemid, N2Mmid, etamid, fPmid = middata
    plotSNRparamContours( likemid, etamid, fPmid, N2Mgrid=N2Mmid, newfig=True, showlabels=True )
    return( middata )

def mkDTDcontoursGround( middata=None, Ngrid=30):
    """ make a figure showing the confidence regions in eta vs fP space, 
    using the all-data constraints
    """
    if middata==None : middata = fitSNR( data='ground', dust='mid', Ngrid=Ngrid, skipN2M=False)
    chi2mid, likemid, N2Mmid, etamid, fPmid = middata
    plotSNRparamContours( likemid, etamid, fPmid, N2Mgrid=N2Mmid, newfig=True, showlabels=True )
    return( middata )


def mkDTDcontoursSystSample( middata=None, highdata=None, lowdata=None, Ngrid=30):
    """ make a figure showing the confidence regions in eta vs fP space, 
    comparing the all-data constraints for three different dust assumptions
    to probe systematic uncertainties
    """

    if middata==None : middata = fitSNR( data='ALL', dust='mid', Ngrid=Ngrid, skipN2M=False)
    if highdata==None: highdata = fitSNR( data='ALL', dust='high', Ngrid=Ngrid, skipN2M=True)
    if lowdata ==None : lowdata = fitSNR( data='ALL', dust='low', Ngrid=Ngrid, skipN2M=True) 

    chi2mid, likemid, N2Mmid, etamid, fPmid = middata
    chi2high, likehigh, N2Mhigh, etahigh, fPhigh = highdata
    chi2low, likelow, N2Mlow, etalow, fPlow = lowdata

    plotSNRparamContours( likemid, etamid, fPmid, N2Mgrid=N2Mmid, newfig=True, showlabels=True )
    plotSNRparamContours( likehigh, etahigh, fPhigh, colors='0.3', linestyles='--' )
    plotSNRparamContours( likelow, etalow, fPlow, colors='0.3', linestyles='--' )

    return( middata, highdata, lowdata )


def mkDTDcontoursHSTvsGround( hstdata=None, grounddata=None,  Ngrid=30):
    """ make a figure showing the confidence regions in eta vs fP space, 
    comparing the ground-only constraints and the HST-only constraints
    """
    if hstdata==None :   hstdata = fitSNR( data='C+C', dust='mid', Ngrid=Ngrid, skipN2M=False )
    if grounddata==None: grounddata = fitSNR( data='GROUND', dust='mid', Ngrid=Ngrid, skipN2M=True )

    likeH, N2MH, etaH, fPH = hstdata
    likeG, N2MG, etaG, fPG = grounddata

    plotSNRparamContours( likeH, etaH, fPH, N2Mgrid=N2MH, newfig=True, showlabels=False, linestyles='--', color='k' )
    plotSNRparamContours( likeG, etaG, fPG, N2Mgrid=None, newfig=False, showlabels=False, linestyles='-', colors='m')

    from matplotlib import pyplot as pl
    ax = pl.gca()
    ax.text( 0.60, 1.6, 'Ground-based', ha='center', va='center', color='m', rotation=-33)
    ax.text( 0.25, 3.25, 'CANDELS+CLASH', ha='center', va='center', color='k', rotation=-55)
    return( hstdata, grounddata )

def scumsum( a ): 
    """ 
    Sorted Cumulative Sum function : 
    Construct an array "sumabove" such that the cell at index i in sumabove 
    is equal to the sum of all cells from the input array "a" that have a
    cell value higher than a[i]  
    (useful for plotting confidence region contours)
    """    
    # Collapse the array into 1 dimension
    sumabove = a.ravel()
    
    # Sort the raveled array by descending cell value
    iravelsorted = sumabove.argsort( axis=0 )[::-1]

    # Reassign each cell to be the cumulative sum of all 
    # input array cells with a higher value :  
    sumabove[iravelsorted] = sumabove[iravelsorted].cumsum() 

    # Now unravel back into shape of original array and return
    return( sumabove.reshape( a.shape ) )



from matplotlib import pyplot as pl
import numpy as np


def mkfigToyModelDTD( presfig=True ):
    """ 
    make a figure showing the DTD toy model with a variable 
    prompt fraction.
    """
    from hstsnpipe.tools.figs import plotsetup, colors
    if presfig : plotsetup.presfig()
    else : plotsetup.halfpaperfig()
    from dtdconvolve import dtd_t1_fprompt as dtd

    t = np.arange( 0, 10, 0.001 )
    dtd1 = dtd( t, eta = 1, f0=0.1 )
    dtd5 = dtd( t, eta = 1, f0=0.5 )
    dtd8 = dtd( t, eta = 1, f0=0.8 )
    pl.clf()
    pl.loglog( t, dtd1, color=colors.red, ls='-.')
    pl.loglog( t, dtd5, color=colors.green, ls='--')
    pl.loglog( t, dtd8, color=colors.blue, ls='-')
    ax = pl.gca()
    ax.set_xlim( 0.01, 10.01 )
    ax.set_ylim( 0.1, 100 )
    ax.set_xticklabels( ['','100 Myr', '1 Gyr', '10  Gyr'] )
    ax.set_yticklabels( ['','1', '10', '100'] )

    #ax.axvline( 0.04, ls=':', color='0.3' )
    #ax.text( 0.038, 0.2, '40 Myr', ha='right', va='bottom', rotation=90, fontsize='small', color='0.3' )
    ax.text( 0.04, 0.15, '40 Myr', ha='center', va='bottom', rotation=90, fontsize='small', color='0.3' )
    ax.text( 0.5, 0.15, '500 Myr', ha='center', va='bottom', rotation=90, fontsize='small', color='0.3' )

    ax.text( 0.055, 40, 'Prompt Component', ha='left', va='bottom', fontsize='small', color='0.3' )
    ax.text( 0.1, 24, 'f$_{P}$=0.8', ha='center', va='top', fontsize='small', color=colors.blue)
    ax.text( 0.1, 6, 'f$_{P}$=0.5', ha='center', va='top', fontsize='small', color=colors.green)
    ax.text( 0.1, 0.6, 'f$_{P}$=0.1', ha='center', va='top', fontsize='small', color=colors.red)
    ax.text( 2.,  1.8, 'Delayed Component', ha='center', va='center', fontsize='small', color='0.3', rotation=-30 )
    ax.text( 2.,  1., '$\propto \eta t^{-1}$', ha='center', va='center', color='0.3', rotation=-30 )

    ax.set_xlabel('Time since formation')
    ax.set_ylabel(r'SN rate~~~[ yr$^{-1}$ M$_{\bigodot}^{-1}$]')
    fig = pl.gcf()
    fig.subplots_adjust( left=0.1, bottom=0.13, right=0.95, top = 0.95 )


def mkfigRatesDTD( zFirst=20, zrange=[0.01,2.8], twopanels=False, 
                   logscale=False, ageaxis=True, label=True, presfig=False, **kwargs ):
    """ 
    Making a figure for the CANDELS rates paper that shows 
    observed rates and t-1 power law models with different values for the 
    fraction of SNIa that are prompt (t<500 Myr).
    """
    from cosmo import agez, zfromt
    import dtdconvolve
    from ratetable import ALLHST, LOWZ, RATELISTGROUND, RATELISTHST
    from hstsnpipe.tools.figs import plotsetup

    if twopanels==1 : 
        mkfigRatesDTD( zFirst=zFirst, zrange=zrange, twopanels=2, 
                       logscale=logscale, ageaxis=ageaxis, label=label, 
                       presfig=1, **kwargs )
        mkfigRatesDTD( zFirst=zFirst, zrange=zrange, twopanels=3, 
                       logscale=logscale, ageaxis=ageaxis, label=label, 
                       presfig=1, **kwargs )

    if logscale : scale = 1e-4
    else : scale = 1

    if presfig :
        if twopanels : 
            plotsetup.presfig( int(twopanels) )
        else : 
            plotsetup.presfig()
    else :
        plotsetup.fullpaperfig(1, figsize=[8,4.5])

    top=0.88
    ax1 = pl.axes( [0.12,0.12, 0.86, 0.76] )
    for R in RATELISTGROUND : 
        R.mfc='0.3'
        R.mec='0.3'
        R.color='0.3'
        R.ms=8
        R.marker='o'
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.errsysplus=0*R.errsysplus
        R.errsysminus=0*R.errsysminus
        R.plot( thicksys=False, zorder=-100, scalerates=scale)

    if twopanels==3 or not twopanels : 
        #ALLHST.z += 0.07
        #ALLHST.zerrminus += 0.02
        #ALLHST.zerrplus -= 0.02
        #ALLHST.plot(zorder=100,scaleunits=scale)
        for R in RATELISTHST : 
            if R.reference in [ 'CANDELS','CLASH'] : 
                    alpha=0.8
            else : 
                    R.mfc='w'
                    R.mec='0.3'
                    R.color='0.3'
                    R.ms=8
                    alpha=1
            R.zerrplus=0*R.zerrplus
            R.zerrminus=0*R.zerrminus
            R.errsysplus=0*R.errsysplus
            R.errsysminus=0*R.errsysminus
            R.plot( thicksys=False, zorder=10000, scalerates=scale, alpha=alpha)

    z = np.arange( zrange[0], zrange[1], 0.1 )
    snrlist = []
    for data,ax,c,ls in zip(['GROUND','HST'],[ax1,ax1],[color1,color2],['-','--']) : 
        if data == 'GROUND' : 
            eta, fp = 1.38e-4, 0.59
        elif data == 'HST' : 
            # eta, fp = 3.66e-4, 0.05
            eta, fp = 2.25e-4, 0.21
            if twopanels==2: continue
        snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=False, 
                                  t0=0.5, tmin=0.04, f0=fp, eta=eta, **kwargs )
        line = ax1.plot( z, snr*scale, color=c, ls=ls, zorder=2000, lw=2 )

    ax1.set_xlim([-0.05, zrange[-1]+0.05])
    ax1.set_xlabel('Redshift')

    if logscale : 
        ax1.set_ylabel('SNIa Rate [ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]', ha='center', va='center' )
        ax1.set_yscale('log')
        ax1.set_ylim([0.101*scale, 2.55*scale])
    else : 
        ax1.set_ylabel('SNIa Rate [ 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]', ha='center', va='center' )
        ax1.set_ylim( 0, 1.49 )

    if logscale  : 
        ax1.text( 2.78,scale*0.85,r'\noindent Fraction of SN Ia\\that are "prompt"\\ ($\tau<$500 Myr) :',ha='right',va='bottom', backgroundcolor='w')
        ax1.text( 2.32,scale*0.57,r'$f_{P}=0.6$',ha='left',va='top',color=darkgreen,backgroundcolor='w',fontsize=14,zorder=10000)
        ax1.text( 2.25,scale*0.27,r'$f_{P}=0.2$',ha='left',va='top',color='teal',backgroundcolor='w',fontsize=14,zorder=10000)
        #ax1.text(0.0,scale*1.8,'Weighted Ave.,\nGround-based',ha='left',va='top',color='saddlebrown')
        #ax1.text(1.4,scale*0.17,'Weighted Ave., HST',ha='right',va='top',color='blue')
        #ax1.plot([0.20,0.34],[scale*1.20,scale*0.47],ls='-',marker=' ',color='saddlebrown',lw=0.5)
        #ax1.plot([1.21,1.58],[scale*0.18,scale*0.31],ls='-',marker=' ',color='blue',lw=0.5)
        ax1.yaxis.set_label_coords( -0.08, 0.5 )

    elif label:                    
        if presfig: 
            ax1.text( 2.29,0.88,r'Fraction of SN Ia',ha='left',va='bottom', fontsize='small')
            ax1.text( 2.29,0.79,r'that are "prompt"',ha='left',va='bottom', fontsize='small')
            ax1.text( 2.29,0.70,r'($\tau<$500 Myr) :',ha='left',va='bottom', fontsize='small')
            ax1.text( 2.32,0.63,r'$f_{P}=0.6$',ha='left',va='top',color=color1,backgroundcolor='w',fontsize=14,zorder=10000, rotation=-19)
            ax1.text( 2.32,0.32,r'$f_{P}=0.2$',ha='left',va='top',color=color2,backgroundcolor='w',fontsize=14,zorder=10000, rotation=-16)
            ax1.set_xlim(-0.05, 2.95)
            #ax1.text( 2.82,0.85,r'\noindent Fraction of SN Ia\\that are "prompt"\\ ($\tau<$500 Myr) :',ha='right',va='bottom', backgroundcolor='w', fontsize='large')
            #if twopanels==2 or not twopanels : 
            #    ax1.text( 2.32,0.63,r'$f_{P}=0.6$',ha='left',va='top',color='m',backgroundcolor='w',zorder=10000, rotation=-19)
            #if twopanels==3 or not twopanels : 
            #    ax1.text( 2.32,0.32,r'$f_{P}=0.2$',ha='left',va='top',color='teal',backgroundcolor='w',zorder=10000, rotation=-16)
        else : 
            ax1.text( 2.3,0.88,r'Fraction of SN Ia',ha='left',va='bottom', fontsize='large')
            ax1.text( 2.3,0.79,r'that are "prompt"',ha='left',va='bottom', fontsize='large')
            ax1.text( 2.3,0.70,r'($\tau<$500 Myr) :',ha='left',va='bottom', fontsize='large')
            ax1.text( 2.32,0.63,r'$f_{P}=0.6$',ha='left',va='top',color=color1,backgroundcolor='w',fontsize=14,zorder=10000, rotation=-19)
            ax1.text( 2.32,0.32,r'$f_{P}=0.2$',ha='left',va='top',color=color2,backgroundcolor='w',fontsize=14,zorder=10000, rotation=-16)

        # if twopanels==2 or not twopanels : 
        #     ax1.text(0.02,1.4,'Weighted Average,\nGround-based',ha='left',va='top',color='saddlebrown')
        #     ax1.plot([0.32,0.55],[1.23,0.65],ls='-',marker=' ',color='saddlebrown',lw=0.5)
        # if twopanels==3 or not twopanels : 
        #     ax1.text(2.5,1.35,'Weighted Average, HST',ha='right',va='top',color='blue')
        #     ax1.plot([1.7,2.],[0.65,1.27],ls='-',marker=' ',color='blue',lw=0.5)
    ax1.yaxis.set_label_coords( -0.08, 0.5 )

    if ageaxis : 
        axtop = ax1.twiny()
        axtop.set_xlim( ax1.get_xlim() )
        ageticks = np.array( [13,8,5,3] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')
    
    fig = pl.gcf()
    fig.subplots_adjust(0.08, 0.15, 0.98, 0.85)

    pl.draw()
    return( ax1 )


def mkfigFFSN( zFirst=20, zrange=[0.01,3.5], twopanels=False,
               logscale=False, ageaxis=True, label=True, paperfig=True, **kwargs ):
    """ 
    Making a figure for Fall 2013 job applications that shows 
    observations and DTD-based curves, including projections for the FFSN sample.
    """
    from cosmo import agez, zfromt
    import dtdconvolve
    from ratetable import GROUND, LOWZ, HIZHST, HIZFF, RATELISTALL, ALL 

    if logscale : scale = 1e-4
    else : scale = 1

    if paperfig :
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.fullpaperfig(figsize=[4.5,4.5])
    else :
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.presfig(figsize=[8,5])

    ax1 = pl.axes( [0.16,0.12, 0.80, 0.76] )
    for R in RATELISTALL : 
        R.mfc='0.7'
        R.mec='0.7'
        R.color='0.7'
        R.ms=8
        R.marker='o'
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100, scalerates=scale)

    LOWZ.plot(zorder=100, scalerates=scale)
    HIZHST.plot(zorder=200, scalerates=scale)
    HIZFF.plot(zorder=300, scalerates=scale)

    z = np.arange( zrange[0], zrange[1], 0.1 )
    snrlist = []
    for data,ax,c,ls in zip(['GROUND','HST'],[ax1,ax1],['m','teal'],['-','--']) : 
        #eta, etaErr, fprompt, fpromptErr = candelsDTD.fitSNR( \
        #    data=data, dust='mid', t0=0.5, tmin=0.04, 
        #    Ngrid=0, skipN2M=True, showplots=False )
        
        if data == 'GROUND' : 
            eta, fp = 1.7e-4, 0.55
            ax = ax1

        elif data == 'HST' : 
            eta, fp = 3.47e-4, 0.046
            eta, fp = 2.8e-4, 0.1
            ax = ax1

        snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=False, 
                                  t0=0.5, tmin=0.04, f0=fp, eta=eta, **kwargs )

        line1 = ax1.plot( z, snr*scale, color=c, ls=ls, zorder=2000 )
        line2 = ax1.plot( z, snr*scale, color=c, ls=ls, zorder=2000 )

    # ax.text( 2.7,1.7,r'DTD$\propto t^{-1}$\\ for $t>500 Myr$',ha='right',va='top')
    ax1.set_xlim([-0.05, 3.75])
    ax1.set_xlabel('Redshift')

    if logscale : 
        ax1.set_ylabel('SNIa Rate [ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]', ha='center', va='center' )
        ax1.set_yscale('log')
        ax1.set_ylim([0.101*scale, 2.55*scale])
    else : 
        ax1.set_ylabel('SNIa Rate [ 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]', ha='center', va='center' )
        ax1.set_ylim( 0, 1.49 )

    if logscale  : 
        # ax1.text( 2.78,scale*0.85,r'\noindent Fraction of SN Ia\\that are "prompt"\\ ($\tau<$500 Myr) :',ha='right',va='bottom', backgroundcolor='w')
        # ax1.text( 2.32,scale*0.57,r'$f_{P}=0.5$',ha='left',va='top',color='m',backgroundcolor='w',fontsize=14,zorder=10000)
        # ax1.text( 2.25,scale*0.27,r'$f_{P}=0.05$',ha='left',va='top',color='teal',backgroundcolor='w',fontsize=14,zorder=10000)
        # ax1.text(0.0,scale*1.8,'Weighted Ave.,\nGround-based',ha='left',va='top',color='saddlebrown')
        # ax1.text(1.4,scale*0.17,'Weighted Ave., HST',ha='right',va='top',color='blue')
        # ax1.plot([0.20,0.34],[scale*1.20,scale*0.47],ls='-',marker=' ',color='saddlebrown',lw=0.5)
        # ax1.plot([1.21,1.58],[scale*0.18,scale*0.31],ls='-',marker=' ',color='blue',lw=0.5)
        ax1.yaxis.set_label_coords( -0.12, 0.5 )

    else  : 
        #ax1.text( 2.72,0.85,r'\noindent Fraction of SN Ia\\that are "prompt"\\ ($\tau<$500 Myr) :',ha='right',va='bottom', backgroundcolor='w')
        #ax1.text( 2.92,0.52, r'$f_{P}=0.5$',ha='left',va='top',color='m',fontsize=14,zorder=10000, rotation=-35)
        #ax1.text( 3.2,0.18,r'0.05',ha='left',va='top',color='teal',fontsize=14,zorder=10000, rotation=-20)

        ax1.text( 3.7,0.42,r'\noindent Prompt\\ SNIa\\ Fraction',ha='right',va='bottom', backgroundcolor='w')
        ax1.text( 3.3,0.37, '0.5',ha='left',va='top',color='m',fontsize=14,zorder=10000,backgroundcolor='w')
        ax1.text( 3.2,0.12,'0.05',ha='left',va='top',color='teal',fontsize=14,zorder=10000,backgroundcolor='w')

        ax1.text(0.05,0.1,'Ground-based',ha='left',va='top',color='saddlebrown',backgroundcolor='w')
        ax1.text(3.5,1.35,'HST (to date)',ha='right',va='top',color='darkgreen',backgroundcolor='w')
        ax1.text(3.65,1.18,'$+$FrontierSN\n(projected)',ha='right',va='top',color='darkblue',backgroundcolor='w')
        #ax1.plot([0.32,0.55],[1.23,0.65],ls='-',marker=' ',color='saddlebrown',lw=0.5)
        ax1.plot([2.3,2.68],[0.9,1.27],ls='-',marker=' ',color='darkgreen',lw=0.5)
        ax1.plot([2.6,3.00],[0.39,1.12],ls='-',marker=' ',color='darkblue',lw=0.5)
        ax1.yaxis.set_label_coords( -0.12, 0.5 )

    if ageaxis : 
        axtop = ax1.twiny()
        axtop.set_xlim( ax1.get_xlim() )
        ageticks = np.array( [13,8,5,3,2] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')
    
    fig = pl.gcf()
    fig.subplots_adjust(0.08, 0.15, 0.98, 0.85)

    pl.draw()
    return( ax1, ax1 )



def testplots():
    from ratetable import RATELISTGROUND, RATELISTHST, R13, G13
    from dtdconvolve import snrate
    from matplotlib import pyplot as pl
    import numpy as np

    pl.clf()
    for R in RATELISTGROUND : 
        #R.mfc='g'
        #R.mec='g'
        #R.marker='o'
        R.errsysplus=np.zeros( len(R.z) )
        R.errsysminus=np.zeros( len(R.z) )
        R.plot(thicksys=False)

    for R in RATELISTHST : 
        R.mfc='g'
        R.mec='g'
        R.marker='^'    
        R.errsysplus=np.zeros( len(R.z) )
        R.errsysminus=np.zeros( len(R.z) )
        R.plot(thicksys=False)

    R13.mfc='r'
    R13.mec='r'
    R13.marker='D'    
    R13.errsysplus=np.zeros( len(R13.z) )
    R13.errsysminus=np.zeros( len(R13.z) )
    R13.plot(thicksys=False)

    G13.mfc='darkorange'
    G13.mec='darkorange'
    G13.marker='s'    
    G13.errsysplus=np.zeros( len(G13.z) )
    G13.errsysminus=np.zeros( len(G13.z) )
    G13.plot(thicksys=False)


    etahst = 3.93609e-04
    fphst = 0.000

    etagnd = 1.37778e-04
    fpgnd = 0.590

    z = np.arange( 0.01, 2.5, 0.05 )
    gndmod = snrate( z, DTDmodel='t-1+p', SFRmodel='B13mid', 
                     eta=etagnd, f0=fpgnd, t0=0.5, tmin=0.04, tmax=13.5)    
    hstmod = snrate( z, DTDmodel='t-1+p', SFRmodel='B13mid', 
                     eta=etahst, f0=fphst, t0=0.5, tmin=0.04, tmax=13.5)    

    ccmod = snrate( z, DTDmodel='t-1+p', SFRmodel='B13mid', 
                    eta=2.25e-4, f0=0.21, t0=0.5, tmin=0.04, tmax=13.5)    

    r13mod = snrate( z, DTDmodel='t-1+p', SFRmodel='B13mid', 
                     eta=1.26e-4, f0=0.53, t0=0.5, tmin=0.04, tmax=13.5)    

    pl.plot( z, gndmod, color='k', ls='-', marker=' ', lw=2 )
    pl.plot( z, hstmod, color='b', ls='--', marker=' ', lw=2 )
    pl.plot( z, ccmod, color='g', ls='--', marker=' ', lw=2 )
    pl.plot( z, r13mod, color='m', ls='--', marker=' ', lw=2 )
