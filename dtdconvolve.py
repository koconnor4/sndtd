import numpy as np

def snrate(z, DTDmodel='t-1+p', SFRmodel = 'B13', snclass='TN', 
           zFirst = 20,  H0=70, Om=0.3, Ode=0.7, tstep=0.01,
           normalize=False, debug=False, **kwargs):
    """Compute the SN1a rate as a function of redshift z.
    Returns a volumetric rate in units of SNuVol = [10-4 yr-1 Mpc-3]

    The available models for SNR(z) are :
       Model      parameters        description
    ----------------------------------------------------       
    'constant'     [rate]        no change with z 
    'sfr'          [eta]         prop. to star form'n rate 
    '2comp'        [alpha=eta,beta=tau]  Scann.&Bildsten 2-component model 
    'WD+He'        [eta]         short delay time, from Wang:2009a
    'WD+MS'        [eta]         interm. delay time, from Wang:2009b
    'WD+MSZ1'      [eta]         from Meng:2009, low metallicity
    'WD+MSZ2'      [eta]         from Meng:2009, high metallicity
    'WD+RG'        [eta]         long delay time, from Wang:2009b
    'SD'           [eta]         all three Single Degen. models from Wang:2009a,b
    'lognormal'    [eta,tau]     log-normal from Forster:2004
    't-1'          [eta,tau]     t-1 with first explosion at t=tau Gyr
    't-1+p'        [eta,tau]     t-1 with tau=fraction of SNe that are prompt

    The available models for the SFR(z) are : 
    'madau' : sum of exponentials (used by S04)
              parameter sets 'corrected' or 'uncorrected' 
              for dust extinction
    'cole'  : polynomial in z
              coleH06a uses modified Salpeter IMF
              coleH06b uses Baldry & Glazebrook IMF
              coleL08 uses the Li:2008 modified fit
    'behroozi' : polynomial in z 
              behroozi uses the post-2006 SFR(z) data
              behrooziHB06 uses Hopkins&Beacom 2006 data

    Normalize=True pins the curve to the volumetric SNIa 
    rate from SNLS at z=0.55:  
         SNR  =  (0.48 +/- 0.06)  * 10^-4 SNe yr-1 Mpc-3

    tstep : the time sampling step size in Gyr, for evaluating 
       the convolution integral.  Should be small enough to 
       accurately sample the structure of the DTD.
    """
    #from pytools import cosmo
    import astropy.units as u
    from astropy.cosmology import FlatLambdaCDM,z_at_value
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
    #import exceptions
    from scipy.interpolate  import interp1d
    import sfrate
    if debug: import pdb; pdb.set_trace()


    # Create a time array, in units of Gyr, sampling the age of the 
    # universe from the epoch of first star formation to the present
    #  * age of univ. at z of First Stars
    tFS = 13.77 - cosmo.lookback_time(zFirst).value
    #tFS = cosmo.agez(zFirst, H0=H0, Om=Om, Ode=Ode, unit='Gyr') 
    #  * age of universe today
    tNow = 13.77
    #tNow = cosmo.agez( 0, H0=H0, Om=Om, Ode=Ode, unit='Gyr') 
    t = np.arange( tFS+tstep/2., tNow, tstep )

    # Fill the arrays for the Star formation history and the DTD 
    # Note that we evaluate the SFR at times in the t array, 
    # relative to t=0 at the big bang. We evaluate the DTD at 
    # times t-tFS, relative to t=0 at the epoch of first stars.
    SFR = sfrate.sfrate( t, model=SFRmodel )
    if DTDmodel == 't-1+p' :
        # kwargs : eta, f0, t0, tmin, tmax
        DTD = dtd_t1_fprompt( t-tFS, **kwargs )
    else : 
        # Various models from Binary evolution simulations
        # kwargs : eta, tau
        DTD = delaytimedist( t-tFS, model=DTDmodel, **kwargs )
    # Do a discrete convolution of the SFR and DTD
    SNRt = tstep * np.convolve( SFR, DTD, mode='full')[:len(t)]
    # Set up an interpolator function to evaluate the SNR at any 
    # time t, relative to the big bang
    tz = 13.77 - cosmo.lookback_time(z).value
    #tz = cosmo.agez( z, H0=H0, Om=Om, Ode=Ode, unit='Gyr') 
    getSNRt = interp1d( t, SNRt, kind='linear', bounds_error=False, fill_value=0 )
    if normalize:
        #if snclass.startswith('TN'): snr0,z0 = 0.45e-4, 0.45
        if snclass.startswith('TN'): snr0,z0 = 0.67e-4, 0.79
        elif snclass=='Ia+': snr0,z0 = 0.2 * 0.45e-4, 0.45
        elif snclass=='Ia-': snr0,z0 = 0.16 * 0.45e-4, 0.45
        elif snclass=='Ia':  snr0,z0 = 0.64 * 0.45e-4, 0.45
        elif snclass=='CC':  snr0,z0 = 1.92e-4, 0.26
        tz0 = 13.77 - cosmo.lookback_time(z0).value
        #tz0 = cosmo.agez( z0, H0=H0, Om=Om, Ode=Ode, unit='Gyr') 
        snr0mod = getSNRt( tz0 )
        normfactor = snr0 / snr0mod
    else : 
        normfactor = 1

    # interpolate to the times (i.e. redshift values) of interest, 
    # apply the normalization, and multiply by 10^4 to cast into
    # SNuVol units
    SNRz = getSNRt( tz ) * normfactor * 1e4
    return( SNRz )



def old_snrate():        
    # Convert the input redshift array to time, in units of Gyr 
    #   Note that we start by sorting z in place, then pass it to the cosmo.agez 
    #   function with the [::-1] slicing, which reverses the sorted z array. Thus
    #   we get time sequences that start at t=0, with the big bang.
    #   The resulting SNR computed  below will therefore be reversed relative 
    #   to the input z array, so we will return it with [::-1] slicing.
    if np.iterable( z ) : 
        z.sort()
        z = z[::-1]
    #  * time since the Big Bang in Gyr:
    tBB = cosmo.agez( z, H0=H0, Om=Om, Ode=Ode, unit='Gyr')
    #  * time since the formation of the first stars
    t = tBB - tFS
   
    # Evaluate SFR(t) at each redshift (time) step:
    if np.iterable(SFR) : 
        if len( SFR ) != len(z) : 
            raise exceptions.RuntimeError("The SFR array must be the same length as the z array!")

    # evaluate the DTD model at each redshift (time) step:
    if model == 'sfr' : 
        # SNR follows the SFR directly 
        SNR = normfactor * eta * SFR

    elif model == '2comp' :
        #  Scan. & Bil. 2005 two-component model
        #  Neill:2006   A=1.4e-14   B=8.0e-4
        #  component A is slow, prop. to total stellar Mass
        #  (i.e. the integrated Star Formation History)
        #  component B is prompt, prop. to the instantaneous SFR
        if ('A' not in kwargs) or ('B' not in kwargs) : 
            raise exceptions.RuntimeError( '2comp model requires keyword arguments A and B')
        A = kwargs['A']
        B = kwargs['B']

        # Total stellar mass at each redshift is the integral of the SFR(t)
        # curve from t=tFS to t.  We adopt the user-provided redshift array 
        # to define the redshift spacing, and add on Mtot0, representing 
        # the total mass of stars formed prior to the highest redshift 
        # in the user's z array
        if z.max() < zFirst-0.01 : 
            z0 = np.arange( z.max(), zFirst, 0.01 )
            tBB0 = cosmo.agez( z0[::-1], H0=H0, Om=Om, Ode=Ode, unit='Gyr')
            SFR0 = sfrate( tBB0, model=SFRmodel )
            Mtot0 = (SFR0 * np.append(0,np.diff(tBB0))).sum()
        Mtot = SFR.cumsum() * np.append( 0, np.diff(t) ) + Mtot0
        SNR = normfactor * ( A * Mtot + B * SFR )

    elif model == 't-1+p' :
        # kwargs : eta, f0, t0, tmin, tmax
        import pdb; pdb.set_trace()
        dtd = dtd_t1_fprompt( t, **kwargs )
        SNR = normfactor * np.convolve( SFR, dtd, mode='full' )

    # Some models from Binary evolution simulations
    else : 
        dtd = delaytimedist( t, model=model, eta=eta, tau=tau )
        SNR = normfactor * np.convolve( SFR, dtd, mode='full' )

    return( SNR[::-1] )
        
foo = """
model.startswith('SD'): 
        # single degenerate model includes 2 or 3 SD progenitor tracks
        # The rate from each one is evaluated separately, then
        # scaled and added together
        snr, abserr = 0, 0 
        #eta = parameters[0]
        if model=='SD':
            # Wang:2009 SD cocktail
            modlist = ['WD+He','WD+MS','WD+RG']
            if iterable(eta) : 
                scalelist = eta
            else : 
                scalelist = [0.1*eta,0.1*eta,1e-4*eta]

        elif model=='SD2':
            # Wang:2009 just He and RG
            modlist = ['WD+MS','WD+RG']
            scalelist = [0.1*eta,1e-4*eta]

        elif model in ['SDZ','SD(Z)','SD-Z'] : 
            # Kobayashi:2009 SD cocktail
            modlist = ['WD+MS.Z','WD+RG.Z']
            scalelist = [6*eta,1*eta]

        elif model in ['SDRGZ','SDRG(Z)','SDRG-Z'] : 
            # Kobayashi:2009 no MS
            modlist = ['WD+RG.Z']
            scalelist = [eta]

        else : 
            raise exceptions.RuntimeError(
                "SNR Model %s is unknown."%model )
            
        for i in range(len( modlist )) : 
            modeli,scalei = modlist[i], scalelist[i]
            # Metallicity (rel. to solar) as a function of redshift.
            # Only used in the Kobayashi SD models
            # Zrelsun = 1.5 * 10**(-gammaZ*z) 
            subgral,abserr = scint.quad( snrgrand, tF, tz, args=(tz, modeli, eta[i], tau, SFR ) )
            if subgral : 
                if abserr/subgral > 0.01 : 
                    print( "Warning1 : large errors (%i percent) in Num. integr'n of SNR for %s."%((abserr/subgral)*100,model)) 
            snr += scalei*subgral 

    elif model.startswith('WD'): 
        # A single progenitor track 
        # e.g.  WD+MS or WD+RG

        # eta sets the scale of the DTD
        #if len(parameters) : eta = parameters[0]
        ## Metallicity (rel. to solar) as a function of redshift.
        ## Only used in the Kobayashi SD models
        #Zrelsun = 1.5 * 10**(-gammaZ*z) 

        # evaluate the SNR convolution integral, fixing eta=1
        snr,abserr = scint.quad( snrgrand, tF, tz, args=(tz, model, eta, tau, SFR ) )
        if snr : 
            if abserr/snr > 0.01 : 
                print( "Warning2 : large errors (%i percent) in Num. integr'n of SNR for %s."%((abserr/snr)*100,model)) 
        # snr *= eta

    elif model in ['exp','lognormal','t-1','gauss','t-1+p'] : 
        # evaluate the SNR convolution integral
        snr,abserr = scint.quad( snrgrand, tF, tz, args=(tz, model, eta, tau, SFR ) )
        if snr : 
            if abserr/snr > 0.01 : 
                print( "Warning3 : large errors (%i pct.) in Num. integr'n of SNR."%(abserr/snr*100) )

    else : 
        raise exceptions.RuntimeError(
            "%s is not a recognized SNRate model"%model )

    snr = snr * normfactor
        
    return(snr)
"""

def getNormFactor( SNRinterp, snclass='TN' ):
    """ get the normalization factor to adjust a SN rate curve 
    to match observations at z~0.5.
    For TNSN, we normalize to the SNLS rate at z=0.45
    using the Li:2001 peculiarity fractions
    For CCSN, we normalize to the Cappellaro rate at z=0.26
    """
    if snclass.startswith('TN'): snr0,z0 = 0.45e-4, 0.45
    elif snclass=='Ia+': snr0,z0 = 0.2 * 0.45e-4, 0.45
    elif snclass=='Ia-': snr0,z0 = 0.16 * 0.45e-4, 0.45
    elif snclass=='Ia':  snr0,z0 = 0.64 * 0.45e-4, 0.45
    elif snclass=='CC':  snr0,z0 = 1.92e-4, 0.26
    snr0mod = snrate( z0, snclass=snclass, DTDmodel=DTDmodel, 
                      SFRmodel=SFRmodel, zFirst=zFirst,  
                      H0=H0, Om=Om, Ode=Ode, 
                      normalize=False, **kwargs)
    return( snr0 / snr0mod )
    

def snrgrand(tp, t, model='SD', eta=None, tau=None, SFR='B13' ):
    """
    The SN1a rate is the convolution of the star
    formation rate SFR(t) with the delay time 
    function:

    SNR(t) =  int(tF,t){ SFR(t') x Phi(t-t') dt' }

    This function defines the integrand of the 
    convolution integral.  

    SFR may be a string, specifying the model 
    and paramerter set to use, or it may be an 
    interpolating function for SFR(t). 

    All times are in units of Gyr, and are relative 
    to the start of star formation in the universe,
    NOT relative to the beginning of the universe.
    """
    from numpy import iterable
    import sfrate

    # the delay time variable : t-t'
    td = t - tp

    # evaluate the delay time distribution at time td 
    phi = delaytimedist( td, model=model, eta=eta, tau=tau )

    # evaluate the instantaneous Star formation rate at time t'
    if isinstance( SFR, str) : 
        SFRnow = sfrate.sfrate( tp, model=SFR )
    else : 
        SFRnow = SFR( tp ) 

    # return value is the evaluated integrand 
    return(  SFRnow * phi )




def dtd_t1_fprompt( t, eta=1e-13, f0=0.2, t0=0.5, tmin=0.04, tmax=13.5  ):
    """ The delay time distribution, assuming 
    a t-1 power law that begins at t0 (~500 Myr). 
    For t<t0, a prompt component is added as a 
    constant rate from t=tmin (the min time needed 
    to become a WD) to t=t0.
    The constant component is scaled such that the 
    total fraction of SNe that are prompt is equal 
    to the user parameter f0.

    INPUTS: 
      t   : an input array of time since 
            onset of star formation [Gyr]
      eta : scaling factor [SNIa yr-1 Msun-1]
      f0  : the fraction of all SNe that are 'prompt'
             i.e. they explode before t=t0
      t0  : start of the t-1 distribution
      tmin: minimum time needed to make the first SNIa
            (default of 40 Myr is basically a guess)
      tmax: maximum time allowed for the last SNIa 
            (13.5 Gyr ~ age of the universe)
    """
    # Compute the constant factor K :
    # fprompt = Nprompt / (Nlate + Nprompt) 
    #         = K*(t0-tmin) / ( eta*ln(tmax/t0) + K*(t0-tmin) )
    K = ( f0/(1-f0) ) * eta * np.log(tmax/t0) / (t0-tmin)
    phi = np.where( t>t0,  eta / t, np.where( t>tmin, K, 0) )
    return( phi )


def delaytimedist( t, model, eta=None, tau=None, Z=None  ):
    """ the delay time distribution (Phi) 
    t : an input array of time since onset of star formation, in units of Gyr. 
    model : a string giving the name of the DTD model (see below)
    eta : efficiency factor  [ SNIa yr-1 Msun-1 ]
    tau : characteristic time [ Gyr ]

    Model choices : 
       'WD+He','WD+MS','WD+RG' : 
          individual single degenerate models from Bo Wang 2009
           parameters = [ eta, ]

       'SD' :
           composite single degenerate model, all of the above from Wang:2009.
           parameters = [ etaHe, etaMS, etaRG, ]

       't-1' = 'DD' = 'WD+WD' = 'ruiter' :
           Double Degenerate model from Ruiter:2012
           parameters = [ eta, ]

       'WD+MS.Z1' (alphaCE=3.0, Z=0.001)
       'WD+MS.Z2' (alphaCE=1.0, Z=0.03)
          fixed-metallicity models from Meng:2009
           parameters = [ eta, ]

       'WD+RG.Z' (long delay time)
       'WD+MS.Z' (short delay time)
          Metallicity dependent models from Kobayashi et al. 2009 
          t1 and slope vary with metallicity Z 
          parameters = [ eta, Zrelsun ]

       't-1' = 'DD' = 'WD+WD' = 'ruiter' :
           Double Degenerate model from Ruiter:2012
           parameters = [ eta, ]

       Other models: 'exp','gauss','2comp','lognormal'
    """
    from numpy import exp, sqrt, pi, log, iterable, array, \
        append, zeros, where
    import exceptions 
    if not iterable( t ): t = array([t]) 
    # parameters = [1,1]
    parameters = [ eta, tau ]

    # exponential model:
    if model=='exp':
        tau = parameters[0]
        return( exp(-t/tau)/tau )

    elif model in ['Ruiter','ruiter']:
        # Double Degenerate model mimicking  Ruiter:2012. 
        # a simple power law, truncated on the low end
        if eta == None : eta = 1e-13 # SNIa yr-1 Msun-1 (see Ruiter:2012, Fig7)
        t0    = 0.09   # time of first SN explosion [Gyr]
        t1    = 0.150  # time for start of t-1 DTD [Gyr]
        phi = where( t>t1, eta / t, where(t>t0, 0.0003*eta/t**4, 0) )

    elif model in ['t-1','DD','WD+WD']:
        # Double Degenerate model 
        # a t-1 power law, truncated on the low end at t=tau
        if eta == None : eta = 1e-13 # SNIa yr-1 Msun-1 (see Ruiter:2012, Fig7)
        if tau == None : tau = 0.1  # time of first SN explosion [Gyr]
        phi = where( t>tau, eta / t, 0 )

    elif model == 't-1+p':
        # a t-1 power law, truncated at t0=500 Myr, 
        # but with a prompt component with t<500 Myr,
        # added as a constant rate from t=tmin (the min
        # time needed to become a WD) to t=t0=500 Myr.
        # The constant component is scaled by a factor K
        # such that the total fraction of SNe that are 
        # prompt is equal to the user parameter tau

        # Compute the constant factor K :
        # fprompt = Nprompt / (Nlate + Nprompt) 
        #         = K*(t0-tmin) / ( eta*ln(tmax/t0) + K*(t0-tmin) )
        if eta == None : eta = 1e-13 # SNIa yr-1 Msun-1 (see Ruiter:2012, Fig7)
        if tau == None : tau = 0.1  # fraction of SNe that are prompt
        fprompt = tau
        t0 = 0.5  # 500 Myr : start of the t-1 distribution
        tmin = 0.01 # 10 Myr : time for an 8 Msun star to become a WD
        tmax = 13.3 # maximum age of a WD in the present universe
        K = ( fprompt/(1-fprompt) ) * eta * np.log(tmax/t0) / (t0-tmin)
        phi = where( t>t0, eta / t, K )

    elif model=='gauss':
        # A single-parameter gaussian delay-time dist'n model, 
        # following Dahlen+ 2008. 
        # For D08, using the Giavalisco+ SFR(z), they found 
        # the best-fit tau = 3.4.
        # When using the Behroozi+ 2012 SFR(z) a better fit is
        # eta = 0.0009, tau=2.2
        if eta==None : eta = 0.0009  # normalization factor
        if tau==None : tau=2.2  # peak of DTD gaussian [Gyr]
        t1    = 0.045    # time of first SN explosion [Gyr]
        sigma = 0.2*tau  # width of DTD gaussian [Gyr]
        phi = where( t>t1, eta * exp( -(t - tau)**2/(2*sigma**2)), 0 )
        
    # Single Degenerate models from Wang:2009a,b
    # for eta values, see Wang et al 2010 Sec.3 (arXiv:1002.4229)
    elif model == 'WD+He' : 
        # WD+He : short delay time model from Wang:2009a
        if eta==None : eta   = 1e-15 # SNIa yr-1 Msun-1
        t1    = 0.045  # time of first SN explosion [Gyr]
        sigma = 0.012   # width of DTD gaussian [Gyr]
        tau   = 0.079   # peak of DTD gaussian [Gyr]
        phi = where( t>t1, eta * exp( -(t - tau)**2/(2*sigma**2)), 0 )

    elif model == 'WD+MS' : 
        # WD+MS : intermediate delay time model from Wang:2009b
        if eta==None : eta = 1e-15 # SNIa yr-1 Msun-1
        t1    = 0.3   # time of first SN explosion [Gyr]
        sigma = 0.08    # width of DTD gaussian [Gyr]
        tau   = 0.5   # peak of DTD gaussian [Gyr]
        phi = where( t<tau, where( t<t1, 0, eta * exp( -(t-tau)**2/(2*sigma**2)) ), eta * exp( -2*(t-tau) ) )

    elif model == 'WD+RG' : 
        # WD+RG : long delay time model from Wang:2009b
        if eta==None : eta   = 1e-18 # SNIa yr-1 Msun-1
        t1    = 3.55  # time of first SN explosion [Gyr]
        sigma = 3.5   # width of DTD gaussian [Gyr]
        tau   = 5.62  # peak of DTD gaussian [Gyr]
        phi = where( t>t1, eta * exp( -(t - tau)**2/(2*sigma**2)), 0 )

    elif model == 'SD' : 
        # Composite SD model from Wang:2009a,b  (WD+He,MS,RG)
        # for eta values, see Wang et al 2010 Sec.3 (arXiv:1002.4229)
        if iterable(eta) : etaHe, etaMS, etaRG = eta
        elif eta==None : etaHe, etaMS, etaRG = [ 1e-14, 1e-14, 1e-16 ]
        else : etaHe, etaMS, etaRG = eta, eta, eta
        phiHe = delaytimedist( t, 'WD+He', etaHe )
        phiMS = delaytimedist( t, 'WD+MS', etaMS )
        phiRG = delaytimedist( t, 'WD+RG', etaRG )
        phi = phiHe + phiMS + phiRG

    elif model == 'WD+MS.Z1' : 
        # intermediate delay time model from Meng:2009
        #  assuming alpha_CE = 3.0, Z=0.001
        if eta==None : eta   = 1
        t1    = 0.8    # time of first SN explosion [Gyr]
        sigma = 0.2    # width of DTD gaussian [Gyr]
        tau   = 1.0   # peak of DTD gaussian [Gyr]
        phi = where( t>t1, 
                     eta * exp( -(t - tau)**2/(2*sigma**2)), 
                     0 )

    elif model == 'WD+MS.Z2' : 
        # intermediate delay time model from Meng:2009
        #  assuming alpha_CE = 1.0, Z=0.03
        if eta==None : eta   = 1
        t1    = 0.32   # time of first SN explosion [Gyr]
        sigma = 0.15    # width of DTD gaussian [Gyr]
        tau   = 0.63   # peak of DTD gaussian [Gyr]
        phi = where( t>t1, 
                     eta * exp( -(t - tau)**2/(2*sigma**2)), 
                     0 )

    elif model == 'WD+RG.Z' : 
        # WD+RG : long delay time model from Kobayashi:2009
        #  t1 and slope vary with metallicity Z 
        eta   = parameters[0]
        Z     = parameters[1]  # Metallicity (rel. to solar)
        # metallicity-dependant parameters:
        t1 = 3.1 * exp( -4 * Z ) + 0.25 # time of first SNIa
        tau = 0.15 * Z - 1.54  # slope of power-law
        t2    = 10.1   # time of last SN explosion [Gyr]
        phi = where( t>t2, 0,  
                     where( t<t1, 0, 
                            eta * 1e-3 * t**(tau) ) )
    elif model == 'WD+MS.Z' : 
        # WD+MS : short delay time model from Kobayashi:2009
        #  t1 and tau vary with metallicity Z 
        eta   = parameters[0]
        Z     = parameters[1] # Metallicity (rel. to solar)
        # metallicity-dependant parameters:
        t1 = 100 * exp( -8.5 * Z**0.2 ) + 0.06 # time of first SNIa
        tau = 0.15 * Z - 1.54  # slope of power-law
        t2    = 1.4   # time of last SN explosion [Gyr]
        phi = where( t>t2, 0,  
                     where( t<t1, 0, 
                            eta * 1e-3 * t**(tau) ) )

    elif model in ['2comp','SFR'] :
        # Scan. & Bil. 2005 two-component model
        # Here we integrate the SFR up to the given redshift
        # to get the first component, which is proportional to
        # the total mass of stars.
        phi = 1
    elif model=='lognormal' : 
        # log-normal model from Forster:2006
        if eta == None : eta = 1 
        if tau == None : tau = 3 # Gyr
        sigma = 0.5*tau
        phi = eta * (1 / sqrt(2*pi*log(sigma)**2) ) * \
            exp( -(log(t/tau)**2/(2*log(sigma)**2)) ) / t
    else : 
        raise exceptions.RuntimeError(
            "%s is not a recognized SNRate model"%model )

    if iterable(phi) :
        if len(phi) == 1 : return( phi[0]) 
    return( phi ) 


