def behroozi2012( x, abscissa='redshift', fitto='B13'):
    """alias for behroozi2013()"""
    behroozi2013( x, abscissa=abscissa, fitto=fitto)

def B13( x, abscissa='redshift', fitto='B13'):
    """alias for behroozi2013()"""
    behroozi2013( x, abscissa=abscissa, fitto=fitto)

def behroozi2013( x, abscissa='redshift', fitto='B13'):
    """ the star formation rate as a function of redshift
    in units of  [ Msun yr-1 Mpc-3 ], following Behroozi et al 2013b

    SFR(z) = C / [ 10^A(z-z0) + 10^B(z-z0) ]  Msun-1 yr-1 Mpc-3 h3

    abscissa : redshift = the x array is redshift
               time or age = the x array is time (age of the universe i Gyr)

    fitto :  B12,B13 = Behroozi 2013b recent compilation
             HB06 = Hopkins & Beacom 2006 older compilation
             B13high = upper systematics limit
             B13low = lower systematics limit
             B13sharp = get a sharp peak, lower systematics limit for z>1.5
             B13flat = get a flatter curve, upper systematics limit for z>1.5
    """
    if fitto == 'HB06' : 
        # Behroozi:2013, fitting old Hopkins+Beacom:2006 data
        (z0,a,b,c) = (0.840,-1.311,0.085,0.143)
    elif fitto == 'B13high' : 
        # Behroozi:2013 compilation, taking upward systematic shift everywhere
        (z0,a,b,c) = (1.293, -0.980,  0.219,  0.269 )
    elif fitto == 'B13low' : 
        # Behroozi:2013 compilation, taking downward systematic shift everywhere
        (z0,a,b,c) = (1.176, -0.980,  0.254,  0.116)
    elif fitto == 'B13sharp' : 
        # Behroozi:2013 compilation, taking downward shift for z>1.5
        (z0,a,b,c) = (1.121, -1.054,  0.287,  0.157)
    elif fitto == 'B13flat' : 
        # Behroozi:2013 compilation, taking upward shift for z>1.5 
        (z0,a,b,c) = (1.464, -0.867,  0.203,  0.218)
    else : 
        # baseline fit to the Behroozi:2013 compilation of recent data w/ revised dust corrections
        (z0,a,b,c) = (1.243,-0.997,0.241,0.180)
        
    if abscissa in ['redshift','z'] : 
        z = x 
    elif abscissa in ['time','age'] : 
        # the Behroozi parameterization is in redshift space,
        # and we are working in age-of-universe space. 
        # To do the conversion, we assume a flat LambdaCDM 
        # cosmology, and use the formula from Peebles:1993
        # p.315, solved for z. 
        from cosmo import zfromt
        z = zfromt( x, H0=70 ) 
    else : 
        raise exceptions.RuntimeError( "abscissa must be 'redshift' or 'time'")

    # now compute the SFR density
    sfr = c / ( 10**(a*(z-z0)) + 10**(b*(z-z0)))
    return( sfr )

def yuksel2008( x, abscissa='redshift', fitto='H11'):
    """ the star formation rate as a function of redshift,
    following Yuksel et al 2008

    SFR(z) = sfr0 * [(1+z)**(a*eta) + ( (1+z)/B )**(b*eta) + ( (1+z)/C )**(c*eta) ] ** (1/eta)
      B = (1+z1)**(1-(a/b))
      C = (1+z1)**((b-a)/c) * (1+z2)**(1-b/c)

    abscissa : redshift = the x array is redshift
               time or age = the x array is time (age of the universe)

    fitto :  HB06 = Hopkins & Beacom 2006 dust-corrected compilation
    """

    # Parameterization from Horiuchi et al 2011
    sfr0 =  0.016 * (73/70.)**3
    a,b,c = 3.4, -0.3, -3.5
    z1, z2 = 1., 4.
    eta = -10. 
    B = (1+z1)**(1-(a/b))
    C = (1+z1)**((b-a)/c) * (1+z2)**(1-(b/c))

    if abscissa in ['redshift','z'] : 
        z = x 
    elif abscissa in ['time','age'] : 
        # this parameterization is in redshift space,
        # and we are working in age-of-universe space. 
        # To do the conversion, we assume a flat LambdaCDM 
        # cosmology, and use the formula from Peebles:1993
        # p.315, solved for z. 
        from cosmo import zfromt
        z = zfromt( x, H0=70 ) 
    else : 
        raise exceptions.RuntimeError( "abscissa must be 'redshift' or 'time'")

    sfr = sfr0 * ( (1+z)**(a*eta) + ( (1+z)/B )**(b*eta) + ( (1+z)/C )**(c*eta) ) ** (1/eta)

    return( sfr )



def sfrate( t, model='B13' ):
    """ the star formation rate as a function of time.  Time t is in Gyr
    since the beginning of the universe.
    """
    from numpy import exp
    #from cosmo import agez # requires scipy.integrate
    #from cosmo import agezAlt as agez
    from cosmo import zfromt
    import exceptions

    # Cole:2001  SFH = (a+bz)/(1+(z/c)**d) h Msun yr-1 Mpc-3
    #   (a,b,c,d) = (0.0166,0.1848,1.9474,2.6316)  (corrected for dust)
    #   (a,b,c,d) = (0.0,0.0798,1.658,3.105)  (uncorrected for dust)
    # Hopkins & Beacom, 2006  (using the Cole:2001 form)
    #   (a,b,c,d) = (0.017,0.13,3.3,5.3)  (Modified Salpeter)
    #   (a,b,c,d) = (0.0118,0.08,3.3,5.2)  (Baldry & Glazebrook)
    # Li,L.-X. 2008 (using the Cole:2001 form)
    #   (a,b,c,d) = (0.0157,0.118,3.23,4.66)

    # redshift for the onset of star formation
    # z = 25   # reference??? Maio:2009, range=12-40
                  
    # The Star Formation Rate as a function of t, the age of the 
    # universe since star formation began.
    
    if model.startswith('B13') :
        sfr = behroozi2013( t, abscissa='time', fitto=model)

    elif model.startswith( 'madau' ) : 
        # following Strolger et al. 2004, using parameters a,b,c,d for 
        # reddening corrected SFR(t) formula from Giavalisco:2004
        #  TODO  : update params with fit to Bouwens:2009 data
        if model.endswith('??') : 
            a,b,c,d,t0 = 0.182, 1.26, 1.865, 0.071, 13.66
        elif model.endswith('S04') : 
            # parameterization of Strolger:2004
            a,b,c,d,t0 = 0.021, 2.12, 1.69, 0.207, 13.66
        else :
            print("%s not a recognized SF model."%model)
            return()
        sfr = a*(t**b*exp(-t/c) + d*exp(d*(t-t0)/c))

    elif model.startswith('cole') :
        # Cole:2001 form:  sfr = (a+b*z)/(1+(z/c)**d)         
        if model=='cole' : 
            # Cole:2001, dust-corrected
            (a,b,c,d) = (0.0166,0.1848,1.9474,2.6316) 
        elif model.endswith('H06a') : 
            # Hopkins:2006, using Modified SalA IMF
            (a,b,c,d) = (0.017,0.13,3.3,5.3) 
        elif model.endswith('H06b') : 
            # Hopkins:2006, using Baldry & Glazebrook IMF
            (a,b,c,d) = (0.0118,0.08,3.3,5.2) 
        elif model.endswith('L08') : 
            # Li:2008, using Salpeter IMF
            (a,b,c,d) = (0.0157,0.118,3.23,4.66) 

        elif model.endswith('mid') : 
            # Li:2008, with hi-z data
            (a,b,c,d) = (0.0157, 0.118, 3.23, 4.66)

        elif model.endswith('steep') : 
            #(a,b,c,d) = (0.0157,0.118,2.5,5.06)    # same as mid for z<2
            (a,b,c,d) = (0.006,0.08,2.7,4.8)

        elif model.endswith('flat') : 
            #(a,b,c,d) = (0.013,0.13,3.63,2.5)   # same as mid for z<2
            (a,b,c,d) = (0.03,0.15,3.43,3.2)   

        else : 
            print("%s not a recognized  SF model."%model)
            return()
        # the Cole parameterization is in redshift space,
        # and we are working in age-of-universe space. 
        # To do the conversion, we assume a flat LambdaCDM 
        # cosmology, and use the formula from Peebles:1993
        # p.315, solved for z. 
        from cosmo import zfromt
        z = zfromt( t, H0=70 ) 
        # now compute the SFR density
        sfr = (a+b*z)/(1+(z/c)**d)         
    else : 
        print("%s not a recognized  SF model."%model)
        return()
    return( sfr ) 
