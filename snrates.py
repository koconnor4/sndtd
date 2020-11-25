# 2013.05.22 S.Rodney
# Computing the CANDELS rates v1 for the paper
import numpy as np
from scipy.stats.distributions import poisson
#from hstsnpipe.tools import snana

from hstsnpipe.tools import cosmo
_zFirst = 20
_H0=70
_Om=0.3
_Ode=0.7

#  * age of univ. at z of First Stars  [Gyr]
_tFS = cosmo.agez(_zFirst, H0=_H0, Om=_Om, Ode=_Ode, unit='Gyr') 
#  * age of universe today [Gyr]
_tNow = cosmo.agez( 0, H0=_H0, Om=_Om, Ode=_Ode, unit='Gyr') 


def getNctrl( survey='HST', field='all', dustmodel='mid',
              Nsim=2000, decliners=False, showplots=False, 
              Nbins=5, dz=0.5, clobber=False, verbose=False, **kwargs):
    """ set up and run a simulation to compute the 'control count':
    the number of SNIa that would be detected in the given survey 
    field if the Ia rate were flat at 1 SNuVol = 1 SNIa 10-4 yr-1 Mpc-3 h3

    field : ['all','all1','all2','gsd','gsw','uds','cos','egs','gnd','gnw']
    dustmodel : ['all','high','mid','low']
    decliners : include SNe that peak before the survey starts? 

    To ~match Dahlen+ 2008 use dz=0.4 and Nbins=4 or 5
    """
    from hstsnpipe.tools.snana import SimTable
    from hstsnpipe.tools.snana.constants import SurveyDataDict
    from hstsnpipe.tools.snana.simulate import mkinput, mksimlib, dosim
    import os

    sndatadir = os.environ['SNDATA_ROOT']

    # scale the survey area by a large factor so that the control 
    # count (an integer) carries sufficient significant digits
    areascale = 1e6

    # set up the redshift bins for computing rates 
    # and for plotting detection efficiency histograms
    HbinsDetEff = np.arange(20,26.5,0.1)

    if decliners : 
        mjdrangelist = [ [-100,20], [-100,20],  [-70,20],  [-40,20],  [0,0]  ]
    else : 
        # When decliners are excluded, we drop the first epoch of every field,
        # and set the PKMJD range so that we simulate only SNe that are 
        # detectable as positive flux deviations relative to that template epoch
        # (see rates.controlTime for the derivation of these points)
        mjdrangelist = [ [-35,19], [-37,21],  [-19,22],  [-15,15],  [-5,5]  ]

    zrangelist = [ [ round( max(a,0.001),3), round( a+dz,3) ] for a in  \
                       np.linspace( 0, Nbins*dz, Nbins+1 )[:-1]  ]

    if Nbins == 6 : 
        mjdrangelist = [ [-35,19], [-37,21],  [-19,22], [-17,18], [-15,15],  [-5,5]  ]
    elif Nbins == 5 : 
        mjdrangelist = [ [-35,19], [-37,21],  [-19,22],  [-15,15],  [-5,5]  ]
    elif Nbins == 4 :
        mjdrangelist = [ [-35,19], [-37,21],  [-19,22],  [-15,15] ]
    elif Nbins==3 : 
        mjdrangelist = [ [-21,37], [-23,20], [-12,12] ]

    surveydata = SurveyDataDict[survey] 
    if field=='all': fieldlist = ['gsd','gsw','cos','uds','egs','gnd','gnw']
    else : fieldlist = [field]
    if dustmodel=='all': dustlist = ['high','mid','low']
    else : dustlist = [dustmodel]

    Nctrl = {}
    Nctrldz05 = {}
    zdz05 = {}
    DetEffH = {}
    DetEffz = {}
    for field in fieldlist : 
        simroot = '%s_%s_Nctrl'%(survey,field)
        simlibfile = '%s.simlib'%simroot

        for dustmodel in dustlist : 
            for zrange,mjdrange in zip( zrangelist, mjdrangelist ) :
                zbinsDetEff = np.arange(zrange[0],zrange[1]+0.05,0.05)
                zstr = 'z%.2f'%round(np.mean(zrange),2)

                # Set up the simulated PEAKMJD range, extending before/beyond
                # the survey dates according to the redshift-dependent  
                # ranges of detectability, as given in mjdrangelist
                if decliners : mjdlist = sorted( surveydata.EPOCHLIST[field] )
                else : mjdlist = sorted( surveydata.EPOCHLIST[field][1:] )
                mjdpk0 = min(mjdlist) + mjdrange[0]
                mjdpk1 = max(mjdlist) + mjdrange[1]

                # check for existing simulation products
                simname = '%s_dust%s_%s'%(simroot,dustmodel,zstr)
                simdatadir = os.path.join( sndatadir,'SIM/%s'%simname )
                simisdone = np.all([ os.path.isfile( os.path.join(simdatadir,'%s_%s.FITS'%(simname,sfx)) ) 
                                     for sfx in ['PHOT','HEAD'] ] )

                if not simisdone or clobber : 
                    if verbose>2: print(" simsdone=%s  clobber=%s ..."%(simisdone,clobber))
                    if verbose>1: print(" Running SNANA simulation for %s ..."%simname)

                    # make the simlib file
                    mjdlist,bandlist,etimelist = [],[],[]
                    for mjd in surveydata.EPOCHLIST[field] : 
                        for band,etime in zip(surveydata.SEARCHBANDS,surveydata.SEARCHETIMES) :
                            mjdlist.append( mjd ) 
                            bandlist.append( band )
                            etimelist.append( etime )
                            continue
                        continue
                    simlibfile = mksimlib( simlibfile, survey=survey, field=field,
                                           mjdlist=mjdlist, bandlist=bandlist, etimelist=etimelist,
                                           perfect=False, clobber=clobber, verbose=verbose )

                    # make the input file
                    inputfile = '%s.input'%simname 
                    inputfile = mkinput( simname,  inputfile=inputfile, simlibfile=simlibfile,
                                         survey=survey, field=field, simtype='Ia', Nsim=Nsim, 
                                         ratemodel='constant',       # SNR(z) : 'constant','powerlaw','flat'
                                         dustmodel=dustmodel,   # p(Av) or p(c) model : 'high','mid','low'
                                         applySearchEff=0,  # Keep all SNe, raising flags if not detected
                                         smear = True, # smear the magnitudes to reflect intrinsic variations ?
                                         clobber=clobber, verbose=verbose,
                                         GENRANGE_REDSHIFT= [ '%.3f %.3f'%(zrange[0],zrange[1]),'# range of simulated redshifts'],
                                         GENRANGE_PEAKMJD= [ '%.1f %.1f'%(mjdpk0,mjdpk1),'# range of simulated peak MJD dates'],
                                         GENFILTERS=['JH','# List of filters to simulate'],
                                         SOLID_ANGLE=[surveydata.SOLIDANGLE[field]*areascale,'# Field size in (scaled) steradians'],
                                         **kwargs )

                    # run the simulation 
                    dosim( inputfile, simname=simname, perfect=False, simargs='', verbose=verbose )

                # read in the simulation results
                sim = SimTable( simname, verbose=verbose )

                # Compute the Control Count : the number of  simulated detections        
                # First, normalize the distributions and scale by NSNSURVEY,
                # which is the number of SN events in the survey volume/time.
                # Also remove the area scaling factor. 
                idet = sim.DUMP['idet']
                Nctrlz = idet.sum() * float(sim.NSNSURVEY) / Nsim / areascale
                Nctrl['%s.%s.%s'%(field,dustmodel,zstr)] = Nctrlz

                # Compute the detection efficiency vs peak H mag
                NdetH,Hbins = np.histogram( sim.Hpk[idet], bins=HbinsDetEff )
                NsimH,Hbins = np.histogram( sim.Hpk, bins=HbinsDetEff )
                DetEffH['%s.%s.%s'%(field,dustmodel,zstr)] = np.ma.masked_invalid( NdetH / NsimH.astype(float) )
  
                # Compute the detection efficiency vs redshift
                Ndetz,zbins = np.histogram( sim.z[idet], bins=zbinsDetEff )
                Nsimz,zbins = np.histogram( sim.z, bins=zbinsDetEff )
                DetEffz['%s.%s.%s'%(field,dustmodel,zstr)] = np.ma.masked_invalid( Ndetz / Nsimz.astype(float) )

                # Record the ctrl count in bins of dz=0.05, 
                # for computing the rate with a rolling window function. 
                Nctrldz05['%s.%s.%s'%(field,dustmodel,zstr)] = Ndetz * float(sim.NSNSURVEY) / Nsim / areascale
                if verbose : print( "%.2f  %.2f   %.3f"%( zrange[0], zrange[1], Nctrlz ) )
            
    if showplots : 
        from matplotlib import pyplot as pl
        # Distribution of all SNe exploding in this volume
        # (the total number here depends on the user-defined Nsim=NGENTOT_LC)
        zbinsplot=np.arange(0,3.,0.1)
        Nall,zall = np.histogram( sim.DUMP['REDSHIFT'], bins=zbinsplot )
        Nallsurvey = Nall * normfactor
        Ndetsurvey = Ndet * normfactor
        pl.clf()
        pl.plot( zbinsplot[1:], Nallsurvey, drawstyle='steps-post', color='b', label=r'%i total SNe in survey volume \& time'%sim.NSNSURVEY )
        pl.plot( zbinsplot[1:], Ndetsurvey, drawstyle='steps-post', color='r', label='%i detectable SNe'%Ndetsurvey.sum() )
        pl.xlabel('redshift')
        pl.ylabel('Number of SNe')
        pl.legend( loc='upper left', frameon=False, numpoints=3, borderpad=0.2)

        printNctrl( Nctrl, dustlist, fieldlist, zrangelist )
    return( Nctrl) #, NctrlRoll, zRoll, DetEffH, DetEffz )


def printNctrl( Nctrl, dustlist, fieldlist, zrangelist ) :
    """ print a summary table of Nctrl values
    and print Nctrl arrays to be planted in ratetable.py """
    NctrlArrays = {}
    for dustmodel in dustlist : 
        print( "\n\n---- Dust = %s ----"%dustmodel )
        header = '   z       ' + '    '.join( fieldlist ) + '      ALL'
        table = ''
        nctrlarray = []
        nctrlarraystr = "Nctrl_%sdust = array([ "%( dustmodel )
        NctrlAllFieldsAllz = 0
        for zrange in zrangelist : 
            zstr = 'z%.2f'%round(np.mean(zrange),2)
            table += '\n%.1f-%.1f '%(tuple(zrange))
            NctrlAllFields = 0
            for field in fieldlist : 
                key = '%s.%s.%s'%(field,dustmodel,zstr) 
                table += ' %6.3f'%Nctrl[key]
                NctrlAllFields += Nctrl[key]
            table += '  %6.3f'%NctrlAllFields
            NctrlAllFieldsAllz += NctrlAllFields
            nctrlarray.append( NctrlAllFields )
            nctrlarraystr += '%.3f, '% NctrlAllFields
        table += '\n total: %.3f'%NctrlAllFieldsAllz
        nctrlarraystr += '])'
        print( header + table )
        print( nctrlarraystr )
        NctrlArrays[ dustmodel ] = np.array( nctrlarray )
    return( NctrlArrays )

def getNctrldz05(NctrlDict, dustlist, fieldlist, zrangelist ):
    """ compute the ctrl count in rolling windows of width dz
    NOT YET FUNCTIONAL
    """

    NctrlRoll = {}
    for dustmodel in dustlist : 
        print( "\n\n---- Dust = %s ----"%dustmodel )
        header = '   z      ' + '     '.join( fieldlist ) + '      ALL'
        table = ''
        NctrlAllFieldsAllz = 0
        #NctrlAllFieldsAllz_dz05 = []
        for zrange in zrangelist : 
            zstr = 'z%.2f'%round(np.mean(zrange),2)
            table += '\n%.1f-%.1f '%(tuple(zrange))
            NctrlAllFields = 0
            #NctrlAllFields_dz05 = []
            for field in fieldlist : 
                key = '%s.%s.%s'%(field,dustmodel,zstr) 
                table += ' %6.3f'%Nctrl[key]
                NctrlAllFields += Nctrl[key]
                #NctrlAllFields_dz05.append( Nctrldz05[key] )
            table += '  %6.3f'%NctrlAllFields
            NctrlAllFieldsAllz += NctrlAllFields
            NctrlAllFieldsAllz_dz05.append( np.sum( NctrlAllFields_dz05, axis=0 ) )
        table += '\n total: %.3f'%NctrlAllFieldsAllz
        print( header + table )


    NctrlAll_dz05 = np.ravel( NctrlAllFieldsAllz_dz05 )
    zdz05 = np.arange( np.round(zrangelist[0][0],2), zrangelist[-1][-1], 0.05 )+0.025
    zRoll = np.arange( np.round(np.mean(zrangelist[0]),2), zrangelist[-1][-1]-np.diff(zrangelist[-1])/2.+0.05, 0.05 )
    NctrlRoll[dustmodel] = []
    for zR in zRoll : 
        iinbin = np.where( (zdz05>=zR-(dz/2.)) & (zdz05<zR+(dz/2.)) )
        NctrlRoll[dustmodel].append( np.ravel(NctrlAllFieldsAllz_dz05)[iinbin].sum() )



def printRates( decliners=False, dz=0.5 ):
    """ count the SN detections, compute the rates, print them out
    keyword args passed to countDetections
    """

    # Control counts computed using snana.simulate.getNctrl 
    # TODO : option to re-run sims and re-compute the control counts 
    # 2013.07.11 : using updated detection efficiency curves 

    if dz==0.4 : 
        # 6 bins dz=0.4 to z=2.4
        # ---- Dust = mid ----
        #    z      gsd     gsw     cos     uds     egs     gnd     gnw      ALL
        # 0.0-0.4   0.724  0.046  0.232  0.244  0.250  0.758  0.109   2.363
        # 0.4-0.8   2.890  0.193  0.962  1.012  1.034  2.986  0.451   9.527
        # 0.8-1.2   3.739  0.181  0.902  0.947  0.972  3.825  0.426  10.992
        # 1.2-1.6   3.170  0.124  0.619  0.649  0.666  3.279  0.292   8.799
        # 1.6-2.0   1.929  0.030  0.148  0.156  0.160  1.990  0.070   4.480
        # 2.0-2.4   0.935  0.015  0.073  0.078  0.079  0.991  0.035   2.206
        #  total: 38.368
        NctrlLowDust  = np.array([ 2.367, 9.862, 12.999, 12.672, 6.451, 3.176 ])
        Nctrl         = np.array([ 2.363, 9.527, 10.992,  8.799, 4.480, 2.206 ])
        NctrlHighDust = np.array([ 2.344, 8.375,  7.624,  5.379, 2.737, 1.347 ])
        zbinedges = np.array([0,0.4,0.8,1.2,1.6,2.0,2.4])

    if dz==0.5 : 
        #    5 bins dz=0.5  to z=2.5
        # Recalculated 2013.10.29 with final areas and cadences from Tomas
        NctrlHighDust = np.array([ 4.064, 12.398, 9.108, 4.681, 1.453, ])
        Nctrl          = np.array([ 4.101, 14.105, 13.161, 7.675, 2.524, ])
        NctrlLowDust  = np.array([ 4.107, 14.636, 15.618, 11.147, 4.246, ])
        zbinedges = np.array([0,0.5,1.0,1.5,2.0,2.5])


    elif dz==0.6 : 
        # 4 bins dz=0.6 to z=2.4
        # ---- Dust = mid ----
        #    z      gsd     gsw     cos     uds     egs     gnd     gnw      ALL
        # 0.0-0.6   1.929  0.122  0.612  0.644  0.659  1.994  0.288   6.249
        # 0.6-1.2   5.543  0.366  1.829  1.908  1.952  5.687  0.849  18.135
        # 1.2-1.8   4.438  0.221  1.103  1.167  1.192  4.500  0.524  13.147
        # 1.8-2.4   1.900  0.077  0.384  0.403  0.414  1.871  0.180   5.230
        #  total: 42.761
        NctrlLowDust  = np.array([ 6.261, 18.773, 15.548, 7.532])
        Nctrl         = np.array([ 6.249, 18.135, 13.147, 5.230])
        NctrlHighDust = np.array([ 6.198, 15.941,  9.119, 3.197])
        zbinedges = np.array([0,0.6,1.2,1.8,2.4])

    zbins = zbinedges[:-1] + np.diff( zbinedges )/2.
    Nbins = len(zbins)

    # Systematic uncertainty estimates:
    dNctrlPlus = NctrlLowDust - Nctrl
    dNctrlMinus = Nctrl - NctrlHighDust

    Nobsz, dNobszSysPlus, dNobszSysMinus, Nobsh, dNobshSysPlus, dNobshSysMinus  = countDetections( zbinedges=zbinedges )

    if decliners : # compute the rate with decliners
        print("WARNING : decliners ctrl counts are out-of-date!!")
        return(0)
        # shortcut: with-decliners systematics based on no-decliner sets defined above
        Nctrl = np.array([2.670, 10.380, 13.323, 9.258])
        Nobs, dNobsSysPlus, dNobsSysMinus = countDetections( decliners=True, **kwargs )

    # compute the rate
    ratez = Nobsz / Nctrl
    rateh = Nobsh / Nctrl

    # Systematic uncertainty: shift the observed and control counts 
    # coherently and compare to the baseline rates
    dratezSysMinus = ratez - ( Nobsz - dNobszSysMinus ) / (Nctrl+dNctrlPlus)
    dratezSysPlus = ( Nobsz + dNobszSysPlus ) / (Nctrl-dNctrlMinus) - ratez

    dratehSysMinus = rateh - ( Nobsh - dNobshSysMinus ) / (Nctrl+dNctrlPlus)
    dratehSysPlus = ( Nobsh + dNobshSysPlus ) / (Nctrl-dNctrlMinus) - rateh

    # Statistical uncertainty: compute the 1-sigma Poisson confidence limits
    NobszMinus, NobszPlus = poissonLimits( Nobsz, 1 ) 
    dNobszStatMinus, dNobszStatPlus  = Nobsz-NobszMinus, NobszPlus-Nobsz
    dratezStatMinus = dNobszStatMinus / Nctrl
    dratezStatPlus = dNobszStatPlus / Nctrl

    NobshMinus, NobshPlus = poissonLimits( Nobsh, 1 ) 
    dNobshStatMinus, dNobshStatPlus  = Nobsh-NobshMinus, NobshPlus-Nobsh
    dratehStatMinus = dNobshStatMinus / Nctrl
    dratehStatPlus = dNobshStatPlus / Nctrl


    print( "%4s  %5s   %10s    %12s   %5s   %10s    %12s"%('z','SNRz','Stat.','Syst.','SNRh','Stat.','Syst.') )
    for i in range(len(zbins)) : 
        print( "%5.2f  %5.2f  +%5.2f -%5.2f   +%5.2f -%5.2f   %5.2f  +%5.2f -%5.2f   +%5.2f -%5.2f"%(
                zbins[i], 
                ratez[i], dratezStatPlus[i], dratezStatMinus[i], 
                dratezSysPlus[i], dratezSysMinus[i],
                rateh[i], dratehStatPlus[i], dratehStatMinus[i], 
                dratehSysPlus[i], dratehSysMinus[i] ))

    fmttuple = ('%.2f,'*len(ratez)).rstrip(',')
    fmtstr =  """
R13z.Ndet        = array([%s])
R13z.rate        = array([%s])
R13z.errstatplus = array([%s])
R13z.errstatminus= array([%s])
R13z.errsysplus  = array([%s])
R13z.errsysminus = array([%s])
"""%( fmttuple,fmttuple,fmttuple,fmttuple,fmttuple,fmttuple) 
    print( fmtstr % ( tuple(Nobsz) + tuple(ratez) + 
                      tuple(dratezStatPlus) + tuple(dratezStatMinus) + 
                      tuple(dratezSysPlus) + tuple(dratezSysMinus) ) )

    fmtstr =  """
R13h.Ndet        = array([%s])
R13h.rate        = array([%s])
R13h.errstatplus = array([%s])
R13h.errstatminus= array([%s])
R13h.errsysplus  = array([%s])
R13h.errsysminus = array([%s])
"""%( fmttuple,fmttuple,fmttuple,fmttuple,fmttuple,fmttuple) 
    print( fmtstr % ( tuple(Nobsh) + tuple(rateh) + 
                      tuple(dratehStatPlus) + tuple(dratehStatMinus) + 
                      tuple(dratehSysPlus) + tuple(dratehSysMinus) ) )


    print( "\n\n USING REDSHIFT DEPENDENT PRIOR")
    for i in range(len(zbins)) : 
        print( "$%.2f$ & %.2f & $^{+%.2f}_{%.2f}$ & $^{+%.2f}_{%.2f}$ & & %.2f & $^{+%.2f}_{%.2f}$ & & %.2f & $^{+%.2f}_{%.2f}$ & $^{+%.2f}_{%.2f}$  \\\\[0.7em]" % ( 
                zbins[i], 
                Nobsz[i], dNobszStatPlus[i], -dNobszStatMinus[i], 
                dNobszSysPlus[i], -dNobszSysMinus[i], 
                Nctrl[i], dNctrlPlus[i], -dNctrlMinus[i], 
                ratez[i], dratezStatPlus[i], -dratezStatMinus[i], 
                dratezSysPlus[i], -dratezSysMinus[i], 
                ) )

    print( "\n\n Change in Counts / Rates when USING HOST GALAXY PRIOR")
    for i in range(len(zbins)) : 
        print( "%.2f & %+.2f & %+.2f  & %i \\\\" % ( 
               zbins[i], Nobsh[i]-Nobsz[i], rateh[i]-ratez[i],
               100* (rateh[i]-ratez[i])/dratezSysPlus[i]  ) )


def countDetections( datfile='candels.csv', 
                     zbinedges=[0,0.5,1.0,1.5,2.0,2.5], 
                     decliners=False, Bgrade=False ):
    """ read in the .dat file (exported from the google doc as a .csv file) 
    and count up the number of SNIa detections as a function of redshift """
    from astropy.io import ascii
    import sys
    import os

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    if not os.path.isfile( datfile ) :
        datfile = os.path.join( os.path.dirname( thisfile ), datfile )
    if not os.path.isfile( datfile ) :
        print("Cannot find %s"%datfile)
        return(-1) 

    c1dat = ascii.read( datfile )

    name = c1dat['name']
    nickname = c1dat['nickname']
    zsn = c1dat['zSN']
    zsnerr = c1dat['dzSN']
    # confidence = c1dat['confidence']
    # decliner = c1dat['decliner']
    
    piaz = c1dat['PIaz']
    dpiazPlus = np.array( [float(dpia.split()[0]) for dpia in c1dat['dPIaz']] )
    dpiazMinus = np.array( [abs(float(dpia.split()[-1])) for dpia in c1dat['dPIaz']] )

    piah = c1dat['PIaHost']
    dpiahPlus = np.array( [float(dpia.split()[0]) for dpia in c1dat['dPIaHost']] )
    dpiahMinus = np.array( [abs(float(dpia.split()[-1])) for dpia in c1dat['dPIaHost']] )

    Ndetz = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]
    dNdetzPlus = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]
    dNdetzMinus = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]
    Ndeth = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]
    dNdethPlus = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]
    dNdethMinus = [ [] for izbin in range( 0,len(zbinedges)-1 ) ]

    for isn in range( len(zsn) ) : 
        # if not decliners and decliner[isn] : continue
        # if not Bgrade and confidence[isn] != 'A' : continue

        zsnmin = zsn[isn]-zsnerr[isn]
        zsnmax = zsn[isn]+zsnerr[isn]
        for izbin in range( 1,len(zbinedges) ) : 
            zbinmin = zbinedges[izbin-1]
            zbinmax = zbinedges[izbin]
            if zsnmin>zbinmax : continue
            if zsnmax<zbinmin : continue
            if zsnmin>=zbinmin and zsnmax<zbinmax : 
                piascale = 1 # wholly inside the bin
            else : 
                zinbinmax = min( zbinmax, zsnmax )
                zinbinmin = max( zbinmin, zsnmin )
                dzsn = zsnmax - zsnmin # full range of the SN z uncertainty
                dzinbin = zinbinmax - zinbinmin # length of z dist'n in this bin
                piascale = dzinbin / dzsn # fraction of  of z dist'n in this bin
           
            Ndetz[izbin-1].append( piaz[isn] * piascale )
            dNdetzPlus[izbin-1].append( dpiazPlus[isn] * piascale )
            dNdetzMinus[izbin-1].append( dpiazMinus[isn] * piascale )

            Ndeth[izbin-1].append( piah[isn] * piascale )
            dNdethPlus[izbin-1].append( dpiahPlus[isn] * piascale )
            dNdethMinus[izbin-1].append( dpiahMinus[isn] * piascale )

    for izbin in range(len(Ndetz) ) : 
        Ndetz[izbin] = np.sum( Ndetz[izbin] )
        dNdetzPlus[izbin] = np.sum( dNdetzPlus[izbin] )
        dNdetzMinus[izbin] = np.sum( dNdetzMinus[izbin] )
        Ndeth[izbin] = np.sum( Ndeth[izbin] )
        dNdethPlus[izbin] = np.sum( dNdethPlus[izbin] )
        dNdethMinus[izbin] = np.sum( dNdethMinus[izbin] )


    return( np.array( Ndetz ), np.array( dNdetzPlus ), np.array( dNdetzMinus ), 
            np.array( Ndeth ), np.array( dNdethPlus ), np.array( dNdethMinus ) )



def rollingWindow( datfile='candels1.csv', zrange=[0.25,1.75], 
                   zstep=0.05, dzwindow=0.5, 
                   decliners=False, Bgrade=False ):
    """ count detections using a rolling window function """

    Ndet = []
    zarray = np.arange( zrange[0], zrange[1]+zstep, zstep )
    for z in zarray :
        Ndet.append( countDetections( datfile=datfile, zbinedges=[z-dzwindow/2.,z+dzwindow/2.])[0] )
    Ndet = np.array( Ndet ).ravel()

    # Nctrl Roll computed using the function snana.simulate.getNctrl()
    # TODO : allow this to be recomputed on-the-fly
    # NctrlRoll = np.array( [1.7719284614999999, 2.2324409864999999, 2.7542166554999996, 3.3049993830000002, 3.8008305465000003, 4.273929109, 4.7894081695000006, 5.2204835935000009, 5.6277937215000016, 6.0745925459999999, 6.4251440540000004, 6.6683275189999991, 6.8075586459999995, 6.9556416840000006, 7.0395220425000007, 7.0873497905000002, 7.0644104125, 7.0049088739999998, 6.9269568655000002, 6.781903733, 6.6434476295000007, 6.490404379000001, 6.355524925500001, 6.1521560600000011, 5.9337193710000014, 5.7562845095000004, 5.5132619064999995, 5.2847828369999998, 5.050937921500001, 4.8283876695000005, 4.5574059724999998] )
    NctrlRoll = np.array( [2.0114413330000001, 2.5301089665000003, 3.1134532309999998, 3.7281001389999999, 4.2792163785000001, 4.8073632490000007, 5.3857863485000008, 5.8624771670000007, 6.3147089325000012, 6.8109430420000008, 7.2041525510000008, 7.4726643535000008, 7.6224325220000004, 7.7788878994999999, 7.8611566015000003, 7.9093654324999996, 7.8644505739999984, 7.7881678259999996, 7.6836261735000004, 7.5071966010000004, 7.3464701195000011, 7.1591129120000012, 7.0016954505000015, 6.7554643545000017, 6.5031607635000013, 6.2920727160000007, 6.0071201245000001, 5.7514208345000002, 5.4916525920000012, 5.2342252325, 4.9054012819999997] )

    return( zarray, Ndet, Ndet / NctrlRoll )


def googledoc2tex( datfile='candels.csv', outdir=None ):
    """ read in an exported .csv version of the google doc table 
    write it out as a .tex table """
    from astropy.io import ascii
    from astropy.table import Table, Column
    try : 
        from astropy.coordinates import ICRS
    except : 
        from astropy.coordinates import ICRSCoordinates as ICRS
    from astropy import cosmology
    import sys
    import os
    import numpy as np

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    if not os.path.isfile( datfile ) :
        datfile = os.path.join( os.path.dirname( thisfile ), datfile )
    if not os.path.isfile( datfile ) :
        print("Cannot find %s"%datfile)
        return(-1) 

    if outdir == None : 
        outdir = os.path.abspath( os.path.expanduser("~/Dropbox/Papers/hstsnrates/TABLE/"))

    c1dat = Table.read(datfile, format='ascii' )

    # reformat the Decl. strings
    hostdec = [ dec.replace('-','$-$').replace('+','$+$') for dec in c1dat['Host DEC'] ]
    sndec = [ dec.replace('-','$-$').replace('+','$+$') for dec in c1dat['DEC'] ]
    snra = c1dat['RA']
    hostra = c1dat['Host RA']

    # make new data arrays to hold  Pz(Ia) and Phost(Ia) with uncertainties 
    piaz = c1dat['PIaz']
    dpiaz = c1dat['dPIaz']
    piaznew = [ '%.2f'%p + dp.replace('+',' $^{+').replace(' -','}_{-') + '}$' 
               for p,dp in zip( piaz, dpiaz) ]

    piaHost = c1dat['PIaHost']
    dpiaHost = c1dat['dPIaHost']
    piaHostnew = [ '%.2f'%p + dp.replace('+',' $^{+').replace(' -','}_{-') + '}$' 
               for p,dp in zip( piaHost, dpiaHost) ]

    zsn = c1dat['zSN']
    dzsn = [ '(%s)'%dz for dz in c1dat['dzSN'] ]
    zhost = c1dat['zHost']
    dzhost = [ '(%s)'%dz for dz in c1dat['dzHost'] ]

    # Compute the host - SN separation in arcsec and kpc
    cosmo = cosmology.FlatLambdaCDM( H0=_H0, Om0=_Om )
    rasn = c1dat['RA']
    decsn = c1dat['DEC']
    rahost = c1dat['Host RA']
    dechost = c1dat['Host DEC']
    zhost = c1dat['zHost']
    dzhost = c1dat['dzHost']
    darcsec, dkpc = [],[]
    for i in range(len(rasn)): 
        try : 
            csn = ICRS( colon2hms( rasn[i] ) + ' ' + colon2dms( decsn[i] ) )
            chost = ICRS( colon2hms( rahost[i] ) + ' ' + colon2dms( dechost[i] ) )
            sep = csn.separation( chost )
            try : 
                separcmin = sep.arcminute
                separcsec = sep.arcsecond
            except : 
                separcmin = sep.arcmins
                separcsec = sep.arcsecs
            dkpc.append( '%.1f'%cosmo.kpc_proper_per_arcmin( separcmin ) )
            darcsec.append( '%.2f'%separcsec )
            # print( '%9s %6.3f %6.1f'%( c1dat['name'][i], separcsec, dkpc ) )
        except : 
            darcsec.append( '\\nodata' )
            dkpc.append( '\\nodata' )
            #print( '%9s %6.3f %6.1f'%( c1dat['name'][i], -0.000, -0.0 ) )

    # icos = np.where( [ n.startswith('COS') for n in c1dat['name'] ] )[0].tolist()
    # iegs = np.where( [ n.startswith('EGS') for n in c1dat['name'] ] )[0].tolist()
    # iuds = np.where( [ n.startswith('UDS') for n in c1dat['name'] ] )[0].tolist()
    # igsa = np.where( [ n.startswith('GS') for n in c1dat['name'] ] )[0].tolist()
    # igna = np.where( [ n.startswith('GN') for n in c1dat['name'] ] )[0].tolist()
    # iwide = icos + iegs + iuds
    # ideep = igsa + igna

    preambleSN = r"""
% DELETE THE \begin{deluxetable} ABOVE THIS LINE
% IT IS REPLACED BY THE FOLLOWING:
\renewcommand{\arraystretch}{1.8}
\begin{small}
\ifms
  \begin{deluxetable}{llllllllp{2.1in}}
\else
  \begin{deluxetable*}{llllllllp{2.1in}}
\fi
\tablecolumns{9}
 """

    preambleHost = r"""
% DELETE THE \begin{deluxetable} ABOVE THIS LINE
% IT IS REPLACED BY THE FOLLOWING:
\renewcommand{\arraystretch}{1.8}
\begin{small}
\ifms
  \begin{deluxetable}{llllrcclll}
\else
  \begin{deluxetable*}{llllrcclll}
\fi
\tablecolumns{10}
 """

    footerSN = r"""
\tablenotetext{a}{ Type Ia SN classification probability from STARDUST, 
  using the redshift-dependent class prior. Uncertainties reflect systematic
  biases due to the class prior and extinction assumptions
  (Sections~\ref{sec:ClassPrior} and  \ref{sec:HostAvDistribution}). }
\tablenotetext{b}{ Type Ia SN classification probability from STARDUST, 
  using the {\em galsnid} host galaxy prior. Uncertainties reflect systematic
  biases due to the class prior and extinction assumptions.}
\tablenotetext{c}{ Posterior redshift and uncertainty, as determined 
  by the STARDUST light curve fit.}
\tablenotetext{d}{ The {\it host / SN} values indicate whether the redshift is derived from the host galaxy, the SN itself, or a combination;  {\it spec-z / phot-z} specify a spectroscopic or photometric redshift. A value of {\it host+SN phot-z} means the redshift is derived from a STARDUST light curve fit, with the host galaxy phot-z used as a prior.}

\ifms  \end{deluxetable}
\else  \end{deluxetable*}
\fi
\end{small}
\renewcommand{\arraystretch}{1.0}

% DELETE THE \end{deluxetable} BELOW THIS LINE
% IT IS REPLACED BY THE PRECEDING

"""

    footerHost = r"""
\tablenotetext{a}{ Physical separation between the SN and center of the host, 
  computed from the measured angular separation in the preceding column, 
  assuming a flat \LCDM\ cosmology with \Ho=70, \Om=0.3}
\tablenotetext{b}{ Visual classifications for host galaxy morphology: 
   s = spheroid, d = disk, i = irregular}
\tablenotetext{c}{ Template-matching classification of host galaxy SED:
   P = Passive, A = Active, SB = Starburst type}
\tablenotetext{d}{ Unpublished spectroscopic observations are given as {\em
    Observatory}+{\em Instrument} ({\em name of PI}). Host galaxy
  photometric redshifts are marked as {\em phot-z} (Dahlen et al. in
  prep).}

\ifms  \end{deluxetable}
\else  \end{deluxetable*}
\fi
\end{small}
\renewcommand{\arraystretch}{1.0}

% DELETE THE \end{deluxetable} BELOW THIS LINE
% IT IS REPLACED BY THE PRECEDING

"""

    latexdictSN={'tabletype':'deluxetable*','preamble': preambleSN,'tablefoot': footerSN } # 'tablehead':header,
    latexdictHost={'tabletype':'deluxetable*','preamble': preambleHost,'tablefoot': footerHost } # 'tablehead':header,

    # Output latex tables are segregated by redshift
    ihiz = np.where( c1dat['zSN'] > 1.5 )
    iloz = np.where( c1dat['zSN'] <= 1.5 )


    # Construct a new table for the SN data and write it out to latex
    caphiz = r'16 Supernovae with $z>1.5$ (see Appendix for the remainder)\label{tab:highzSN}'
    caploz = r'49 Supernovae with $z<1.5$\label{tab:lowzSN}'
    snData = [ c1dat['name'], snra, sndec, piaznew, piaHostnew, zsn, dzsn, c1dat['zSN Source'] ]
    snNames = [ 'Name','R.A. (J2000)','Decl. (J2000)','P(Ia$|$D$_{z}$)\\tablenotemark{a}',
                  'P(Ia$|$D$_{host}$)\\tablenotemark{b}','z$_{\\mbox{\\scriptsize SN}}$\\tablenotemark{c}','($\pm$)',
                  '$z$ Source\\tablenotemark{d}' ] 
    snTable = Table( data=snData, names=snNames  )

    hiztexfile = os.path.join( outdir, 'highzSNTable.tex')
    loztexfile = os.path.join( outdir, 'lowzSNTable.tex')
    ascii.write( snTable[ihiz], hiztexfile, Writer=ascii.AASTex,
                 caption=caphiz, latexdict=latexdictSN, col_align='lllllllll' )
    ascii.write( snTable[iloz], loztexfile, Writer=ascii.AASTex,
                 caption=caploz, latexdict=latexdictSN, col_align='lllllllll' )
    for texfile in [hiztexfile, loztexfile] : 
        fin = open( texfile, 'r' )
        inlines = fin.readlines()
        fin.close()
        fout = open( texfile, 'w' )
        fout.writelines( inlines[4:-5] )
        fout.close()


    # Construct a new table for the Host data and write it out to latex
    caphiz = r'Host galaxies of 16 Supernovae with $z>1.5$ (see Appendix \ref{app:lowz} for the remainder)\label{tab:highzHosts}'
    caploz = r'Host galaxies of 49 Supernovae with $z<1.5$\label{tab:lowzHosts}'
    hostData = [ c1dat['name'], hostra, hostdec, darcsec, dkpc, c1dat['Host Morphology'], c1dat['Host SED Type'], zhost, dzhost, c1dat['zHost Source'] ]
    hostNames = [ 'SN','R.A. (J2000)','Decl. (J2000)', 'd[\\arcsec]','d[kpc]\\tablenotemark{a}',
                  'Morph.\\tablenotemark{b}','SED\\tablenotemark{c}',
                  'z$_{\\mbox{host}}$','($\pm$)','$z$ Reference\\tablenotemark{d}' ] 
    hostTable = Table( data=hostData, names=hostNames  )

    hiztexfile = os.path.join( outdir, 'highzHostTable.tex')
    loztexfile = os.path.join( outdir, 'lowzHostTable.tex')
    ascii.write( hostTable[ihiz], hiztexfile, Writer=ascii.AASTex,
                 caption=caphiz, latexdict=latexdictHost, col_align='llllrcclll' )
    ascii.write( hostTable[iloz], loztexfile, Writer=ascii.AASTex,
                 caption=caploz, latexdict=latexdictHost, col_align='llllrcclll' )
    for texfile in [hiztexfile, loztexfile] : 
        fin = open( texfile, 'r' )
        inlines = fin.readlines()
        fin.close()
        fout = open( texfile, 'w' )
        fout.writelines( inlines[4:-5] )
        fout.close()
    
    return(snTable, hostTable)

    

def poissonLimits( N, confidence=1 ):
    """   
    Adapted from P.K.G.Williams : 
    http://newton.cx/~peter/2012/06/poisson-distribution-confidence-intervals/

    Let's say you observe n events in a period and want to compute the k
    confidence interval on the true rate - that is, 0 < k <= 1, and k =
    0.95 would be the equivalent of 2sigma. Let a = 1 - k, i.e. 0.05. The
    lower bound of the confidence interval, expressed as a potential
    number of events, is
       scipy.special.gammaincinv (n, 0.5 * a)
    and the upper bound is
       scipy.special.gammaincinv (n + 1, 1 - 0.5 * a)
    
    The halving of a is just because the 95% confidence interval is made
    up of two tails of 2.5% each, so the gammaincinv function is really,
    once you chop through the obscurity, exactly what you want.

    INPUTS : 
      N : the number of observed events

      confidence : may either be a float <1, giving the exact 
          confidence limit desired (e.g.  0.95 or 0.99)
          or it can be an integer in [1,2,3], in which case 
          we set the desired confidence interval to match 
          the 1-, 2- or 3-sigma gaussian confidence limits
             confidence=1 gives the 1-sigma (68.3%) confidence limits
             confidence=2  ==>  95.44% 
             confidence=3  ==>  99.74% 
    """
    from scipy.special import gammaincinv as ginv
    if confidence<1 : k = confidence
    elif confidence==1 : k = 0.6826
    elif confidence==2 : k = 0.9544
    elif confidence==3 : k = 0.9974
    else :
        print( "ERROR : you must choose nsigma from [1,2,3]")
        return( None )
    lower = ginv( N, 0.5 * (1-k) )
    upper = ginv( N+1, 1-0.5*(1-k) )
    return( lower, upper )
 

        

def colon2hms( coordstr ) :
    """ convert an R.A. string containing colons
    into one using 'h' 'm' 's' """
    coordlist = coordstr.split(':')
    return( '%sh%sm%ss'%tuple(coordlist) )

def colon2dms( coordstr ) :
    """ convert a DEC string containing colons
    into one using 'd' 'm' 's' """
    coordlist = coordstr.split(':')
    return( '%sd%sm%ss'%tuple(coordlist) )


def mkRatePointsFig( rgb=True, outfile=None, ageaxis=True, logscale=False, presfig=False, **kwarg):
    """ 
    Make a plot of measured SN rate points vs z  
    for the Rodney et al 2013 rates paper.
    """
    from cosmo import agez, zfromt
    from hstsnpipe.tools.rates.ratetable import RATELISTALL, RATELISTLOWZ, RATELISTHIGHZ, R13,R13h
    from hstsnpipe.tools.figs import plotsetup
    from matplotlib import pyplot as pl

    if presfig : 
        plotsetup.presfig()
    else :
        plotsetup.fullpaperfig(figsize=[8,5])

    if logscale : scale = 1e-4
    else : scale = 1

    ax1 = pl.axes( [0.1,0.12,0.87,0.76] )
    for i in range(len(RATELISTLOWZ)):
        R = RATELISTLOWZ[i]
        R.mfc='0.8'
        R.mec='0.8'
        R.color='0.8'
        R.mew=0
        R.marker='o'
        R.zerrplus = np.zeros( len( R.zerrplus ) )
        R.zerrminus = np.zeros( len( R.zerrplus ) )
        R.plot( thicksys=False, zorder=-10*i, scalerates=scale, **kwarg )
    for i in range(len(RATELISTHIGHZ)-1):
        R = RATELISTHIGHZ[i]
        R.ms=12
        R.zerrplus = np.zeros( len( R.zerrplus ) )
        R.zerrminus = np.zeros( len( R.zerrplus ) )
        pl.plot( R.z, R.rate, marker=R.marker, ms=12, mew=0, ls=' ',
                color='w', alpha=1, zorder=10*i )
        R.plot( thicksys=False, zorder=10*(i+1), scalerates=scale, alpha=0.3, **kwarg )

    R13.mfc='darkorange'
    R13.mec='darkred'
    R13.color='darkred'
    R13.ms=15
    R13.plot( thicksys=True, zorder=1000, scalerates=scale )

    #from matplotlib.font_manager import FontProperties
    #fs=FontProperties(size='9')

    ax = pl.gca()
    ax.set_xlim([-0.05, 2.81])
    pl.xlabel("Redshift")

    if logscale : 
        ax.set_yscale("log",nonposy="clip")
        ax.set_ylim([0.011*scale, 5.55*scale])

        pl.ylabel(r"SNR ~$\left[{\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")
        ax.yaxis.set_ticklabels( ['','$10^{-5}$', '$10^{-4}$', ''] )

        #ax1.text( 2.29, 0.32*scale, 'CANDELS', ha='left', va='top', color='darkred', rotation=-90 )
        ax1.text( 2.16, 0.32*scale, 'CANDELS', ha='left', va='top', color='darkred', rotation=-90 )
        ax1.text( 1.64, 0.28*scale, 'CLASH', ha='right', va='top', color='darkmagenta', rotation=-90 )
        ax1.text( 1.61, 1.02*scale, 'GOODS', ha='center', va='bottom', color='darkgreen', rotation=-90 )
        ax1.text( 1.54, 2.53*scale, 'CSS', ha='center', va='bottom', color='darkcyan', rotation=-90)
        ax1.text( 1.70, 1.75*scale, 'SDF', ha='center', va='bottom', color='0.2', rotation=-90 )

        ax1.text( 2.45, 3.0*scale, r'\noindent Stat.+Syst.\\ Error',ha='left',va='center', color='darkred')
        ax1.text( 2.45, 1.7*scale, 'Syst. Only',ha='left',va='center', color='darkorange')
        ax1.plot( [2.27, 2.42], scale*np.array([ 1.37, 2.80]), marker=' ', ls='-', color='darkred', lw=0.8 )
        ax1.plot( [2.28, 2.42], scale*np.array([ 0.91, 1.71]), marker=' ', ls='-', color='darkorange', lw=0.8 )

    else : 
        ax.set_ylim([0.0, 2.55])
        pl.ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")

        if presfig: txtrot,txtshift=90,0.3
        else : txtrot,txtshift=0,0
        ax1.text( 2.33, 0.72, 'CANDELS', ha='left', va='top', color='darkred' )
        ax1.text( 1.64, 0.21, 'CLASH', ha='left', va='top', color='darkmagenta' )
        ax1.text( 1.27, 1.31+2*txtshift, 'GOODS', ha='left', va='top', color='darkgreen',rotation=txtrot )
        ax1.text( 1.13, 1.44+txtshift, 'CSS', ha='right', va='top', color='darkcyan',rotation=txtrot )
        ax1.text( 1.67, 1.15+txtshift, 'SDF', ha='right', va='top', color='0.2',rotation=txtrot )
        
        ax1.text( 1.9, 2.22, 'Stat. + Syst. Error',ha='left',va='center', color='darkred')
        ax1.text( 1.95, 2.0, 'Systematic Only',ha='left',va='center', color='darkorange')
        ax1.plot( [1.77, 1.92], [ 1.45, 2.15], marker=' ', ls='-', color='darkred', lw=0.8 )
        ax1.plot( [1.8, 2.02], [ 1.02, 1.93], marker=' ', ls='-', color='darkorange', lw=0.8 )

    axtop = ax1.twiny()
    axtop.set_xlim( ax1.get_xlim() )
    if ageaxis: 
        ageticks = np.array( [13,8,5,3] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')


    pl.draw()
    if outfile: 
        pl.savefig(outfile)
        return(outfile)
    else:
        return( ax )



def mkPreCANDELSRatesFig( outfile=None, ageaxis=True,  **kwarg):
    """ 
    Make a plot of measured pre-CANDELS SN rate points vs z  
    for the Rodney et al 2014 rates paper.
    """
    from cosmo import agez, zfromt
    from hstsnpipe.tools.rates.ratetable import RATELISTPRECC
    from hstsnpipe.tools.figs import plotsetup
    from matplotlib import pyplot as pl
    plotsetup.halfpaperfig()

    ax1 = pl.axes( [0.15,0.12,0.8,0.76] )
    for R in RATELISTPRECC:
        if R.ref.startswith("Barbary") : R.z[-1] = 1.515
        if R.ref.startswith("Dahlen") : R.z[-1] = 1.64
        R.zerrplus = np.zeros( len( R.zerrplus ) )
        R.zerrminus = np.zeros( len( R.zerrplus ) )
        R.plot( thicksys=False, **kwarg )

    ax = pl.gca()
    ax.set_xlim([-0.05, 1.81])
    pl.xlabel("Redshift")

    ax.set_ylim([0.0, 2.55])
    pl.ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")

    ax.text( 1.62, 0.53, 'GOODS', ha='right', va='bottom', color='darkgreen', rotation=90 )
    ax.text( 1.50, 0.98, 'CSS', ha='right', va='bottom', color='darkcyan', rotation=90 )
    ax.text( 1.67, 1.15, 'SDF', ha='right', va='bottom', color='0.2', rotation=90 )

    axtop = ax1.twiny()
    axtop.set_xlim( ax1.get_xlim() )
    if ageaxis: 
        ageticks = np.array( [13,8,6,4] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')

    pl.draw()
    if outfile: 
        pl.savefig(outfile)
        return(outfile)
    else:
        return( ax )




def mkGalsnidTestFig():
    """ make a figure showing
    the results of the galsnid systematics check,
    for the appendix of the CANDELS rates paper """
    from cosmo import agez, zfromt
    from hstsnpipe.tools.rates.ratetable import RATELISTALL, RATELISTLOWZ, RATELISTHIGHZ, R13,R13h
    from hstsnpipe.tools.figs import plotsetup
    from matplotlib import pyplot as pl
    plotsetup.fullpaperfig(figsize=[8,5])


    # lower axes, showing the systematics check
    # from adopting the host-based prior
    ax  = pl.axes( [0.12,0.1,0.84,0.2] )

    # First, plot just the points, no error bars 
    # showing the % change relative to the baseline rate
    R13h.rate = (R13h.rate-R13.rate)/R13.rate * 100
    R13h.errsysplus = R13.errsysplus * 0
    R13h.errsysminus = R13.errsysminus * 0
    R13h.errstatplus = R13.errstatplus * 0
    R13h.errstatminus = R13.errstatminus * 0
    R13h.zerrplus = R13.zerrplus * 0
    R13h.zerrminus = R13.zerrminus * 0
    R13h.mfc='None'
    R13h.ms=10
    R13h.alpha=0.0
    R13h.plot( thicksys=False )

    # Now plot just the true R13 systematic errors
    # centered on the 0 line, showing the 
    # systematic error range from the baseline analysis
    R13h.rate = R13h.rate * 0
    R13h.mfc='darkorange'
    R13h.alpha=1.0
    R13h.errsysplus = R13.errsysplus/R13.rate  * 100
    R13h.errsysminus = R13.errsysminus/R13.rate  * 100
    R13h.marker= '-'
    R13h.ms=8
    R13h.plotSysErr()

    ax2.set_xlim( ax1.get_xlim() )
    ax2.set_ylim( -105, 125 )
    #ax2.xaxis.set_ticks_position('bottom')
    ax2.axhline( 0, ls='--', color='0.5', lw=0.7, zorder=-1000 )
    ax2.text( 0.05, 0.85, 'Systematics Check: adopting the host-based prior',
              transform=ax2.transAxes, ha='left', va='top', fontsize='small' )
    ax2.set_ylabel( '\%\ Change' )
    ax2.set_xlabel( 'Redshift' )
    ax2.yaxis.set_ticks( [-50, 0, 50 ] )

    #ax2.text( 1.9, 2.22, 'Baseline\n Syst. Err.',ha='left',va='center', color='darkred')
    #ax2.text( 1.95, 2.0, 'Systematic Only',ha='left',va='center', color='darkorange')
    #ax1.plot( [1.77, 1.92], [ 1.45, 2.15], marker=' ', ls='-', color='darkred', lw=0.8 )
    #ax1.plot( [1.8, 2.02], [ 1.02, 1.93], marker=' ', ls='-', color='darkorange', lw=0.8 )
