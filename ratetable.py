"""
090514  S.Rodney
Table of TNSN Rates compiled from the literature (see esp. Blanc:2008, Valiante:2009)
Redshift bin numbers are zmin, zeff, zmax 
where zeff is the volume-weighted effective redshift.
Rates are given in units of 10-4 yr-1 Mpc-3 h_70-3

# z-bin            rate       error(stat)    error(syst)           Reference
 0.2,0.4,0.6       0.53      +0.39,-0.17     +0.00,-0.00          Kuznetsova:2008
 0.6,0.8,1.0       0.93      +0.25,-0.25     +0.00,-0.00          Kuznetsova:2008
 1.0,1.2,1.4       0.75      +0.35,-0.30     +0.00,-0.00          Kuznetsova:2008
 1.4,1.55,1.7      0.12      +0.58,-0.12     +0.00,-0.00          Kuznetsova:2008

 0.2,0.47,0.6      0.80      +0.37,-0.27     +1.66,-0.26          Dahlen:2008
 0.6,0.83,1.0      1.30      +0.33,-0.27     +0.73,-0.51          Dahlen:2008
 1.0,1.21,1.4      1.32      +0.36,-0.29     +0.38,-0.32          Dahlen:2008
 1.4,1.61,1.7      0.42      +0.39,-0.23     +0.19,-0.14          Dahlen:2008

 0.2,0.4,0.6       0.69      +0.20,-0.40     +0.00,-0.00          Dahlen:2004
 0.4,0.8,1.0       1.57      +0.50,-0.40     +0.00,-0.00          Dahlen:2004
 1.0,1.2,1.4       1.15      +0.70,-0.50     +0.00,-0.00          Dahlen:2004
 1.4,1.6,1.8       0.44      +0.50,-0.3      +0.00,-0.00          Dahlen:2004

 0.5,0.80,1.0      0.43      +0.36,-0.32     +0.00,-0.00          Poznanski:2007
 1.0,1.30,1.5      1.05      +0.45,-0.56     +0.00,-0.00          Poznanski:2007
 1.5,1.75,2.0      0.81      +0.79,-0.60     +0.00,-0.00          Poznanski:2007

 0.2,0.25,0.3      0.17      +0.17,-0.17     +0.00,-0.00          Barris:2006
 0.3,0.35,0.4      0.53      +0.24,-0.24     +0.00,-0.00          Barris:2006
 0.4,0.45,0.5      0.73      +0.24,-0.24     +0.00,-0.00          Barris:2006
 0.5,0.55,0.6      2.04      +0.38,-0.38     +0.00,-0.00          Barris:2006
 0.6,0.65,0.7      1.49      +0.31,-0.31     +0.00,-0.00          Barris:2006
 0.7,0.75,0.8      1.78      +0.34,-0.34     +0.00,-0.00          Barris:2006
 
 0.3,0.46,0.5      0.48      +0.17,-0.17     +0.00,-0.00          Tonry:2003
 0.2,0.47,0.6      0.42      +0.06,-0.06     +0.13,-0.09          Neill:2006      

 0.0,0.13,0.2      0.20      +0.07,-0.05     +0.05,-0.05          Blanc:2004
# 0.0,0.14,0.2      0.34      +0.27,-0.16     +0.11,-0.06          Hardin:2000 (superceded by blanc)

 0.0,0.01,0.05     0.28      +0.09,-0.09     +0.00,-0.00          Cappellaro:1999 
 0.25,0.55,0.85    0.52      +0.10,-0.09     +0.00,-0.00          Pain:2002
 0.02,0.03,0.04    0.28      +0.11,-0.11     +0.00,-0.00          Mannucci:2005
 0.05,0.10,0.15    0.32      +0.15,-0.15     +0.00,-0.00          Madgwick:2003
 0.06,0.11,0.16    0.37      +0.10,-0.10     +0.00,-0.00          Strolger:2003
 0.15,0.20,0.25    0.20      +0.08,-0.08     +0.00,-0.00          Horesh:2008
 0.1,0.30,0.6      0.34      +0.16,-0.15     +0.21,-0.22          Botticella:2008

"""
from __init__ import rate
from numpy import array, arange, zeros, ones, digitize, average, sqrt
from scipy.interpolate import interp1d
from copy import deepcopy

def rate_prediction( fp=0.5 ): 
    zvals = arange( 0.0, 3.5, 0.1 )
    if fp == 0.5 : 
        snrt1vals = array([ 0.        ,  0.25715084,  0.29598062,  0.33910441,  0.38636852,
                            0.437304  ,  0.49103477,  0.54620207,  0.60093896,  0.6529356 ,
                            0.69962025,  0.7384616 ,  0.76733809,  0.78487899,  0.79067015,
                            0.78526867,  0.77000951,  0.74671747,  0.71738788,  0.68393567,
                            0.64803069,  0.61103073,  0.57397242,  0.53760484,  0.50244008,
                            0.46880537,  0.43688996,  0.40678193,  0.37850176,  0.35202232,
                            0.32728607,  0.30421701,  0.28272958,  0.26273183,  0.24413385])
    snrate_tminus1 = interp1d( _zvals, _snrt1vals )




def poissonerrExact( N, confidence=1 ):
    """   
    The poisson error for a signal with N counts.
    N may be a scalar or a numpy array.
    This returns the exact solution, accurate even for low N or N=0.

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
    from numpy import mean
    if confidence<1 : k = confidence
    elif confidence==1 : k = 0.6826
    elif confidence==2 : k = 0.9544
    elif confidence==3 : k = 0.9974
    else :
        print( "ERROR : you must choose nsigma from [1,2,3]")
        return( None )

    lower = ginv( N, 0.5 * (1-k) )
    upper = ginv( N+1, 1-0.5*(1-k) )
    meanerr = (upper-lower) / 2.
    return( upper-N, N-lower  )


# define a rectangle shape for scatter points
rectangle1 = [ (-.1,.5),(.1,.5),(.1,-.5),(-.1,-.5) ]
rectangle2 = [ (-.2,.5),(.2,.5),(.2,-.5),(-.2,-.5) ]


BT06=rate()
BT06.z           = array([0.27,0.37,0.47,0.57,0.67,0.77])
BT06.zerrplus    = array([0.05,0.05,0.05,0.05,0.05,0.05]) 
BT06.zerrminus   = array([0.05,0.05,0.05,0.05,0.05,0.05])
BT06.rate        =array([0.17,0.53,0.73,2.04,1.49,1.78])
BT06.errstatplus =array([0.17,0.24,0.24,0.38,0.31,0.34])
BT06.errstatminus=array([0.169,0.24,0.24,0.38,0.31,0.34])
BT06.errsysplus  =array([0.00,0.00,0.00,0.00,0.00,0.00])
BT06.errsysminus =array([0.00,0.00,0.00,0.00,0.00,0.00])
BT06.reference="Barris and Tonry 2006"
BT06.ref="Barris \& Tonry 2006 (IfA)"
BT06.marker='o'
BT06.color='green'
BT06.mec='darkgreen'
BT06.mfc='lightgreen'
BT06.ms=10

B04 = rate()
B04.ndet        =array([14])
B04.z           =array([0.13])
B04.zerrplus    =array([0.07])
B04.zerrminus   =array([0.13])
B04.rate        =array([0.2])
B04.errstatplus =array([0.07])
B04.errstatminus=array([0.05])
B04.errsysplus  =array([0.05 ])
B04.errsysminus =array([0.05 ])
B04.reference = "Blanc et al. 2004"
B04.ref = "Blanc+ 2004"
B04.marker='o'
B04.color='gray'
B04.mec='gray'
B04.mfc='white'
B04.ms=8

B08 = rate()
B08.ndet        =array([31.05])
B08.z           =array([0.30])
B08.zerrplus    =array([0.3])
B08.zerrminus   =array([0.2])
B08.rate        =array([0.34])
B08.errstatplus =array([0.16])
B08.errstatminus=array([0.15])
B08.errsysplus  =array([0.21 ])
B08.errsysminus =array([0.22 ])
B08.reference="Botticella et al. 2008"
B08.ref="Botticella+ 2008"
B08.marker='o'
B08.color='gray'
B08.mec='gray'
B08.mfc='white'
B08.ms=8

B12=rate()
B12.ndet = array([ 1.18, 5.63, 1.12 ])
B12.z         = array([0.807,1.187,1.535])
B12.zerrplus  = array([0.193,0.213,0.265])
B12.zerrminus = array([0.207,0.187,0.135])
B12.rate      = array([1.18, 1.33, 0.77])
B12.errstatplus  = array([0.60, 0.65, 1.07])
B12.errstatminus = array([0.45, 0.49, 0.54])
B12.errsysplus   = array([0.44, 0.69, 0.44])
B12.errsysminus  = array([0.28, 0.26, 0.76])
# If you include their non-detection at low-z.
#B12.z         = array([0.442,0.807,1.187,1.535])
#B12.zerrplus  = array([0.158,0.193,0.213,0.265])
#B12.zerrminus = array([0.242,0.207,0.187,0.135])
#B12.rate      = array([0.00, 1.18, 1.33, 0.77])
#B12.errstatplus  = array([0.50, 0.60, 0.65, 1.07])
#B12.errstatminus = array([0.00, 0.45, 0.49, 0.54])
#B12.errsysplus   = array([0.00, 0.44, 0.69, 0.44])
#B12.errsysminus  = array([0.00, 0.28, 0.26, 0.77])
B12.reference="Barbary et al. 2012"
B12.ref="Barbary+ 2012 [HST+ACS]"
B12.marker='v'
B12.color='c'
B12.mec='darkblue'
B12.mfc='c'
B12.ms=8


C99 = rate()
C99.ndet        =array([70])
C99.z           =array([0.01])
C99.zerrplus    =array([0.01])
C99.zerrminus   =array([0.01])
C99.rate        =array([0.28])
C99.errstatplus =array([0.09])
C99.errstatminus=array([0.09])
C99.errsysplus  =array([0.00 ])
C99.errsysminus =array([0.00 ])
C99.reference = "Cappellaro et al. 1999"
C99.ref = "Cappellaro+ 1999"
C99.marker='o'
C99.color='gray'
C99.mec='gray'
C99.mfc='white'
C99.ms=8


D04=rate()
D04.z=array([0.4,0.8,1.2,1.6])
D04.zerrplus = array([0.2, 0.2, 0.2, 0.2])
D04.zerrminus = array([0.2, 0.2, 0.2, 0.2])
D04.rate        =array([0.69,1.57,1.15,0.44])
D04.errstatplus =array([0.2, 0.5, 0.7, 0.5])
D04.errstatminus=array([0.4, 0.4, 0.5, 0.3])
D04.errsysplus=array([0,0,0,0])
D04.errsysminus=array([0,0,0,0])
D04.reference="Dahlen et al. 2004"
D04.ref = "Dahlen+ 2004"
D04.marker='o'
D04.color='gray'
D04.mec='gray'
D04.mfc='gray'
D04.ms=8

D08=rate()
D08.z        =array([0.47,0.83,1.21,1.61])
D08.zerrplus = array([0.2, 0.2, 0.2, 0.2])
D08.zerrminus = array([0.2, 0.2, 0.2, 0.2])
D08.rate        =array([0.80,1.30,1.32,0.42])
D08.errstatplus =array([0.37,0.33,0.36,0.39])
D08.errstatminus=array([0.27,0.27,0.29,0.23])
D08.errsysplus  =array([1.66,0.73,0.38,0.19])
D08.errsysminus =array([0.26,0.51,0.32,0.14])
D08.reference="Dahlen et al. 2008"
D08.ref="Dahlen+ 2008 [HST+ACS]"
D08.marker='^'
D08.color='darkgreen'
D08.mec='darkgreen'
D08.mfc='lightgreen'
D08.ms=8
D08.ndet = array([ 3, 8, 20, 3 ])

D10=rate()
D10.z        =array([0.0375,   0.1,  0.15,  0.2, 0.25, 0.3])
D10.ndet     =array([     5,    10,    64,  111,  160, 166])
D10.zerrplus = array([0.0125, 0.025, 0.025, 0.025, 0.025, 0.025])
D10.zerrminus = array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
D10.rate=array([0.278, 0.259, 0.307, 0.348, 0.365, 0.434])
D10.errstatplus=array([0.112,0.052,0.038,0.032,0.031,0.037])
D10.errstatminus=array([0.083,0.044,0.034,0.030,0.028,0.034])
D10.errsysplus=array([0.015,0.028,0.035,0.082,0.182,0.396])
D10.errsysminus=array([0.00,0.001,0.005,0.007,0.012,0.016])
D10.reference="Dilday et al. 2010"
D10.ref="Dilday+ 2010 [SDSS]"
D10.marker='o'
D10.color='gray'
D10.mec='gray'
D10.mfc='white'
D10.ms=8


G11 = rate()
G11.ndet        =array([20.3,28.0,10.0])
G11.z           =array([0.73,1.23,1.69])
G11.zerrplus    =array([0.27,0.23,0.23])
G11.zerrminus   =array([0.23,0.24,0.23])
G11.rate        =array([0.79,0.84,1.02])
G11.count       =array([20,28,10])
G11.errtotplus =array([0.33,0.25,0.54])
G11.errtotminus=array([0.41,0.28,0.37])
G11.errstatplus = G11.rate * poissonerrExact( G11.ndet )[0] / G11.ndet
G11.errstatminus= G11.rate * poissonerrExact( G11.ndet )[1] / G11.ndet
G11.errsysplus  = sqrt( G11.errtotplus**2 - G11.errstatplus**2 )
G11.errsysminus = sqrt( G11.errtotminus**2 - G11.errstatminus**2 )
G11.reference="Graur et al. 2011"
G11.ref="Graur+ 2011 [SDF]"
G11.marker='o'
G11.color='k'
G11.mec='k'
G11.mfc='0.3'
G11.ms=8


GM13 = rate()
GM13.ndet        =array([90])
GM13.z           =array([0.11])
GM13.zerrplus    =array([0.10])
GM13.zerrminus   =array([0.10])
GM13.rate        =array([0.247])
GM13.count       =array([90])
GM13.errstatplus =array([0.029]) 
GM13.errstatminus=array([0.026]) 
GM13.errsysplus  =array([0.016]) 
GM13.errsysminus =array([0.031]) 
GM13.reference="Graur & Maoz 2013"
GM13.ref="Graur+Maoz 2013"
GM13.marker='o'
GM13.color='k'
GM13.mec='k'
GM13.mfc='0.3'
GM13.ms=8


H00 = rate()
H00.z           =array([0.14])
H00.ndet        =array([4])
H00.zerrplus    =array([0.06])
H00.zerrminus   =array([0.12])
H00.rate        =array([0.28])
H00.errstatplus =array([0.22])
H00.errstatminus=array([0.13])
H00.errsysplus  =array([0.07 ])
H00.errsysminus =array([0.04 ])
H00.reference="Hardin et al. 2000"
H00.ref="Hardin+ 2000"
H00.marker='o'
H00.color='gray'
H00.mec='gray'
H00.mfc='white'
H00.ms=8

H08 = rate()
H08.ndet        =array([17])
H08.z           =array([0.2])
H08.zerrplus    =array([0.05])
H08.zerrminus   =array([0.05])
H08.rate        =array([0.189])
H08.errstatplus =array([0.042])
H08.errstatminus=array([0.042])
H08.errsysplus  =array([0.00 ])
H08.errsysminus =array([0.00 ])
H08.reference="Horesh et al. 2008"
H08.ref="Horesh+ 2008"
H08.marker='o'
H08.color='gray'
H08.mec='gray'
H08.mfc='white'
H08.ms=8


K08 = rate()
K08.z = array([0.4, 0.8, 1.2, 1.55])
K08.zerrplus = array([0.2, 0.2, 0.2, 0.15])
K08.zerrminus = array([0.2, 0.2, 0.2, 0.15])
K08.rate = array([0.53,0.93,0.75,0.12])
K08.errstatplus = array([0.39,0.25,0.35,0.58])
K08.errstatminus= array([0.17,0.25,0.30,0.11])
K08.errsysplus = array([0,0,0,0])
K08.errsysminus= array([0,0,0,0])
K08.reference="Kuznetsova et al. 2008"
K08.ref="Kuznetsova+ 2008 (HST)"
K08.marker='o'
K08.color='gray'
K08.mec='gray'
K08.mfc='white'
K08.ms=8

L11 = rate()
L11.ndet        =array([274])
L11.z           =array([0.01])
L11.zerrplus    =array([0.009])
L11.zerrminus   =array([0.009])
L11.rate        =array([0.265])
L11.errstatplus =array([0.034])
L11.errstatminus=array([0.033])
L11.errsysplus  =array([0.043])
L11.errsysminus =array([0.043])
L11.reference="Li et al. 2011"
L11.ref="Li+ 2011"
L11.marker='o'
L11.color='gray'
L11.mec='gray'
L11.mfc='white'
L11.ms=8


M03 = rate()
M03.z           =array([0.1])
M03.zerrplus    =array([0.03])
M03.zerrminus   =array([0.03])
M03.rate        =array([0.32])
M03.errstatplus =array([0.15])
M03.errstatminus=array([0.15])
M03.errsysplus  =array([0.00 ])
M03.errsysminus =array([0.00 ])
M03.reference="Madgwick et al. 2003"
M03.ref="Madgwick+ 2003"
M03.marker='o'
M03.color='gray'
M03.mec='gray'
M03.mfc='white'
M03.ms=8


M12 = rate()
M12.ndet        = array([7])
M12.z           = array([0.62 ])
M12.zerrplus    = array([0.38 ])
M12.zerrminus   = array([0.52 ])
M12.rate        = array([1.29 ])
M12.errstatplus = array([0.88 ])
M12.errsysplus  = array([0.27 ])
M12.errstatminus= array([0.57 ])
M12.errsysminus = array([0.28 ])
M12.reference="Melinder et al. 2012"
M12.ref="Melinder+ 2012"
M12.marker='o'
M12.color='gray'
M12.mec='gray'
M12.mfc='white'
M12.ms=8

N06 = rate()
N06.z           =array([0.47])
N06.zerrplus    =array([0.13])
N06.zerrminus   =array([0.27])
N06.rate        =array([0.42])      
N06.errstatplus =array([0.06])
N06.errstatminus=array([0.06])
N06.errsysplus  =array([0.13])
N06.errsysminus =array([0.09])
N06.reference='Neill et al. 2006'
N06.ref='Neill+ 2006'
N06.marker='o'
N06.color='gray'
N06.mec='gray'
N06.mfc='w'
N06.ms=8

P02 = rate()
P02.ndet        =array([38])
P02.z           =array([0.55])
P02.zerrplus    =array([0.3])
P02.zerrminus   =array([0.3])
P02.rate        =array([0.568])
P02.errstatplus =array([0.10])
P02.errstatminus=array([0.09])
P02.errsysplus  =array([0.10 ])
P02.errsysminus =array([0.09 ])
P02.reference = "Pain et al. 2002"
P02.ref = "Pain+ 2002"
P02.marker='o'
P02.color='gray'
P02.mec='gray'
P02.mfc='w'
P02.ms=8

P07=rate()
P07.z=array([0.80,1.30,1.75])
P07.zerrplus = array([0.2, 0.2, 0.25])
P07.zerrminus = array([0.3, 0.3, 0.25])
P07.rate=array([0.43,1.05,0.81])
P07.errstatplus=array([0.36,0.45,0.79])
P07.errstatminus=array([0.32,0.56,0.60])
P07.errsysplus=array([0.00,0.00,0.00])
P07.errsysminus=array([0.00,0.00,0.00])
P07.reference="Poznanski et al. 2007"
P07.ref="Poznanski+ 2007"
P07.marker='o'
P07.color='gray'
P07.mec='gray'
P07.mfc='w'
P07.ms=8

P12=rate()
P12.z=arange(0.15,1.10,0.1)
P12.zerrplus=zeros(len(P12.z))+0.05
P12.zerrminus=zeros(len(P12.z))+0.05
P12.ndet = array([ 4, 16, 31, 42, 72, 91, 110, 128, 141, 50 ])
P12.rate        = array([0.14,0.28,0.36,0.36,0.48,0.48,0.58,0.57,0.77,0.74])
# correct the rates for missing 91bg-like SNe
P12.rate *= 100/85.
P12.errstatplus = array([0.09,0.07,0.06,0.06,0.06,0.05,0.06,0.05,0.08,0.12])
P12.errstatminus= array([0.09,0.07,0.06,0.06,0.06,0.05,0.06,0.05,0.08,0.12])
P12.errsysplus  = array([0.06,0.06,0.05,0.04,0.04,0.04,0.05,0.06,0.10,0.10])
P12.errsysminus = array([0.12,0.07,0.06,0.05,0.05,0.06,0.07,0.07,0.12,0.13])
P12.reference="Perrett et al. 2012"
P12.ref="Perrett+ 2012 [SNLS]"
P12.marker='o'
P12.color='gray'
P12.mec='gray'
P12.mfc='w'
P12.ms=8


R14=rate()  # Full CANDELS sample, No Decliners
R14.ndet        = array([1.46,7.19,8.47,5.54,1.24])
R14.rate        = array([0.36,0.51,0.64,0.72,0.49])
R14.errstatplus = array([0.60,0.27,0.31,0.45,0.95])
R14.errstatminus= array([0.26,0.19,0.22,0.30,0.38])
R14.errsysplus  = array([0.12,0.23,0.34,0.50,0.45])
R14.errsysminus = array([0.35,0.19,0.23,0.28,0.24])
# Pre-bugfix : (as originally submitted)
#R14.ndet        = array([1.33,6.53,6.55,5.44,1.40])
#R14.rate        = array([0.32,0.46,0.50,0.71,0.55])
#R14.errstatplus = array([0.59,0.26,0.28,0.45,0.97])
#R14.errstatminus= array([0.24,0.18,0.19,0.29,0.41])
#R14.errsysplus  = array([0.12,0.21,0.31,0.50,0.53])
#R14.errsysminus = array([0.31,0.16,0.18,0.30,0.29])
R14.z           = array([0.25,0.75,1.25,1.75,2.25])
R14.zerrplus    = array([0.25,0.25,0.25,0.25,0.25])     
R14.zerrminus    = array([0.25,0.25,0.25,0.25,0.25])    
R14.reference="CANDELS"
R14.ref="This Work [HST+WFC3]"
R14.marker= 's'
R14.mfc='darkorange'
R14.mec='darkred'
R14.color='darkred'
R14.ms= 13


R14h=rate()  # Full CANDELS sample, No Decliners, using galsnid for priors

R14h.ndet        = array([1.59,9.18,9.40,5.95,1.57])
R14h.rate        = array([0.39,0.65,0.71,0.77,0.62])
R14h.errstatplus = array([0.61,0.29,0.32,0.47,0.99])
R14h.errstatminus= array([0.27,0.21,0.23,0.31,0.44])
R14h.errsysplus  = array([0.07,0.20,0.36,0.50,0.48])
R14h.errsysminus = array([0.39,0.23,0.23,0.25,0.26])
# Pre-bugfix : (as originally submitted)
#R14h.ndet        = array([1.43,8.56,7.00,6.04,1.75])
#R14h.rate        = array([0.35,0.61,0.53,0.79,0.69])
#R14h.errstatplus = array([0.60,0.29,0.29,0.47,1.01])
#R14h.errstatminus= array([0.26,0.20,0.20,0.31,0.47])
#R14h.errsysplus  = array([0.08,0.19,0.28,0.50,0.52])
#R14h.errsysminus = array([0.34,0.20,0.16,0.25,0.30])
R14h.z           = array([0.25,0.75,1.25,1.75,2.25])
R14h.zerrplus    = array([0.25,0.25,0.25,0.25,0.25])     
R14h.zerrminus    = array([0.25,0.25,0.25,0.25,0.25])    
R14h.reference="CANDELS"
R14h.ref="This Work [HST+WFC3]"
R14h.marker= 's'
R14h.mfc='darkorange'
R14h.mec='darkred'
R14h.color='darkred'
R14h.ms= 8

FF=rate()  # HST FrontierSN sample (including GLASS and ground imaging
#FF.rate        = array([0.34,0.67,0.75,0.52,0.23,0.16,0.10])  # assuming a decline after z=1.5
FF.rate        = array([0.34,0.67,0.96,0.97,0.77,0.58,0.44])  # assuming a simple t-1 power law
FF.z           = array([0.25,0.75,1.25,1.75,2.25,2.75,3.25])
FF.zerrplus    = array([0.25,0.25,0.25,0.25,0.25,0.25,0.25])     
FF.zerrminus   = array([0.25,0.25,0.25,0.25,0.25,0.25,0.25])    
# Anticipated SNIa yields from the 3 FFSN components : 
Ndet_HFF   = array([  0.7, 2.3, 2.5, 2.5, 1.0, 0.3, 0.1 ] ) # from SNANA sims
Ndet_GLASS = array([  0.6, 3.8, 4.2, 1.3, 0.1, 0.0, 0.0 ] ) # back-of-the-envelope
Ndet_GROUND= array([  1.3, 4.6, 6.0, 2.8, 0.0, 0.0, 0.0 ] ) # basically an educated guess
FF.ndet = Ndet_HFF + Ndet_GLASS + Ndet_GROUND 
FF.errstatplus = FF.rate * poissonerrExact( FF.ndet )[0] / FF.ndet
FF.errstatminus= FF.rate * poissonerrExact( FF.ndet )[1] / FF.ndet
FF.derrstat    = (FF.errstatplus+FF.errstatminus) /2.
FF.errsysplus  = FF.rate * 0.3
FF.errsysminus = FF.rate * 0.3
FF.reference="FrontierSN"
FF.ref="FFSN"
FF.marker='D'
FF.color='darkcyan'
FF.mec='darkcyan'
FF.mfc='darkcyan'
FF.ms= 10

RT10=rate()              
RT10.ndet = array([1.95,4.01,5.11,6.49,10.09,14.29,15.4,13.2,11.01])
RT10.z = array([ 0.15, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05])
RT10.zerrplus  = array([ 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ])
RT10.zerrminus = array([ 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ])
RT10.rate = array([0.32,0.34,0.31,0.32,0.49,0.68,0.78,0.76,0.79])
RT10.errstatplus = array([0.23,0.19,0.15,0.14,0.17,0.21,0.22,0.25,0.28])
RT10.errsysplus = array([0.07,0.07,0.12,0.07,0.14,0.23,0.31,0.32,0.36])
RT10.errstatminus = array([0.23,0.19,0.15,0.14,0.17,0.21,0.22,0.25,0.28])
RT10.errsysminus = array([0.06,0.03,0.04,0.07,0.08,0.14,0.16,0.26,0.41])
RT10.count = array([1.95 ,4.01 ,5.11 ,6.49 ,10.09,14.29,15.43,13.21,11.01])
RT10.reference="Rodney and Tonry 2010"
RT10.ref="Rodney \& Tonry 2010 (IfA)"
RT10.marker='o'
RT10.color='gray'
RT10.mec='gray'           
RT10.mfc='w'
RT10.ms=8

S03 = rate()
S03.z           =array([0.11])
S03.zerrplus    =array([0.05])
S03.zerrminus   =array([0.05])
S03.rate        =array([0.37])
S03.errstatplus =array([0.1])
S03.errstatminus=array([0.1])
S03.errsysplus  =array([0.00 ])
S03.errsysminus =array([0.00 ])
S03.reference="Strolger et al. 2003"
S03.ref="Strolger+ 2003"
S03.marker='o'
S03.color='gray'
S03.mec='gray'
S03.mfc='w'
S03.ms=8

T03 = rate()
T03.ndet        =array([8])
T03.z           =array([0.46])
T03.zerrplus    =array([0.14])
T03.zerrminus   =array([0.16])
T03.rate        =array([0.48])
T03.errstatplus =array([0.17])
T03.errstatminus=array([0.17])
T03.errsysplus  =array([0.0 ])
T03.errsysminus =array([0.0 ])
T03.reference="Tonry et al. 2003"
T03.ref="Tonry+ 2003"
T03.marker='o'
T03.color='gray'
T03.mec='gray'
T03.mfc='white'
T03.ms=8



G14=rate()  # Full CLASH sample, parallel fields
# # 4 bins of width dz=0.6
G14.ndet        = array([2.0, 5.1, 3.9])
G14.z           = array([0.42, 0.94, 1.59] )
G14.zerrplus    = array([0.18, 0.26, 0.21] )     
G14.zerrminus   = array([0.42, 0.34, 0.39] )     
G14.rate        = array([0.46,0.45,0.45])
G14.errstatplus = array([0.42,0.22,0.34])
G14.errstatminus= array([0.32,0.19,0.22])
G14.errsysplus  = array([0.10,0.13,0.05])
G14.errsysminus = array([0.13,0.06,0.09])
G14.reference="CLASH"
G14.ref="Graur et al. 2013"
G14.marker= 'D'
G14.color='orchid'
G14.mec='darkorchid'
G14.mfc='orchid'
G14.ms= 13


W30=rate()  # Full WFIRST sample. 2750 SN to z=1.7, c.2030
W30.ndet        = array([ 70,  210, 400, 220, 320, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140])
W30.z           = array([0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65])
W30.zerrplus    = zeros(len(W30.z))+0.05
W30.zerrminus   = zeros(len(W30.z))+0.05
W30.rate        = ones( len(W30.z) )

W30.errstatplus = W30.rate * poissonerrExact( W30.ndet )[0] / W30.ndet
W30.errstatminus= W30.rate * poissonerrExact( W30.ndet )[1] / W30.ndet
#W30.errstatplus = sqrt( W30.ndet ) * W30.rate 
#W30.errstatminus= sqrt( W30.ndet ) * W30.rate 
W30.errsysplus  = W30.rate * 0.1
W30.errsysminus = W30.rate * 0.1
W30.reference="WFIRST-AFTA"
W30.ref="c. 2030"
W30.marker= 'D'
W30.color='b'
W30.mec='b'
W30.mfc='teal'
W30.ms= 13


W30GO=rate()  # additional WFIRST GO sample. 170 SN from z=1.5 to z=2.5, c.2030
W30GO.ndet        = array([100, 100, 100,  80,  60 ])
W30GO.z           = array([1.6, 1.8, 2.0, 2.2, 2.4 ])
W30GO.zerrplus    = zeros(len(W30GO.z))+0.1
W30GO.zerrminus   = zeros(len(W30GO.z))+0.1
W30GO.rate        = array([ 0.9, 0.9, 0.8, 0.7, 0.6 ])
W30GO.errstatplus = W30GO.rate * poissonerrExact( W30GO.ndet )[0] / W30GO.ndet
W30GO.errstatminus= W30GO.rate * poissonerrExact( W30GO.ndet )[1] / W30GO.ndet
W30GO.errsysplus  = W30GO.rate * 0.1
W30GO.errsysminus = W30GO.rate * 0.1
W30GO.reference="WFIRST-AFTA-GO"
W30GO.ref="c. 2030"
W30GO.marker= 'D'
W30GO.color='k'
W30GO.mec='k'
W30GO.mfc='0.7'
W30GO.ms= 13




#------------------------------------------------------------
# CC Rates : 

S09cc = rate()
S09cc.z           = array([0.0033 ])
S09cc.zerrplus    = array([0.0033 ])
S09cc.zerrminus   = array([0.0033 ])
S09cc.rate        = array([0.96 ]) # lower limit
S09cc.errstatplus = array([1.00 ])
S09cc.errsysplus  = array([1.00 ])
S09cc.errstatminus= array([0.1 ])
S09cc.errsysminus = array([0.2 ])
S09cc.reference="Smartt et al. 2009"
S09cc.ref="Smartt+ 2009 [CC]"
S09cc.marker='d'
S09cc.color='gray'
S09cc.mec='gray'
S09cc.mfc='white'
S09cc.ms=8

C99cc = rate()
C99cc.z           = array([0.01 ])
C99cc.zerrplus    = array([0.005 ])
C99cc.zerrminus   = array([0.005 ])
C99cc.rate        = array([0.43 ])
C99cc.errstatplus = array([0.17 ])
C99cc.errsysplus  = array([0.00 ])
C99cc.errstatminus= array([0.17 ])
C99cc.errsysminus = array([0.0 ])
C99cc.reference="Cappellaro et al. 1999"
C99cc.ref="Cappellaro+ 2009 [CC]"
C99cc.marker='s'
C99cc.color='gray'
C99cc.mec='gray'
C99cc.mfc='white'
C99cc.ms=8

L11cc = rate()
L11cc.z           = array([0.007 ])
L11cc.zerrplus    = array([0.007 ])
L11cc.zerrminus   = array([0.007 ])
L11cc.rate        = array([0.62 ])
L11cc.errstatplus = array([0.07 ])
L11cc.errsysplus  = array([0.17 ])
L11cc.errstatminus= array([0.07 ])
L11cc.errsysminus = array([0.15 ])
L11cc.reference="Li et al. 2011"
L11cc.ref="Li+ 2011 [CC]"
L11cc.marker='^'
L11cc.color='gray'
L11cc.mec='gray'
L11cc.mfc='white'
L11cc.ms=8

B08cc = rate()
B08cc.z           = array([0.21 ])
B08cc.zerrplus    = array([0.05 ])
B08cc.zerrminus   = array([0.05 ])
B08cc.rate        = array([1.15 ])
B08cc.errstatplus = array([0.43 ])
B08cc.errsysplus  = array([0.42 ])
B08cc.errstatminus= array([0.33 ])
B08cc.errsysminus = array([0.36 ])
B08cc.reference="Botticella et al. 2008"
B08cc.ref="Botticella+ 2008 [CC]"
B08cc.marker='>'
B08cc.color='gray'
B08cc.mec='gray'
B08cc.mfc='white'
B08cc.ms=8

C05cc = rate()
C05cc.z           = array([0.26 ])
C05cc.zerrplus    = array([0.05 ])
C05cc.zerrminus   = array([0.05 ])
C05cc.rate        = array([1.88 ])
C05cc.errstatplus = array([0.71 ])
C05cc.errsysplus  = array([0.00 ])
C05cc.errstatminus= array([0.58 ])
C05cc.errsysminus = array([0.00 ])
C05cc.reference="Cappellaro et al. 2005"
C05cc.ref="Cappellaro+ 2005 [CC]"
C05cc.marker='o'
C05cc.color='gray'
C05cc.mec='gray'
C05cc.mfc='white'
C05cc.ms=8

B09cc = rate()
B09cc.z           = array([0.30 ])
B09cc.zerrplus    = array([0.05 ])
B09cc.zerrminus   = array([0.05 ])
B09cc.rate        = array([1.63 ])
B09cc.errstatplus = array([0.34 ])
B09cc.errsysplus  = array([0.37 ])
B09cc.errstatminus= array([0.34 ])
B09cc.errsysminus = array([0.28 ])
B09cc.reference="Bazin et al. 2009"
B09cc.ref="Bazin+ 2009 [CC]"
B09cc.marker='<'
B09cc.color='gray'
B09cc.mec='gray'
B09cc.mfc='white'
B09cc.ms=8

G11cc = rate()
G11cc.z           = array([0.66 ])
G11cc.zerrplus    = array([0.05 ])
G11cc.zerrminus   = array([0.05 ])
G11cc.rate        = array([6.9 ])
G11cc.errstatplus = array([9.9 ])
G11cc.errsysplus  = array([0.00 ])
G11cc.errstatminus= array([5.4 ])
G11cc.errsysminus = array([0.00 ])
G11cc.reference="Graur et al. 2011"
G11cc.ref="Graur+ 2011 [CC]"
G11cc.marker='d'
G11cc.color='gray'
G11cc.mec='gray'
G11cc.mfc='white'
G11cc.ms=8

M12cc = rate()
M12cc.z           = array([0.39, 0.73 ])
M12cc.zerrplus    = array([0.11, 0.17 ])
M12cc.zerrminus   = array([0.29, 0.23 ])
M12cc.rate        = array([3.29, 6.40 ])
M12cc.errstatplus = array([3.08, 5.30 ])
M12cc.errsysplus  = array([1.98, 3.65 ])
M12cc.errstatminus= array([1.78, 3.12 ])
M12cc.errsysminus = array([1.45, 2.11 ])
M12cc.reference="Melinder et al. 2012"
M12cc.ref="Melinder+ 2012 [CC]"
M12cc.marker='s'
M12cc.color='gray'
M12cc.mec='gray'
M12cc.mfc='white'
M12cc.ms=8

B12cc = rate()
B12cc.z           = array([0.001 ])
B12cc.zerrplus    = array([0.002 ])
B12cc.zerrminus   = array([0.001 ])
B12cc.rate        = array([1.1 ])
B12cc.errstatplus = array([0.4 ])
B12cc.errsysplus  = array([0.3 ])
B12cc.errstatminus= array([0.0 ])
B12cc.errsysminus = array([0.0 ])
B12cc.reference="Botticella et al. 2012"
B12cc.ref="Botticella+ 2012 [CC]"
B12cc.marker='o'
B12cc.color='gray'
B12cc.mec='gray'
B12cc.mfc='white'
B12cc.ms=8

D12cc = rate()
D12cc.z         = array([0.39, 0.73, 1.11 ])
D12cc.zerrplus  = array([0.11, 0.17, 0.19 ])
D12cc.zerrminus = array([0.29, 0.23, 0.21 ])
#D12cc.rate    = array([2.24, 4.86, 5.95 ]) # without dust correction
D12cc.rate      = array([3.00, 7.39, 9.57 ]) # with dust correction
D12cc.errstatplus =array([1.28,1.86, 3.76 ])
D12cc.errsysplus  =array([1.04,3.20, 4.96 ])
D12cc.errstatminus=array([0.94,1.52, 2.80 ])
D12cc.errsysminus =array([0.57,1.60, 2.80 ])
D12cc.reference="Dahlen et al. 2012"
D12cc.ref="Dahlen+ 2012 [CC]"
D12cc.marker='o'
D12cc.color='gray'
D12cc.mec='gray'
D12cc.mfc='white'
D12cc.ms=8



CCRATELIST = [ C99cc, S09cc, C05cc, B08cc, B09cc, L11cc, G11cc, D12cc, M12cc, B12cc ]


# These rate lists are composed of all the non-redundant rates 
# i.e. it uses Rodney&Tonry2010 but not Barris&Tonry2006
RATELISTLOWZ = [ B04, B08, C99, D10, H00, H08, GM13, L11, P02, P12, RT10, T03, M12 ]
RATELISTHIGHZ = [ D08, G11, B12, G14, R14 ]
RATELISTGROUND = RATELISTLOWZ + [ G11 ]
RATELISTPRECC = RATELISTLOWZ + [ D08, G11, B12 ]
RATELISTALL = RATELISTLOWZ + RATELISTHIGHZ
RATELISTHST = [ D08, B12, G14, R14 ]

RATELISTCANDC = [ G14, R14 ]

def mkCOMPOSITE( ratelist=RATELISTALL, zbinedges = [0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3], 
                 debug=False ):
    """ combine all the rates into a composite rate estimate
    binned by redshift """
    if debug: import pdb; pdb.set_trace()
    COMP = rate()

    COMP.z         = [  ]
    COMP.zerrplus  = [  ]
    COMP.zerrminus = [  ]
    COMP.ndet = []
                 
    for i in range(len(zbinedges)-1) :
        COMP.z.append( (zbinedges[i] + zbinedges[i+1])/2. )
        COMP.zerrplus.append( (zbinedges[i+1] - COMP.z[-1]) )
        COMP.zerrminus.append( (COMP.z[-1]-zbinedges[i]) )

    # separate the rates, errors, and counts into redshift bins 
    COMPrbinned = [ [] for i in range(len(COMP.z)) ]
    COMPrbinnedErrPlus = [ [] for i in range(len(COMP.z)) ]
    COMPrbinnedErrMinus =[ [] for i in range(len(COMP.z)) ]
    COMPrbinnedSysErrPlus = [ [] for i in range(len(COMP.z)) ]
    COMPrbinnedSysErrMinus =[ [] for i in range(len(COMP.z)) ]
    COMPrbinnedNdet = [ [] for i in range(len(COMP.z)) ]
    for R in ratelist :
        iRz = digitize( R.z, zbinedges )
        for ibin,iR in zip(iRz,range(len(R.rate)) ) : 
            if ibin>len(COMP.z) : continue
            if ibin<1 : continue
            COMPrbinned[ibin-1].append( R.rate[iR] )
            COMPrbinnedErrPlus[ibin-1].append( R.errplus[iR] )
            COMPrbinnedErrMinus[ibin-1].append( R.errminus[iR] )
            COMPrbinnedSysErrPlus[ibin-1].append( R.errsysplus[iR] )
            COMPrbinnedSysErrMinus[ibin-1].append( R.errsysminus[iR] )
            if 'ndet' in R.__dict__ : 
                COMPrbinnedNdet[ibin-1].append( R.ndet[iR] ) 
            elif 'Ndet' in R.__dict__ : 
                COMPrbinnedNdet[ibin-1].append( R.Ndet[iR] ) 
            else : 
                # number of detections not yet defined in ratetable.py 
                # for ground-based surveys
                COMPrbinnedNdet[ibin-1].append( 10. )

    # In each z-bin, compute the mean rates, composite poisson errors, 
    # and count-weighted systematic error sums
    COMP.rate = []
    COMP.errstatplus = []
    COMP.errstatminus = []
    COMP.errsysplus = []
    COMP.errsysminus = []
    for irbin in range(len(COMPrbinned)):
        Ninbin = len( COMPrbinned[irbin] )
        if not Ninbin : continue
        weights = 1. / ( array(COMPrbinnedErrPlus[irbin])**2 + array(COMPrbinnedErrMinus[irbin])**2 )
        COMP.rate.append( average( COMPrbinned[irbin], weights=weights ) )
        # standard err on the mean (NOTE: This is not statistically proper!! 
        #  a better approach is to approximate the asymmetric log-likelihood,
        #  see Roger Barlow's paper : http://arxiv.org/pdf/physics/0306138.pdf
        # COMP.errstatplus.append( 1 / sqrt( (1./array(COMPrbinnederrPlus[irbin])**2).sum() ) )
        # COMP.errstatminus.append( 1 / sqrt( (1./array(COMPrbinnederrMinus[irbin])**2).sum() ) )

        # recompute the statistical error as the poisson error of the combined sample
        Ndet = sum(COMPrbinnedNdet[irbin])
        COMP.ndet.append( Ndet )
        pErrPlus, pErrMinus = poissonerrExact( Ndet )
        COMP.errstatplus.append( COMP.rate[-1] * pErrPlus / Ndet )
        COMP.errstatminus.append( COMP.rate[-1] * pErrMinus / Ndet  )

        # compute the systematic error as the count-weighted sum of systematic errors
        weights = array( COMPrbinnedNdet[irbin] ) / sum(COMPrbinnedNdet[irbin])
        COMP.errsysplus.append( average( COMPrbinnedSysErrPlus[irbin], weights=weights ) ) 
        COMP.errsysminus.append( average( COMPrbinnedSysErrMinus[irbin], weights=weights ) ) 
    COMP.z = array( COMP.z )
    COMP.zerrplus = array( COMP.zerrplus )
    COMP.zerrminus = array( COMP.zerrminus )
    COMP.ndet = array( COMP.ndet )
    COMP.rate = array( COMP.rate )
    COMP.errstatplus = array( COMP.errstatplus )
    COMP.errstatminus = array( COMP.errstatminus )
    COMP.errsysplus = array( COMP.errsysplus )
    COMP.errsysminus = array( COMP.errsysminus )
    COMP.reference="All surveys"
    COMP.ref="All"
    COMP.marker='o'
    COMP.color='0.8'
    COMP.mec='0.8'
    COMP.mfc='w'
    COMP.ms=12
    return( COMP )

GROUND = mkCOMPOSITE(ratelist=RATELISTGROUND, zbinedges=[0.1,0.3,0.5,0.7,0.9,1.1,1.4,2.0])
GROUND.reference="Ground-based surveys"
GROUND.ref="Ground Surveys"
GROUND.color='saddlebrown'
GROUND.mfc='gold'
GROUND.mec='saddlebrown'
GROUND.reference="All completed surveys (ground+HST)"

RATELISTLOWZ2 = [ B04, B08, C99, D10, H00, H08, GM13, L11, M12, P02, P12, RT10, T03 ]
#RATELISTLOWZ2 = [ D10, P12 ]
LOWZ = mkCOMPOSITE(ratelist=RATELISTLOWZ2, zbinedges=[0.1,0.3,0.5,0.7,0.9,1.1])
LOWZ.reference="Ground-based surveys"
LOWZ.ref="Ground Surveys"
LOWZ.color='darkmagenta'
LOWZ.mfc='gold'
LOWZ.mec='darkmagenta'
LOWZ.reference="All completed surveys (ground+HST)"


ALL = mkCOMPOSITE(ratelist=RATELISTALL, zbinedges = [0.1,0.3,0.5,0.7,0.9,1.1,1.4,2.0])
ALL.color='saddlebrown'
ALL.mfc='gold'
ALL.mec='saddlebrown'
ALL.reference="All completed surveys (ground+HST)"

ALLPRECC = mkCOMPOSITE(ratelist=RATELISTPRECC, zbinedges = [0.1,0.3,0.5,0.7,0.9,1.1,1.4,1.8])
ALLPRECC.color='darkgreen'
ALLPRECC.reference="All completed surveys (ground+HST) prior to C+C"
   

ALLHST = mkCOMPOSITE(ratelist=[D08,B12,R14,G14], zbinedges = [0.2,0.6,1.0,1.4,1.8])
ALLHST.color='darkblue'
ALLHST.mfc='cyan'
ALLHST.mec='darkblue'
ALLHST.reference="All HST surveys (D08,B12,G14,R14)"

ALLCC = mkCOMPOSITE(ratelist=[R14,G14], zbinedges = [0.2,0.6,1.0,1.4,1.8])
ALLHST.color='darkmagenta'
ALLHST.mfc='darkcyan'
ALLHST.mec='darkmagenta'
ALLHST.reference="CANDELS+CLASH (R14,G14)"


#------------------------------------------------------------
# Projections for the final constraints including frontier fields

HIZHST = mkCOMPOSITE(ratelist=RATELISTHST, zbinedges = [0.0,0.5,1.0,1.5,2.0,2.5])
HIZHST.reference="HST surveys at z>1"
HIZHST.mec='k'  
HIZHST.mfc='g'  
HIZHST.color='g'
HIZHST.marker='s'

HIZFF = mkCOMPOSITE(ratelist=RATELISTHST+[FF], zbinedges = [0.0,0.5,1.0,1.5,2.0,2.5])
HIZFF.z = array([ 0.25, 0.75, 1.25, 1.75, 2.25])
HIZFF.ndet = array( [3 + 2+2+3, 8 + 3+7+5, 20+4+8+8, 3+2+6+4, 1.5+2])
HIZFF.errstatplus = HIZFF.rate * poissonerrExact( HIZFF.ndet )[0] / HIZFF.ndet
HIZFF.errstatminus= HIZFF.rate * poissonerrExact( HIZFF.ndet )[1] / HIZFF.ndet
HIZFF.errsysplus =  array(HIZFF.rate) * array( [ 0.2, 0.2, 0.3, 0.5, 0.7 ] )
HIZFF.errsysminus = array(HIZFF.rate) * array( [ 0.2, 0.2, 0.3, 0.4, 0.5 ] )
#HIZFF.rate = array([0.67, 0.37, 0.2])
HIZFF.zerrminus = array([0.25, 0.25, 0.25, 0.25, 0.25])
HIZFF.zerrplus = array([0.25, 0.25, 0.25, 0.25, 0.25])
HIZFF.reference="High-z HST surveys (incl. FFSN projection) at z>1"
HIZFF.marker='D'
HIZFF.color='g'    
HIZFF.mec='forestgreen'  
HIZFF.mfc='teal'      


HIZFF20 = deepcopy( HIZFF )
HIZFF20.rate = array( [ 0.35, 0.57, 0.66, 0.54, 0.28 ] ) 
HIZFF20.errstatminus[-1] = HIZFF20.rate[-1]*0.9
HIZFF20.elw=2
HIZFF20.color='r'    
HIZFF20.mec='darkred'  
HIZFF20.mfc='palevioletred'      

HIZFF80 = deepcopy( HIZFF )
HIZFF80.rate = array( [ 0.35, 0.65, 1.1, 1.2, 0.93 ] ) 
HIZFF80.elw=2
HIZFF80.color='cadetblue'    
HIZFF80.mec='darkblue'  
HIZFF80.mfc='cadetblue'      
