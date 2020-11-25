from matplotlib import ticker
from matplotlib import pyplot as pl
from pytools import colorpalette, plotsetup
import numpy as np

#import astropy.units as u
#from astropy.cosmology import FlatLambdaCDM
#from scipy.interpolate import interp1d
#LCDM = FlatLambdaCDM( 70, 0.3 )
#zsamples = np.logspace( -3, 0.7, num=50, base=10 )
#tsamples = LCDM.age( zsamples )
#getzfromt = interp1d( tsamples, zsamples, bounds_error=False, fill_value=tsamples[-1], assume_sorted=True )
#def zfromt( t ):
#    return( getzfromt( t ) )

from pytools import cosmo
zfromt = cosmo.zfromt

import os
import sys

thisfile = sys.argv[0]
if 'ipython' in thisfile : thisfile = __file__
_thisdir = os.path.dirname( os.path.abspath( thisfile ) )


from pytools import colorpalette as cp

def ratepointsSeq1( zmax=2.55,  ageaxis=True, slide=False ):
    """
    Making a figure for a review talk showing how the observed rates vs z
    have improved since c. 2009
    """
    from ratetable import B04,B08,BT06,C99,D08,H00,H08,M03,N06,P02,P07,S03,T03
    plotsetup.presfig(figsize=[15,8])
    pl.clf()
    ax = pl.gca()

    ratelist2009 = [ B04,B08,BT06,C99,D08,H00,H08,M03,N06,P02,P07,S03,T03 ]

    for R in ratelist2009 :
        c,m,ms='0.5','o',8
        if slide and R is BT06:
            c,m,ms=cp.red,'^',10
            R.z -=0.02
        if slide and R is D08:
            c,m,ms=cp.green,'D',10
        R.mfc=c
        R.mec=c
        R.color=c
        R.ms=ms
        R.marker=m
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100, label=R.ref)

    ax.text(0.95,0.9,'SN Ia Rate measures c.2009',fontsize='x-large',ha='right',va='top', transform=ax.transAxes)
    # ax.legend( numpoints=1, handlelength=0.2, loc='upper right', title='Rate Measures c.2009')

    if slide==1 :
        ax.text( 1.35, 1.7, 'HST GOODS\n Dahlen+ 2008', color=cp.green, fontsize='x-large')
        ax.text( 0.6, 2.2, 'IfA Deep\n Barris+Tonry 2006', color=cp.red, fontsize='x-large')

    if slide>1 :
        import dtdconvolve
        z = np.arange( 0.01, zmax+0.05, 0.05 )
        snr = dtdconvolve.snrate( z, DTDmodel='gauss', normalize=False,
                                  eta=0.0005, tau=3.4,  SFRmodel='madau??')
        ax.plot( z, snr, color=cp.green, ls='-', marker=' ' )
        ax.text( 0.6, 2.2, 'IfA Deep\n Barris+Tonry 2006', color=cp.red, fontsize='x-large')
        ax.text( 1.35, 1.7, 'Gaussian DTD, $\\tau \sim 3$ Gyr\nDahlen+ 2008', color=cp.green, fontsize='x-large')

    addaxislabels( zmax=zmax )
    if ageaxis : ax = addageaxis()



def ratepointsSeq2( zmax=2.55,  ageaxis=True, slide=0, lines=False, **kwargs ):
    """
    Making a sequence of figures for a review talk showing how the observed rates vs z
    have improved since c. 2009
    """
    from ratetable import B04,B08,B12,BT06,C99,D04,D08,D10,G11,G14,GM13,H00,H08,K08,L11,LOWZ,M03,M12,N06,P02,P07,P12,R14,RT10,S03,T03
    plotsetup.presfig(figsize=[15,8])
    pl.clf()
    ax = pl.gca()


    ratelistOld = [ B04,B08,BT06,C99,D08,H00,H08,M03,N06,P02,P07,S03,T03 ]
    if slide==2 :
        ratelistOld = [ B04,B08,C99,D08,H00,H08,M03,N06,P02,P07,RT10,S03,T03 ]

    if slide==3 :
        ratelistOld = [ B04,B08,C99,D08,D10,H00,G11,M03,P02,P07,RT10,S03,T03 ]
    if slide==4 :
        ratelistOld = [ B04,B08,C99,D08,D10,H00,G11,M03,P02,RT10,S03,T03 ]
    if slide==5 :
        P12.z-=0.01
        ratelistOld = [ B04,B08,C99,D08,D10,H00,G11,GM13,P02,P12,RT10,S03,T03 ]
    if slide==6 :
        P12.z-=0.01
        ratelistOld = [ B04,B08,C99,D08,D10,H00,G11,GM13,G14,P02,P12,RT10,R14,S03,T03 ]



    for R in ratelistOld :
        fc,ec,m,ms='0.5','0.5','o',8
        if R is BT06:
            fc,m,ms='w','^',10
            R.z -=0.02
        if R is D08 and slide==0:
            fc,m,ms='w','D',10
        if R is P07 and slide==3:
            fc,m,ms='w','d',10
        R.mfc=fc
        R.mec=ec
        R.color=ec
        R.ms=ms
        R.marker=m
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100, label=R.ref)


    ratelistNew=[]
    if slide==1 :
        ratelistNew = [RT10]
        ec,fc,m,ms=cp.red,cp.red,'^',10
        ax.text( 1.0, 2.05, 'Rodney+Tonry 2010', color=cp.red, fontsize='x-large')

    if slide==2 :
        ratelistNew = [D10]
        ec,fc,m,ms=cp.bluegray,cp.bluegray,'v',10
        ax.text( 0.2, 1.5, 'SDSS\n Dilday+ 2010', color=cp.bluegray, fontsize='x-large')

    if slide==3 :
        ratelistNew = [G11]
        ec,fc,m,ms=cp.teal,cp.teal,'d',10
        ax.text( 1.3, 2.2, 'SDF\n Graur+ 2011', color=cp.teal, fontsize='x-large')


    if slide==4 :
        P12.z-=0.01
        ratelistNew = [P12]
        ec,fc,m,ms=cp.purple,cp.purple,'s',10
        ax.text( 0.6, 2.2, 'SNLS\n Perrett+ 2012', color=cp.purple, fontsize='x-large')

    if slide==5 :
        ratelistNew = [G14,R14]
        ec,fc,m,ms=cp.darkgold,cp.darkgold,'D',10
        ax.text( 0.6, 2.2, 'CLASH+CANDELS\n Graur+ 2014, Rodney+ 2014', color=cp.gold, fontsize='x-large')

    if slide==6 :
        ratelistNew = []
        ax.text(0.95,0.9,'SN Ia Rate measures c.2014',fontsize='x-large',ha='right',va='top', transform=ax.transAxes)


    for R in ratelistNew :
        R.mfc=fc
        R.mec=ec
        R.color=ec
        R.ms=ms
        R.marker=m
        R.zerrplus=0*R.zerrplus
        R.zerrminus=0*R.zerrminus
        R.plot( thicksys=False, zorder=-100, label=R.ref)

    # ax.legend( numpoints=1, handlelength=0.2, loc='upper right', title='Rate Measures c.2009')

    if slide==6 :
        import dtdconvolve
        z = np.arange( 0.01, zmax+0.05, 0.05 )
        snr = dtdconvolve.snrate( z, DTDmodel='t-1', normalize=True,
                                  SFRmodel='B13')
        ax.plot( z, snr, color=cp.green, ls='-', marker=' ' )
        ax.text( 2., 1., 'DTD $\propto$ \ $t^{-1}$', color=cp.green, fontsize='x-large')

    addaxislabels( zmax=zmax )
    if ageaxis : ax = addageaxis()

    if slide==6 :
        ax.set_ylim(0,1.98)
        pl.draw()


# best-fit eta and eta err
etaBest = {
    't-1': [1.45e-13,3.27157e-15],
    't-1fp': [1.21770e-13,1.50598e-14],

    'BPS.DD.Ruiter'          : [   1.952 ,   0.044 ],
    'BPS.DD.Bours/Toonen'    : [   5.317 ,   0.120 ],
    'BPS.DD.Wang/Han'        : [   1.884 ,   0.042 ],
    'BPS.DD.Claeys'          : [   2.314 ,   0.052 ],
    'BPS.DD.Mennekens'       : [   4.704 ,   0.106 ],
    'BPS.DD.Yungelson'       : [   2.794 ,   0.063 ],

    'BPS.SD.Ruiter'          : [  46.200 ,   1.041 ],
    'BPS.SD.Bours/Toonen'    : [  13.028 ,   0.294 ],
    'BPS.SD.Wang/Han'        : [   2.938 ,   0.066 ],
    'BPS.SD.Claeys'          : [  67.186 ,   1.515 ],
    'BPS.SD.Mennekens'       : [   2.722 ,   0.061 ],
    'BPS.SD.Yungelson'       : [ 1746.575 ,  39.371 ],

    'BPS.SD+DD.Ruiter'       : [   1.873 ,   0.042 ],
    'BPS.SD+DD.Bours/Toonen' : [   3.776 ,   0.085 ],
    'BPS.SD+DD.Wang/Han'     : [   1.148 ,   0.026 ],
    'BPS.SD+DD.Claeys'       : [   2.237 ,   0.050 ],
    'BPS.SD+DD.Mennekens'    : [   1.724 ,   0.039 ],
    'BPS.SD+DD.Yungelson'    : [   2.790 ,   0.063 ],

    }

fpBest = {'t-1fp': [0.554,0.074]}



def plot_rate_points( dataset='2014', zmax=2.55,  ageaxis=True, showerr=True,
                      logscale=False, systematicshift=None ):
    """
    Making a figure for a review talk showing all observed SNIa rates c.2014
    """
    import dtdconvolve
    from ratetable import B04,B08,B12,BT06,C99,D04,D08,D10,G11,G14,GM13,H00,H08,K08,L11,LOWZ,M03,M12,N06,P02,P07,P12,R14,RT10,S03,T03
    plotsetup.presfig(figsize=[15,8])
    pl.clf()
    ax = pl.gca()

    if dataset.lower() in ['2014','all'] :
        ratelist = [ B04,B08,B12,C99,D08,D10,H00,G11,GM13,G14,L11,M12,P02,P12,RT10,R14,S03,T03 ]
    elif dataset.lower() in ['2012'] :
        ratelist = [ B04,B08,C99,B12,D08,D10,H00,G11,GM13,L11,M12,P02,P12,RT10,S03,T03 ]
    elif dataset.lower() in ['2012ground'] :
        ratelist = [ B04,B08,C99,D10,H00,G11,GM13,L11,M12,P02,P12,RT10,S03,T03 ]

    if logscale : scale=1e-4
    else : scale = 1

    zlist = np.array([])
    for R in ratelist :
        fc,ec,m,ms='0.5','0.5','o',10
        R.mfc='0.8'
        R.mec='0.5'
        R.color=ec
        R.ms=ms
        R.marker=m
        R.elw=2.5

        zshift = np.zeros( len(R.z) )
        for i in range(len(R.z)) :
            while np.any( [ np.abs(R.z[i]+zshift[i]-z)<0.01 for z in zlist] ) :
                zshift[i]+=0.008
        #if R == P12 :
        #    zshift = -0.02
        #elif R==RT10 :
        #    zshift= +0.02
        zlist = np.append( zlist, R.z+zshift )
        R.plot( showerr=showerr, zorder=-100, label=R.ref,
                scalerates=scale, showzerr=False, zshift=zshift,
                systematicshift=systematicshift, alpha=1)

    if logscale :
        ax.set_yscale( 'log', basey=10, nonposy='clip'  )
    else :
        ax.set_ylim(0,1.5)
    addaxislabels( zmax=zmax, logscale=logscale )
    if ageaxis : ax = addageaxis()
    pl.draw()
    return( ax )


grouplist = [ 'Ruiter', 'Bours/Toonen','Wang/Han',  'Claeys', 'Mennekens',  'Yungelson'  ]
colorlist = [  cp.red,    cp.purple,   cp.darkblue,   cp.gold,  cp.green,    cp.lightblue ]
reflist =   [ 'Ruiter, Belcynski \& Fryer 2009', 'Bours, Toonen \& Nelemans 2013','Wang \& Han 2012', 'Claeys et al. 2014','Mennekens et al. 2010','Yungelson 2010' ]


def bpsdtd_fits_for_chicago( ):
    outlinelist = []
    # for modeltype in ['DD','SD','SD+DD'] :
    for modeltype in ['SD+DD'] :
        for group in grouplist :
            dtdmodel='BPS.%s.%s'%(modeltype,group)
            etabest, etaerr = fitSNR( data='ALL', dtdmodel=dtdmodel,
                                      systematicshift='mid', sfrmodel='B13' )

            modelkey = "%-25s:"% ( "'%s'"%dtdmodel )
            outline = "%s [ %7.3f ,  %6.3f ],"%( modelkey, etabest, etaerr )
            print( outline )
            outlinelist.append( outline )
    for outline in outlinelist:
        print( outline )
    return( outlinelist )


def rateplots_for_chicago( figlist=[1], zmax=2.5, dataset='2014' ):
    from pytools import plotsetup
    z = np.arange( 0.01, zmax+0.05, 0.05 )

    if 1 in figlist :
        # just the data points, showing composite sys+stat error
        plotsetup.presfig( 1, figsize=[15,8])
        ax = plot_rate_points( dataset=dataset, zmax=zmax, ageaxis=True, showerr=True, logscale=False)
        ax.set_ylim(0,1.5)
        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_%s.png'%dataset))


    if 2 in figlist :
        # t-1 DTD fit, using statistical errors
        plotsetup.presfig( 2, figsize=[15,8])
        ax = plot_rate_points( zmax=zmax, ageaxis=True, showerr='stat', logscale=False)
        eta, etaerr  = etaBest['t-1']
        snrt1 = snrate( z, dtdmodel='t-1', sfrmodel='B13', fitpar=eta )

        ax.plot( z, snrt1, color=cp.red, ls='--', marker=' ' )
        ax.text( 2., 1., 'DTD $\propto$ \ $t^{-1}$', ha='center',color=cp.red, fontsize='x-large')
        snrt1a = snrate( z, dtdmodel='t-1', sfrmodel='B13', fitpar=eta+etaerr )
        snrt1b = snrate( z, dtdmodel='t-1', sfrmodel='B13', fitpar=eta-etaerr )
        ax.fill_between( z, snrt1b, snrt1a, color=cp.red, alpha=0.3 )
        ax.set_ylim(0,1.5)
        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_dtd_t-1.png'))

    for fig in [3,4,5,6,7,8,9]:
        # scaled fits to the DD, SD and SD+DD models
        if fig not in figlist : continue
        plotsetup.presfig( fig, figsize=[15,8])
        ax = plot_rate_points( zmax=zmax, ageaxis=True, showerr='stat', logscale=False)

        for group, color, ref in zip( grouplist, colorlist, reflist ) :
            if fig in [3,6] : bpsmod = 'BPS.DD.' + group
            elif fig in [4,7] : bpsmod = 'BPS.SD.' + group
            elif fig in [5,8,9] : bpsmod = 'BPS.SD+DD.' + group
            if fig < 6 : eta, etaerr = 1,0
            else : eta,etaerr = etaBest[ bpsmod ]

            if etaerr :
                snrbpsa = snrate( z, dtdmodel=bpsmod, fitpar=eta+etaerr  )
                snrbpsb = snrate( z, dtdmodel=bpsmod, fitpar=eta-etaerr  )
                ax.fill_between( z, snrbpsb, snrbpsa, color=color, alpha=0.3 )
            else :
                snrbps = snrate( z, dtdmodel=bpsmod, fitpar=1  )
                ax.plot( z, snrbps, color=color, ls='-', marker=' ', label=ref)

            modfam = bpsmod.split('.')[1]
            if fig!=9 :
                ax.text(2., 1., modfam + ' Models', ha='center',color='k', fontsize='x-large')
        ax.set_ylim(0,1.5)
        pl.draw()
        if fig in [3,4,5] : pngfile = 'snr_dtd_bps_%s.png'%bpsmod.split('.')[1].lower()
        if fig in [6,7,8] : pngfile = 'snr_dtd_bps_%s_fit.png'%bpsmod.split('.')[1].lower()
        if fig in [9] : pngfile = 'snr_dtd_bps_%s_fit_nolabel.png'%bpsmod.split('.')[1].lower()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/' + pngfile ) )


    if 10 in figlist :
        # Showing different fp tracks
        plot_rate_points( '2014' )
        plotsetup.presfig( 10, figsize=[15,8])
        ax = plot_rate_points( zmax=zmax, ageaxis=True, showerr=True, logscale=False)
        eta, etaerr  = etaBest['t-1fp']
        fp, fperr  = fpBest['t-1fp']
        snr02 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*1.6,0.2] )
        snr05 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta,fp] )
        snr08 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*0.5,0.8] )
        ax.fill_between( z, snr08*0.98, snr08*1.02, color=cp.lightblue, alpha=0.5 )
        ax.fill_between( z, snr05*0.98, snr05*1.02, color=cp.green, alpha=0.5 )
        ax.fill_between( z, snr02*0.98, snr02*1.02, color=cp.red, alpha=0.5 )
        ax.set_ylim(0,1.5)

        ax.text(2., 1.4, 'Prompt\nFraction', ha='center',va='top', color='k', fontsize='x-large')
        ax.text(2., 1.1, '0.8', ha='center',va='center', color=cp.darkbluegray, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.76, '0.5', ha='center',va='center', color=cp.darkgreen, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.43, '0.2', ha='center',va='center', color=cp.darkred, backgroundcolor='w', fontsize='x-large')

        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_dtd_fprompt.png' ) )

    if 11 in figlist :
        # Showing future constraints from HST with FFSN
        plotsetup.presfig( 11, figsize=[15,8])
        ax = plot_rate_points( '2012ground', zmax=zmax, ageaxis=True, showerr='stat', logscale=False)
        eta, etaerr  = etaBest['t-1fp']
        fp, fperr  = fpBest['t-1fp']
        snr02 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*1.6,0.2] )
        snr05 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta,fp] )
        snr08 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*0.5,0.8] )
        ax.fill_between( z, snr08*0.98, snr08*1.02, color=cp.lightblue, alpha=0.5 )
        ax.fill_between( z, snr05*0.98, snr05*1.02, color=cp.green, alpha=0.5 )
        ax.fill_between( z, snr02*0.98, snr02*1.02, color=cp.red, alpha=0.5 )
        ax.text(2., 1.4, 'Prompt\nFraction', ha='center',va='top', color='k', fontsize='x-large')
        ax.text(2., 1.1, '0.8', ha='center',va='center', color=cp.darkbluegray, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.76, '0.5', ha='center',va='center', color=cp.darkgreen, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.43, '0.2', ha='center',va='center', color=cp.darkred, backgroundcolor='w', fontsize='x-large')

        from ratetable import HIZFF80, HIZFF20
        HIZFF80.plot( zshift=-0.03, showerr='stat' )
        HIZFF20.plot( zshift=-0.01, showerr='stat' )
        ax.set_ylim(0,1.5)
        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_FFSN_stat.png' ) )


    if 12 in figlist :
        plotsetup.presfig( 11, figsize=[15,8])
        ax = plot_rate_points( '2012ground', zmax=zmax, ageaxis=True, showerr=True, logscale=False)
        eta, etaerr  = etaBest['t-1fp']
        fp, fperr  = fpBest['t-1fp']
        snr02 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*1.6,0.2] )
        snr05 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta,fp] )
        snr08 = snrate( z, dtdmodel='t-1fp', sfrmodel='B13', fitpar= [eta*0.5,0.8] )
        ax.fill_between( z, snr02*0.98, snr02*1.02, color=cp.lightblue, alpha=0.5 )
        ax.fill_between( z, snr05*0.98, snr05*1.02, color=cp.green, alpha=0.3 )
        ax.fill_between( z, snr08*0.98, snr08*1.02, color=cp.red, alpha=0.3 )
        ax.set_ylim(0,1.5)
        pl.draw()

        # same as 11, but with systematic errors included
        from ratetable import HIZFF80, HIZFF20
        HIZFF80.plot( zshift=-0.03, showerr='true' )
        HIZFF20.plot( zshift=-0.01, showerr='true' )
        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_FFSN_stat+syst.png' ) )

    if 13 in figlist :
        # PS1 and DES contributions at low-z
        print('no fig 10 yet')

    if 14 in figlist :
        # HSC contributions at mid-z
        print('no fig 11 yet')

    if 15 in figlist or 16 in figlist :
        plotsetup.presfig( 16, figsize=[15,8])
        # WFIRST contributions at high z
        from ratetable import W30, HIZFF, W30GO

        for sample, color, ls in zip( [HIZFF], [cp.darkgreen], ['-','-','-','-']) :
            z = np.append( np.append( [sample.z[0]-0.1], sample.z), [sample.z[-1]+2*sample.zerrplus[-1]] )
            ndet = np.append( np.append( [0], sample.ndet), [0] )
            pl.plot( z, ndet, drawstyle='steps-post', color=color, ls=ls )

        ax = pl.gca()
        zjwst = [ 2.0, 3.8 ]
        njwst0 = [ 0, 0 ]
        njwst = [ 8, 8 ]
        ax.fill_between( zjwst, njwst0, njwst, color=cp.darkgold, zorder=-100 )
        ax.set_ylim( 0, 50 )
        ax.set_xlim( 0.01, 3.6 )
        ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.5 ) )
        ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.25 ) )
        ax.yaxis.set_major_locator( ticker.MultipleLocator( 10 ) )
        ax.yaxis.set_minor_locator( ticker.MultipleLocator( 2 ) )
        ax.set_xlabel('redshift')
        ax.set_ylabel('Number of SN Detections')

        if 16 not in figlist :
            ax.text( 2.0, 30, 'All HST Surveys\n c.2016', ha='left', va='bottom', color=cp.darkgreen, fontsize='x-large', fontweight='bold' )
            ax.text( 2.75, 12, 'JWST range', ha='center', va='bottom', color=cp.darkgold, fontsize='x-large', fontweight='bold'  )
            pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/jwst_bar.png' ) )

        if 16 in figlist :
            for sample, color, ls in zip( [W30, HIZFF, W30GO], [cp.darkblue,cp.darkgreen,cp.darkred,cp.gold], ['-','-','-','-']) :
                z = np.append( np.append( [sample.z[0]-0.1], sample.z), [sample.z[-1]+0.1] )
                ndet = np.append( np.append( [0], sample.ndet), [0] )
                pl.plot( z, ndet, drawstyle='steps-post', color=color, ls=ls )

            ax.text( 1.25, 60, 'HST', ha='right', va='bottom', color=cp.darkgreen, fontsize='x-large', fontweight='bold' )

            ax.fill_between( zjwst, njwst0, [40,40], color=cp.darkgold, zorder=-100 )
            ax.text( 3.15, 50, 'JWST', ha='center', va='bottom', color=cp.darkgold, fontsize='x-large', fontweight='bold'  )

            ax.text( 0.9,  180, 'WFIRST Baseline', ha='left', va='bottom', color=cp.darkblue, fontsize='x-large', fontweight='bold'  )
            ax.text( 2.0, 120, 'WFIRST GO', ha='left', va='bottom', color=cp.darkred, fontsize='x-large', fontweight='bold'  )

            ax.set_ylim( 0, 420 )
            ax.yaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
            ax.yaxis.set_minor_locator( ticker.MultipleLocator( 20 ) )
            pl.draw()
            pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/wfirst_histogram.png' ) )



    if 20 in figlist :
        # Showing different fp tracks
        plot_rate_points( '2014' )
        plotsetup.presfig( 10, figsize=[15,8])
        ax = plot_rate_points( zmax=zmax, ageaxis=True, showerr=True, logscale=False)
        eta, etaerr  = etaBest['t-1fp']
        fp, fperr  = fpBest['t-1fp']
        snr02 = snrate( z, dtdmodel='t-1fp', sfrmodel='HB06', fitpar= [eta*1.6,0.2] )
        snr05 = snrate( z, dtdmodel='t-1fp', sfrmodel='HB06', fitpar= [eta,fp] )
        snr08 = snrate( z, dtdmodel='t-1fp', sfrmodel='HB06', fitpar= [eta*0.5,0.8] )
        ax.fill_between( z, snr08*0.98, snr08*1.02, color=cp.lightblue, alpha=0.5 )
        ax.fill_between( z, snr05*0.98, snr05*1.02, color=cp.green, alpha=0.5 )
        ax.fill_between( z, snr02*0.98, snr02*1.02, color=cp.red, alpha=0.5 )
        ax.set_ylim(0,1.5)

        ax.text(2., 1.4, 'Prompt\nFraction', ha='center',va='top', color='k', fontsize='x-large')
        ax.text(2., 1.1, '0.8', ha='center',va='center', color=cp.darkbluegray, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.76, '0.5', ha='center',va='center', color=cp.darkgreen, backgroundcolor='w', fontsize='x-large')
        ax.text(2., 0.43, '0.2', ha='center',va='center', color=cp.darkred, backgroundcolor='w', fontsize='x-large')

        pl.draw()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/snr_dtd_fprompt_HB06.png' ) )


    for fig in [23,24,25,26,27,28,29]:
        # scaled fits to the DD, SD and SD+DD models
        if fig not in figlist : continue
        plotsetup.presfig( fig, figsize=[15,8])
        ax = plot_rate_points( zmax=zmax, ageaxis=True, showerr='stat', logscale=False)

        for group, color, ref in zip( grouplist, colorlist, reflist ) :
            if fig in [23,26] : bpsmod = 'BPS.DD.' + group
            elif fig in [24,27] : bpsmod = 'BPS.SD.' + group
            elif fig in [25,28,29] : bpsmod = 'BPS.SD+DD.' + group
            if fig < 26 : eta, etaerr = 1,0
            else : eta,etaerr = etaBest[ bpsmod ]

            if etaerr :
                snrbpsa = snrate( z, dtdmodel=bpsmod, fitpar=eta+etaerr  )
                snrbpsb = snrate( z, dtdmodel=bpsmod, fitpar=eta-etaerr  )
                ax.fill_between( z, snrbpsb, snrbpsa, color=color, alpha=0.3 )
            else :
                snrbps = snrate( z, dtdmodel=bpsmod, sfrmodel='HB06', fitpar=1  )
                ax.plot( z, snrbps, color=color, ls='-', marker=' ', label=ref)

            modfam = bpsmod.split('.')[1]
            if fig!=9 :
                ax.text(2., 1., modfam + ' Models', ha='center',color='k', fontsize='x-large')
        ax.set_ylim(0,1.5)
        pl.draw()
        if fig in [23,24,25] : pngfile = 'snr_dtd_bps_%s_HB06.png'%bpsmod.split('.')[1].lower()
        if fig in [26,27,28] : pngfile = 'snr_dtd_bps_%s_HB06_fit.png'%bpsmod.split('.')[1].lower()
        if fig in [29] : pngfile = 'snr_dtd_bps_%s_fit_HB06_nolabel.png'%bpsmod.split('.')[1].lower()
        pl.savefig(os.path.expanduser('~/Dropbox/Presentations/1409.Chicago/' + pngfile ) )




def wfirsthighz():
    """  how many SN Ia could we find at z<2.5 from the WFIRST SN survey ?
    :return:
    """
    from pytools import cosmo

    areaD = 5.81 # deg2
    areaM = 8.96 # deg2
    areaW = 27.44 # deg2
    areaSky = 41253. # deg2 all sky

    tsurvey = 2 # years (total survey time in the rest frame)

    zbins = np.array( [ 1.6, 1.8, 2.0, 2.2, 2.4 ] )

    # survey volume in each redshift bin  (Mpc^3) :
    volDz = (areaD/areaSky) * np.array( [ cosmo.volume( z+0.1 ) - cosmo.volume(z-0.1) for z in zbins ] )
    volMz = (areaM/areaSky) * np.array( [ cosmo.volume( z+0.1 ) - cosmo.volume(z-0.1) for z in zbins ] )
    volWz = (areaW/areaSky) * np.array( [ cosmo.volume( z+0.1 ) - cosmo.volume(z-0.1) for z in zbins ] )

    # SN rate in each redshift bin
    snrz = np.array( [ 0.92, 0.9, 0.8, 0.7, 0.6 ] ) *1e-4 # SN yr-1 Mpc^3

    # total number of explosions within the survey volume at each z
    NexpDz =  snrz * tsurvey * volDz / (1+zbins)
    NexpMz =  snrz * tsurvey * volMz / (1+zbins)
    NexpWz =  snrz * tsurvey * volWz / (1+zbins)

    # approximate detection efficiencies in each z bin
    deteffDz = np.array( [1.0, 1.0, 1.0, 0.9, 0.8] )
    deteffMz = np.array( [1.0, 1.0, 1.0, 0.8, 0.7] )
    deteffWz = np.array( [1.0, 1.0, 1.0, 0.6, 0.3] )

    # Number detected in each z bin
    NdetDz = NexpDz * deteffDz
    NdetMz = NexpMz * deteffMz
    NdetWz = NexpWz * deteffWz

    NdetTot = NdetDz + NdetMz + NdetWz
    print( NdetTot )




def addaxislabels( zmax=2.55, logscale=False ):
    fig = pl.gcf()
    fig.subplots_adjust( left=0.12, bottom=0.12, top=0.86, right=0.95 )
    ax = pl.gca()
    ax.set_xlabel('Redshift')
    if logscale :
        ax.set_ylabel('SNR ( yr$^{-1}$ Mpc$^{-3}$ h$_{70}}^{3}$ )')
    else :
        ax.set_ylabel('SNR ( 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^{3}$ )')
    ax.set_xlim([-0.05, zmax] )
    #ax.set_ylim([0, 2.98])
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.5 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.25 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.5 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )


def addageaxis():
    ax = pl.gca()
    axtop = ax.twiny()
    axtop.set_xlim( ax.get_xlim() )
    ageticks = np.array( [13,8,5,3] )
    zageticks = zfromt( ageticks )
    axtop.set_xticks( zageticks )
    axtop.set_xticklabels( ageticks )
    axtop.set_xlabel('Age of Universe (Gyr)')
    pl.draw()
    return( ax )



def dtd_bps( dtdmodel='DD', bpsgroup='Yungelson', datstyle='edges' ) :
    """ returns an interpolator that gives the DTD(t) value at any time t [Gyr]
    within a hubble time.
    :param dtdmodel:  Provide a BPS model as e.g. :  'bps.dd.yungelson'
    :param bpsgroup:  The BPS modeling group name (.dat file column name)
    :param datstyle:  Use the 'edges' data file or the 'mid' data file
    :return:
    """
    from astropy.io import ascii
    from scipy.interpolate import interp1d

    bpsdatfile = os.path.join( _thisdir, 'bps_dtd_%s.%s.dat'%( dtdmodel.lower(), datstyle) )
    bpsdat = ascii.read( bpsdatfile,  format='commented_header',
                         header_start=-1, data_start=0 )

    #   following Reddy+ 2009 we convert from Kroupa to Salpeter
    # by multiplying a factor 1.7, then we multiply by 0.7 to get
    # from Salpeter to diet-Salpeter :
    snr = 1.7 * 0.7 * bpsdat[bpsgroup]  # specific SN rate [ yr-1 Msun-1 ]

    # NOTE : in G14 did we use a factor 0.7 to convert from Kroupa to diet-Salpeter ?
    # snr = 0.7 * bpsdat[bpsgroup]  # specific SN rate [ yr-1 Msun-1 ]

    # The time array in the .dat files from Gijs gives the right-hand side
    # of each time bin (i.e. the long-delay edge of each bin).
    # Note the factor of 10^-3 which converts the time  variable
    # from Myr (in the .dat files) into Gyr (as expected by the
    # convolution integral in snrate()
    tGyr = bpsdat['t'] * 1e-3  # time [ Gyr ]

    dtdinterp = interp1d( tGyr, snr, bounds_error=False, fill_value=0 )

    return( dtdinterp )


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


def snrate( z, dtdmodel,  sfrmodel='B13', fitpar=1, tstep=0.01, zFirst = 20, H0=70, Om=0.3, bpsdatstyle='edges' ):
    """  Returns the predicted SN Ia rate at the given redshift for the
    given DTD model in units of SNuVol =  10-4 yr-1 Mpc-3

    :param z:  redshift array
    :param dtdmodel:  Provide a BPS model as e.g. :  'bps.dd.yungelson'
    :param fitpar : SN creation efficiency factor, to scale the given DTD
    :param tstep:  time step in Gyr
    :param zFirst:  redshift of first star formation
    :param H0: hubble constant
    :param Om: Omega_matter
    :return:  np.NDArray with the value of the SNR at each redshift in z.
    """
    import numpy as np
    from scipy.interpolate import interp1d
    import csfr
    Ode=1-Om


    # Create a time array, in units of Gyr, sampling the age of the
    # universe from the epoch of first star formation to the present
    #  * age of univ. at z of First Stars
    tFS = cosmo.agez(zFirst, H0=H0, Om=Om, Ode=Ode, unit='Gyr')
    tFS = 0
    #  * age of universe today
    tNow = cosmo.agez( 0, H0=H0, Om=Om, Ode=Ode, unit='Gyr')

    # Array giving time since big bang
    tBB = np.arange( tFS+tstep/2., tNow, tstep )

    # Fill the arrays for the Star formation history and the DTD
    # Note that we evaluate the SFR at times in the tBB array,
    # relative to t=0 at the big bang. We evaluate the DTD at
    # times tBB-tFS, i.e. relative to tFS at the epoch of first stars.
    zvals_from_tarray = np.array( [ cosmo.zfromt( tt, H0=H0, Om=Om) for tt in tBB ]  )

    # define the star formation rate density [ Msun yr-1 Mpc-3 ]
    if sfrmodel.startswith('B13'):
        if sfrmodel.endswith('high') :
            SFRt = csfr.sfrB13( zvals_from_tarray, dust='high')
        elif sfrmodel.endswith('low') :
            SFRt = csfr.sfrB13( zvals_from_tarray, dust='low')
        else :
            SFRt = csfr.sfrB13( zvals_from_tarray, dust='mid')
        getSFRz = csfr.sfrB13
    elif sfrmodel.startswith('HB06') :
        SFRt = csfr.sfrHB06( zvals_from_tarray )
    elif sfrmodel=='tophat' :
        tophat = lambda z : np.where( (1.<z) & (z<1.2), 1., 0. )
        SFRt = tophat( zvals_from_tarray)
        getSFRz = tophat


    # NOTE : the DTD function should return values of the SN rate in units of ( yr-1 Msun-1 )
    #   as a function of time t in (Gyr).

    if np.iterable(fitpar) :
        eta = fitpar[0]
        if dtdmodel.lower() == 't-1fp' :
            # t-1 delay time distribution for t>0.5 Gyr, flat below.
            fp = fitpar[1]
    else :
        eta = fitpar

    if dtdmodel.lower().startswith( 'bps' ) :
        # Note the factor of 1e4, since these BPS DTD models are in units
        #  of (yr-1 Msun-1) and we are working in units of (10-4 yr-1 Msun-1)
        dtdmodparts = dtdmodel.split('.')
        dtdfunc = dtd_bps( dtdmodparts[1], dtdmodparts[2], datstyle=bpsdatstyle  )
    elif dtdmodel.lower() == 't-1' :
        # t-1 delay time distribution.
        dtdfunc = lambda tGyr : np.where( (0.04<tGyr) & (tGyr<13.7), 1./tGyr, 0.0 )
    elif dtdmodel.lower() == 't-1fp' :
        # t-1 delay time distribution for t>0.5 Gyr, flat below.
        t0=0.5
        tmin=0.04
        tmax=13.5
        K = ( fp/(1-fp) ) * np.log(tmax/t0) / (t0-tmin)
        dtdfunc  = lambda tGyr : np.where( tGyr>t0,  1. / tGyr, np.where( tGyr>tmin, K, 0) )

    elif dtdmodel.lower().startswith('delta') :
        # deltaGyr = float( dtdmodel.lstrip('delta') )
        dtdfunc = lambda tGyr : np.where( np.abs(tGyr-2.)<tstep, 1., 0.0 )


    DTDt = eta * dtdfunc( tBB - tFS )  # [ yr-1 Msun-1 ]

    # Do a discrete convolution of the SFR and DTD.
    # Note that we multiply the time step (in Gyr) by 1e9, because the convolution
    # integral is done in yr (since both SFR and DTD carry units of yr-1, not Gyr-1)
    # and we also multiply the output rate by 1e4 to get units of ( 10-4 yr-1 Mpc-3 )
    SNRt = 1e4 * tstep * 1e9 * np.convolve( SFRt, DTDt, mode='full')[:len(tBB)]

    # Set up an interpolator function to evaluate the SNR at any
    # time t, relative to the big bang
    tz = cosmo.agez( z, H0=H0, Om=Om, Ode=Ode, unit='Gyr')
    getSNRt = interp1d( tBB, SNRt, kind='linear', bounds_error=False, fill_value=0 )

    return( getSNRt(tz) )



def fitSNR( data='ALL', dtdmodel='t-1', systematicshift='mid', sfrmodel='B13',
            tstep=0.01, Ngrid=0, skipN2M=True, showplots=False, sneakpeek=False ):
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
    import numdifftools as ndt
    from ratetable import ALL, RATELISTALL
    start = time.time()

    ratelist = RATELISTALL
    composite = ALL

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
    from scipy import optimize as scopt
    # sig2 = np.max([robsStatP,robsStatM],axis=0)**2
    # Ndof = len(robs) - 3

    # initial guesses for the DTD model parameters
    # Note that we use scaling factors so that the values
    # handled by the scipy minimizer are close to unity.
    eta0, etascale = 2.,  1.
    fprompt0, fpscale = 0.5, 1.

    if dtdmodel.startswith('t-1'):
        etascale=1e-13


    # Set the CSFR(z) curve and the observed points based on the
    # systematic error test case defined by the dust and csfr parameters
    sfrmodel = sfrmodel + systematicshift
    if systematicshift=='high' : robsfit = robs+robsSysP
    elif systematicshift=='mid' : robsfit = robs
    elif systematicshift=='low' : robsfit = robs-robsSysM

    def neglikefunc( fitpar ) :
        """ negative likelihood function to minimize.

        model='fprompt' :
        fitpar[0] = eta : the efficiency factor for Ia production
        fitpar[1] = fprompt : the prompt Ia fraction

        model='t-1'
        fitpar = eta

        model= (any BPS DTD model)
        fitpar = eta

        """

        eta = fitpar[0] * etascale
        if len(fitpar) == 2 :
            fprompt = fitpar[1] * fpscale
            fitparscaled = [eta,fprompt]
        else :
            fitparscaled = [eta]

        # define the SNR model using a given prompt fraction
        # evaluated at each z bin center

        # TODO : add in the fprompt component
        rmod = snrate( zobs, dtdmodel,  sfrmodel=sfrmodel, fitpar=fitparscaled, tstep=tstep, zFirst = 20, H0=70, Om=0.3 )

        # TODO : integrate the SNR model over the redshift bin to match the observed bin
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

    # maximize the likelihood (i.e. minimize the negative log likelihood)
    if dtdmodel.endswith('fp'):
        output = scopt.fmin_l_bfgs_b( neglikefunc, [eta0,fprompt0],
                                      fprime=None, args=(), approx_grad=True,
                                      bounds=[ (0,None), (0,0.99) ] )
    else :
        output = scopt.fmin_l_bfgs_b( neglikefunc, [eta0],
                                      fprime=None, args=(), approx_grad=True,
                                      bounds=[ (0,None) ] )
    fitparBest, likemax, infodict = output
    hessian = ndt.Hessian( fun=neglikefunc )
    invhess = np.linalg.inv( hessian( fitparBest ) )

    # estimate the uncertainties from the sqrt of the diagonal elements
    # of the inverse hessian matrix and then
    # rescale eta and fprompt to their natural units
    fitparErr = np.sqrt( np.diagonal( invhess ) )
    if dtdmodel.endswith('fp'):
        etaErr, fpromptErr = fitparErr
        fpromptBest = fitparBest[1]
        fpromptBest *= fpscale
        fpromptErr *= fpscale
    else :
        etaErr = fitparErr
    etaBest = fitparBest[0]
    etaBest *= etascale
    etaErr *= etascale

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
    label='%s data, %s dtd, %s dust, %s CSFR(z)'%(data, dtdmodel, systematicshift, sfrmodel)
    print( label )
    print( "  likemax = %.3e "%(likemax ) )
    print( "  eta = %.5e  +-%.5e"%(etaBest, etaErr ) )
    if dtdmodel.endswith('fp') :
        print( "  fprompt = %.3f  +-%.3f\n"%(fpromptBest, fpromptErr ) )

    if showplots :
        from matplotlib import pyplot as pl
        plot_rate_points( slide=1, showerr='stat', eta=etaBest, etaerr=etaErr )

        pl.draw()

    end = time.time()
    print( "total execution time = %.1f seconds"%(end-start) )
    if Ngrid : return( likegrid, N2Mgrid, etagrid, fPgrid )
    elif dtdmodel.endswith('fp') : return( etaBest, etaErr, fpromptBest, fpromptErr )
    else : return( etaBest, etaErr )



def dtdtestfig( dtdmodel='BPS.SD.Delta', datstyle='mid'):
    from pytools import cosmo
    import csfr

    z = np.arange( 0.01, 2.2, 0.01 )
    snr,sfr = snrate( z, dtdmodel=dtdmodel, fitpar=1, sfrmodel='B13', tstep=0.01, bpsdatstyle=datstyle)

    # tophat = lambda z : np.where( (1.<z) & (z<1.2), 1., 0. )
    # sfr = tophat( z )
    # sfr2  = csfr.sfrB13(z)

    fig1 = pl.figure(1)

    tz = cosmo.agez( z, H0=70, Om=0.3, Ode=0.7 )
    tnow = cosmo.agez( 0, H0=70, Om=0.3, Ode=0.7 )
    tlookback = tnow - tz
    #pl.plot( tlookback, snrWangSD, 'b--'  )
    #pl.plot( tlookback, snrYungSD, 'k-'  )
    # pl.plot( tz, sfr2/sfr2.max(), 'k--'  )
    dtdfunc = dtd_bps('SD','Yungelson',datstyle=datstyle)
    dtdt = dtdfunc(tlookback)

    ax1 = pl.subplot( 1 , 2, 1)
    pl.plot( z, sfr/sfr.max(), 'b-'  )
    pl.plot( z, snr/snr.max(), 'r:'  )
    ax1.set_xlabel('redshift')
    ax1.set_ylabel( 'SN or SF rate')

    ax2 = pl.subplot( 1 , 2,  2)
    pl.plot( tlookback, sfr/sfr.max(), 'b-'  )
    pl.plot( tlookback, snr/snr.max(), 'r:'  )
    pl.plot( tlookback, dtdt/dtdt.max(), 'k--'  )
    ax2.set_xlabel('lookback time (Gyr)')
    ax2.set_ylabel( 'DTD or SN or SF rate')

    ax2.axvline( tlookback[snr.argmax()], color='r')
    ax2.axvline( tlookback[sfr.argmax()], color='b')
    print( tlookback[snr.argmax()], tlookback[sfr.argmax()], tlookback[snr.argmax()] - tlookback[sfr.argmax()])

    ax1.set_ylim( -0.1, 1.2 )
    ax2.set_ylim( -0.1, 1.2 )
    pl.draw()

