"""
2012.09.10 S.Rodney
Some figures showing rates vs z for Garching Pres.
"""
from ratetable import *
from matplotlib import pyplot as p
import numpy as np
from hstsnpipe.tools.snana import extinction 

def groundrates( outfile='groundRates.pdf', label=True, **kwarg):
    import rcpar
    rcpar.usetex()
    rcpar.fullpaperfig(2,figsize=[8,4])
    rcParams['lines.linewidth']=1.5
    rcParams['font.size']=12
    rcParams['axes.labelsize']=14
    p.ioff()
    ax = p.gca()

    D10.marker='<'; D10.mfc='w'; D10.mec='0.1'; 
    G11.marker='o'; G11.mfc='w'; G11.mec='0.1'; 
    RT10.marker='v'; RT10.mfc='w'; RT10.mec='0.1'; 
    P12.marker='s'; P12.mfc='0.6'; P12.mec='0.4'; 

    ratelist = [RT10,P12,D10,G11]  

    for i in range(len(ratelist)):
        R = ratelist[i]
        R.plot( **kwarg )
    p.xlabel("Redshift")
    p.ylabel("SNR $\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")

    from matplotlib.font_manager import FontProperties
    fs=FontProperties(size='9')

    ax = p.gca()
    xlim = [-0.05, 2.61]
    ylim = [-0.05, 2.01]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #ax.set_yticklabels(['0.1','1','10'])
    p.subplots_adjust(0.11, 0.125, 0.98, 0.98)

    if label : 
        labels = [ 
            (D10,(-0.25,-0.4), {'lw':0,'headwidth':0,'width':0,'color':D10.mec}),
            (P12,(-0.4,-0.55), {'lw':0,'headwidth':0,'width':0,'color':P12.mec}),
            (G11,(0.1,0.15), {'lw':0,'headwidth':0,'width':0,'color':G11.mec}),
            (RT10,(-0.6,0.35), {'lw':0,'headwidth':0,'width':0,'color':RT10.mec}),
            ]

        for rt,xyoff,aprop in labels : 
            p.annotate(rt.ref, 
                     xy=[rt.z[-1],rt.rate[-1]], 
                     xytext=[rt.z[-1]+xyoff[0],rt.rate[-1]+xyoff[1]], 
                     xycoords='data', textcoords='data', 
                     arrowprops=aprop, color=rt.mec)

    p.savefig(outfile)
    p.draw()
    p.ion()
    return(outfile)

def goodsRates(outfile='goodsRates.pdf'):
    p.clf()
    groundrates( label=False )
    D08.marker='^'; D08.mfc='0.6'; D08.mec='r'; 
    D08.plot()

    aprop=  {'lw':0,'shrink':0.05,'headwidth':0,'width':0,'color':D08.mec}
    p.annotate(D08.ref, 
             xy=[D08.z[-1],D08.rate[-1]], 
             xytext=[D08.z[-1]+0.2, D08.rate[-1]-0.2],
             xycoords='data', textcoords='data', 
             arrowprops=aprop, color=D08.mec )
    p.savefig(outfile)
    p.draw()
    p.ion()
    return(outfile)


def ratepoints( dataset='R13fig', rgb=True, outfile=None,
                ageaxis=True, **kwarg):
    """ 
    simple SN rate plot. 

    R13fig : make the figure for the Rodney et al 2013 rates paper
    """
    import candelsRates
    from cosmo import agez, zfromt

    ratelist = RATELISTALL

    if dataset=='R13fig' :
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.fullpaperfig(figsize=[8,4])
        ax = p.gca()
        for i in range(len(RATELISTLOWZ)):
            R = RATELISTLOWZ[i]
            R.mfc='0.8'
            R.mec='0.8'
            R.color='0.8'
            R.mew=0
            R.marker='o'
            R.zerrplus = np.zeros( len( R.zerrplus ) )
            R.zerrminus = np.zeros( len( R.zerrplus ) )
            R.plot( thicksys=False, zorder=-10*i, **kwarg )
        for i in range(len(RATELISTHIGHZ)-1):
            R = RATELISTHIGHZ[i]
            R.ms=12
            R.zerrplus = np.zeros( len( R.zerrplus ) )
            R.zerrminus = np.zeros( len( R.zerrplus ) )
            p.plot( R.z, R.rate, marker=R.marker, ms=12, mew=0, ls=' ',
                    color='w', alpha=1, zorder=10*i )
            R.plot( thicksys=False, zorder=10*(i+1), alpha=0.3, **kwarg )

        R13.mfc='darkorange'
        R13.mec='darkred'
        R13.color='darkred'
        R13.ms=15
        R13.plot( thicksys=True, zorder=1000 )
        ax.text( 2.33, 0.72, 'CANDELS', ha='left', va='top', color='darkred' )
        ax.text( 1.49, 0.21, 'CLASH', ha='right', va='top', color='darkmagenta' )
        ax.text( 1.27, 1.31, 'GOODS', ha='left', va='top', color='darkgreen' )
        ax.text( 1.13, 1.44, 'CSS', ha='right', va='top', color='darkcyan' )
        ax.text( 1.67, 1.15, 'SDF', ha='right', va='top', color='0.2' )

    elif dataset =='R13figB' : 
        ratelist = RATELISTALL + [ ALL ] # B12, D08, G13, R13,
        for i in range(len(ratelist)-2):
            R = ratelist[i]
            R.mfc='0.8'
            R.mec='0.8'
            R.color='0.8'
            R.mew=0
            R.ms=8
            R.zerrplus = np.zeros( len( R.zerrplus ) )
            R.zerrminus = np.zeros( len( R.zerrplus ) )
            R.plot( thicksys=False, zorder=-i, **kwarg )
        R13.mfc='0.7'
        R13.mec='0.7'
        R13.color='0.7'
        R13.ms=15
        R13.plot( thicksys=False , zorder=1000 )
        ALL.plot( thicksys=True , zorder=10000 )

    elif dataset =='GROUND' : 
        ratelist = RATELISTGROUND + [GROUND]
        for i in range(len(ratelist)-1):
            R = ratelist[i]
            R.mfc='0.8'
            R.mec='0.8'
            R.color='0.8'
            R.mew=0
            R.ms=8
            R.zerrplus = np.zeros( len( R.zerrplus ) )
            R.zerrminus = np.zeros( len( R.zerrplus ) )
            R.plot( thicksys=False, zorder=-i, **kwarg )
        R13.mfc='0.7'
        R13.mec='0.7'
        R13.color='0.7'
        R13.ms=15
        R13.plot( thicksys=False , zorder=1000 )
        GROUND.plot( thicksys=True , zorder=10000 )

    elif dataset =='ALLHST' : 
        ALLHST.plot( thicksys=True, **kwarg )

    elif dataset =='HST' : 
        ratelist = [ B12, D08, G13, R13 ]
        for R in ratelist :
            R.plot( thicksys=True, **kwarg )

    elif dataset =='COMPOSITE' : 
        ratelist = [ ALL ]
        ratelist[-1].mfc='darkorange'
        ratelist[-1].mec='darkred'
        ratelist[-1].color='darkred'
        for i in range(len(ratelist)-1):
            R = ratelist[i]
            R.plot( thicksys=False, zorder=-10*i, **kwarg )
        ratelist[-1].plot( thicksys=True, zorder=1000 )
    elif dataset =='ALLPRECC' : 
        ratelist = RATELISTPRECC + [ALLPRECC]
        for r in ratelist[:-1] :
            r.mfc='0.8'
            r.mec='0.8'
            r.color='0.8'
            r.mew=0
        ratelist[-1].mfc='darkcyan'
        ratelist[-1].mec='darkgreen'
        ratelist[-1].color='darkgreen'
        for i in range(len(ratelist)-1):
            R = ratelist[i]
            R.zerrplus = np.zeros( len( R.zerrplus ) )
            R.zerrminus = np.zeros( len( R.zerrplus ) )
            R.plot( thicksys=False, zorder=-10*i, **kwarg )
        ratelist[-1].plot( thicksys=True, zorder=1000 )


    p.xlabel("Redshift")
    p.ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")

    from matplotlib.font_manager import FontProperties
    fs=FontProperties(size='9')

    ax = p.gca()
    ax.set_xlim([-0.05, 2.81])
    ax.set_ylim([0.0, 2.55])
    p.subplots_adjust(0.11, 0.125, 0.98, 0.87)

    # if label==1 : 
    #     ax.text(0.55,0.15,'Ground-based Surveys',ha='left',va='top',fontsize='small',color='0.3')
    #     ax.text(1.78,2.15,'CANDELS',ha='left',va='center',fontsize='medium',color='darkred')
    #     ax.text(1.65,0.25,'GOODS+PANS',ha='left',va='top',fontsize='small',color='m')
    #     ax.text(1.52,1.65,'SCP',ha='right',va='top',fontsize='small',color='darkcyan')
    #     ax.text(1.67,1.35,'SDF',ha='right',va='top',fontsize='small',color='darkgreen')
    #     ax.text( 2.0, 1.8, 'Stat. + Syst. Error',ha='left',va='center', color='darkred')
    #     ax.text( 2.2, 1.1, 'Systematic Only',ha='left',va='center', color='darkred')
    #     ax.plot( [1.77, 1.97], [ 1.45, 1.76], marker=' ', ls='-', color='darkred', lw=0.8 )
    #     ax.plot( [1.8, 2.15], [ 0.85, 1.1], marker=' ', ls='-', color='darkred', lw=0.8 )

    if dataset=='R13' : 
        # ax.text(0.05,1.75,'Projected rates from the full CANDELS sample',ha='left',va='center',fontsize='medium',color='darkgreen')
        ax.text(2.18,1.75,'SN Ia rates from the\n full CANDELS sample',ha='left',va='top',fontsize='medium',color='darkgreen')
        ax.text(2.22,1.2,'pre-CANDELS avg',ha='left',va='top',fontsize='medium',color='0.4')
        ax.plot( [1.75, 2.17], [ 0.65,1.65], marker=' ', ls='-', color='darkgreen', lw=0.8 )
        ax.plot( [1.8, 2.2], [ 0.65,1.12], marker=' ', ls='-', color='k', lw=0.8 )

        ax.text( 0.55, 1.55, 'All Error',ha='right',va='center', color='k')
        ax.text( 0.45, 1.3, 'Systematic\nError Only',ha='right',va='top', color='darkgreen')
        ax.plot( [0.57, 0.73], [ 1.5, 0.95], marker=' ', ls='-', color='k', lw=0.8 )
        ax.plot( [0.47, 0.73], [ 1.2, 0.7], marker=' ', ls='-', color='darkgreen', lw=0.8 )

        pcterr = ( np.sqrt( R13.errstatplus**2 + R13.errstatminus**2 ) / R13.rate ) * 100
        ax.text( 0.25, 2.3, '%d'%pcterr[0]+' \%',ha='center',va='top', color='0.3')
        ax.text( 0.75, 2.3, '%d'%pcterr[1]+' \%',ha='center',va='top', color='0.3')
        ax.text( 1.25, 2.3, '%d'%pcterr[2]+' \%',ha='center',va='top', color='0.3')
        ax.text( 1.75, 2.3, '%d'%pcterr[3]+' \%',ha='center',va='top', color='0.3')
        ax.text( 2.15, 2.3, '%d'%pcterr[4]+' \%',ha='center',va='top', color='0.3')
        ax.text( 2.3,  2.3, '(stat. err as pct.)',ha='left',va='top', color='0.3')

        zbinedges=np.append( R13.z-R13.zerrminus, R13.z[-1]+R13.zerrplus[-1])
        Ndet, dNdetPlus, dNdetMinus = candelsRates.countDetections( zbinedges=zbinedges )

        ax.text( 0.25, 2.1, r'%.1f$_{-%.1f}^{+%.1f}$'%(round(Ndet[0],1), round(dNdetPlus[0],1), round(dNdetMinus[0],1)),ha='center',va='top', color='0.3', backgroundcolor='w')
        ax.text( 0.75, 2.1, r'%.1f$_{-%.1f}^{+%.1f}$'%(round(Ndet[1],1), round(dNdetPlus[1],1), round(dNdetMinus[1],1)),ha='center',va='top', color='0.3', backgroundcolor='w')
        ax.text( 1.25, 2.1, r'%.1f$_{-%.1f}^{+%.1f}$'%(round(Ndet[2],1), round(dNdetPlus[2],1), round(dNdetMinus[2],1)),ha='center',va='top', color='0.3', backgroundcolor='w')
        ax.text( 1.75, 2.1, r'%.1f$_{-%.1f}^{+%.1f}$'%(round(Ndet[3],1), round(dNdetPlus[3],1), round(dNdetMinus[3],1)),ha='center',va='top', color='0.3', backgroundcolor='w')
        ax.text( 2.15, 2.1, r'%.1f$_{-%.1f}^{+%.1f}$'%(round(Ndet[4],1), round(dNdetPlus[4],1), round(dNdetMinus[4],1)),ha='center',va='top', color='0.3', backgroundcolor='w')
        ax.text( 2.3, 2.1, '(\# detected)',ha='left',va='top', color='0.3')

        ax.set_ylim( 0.,2.6 )

    if dataset=='HST' : 
        ax.text(0.5,0.95,'SN Ia rates from a decade of HST surveys',transform=ax.transAxes,ha='center',va='top',fontsize='medium',color='k')

        ax.text( 0.38, 2.15, 'SCP (ACS)',ha='right',va='center', color='blue')
        ax.text( 0.38, 1.85, 'GOODS (ACS)',ha='right',va='top', color='darkgreen')
        ax.text( 1.8, 2.15, 'CLASH (WFC3+ACS)',ha='left',va='top', color='red')
        ax.text( 1.9, 1.95, 'CANDELS (WFC3)',ha='left',va='top', color='m')

        ax.plot( [0.4, 0.78], [ 2.12, 1.21], marker=' ', ls='-', color='b', lw=0.8 )
        ax.plot( [0.18, 0.45], [ 1.72, 0.9], marker=' ', ls='-', color='darkgreen', lw=0.8 )

        ax.plot( [1.85, 1.52], [ 1.98, 0.5], marker=' ', ls='-', color='r', lw=0.8 )
        ax.plot( [2.05, 1.77], [ 1.83, 0.72], marker=' ', ls='-', color='m', lw=0.8 )

        ax.set_ylim( 0.,2.6 )

    #if dataset=='ALL' : 
    # ax.text(0.5,0.95,'Composite SN Ia rates from all surveys',transform=ax.transAxes,ha='center',va='top',fontsize='medium',color='darkred')

    if ageaxis :              
        axtop = ax.twiny()
        axtop.set_xlim( ax.get_xlim() )
        ageticks = np.array( [13,8,5,3] )
        zageticks = zfromt( ageticks )
        axtop.set_xticks( zageticks )
        axtop.set_xticklabels( ageticks )
        axtop.set_xlabel('Age of Universe [Gyr]')

    p.draw()
    #p.ion()
    if outfile: 
        p.savefig(outfile)
        return(outfile)
    else:
        return( ax )


def dtd3( SFRmodel='behroozi', zFirst=20, 
          zRange=[0.01,2.8], logz=False, normalize=True, 
          label=True, **kwargs ):
    """ plot 3 DTD models using the given SFR(z) model. 
    "t-1" : the t^-1 power-law 
    "gauss" : a 1-component gaussian, following Dahlen+ 2008 
    "2comp" : the Scannapieco+Bildsten 2-component model
    """
    import rates
    p.clf()
    ax = p.gca()
    ratepoints(showR13=True, label=False)

    dtdmodel(model='t-1', parameters=[0.007],
             SFRmodel=SFRmodel, zFirst=zFirst,
             zRange=zRange, gammaZ=-0.15,
             logz=logz, 
             color='b', **kwargs)

    dtdmodel(model='2comp', parameters=[0.044,2.0],
             SFRmodel=SFRmodel, zFirst=zFirst,
             zRange=zRange, gammaZ=-0.15,
             logz=logz, 
             color='g', **kwargs )

    dtdmodel(model='gauss', parameters=[0.0009,1.5], 
             SFRmodel=SFRmodel, zFirst=zFirst,
             zRange=zRange, gammaZ=-0.15,
             logz=logz,  
             color='r', **kwargs )

    ax = p.gca()
    ax.set_ylim( [0.007, 2.01] )
    ax.set_xlim( [0.0, 2.81] )
    #setp( ax.get_xticklabels(), visible=False)
    ax.set_ylabel(r"$10^{-4}$ yr$^{-1}$ Mpc$^{-3}$ ${h_{70}}^3$")
    if logz:
        ax.set_xlabel("log( 1 + z )")
    else : 
        ax.set_xlabel("Redshift")

    if label :
        ax.text( 2.06, 0.67, r'DTD $\propto$ t$^{-1}$', color='b', ha='left',va='bottom')
        ax.text( 2.05, 1.26, r'SNR $\propto$ A M$_*$ + B \.M$_*$', color='g', ha='left',va='bottom')
        ax.text( 2.42, 0.1,  r'gaussian, $\tau$=1.5 Gyr', color='r', ha='right',va='bottom')


def dtd2( SFRmodel='behroozi', zFirst=20, 
          zRange=[0.01,2.8], logz=False, 
          showR13=True, normalize=True, 
          label=True, **kwargs ):
    """ plot 2 DTD models using the given SFR(z) model. 
    "t-1" : the t^-1 power-law (i.e. the DD scenario)
    "SD"  : a cocktail of SD models from Wang:2009 (WD+MS,RG,He)
    """
    import rates
    import rcpar
    rcpar.usetex()
    rcpar.presfig(figsize=[10,5])
    p.ioff()
    p.clf()
    ax = p.gca()
    ratepoints(showR13=showR13, label=False)

    ratelist = [P12,D10,G11,D08]  
    if showR13 : ratelist.append( R13 )

    z = np.arange( zRange[0], zRange[1], 0.1 )
    dtdmodel(model='WD+WD', eta=1,
             SFRmodel=SFRmodel, zFirst=zFirst,
             z=z, logz=logz, scale2fit=ratelist, 
             color='b', lw=2, **kwargs)

    dtdmodel(model='SD', eta=1,
             SFRmodel=SFRmodel, zFirst=zFirst,
             z=z, logz=logz, scale2fit=ratelist, 
             color='g', lw=2, ls='--', **kwargs )

    ax = p.gca()
    #ax.set_ylim( [0.007, 2.01] )
    #ax.set_xlim( [0.0, 2.81] )
    #setp( ax.get_xticklabels(), visible=False)
    ax.set_ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")
    if logz:
        ax.set_xlabel("log( 1 + z )")
    else : 
        ax.set_xlabel("Redshift")

    if label :
        if SFRmodel=='behroozi':
            ax.text( 2.25, 0.3, r'DD (t$^{-1})$', color='b', ha='left',va='top')
            ax.text( 2.13, 0.8, r'SD (Wang 2009)', color='g', ha='left',va='bottom')
        elif SFRmodel.startswith('cole'):
            ax.text( 2.25, 0.5, r'DD (t$^{-1})$', color='b', ha='left',va='top')
            ax.text( 2.16, 1.13, r'SD (Wang 2009)', color='g', ha='left',va='bottom')


def powerlawplot_tonset( SFRmodel='behroozi', zFirst=20, 
                         zRange=[0.01,2.8], logz=False, 
                         showR13=True, normalize=True, 
                         label=False, **kwargs ):
    """ plot t-1 power law models with different values for the 
    onset of the first SNIa progenitors
    """
    import rates
    import rcpar
    rcpar.usetex()
    rcpar.presfig(figsize=[10,5])
    p.ioff()
    p.clf()
    ax = p.gca()
    ratepoints(showR13=showR13, label=False)

    ratelist = [P12,D10,G11,D08,R13]  
    #if showR13 : ratelist.append( R13 )

    z = np.arange( zRange[0], zRange[1], 0.1 )
    zlist, snrlist = [],[]
    for tau,c,ls in zip([0.001,0.05,0.5],['b','g','r'],['-','--','-.']) : 
        ztau,snrtau = dtdmodel(model='t-1', eta=1, tau=tau, 
                             SFRmodel=SFRmodel, zFirst=zFirst,
                             z=z, logz=logz, scale2fit= 'normalize', #ratelist, 
                             lw=2, color=c, ls=ls, label='tF=%.2f'%tau, 
                             **kwargs)
        zlist.append( ztau ) 
        snrlist.append( snrtau )

    ax = p.gca()
    ax.text( 2.7,1.7,r'DTD$\propto t^{-1}$ for $t>t_0$',ha='right',va='top')
    ax.text( 2.3,1.1,r'$t_0=1 Myr$',ha='left',va='top',color='b')
    ax.text( 2.4,0.7,r'$50 Myr$',ha='left',va='top',color='g')
    ax.text( 2.2,0.4,r'$500 Myr$',ha='left',va='top',color='r')
    #ax.set_ylim( [0.007, 2.01] )
    #ax.set_xlim( [0.0, 2.81] )
    #setp( ax.get_xticklabels(), visible=False)
    ax.set_ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")
    if logz:
        ax.set_xlabel("log( 1 + z )")
    else : 
        ax.set_xlabel("Redshift")

    if label :
        ax.legend(loc='upper right',handlelength=1.5, frameon=False, handletextpad=0.5 )

def snrFprompt( dataset='R13figB', SFRmodel='behroozi', zFirst=20, showR15=False,
                zrange=[0.01,2.8], logz=False, 
                normalize=False, label=True, savetable=False, 
                paperfig=True, **kwargs ):
    """ plot t-1 power law models with different values for the 
    fraction of SNIa that are prompt (t<500 Myr).
    """
    from cosmo import agez, zfromt
    import dtdconvolve
    import candelsRates

    if paperfig :
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.fullpaperfig(figsize=[8,4])
    else :
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.presfig(figsize=[8,5])

    ax = p.gca()
    z = np.arange( zrange[0], zrange[1], 0.1 )
    snrlist = []
    tau = 0.5
    for fprompt,c,ls in zip([0.8,0.45,0.1],['b','g','r'],['-','--','-.']) : 
        if fprompt == 0 : 
            snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=normalize, 
                                      t0=0.01, f0=fprompt, **kwargs)
        else : 
            snr = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=normalize, 
                                      f0=fprompt, **kwargs)
        snrlist.append( snr )
        if logz: z = np.log10( 1+z )
        line = p.plot( z, snr, color=c, ls=ls, zorder=-10 )

    ax = ratepoints( dataset=dataset )

    # ax.text( 2.7,1.7,r'DTD$\propto t^{-1}$\\ for $t>500 Myr$',ha='right',va='top')
    if label : 
        ax.text( 2.72,1.4,r'\noindent Fraction of SN Ia\\that are "prompt"\\ ($\tau<$500 Myr)',ha='right',va='bottom', backgroundcolor='w')
        ax.text( 2.4,1.04,r'$f_{P}=0.8$',ha='left',va='top',color='b',backgroundcolor='w',fontsize=14)
        ax.text( 2.45,0.59,r'0.45',ha='left',va='top',color='g',backgroundcolor='w',fontsize=14)
        ax.text( 2.4,0.32,r'0.1',ha='left',va='top',color='r',backgroundcolor='w',fontsize=14)
    ax.text(0.4,1.25,'Weighted\nAverage',ha='right',va='bottom',color='saddlebrown')
    ax.plot([0.42,0.55],[1.22,0.65],ls='-',marker=' ',color='saddlebrown',lw=0.5)

    if savetable : 
        import sfrate
        t = agez( z )
        sfr = sfrate.behroozi2012( z, fitto='B12' )
        fout=open(savetable,'w')
        print >> fout, "%7s  %7s  %8s  %8s  %8s  %8s"%('redshift','age(Gyr)','SFR','SNR(0.8)','SNR(0.5)','SNR(0.1)')
        for zz,tt,sf,s08,s05,s01 in zip( z, t, sfr, snrlist[0], snrlist[1], snrlist[2] ):
            print >> fout, "%7.4f  %7.3f  %8.5e  %8.5e  %8.5e  %8.5e"%(
                zz, tt, sf, s08, s05, s01 )
        fout.close()

    ax.set_xlim([-0.05, zrange[-1]-0.05])
    p.draw()
    return( ax )

    junk="""
    ax.set_xlim([-0.05, 2.81])
       
    axtop = ax.twiny()
    axtop.set_xlim( ax.get_xlim() )
    ageticks = np.array( [13,8,5,3] )
    zageticks = zfromt( ageticks )
    axtop.set_xticks( zageticks )
    axtop.set_xticklabels( ageticks )
    axtop.set_xlabel('Age of Universe [Gyr]')
    p.subplots_adjust( top=0.87 )

    ax.set_ylabel(r"SNR ~$\left[ 10^{{-4}}~ {\rm yr}^{-1}~ {\rm Mpc}^{{-3}}~ {{\rm h}_{70}}^{3}\right]$")
    if logz:
        ax.set_xlabel("log( 1 + z )")
    else : 
        ax.set_xlabel("Redshift")
    ax.set_ylim([-0.15, 2.05])
"""



def dtdFprompt(  ):
    """ 
    Plot theoretical t-1 DTD curves with varying fraction of prompt SNe
    """
    from cosmo import agez, zfromt
    import dtdconvolve 
    import rcpar 
    rcpar.presfig(1)
    rcpar.usetex()
    p.clf()
    ax = p.gca()

    t = 10**( np.arange( 6, 10.5, 0.05) ) / 1e9  # time in Gyr
    dtd10 = dtdconvolve.dtd_t1_fprompt( t, eta=1e-13, f0=0.1, t0=0.5, tmin=0.01, tmax=13.5  )
    dtd50 = dtdconvolve.dtd_t1_fprompt( t, eta=1e-13, f0=0.5, t0=0.5, tmin=0.01, tmax=13.5  )
    dtd80 = dtdconvolve.dtd_t1_fprompt( t, eta=1e-13, f0=0.8, t0=0.5, tmin=0.01, tmax=13.5  )
    #p.loglog( t, dtd10, lw=2, color='r', ls='-.', label='0.1')
    #p.loglog( t, dtd50, lw=2, color='g', ls='--', label='0.5')
    #p.loglog( t, dtd80, lw=2, color='b', ls='-', label='0.8')

    p.semilogy( t, dtd80, lw=2, color='b', ls='-', label='0.8')
    p.semilogy( t, dtd50, lw=1.8, color='g', ls='--', label='0.5')
    p.semilogy( t, dtd10, lw=1.5, color='r', ls='-.', label='0.1')

    ax.set_ylim( 2e-14,5e-12)
    ax.set_xlim( -0.1, 3.1 )

    ax.set_ylabel( r'SNIa ~~yr$^{-1}$ Msun$^{-1}$')
    ax.set_xlabel( r'Delay Time [ Gyr ]')

    ax.text( 0.51, 2.5e-12, 'f$_{prompt}$=0.8', color='b', ha='left',va='center')
    ax.text( 0.52, 7e-13, '0.5', color='g', ha='left',va='center')
    ax.text( 0.51, 7e-14, '0.1', color='r', ha='left',va='center')
    p.draw()



def futureHSTrates( zrange=[0.01,2.8], **kwargs):
    """ plot the improving precision of the HST rates
    Now = first half of candels
    End = CANDELS+CLASH+GOODS+PANS
   
    """
    import dtdconvolve
    from ratetable import R13
    ax = p.gca()
    
    #pcterrNow = np.array( [ [100,370], [58,103],[67,136], [72,158] ]  )
    pcterrNow = np.array( [ R13.errstatminus / R13.rate, 
                            R13.errstatplus / R13.rate ] ) * 100.
    pcterrEnd = np.array( [ [ 39, 21, 24, 41 ], 
                            [ 59, 27, 30, 59 ] ] )

    z = np.arange( zrange[0], zrange[1], 0.1 )
    zlist, snrlist = [],[]
    tau = 0.5
    snrlist = []

    snr80 = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=True, f0=0.8, **kwargs)
    snr50 = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=True, f0=0.5, **kwargs)
    snr10 = dtdconvolve.snrate( z, DTDmodel='t-1+p', normalize=True, f0=0.1, **kwargs)
    
    snr80 = (snr80 - snr50)/snr50 * 100
    snr10 = (snr10 - snr50)/snr50 * 100

    ax.plot( z, snr80, color='b', ls='-' , lw=2)
    ax.axhline( 0, color='g', ls='--' , lw=2)
    ax.plot( z, snr10, color='r', ls='-.', lw=2)

    ax.errorbar( [0.25,0.75,1.25,1.75], [0,0,0,0], pcterrNow, marker='D', ls=' ',  capsize=0, elinewidth=2, ecolor='k', color='k' ) 
    ax.errorbar( [0.3,0.8,1.3,1.8], [0,0,0,0], pcterrEnd, marker='D', ls=' ',  capsize=0, elinewidth=2, ecolor='m', color='m' ) 

    ax.set_ylabel('Percent Error')
    ax.set_xlabel('Redshift')

    ax.text( 1.95, 60, 'All HST Field Surveys', color='m' )
    ax.text( 1.8, 100, 'First half of CANDELS', color='k' )
    ax.set_xlim( -0.1, 3.1 )

    p.draw()



def dtdmodel(model='WD+WD', eta=1, tau=1,
             SFRmodel='behroozi', zFirst=11, 
             z=np.arange( 0.01, 2.8, 0.1 ), gammaZ=-0.15,
             logz=False, scale2fit=[P12,D10,G11,D08], 
             debug=False, verbose=True, **kwarg):
    """ 
    NOTE : THIS IS RATHER BROKEN RIGHT NOW

    Using the given DTD parameters, convolve the DTD model 
    with the SFR(z) curve named by SFRmodel to get SNR(z).
    Then plot it.
    normalize : normalize the model so that it agrees with 
    the observed rate at z=0.55 from Perrett+ 2012.
    """
    if debug: import pdb; pdb.set_trace()
    import dtdconvolve

    scale = 1
    if scale2fit in ['norm','normalize','Normalize']:
        snr = dtdconvolve.snrate( z, model=model, eta=eta, tau=tau, 
                                  SFR=SFRmodel, zFirst=zFirst, normalize=True)

    elif len(scale2fit) : 
        from scipy import optimize as scopt
        zobs, snrobs, snrerrplus, snrerrminus = [], [], [], []
        for rate in scale2fit : 
            zobs += rate.z.tolist()
            snrobs += rate.rate.tolist()
            snrerrplus += rate.errplus.tolist()
            snrerrminus += rate.errminus.tolist()

        isort = np.argsort(zobs)
        zobs = np.array( zobs )[ isort ]
        snrobs = np.array( snrobs )[ isort ]
        snrerrplus = np.array( snrerrplus )[ isort ]
        snrerrminus = np.array( snrerrminus )[ isort ]

        # evaluate the model SNR at the observed redshifts,
        # normalized to the Perrett:2012 z=0.55 rate measurement
        # (so that we start off pretty close)
        snrmodel = dtdconvolve.snrate( zobs, model=model, eta=eta, tau=tau,
                                       SFR=SFRmodel, zFirst=zFirst, normalize=False )
        # SNR is reported in NIa yr-1 Mpc-3 h70^3
        # change units to  NIa 10^-4 yr-1 Mpc-3 h70^3
        snrmodel *= 1e4

        def chi2( scalefactor ):
            snrfit = snrmodel * scalefactor 
            snrerr = np.where( snrfit>=snrobs, snrerrplus, snrerrminus )
            return np.sum((snrobs - snrfit)**2 / snrerr**2)

        scale,chi2min,niter,ncalls,warnflag  = scopt.fmin( 
            chi2, 1.0, disp=0, full_output=True,retall=False )
        Ndof = len(zobs)-1
        if verbose: 
            print( 'Chi2 fit results for %s:'%model )
            print( '  scale = %.3e:'%scale )
            print( '  chi2 = %.3f:'%chi2min )
            print( '  Ndof = %i:'%Ndof )
            print( '  chi2_red = %.3f'%(chi2min / Ndof) )
            print( '  Niter = %i'%niter )
            print( '  Ncalls = %i'%ncalls )
            print( '  warning = %i'%warnflag )
    else : 
        snr = dtdconvolve.snrate( z, model=model, eta=eta, tau=tau, 
                                  SFR=SFRmodel, zFirst=zFirst, normalize=False)

    # SNR is reported in NIa yr-1 Mpc-3 h70^3
    # change units to  NIa 10^-4 yr-1 Mpc-3 h70^3
    snr *= 1e4

    # multiply by the scale factor 
    # (this is optimized if scale2fit list was provided, 
    #  otherwise its unity)
    snr *= scale

    if logz: z = np.log10( 1+z)
    line = p.plot( z, snr, **kwarg)
    return( z, snr ) 


def dtdchecks( ):
    """ plots to verify that the primary DTD models 
    match the published figures of reference.  
    Compare Figure 1 to Ruiter et al 2012, Fig 7
    Compare Figure 2 to Wang et al 2010, Fig 2
    """
    import dtdconvolve
    
    t = 10**( np.arange( 7, 10.5, 0.05) ) / 1e9  # time in Gyr

    p.figure(1)
    p.clf()
    phiDD = dtdconvolve.delaytimedist( t, 'DD' )
    p.loglog( t*1000, phiDD*1e14 )
    ax = p.gca()
    ax.set_xlim( 1e2, 1.2e4 )
    ax.set_ylim( 1e-1, 1.05e2 )
    ax.text( 0.95, 0.95, "WD+WD model Delay Time Dist'n \nCompare to Ruiter et al 2012, Fig 7.", ha='right',va='top', transform=ax.transAxes )
    ax.set_xlabel( r'Time [Myr]')
    ax.set_ylabel( r'SNIa yr$^{-1}$ Msun$^{-1}$')
    p.grid()

    p.figure(2)
    p.clf()
    phiHe = dtdconvolve.delaytimedist( t, 'WD+He' )
    phiMS = dtdconvolve.delaytimedist( t, 'WD+MS' )
    phiRG = dtdconvolve.delaytimedist( t, 'WD+RG' )
    p.semilogy( np.log10(t*1e9), phiHe*1e14 )
    p.semilogy( np.log10(t*1e9), phiMS*1e14 )
    p.semilogy( np.log10(t*1e9), phiRG*1e14 )
    ax = p.gca()
    ax.set_ylim(1e-5,2.2)
    ax.set_xlim( 7.38, 10.7 )
    ax.text( 0.95, 0.95, "SD model Delay Time Dist'ns \nCompare to to Bo Wang et al 2010, Fig 2.", ha='right',va='top', transform=ax.transAxes )
    ax.set_xlabel( r'Time [log(t/1e9 yr)]')
    ax.set_ylabel( r'SNIa yr$^{-1}$ Msun$^{-1}$')
    p.grid()
    p.draw()


def dtdplot( wang=True, ruiter=True, normalize=False, showredshift=False ):
    """ Plot theoretical DTD curves, with t=0 on the left
    and time increasing to the right.  
    normalize :  
        True = re-scale so that the DTD curves integrate to unity
        False = use the true theoretical efficiency factors to 
            plot in actual units of SN yr-1 Msun-1
    showredshift :  add a redshift axis along the top, with a
        t=0 line fixed to z=2.31, the peak of the cosmic SFR curve
    """
    from cosmo import agez, zfromt
    import dtdconvolve 
    import rcpar 
    rcpar.presfig(1)
    rcpar.usetex()
    p.clf()

    t = 10**( np.arange( 7, 10.5, 0.05) ) / 1e9  # time in Gyr
    #logt = np.log10(t*1e9) # log of time in Gyr

    #z0 = 2.31  # z for peak of SFR(z) curve
    #z = zfromt( agez( z0 ) + t ) # time elapsed since peak of SFR(z)
    #logz = np.log10( 1+z )

    ax = p.gca()

    # Bo Wang SD Models
    wdhe = 1e14*dtdconvolve.delaytimedist( t, 'WD+He' )
    wdms = 1e14*dtdconvolve.delaytimedist( t, 'WD+MS' )
    wdrg = 1e14*dtdconvolve.delaytimedist( t, 'WD+RG' )
    wdhe = np.where( wdhe>1e-5, wdhe, 1e-5 )
    wdms = np.where( wdms>1e-5, wdms, 1e-5 )
    wdrg = np.where( wdrg>1e-5, wdrg, 1e-5 )

    lw = p.rcParams['lines.linewidth']+1

    p.loglog( t, wdhe, lw=lw, color='green', ls=':', label='WD+He')
    p.loglog( t, wdms, lw=lw, color='green', ls='--', label='WD+MS')
    p.loglog( t, wdrg, lw=lw, color='green', ls='-', label='WD+RG')

    # Ruiter DD model
    wdwd = 1e14*dtdconvolve.delaytimedist( t, 'DD' )
    wdwd = np.where( wdwd>1e-5, wdwd, 1e-5 )
    p.loglog( t, wdwd, lw=lw, color='b', ls='-', label='WD+WD')

    # ax.invert_xaxis()
    #print ax.get_xlim()
    ax.set_ylim(1e-5,100)
    ax.set_xlim( 15, -0.1  )

    ax.set_ylabel( r'10$^{-14}$ ~~SNIa ~~yr$^{-1}$ Msun$^{-1}$')
    ax.set_xlabel( r'Time [ Gyr ]')
    #ax.set_ylim([1.1e-5,4e-2])
        #ax.set_xticks( [8, 9, 10] )
    #ax.set_xlim([zfromt( agez(9.0) , ] )  1.1e-5,4e-2])
    p.draw()



def ccrates( logplot=False, labelobs=False ):
    """ plot the CC rates measurements 
    overplot a few SFR(z) lines 
    Scaled for the CCSN generation efficiency : 
      SNR(z) = k * h^2 * SFR(z)
      k = 0.007 Msun-1  (Dahlen+ 2012)
    """
    from ratetable import CCRATELIST
    import sfrate
    
    for ccr in CCRATELIST : ccr.plot()
    ax = p.gca()
    ax.set_xlabel('Redshift')
    ax.set_ylabel('CC SN Rate [ 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]' )
    
    k = 0.0070 
    z = np.arange( 0.01, 2.5, 0.05 )
    sfrB12 = sfrate.behroozi2012( z, fitto='B12' )
    sfrHB06 = sfrate.behroozi2012( z, fitto='HB06' )
    sfrY08 = sfrate.yuksel2008( z, fitto='H11' )

    r1= ax.plot( z, 0.9*k*sfrY08*1e4, 'b:', label='Horiuchi+ 2011')
    r2= ax.plot( z, 1.5*k*sfrHB06*1e4, 'g--', label='Hopkins \& Beacom 2006')
    r3= ax.plot( z, 1.6*k*sfrB12*1e4, 'r-', label='Behroozi+ 2012')
    if logplot: ax.semilogy()
    if labelobs : ax.legend(loc='lower right')
    else : ax.legend([r1,r2,r3], ['Horiuchi+ 2011','Hopkins \& Beacom 2006','Behroozi+ 2012'], loc='lower right')

    ax.grid()
    ax.set_ylim(0.2, 18)
    ax.set_xlim(-0.1, 2.6)
    p.draw()


def iarates( logplot=False, labelobs=False ):
    """ plot the SNIa rates measurements 
    overplot a few SNR(z) lines 
    """
    from ratetable import RATELISTBEST
    import sfrate
    from dtdconvolve import snrate
    
    for tnr in RATELISTBEST : 
        tnr.plot()
    ax = p.gca()
    ax.set_xlabel('Redshift')
    ax.set_ylabel('SNIa Rate [ 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]' )
    
    #z = arange( 0.01, 2.6, 0.05 )
    z = np.arange( 0.01, 2.1, 0.1 )

    snrPrompt = snrate( z, model='SD', eta=[8,0.1,0.05], SFR='behroozi', normalize=True )
    #snrMixed = snrate( z, model='WD+WD', eta=1, SFR='behroozi', normalize=True )
    snrMixed = snrate( z, model='SD', eta=[5,1,2], SFR='behroozi', normalize=True )
    #snrDelayed = snrate( z, model='gauss', eta=1, tau=2, SFR='behroozi', normalize=True )
    snrDelayed = snrate( z, model='lognormal', eta=1, tau=0.1, SFR='behroozi', normalize=True )

    r1= ax.plot( z, snrPrompt*1e4, 'b:', label='Prompt')
    r2= ax.plot( z, snrMixed*1e4, 'g--', label='Mixed')
    r3= ax.plot( z, snrDelayed*1e4, 'r-', label='Delayed')
    if logplot : ax.semilogy()
    if labelobs: ax.legend(loc='lower right')
    else : ax.legend([r1,r2,r3], ['Prompt','Mixed','Delayed'], loc='lower right')
    ax.grid()
    ax.set_ylim(0., 2.1)
    ax.set_xlim(-0.1, 2.6)
    p.draw()



def plotExtinctionDistributions( paperfig=False ):
    """ plot the distributions of SALT2 C and host Av values:
    Type Ia on top, CC on bottom
    Dashed = high dust
    Solid = mid dust
    Dash-dot = low dust 

    paperfig : set to True to generate a plot suitable for inclusion
     as a paper figure (thicker lines, sized to fit in a text column)
    """
    from matplotlib import pylab as p
    import numpy as np
    if paperfig : 
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.halfpaperfig()

    meanIa, meanCC = extinction.meanAv()
    
    # Background : whole-frame ylabel
    fig = p.gcf()
    axbg = fig.add_subplot( 1,1,1, frameon=False)
    axbg.set_xticks([])
    axbg.set_yticks([])
    axbg.set_ylabel('Relative Probability', labelpad=20)

    # Top : the Ia extinction distributions
    ax1 = fig.add_subplot( 2,1,1 )
    Cin, Cobs = extinction.convolveCAv( 'highIa' );
    Avin = Cin + 0.1
    ax1.plot( Avin, extinction.highIa_c( Cin ),'r--' , 
              label='high: Gauss($\sigma$=0.62)')

    Cin, Cobs = extinction.convolveCAv( 'midIa');
    Avin = Cin + 0.1
    ax1.plot( Avin, extinction.midIa_c( Cin ),'g-', 
              label='mid: Exp.($\\tau$=0.33)' )

    Cin, Cobs = extinction.convolveCAv( 'lowIa' );
    Avin = Cin + 0.1
    ax1.plot( Avin, extinction.lowIa_c( Cin ),'b-.', 
              label='low: Gauss($\sigma$=0.15)+Exp.($\\tau$=0.15)')

    ax1.text( meanIa[0]+0.1, 0.7, '$\langle A_{V} \\rangle$=%.1f'%meanIa[0],ha='center',va='center',color='r',backgroundcolor='w')
    ax1.text( meanIa[1], 0.5, '%.1f'%meanIa[1],ha='center',va='center',color='g',backgroundcolor='w')
    ax1.text( meanIa[2]+0.1, 0.2, '%.1f'%meanIa[2],ha='center',va='center',color='b',backgroundcolor='w')
    legend1 = ax1.legend( loc='upper right', frameon=False, title='SN Ia Dust Priors : p(A$_{V}|$Ia)',
                          handlelength=1.8, handletextpad=0.1 )
    p.setp( legend1.get_title(), fontsize='medium' )
    p.setp( legend1.get_texts(), fontsize='small' )
    ax1top = ax1.twiny()
    ax1top.set_xlim( ax1.get_xlim()[0]-0.1, ax1.get_xlim()[1]-0.1 )
    ax1top.set_xlabel('SALT2 $\mathscript{C}')

    # bottom : the CC extinction distributions
    ax2 = fig.add_subplot( 2,1,2, sharex=ax1, sharey=ax1 )
    Avin = np.arange( -0.1, 3.0, 0.01 )
    ax2.plot( Avin, extinction.highCC( Avin )/extinction.highCC( Avin ).max(),
              'r--', label='high: $\sigma$=0.8, $\\tau$=2.8, $R_0$=3')
    ax2.plot( Avin, extinction.midCC( Avin )/extinction.midCC( Avin ).max(),
              'g-', label='mid: $\sigma$=0.6, $\\tau$=1.7, $R_0$=4')
    ax2.plot( Avin, extinction.lowCC( Avin )/extinction.lowCC( Avin ).max(),
              'b-.', label='low: $\sigma$=0.15, $\\tau$=0.5, $R_0$=1' )
    ax2.set_xlabel('host extinction (A$_{V}$)')


    legend2 = ax2.legend( loc='upper right', frameon=False, title='CC SN Dust Priors : p(A$_{V}|$CC)')
    p.setp( legend2.get_title(), fontsize='medium' )
    p.setp( legend2.get_texts(), fontsize='small' )

    ax2.text( 0.72, 0.8, '$\langle A_{V} \\rangle$=%.1f'%meanCC[0],ha='center',va='center',color='r',backgroundcolor='w')
    ax2.text( 0.67, 0.5, '%.1f'%meanCC[1],ha='center',va='center',color='g',backgroundcolor='w')
    ax2.text( 0.4, 0.25, '%.1f'%meanCC[2],ha='center',va='center',color='b',backgroundcolor='w')

    ax1.set_xlim( -0.3, 2.8 )
    ax1.set_ylim( -0.02, 1.08 )
    fig.subplots_adjust( wspace=0, hspace=0, left=0.1, right=0.97, bottom=0.1, top=0.9)






#------------------------------------------------------------
#  old functions follow: may not work 
#------------------------------------------------------------


def likemap( ratelist=RATELISTALL,
             Amin=0, Amax=1, nA = 21, Ascale=1e-7, 
             Bmin=0, Bmax=10, nB = 21, Bscale=1e-4, 
             logA = False, 
             showimage=False, showcontour=False, 
             verbose=False, debug=False):
    """ 
    construct a 2-D map of likelihoods
    using the 2-component model.  The parameter
    grids are defined as [min, max, Nstep]
    and are scaled by 1e-4 for fitting
    against the data. 
    """
    from numpy import array, zeros, linspace
    if debug: import pdb; pdb.set_trace()

    like = zeros( [nA, nB] )
    Agrid, Astep = linspace( Amin, Amax, nA, retstep=True )    
    Bgrid, Bstep = linspace( Bmin, Bmax, nB, retstep=True )

    for i in xrange( nA ) :
        A = Agrid[i] * Ascale
        if logA : A = 10**(Agrid[i] * Ascale)
        for j in xrange( nB ):
            B = Bgrid[j] * Bscale
            
            chi2 = fitmodel( ratelist=ratelist, parameters=[A, B], 
                             model='2comp', SFRmodel='colemid', 
                             showfit=False, showpoints=False, 
                             dofit = False, debug=False, 
                             verbose=verbose )
            like[i,j] = exp( -0.5 * chi2 )

    from pylab import imshow, contour, clf, gca
    if showimage: 
        imshow( like, interpolation='nearest',aspect='auto', 
                cmap=cm.gray_r, origin='lower', 
                extent=[Bmin,Bmax,Amin,Amax] )
    if showcontour: 
        # draw 68 and 95% contours
        # Set up an array 'Labove' such that  the cell value in Labove[i,j] 
        # is equal to the sum of likelihood values from all cells that have a 
        # likelihood value higher than like[i,j]  
        like1d = like.ravel() / like.sum()
        Labove = array( [ (like1d[where( like1d > like1d[i]) ]).sum() for i in range(len(like1d)) ] )

        # Now unravel the Labove array back into 2-dimensions,
        Labove = Labove.reshape( like.shape ) 
        # find the right and left edges of the 95% contour:
        B95max = where(Labove<=0.68)[0].max()  * Bstep + Bmin
        B95min = where(Labove<=0.68)[0].min()  * Bstep + Bmin
        # find the top and bottom edges of the 95% contour:
        A95max = where(Labove<=0.68)[1].max()  * Astep + Amin
        A95min = where(Labove<=0.68)[1].min()  * Astep + Amin

        Apk = int(like.argmax() / nB) * Astep + Amin
        Bpk = (like.argmax() % nB) * Bstep + Bmin

        import colors
        contour( Labove, [0.68, 0.95], aspect='auto', origin='lower',  
                 colors=[colors.orange,colors.blue], 
                 extent=[Bmin,Bmax,Amin,Amax])
        print("A = %.2f +%.2f -%.2f"%(Apk, A95max-Apk, Apk-A95min) )
        print("B = %.2f +%.2f -%.2f"%(Bpk, B95max-Bpk, Bpk-B95min) )
        
        ax = gca()
        ax.set_xlabel("B")
        ax.set_ylabel("A")

    return( like )


def fitmodel( ratelist=RATELISTALL, 
              model='2comp', parameters=[0.1e-4,3.6e-4], 
              SFRmodel='colemid', gammaZ=-0.15, 
              showfit=False, showpoints=False, 
              dofit = True, 
              debug=False, verbose=True,
              **kwarg) : 
    """ 
    find the min chi2 fit to the data.
    dofit=False : don't minimize, just return chi2.
    """
    from scipy import optimize as scopt
    from iceberg import snrate 
    if debug: import pdb; pdb.set_trace()

    snr, dsnr, z  = [], [], []
    for R in ratelist : 
        snr +=  R.rate.tolist() 
        dsnr += R.meanerr.tolist() 
        z += R.z.tolist() 
    snr = array(snr)
    dsnr = array(dsnr) 
    z = array(z)
    
    if model=='2comp' : 
        chi2 = lambda param : ( ( snr - 1e4*snrate( z, model=model, parameters=param, 
                                                    SFR=SFRmodel, normalize=False))**2 
                                / dsnr**2 ).sum()
    #elif model=='SD(Z)' : 
    #    chi2 = lambda param : ( ( 
    #            snr - 1e4*param[0]*snrate( z, model=model, 
    #                                       parameters=[1,param[1]], 
    #                                       SFR=SFRmodel, 
    #                                       normalize=False))**2 
    #                            / dsnr**2 ).sum()

    else : 
        modelrate = snrate( z, model=model, parameters=[1], 
                            gammaZ=gammaZ, SFR=SFRmodel, 
                            normalize=False) 
        chi2 = lambda param : ( ( snr - 1e4*modelrate*param[0])**2 
                                / dsnr**2 ).sum() 
        
    if dofit: fitparam = scopt.fmin( chi2, parameters, disp=0 )
    else : fitparam = parameters

    chi2min = chi2(fitparam)
    if verbose: 
        print( "---------------------")
        print( "SNR Model : %s"% model ) 
        print( "SFR Model : %s"% SFRmodel ) 
        print( "Min chi2 : %.3f"% chi2min ) 
        for i in range(len(fitparam) ) : 
            print( "param %i : %.3e"%(i,fitparam[i]) )
 
    if showpoints : 
        plotall(ratelist, showlegend=False, usecolors=True)
    if showfit : 
        z,snr = plotmodel( model=model, parameters=fitparam,
                          SFRmodel=SFRmodel, gammaZ=gammaZ,
                          **kwarg)
        #return( z,snr ) 
    return( chi2min ) 




def plotall(ratelist, showlegend=True, usecolors=True, 
            **kwarg):
    """ plot all the TNSN rate measurements in ratelist"""
    import pylab
    import rcpar
    rcpar.usetex()
    pylab.ioff()
    #ax = pylab.axes()
    ax = pylab.gca()
    ax.semilogy()
        
    for i in range(len(ratelist)):
        R = ratelist[i]
        R.plot( **kwarg )
    pylab.xlabel("redshift")
    pylab.ylabel(r"SNR ~~~$\left[ ~10^{-4}~ {\rm yr}^{-1}~ {\rm Mpc}^{-3}~ {{\rm h}_{70}}^{3}~\right]$")

    from matplotlib.font_manager import FontProperties
    fs=FontProperties(size='10')
    if showlegend : 
        pylab.legend(loc="upper left",numpoints=1,handlelen=0,
                     pad=0.05, axespad=0.02, prop=fs)
    #pylab.title("Type Ia SN Rates")
    ax = pylab.gca()
    xlim = [-0.35, 1.9]
    ylim = [ 0.088,  4.1]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_autoscale_on(False)
    
    #pylab.ion()
    #pylab.draw()
    pylab.show()

    #from pylab import twiny, gca
    #ax2 = twiny(ax)
    ##ax2.set_yscale('log')
    #ax2.set_xlim(xlim)
    #ax2.set_ylim(ylim)
    ##ax2.xaxis.set_ticks_position('bottom')
    ##ax2.xaxis.set_ticks_position('both')
    #yax = ax2.yaxis
    #yax.set_ticks_position('right') # move tick labels to right
    ##yax.set_label_position('right') # move yaxis label to right
    #yax.set_ticks_position('both')  # restore left-side tick marks
    ##label = yax.get_label()
    ##setp(label, rotation=-90)
       


def ratePriors( logplot=False, labelobs=False, paperfig=True ):
    """ plotting simple linear + exponential SN rate models
    for use as SN rate priors  """
    from ratetable import ALL, CCRATELIST
    import sfrate
    from models import linexp as snr  # our simple linear+exp model
    import numpy as np
    from matplotlib import pyplot as p
    from priors import mkpriors
    if paperfig : 
        from hstsnpipe.tools.figs import plotsetup
        plotsetup.halfpaperfig(1, figsize=[4.5,6])
    from cosmo import zfromt

    z, iar,ccr,iafrac = mkpriors()
    iar1, iar2, iar3 = iar
    ccr1,ccr2,ccr3 = ccr
    iafracMin, iafracBase, iafracMax = iafrac

    fig = p.gcf()
    ax12bg = fig.add_subplot( 2,1,1, frameon=False )
    ax12bg.set_xticks([])
    ax12bg.set_yticks([])
    ax12bg.set_ylabel('SN Rate [SNuVol]', labelpad=20)
    ax1 = fig.add_subplot( 3,1,1 )
    ax2 = fig.add_subplot( 3,1,2, sharex=ax1 )
    ax3 = fig.add_subplot( 3,1,3, sharex=ax1 ) 

    p.axes( ax1 )
    ALL.plot()

    r1= ax1.plot( z, iar1, 'r-', label='Prompt')
    r2= ax1.plot( z, iar2, 'g-', label='Mixed')
    r3= ax1.plot( z, iar3, 'b-', label='Delayed')
    if logplot : ax1.semilogy()
    #if labelobs: ax1.legend(loc='upper right')
    #else : ax1.legend([r1,r2,r3], ['Prompt','Mixed','Delayed'], loc='upper right',frameon=False, handlelength=1.5, handletextpad=0.5)
    ax1.set_ylim(0.01, 1.6)
    ax1.set_xlim(-0.1, 2.6)
    ax1.text( 0.1, 1.4, 'SNe Ia', ha='left', va='top', fontsize='x-large' )
    
    # ------------------------------------------------------------
    # CC SN Rates follow the cosmic SFR(z) 
    p.axes( ax2 )
    for ccr in CCRATELIST : ccr.plot()
    ax2.set_xlabel('Redshift')
    #ax2.set_ylabel('CC SN Rate [ 10$^{-4}$ yr$^{-1}$ Mpc$^{-3}$ h$_{70}^3$ ]' )
    # ax2.set_ylabel('CC SN Rate [ SNuVol ]' )
    
    k = 0.0070 
    sfrB12 = sfrate.behroozi2012( z, fitto='B12' )
    ccr0 = 1.6*k*sfrB12*1e4

    r1= ax2.plot( z, ccr1, 'b-', label='High', zorder=30)
    r2= ax2.plot( z, ccr2, 'g-', label='Mid', zorder=30)
    r3= ax2.plot( z, ccr3, 'r-', label='Low', zorder=30 )
    r0= ax2.plot( z, ccr0, 'm--', label='Behroozi+ 2012', lw=p.rcParams['lines.linewidth']*2, zorder=20)
    if logplot: ax2.semilogy()
    #if labelobs : ax2.legend(loc='upper right')
    #else : ax2.legend([r1,r2,r3,r0], ['High','Mid','Low','$\propto$SFR(z)'], loc='lower right', frameon=False, handlelength=1.5, handletextpad=0.5)

    #ax2.grid()
    ax2.set_ylim(0.2, 17.8)
    ax2.set_xlim(-0.1, 2.6)
    ax2.text( 0.1, 16, 'CC SNe', ha='left', va='top', fontsize='x-large' )

    # ------------------------------------------------------------
    # Plot the fraction of SNe that are Type Ia as a function of z
    p.axes( ax3 )

    ax3.plot( z, iafracMax, color='r' )
    ax3.plot( z, iafracBase, color='g' )
    ax3.plot( z, iafracMin, color='b' )
    #ax3.grid()
    ax3.set_xlabel('Redshift')
    ax3.text( 2.5, 0.32, 'Fraction of all SNe that are SN Ia', ha='right', va='top',fontsize='large', color='k')
    ax3.set_xlim(0., 2.6)
    ax3.set_ylim(-0.03, 0.37)
    ax3.text( 2.2, 0.155, 'high', color='r', backgroundcolor='w', ha='left', va='center')
    ax3.text( 2.0, 0.078, 'mid', color='g', backgroundcolor='w' , ha='left', va='center')
    ax3.text( 1.7, 0.023, 'low', color='b', backgroundcolor='w' , ha='left', va='center')

    axtop = ax1.twiny()
    axtop.set_xlim( ax1.get_xlim() )
    ageticks = np.array( [13,8,5,3] )
    zageticks = zfromt( ageticks )
    axtop.set_xticks( zageticks )
    axtop.set_xticklabels( ageticks )
    axtop.set_xlabel('Age of Universe [Gyr]')

    p.setp( ax1.get_xticklabels(), visible=False )
    p.setp( ax2.get_xticklabels(), visible=False )

    fig.subplots_adjust( left=0.12, bottom=0.08, right=0.97, top=0.92, hspace=0, wspace=0)


