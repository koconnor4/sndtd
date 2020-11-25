# 2013.05.24
# S.Rodney
from matplotlib.pyplot import *
from numpy import *
import pprint

from hstsnpipe.tools import snana

from glob import glob
iadatfiles = glob('/Users/rodney/Dropbox/GOODS-ACS/HST_GOODS_*.dat')
ccdatfiles = glob('/Users/rodney/Dropbox/GOODS-ACS/CC/HST_GOODS_*.dat')

goldIalist = ['Yowie','Qiqirn','Anguta',
              'Frodo','Rakke','Spock','Elvis','Manipogo','Mcenroe','Thames',
              'Patuxent','Ramone','Ombo','Ferguson','Eagle','Strolger',
              'Gabi','Greenburg','Nanna','Koekemoer','Lancaster',
              'Aphrodite','Thoth','Borg','Mcguire','Sasquatch',]
silverIalist=['Bilbo','Inanna','Isis','Hawk','Dickinson','Vilas',
              'Athena','Gilgamesh','Zwicky','Kurage','Torngasak',]
bronzeIalist=['Inanna',]
goldCClist = ['denethor',]
silverCClist=['agugux','michabo','hades','artemis','raven','jimi',]
bronzeCClist=['bast','loki','heineken','re','jagger','connors','jefferson','apollo',]

candelsdatfiles = glob('/Users/rodney/Dropbox/CANDELS/CANDELS1_classify/HST_CANDELS1_*.dat')

# RISING SNE FROM CANDELS1 : 
gsadatfiles = [ '/Users/rodney/Dropbox/CANDELS/CANDELS1_classify_v2/HST_CANDELS1_%s.dat'%sn 
                for sn in ['primo','tumbleweed','adams','buchanan','bush','ford',
                           'harrison','lincoln','madison','quincy','roosevelt',
                           'taylor','vanburen','washington','workman','agnew',
                           'rockefeller','jackson'] ]
widedatfiles = [ '/Users/rodney/Dropbox/CANDELS/CANDELS1_classify_v2/HST_CANDELS1_%s.dat'%sn 
                 for sn in ['aidan','carter','clinton','eisenhower','herbert','kennedy',
                            'mikulski','mondale','reagan','truman','wilson','garfield','hughes'] ]

def colorMagClassTest( testset='Ia1', clobber=False, verbose=False ):
    """ compute color mag classifications for all the gold SNIa 
    and the gold/silver/bronze CC SNe. Print a table of results """

    # does increasing the model error improve the color-mag classifications?
    snlist = []
    print('Ia.mod.err= 0.0   0.08   0.24')

    if testset=='Ia1' : datfilelist = iadatfiles[:19]
    elif testset=='Ia2' : datfilelist = iadatfiles[19:]
    elif testset=='cc' : datfilelist = ccdatfiles

    if testset=='c1' : datfilelist = candelsdatfiles[:13]
    elif testset=='c2' : datfilelist = candelsdatfiles[13:26]
    elif testset=='c3' : datfilelist = candelsdatfiles[26:]

    for datfile in datfilelist : 
        sn0 = snana.SuperNova( datfile )
        snana.classify.colorMagClassify( sn0, classfractions='mid', dustmodel='mid', bands='all', Nsim=2000, clobber=clobber, verbose=verbose, modelerror=[0.0,0.,0.] )

        sn1 = snana.SuperNova( datfile )
        snana.classify.colorMagClassify( sn1, classfractions='mid', dustmodel='mid', bands='all', Nsim=2000, clobber=False, verbose=verbose, modelerror=[0.08,0.1,0.1] )

        sn2 = snana.SuperNova( datfile )
        snana.classify.colorMagClassify( sn2, classfractions='mid', dustmodel='mid', bands='all', Nsim=2000, clobber=False, verbose=verbose, modelerror=[0.24,0.3,0.3] )

        print('%5s      %5.2f  %5.2f  %5.2f  '%( sn0.nickname[:3], sn0.PIaColor, sn1.PIaColor, sn2.PIaColor ) )
        
        snlist.append( [sn0,sn1,sn2] )
    return( snlist ) 


def systematicsClassifications( testset='gsa', interactive=False, clobber=False, verbose=False ):
    """ compute LCgrid classifications for all the gold SNIa 
    and the gold/silver/bronze CC SNe. Print a table of results """

    import gc # garbage collector
    #snlist = []

    if testset=='Ia1' : datfilelist = iadatfiles[:19]
    elif testset=='Ia2' : datfilelist = iadatfiles[19:]
    elif testset=='cc' : datfilelist = ccdatfiles

    elif testset=='c1' : datfilelist = candelsdatfiles[:15]
    elif testset=='c2' : datfilelist = candelsdatfiles[15:30]
    elif testset=='c3' : datfilelist = candelsdatfiles[30:]
    elif testset=='c4' : datfilelist = candelsdatfiles[15:]

    # ONLY RISING SNE: 
    elif testset=='gsa1' : datfilelist = gsadatfiles[:11]
    elif testset=='gsa2' : datfilelist = gsadatfiles[11:]
    elif testset=='wide1' : datfilelist = widedatfiles[:7]
    elif testset=='wide2' : datfilelist = widedatfiles[7:]

    else : datfilelist = testset
    
    cprior = {'high':snana.extinction.highIa_c, 
              'mid':snana.extinction.midIa_c, 
              'low':snana.extinction.lowIa_c }
    avprior = {'high':snana.extinction.highCC, 
               'mid':snana.extinction.midCC, 
               'low':snana.extinction.lowCC }

    pmatrices = []
    for datfile in datfilelist : 
        pmatrix = []
        for dustmodel in ['high','mid','low'] :
            plist = []
            for ratemodel in ['high','mid','low'] :
                try: 
                    sn = snana.SuperNova( datfile )

                    print( " LCgrid classification for %s : dust %s  : rate %s "%( sn.nickname.upper(), dustmodel, ratemodel ) )
                    if sn.zerr<0.02 : 
                        nlogz=1
                        ncolorpar=20
                        nlumipar=30
                        npkmjd=40
                    elif sn.zerr<0.11 : 
                        nlogz=5
                        ncolorpar=15
                        nlumipar=20
                        npkmjd=30
                    else : 
                        nlogz=10
                        ncolorpar=10
                        nlumipar=15
                        npkmjd=30
                    sn.doGridClassify( clobber=clobber, verbose=verbose, 
                                       nlogz=nlogz, ncolorpar=ncolorpar, nlumipar=nlumipar, ncolorlaw=1, npkmjd=npkmjd, 
                                       cprior=cprior[dustmodel], avprior=avprior[dustmodel],classfractions=ratemodel, 
                                       useLuminosityPrior=True, trestrange=[-20,60], modelerror=[0.08,0.1,0.1] )
                    figure(3)
                    clf()
                    sn.plotMaxLikeModels( 'flux','mjd', showchi2=False )
                    if interactive: 
                        userin = raw_input("Adjust plot spacing. Hit <return> to save and continue.")
                    savefig('/Users/rodney/Dropbox/CANDELS/LCgridClass/%s.dust%s.rate%s.LCfit.pdf'%(sn.nickname,dustmodel,ratemodel))                
                    figure(4)
                    close()
                    sn.plotClassStatsGrid()
                    savefig('/Users/rodney/Dropbox/CANDELS/LCgridClass/%s.dust%s.rate%s.stats.pdf'%(sn.nickname,dustmodel,ratemodel))

                    print('%5s dust%s rate%s :  %5.2f  '%( sn.nickname[:3], dustmodel, ratemodel, sn.PIa ) )

                    #snlist.append( sn )
                    plist.append( '%.2f'%sn.PIa )
                except : 
                    #snlist.append( sn ) 
                    plist.append( '-1.0')
                finally: 
                    print(sn.nickname.upper())
                    del sn
                    gc.collect()
            pmatrix.append( plist )
            pprint.pprint( array(pmatrix,dtype=float) ) 
        pmatrices.append( array(pmatrix,dtype=float) ) 

    for datfile,pmatrix in zip( datfilelist, pmatrices) : 
        sn = snana.SuperNova( datfile )
        pmatrix = array( pmatrix )
        print( '\n-------------------\n%s : %.2f +%.2f -%.2f\n'%(
                sn.nickname, pmatrix[1][1], pmatrix.max()-pmatrix[1][1], pmatrix[1][1]-pmatrix.min()) )
        pprint.pprint( pmatrix )
        
    return( ) 
    



def doClassify( name, survey='GOODS', interactive=True, ColorMag=True, LCgrid=True, 
             clobber=3, verbose=1 ) : 
    """ interactively classify a GOODS SN and make plots """
    print("Working on %s"%name)

    sn = snana.SuperNova( 'HST_%s_%s.dat'%(survey,name) )
    if clobber: 
        sn.writedatfile( 'HST_%s_%s.dat'%(survey,name), mag2fluxcal=True )
        sn = snana.SuperNova( 'HST_%s_%s.dat'%(survey,name) )

    if interactive : 
        print(" %s: initial LC plot and .dat file adjustment..."%name)
        while True : 
            figure(1)
            clf() 
            sn = snana.SuperNova( 'HST_%s_%s.dat'%(survey,name) )
            sn.plotLightCurve( 'flux','trest', showpkmjdrange=True, autozoom=False )
            userin = raw_input("Adjust .dat file. Type <return> to continue, 'r' to re-plot:  ")
            if userin.strip().lower() != 'r' : break
        close(1)

    if ColorMag : 
        print(" %s: making a color-mag plot..."%name)
        snCM = snana.SuperNova( 'HST_%s_%s.dat'%(survey,name) )
        figure(2)
        clf()
        figure(2)
        snCM.plotColorMagAll( mjd='peak', linelevels=[0.95,0.68,0.0], showpia=True, clobber=clobber, verbose=verbose )
        if interactive: 
            userin = raw_input("Adjust plot spacing. Hit <return> to save and continue.")
        savefig('/Users/rodney/Desktop/%s.colormag.pdf'%name)
    else : snCM=None

    sn = snana.SuperNova( 'HST_%s_%s.dat'%(survey,name) )
    if LCgrid and sn.zerr > 0.2 : 
        print(" zerr = %.2f... skipping Grid classification"%(sn.zerr))

    elif LCgrid : 
        print( " %s : doing the grid classification and plots... "%name )
        if sn.zerr<0.02 : nlogz=1
        elif sn.zerr<0.11 : nlogz=5
        else : nlogz=10
        sn.doGridClassify( clobber=clobber, verbose=verbose, nlogz=nlogz, ncolorpar=20, nlumipar=30, ncolorlaw=1, npkmjd=40, useLuminosityPrior=True, trestrange=[-20,60], modelerror=[0.08,0.1,0.1] )
        figure(3)
        clf()
        figure(3)
        sn.plotMaxLikeModels( 'flux','mjd', showchi2=False )
        if interactive: 
            userin = raw_input("Adjust plot spacing. Hit <return> to save and continue.")
        savefig('/Users/rodney/Desktop/%s.LCfit.pdf'%name)
        figure(4)
        close()
        sn.plotClassStatsGrid()
        savefig('/Users/rodney/Desktop/%s.stats.pdf'%name)
    else : sn=None

    print( "All done with %s "%name)
    return(snCM, sn)

def athenaPlots( sn ):
    """ plotting SN athena classification deep diagnostics to figure out why she's called a Ibc """
    from hstsnpipe.tools import snana
    if isinstance( sn, basestring ) : 
        sn = snana.SuperNova( sn ) 
        sn.doGridClassify( clobber=0, verbose=3, nlogz=1, ncolorpar=20, nlumipar=20, ncolorlaw=1, npkmjd=30, 
                           useLuminosityPrior=True, trestrange=[-20,80], modelerror=[0.08,0.1,0.1] )
    simIa = sn.ClassSim.Ia
    simIbc = sn.ClassSim.Ibc


    # plot Like vs Magoff, x1, c, Av, 
    figure(1)
    clf()
    ax1 = subplot( 3, 2, 1 )
    semilogy( simIa.priorMagoff, simIa.LIKE, 'ro', alpha=0.1 )
    xlabel('prior( Magoff )' )
    ax2 = subplot( 3, 2, 2, sharex=ax1, sharey=ax1 )
    semilogy( simIbc.priorMagoff, simIbc.LIKE, 'go', alpha=0.1 )
    xlabel('prior( Magoff )' )

    ax1 = subplot( 3, 2, 3 )
    semilogy( simIa.pc, simIa.LIKE, 'ro', alpha=0.1 )
    xlabel('prior( c )' )

    ax2 = subplot( 3, 2, 4, sharex=ax1, sharey=ax1 )
    semilogy( simIbc.pAv, simIbc.LIKE, 'go', alpha=0.1 )
    xlabel('prior( Av )' )

    ax1 = subplot( 3, 2, 5 )
    semilogy( simIa.px1, simIa.LIKE, 'ro', alpha=0.1 )
    xlabel('prior( x1 )' )

    figure(2)
    clf()
    ax1 = subplot( 3, 2, 1 )
    semilogy( simIa.priorMagoff, simIa.postProb, 'ro', alpha=0.1 )
    ax2 = subplot( 3, 2, 2, sharex=ax1, sharey=ax1 )
    semilogy( simIbc.priorMagoff, simIbc.postProb, 'go', alpha=0.1 )

    ax1 = subplot( 3, 2, 3 )
    semilogy( simIa.pc, simIa.postProb, 'ro', alpha=0.1 )
    ax2 = subplot( 3, 2, 4, sharex=ax1, sharey=ax1 )
    semilogy( simIbc.pAv, simIbc.postProb, 'go', alpha=0.1 )

    ax1 = subplot( 3, 2, 5 )
    semilogy( simIa.px1, simIa.postProb, 'ro', alpha=0.1 )


def mkColorMagPlot( Nsim=2000, dustmodel='mid', ratemodel='mid', 
                    namelabels=False, clobber=False, verbose=True, **kwargs ):
    """ plot color-mag diagrams for 3-4 redshift bins, 
    overplotting the near-peak points for all the GOODS SNe"""
    from hstsnpipe.tools import snana

    from hstsnpipe.tools.figs import plotsetup
    plotsetup.fullpaperfig( 1, figsize=[8,4] )
    fig = gcf()
    fig.subplots_adjust( left=0.08, bottom=0.13, right=0.98, top=0.9, wspace=0.27, hspace=0.1 )
   
    # run a simulation for this redshift bin
    allbands='XVZEF'
    iax=0
    for zbin,color,mag in zip( [0.25,0.75,1.25],
                               ['X-Z','X-Z','X-E'],
                               ['X',  'Z',  'E' ] ): 
        iax +=1
        zrange = [max(0.001,zbin-0.25),zbin+0.25]
        simroot = 'goodsColorMag_z%.2f'%zbin
        
        dosim=False
        if not clobber : 
            try : 
                simset = snana.classify.rdClassSim( simroot, verbose=verbose )
                dosim=False
            except : 
                dosim=True
        if dosim or clobber : 
            if verbose: print("running snana.simulate")
            if verbose>2 : print( 'dosim=%s  clobber=%s'%(dosim,clobber))
            simroot = snana.classify.doMonteCarloSim( 
                simroot=simroot, Nsim=Nsim, bands=allbands, 
                zrange=zrange, avrange=[0,3], 
                mjd0=55130, Nepoch=1, cadence=1, pkmjdrange=[55100,55160],
                dustmodel=dustmodel, ratemodel=ratemodel, 
                clobber=clobber, verbose=verbose, perfect=True )
            simset = snana.classify.rdClassSim( simroot, verbose=verbose)

        if verbose: print( "plotting panel %i"%iax)
        subplot( 1,3, iax)
        snana.simplot.plotSimClouds( simset.II,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ibc,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ia,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        ax = gca()
        if iax==1 : 
            ax.set_xlim( -1.4,1.4 )
            ax.set_ylim( 28.5,16.5 )
        elif iax==2 : 
            ax.set_xlim( -1.2,1.7 )
            ax.set_ylim( 28.5,20.55 )
        elif iax==3 : 
            ax.set_xlim( -0.9,4.5 )
            ax.set_ylim( 28.5,21.5 )

        Nsn = plotColorMagObs(color,mag, zrange=zrange, verbose=verbose, namelabels=namelabels )
        title( "%i SNe with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )


def mkColorMagPlotSAMPLER( Nsim=2000, dustmodel='mid', ratemodel='mid', namelabels=True, 
                           clobber=False, verbose=True, **kwargs ):
    """ plot color-mag diagrams for 3-4 redshift bins, 
    overplotting the near-peak points for all the GOODS SNe"""
    from hstsnpipe.tools import snana

    from hstsnpipe.tools.figs import plotsetup
    plotsetup.fullpaperfig( 1, figsize=[12,10] )
    fig = gcf()
    fig.subplots_adjust( left=0.12, bottom=0.12, right=0.95, top=0.9, wspace=0.1, hspace=0.1 )
   
    
    # run a simulation for this redshift bin
    allbands='XVZEF'
    iax=0
    for zbin in [0.25,0.75,1.25] : 
        iax +=1
        zrange = [max(0.001,zbin-0.25),zbin+0.25]
        simroot = 'goodsColorMag_z%.2f'%zbin
        
        dosim=False
        if not clobber : 
            try : 
                simset = snana.classify.rdClassSim( simroot, verbose=verbose )
                dosim=False
            except : 
                dosim=True
        if dosim or clobber : 
            if verbose: print("running snana.simulate")
            if verbose>2 : print( 'dosim=%s  clobber=%s'%(dosim,clobber))
            simroot = snana.classify.doMonteCarloSim( 
                simroot=simroot, Nsim=Nsim, bands=allbands, 
                zrange=zrange, avrange=[0,3], 
                mjd0=55130, Nepoch=1, cadence=1, pkmjdrange=[55100,55160],
                dustmodel=dustmodel, ratemodel=ratemodel, 
                clobber=clobber, verbose=verbose, perfect=True )
            simset = snana.classify.rdClassSim( simroot, verbose=verbose)

        if verbose: print( "plotting panel %i"%iax)
        subplot( 2,4, iax)
        snana.simplot.plotSimClouds( simset.II,  xaxis='X-Z', yaxis='Z', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], **kwargs )
        snana.simplot.plotSimClouds( simset.Ibc,  xaxis='X-Z', yaxis='Z', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], **kwargs )
        snana.simplot.plotSimClouds( simset.Ia,  xaxis='X-Z', yaxis='Z', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], **kwargs )

        ax = gca()
        ax.set_xlim( -1,2 )
        ax.set_ylim( 28,19 )

        Nsn = plotColorMagObs('X-Z','Z', zrange=zrange, verbose=verbose, namelabels=namelabels )
        title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )

    subplot( 2,4, 4 )
    snana.simplot.plotSimClouds( simset.II,  xaxis='Z-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ibc,  xaxis='Z-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ia,  xaxis='Z-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    Nsn = plotColorMagObs('Z-E','E', zrange=zrange, verbose=verbose, namelabels=namelabels )
    title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )

    subplot( 2,4, 5 )
    snana.simplot.plotSimClouds( simset.II,  xaxis='X-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ibc,  xaxis='X-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ia,  xaxis='X-E', yaxis='E', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    Nsn = plotColorMagObs('X-E','E', zrange=zrange, verbose=verbose, namelabels=namelabels )
    title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )

    subplot( 2,4, 6 )
    snana.simplot.plotSimClouds( simset.II,  xaxis='Z-E', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ibc,  xaxis='Z-E', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ia,  xaxis='Z-E', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    Nsn = plotColorMagObs('X-E','F', zrange=zrange, verbose=verbose, namelabels=namelabels )
    title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )

    subplot( 2,4, 7 )
    snana.simplot.plotSimClouds( simset.II,  xaxis='Z-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ibc,  xaxis='Z-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ia,  xaxis='Z-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    Nsn = plotColorMagObs('X-E','F', zrange=zrange, verbose=verbose, namelabels=namelabels )
    title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )

    subplot( 2,4, 8 )
    snana.simplot.plotSimClouds( simset.II,  xaxis='X-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ibc,  xaxis='X-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    snana.simplot.plotSimClouds( simset.Ia,  xaxis='X-F', yaxis='F', snmags={}, 
                                 linelevels=[0.95,0.68,0], sidehist=False, 
                                 mjdrange=[55130,55130], **kwargs )
    Nsn = plotColorMagObs('X-E','F', zrange=zrange, verbose=verbose, namelabels=namelabels )
    title( "%i SNIa with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )




def mkColorMagPlotCANDELS( Nsim=2000, dustmodel='mid', ratemodel='mid', 
                           namelabels=False, clobber=False, verbose=True, **kwargs ):
    """ plot color-mag diagrams for 3-4 redshift bins, 
    overplotting the near-peak points for all the CANDELS SNe"""
    from hstsnpipe.tools import snana

    from hstsnpipe.tools.figs import plotsetup
    plotsetup.fullpaperfig( 1, figsize=[8,4] )
    fig = gcf()
    fig.subplots_adjust( left=0.08, bottom=0.13, right=0.98, top=0.9, wspace=0.27, hspace=0.1 )
   
    # run a simulation for this redshift bin
    allbands='VIWJH'
    iax=0
    for zbin,color,mag in zip( [0.25,0.75,1.25,1.75],
                               ['W-J','W-J','W-J','W-J'],
                               ['J',  'J',  'J' , 'J'] ): 
        iax +=1
        zrange = [max(0.001,zbin-0.25),zbin+0.25]
        simroot = 'candelsColorMag_z%.2f'%zbin
        
        dosim=False
        if not clobber : 
            try : 
                simset = snana.classify.rdClassSim( simroot, verbose=verbose )
                dosim=False
            except : 
                dosim=True
        if dosim or clobber : 
            if verbose: print("running snana.simulate")
            if verbose>2 : print( 'dosim=%s  clobber=%s'%(dosim,clobber))
            simroot = snana.classify.doMonteCarloSim( 
                simroot=simroot, Nsim=Nsim, bands=allbands, 
                zrange=zrange, avrange=[0,3], 
                mjd0=55130, Nepoch=1, cadence=1, pkmjdrange=[55100,55160],
                dustmodel=dustmodel, ratemodel=ratemodel, 
                clobber=clobber, verbose=verbose, perfect=True )
            simset = snana.classify.rdClassSim( simroot, verbose=verbose)

        if verbose: print( "plotting panel %i"%iax)
        subplot( 2,4, iax)
        snana.simplot.plotSimClouds( simset.II,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ibc,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ia,  xaxis=color, yaxis=mag, snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        ax = gca()
        if iax==1 : 
            ax.set_xlim( -1.4,1.4 )
            ax.set_ylim( 28.5,16.5 )
        elif iax==2 : 
            ax.set_xlim( -1.2,1.7 )
            ax.set_ylim( 28.5,20.55 )
        elif iax==3 : 
            ax.set_xlim( -0.9,4.5 )
            ax.set_ylim( 28.5,21.5 )

        Nsn = plotColorMagObs(color,mag, sample='candels', zrange=zrange, verbose=verbose, namelabels=namelabels )
        title( "%i SNe with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )
        
        subplot( 2,4, iax+4)
        snana.simplot.plotSimClouds( simset.II,  xaxis='I-J', yaxis='J', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ibc,  xaxis='I-J', yaxis='J', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        snana.simplot.plotSimClouds( simset.Ia,  xaxis='I-J', yaxis='J', snmags={}, 
                                     linelevels=[0.95,0.68,0], sidehist=False, 
                                     mjdrange=[55130,55130], binrange=[[-3,8],[16,32]], bins=50,
                                     **kwargs )
        ax = gca()
        if iax==1 : 
            ax.set_xlim( -1.4,1.4 )
            ax.set_ylim( 28.5,16.5 )
        elif iax==2 : 
            ax.set_xlim( -1.2,1.7 )
            ax.set_ylim( 28.5,20.55 )
        elif iax==3 : 
            ax.set_xlim( -0.9,4.5 )
            ax.set_ylim( 28.5,21.5 )

        Nsn = plotColorMagObs('I-J','J', sample='candels', zrange=zrange, verbose=verbose, namelabels=namelabels )
        title( "%i SNe with $%.1f < z < %.1f$"%( Nsn, zrange[0], zrange[1] ) )



def plotColorMagObs( color='X-Z', mag='Z', zrange=[0.001, 0.5], namelabels=False, 
                     sample='goods', showcc=False, showlc=False, verbose=False):
    """ plot observed color-mag points from real SNe with observations within 30 days 
    of their nominal peak brightness """

    band1 = color[0]
    band2 = color[-1]
    band3 = mag
    colormagbands = unique( [band1,band2,band3] )

    pkmag = {}
    pkmagerr = {}
    z = {}

    pkcolor= {}
    pkcolerr={}
    pkmag={}
    pkmagerr={}
    pkmjd={}
    z={}

    if sample=='goods' : 
        datfilelist = iadatfiles + ccdatfiles
    elif sample=='candels' : 
        datfilelist = candelsdatfiles

    Nsn = 0
    for datfile in datfilelist : 
        sn = snana.SuperNova( datfile ) 
        if not sn.z > zrange[0] : continue
        if not sn.z <= zrange[1] : continue
        
        name = sn.nickname
        if name in goldIalist : 
            color='goldenrod'
            marker='D'
        elif name in silverIalist : 
            color='0.8'
            marker='D'
        elif name in goldCClist + silverCClist + bronzeCClist : 
            color='c'
            marker='o'
        elif sample=='candels': 
            color='k'
            marker='o'
        else : 
            if verbose>5 : print("skipping %s  (not in gold/silver list)"%name) 
            continue
        
        # locate coincident observation sets (i.e. w/in 5 days)
        # that have all three bands of interest
        mjdlist, colorlist, maglist, colerrlist, magerrlist = [], [], [], [], []
        iband1 = where( sn.FLT == band1 )[0] 
        iband2 = where( sn.FLT == band2 )[0] 
        iband3 = where( sn.FLT == band3 )[0] 
        if any( [ len(iband1)==0,  len(iband2)==0, len(iband3)==0 ] ): 
            if verbose>5 : print("skipping %s  (no good color-mag epochs)"%name) 
            continue
        for i1 in iband1 : 
            if abs( sn.MJD[i1] - sn.pkmjd ) > 30 : continue
            dmjd2 = abs( sn.MJD[iband2] - sn.MJD[i1] ) 
            dmjd3 = abs( sn.MJD[iband3] - sn.MJD[i1] ) 
            if dmjd2.min() > 5 or dmjd3.min() > 5 : continue
            imjd2min = dmjd2.argmin()
            imjd3min = dmjd3.argmin()
            mjd1 = sn.MJD[i1]
            mjd2 = sn.MJD[iband2][imjd2min] 
            mjd3 = sn.MJD[iband3][imjd3min] 
            thismjd = mean([mjd1,mjd2,mjd3])
            if thismjd in mjdlist : continue
            thiscolor = sn.MAG[i1]- sn.MAG[iband2][imjd2min] 
            thiscolerr = sqrt(sn.MAGERR[i1]**2 + sn.MAGERR[iband2][imjd2min]**2)
            thismag = sn.MAG[iband3][imjd3min] 
            thismagerr = sn.MAGERR[iband3][imjd3min] 

            colorlist.append( thiscolor )
            colerrlist.append( thiscolerr )
            maglist.append( thismag ) 
            magerrlist.append( thismagerr ) 
            mjdlist.append( thismjd )

        if showlc : 
            figure(100)
            clf()
            sn.plotLightCurve( 'mag','mjd', showpkmjdrange=True, autozoom=False )
            plot( mjdlist, maglist, 
                  marker='D', ms=14, alpha=0.3, ls=' ', color='yellow', mew=0.5,mec='k')
            draw()
            userin=raw_input("Showing %s. \n   hit <return> to continue"%sn.nickname)
            if userin=='pdb' : import pdb; pdb.set_trace()
            close(100)

        if not len(colorlist) : 
            if verbose>5 : print("skipping %s  (no good color-mag epochs)"%name) 
            continue

        if len(colorlist): 
            if namelabels : 
                errorbar( colorlist, maglist, magerrlist, colerrlist, 
                          ls=' ', marker=' ', mfc=color, mec='w', mew=2, ms=8, color='k' )
                for i in range(len(colorlist)) :
                    text( colorlist[i], maglist[i], name[:3], 
                          color='k', backgroundcolor=color, ha='center',va='center', alpha=0.4 )
            else : 
                errorbar( colorlist, maglist, magerrlist, colerrlist, 
                          ls=' ', marker=marker, mfc=color, mec='w', mew=2, ms=8, color=color )

        Nsn += 1
                
    return( Nsn  ) 



def plotMiniMagLC( sn, dustmodel='mid', inflateErrors=1, 
                   interactive=True, maxtype='probability', showclasstable=True,
                   bands='all', plotPosteriors=False, clobber=False, verbose=3 ) : 
    """ For the given supernova, make a "miniature" light curve plot
    showing magnitudes vs mjd, overlaid with the best-fit light curve
    from the sub-class that has the largest posterior probability.
    E.g. if the given SN has P(Ia) = 0.7, then the best-fit Type Ia
    model is shown. 

    If needed, run doGridClassify first to get classification
    probabilities, using reasonable default values. If the default
    doGridClassify settings don't provide enough sampling, then run
    doGridClassify separately and pass this function a SuperNova
    object that already has posterior probabilities defined.

    sn : Name of the supernova .dat file, or a snana.SuperNova object

    dustmodel : 'mid', 'high', 'low', or a function from snana.extinction
    errscale : factor by which to scale the modelerror values used in 
               computing the posterior probabilities (if the best-fit chi2 
               values are large, you can use errscale>1 to get more realistic 
               posterior probabilities)
    interactive : if True, prompt the user for adjustments of the plot window 
               axes, and ask where to save the final figure.  When the
               prompt asks for X/Y limits, enter two numbers separated
               by a space to give the MJD and mag plot limits,
               respectively.
    maxtype  : ['probability'] plot light curves that maximize the posterior 
                   probability (i.e. include priors)
               ['likelihood'] plot max-likelihood (min-chi2) curves 
    clobber : [1] force a re-calculation of the posterior probabilities, and 
                  re-make the light curve plots
              [2] re-read the SNANA simulation results, recompute probs, 
                  and re-make plots
              [3+] force a re-run of the SNANA simulations, recompute probs, 
                  and re-make plots
    """
    from  matplotlib import pyplot as pl
    from  matplotlib import ticker
    import numpy as np
    from hstsnpipe.tools.figs import plotsetup
    import os 
    import candelsClassify
    

    if showclasstable : 
        plotsetup.halfpaperfig( 1, figsize=[4.2,2.5] )
    else : 
        plotsetup.halfpaperfig( 1, figsize=[2.5,2] )
       
    if isinstance( sn, str ):
        if '.dat' not in sn : sn = "HST_CANDELS1_%s.dat"%sn
        sn = snana.SuperNova( sn )

    if bands=='all' : bands = sn.bands
    if 'PIa' not in sn.__dict__ :
        sn = candelsClassify.doOneSN( sn, inflateErrors=inflateErrors, returnsn=True, bands=bands, clobber=clobber, verbose=verbose-1 )

    if sn.PIa > 0.5 : 
        sn.plotLightCurve("mag","mjd", showclassfit='Ia.max%s'%maxtype[:4], showlegend=True, showclasstable=showclasstable, bands=bands );
    elif sn.PIbc > sn.PII : 
        sn.plotLightCurve("mag","mjd", showclassfit='Ibc.max%s'%maxtype[:4], showlegend=True, showclasstable=showclasstable, bands=bands );
    else : 
        sn.plotLightCurve("mag","mjd", showclassfit='II.max%s'%maxtype[:4], showlegend=True, showclasstable=showclasstable, bands=bands );

    fig = pl.gcf()

    if len(bands)==1 : ncol=1
    if len(bands)==2 : ncol=2
    if len(bands)==3 : ncol=3
    if len(bands)==4 : ncol=4
    if len(bands)==5 : ncol=3
    if len(bands)==6 : ncol=3
    if len(bands)==7 : ncol=4
    if len(bands)>7 : ncol=4
    ax = pl.gca()
    leg = ax.legend( loc='upper left', title=r'%s'%(sn.name),
                     ncol=ncol, columnspacing=0.25, 
                     frameon=False, bbox_to_anchor=(0.1,1.3), borderaxespad=0, 
                     numpoints=1, handlelength=0.2, handletextpad=0.4, labelspacing=0.2  )
    pl.setp( leg.get_title(), fontsize='large' )
    pl.setp( leg.get_texts(), fontsize='small' )

    ax.set_xlabel('MJD :', ha='right')
    ax.xaxis.set_label_coords( -0.005, -0.025)
    ax.set_ylabel('Vega mag')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax.xaxis.set_major_locator( ticker.MultipleLocator( 30 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 15 ) )
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )
    if showclasstable : 
        fig.subplots_adjust( left=0.13, right=0.73, bottom=0.12, top=0.77 )
        ax.texts[-2].set_transform( ax.transAxes )
        ax.texts[-2].set_verticalalignment('top')
        ax.texts[-2].set_horizontalalignment('left')
        ax.texts[-2].set_position( [1.02, 1.3] )
        ax.texts[-1].set_transform( ax.transAxes )
        ax.texts[-1].set_verticalalignment('top')
        ax.texts[-1].set_horizontalalignment('left')
        ax.texts[-1].set_position( [1.02, 0.85] )

    else : 
        fig.subplots_adjust( left=0.1, right=0.95, bottom=0.13, top=0.78 )

    ax.set_ylim( [28.9, 23.0] )
    ax.set_xlim( [ sn.maxProbIaModel.pkmjd- 15*(1+sn.z), sn.maxProbIaModel.pkmjd+35*(1+sn.z)] )


    savedir = os.path.abspath( os.path.expanduser("~/Dropbox/Papers/hstsnrates/FIG/LCFIT/"))
    if not os.path.isdir( savedir ) : 
        savedir = os.path.abspath( '.' ) 

    if plotPosteriors : 
        sn.plotClassStatsGrid()
        pl.savefig(os.path.join( savedir, '%s.posteriors.pdf'%( sn.name)) )
    pl.draw()

    fig = pl.figure(1)
    savefile = os.path.join( savedir, "%sLC.pdf"%sn.nickname.lower())
    if interactive: 
        loop = True
        xticker=30
        yticker=1
        while loop : 
            xlim = [ float(s) for s in raw_input("Set X lim :  ").split() ]
            if len(xlim) : 
                ax.set_xlim( xlim ) 
                pl.draw()
            ylim = [ float(s) for s in raw_input("Set Y lim :  ").split() ]
            if len(ylim) : 
                ax.set_ylim( ylim ) 
                pl.draw()
            xtickerin = raw_input("Set X tick mark spacing [%i] : "%xticker )
            if xtickerin :
                xticker= float(xtickerin)
                ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
                ax.xaxis.set_major_locator( ticker.MultipleLocator( xticker) )
                ax.xaxis.set_minor_locator( ticker.MultipleLocator( xticker/2 ) )
                pl.draw()
            ytickerin = raw_input("Set Y tick mark spacing [%i] : "%yticker)
            if ytickerin :
                yticker=float(ytickerin)
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
                ax.yaxis.set_major_locator( ticker.MultipleLocator( yticker ) )
                ax.yaxis.set_minor_locator( ticker.MultipleLocator( yticker/2 ) )
                pl.draw()
            if not len(xlim) and not len(ylim) and not xtickerin and not ytickerin : break
        newsavefile = raw_input("Filename for saving [%s] :  "%savefile )
        if len(newsavefile.strip())!=0 : savefile = newsavefile
    savefig( savefile )
    if savefile.endswith('.pdf') : savefig( savefile.replace('.pdf','.eps') )
    return( sn ) 


# STARDUST Classifications of the GOODS+PANS SNe 
# Updated 2013.08.02 after major bug-fix in classifer
#BronzeCC = np.array([ 0.0, 0.0, 0.0, 0.1, 0.12, 0.54, 0.9 ])
SilverCC = np.array([0.00,0.00,0.00,0.00,0.00,0.02,0.15,0.06,0.93])
SilverCCpk = np.array([0.10,0.10,0.00,0.16,0.03,0.02,0.42,0.41,0.93])

GoldIa = np.array([0.00,1.00,0.01,0.55,0.76,0.86,0.95,0.95,0.98,0.99,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00])
GoldIapk = np.array([0.36,1.00,0.33,0.99,0.68,0.97,0.83,0.98,0.54,0.97,1.00,0.92,0.97,0.72,0.96,1.00,1.00,0.88,1.00,0.98,0.99,0.98,0.97,0.94,0.91,0.88,1.00])

SilverIa= np.array([0.58,0.00,0.01,0.04,0.31,0.43,0.91,0.99,1.00,1.00,1.00,1.00])
SilverIapk = np.array([0.34,0.00,0.34,0.47,0.91,0.43,0.86,0.51,0.07,1.00,1.00,0.57])



def plotGOODSPANShist():
    """ make a simple histogram showing classification validation test results """
    
    from matplotlib import pyplot as pl
    import numpy as np
    import plotsetup

    pbins = np.arange( 0, 1.01, 0.1 )
    gold, pbins = np.histogram( GoldIa, bins=pbins )
    silver, pbins = np.histogram( SilverIa, bins=pbins )
    # bronze, pbins = np.histogram( BronzeIa, bins=pbins )

    cc, pbins= np.histogram( SilverCC, bins=pbins)

    fig = plotsetup.halfpaperfig(1, figsize=[8,3])
    ax1 = fig.add_subplot(121)
    ax1.plot( pbins, np.append(gold,gold[-1]), color='darkorange', drawstyle='steps-post', lw=3, ls='-' )
    ax1.plot( pbins, np.append(silver,silver[-1]), color='0.3', drawstyle='steps-post', lw=3, ls='--' )
    ax1.plot( pbins, np.append(cc,cc[-1]), color='darkcyan', drawstyle='steps-post', lw=3, ls=':' )

    ax1.set_xlabel('P(Ia) from Light Curve Fitting')
    ax1.set_ylabel('Number of SNe')
    ax1.set_xlim( -0.05, 1.05 )
    ax1.text(0.88,15,'Gold SNe Ia', color='darkorange',fontweight='bold',fontsize='x-large', ha='right', va='top')
    ax1.text(0.85,9,'Silver SNe Ia', color='0.3',fontweight='bold',fontsize='x-large', ha='right', va='top')
    ax1.text(0.15,5,'Silver CC SNe', color='darkcyan',fontweight='bold',fontsize='x-large', ha='left', va='top')
    ax1.text( 0.08, 0.93, 'Classifier Validation Test:\n 48 SNe from HST+ACS Surveys\nUsing Full Light Curve', fontsize='large',
             transform=ax1.transAxes, va='top', ha='left')
    fig.subplots_adjust(left=0.1,right=0.96, top=0.96, bottom=0.17)


    gold, pbins = np.histogram( GoldIapk, bins=pbins )
    silver, pbins = np.histogram( SilverIapk, bins=pbins )
    cc, pbins= np.histogram( SilverCCpk, bins=pbins)
  
    ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)
    ax2.plot( pbins, np.append(gold,gold[-1]), color='darkorange', drawstyle='steps-post', lw=3, ls='-' )
    ax2.plot( pbins, np.append(silver,silver[-1]), color='0.3', drawstyle='steps-post', lw=3, ls='--' )
    ax2.plot( pbins, np.append(cc,cc[-1]), color='darkcyan', drawstyle='steps-post', lw=3, ls=':' )

    ax2.set_xlabel('P(Ia) from Light Curve Fitting')
    ax2.set_ylabel('Number of SNe')
    ax2.set_xlim( -0.05, 1.05 )
    ax2.set_ylim( -0.05, 26.5 )
    #ax2.text(0.88,15,'Gold SNe Ia', color='darkorange',fontweight='bold',fontsize='x-large', ha='right', va='top')
    #ax2.text(0.85,9,'Silver SNe Ia', color='0.3',fontweight='bold',fontsize='x-large', ha='right', va='top')
    #ax2.text(0.15,5,'Silver CC SNe', color='darkcyan',fontweight='bold',fontsize='x-large', ha='left', va='top')
    ax2.text( 0.08, 0.93, 'Using Only Peak Data', fontsize='large',
             transform=ax2.transAxes, va='top', ha='left')
    ax2.yaxis.set_ticks_position( 'right' )
    ax2.yaxis.set_ticks_position( 'both' )
    ax2.yaxis.set_label_position( 'right' )

    fig.subplots_adjust(left=0.1,right=0.96, top=0.96, bottom=0.17, wspace=0. )
    pl.draw()
