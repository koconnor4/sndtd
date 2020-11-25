__all__=[ "csfr", "snratefigs" ]
import snratefigs
import csfr

from numpy import arange,array, zeros, sqrt, abs, max, min



class rate(object):
    def __init__(self):
        self.reference=''
        self.z=array([])
        self.zerrplus=array([])
        self.zerrminus=array([])

        self.rate=array([])
        self.errstatplus=array([])
        self.errsysplus=array([])

        self.errstatminus=array([])
        self.errsysminus=array([])
        
        self.elw=1.5

    @property
    def tex(self):
        """ print a tex string for a data table """
        texset = ""
        for i in range(len(self.z)) :
            ztex = "%.2f"%self.z[i]
            if self.zerrplus[i] == self.zerrminus[i] : 
                zerrtex = "$\pm$%.2f"%self.zerrplus[i] 
            else : 
                zerrtex = "$\pm^{%.2f}_{%.2f}$"%(self.zerrplus[i],self.zerrminus[i])
            snrtex = "%.2f"%self.rate[i]
            if self.errplus[i] == self.errminus[i] : 
                snrerrtex = "$\pm$%.2f"%self.errplus[i] 
            else : 
                snrerrtex = "$\pm^{%.2f}_{%.2f}$"%(self.errplus[i],self.errminus[i])
            reftex = "\citealt{%s:%s}"%(self.reference.split()[0],self.reference.split()[-1])
            texset += "%s %s & %s %s & %s \\\\ \n"%( 
                ztex, zerrtex, snrtex, snrerrtex, reftex )
        texset = texset[:-2]
        return(texset)

    @property
    def errplus(self):
        from numpy import sqrt
        return( sqrt( self.errstatplus**2 + self.errsysplus**2 )  )
        # return( self.errstatplus + self.errsysplus )
    
    @property
    def errminus(self):
        from numpy import sqrt, where
        errminus1 = sqrt( self.errstatminus**2 + self.errsysminus**2 )
        return( where( self.rate-errminus1>0, errminus1, self.rate*0.99 ) )
        # return( self.errstatminus + self.errsysminus )
        # return( array( [ min([errvec[i],self.rate[i]-0.1]) for i in range(len(errvec))] ) )

    @property
    def err(self):
        from numpy import array
        errvec = array([self.errminus, self.errplus])
        return(errvec)

    @property
    def zerr(self):
        from numpy import array
        zerrvec = array([self.zerrminus, self.zerrplus])
        return(zerrvec)

    @property
    def meanerr(self):
        from numpy import abs
        mnerr = (abs(self.errminus) + abs(self.errplus)) / 2
        return(mnerr)
        
    def plot(self,  scalerates=1, showerr=True, showzerr=False,
             zshift=0, systematicshift=None, zorder=0,  **kwargs):
        """ plot with double error bars : 
        thicker lighter systematic errs behind thinner darker composite errs

        Set zorder to a positive value to plot in the foreground, and set 
        it to a negative value for background plotting. 

        showerr :  [ True, False, 'sys', 'stat', 'thicksys' ]

        Default units are SNuVol = 10^-4 yr-1 Mpc-3
        Set scalerates to 1e-4 in order to plot in units of yr-1 Mpc-3
        """
        from pylab import plot,errorbar

        if zorder<0 : 
            ezorder=zorder*10
            zorder=-2
        elif zorder>0 : 
            ezorder=10
        elif zorder==0 : 
            ezorder=-1
            zorder=1

        redshift = self.z + zshift
        rates = self.rate * scalerates

        if systematicshift in ['high','up','+'] :
            rates = rates + self.errsysplus * scalerates
        elif systematicshift in ['low','down','-'] :
            rates = rates - self.errsysminus * scalerates

        if showzerr : zerr = self.zerr
        else : zerr= self.zerr * 0


        if not showerr :
            err = zeros(len(self.rate))
        else :
            err = array( [self.errminus, self.errplus] ) * scalerates
            if isinstance( showerr, str ) :
                if showerr.lower().startswith('stat') :
                    err = array([self.errstatminus, self.errstatplus]) * scalerates
                elif showerr.lower().find('sys')>=0 :
                    err = array([self.errsysminus, self.errsysplus]) * scalerates
                    if showerr.lower().startswith('thick') :
                        # thick lines for Systematic errors
                        errorbar( redshift, rates, err,
                                  ls=' ', marker=' ', ecolor=self.mfc, color=self.color,
                                  barsabove=False, elinewidth=self.elw*5, capsize=0,
                                  zorder=zorder, ms=self.ms, **kwargs)
                        err = array( [self.errminus, self.errplus] ) * scalerates

        errorbar( redshift, rates, err, zerr,
                 ls=' ', marker=' ', color=self.color,
                 barsabove=False, elinewidth=self.elw, ecolor=self.mec, capsize=0,
                 zorder=ezorder, **kwargs)
        plot( redshift, rates, ls=' ',
             marker=self.marker, color=self.color, mfc=self.mfc, mec=self.mec, 
             zorder=zorder+1, ms=self.ms,  **kwargs)


    def plotSysErr(self, thicksys=True, zorder=0, scaleunits=1, **kwargs):
        """ plot showing only the systematic error bars 

        Set zorder to a positive value to plot in the foreground, and set 
        it to a negative value for background plotting. 
        """
        from pylab import plot,errorbar,scatter

        if zorder<0 : 
            ezorder=zorder*10
            zorder=-2
        elif zorder>0 : 
            ezorder=10
        elif zorder==0 : 
            ezorder=-1
            zorder=1

        if thicksys : 
            # thick lines for Systematic errors 
            errorbar(self.z, scaleunits*self.rate, 
                     scaleunits*array([self.errsysminus, self.errsysplus]), 
                     ls=' ', marker=' ', ecolor=self.mfc, color=self.color,
                     barsabove=False, elinewidth=self.elw*5, capsize=0,
                     zorder=zorder, ms=self.ms, **kwargs)
        else : 
            # thin lines for systematic errors
            errorbar(self.z, scaleunits*self.rate, 
                     scaleunits*array([self.errstatminus, self.errstatplus]), 
                     self.zerr, ls=' ', marker=' ', color=self.color, 
                     barsabove=False, elinewidth=self.elw, ecolor=self.mec, capsize=0,
                     zorder=ezorder, mfc=self.mfc, mec=self.mec, ms=self.ms,  **kwargs)
        plot(self.z, scaleunits*self.rate, ls=' ', marker=self.marker, color=self.color, 
             zorder=zorder+1, mfc=self.mfc, mec=self.mec, ms=self.ms,  **kwargs)


    def plotStatErr(self, zorder=0, scaleunits=1, **kwargs):
        """ plot showing only the statistical error bars 

        Set zorder to a positive value to plot in the foreground, and set 
        it to a negative value for background plotting. 
        """
        from pylab import plot,errorbar,scatter

        if zorder<0 : 
            ezorder=zorder*10
            zorder=-2
        elif zorder>0 : 
            ezorder=10
        elif zorder==0 : 
            ezorder=-1
            zorder=1

        # statistical errors
        errorbar(self.z, scaleunits*self.rate, 
                 scaleunits*array([self.errstatminus, self.errstatplus]), 
                 self.zerr, ls=' ', marker=' ', color=self.color, 
                 barsabove=False, elinewidth=self.elw, ecolor=self.mec, capsize=0,
                 zorder=ezorder, mfc=self.mfc, mec=self.mec, ms=self.ms,  **kwargs)
        plot(self.z, scaleunits*self.rate, ls=' ', marker=self.marker, color=self.color, 
             zorder=zorder+1, mfc=self.mfc, mec=self.mec, ms=self.ms,  **kwargs)


    def logzplot(self,**kwargs):
        from pylab import plot,errorbar
        from numpy import log10,log
        errorbar(log10(1+self.z), self.rate, 
                 self.err, self.zerr/(1+self.z)/log(10), 
                 ls=' ', marker=self.marker, mfc=self.mfc,
                 mec=self.mec, color=self.mec, 
                 barsabove=False,label=self.reference,
                 ms=self.ms, 
                 capsize=0, **kwargs)

    def printtable( self ):
        print( self.reference + ' : ')
        print( '    z   Ndet     SNR      stat          syst ' )
        for i in range(len(self.z) ): 
            print( "%5.2f  %5.2f   %5.2f  +%4.2f -%4.2f   +%4.2f -%4.2f"%(
                    self.z[i], self.ndet[i], self.rate[i], 
                    self.errstatplus[i], self.errstatminus[i],
                    self.errsysplus[i], self.errsysminus[i] ) )
