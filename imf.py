import numpy as np

def kroupa93(m) : 
    if not np.iterable( m ) : m = [m]
    i0 = np.where( m < 0.08 ) 
    i1 = np.where( (0.08<=m) & (m<0.5) )
    i2 = np.where( (0.5<=m) & (m<1.0) )
    i3 = np.where( 1 <= m )
              
    xi =  m[i0] * 0
    xi = np.append( xi, 0.579 * m[i1] ** -1.3 )
    xi = np.append( xi, 0.310 * m[i2] ** -2.2 )
    xi = np.append( xi, 0.310 * m[i3] ** -2.7 )

    return( xi )


def dietSalpeter(m) : 
    if not np.iterable( m ) : m = [m]
    i0 = np.where( m < 0.188 ) 
    i1 = np.where( 0.188 <= m )
              
    xi =  m[i0] * 0
    xi = np.append( xi, 0.218 * m[i1] ** -2.35 )

    return( xi )


def logplot( ):
    from matplotlib import pyplot as pl
    from scipy import integrate as scint

    

    mall = np.arange( 0.08, 125, 0.01 )


    norm_k93 = scint.simps( kroupa93( mall ), mall )
    norm_ds = scint.simps( dietSalpeter( mall ), mall )

    mwd = np.arange( 3, 8, 0.001 )

    ax = pl.gca()

    floor = np.ones( len(mwd) ) * 1e-6
    k93all = kroupa93( mall ) /  norm_k93 
    dsall = dietSalpeter( mall )  /  norm_ds 

    k93wd = kroupa93( mwd ) /  norm_k93 
    dswd = dietSalpeter( mwd )  /  norm_ds 

    ax.plot( np.log10( mall ) , np.log10( k93all) , 'b--', label='Kroupa+ 93' )
    ax.plot( np.log10( mall ) , np.log10( dsall)  , 'r-' , label='diet Salpeter' )

    ax.fill_between( np.log10( mwd ), np.log10( floor ) , np.log10( dswd ), color='r', alpha=1 )
    ax.fill_between( np.log10( mwd ), np.log10( floor ) , np.log10( k93wd), color='b', alpha=1 )

    ax.set_xlim( -1.5, 2.1 )
    ax.set_ylim( -5, 2 )
    ax.set_xlabel( 'log$_{10}$( M )' )
    ax.set_ylabel( 'log$_{10}$( N )' )

    ax.legend( loc='upper right' )
    ax.text( 0.7, -2, 'CO WDs')


def linplot( ):
    from matplotlib import pyplot as pl

    mall = np.arange( 0.08, 10, 0.01 )
    mwd = np.arange( 3, 8, 0.001 )

    ax = pl.gca()

    floor = np.ones( len(mwd) ) * 1e-6
    k93all = kroupa93( mall )
    dSall = dietSalpeter( mall )

    ax.plot( mall ,  k93all , 'b--' )
    ax.plot( mall ,  dSall  , 'r-'  )

    # ax.fill_between( np.log10( mwd ), np.log10( floor ) , np.log10( dietSalpeter( mwd ) ), color='r', alpha=1 )
    # ax.fill_between( np.log10( mwd ), np.log10( floor ) , np.log10( kroupa93( mwd ) ), color='b', alpha=1 )

    # ax.set_xlim( -1.5, 2.1 )
    # ax.set_ylim( -5, 2 )
    ax.set_xlabel( 'M' )
    ax.set_ylabel( 'N' )




def dtd_kroupa_to_diet_conversion():
    from scipy import integrate as scint


    mall = np.arange( 0.08, 125, 0.001 )
    mwd = np.arange( 3, 8, 0.001 )


    norm_k93 = scint.simps( kroupa93( mall ), mall )
    norm_ds = scint.simps( dietSalpeter( mall ), mall )

    floor = np.ones( len(mwd) ) * 1e-6
    k93all = kroupa93( mall )    / norm_k93 
    dSall = dietSalpeter( mall ) / norm_ds  

    k93wd = kroupa93( mwd )    / norm_k93 
    dSwd = dietSalpeter( mwd ) / norm_ds  

    dtdnorm_k93 = scint.simps( k93wd, mwd ) / scint.simps( k93all * mall, mall )
    dtdnorm_dS = scint.simps( dSwd, mwd ) / scint.simps( dSall * mall, mall )

    print( '%.3f   %.3f    %.3f'% ( dtdnorm_k93, dtdnorm_dS, dtdnorm_dS / dtdnorm_k93 ) )
