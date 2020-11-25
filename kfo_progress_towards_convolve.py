def csfr(z):
    # behroozi et al. 2013; table 6 https://iopscience.iop.org/article/10.1088/0004-637X/770/1/57
    z0 = 1.243
    A = -0.997
    B = 0.241
    C = 0.180 # C ~ [Msol yr^-1 comoving Mpc^-3]
    return C/(10**(A*(z-z0)) + 10**(B*(z-z0)))

z = np.linspace(0,2.5)

fig,ax = plt.subplots(figsize=(16,8))
ax.plot(z,[csfr(i) for i in z],label="Behroozi 13")
ax.set_xlabel("Redshift")
ax.set_ylabel("CSFH [$M_{\odot} yr^{-1} Mpc^{-3}$]")
ax.legend()

promptDTD = ascii.read(open("PromptDTD_Rodney2014.csv","rb"))
print(promptDTD.meta)
Hubble = promptDTD[0] # uses clash+candels rates
Ground = promptDTD[1] # uses various ground based dominated by z < 1
All = promptDTD[2] # combines all data
def promptDTD_Params(table):
    eta = table['$\eta_{Ia}$'] # 10^-4 Msol^-1 yr^-1
    NtoM_assumeDTD = table['$N_{Ia}/M_*^{b,c}$'] # 10^-3 Msol^-1 
    NtoM_noDTD = table['$N_{Ia}/M_*^{c,d}$'] # 10^-3 Msol^-1 
    fp = table['f_P']
    return eta,fp
HubbleParams = promptDTD_Params(Hubble)
HubbleEta,HubbleFp = float(HubbleParams[0].split('_')[0]),float(HubbleParams[1].split('_')[0])
GroundParams = promptDTD_Params(Ground)
GroundEta,GroundFp = float(GroundParams[0].split('_')[0]),float(GroundParams[1].split('_')[0])
def SNR(t,eta,fp):
    # Rodney 2014 prompt model to Ia delay times, breaks model at time 500 Myr, 40 Myr for CO WD to form
    # t ~ Gyr 
    # eta ~ SN Ia yr^-1 Msol^-1
    K = 7.132 # defined by the times defining the model
    if t <= 0.04:
        snr = 0
    elif 0.5 > t > 0.04:
        snr = K*eta*(fp/(1-fp))
    elif t >= 0.5:
        snr = eta*t**(-1)
    return snr

t = np.linspace(0.01,13,1000) # Gyr
z = [z_at_value(cosmo.lookback_time, i * u.Gyr) for i in t]

HubbleSNR = [SNR(i,HubbleEta*10**(-4),HubbleFp) for i in t]
GroundSNR = [SNR(i,GroundEta*10**(-4),GroundFp) for i in t]

fig,ax = plt.subplots(figsize=(16,8))
ax.plot(t,HubbleSNR,label="Hubble")
ax.plot(t,GroundSNR,label="Ground")
ax.set_xlabel("t [$Gyr$]")
ax.set_ylabel("SN Ia [$yr^{-1} M_{\odot}^{-1}$]")
ax.set_xlim(0,4)
ax.set_yscale('log')
ax.legend()

#z = np.linspace(0,20,1000)

# the DTD for times from hypothetical burst up to 13 Gyr
t = np.linspace(0.01,13,1000) # Gyr 
HubbleSNR = [SNR(i,HubbleEta*10**(-4),HubbleFp) for i in t]

z = [z_at_value(cosmo.lookback_time, i * u.Gyr) for i in t]
#z = np.linspace(0,20,1000)
#z = z[::-1] # want to convolve the delay times over starting point of star formation at high-z 
CSFH = [csfr(i) for i in z]

tmp = np.convolve(HubbleSNR,CSFH)

tmp = np.convolve(HubbleSNR,CSFH,'same')
tmp2 = np.convolve(CSFH,HubbleSNR,'same')

fig,ax = plt.subplots(figsize=(16,8))
ax.plot(z,tmp/10**(-4))
ax.set_ylabel(r"RIa [$ \times 10^{-4} yr^{-1} Mpc^{-3} h^3_{70}$]")
ax.set_xlim(0,2.5)