import numpy as np
import matplotlib.pyplot as plt
from astroquery.ned import Ned
from astroquery.irsa import Irsa
import astropy.units as u
snname,hostnames = np.loadtxt('selected_sne.dat',usecols=(0,2),unpack=True,delimiter=',',dtype='str')


unique_hosts = np.unique(hostnames)
has_spectra = 0
'''for i,unique_host in enumerate(unique_hosts[0:110]):
    print i
    result_table = Ned.get_image_list(unique_host,item="spectra")
    if len(result_table) > 0:
        has_spectra = has_spectra + 1
print has_spectra
'''
'''
has_photometry = 0
for i,unique_host in enumerate(unique_hosts):
    print i
    try:
        result_table = Ned.get_table(unique_host,table="photometry")
    except:
        continue
    if len(result_table) > 0:
        has_photometry = has_photometry + 1

print has_photometry
'''

# sne	w1mag	w1sigmn	w2mag	w2sigm	w3mag	w3sigm	w4mag	w4sigm	freq1, freq2, freq3, ...... freqn

hasboth = 0


for i,unique_host in enumerate(unique_hosts):
    print i
    isWise = 0
    isUV = 0
    try:
        photometry = Ned.get_table(unique_host,table="photometry")
        freq = photometry["Frequency"]
        isUV = np.where(freq>1280000000000000.0)
        if isUV > 0:
            freq_u = np.unique(freq)
            mag = photometry["Photometry Measurement"]
            sigs = photometry["Uncertainty"]
            signew = np.zeros((len(sigs),1))
            neddata = np.zeros((len(freq_u),3))
            neddata[:,0] = freq_u

            for k, sig in enumerate(sigs[:-2]):
                sigval = sig[3:]
                signew[k] = 0
                if len(sigval) > 0:
                    signew[k] = float(sig[3:])
            for j,freq_u in enumerate(freq_u):
                ind = np.where(freq==freq_u)
                magmed = np.median(mag[ind])
                sigmed = np.median(signew[ind])
                neddata[j,1] = magmed
                neddata[j,2] = sigmed
    except:
        continue


    try:
        wise = Irsa.query_region(unique_host,catalog="allsky_4band_p3as_psd", spatial="Cone",radius=5 * u.arcsec)

        wisedat = np.zeros((4,3))
        wisedat[:,0] = np.asarray([1,2,3,4])
        wisedat[0,1] = float(wise["w1m"][0])
        wisedat[1,1] = float(wise["w2m"][0])
        wisedat[2,1] = float(wise["w3m"][0])
        wisedat[3,1] = float(wise["w4m"][0])

        if float(wise["w1sigm"][0]) > 0:
            wisedat[0,2] = float(wise["w1sigm"][0])
        if float(wise["w2sigm"][0]) > 0:
            wisedat[1,2] = float(wise["w2sigm"][0])
        if float(wise["w3sigm"][0]) > 0:
            wisedat[2,2] = float(wise["w3sigm"][0])
        if float(wise["w4sigm"][0]) > 0:
            wisedat[3,2] = float(wise["w4sigm"][0])
        wisedat = wisedat[~np.isnan(wisedat).any(axis=1)]
        isWise = 1
    except:
        continue
 
    if (isWise>0) & (isUV>0):
        print "whee, both"
        combined = np.vstack((neddata,wisedat))
        hasboth = hasboth+1
        #print combined
        #print data
        #print wisedat

print hasboth
    
'''
print has_uv,has_ir,has_both

blah = Irsa.query_region(unique_hosts[0],catalog="allsky_4band_p3as_psd",
                         spatial="Cone",radius=5 * u.arcsec)
print blah['w3mag'],blah['w3sigm']
'''

