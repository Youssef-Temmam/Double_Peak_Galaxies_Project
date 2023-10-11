#libraries
import os
import pandas
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
import sys 

#defining gaussian function
def gaussian(amplitude, fwhm, mean):
    return lambda x: amplitude * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)

#opening redshift file and reading it
#df = pandas.read_csv("redshift_updated.csv", sep = "\t", index_col = "ProcessingID")
df = pandas.read_csv("redshift.csv", sep = "\t", index_col = "ProcessingID")
index_array = df.index.to_numpy()

#parameters
galaxy = '0894-52615-0339' #np.random.choice(index_array) #'spec-0630-52050-0245'#'spec-0405-51816-0035'# spec-0894-52615-0339
line = 'Hbeta'
#line_skewness = df.loc[galaxy, line+'.Skewness']
#line_kurtosis = df.loc[galaxy, line+'.Kurtosis']
line_amplitude = df.loc[galaxy, 'galaxy.linemeas.' + line + '.LinemeasAmplitude']*1e17
line_mean = df.loc[galaxy, 'galaxy.linemeas.' + line + '.LinemeasLineLambda']
line_deviation = df.loc[galaxy, 'galaxy.linemeas.' + line + '.LinemeasLineWidth']
line_continuum = df.loc[galaxy, 'galaxy.linemeas.' + line + '.LinemeasLineCenterContinuumFlux']*1e17
line_fwhm = line_deviation*2.355

#opening galaxy spectrum .fits file
fits_file_path = os.path.join('.', 'spec-' + galaxy + '.fits')
with fits.open(fits_file_path) as hdul:
            spectrum_hdu = hdul[1]
            spectrum_header = spectrum_hdu.header
            spectrum_table = Table(spectrum_hdu.data)
            spectrum_table.add_index('loglam')

            #converting log wavelength 
            waves = np.array([ 10**x for x in spectrum_table['loglam']])
            spectrum_table['loglam'] = waves
            flux = spectrum_table['flux']
            

            #cutting line spectrum at the interval [lambda - 10 sigma, lambda + 10 sigma]
            line_spectrum_table = spectrum_table.loc[line_mean - 10*line_deviation: line_mean + 10*line_deviation]
            line_fluxes = line_spectrum_table['flux']
            line_waves = line_spectrum_table['loglam']

fitted_line_fluxes = gaussian(line_amplitude, line_fwhm, line_mean)(line_waves) + line_continuum

#plotting line spectra with fit
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
ax.plot(line_waves, line_fluxes, label='Raw flux', color='black', linewidth=1.5)
ax.plot(line_waves, fitted_line_fluxes, label = 'Fitted flux', color='purple', linewidth=2.)
plt.ylabel("Flux" + r"$(erg/s/cm^2/\AA)$")
plt.xlabel("Wavelength" + r"$(\AA)$")
plt.title("Raw and fitted fluxes for H" + r"$\beta$" + " emission line in terms of wavelength for galaxy:spec-" + galaxy)
#line_skewness_str = str(line_skewness)
#line_kurtosis_str = str(line_kurtosis)
#plt.plot([], [], ' ', label='Skewness =' + line_skewness_str)
#plt.plot([], [], ' ', label='Kurtosis =' + line_kurtosis_str)
plt.legend()
plt.show()

