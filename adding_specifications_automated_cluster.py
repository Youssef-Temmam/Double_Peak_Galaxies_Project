#libraries
import csv
import os
import pandas
import numpy as np
from astropy.io import fits
from astropy.table import Table

#defining function that adds new columns to the file and close it
def add_column_to_csv(input_file, new_column_name, new_column_values):
    #creating a temporary file
    temp_file = 'temp.csv'
    with open(input_file, 'r') as input_csv, open(temp_file, 'w', newline = '') as temp_csv:
        #loading data from input fi le
        reader = csv.DictReader(input_csv, delimiter='\t')
        fieldnames = reader.fieldnames + [new_column_name]

        #writting in the temp file
        writer = csv.DictWriter(temp_file, fieldnames=fieldnames)
        for row, new_value in zip(reader, new_column_values):
            row[new_column_name] = new_value
            writer.writerow(row)

    #replacing input file with temp file
    os.replace(temp_file, input_file)

#defining SNR function
def snr(flux, continuum):
    return flux/continuum

#defining second degree polynom function
def degree_2_polynom(coefficient0, coefficient1, coefficient2):
    return lambda x: coefficient2*x**2 + coefficient1*x + coefficient0

#defining mean of a signal as an array
def mean_line(array):
    n = len(array)
    return np.sum(array)/n

#defining variance of a signal as an array
def variance_line(array):
    n = len(array)
    return np.sum((array - mean_line(array))**2)/n

#defining skewness of a signal as an array
def skewness_line(array):
    n = len(array)
    return np.sum((array - mean_line(array))**3)/(n*(variance_line(array)**(3/2)))

#defining kurtosis of a signal as an array
def kurtosis_line(array):
    n = len(array)
    return np.sum((array - mean_line(array))**4)/(n*(variance_line(array)**2))

line_names  = np.array(['Halpha','Hbeta', 'Hgamma', '[OIII](doublet-1/3)', '[OIII](doublet-1)', 
                   '[NII](doublet-1)', '[SII]6716', '[SII]6731'])

#opening redshift file and reading it
df = pandas.read_csv("output_100_000_DP_galaxies/redshift.csv", sep = "\t", index_col = "ProcessingID")

#cleaning data frame from nan values in redshift
df.dropna(subset = ['galaxy.Redshift'])

#cleaning data frame from negative fluxes and fits 
#for line in line_names:
#    df = df[(df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] > 0) 
#                    & (df["galaxy.linemeas." + line + ".LinemeasLineFlux"] > 0)]


for line in line_names:
    #reading data frame for line
    line_wavelengths = df['galaxy.linemeas.' + line + '.LinemeasLineLambda'].values
    line_fitted_flux = df['galaxy.linemeas.' + line + '.LinemeasLineFlux'].values*1e17
    line_fitted_continuum = df['galaxy.linemeas.' + line + '.LinemeasLineCenterContinuumFlux'].values*1e17
    line_deviations = df['galaxy.linemeas.' + line + '.LinemeasLineWidth'].values
    line_continuum_coefficient0 = df['galaxy.linemeas.' + line + '.LinemeasContinuumPCoeff0'].values*1e17
    line_continuum_coefficient1 = df['galaxy.linemeas.' + line + '.LinemeasContinuumPCoeff1'].values*1e17
    line_continuum_coefficient2 = df['galaxy.linemeas.' + line + '.LinemeasContinuumPCoeff2'].values*1e17

    #computing SNR and adding it to the data frame
    df[line + '.SNR'] = snr(line_fitted_flux, line_fitted_continuum)

    #computing skewneess and kurtosis 
    line_Skewness = np.zeros(len(df))*np.nan
    line_Kurtosis = np.zeros(len(df))*np.nan

    for galaxy, index in zip(df.index, range(len(df))):
        #taking deviation and wavelength 
        sigma_line = line_deviations[index]
        lambda_line = line_wavelengths[index]
        coeff0_line = line_continuum_coefficient0[index]
        coeff1_line = line_continuum_coefficient1[index]
        coeff2_line = line_continuum_coefficient2[index]

        #print('before opening', '>>gaalxy ID:', galaxy, '>>line:', line, '>>deviation:', sigma_line, '>>wavelength:', lambda_line)
        
        #computing Skewness and Kurtosis only for fitted lines
        if np.isnan(lambda_line) == False and lambda_line < 9000:
            #opening galaxy spectrum .fits file
            fits_file_path = os.path.join('/net/CESAM/amazed/dataset/SDSS/dr17/sdss/spectro/redux/26/spectra/lite/' + galaxy[: 4], 'spec-' + galaxy + '.fits')
            with fits.open(fits_file_path) as hdul:
                spectrum_hdu = hdul[1]
                spectrum_header = spectrum_hdu.header
                spectrum_table = Table(spectrum_hdu.data)
                #converting wavelength from log to normal
                waves = np.array([ 10**x for x in spectrum_table['loglam']])
                spectrum_table['loglam'] = waves
                spectrum_table.add_index('loglam')

                #print('>>line:', line, '>>deviation', sigma_line, '>>wavelength', lambda_line)
                #cutting the spectrum only around 2.5 deviation of the lines, we can go up to 3
                if line == 'Halpha':
                    spectrum_table_line = spectrum_table.loc[lambda_line - 2*sigma_line: lambda_line + 2*sigma_line]
                if line == '[NII](doublet-1)':
                    spectrum_table_line = spectrum_table.loc[lambda_line - 2*sigma_line: lambda_line + 2*sigma_line]
                #if line == '[NII](doublet-1/3)':
                #    spectrum_table_line = spectrum_table.loc[lambda_line - 2.5*sigma_line: lambda_line + 2.5*sigma_line]
                if line == '[SII]6716]':
                    spectrum_table_line = spectrum_table.loc[lambda_line - 2.5*sigma_line: lambda_line + 2.5*sigma_line]
                if line == '[SII]6731':
                    spectrum_table_line = spectrum_table.loc[lambda_line - 2.5*sigma_line: lambda_line + 2.5*sigma_line]                   
                else:
                    spectrum_table_line = spectrum_table.loc[lambda_line - 4*sigma_line: lambda_line + 4*sigma_line]
                line_spectra = spectrum_table_line['flux']
                line_continuum = degree_2_polynom(coeff0_line, coeff1_line, coeff2_line)(spectrum_table_line['loglam'])

                #computing Skewness and Kurtosis
                line_Skewness[index] = skewness_line(line_spectra - line_continuum)
                line_Kurtosis[index] = kurtosis_line(line_spectra - line_continuum)

    #computing Skewness and Kurtosis in the data frame
    df[line + '.Skewness'] = line_Skewness
    df[line + '.Kurtosis'] = line_Kurtosis

#saving data frame in a new file
df.to_csv('redshift_updated.csv', sep = '\t')
#print(df)
