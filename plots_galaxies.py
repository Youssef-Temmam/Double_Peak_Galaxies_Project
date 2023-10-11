#libraries 
import numpy as np 
import math  
import matplotlib.pyplot as plt  
import scipy.constants as cte  
from mpl_toolkits import mplot3d  
from sympy import *  
from matplotlib import cm  
from scipy.integrate import quad  
import matplotlib.colors as colors  
import matplotlib.cbook as cbook  
from matplotlib import cm  
from sklearn.linear_model import LinearRegression  
import random
import pandas
import sys

#reading the file    
df = pandas.read_csv("./redshift.csv", sep = "\t", index_col = "ProcessingID")
initial_redshift = df["galaxy.Redshift"].values

#creating new data frame with needed columns
#selected_columns = ["galaxy.Redshift", "galaxy.linemeas.Hbeta.LinemeasLineFluxDirectIntegration",
#                    "galaxy.linemeas.Hbeta.LinemeasLineFlux","galaxy.linemeas.Hbeta.LinemeasOffset", 
#                    "galaxy.linemeas.Hbeta.LinemeasVelocity", "galaxy.linemeas.Hbeta.LinemeasAmplitude",
#                    "galaxy.linemeas.Hbeta.LinemeasLineWidth"]

line = "Hbeta"

selected_columns = ["galaxy.Redshift", "galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration",
                    "galaxy.linemeas." + line + ".LinemeasLineFlux","galaxy.linemeas." + line + ".LinemeasOffset", 
                    "galaxy.linemeas." + line + ".LinemeasVelocity", "galaxy.linemeas." + line + ".LinemeasAmplitude",
                    "galaxy.linemeas." + line + ".LinemeasLineWidth"]

new_df = df[selected_columns]
#cleaning data frame from nan values
new_df = new_df.dropna()
#cleaned_df = new_df
#removing null and negative values in fluxes and amplitudes
cleaned_df = new_df[(new_df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] > 0) 
                    & (new_df["galaxy.linemeas." + line + ".LinemeasLineFlux"] > 0)
                    & (new_df["galaxy.linemeas." + line + ".LinemeasAmplitude"] > 0)]

galaxies_id = cleaned_df.index.to_numpy()
redshift = cleaned_df["galaxy.Redshift"].values

#cleaned_df['Flux.relative.difference'] = (cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] 
#                                          - cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFlux"])/cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] 


#Hbeta
fluxdirectInteg = cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"].values*1e17
fluxfitted = cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFlux"].values*1e17
offset_column = cleaned_df["galaxy.linemeas." + line + ".LinemeasOffset"].values
velocity = cleaned_df["galaxy.linemeas." + line + ".LinemeasVelocity"].values
amplitude = cleaned_df["galaxy.linemeas." + line + ".LinemeasAmplitude"].values*1e17
fwhm = cleaned_df["galaxy.linemeas." + line + ".LinemeasLineWidth"].values
deviation = fwhm#/2.355

#flux_diff = cleaned_df['Flux.relative.difference'].values
#print("negative galaxies:", len(cleaned_df[(cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] <= 0) 
#                    & (cleaned_df["galaxy.linemeas." + line + ".LinemeasLineFlux"] <= 0)]))

#arrays parameters
flux_difference = np.abs(fluxfitted - fluxdirectInteg)/fluxdirectInteg
flux_difference_s = (fluxfitted - fluxdirectInteg)/fluxdirectInteg
flux_difference_er = flux_difference[flux_difference != 0]
flux_difference_above10pct_ind = np.where(flux_difference > 0.1)
num_gal_under10pct = np.sum(flux_difference_er<=0.1)

area = amplitude*deviation*np.sqrt(2*np.pi)

max_value_index = np.argmax(flux_difference)
print(flux_difference[max_value_index])
print(cleaned_df.iloc[max_value_index, :])

#plotting fitted flux in terms of direct integration flux
plt.figure(figsize=(12,6))
plt.plot(fluxdirectInteg, fluxfitted, "+")
plt.plot([0, 1400], [0, 1400], color="red")
plt.xlabel("Direc integrated flux"+r"$(erg/s/cm^2)$")
plt.ylabel("Fitted flux"+r"$(erg/s/cm^2)$")
plt.title("Fitted flux in terms of Direc integrated flux in H" + r"$\beta$"  + " for the 1000 galaxies sample")
plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
plt.show()
plt.close()

sys.exit()

#plotting module of flux difference in terms of velocity
plt.figure(figsize=(12,7))
plt.plot(velocity, flux_difference, "+")
plt.ylabel("Module of flux relative difference")
plt.xlabel("Velocity"+r"$(km/s)$")
plt.title("Module of flux relative difference in terms of velocity in H" + r"$\beta$"  + " for the 1000 galaxies sample")
plt.show()
plt.close()

#sys.exit()

#plotting module of flux difference in terms of redshift
plt.figure(figsize=(12,7))
plt.plot(redshift, flux_difference, "+")
plt.ylabel("Module of flux relative difference")
plt.xlabel("Redshift")
plt.title("Module of flux relative difference in H" + r"$\beta$"  + " in terms of redshift for the 1000 galaxies sample")
plt.show()
plt.close()

#plotting module of flux difference in terms of Direct integrated flux
plt.figure(figsize=(12,7))
plt.plot(fluxdirectInteg, flux_difference, "+")
plt.ylabel("Module of flux relative difference")
plt.xlabel("Direct integrated flux"+r"$(erg/s/cm^2)$")
plt.title("Module of flux relative difference in terms of Direct integrated flux in H" + r"$\beta$"  + " for the 1000 galaxies sample")
plt.show()
plt.close()

#plotting Fitted flux in terms of the surface area under the guassian
plt.figure(figsize=(12,7))
plt.plot(area, fluxfitted, "+")
plt.plot([0, 1400], [0, 1400], color="red")
plt.ylabel("Fitted flux"+r"$(erg/s/cm^2)$")
plt.xlabel("Flux computed from the area under the gaussian"+r"$(erg/s/cm^2)$")
plt.title("Fitted flux in terms of flux computed from the area under the gaussian in H" + r"$\beta$"  + " for the 1000 galaxies sample")
plt.show()
plt.close()

#printing useful information about amazed exactitude
print("number of galaxies with calculated redshift:",1000 - np.isnan(initial_redshift).sum(), "/1000")
print("number of galaxies with fitted Halpa line:", len(fluxfitted), "/1000")
print("number of fitted Halpha galaxies with fit error under 10% :", num_gal_under10pct, "/1000")
#print("galaxies with fitted Halpha fit error above 10% :", [galaxies_id[ind] for ind in flux_difference_above10pct_ind ])

