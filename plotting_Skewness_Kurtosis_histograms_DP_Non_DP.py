#libraries
import os
import pandas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sys

#opening redshift file and reading it
df_DP = pandas.read_csv("redshift_updated_DP.csv", sep = "\t", index_col = "ProcessingID")
df_Non_DP = pandas.read_csv("redshift_updated_Non_DP.csv", sep = "\t", index_col = "ProcessingID")

line_names  = np.array(['Halpha','Hbeta', 'Hgamma', '[OIII](doublet-1/3)', '[OIII](doublet-1)',
                   '[NII](doublet-1)', '[SII]6716', '[SII]6731'])

#cleaning nan, negative fluxes and snr from data frame:
for line in line_names:
    df_DP = df_DP.dropna(subset = ["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration",
                                   "galaxy.linemeas." + line + ".LinemeasLineFlux", line + '.SNR'])
    df_DP = df_DP[(df_DP["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] > 0)
                    & (df_DP["galaxy.linemeas." + line + ".LinemeasLineFlux"] > 0)
                    & (df_DP[line + '.SNR'] > 0)]
    df_Non_DP = df_Non_DP.dropna(subset = ["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration",
                                   "galaxy.linemeas." + line + ".LinemeasLineFlux", line + '.SNR'])
    df_Non_DP = df_Non_DP[(df_Non_DP["galaxy.linemeas." + line + ".LinemeasLineFluxDirectIntegration"] > 0)
                    & (df_Non_DP["galaxy.linemeas." + line + ".LinemeasLineFlux"] > 0)
                    & (df_Non_DP[line + '.SNR'] > 0)]

n_bins = 100

line = 'Hbeta'
line_Skewness_DP = df_DP[line + '.Skewness']
line_Kurtosis_DP = df_DP[line + '.Kurtosis']
line_Skewness_Non_DP = df_Non_DP[line + '.Skewness']
line_Kurtosis_Non_DP = df_Non_DP[line + '.Kurtosis']

#plottiNon_ng histograms for Skewness for line
fig, axs = plt.subplots(1, 1, figsize =(6, 7), tight_layout = True)
axs.hist(line_Skewness_DP, bins = n_bins, density = True, alpha=0.5, color = 'red', label = 'DP galaxies sample')
axs.hist(line_Skewness_Non_DP, bins = n_bins, density = True, alpha=0.5, color = 'blue', label = 'Non DP galaxies sample')
axs.set_xlim(-2, 2)
plt.xlabel('Skewness at H' + r"$\beta$")
plt.ylabel('Frequency')
plt.title('Frequency of Galaxies in terms of their Skewness at H' + r"$\beta$")
plt.legend()
plt.show()
plt.close()


#plotting histograms for kurtosis for line
fig, axs = plt.subplots(1, 1, figsize =(6, 7), tight_layout = True)
axs.hist(line_Kurtosis_DP, bins = n_bins, density = True, alpha=0.5, color = 'red', label = 'DP galaxies sample')
axs.hist(line_Kurtosis_Non_DP, bins = n_bins, density = True, alpha=0.5, color = 'blue', label = 'Non DP galaxies sample')
axs.set_xlim(-4, 4)
plt.xlabel('Kurtosis at H' + r"$\beta$")
plt.ylabel('Frequency')
plt.title('Frequency of Galaxies in terms of their Kurtosis at H' + r"$\beta$")
plt.legend()
plt.show()
plt.close()

#sys.exit()

line = 'Halpha'
line_Skewness_DP = df_DP[line + '.Skewness']
line_Kurtosis_DP = df_DP[line + '.Kurtosis']
line_Skewness_Non_DP = df_Non_DP[line + '.Skewness']
line_Kurtosis_Non_DP = df_Non_DP[line + '.Kurtosis']

print(len(df_DP))
print(len(df_Non_DP))
print('Skewness mean for DP:', np.mean(line_Skewness_DP))
print('Skewness mean for Non DP:', np.mean(line_Skewness_Non_DP))

#plottiNon_ng histograms for Skewness for line
fig, axs = plt.subplots(1, 1, figsize =(6, 7), tight_layout = True)
axs.hist(line_Skewness_DP, bins = 100, density = True, alpha=0.5, color = 'red', label = 'DP galaxies sample')
axs.hist(line_Skewness_Non_DP, bins = 100, density = True, alpha=0.5, color = 'blue', label = 'Non DP galaxies sample')
axs.set_xlim(-2, 2)
plt.xlabel('Skewness at H' + r"$\alpha$")
plt.ylabel('Frequency')
plt.title('Frequency of Galaxies in terms of their Skewness at H' + r"$\alpha$")
plt.legend()
plt.show()
plt.close()


#plotting histograms for kurtosis for line
fig, axs = plt.subplots(1, 1, figsize =(6, 7), tight_layout = True)
axs.hist(line_Kurtosis_DP, bins = 100, density = True, alpha=0.5, color= 'red', label = 'DP galaxies sample')
axs.hist(line_Kurtosis_Non_DP, bins = 100, density = True, alpha=0.5, color = 'blue', label = 'Non DP galaxies sample')
axs.set_xlim(0, 6)
plt.xlabel('Kurtosis at H' + r"$\alpha$")
plt.ylabel('Frequency')
plt.title('Frequency of Galaxies in terms of their Kurtosis at H' + r"$\alpha$")
plt.legend()
plt.show()
plt.close()
