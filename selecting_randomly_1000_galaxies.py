#libraries
import numpy as np
import pandas
import os
import sys

#making list of cluster SDSS spectra
plates_directory = "/net/CESAM/amazed/dataset/SDSS/dr17/sdss/spectro/redux/26/spectra/lite"
#array of plates ids
plates_id = np.array(os.listdir(plates_directory))

spectra_id_list = []
for plate in plates_id:
    spectra_id_list += os.listdir(plates_directory + "/" + plate)

#array of SDSS spectra
spectra_SDSS = np.array(spectra_id_list)

#choosing 1000 random galaxies
chosen_galaxies_number = 1000
array_1000 = np.random.choice(spectra_SDSS, chosen_galaxies_number)

#making the input spectra file for amazed
df_1000 = pandas.DataFrame(array_1000)
df_1000.to_csv("nom_1000_galaxies.txt", header=None, index=None)
