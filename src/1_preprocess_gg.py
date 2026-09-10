import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import pygam as gm
import copy as cp

import os
from pathlib import Path

database = 'GG'
dirdatain = './data/'
dirdataout = './output/processed_data/'
# create directories if they dont exist
out_dir = Path(f"./{dirdataout}/GAM_species/{database}")
out_dir.mkdir(parents=True, exist_ok=True)

species = {} # dictionary with the IDs of the species (as they change across networks)

mydf = pd.read_csv(dirdatain+f'/{database}/bale_pollen.csv')
mydf = mydf.fillna(0)
mydf = mydf.rename(columns={'kaBP':'age'})
mydf['age'] = -mydf['age']
mydf = mydf.sort_values(by=['age'])

myspecies = 'Alchemilla,Chenopodiaceae,Artemisia,Cichoroideae,Asteroideae,Carduus,Asteraceae,Vernonia,Polygonum,Tribulus,Galium,Brassicaceae,Caryophylaceae,Swertia,Senecio,Erica,Hypericum,Anthospermum,Cerastium,Apiaceae,Podocarpus,Euclea,Hagenia,Juniperus,Myrsine,Rapanea,Rosa,Schefflera,Brucea,Sapotaceae,Iridaceae,Anthyllis,Maytenus,Buxus,Psydrax,Celastraceae,Celtis,Dodonaea,Ekebergia,Ephedra,Clematis,Myrica,Hypoestes,Ziziphus,Acanthus,Allophylus,Acacia,Combretum,Commiphora,Cussonia,Dobera,Lannea,Macaranga,Olea,Jasminum,Rhus,Securinega,Teclea,Zanthoxylum,Cassia,Capparidaceae,Acalypha,Euphorbia,Euphorbiacea,Phyllantus,Tamarindus,Alchornea,Aloe,Asphodellus,Leguminosae,Rubiaceae,Blepharis,Justicia,Heliotropium,Kohautia,Indigofera,Lamiaceae,Cerealia,Plantago,Planceolata,Poaceae,Rumex,Solanum,Urticaceae,Ricinus'.replace(' ','_').split(',')

mydf_GG_original = pd.read_csv(dirdatain + f'/{database}/bale_pollen.csv')

mydf_GG_original['kaBP'] = - mydf_GG_original['kaBP']
mydf_GG_original = mydf_GG_original.fillna(0)
mydf_GG_original = mydf_GG_original.sort_values(by=['kaBP'])

all_species = list(mydf_GG_original.columns)
metrics_list = ['accrate', 'volume', 'Lycopodium', 'Lyc.conc']

mydf_GG_original['kaBP'] = mydf_GG_original['kaBP'] * 1000  # Convert from ka to years
mydf_GG_original.rename(columns={'kaBP': 'age'}, inplace=True)


multipliers = mydf_GG_original['Lyc.conc'] * mydf_GG_original['accrate'] / (mydf_GG_original['volume'] * mydf_GG_original['Lycopodium'])

mydf_pollen = mydf_GG_original[['age'] + myspecies]
mydf_pollen = mydf_pollen.set_index('age')

for spec in myspecies[:]: # ITERATE OVER A COPY
    if mydf_pollen[spec].sum() == 0:
        print('Species %s has all zeros!' %(spec))
        myspecies.remove(spec)

mydf = mydf_pollen[myspecies]

n = 300
fit_type = f'fit_{n}_0.01'
species = {}
vx = gm.utils.make_2d(mydf.index, verbose=False).astype('float')
for spec in myspecies:
    print('\tIterating on species %s ...' %(spec))
    vy = mydf[spec]

    mygam = gm.GAM(terms='auto', n_splines=n+1, lam=0.01).fit(vx, vy)
    # mygam = gm.GAM(terms='auto').gridsearch(vx, vy)

    species[spec] = {}
    species[spec]['prev_x'] = mydf.index.to_numpy()
    species[spec]['prev_y'] = mydf[spec].to_numpy()

    XX = mygam.generate_X_grid(term=0, n=n)
    YY = mygam.predict(XX)
    species[spec]['x'] = cp.copy(XX)[:,0]
    species[spec]['y'] = cp.copy(YY)

species_df = pd.DataFrame(species)
species_df.to_pickle(f'{dirdataout}/GAM_species/{database}/species_{fit_type}.pkl')

species_index = {}
for i, spec in enumerate(myspecies):
    species_index[spec] = i
pd.to_pickle(species_index, f'{dirdataout}/GAM_species/{database}/species_index.pkl')

# GAM over multipliers
n_splines = n
fit_type = 'fit_%d_0.01' %(n_splines)
vx = gm.utils.make_2d(multipliers.index, verbose=False).astype('float')
vy = multipliers.values

gam_multipliers = gm.GAM(terms='auto', n_splines=n_splines, lam=0.01).fit(vx, vy)

XX = gam_multipliers.generate_X_grid(term=0, n=n_splines)
YY = gam_multipliers.predict(XX)

gam_multipliers_df = pd.DataFrame({
    'cal BP': XX[:, 0],
    'multipliers': YY
})

gam_multipliers_df.to_pickle(f'{dirdataout}/GAM_species/{database}/gam_multipliers_{fit_type}.pkl')