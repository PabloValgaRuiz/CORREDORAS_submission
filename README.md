# CORREDORAS_submission

Code for reproducibility of the manuscript *Fire modulates millennial plant community assembly*.
Pablo Valgañon-Ruiz, Alessio Cardillo, Ana Cano-Herranz, Hugo Saiz,
Adam T. Clark, David Garcia-Callejas, Lena Neuenkamp, Cristina
Ramos-Capon, Penelope Gonzalez-Samperiz, and Graciela Gil-Romera
## HOW TO RUN

It is recommended to run using conda (installed through miniconda). Open a terminal on the repository folder.

```
conda create -p ./.conda python=3.11
conda activate ./.conda
pip install -r requirements.txt
```

From here, do the processing of the data of both sites.

```
python src/1_preprocess_basa.py
python src/1_preprocess_gg.py
```
From there, open the *2_stationary_bootstrap_GC.py* file from the folder *src/* and choose whether to generate the results of 'basa' (Basa de la Mora) or 'GG' (Garba Guracha).

Run the simulations for the following windows in 'basa':
```
python src/2_stationary_bootstrap_GC.py 0 51
python src/2_stationary_bootstrap_GC.py 51 85
python src/2_stationary_bootstrap_GC.py 85 140
```
Or the following in 'GG':
```
python src/2_stationary_bootstrap_GC.py 0 57
python src/2_stationary_bootstrap_GC.py 57 150
python src/2_stationary_bootstrap_GC.py 150 250
python src/2_stationary_bootstrap_GC.py 250 300
```
These take a long time to finish. There are 4000 iterations of the bootstrap, that can be lowered in the file. You can start plotting the results even before all the bootstraps finish.

For the plots, additional packages must be installed. If using a jupyter notebook (in vscode, for example), install ipykernel.
```
conda install -p ./.conda ipykernel --update-deps --force-reinstall
```
Then additional packages.
```
pip install networkx leidenalg iplotx
```

Then run the 3_plot_results_BS_GC.ipynb cells in the *src/* folder.
