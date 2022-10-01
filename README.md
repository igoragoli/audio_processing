# Répresentation Parcimonieuses des Signaux

## `filter_banks/`
The programs in this folder consist in using filter banks and the STFT(Short-Term Fourier Transform) to separate 4 sinusoidal signals.
- `filter_banks.m`
- `parameters.m` (pour la configuration des paramètres pour la TFCT)
- `get_frequencies_fft.m` 
- `filter_banks.ipynb`
- `filter_banks.ipynb`

## `wavelets/`
Denoising and signal separation is performed by using wavelet transforms. The package `ltfat-master` is used in this part.
- `denoising.m`
- `separation.m`

## `supervised_learning``
- `generate_base.m`
- `main.m`
- `train_tl3.mat` 
- `test_tl3.mat`
- `val_tl3.mat`
- `pr_logo3.m`
- `uncontretousoutil.m`
- `lagis_rbf_gaussien.m`
