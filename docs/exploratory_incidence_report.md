# Exploratory: incidence-detrended correlation (NOT pre-registered)

Held-out tiles kept by the registered evaluation (85), class-pure 90 m blocks as in R1-R3. Local incidence from DEM normals and the zero-Doppler look vector. Trend = mean dB in 2 deg local-incidence bins (at least 30 blocks), interpolated; blocks outside the binned range are dropped.

## ascending

Ellipsoid incidence over the tiles: 38.2-45.7 deg.

| Class | Blocks | r (raw) | r (incidence removed) | Variance explained by incidence: real / sim | Slope, dB per 10 deg (20-60 deg): real / sim |
|---|---|---|---|---|---|
| tree cover | 319096 | 0.871 | 0.443 | 0.78 / 0.81 | -1.97 / -3.11 |
| grassland | 264940 | 0.860 | 0.564 | 0.67 / 0.83 | -2.21 / -4.83 |
| cropland | 5975 | 0.277 | 0.146 | 0.07 / 0.79 | -1.59 / -4.53 |
| built-up | 8628 | 0.240 | 0.180 | 0.04 / 0.62 | -1.74 / -1.93 |
| bare / sparse vegetation | 53131 | 0.876 | 0.542 | 0.73 / 0.80 | -2.68 / -4.05 |
| permanent water bodies | 19445 | 0.707 | 0.620 | 0.22 / 0.32 | -5.87 / -6.76 |
| moss and lichen | 10022 | 0.900 | 0.481 | 0.79 / 0.88 | -2.59 / -3.81 |

Median dB per 10 deg of local incidence (real / sim):

| Class | 0-10 | 10-20 | 20-30 | 30-40 | 40-50 | 50-60 | 60-70 | 70-80 | 80-90 |
|---|---|---|---|---|---|---|---|---|---|
| tree cover | 1.9 / 8.0 | -0.1 / 4.4 | -2.3 / 1.3 | -4.4 / -1.4 | -6.3 / -4.2 | -8.1 / -7.4 | -10.0 / -10.7 | -12.4 / -12.1 | -15.7 / -13.5 |
| grassland | 2.3 / 13.3 | -1.0 / 8.8 | -4.9 / 2.8 | -7.8 / -4.1 | -9.7 / -9.7 | -11.1 / -11.9 | -12.6 / -13.2 | -14.7 / -14.4 | -18.1 / -17.1 |
| cropland |  |  |  | -8.1 / -7.4 | -8.8 / -9.0 |  |  |  |  |
| built-up |  |  | -2.4 / -0.5 | -3.1 / -3.0 | -3.4 / -3.8 | -6.9 / -5.7 |  |  |  |
| bare / sparse vegetation | 5.7 / 12.7 | 3.6 / 9.4 | 0.7 / 4.6 | -2.0 / -0.6 | -4.7 / -5.1 | -7.5 / -8.1 | -10.5 / -10.5 | -13.5 / -12.9 | -17.1 / -17.5 |
| permanent water bodies |  |  | -12.1 / 0.8 | -15.1 / -14.7 | -21.1 / -20.7 | -18.8 / -16.5 |  |  |  |
| moss and lichen | 4.1 / 12.2 | 1.5 / 7.3 | -1.8 / 1.7 | -4.5 / -3.9 | -7.0 / -7.8 | -9.6 / -10.2 | -12.0 / -12.4 | -14.8 / -14.8 | -18.0 / -18.5 |

## descending

Ellipsoid incidence over the tiles: 36.9-45.2 deg.

| Class | Blocks | r (raw) | r (incidence removed) | Variance explained by incidence: real / sim | Slope, dB per 10 deg (20-60 deg): real / sim |
|---|---|---|---|---|---|
| tree cover | 319859 | 0.917 | 0.628 | 0.79 / 0.83 | -2.12 / -3.15 |
| grassland | 265592 | 0.861 | 0.562 | 0.66 / 0.85 | -2.27 / -4.94 |
| cropland | 5982 | 0.318 | 0.174 | 0.10 / 0.67 | -1.95 / -3.28 |
| built-up | 8633 | 0.253 | 0.225 | 0.03 / 0.59 | -1.42 / -1.97 |
| bare / sparse vegetation | 54070 | 0.883 | 0.553 | 0.73 / 0.80 | -2.74 / -4.08 |
| permanent water bodies | 19370 | 0.717 | 0.747 | 0.07 / 0.08 | +1.38 / -2.02 |
| moss and lichen | 9989 | 0.894 | 0.524 | 0.75 / 0.86 | -2.86 / -4.21 |

Median dB per 10 deg of local incidence (real / sim):

| Class | 0-10 | 10-20 | 20-30 | 30-40 | 40-50 | 50-60 | 60-70 | 70-80 | 80-90 |
|---|---|---|---|---|---|---|---|---|---|
| tree cover | 3.8 / 8.6 | 1.0 / 5.1 | -2.0 / 1.4 | -4.3 / -1.3 | -6.3 / -4.1 | -8.1 / -7.2 | -10.1 / -10.7 | -12.3 / -12.2 | -15.3 / -13.4 |
| grassland | 3.2 / 13.7 | -0.3 / 9.2 | -4.7 / 2.9 | -7.8 / -3.7 | -9.7 / -9.5 | -11.1 / -11.8 | -12.6 / -13.2 | -14.7 / -14.4 | -17.7 / -16.7 |
| cropland |  |  |  | -8.8 / -7.0 | -9.9 / -9.9 | -11.8 / -11.3 |  |  |  |
| built-up |  |  | -1.6 / -0.2 | -3.4 / -3.0 | -3.9 / -3.9 | -7.0 / -5.8 |  |  |  |
| bare / sparse vegetation | 5.1 / 12.6 | 3.5 / 9.6 | 0.8 / 5.1 | -2.2 / -0.2 | -5.0 / -5.0 | -7.8 / -8.2 | -10.3 / -10.5 | -13.2 / -13.0 | -16.3 / -16.7 |
| permanent water bodies |  |  |  | -21.0 / -19.4 | -20.6 / -20.3 | -17.1 / -19.4 | -16.7 / -18.3 | -20.1 / -21.4 |  |
| moss and lichen | 3.9 / 12.6 | 1.9 / 8.6 | -1.4 / 2.8 | -4.3 / -3.1 | -7.0 / -7.5 | -9.8 / -10.1 | -12.4 / -12.4 | -14.7 / -14.8 | -17.3 / -18.2 |

