# Genome wide association for growth rate at fluctuating and extreme temperatures in _Neurospora crassa_
Genome-wide association for fluctuating and extreme temperatures

This repository describes phenotypic and data and analysis scripts related to the manuscript: [Genome-wide association mapping for growth rate at fluctuating and extreme temperatures](https://www.biorxiv.org/content/10.1101/2024.04.29.591604v1)

## Genotypic data

Note that there is another repository for the nested association mapping population itself: [Neurospora_NAM_population](https://github.com/ikron/Neurospora_NAM_population) and the genotype file is available at [zenodo](https://zenodo.org/records/11120317)

## Phenotypic data

Phenotypic data analysed in this manuscript is in the folder /data/

1. growth_rate_data.csv
   - A csv file that contains the phenotypic measurements for each growth rate assay. The columns are: Genotype = name of the strain, Temp = temperature, Frequency = frequency of the temperature fluctuation, tindex = index for time points, t1 to t8 = values of many mm the mycelium has grown (the first time point should have a value of zero), growthrate = growth rate for the replicate calculting by fitting a linear regression through the points, R2 = R-square value for the fit, Genot = name of the strain that corresponds to name in the genotype file. This data is used to calculate growth rates can generate the 'genomeans.csv' file.
2. mean_growth_rate_data.csv
   - A csv file that contains mean growth rates for all genotypes. Used as input data in GWAS. The columns are: Genot = name of the strain that corresponds to name in the genotype file, Temp = temperature, Frequency = frequency of the temperature fluctuation, meangr = mean growth rate of the three replicates for each strain

## Scripts

Scripts are in folder /scripts/

1. manuscript_codes.R
 - R codes for performing the analyses presented in the manuscript
