# SpatialQPFs
<img width="387" alt="SpatialQPFs_logo" src="https://github.com/Genentech/SpatialQPFs/assets/20170034/d79d1795-5cec-4b31-8a86-b2ab6e595caa">

This repo provides the source code related to the *SpatialQPFs* R library that is reported in the manuscript "SpatialQPFs: An R package for deciphering cell-cell spatial relationship" by Xiao Li. 

## Primary SpatailQPFs functions

| Function                   | Description |
|----------------------------|-------------|
| `Data_Vis()`               | Visualizes the input cell level data, including raw spatial map of individual cells and smoothed spatial density map |
| `Point_pattern_data_uni()` | Main function to generate spatial features using point process data (single cell type). |
| `Point_pattern_data_bi()`  | Main function to generate spatial features using point process data (pairwise cell types). |
| `Point_pattern_data_ITLR()`| This is a function that discriminates a target cell population into 2 subgroups, in spatial relationships with the reference cell population. In case of lymphocytes and tumor cells, the lymphocytes are discriminated into intra-tumoral lymphocytes and adjacent-tumoral lymphocytes. |
| `Areal_data()`             | Main function to generate spatial features using areal data. |
| `Geostatistics_data()`     | Main function to generate spatial features using geostatistical data. |



## Example workflow

### 1. installation

For first-time use, installation can be completed in the user's local computing environment.

```{R}
devtools::install("your_path/SpatialQPFs", dependencies = TRUE, repos="https://cloud.r-project.org/")
```


### 2. Input data visualization

Utilization of *SpatialQPFs*  functionality begins with importing input data, followed by preliminary steps such as identifying tumor ROIs, segmenting and classifying cells, before invoking *SpatialQPFs*.

<img src="https://github.com/user-attachments/assets/e4ce6bf7-a497-4c83-a896-402c898aeb3c" alt="geostat_data" width="800" height="450"/>


Visualization of the input tabular data in original tissue space and its derived spatial density map, users can call ```Data_Vis()``` function.

<img src="https://github.com/user-attachments/assets/1bff6ee3-3f80-4547-8149-5a321e93f587" alt="geostat_data" width="900" height="550"/>



### 3. Point pattern analysis

Point pattern data analysis methods can be conducted by implementing three functions:
- ```Point_pattern_data_uni()```: This function generates spatial features for a single cell type population.
- ```Point_pattern_data_bi()```: It computes the spatial interaction for a pair of cell types.
- ```Point_pattern_data_ITLR()```: This function discriminates a target cell population into subgroups based on their spatial relationships with the reference cell population
<img src="https://github.com/user-attachments/assets/056f9515-7ed4-4e5c-8e87-e9b2f9275227" alt="geostat_data" width="800" height="1000"/>



### 4. Areal data analysis
The analysis of areal data can be performed using the ```Areal_data()``` function. Note that, this function conducts the transformation from point pattern data to areal data behind the scenes, by partitioning the underlying space into square lattices.
<img src="https://github.com/user-attachments/assets/8ab5b844-d6e5-4237-bb76-d7e1097716b2" alt="geostat_data" width="800" height="500"/>



### 5. Geostatistical data analysis
The analysis of geostatistical data can be conducted using the ```Geostatistics_data()``` function. In the current version of *SpatialQPFs* , both semi-variogram and cross-variogram are implemented to illustrate the spatial correlation of either a single cell type or a pair of different cell types, respectively.

<img src="https://github.com/user-attachments/assets/c760f948-5ad8-4aed-8f54-ce7bb3044d5f" alt="geostat_data" width="800" height="600"/>



