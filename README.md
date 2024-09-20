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

```R
## the folder contain the source code names "SpatialQPFs" in the user's local directory
devtools::install("your_path/SpatialQPFs", dependencies = TRUE, repos="https://cloud.r-project.org/")
```
Then, the package can be loaded:
```R
library("SpatialQPFs")
```

### 2. Input data visualization

The utilization of the $SpatialQPFs$ functionality begins by importing cell-level input data (currently only support $.csv$ file, more input format will be supported in the future release). Preliminary steps, including identifying tumor regions of interest (ROIs), segmenting and classifying cells, are performed before invoking *SpatialQPFs*. 

<img src="https://github.com/user-attachments/assets/e4ce6bf7-a497-4c83-a896-402c898aeb3c" alt="geostat_data" width="800" height="450"/>


To visualize the spatial distribution of the input tabular data within the original tissue space and its derived spatial density map, users can call ```Data_Vis()``` function.
For example, by specifying ```cell_type = "Lymphocyte"``` users can plot the lymphocyte population. 
The argument ``path`` and ```file``` accept strings that indicate the directory address where the file is saved and the name of the CSV file. Additionally, ```cell_class_var```, ```x_var``` and ```y_var``` specify the column names from the input CSV file that contain the X, Y-coordinates of cell centroids and the corresponding cell type. See function manual for further details. 

```R
Data_Vis(path = folder_path, file = file_name, cell_class_var = "cell_class",
         x_var = "X0", y_var = "X1", cell_type = "Lymphocyte")  
```

<img src="https://github.com/user-attachments/assets/efcf11b0-7206-490a-b5de-508ea8610c87" alt="geostat_data" width="800" height="450"/>



### 3. Point pattern analysis

Point process data analysis methods can be conducted by implementing three functions:
- ```Point_pattern_data_uni()```: This function generates spatial features for a single cell type population.
- ```Point_pattern_data_bi()```: It computes the spatial interaction for a pair of cell types.
- ```Point_pattern_data_ITLR()```: This function discriminates a target cell population into subgroups based on their spatial relationships with the reference cell population
Additional plots for other point process data analysis methods are presented at ```Tutorial-of-SpatialQPFs.html```.

```R
Point_pattern_data_uni(path = folder_path, file = file_name, 
                       cell_class_var = "cell_class", x_var = "X0", y_var = "X1", 
                       cell_type = "Lymphocyte", scale = 200, myplot = T)

Point_pattern_data_bi(path = folder_path, file = file_name,  
                      cell_class_var = "cell_class", x_var = "X0", y_var = "X1", 
                      from_type = "Lymphocyte", 
                      to_type = "Tumor", 
                      scale = 200, myplot = T)

Point_pattern_data_ITLR(path = folder_path, file = file_name, 
                        cell_class_var = "cell_class", x_var = "X0", y_var = "X1", 
                        from_type = "Lymphocyte", 
                        to_type = "Tumor",  
                        micron_per_pixel = 0.5, 
                        myplot = T)
```                       

<img src="https://github.com/user-attachments/assets/056f9515-7ed4-4e5c-8e87-e9b2f9275227" alt="geostat_data" width="700" height="900"/>



### 4. Areal data analysis
The analysis of areal data can be performed using the ```Areal_data()``` function. 
Note that, this function conducts the transformation from point process data to areal data behind the scenes, by partitioning the underlying space into square lattices. Users can specify the desired side length of the square lattices as ```2*scale```. 


```R
Areal_data(path = folder_path, file = file_name,
           cell_class_var = "cell_class", x_var = "X0", y_var = "X1", 
           from_type = "Lymphocyte", 
           to_type = "Tumor",
           scale = 200, 
           myplot = T)  
```

<img src="https://github.com/user-attachments/assets/954194fa-f9ad-4a91-912e-f02988811cee" alt="geostat_data" width="700" height="900"/>



### 5. Geostatistical data analysis
The analysis of geostatistical data can be conducted using the ```Geostatistics_data()``` function. In the current version of *SpatialQPFs* , both semi-variogram and cross-variogram are implemented to illustrate the spatial correlation of either a single cell type or a pair of different cell types, respectively.

```R
Geostatistics_data(path = folder_path, file = file_name, 
                   cell_class_var = "cell_class", x_var = "X0", y_var = "X1", 
                   from_type = "Lymphocyte", 
                   to_type = "Tumor", 
                   scale = 200, 
                   myplot = T) 
```

<img src="https://github.com/user-attachments/assets/c760f948-5ad8-4aed-8f54-ce7bb3044d5f" alt="geostat_data" width="700" height="500"/>



