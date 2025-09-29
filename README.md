

# SpatialQPFs
<img width="387" alt="SpatialQPFs_logo" src="https://github.com/Genentech/SpatialQPFs/assets/20170034/d79d1795-5cec-4b31-8a86-b2ab6e595caa">

This repo provides the source code related to the *SpatialQPFs* R library that is reported in the manuscript "Deciphering cell to cell spatial relationship for pathology images using SpatialQPFs" by Xiao Li. 

## Citation
If you find this work useful in your research or if you use parts of this code please consider citing this paper:

**Li, X. Deciphering cell to cell spatial relationship for pathology images using SpatialQPFs. Sci Rep 14, 29585 (2024).** https://doi.org/10.1038/s41598-024-81383-1

```
@article{li2024deciphering,
  title={Deciphering cell to cell spatial relationship for pathology images using SpatialQPFs},
  author={Li, Xiao},
  journal={Scientific Reports},
  volume={14},
  number={1},
  pages={1--13},
  year={2024},
  publisher={Nature Publishing Group}
}
```

## What’s new in v2.0.0 (Upcoming)
The upcoming version 2.0.0 of SpatialQPFs will introduce a major new feature: 
- **Spatial entropy measures**
These measures quantify cell type distributional heterogeneity, a key aspect in characterizing tissue organization. Spatial entropy has gained wide adoption in multiplexed imaging and spatial omics data analysis, and this update integrates it directly into the SpatialQPFs framework.
- **Graph-based features**
Graph-derived metrics leverage cell–cell graphs to describe spatial organization, providing a more comprehensive toolkit for cell–cell interaction and tissue architecture analysis. 


The relevant methods has been reviewed and discussed here:

**Li X, Ren X, Venugopal R. Entropy measures for quantifying complexity in digital pathology and spatial omics. iScience. 2025 Jun 20;28(6).** https://www.cell.com/iscience/fulltext/S2589-0042(25)01026-0

```
@article{li2025entropy,
  title={Entropy measures for quantifying complexity in digital pathology and spatial omics},
  author={Li, Xiao and Ren, Xuehan and Venugopal, Raghavan},
  journal={iScience},
  volume={28},
  number={6},
  year={2025},
  publisher={Elsevier}
}
```





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

The utilization of the *SpatialQPFs* functionality begins by importing cell-level input data (currently only support $.csv$ file, more input format will be supported in the future release). Preliminary steps, including identifying tumor regions of interest (ROIs), segmenting and classifying cells, are performed before invoking *SpatialQPFs*. 

<img src="https://github.com/user-attachments/assets/e4ce6bf7-a497-4c83-a896-402c898aeb3c" alt="geostat_data" width="800" height="450"/>


To visualize the spatial distribution of the input tabular data within the original tissue space and its derived spatial density map, users can call ```Data_Vis()``` function.
For example, by specifying ```cell_type = "Lymphocyte"``` users can plot the lymphocyte population. 
The argument ``path`` and ```file``` accept strings that indicate the directory address where the file is saved and the name of the CSV file. Additionally, ```cell_class_var```, ```x_var``` and ```y_var``` specify the column names from the input CSV file that contain the X, Y-coordinates of cell centroids and the corresponding cell type. See function manual for further details. 

```R
Data_Vis(path = folder_path, file = file_name, cell_class_var = "cell_class", 
         x_var = "coordinate_X", y_var = "coordinate_Y", cell_type = "Lymphocyte")  
```

<img src="https://github.com/user-attachments/assets/d1395016-2ab3-48a2-9ca8-12e0ca12fc3a" alt="data_vis" width="800" height="450"/>



### 3. Point pattern analysis

Point process data analysis methods can be conducted by implementing three functions:
- ```Point_pattern_data_uni()```: This function generates spatial features for a single cell type population.
- ```Point_pattern_data_bi()```: It computes the spatial interaction for a pair of cell types.
- ```Point_pattern_data_ITLR()```: This function discriminates a target cell population into subgroups based on their spatial relationships with the reference cell population
  
Additional plots for other point process data analysis methods are presented at ```Tutorial-of-SpatialQPFs.html```.

```R
Point_pattern_data_uni(path = folder_path, file = file_name, 
                       cell_class_var = "cell_class", x_var = "coordinate_X", 
                       y_var = "coordinate_Y", cell_type = "Lymphocyte", 
                       scale = 200, myplot = T)

Point_pattern_data_bi(path = folder_path, file = file_name,  
                      cell_class_var = "cell_class", x_var = "coordinate_X", 
                      y_var = "coordinate_Y", from_type = "Lymphocyte", 
                      to_type = "Tumor", scale = 200, myplot = T)

Point_pattern_data_ITLR(path = folder_path, file = file_name, 
                        cell_class_var = "cell_class", x_var = "coordinate_X", 
                        y_var = "coordinate_Y", from_type = "Lymphocyte", 
                        to_type = "Tumor", micron_per_pixel = 0.5, myplot = T)
```                       

<img src="https://github.com/user-attachments/assets/b68bb32d-918a-490a-b9b2-d14864b5160f" alt="point_data" width="700" height="900"/>



### 4. Areal data analysis
The analysis of areal data can be performed using the ```Areal_data()``` function. 
Note that, this function conducts the transformation from point process data to areal data behind the scenes, by partitioning the underlying space into square lattices. Users can specify the desired side length of the square lattices as ```2*scale```. 

Other areal data analysis plots can be found in ```Tutorial-of-SpatialQPFs.html```.

```R
Areal_data(path = folder_path, file = file_name,
           cell_class_var = "cell_class", x_var = "coordinate_X", 
           y_var = "coordinate_Y", from_type = "Lymphocyte", 
           to_type = "Tumor", scale = 200, myplot = T)     
```

<img src="https://github.com/user-attachments/assets/5dee5e4f-d2b4-43ae-ac1f-9bad86dc6d8c" alt="areal_data" width="700" height="800"/>



### 5. Geostatistical data analysis
The analysis of geostatistical data can be conducted using the ```Geostatistics_data()``` function. In the current version, both semi-variogram and cross-variogram are implemented to illustrate the spatial correlation of either a single cell type or a pair of different cell types, respectively.

Further plots for alternative geostatistical data analysis methods are demonstrated at ```Tutorial-of-SpatialQPFs.html```.

```R
Geostatistics_data(path = folder_path, file = file_name, 
                   cell_class_var = "cell_class", x_var = "coordinate_X", 
                   y_var = "coordinate_Y",  from_type = "Lymphocyte", 
                   to_type = "Tumor", scale = 200, myplot = T)  
```

<img src="https://github.com/user-attachments/assets/26276fb8-a938-4e4c-b5b4-948dd55a003d" alt="geostat_data" width="650" height="550"/>



