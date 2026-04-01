

# SpatialQPFs
<img width="387" alt="SpatialQPFs_logo" src="https://github.com/Genentech/SpatialQPFs/assets/20170034/d79d1795-5cec-4b31-8a86-b2ab6e595caa">

This repo provides the source code related to the **SpatialQPFs** R library that is reported in the manuscript "Deciphering cell to cell spatial relationship for pathology images using SpatialQPFs" by Xiao Li. 

## 📚 Citation
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

## 🚀 What’s new in v2.0.0
The version 2.0.0 of SpatialQPFs introduces a major new update: 
- *Spatial entropy features* ：These measures quantify cell type distributional heterogeneity, a key aspect in characterizing tissue organization. Spatial entropy has gained wide adoption in multiplexed imaging and spatial omics data analysis, and this update integrates it directly into the SpatialQPFs framework.
- *Graph-based features*： Graph-derived metrics leverage cell–cell graphs to describe spatial organization, providing a more comprehensive toolkit for cell–cell interaction and tissue architecture analysis. 

👉 *For a quick overview of these new capabilities, check the* [`Tutorial-of-SpatialQPFs2.0.0.html`](Tutorial-of-SpatialQPFs2.0.0.html) *for an example workflow.*


## 📚 Citation
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
| `DT_graph_uni_subregion_random()`     | Main function to generate graph features using Delaunay triangulation graph (single cell type). |
| `DT_graph_cross_subregion_random()`     | Main function to generate graph features using Delaunay triangulation graph (multiple cell types). |
| `spat_entropy_alltype()`     | Main function to generate spatial entropy for all cell types. |
| `spat_entropy_onetype()`     | Main function to generate spatial entropy for one cell type. |



## Use SpatialQPFs in Python
`SpatialQPFs_rpy2.ipynb` is a tutorial notebook to illustrate how to use this R package in Python.


## Example workflow
Please check the* [`Tutorial-of-SpatialQPFs2.0.0.html`](Tutorial-of-SpatialQPFs2.0.0.html) *for an example workflow.*
