% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Geostatistics_data.R
\name{Geostatistics_data}
\alias{Geostatistics_data}
\title{Main function to generate spatial features using geostatistics data}
\usage{
Geostatistics_data(
  path = "/Users/lix233/Haystack/5862_cell_centers/",
  file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
  cell_class_var = "cell_class",
  x_var = "X0",
  y_var = "X1",
  from_type = "Lymphocyte",
  to_type = "Tumor",
  scale,
  myplot = FALSE
)
}
\arguments{
\item{path}{The path for the directory that contains the data}

\item{file}{The file name in the directory}

\item{cell_class_var}{The column name in each file that indicates the cell class variable}

\item{x_var}{The column name in each file that indicates the cell location x variable}

\item{y_var}{The column name in each file that indicates the cell location y variable}

\item{from_type}{The cell type one wants to use as the "from" cell type}

\item{to_type}{The cell type one wants to use as the "to" cell type}

\item{scale}{The spatial range that user wants to investigate}

\item{myplot}{Whether to plot the results, if available, by default it is set as FALSE}
}
\value{
This function returns the features for point referenced data, using geostatistics methods
}
\description{
This is a function that calculates the features for point referenced data, using geostatistics methods. 
Semi-variogram and crossvariogram methods are implemented.
}
\author{
Xiao Li, \email{li.xiao@gene.com}
}
