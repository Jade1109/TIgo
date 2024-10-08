% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TI.R
\name{runslingshotTrajectory}
\alias{runslingshotTrajectory}
\title{Run Slingshot Trajectory Inference}
\usage{
runslingshotTrajectory(
  seurat_object_file,
  start_cluster,
  end_cluster,
  seurat_cluster,
  color_cells,
  session
)
}
\arguments{
\item{seurat_object_file}{The path to the Seurat object file.}

\item{start_cluster}{The starting cluster for the trajectory.}

\item{end_cluster}{The ending cluster for the trajectory.}

\item{seurat_cluster}{The name of the cluster column in the Seurat object.}

\item{color_cells}{A vector of colors for cells.}

\item{session}{The Shiny session object.}
}
\value{
A list containing the Slingshot output, plot, reduced dimension data, and pseudotime values.
}
\description{
This function runs the Slingshot trajectory inference on a Seurat object. It generates pseudotime values and integrates them into the Seurat metadata.
}
