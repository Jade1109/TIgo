% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TI.R
\name{runmonocle3Trajectory}
\alias{runmonocle3Trajectory}
\title{Run Monocle3 Trajectory Inference}
\usage{
runmonocle3Trajectory(
  gene_expression,
  cell_metadata,
  gene_annotation,
  seurat_cluster,
  start_cluster,
  label_groups_by_cluster,
  color_cells_by,
  label_branch_points_input,
  session,
  label_leaves_input,
  group_input
)
}
\arguments{
\item{gene_expression}{A gene expression matrix.}

\item{cell_metadata}{A data frame containing cell metadata.}

\item{gene_annotation}{A data frame containing gene annotations.}

\item{seurat_cluster}{The name of the cluster column in the Seurat object.}

\item{start_cluster}{The starting cluster for the trajectory.}

\item{label_groups_by_cluster}{A logical value indicating whether to label groups by cluster.}

\item{color_cells_by}{The method used to color cells in the plot.}

\item{label_branch_points_input}{A logical value indicating whether to label branch points.}

\item{session}{The Shiny session object.}

\item{label_leaves_input}{A logical value indicating whether to label leaves.}

\item{group_input}{A grouping variable for alignment.}
}
\value{
A list containing the Monocle3 CDS object and a trajectory plot.
}
\description{
This function runs the Monocle3 trajectory inference on processed input data. It creates a cell data set (CDS) and performs preprocessing, dimension reduction, clustering, and graph learning.
}
