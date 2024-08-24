# TIgo
TIgo - A user-friendly R package for Single-cell RNA-Seq data trajectory inference analysis.

### Installation

TIgo can be installed directly from GitHub using the following code:

```
library("devtools")
install_github("Jade1109/TIgo")

# Start the application
library(TIgo)
startTIgo()
```
### Prerequisites

TIgo requires R>= 4.4.0

## About
TIgo provides a comprehesnive pipline of scRNA-seq analysis, espacially for trajctory inference. 
![TIgo workflow](/figures/workflow.png)

### Object Creation and Clustering

Input file : .h5

Process: 
1. quality control (‘nFeature_RNA > 200’, ‘nFeature_RNA < 2500’, and ‘percent.mt < 5’)
2. Log normalization
3. Identification of variable features
4. Data scaling
5. Dimensionality reduction
6. Batch effect correction
7. Neighbour identification
8. Clustering
9. Visulization
10. Subset

Output：UMAP, .rds file

![TIgo Object Creation and Clustering](/figures/OC.png)

### Trajectory Inference

Input file : .rds file

Output: Trajctory Inference Plot

#### Monocle3

1. Data preprocessing
2. Batch effect correction
3. Dimensionality reduction
4. Clustering
5. Trajectory construction
6. Pseudotime Calculation
7. Visualization
   
![TIgo Trajectory Inference (Monocle3)](/figures/TIM.png)

#### Slingshot

1. TI Analysis
2. Construct a Network of Connections Between Clusters
3. Assign Each Cell to the Most Likely Lineage
4. Generate a Data Frame of Cell Progressions
5. Visualization

![TIgo Trajectory Inference (Slingshot)](/figures/TIS.png)

### Differential expression






#### Input 

TIgo accept input file for each stage:
Object creation: .h5 file 
Tracjectory Inference: .rds file
Differenctial Experince: .csv file





