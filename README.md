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
![TIgo Object Creation and Clustering](/figures/OC.png)
Input file : .h5
Process: 
1. quality control (‘nFeature_RNA > 200’, ‘nFeature_RNA < 2500’, and ‘percent.mt < 5’)
2. Log normalization![image](https://github.com/user-attachments/assets/e69a6b67-b9cb-4ff2-ab2b-7160ab2c93de)





#### Input 

TIgo accept input file for each stage:
Object creation: .h5 file 
Tracjectory Inference: .rds file
Differenctial Experince: .csv file





