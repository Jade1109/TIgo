![image](https://github.com/user-attachments/assets/31926280-42f3-4857-ba52-0eb07dc047a9)![image](https://github.com/user-attachments/assets/a795f9d0-10db-436b-b5a2-f87a1dcfc0e4)![image](https://github.com/user-attachments/assets/f04a7b1d-ed29-4050-88fa-9714737e48ab)![image](https://github.com/user-attachments/assets/7574a9a0-9e86-484a-bd01-d6c1a89250cb)![image](https://github.com/user-attachments/assets/5fafb0e5-52d3-4376-a85a-a0f71b6b992a)# TIgo
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
TIgo provides a comprehensive pipeline from scratch, compares and combines Monocle3 and Slingshot, and supplies the Trajectory information on the key interferon-stimulated gene and immune cell markers for SARS2 infection (Saelens et al., 2019, Street et al., 2018)

![TIgo workflow](/figures/workflow.png)

### Object Creation and Clustering

Input file : .h5

Output：UMAP, .rds file

#### Workflow: 
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

### Differential Expression Analysis

Input file: Output data from tranjectort inference

Output: file containing differential expression genes, heatmap, scatterplot (Monocle3 only)

#### Monocle3
1. DE analysis
2. Filtration of significantly differential expressed genes
3. Creation of module
4. Aggregate gene expression data
5. Visulization

![TIgo Differential Expression Analysis (Monocle3) a. 10 ISGs expression along with pseudotime. b. heatmap of differential expression.
](/figures/DEM.png)

#### Slingshot

1. Fit GAM model for each gene and extract smooth terms
2. Identification of Significant Genes 
3. Visualization

![TIgo Differential Expression Analysis (Slingshot)](/figures/DES.png)

### GO Enrichment

Input files: .csv file
Output: barplot and file of GO enrichment

![TIgo GO Enrichment](/figures/GO.png)

## Test Data
GSE145926, which includes human bronchoalveolar lavage fluid (BALF) cells, was utilized to test the application of this package in COVID-19 studies. This dataset includes samples from 3 healthy controls, 3 patients with moderate COVID-19, and 6 patients with severe COVID-19 (Liao et al., 2020b).


## Built With

* [shiny](https://shiny.rstudio.com/)
* [RStudio](https://www.rstudio.com/)


## Reference
SAELENS, W., CANNOODT, R., TODOROV, H. & SAEYS, Y. 2019. A comparison of single-cell trajectory inference methods. Nature biotechnology, 37, 547-554.
STREET, K., RISSO, D., FLETCHER, R. B., DAS, D., NGAI, J., YOSEF, N., PURDOM, E. & DUDOIT, S. 2018. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC genomics, 19, 1-16.
LIAO, M., LIU, Y., YUAN, J., WEN, Y., XU, G., ZHAO, J., CHENG, L., LI, J., WANG, X., WANG, F., LIU, L., AMIT, I., ZHANG, S. & ZHANG, Z. 2020b. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. Nat Med, 26, 842-844.












