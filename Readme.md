![image](https://img.shields.io/badge/Code-R-blue)
![image](https://img.shields.io/badge/Package-R-blue)
![image](https://img.shields.io/badge/ABI-V%200.4-blue)
## Install Package

```{r}
install.packages("devtools")
devtools::install_github("slphyx/ABI")
```

or

```{r}
install.packages("devtools")
library(devtools)
install_github("slphyx/ABI")
```


## Use Package

```{r}
ABI::run_app_ABI()
```

or

```{r}
library(ABI)
run_app_ABI()
```

## Use function
```{r}
ABI_Helminth()
ABI_Helminth(0.06)
ABI_Helminth(0.02,"NS","18S rRNA")
ABI_Helminth(distance = 0.5,group = "CE",marker = "ITS2")
```

### Use function With Fasta file
```{r}
library(ape)
sequences <- read.dna("file.fasta",format = "fasta")
# Convert sequences to distance matrix
dist_matrix <- dist.dna(sequences, model = "raw")

# Construct neighbor-joining tree
nj_tree <- nj(dist_matrix)
genetic_distance <- cophenetic(nj_tree)

# set column and row name
colM <- colnames(genetic_distance)
rowM <- rownames(genetic_distance)

# select taxa
distance_selected  <- genetic_distance[rowM[1],colM[2]]
# Use ABI_Helminth()
ABI_Helminth(distance = distance_selected,
group = "CE",marker = "ITS2")
```
