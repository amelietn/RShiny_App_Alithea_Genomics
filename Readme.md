# Assessment for Bioinformatics Intern

Two count matrices generated using BRB-seq and DRUG-seq technologies are given. The goal is to perform a comprehensive analysis, and to identify which of the matrices, tech1 or tech2, corresponds to BRB-seq or DRUG-seq. Additionally, we would like to determine which condition, referred to as ‘ConditionA’ or ‘ConditionB’, corresponds to the Huh7 cells treated with a TGFB1 activator.

Tech1 ([Count Matrix1 ](Application/count_mat_tech1.txt)) contains 4 replicates per condition (A, B) and Tech2 ([Count Matrix2 ](Application/count_mat_tech2.txt)) contains 24 replicates per condition. For both of the matrices, the row names are ENSEMBL gene IDs, and the columns are the sample names (replicate number and condition). This R shiny app allows to display the results obtained from the analysis of these two count matrices.

## Installation

This app was coded in R version 4.3.3.

Before using this app, ensure that you have installed the required libraries. You can install the necessary libraries using the `install.packages()` function from the R console.

```R
install.packages(c("shiny","shinydashboard", "reshape2", "ggplot2","pheatmap","RColorBrewer","DESeq2","biomaRt","DT"))
````

These libraries were used in version :

shiny 1.8.1.1 

shinydashboard 0.7.2

reshape2 1.4.4

ggplot2 3.5.0

pheatmap 1.0.12

RColorBrewer 1.1.3

DESeq2 1.42.1

biomaRt 2.58.2

DT 0.33


To download both count matrices tech1 and tech 2 : 

[Count Matrix1 ](Application/count_mat_tech1.txt)

[Count Matrix2 ](Application/count_mat_tech2.txt)

## Usage

The easiest way to open the app is to "Run App" (for app.R) in RStudio. Make sure to download the folder [Application](Application) containing [Count Matrix1 ](Application/count_mat_tech1.txt), [Count Matrix2 ](Application/count_mat_tech2.txt) and [app.R](Application/app.R) (you can directly download the Application.zip). If the files are not in the same folder, you might need to change the path to the files when loading the data :

```R
data1 <- as.matrix(read.csv("/path/to/your/directory/with/count_mat_tech1.txt", header = TRUE, row.names = 1))
data2 <- as.matrix(read.csv("/path/to/your/directory/with/count_mat_tech2.txt", header = TRUE, row.names = 1))
```

There is a total of four tabs on the left, each displaying different results. 

The first tab can be used to explore the two count matrices tech1 and tech2, by selecting a gene (row) it displays the reads for all the samples (replicates and conditions) for this particular gene. For example, it is possible to write "ConditionA" in the Search box, and it will display the reads for each replicate treated with condition A for the selected gene.

The second tab, named "Technology Analysis", contains information about the sequencing depth and number of detected genes for both tech1 and tech2.

The third one is "Correlation Heatmaps" and contains the second part of the analysis, as well as the hypothesis of which matrix corresponds to which technology (BRB-seq and DRUG-seq). 

Finally, the last tab displays all the plots and results of the Differential Expression analysis, the PCA and the hypothesis of which condition (A or B) corresponds to the TGFB-treated Huh7 cells. **The plots usually take some time (~ one minute) to appear (the time it takes R to compute the differential expression analysis).**


## Supplementary material : Expected Plots

Here are the plots you are expected to see on the app :

[Correlation_Heatmaps.pdf](https://github.com/amelietn/RShiny_App_Alithea_Genomics/files/14945230/Correlation_Heatmaps.pdf)

[Heatmaps_top50.pdf](https://github.com/amelietn/RShiny_App_Alithea_Genomics/files/14945234/Heatmaps_top50.pdf)

[PCA_Tech1.pdf](https://github.com/amelietn/RShiny_App_Alithea_Genomics/files/14945261/PCA_Tech1.pdf)

[PCA_Tech2.pdf](https://github.com/amelietn/RShiny_App_Alithea_Genomics/files/14945264/PCA_Tech2.pdf)

[PCA_Tech1_Tech2.pdf](https://github.com/amelietn/RShiny_App_Alithea_Genomics/files/14945268/PCA_Tech1_Tech2.pdf)



To see the initial R Notebook containing all the details, calculations and plots, you can download the Supplementary Material.zip. However, all the plots are already displayed on the R Shiny Application.

