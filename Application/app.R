library(shiny)
library(shinydashboard)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(DT)

# Loading data from count matrices tech1 and tech2
data1 <- as.matrix(read.csv("count_mat_tech1.txt", header = TRUE, row.names = 1))
data2 <- as.matrix(read.csv("count_mat_tech2.txt", header = TRUE, row.names = 1))

# Function to calculate Sequencing Depth
# Sum the reads of all the genes (rows) for each sample (column) and compute the average of these sums
seq_depth <- function(count_matrix){
  total_reads = c()
  for(sample in 1:ncol(count_matrix)){
    reads_per_sample = 0
    for(gene in 1:nrow(count_matrix)){
      reads_per_sample = reads_per_sample + count_matrix[gene,sample]
    }
    total_reads = append(total_reads, reads_per_sample)
  }
  
  mean_reads = mean(total_reads)
  return(mean_reads)
}

# Function to calculate NODG
# Number of genes (rows) with non-zero read counts in each matrix 
NODG <- function(count_matrix){
  NODG = 0
  for(gene in 1:nrow(count_matrix)){
    if (sum(count_matrix[gene, ] > 0) > 0) {
      NODG <- NODG + 1
    }
  }
  return(NODG)
}

# Function for correlation heatmaps
corr_heatmap <- function(count_matrix,tech_nbr){
  log_norm <- log1p(count_matrix)
  correlation <- cor(log_norm) # Correlation matrix
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  par(mfrow=c(1,2))
  pheatmap(correlation, main=paste(tech_nbr," Replicates Correlation"), cluster_rows = FALSE, cluster_cols = FALSE, color=color_palette)
}

# Function for Differential Expression Analysis between condition A and B
DE_analysis <- function(count_matrix,tech_nbr){
  # Creating a list of conditions and renaming them 'A' and 'B'
  conditions <- colnames(count_matrix)
  for (i in 1:length(conditions)){
    if(grepl("ConditionA", conditions[i])){
      conditions[i] = "A"
    }else if(grepl("ConditionB", conditions[i])){
      conditions[i] = "B"
    }
  }
  sample_dataframe <- data.frame(conditions)
  # Creating DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_dataframe, design = ~ conditions)
  # Differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds)
  # Getting top 50 DE genes
  top50 <- head(res[order(res$pvalue), ], 50)
  # Translating Ensembl IDs to gene symbols 
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
  top50_symbols <- ensembl_to_symbol[ensembl_to_symbol$ensembl_gene_id %in% rownames(top50), "external_gene_name"]
  # Gene symbols for TGFB1
  tgfb1_symbol <- ensembl_to_symbol[ensembl_to_symbol$external_gene_name == "TGFB1", "ensembl_gene_id"]
  
  # Adding TGFB1 to the top 50 genes list
  top50_symbols_with_tgfb1 <- c(top50_symbols, "TGFB1")
  top_counts <- count_matrix[c(rownames(top50), tgfb1_symbol), ]
  # Heatmaps
  heatmap <- log2(top_counts + 1) 
  rownames(heatmap) <- c(top50_symbols, "TGFB1")
  pheatmap(heatmap, cluster_rows = TRUE, cluster_cols = FALSE, main = paste("Top 50 DE Genes for ", tech_nbr),fontsize_row = 6)
}

# Function for comparing DE analysis results obtained from tech1 and tech2
comparing_DE <- function(top50_tech1_symbols, top50_tech2_symbols){
  common_DEgenes <- intersect(top50_tech1_symbols, top50_tech2_symbols)
  unique_DEgenes1 <- setdiff(top50_tech1_symbols, top50_tech2_symbols)
  unique_DEgenes2 <- setdiff(top50_tech1_symbols, top50_tech2_symbols)
  return(common_DEgenes)
}

# Function for PCA and scatter plot
pca <- function(tech1, tech2, tech_nbr){
  if(tech_nbr=="Tech 1 and Tech 2"){
    #PCA for tech1 and tech2 combined
    library(stats)
    combined_data <- cbind(data1, data2)
    pca <- prcomp(t(combined_data))
    pc_scores <- as.data.frame(pca$x)
    pc_scores$Technology <- c(rep("tech1", ncol(data1)), rep("tech2", ncol(data2)))
    pc_scores$Condition <- c(rep("conditionA", 4), rep("conditionB", 4), rep(rep(c("conditionA", "conditionB"), each = 3), times = 8))
    ggplot(pc_scores, aes(x = PC1, y = PC2, color = Condition, shape = Technology)) +
      geom_point(size = 3) +
      labs(title = "PCA Plot of Replicates by Condition for Tech1 and Tech2", x = "PC1", y = "PC2") +
      theme_minimal()
  }else if(tech_nbr=="Tech 1"){
    #PCA for only tech1
    pca_tech1 <- prcomp(t(tech1))
    pc_scores1 <- as.data.frame(pca_tech1$x)
    pc_scores1$Condition <- rep(c("conditionA", "conditionB"), each = 4)
    return(
    ggplot(pc_scores1, aes(x = PC1, y = PC2, color = Condition)) + geom_point() + labs(title = "PCA Plot of Replicates by Condition for Tech1", x = "PC1", y = "PC2") + theme_minimal()
    )
  }else if(tech_nbr=="Tech 2"){
    #PCA for only tech2
    pca_tech2 <- prcomp(t(tech2))
    pc_scores2 <- as.data.frame(pca_tech2$x)
    pc_scores2$Condition <- rep(rep(c("conditionA", "conditionB"), each = 3), times = 8)
    return(
      ggplot(pc_scores2, aes(x = PC1, y = PC2, color = Condition)) + geom_point() + labs(title = "PCA Plot of Replicates by Condition for Tech2", x = "PC1", y = "PC2") + theme_minimal()
    )
  }else{
    return("Invalid tech_nbr value. Please provide 'Tech 1 and Tech 2', 'Tech 1', or 'Tech 2'.")
    }
}

# Function for PC loadings
pc <- function(tech1, tech2){
  combined_data <- cbind(tech1, tech2)
  pca <- prcomp(t(combined_data))
  loadings_pc1 <- abs(pca$rotation[, 1])
  rotation_row_names <- rownames(pca$rotation)
  top_genes_pc1 <- rotation_row_names[order(loadings_pc1, decreasing = TRUE)] # Identifying top genes contributing to PC1
  N <- 10  # Number of top genes to select
  top_N_genes_pc1 <- top_genes_pc1[1:N]
  
  #print(top_N_genes_pc1_symbols)
  
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
  top_N_genes_pc1_symbols <- ensembl_to_symbol[ensembl_to_symbol$ensembl_gene_id %in% top_N_genes_pc1, "external_gene_name"]
  
  loadings_pc2 <- abs(pca$rotation[, 2])
  rotation_row_names2 <- rownames(pca$rotation)
  top_genes_pc2 <- rotation_row_names2[order(loadings_pc2, decreasing = TRUE)] # Identifying top genes contributing to PC2
  N <- 10  # Number of top genes to select
  top_N_genes_pc2 <- top_genes_pc2[1:N]
  
  top_N_genes_pc2_symbols <- ensembl_to_symbol[ensembl_to_symbol$ensembl_gene_id %in% top_N_genes_pc2, "external_gene_name"]

  print(top_N_genes_pc2_symbols)
}

# R SHINY APP CODE

ui <- dashboardPage(
  dashboardHeader(title = "Alithea Genomics"),
  dashboardSidebar(
    div(
      style = "margin: 10px;",
      h3("Task"),
      HTML("<p class='small-text'>We are given two count matrices, tech1 and tech2. The objective is to identify which of the matrices, tech1 or tech2, corresponds to BRB-seq or DRUG-seq. Additionally, we would like to determine which condition, referred to as ‘ConditionA’ or ‘ConditionB’, corresponds to the Huh7 cells treated with a TGFB1 activator.</p>")
    ),
    
    tags$hr(),
    
    sidebarMenu(
      menuItem("Count Matrices", tabName = "count_matrices", icon = icon("info")),
      menuItem("Technology Analysis", tabName = "tech_analysis_info", icon = icon("info")),
      menuItem("Replicate Correlation Heatmaps", tabName = "correlation_heatmaps", icon = icon("chart-line")),
      menuItem("Differential Expression Analysis", tabName = "DE_analysis", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tags$script("
      $(document).ready(function() {
        $('#link-to-tab').click(function() {
          $('a[data-value=correlation_heatmaps]').tab('show');
        });
      });
    "),
    tabItems(
      tabItem("count_matrices",
              fluidRow(
                column(6,
                       selectInput("selected_gene1", "Select Gene (Tech 1) :", choices = rownames(data1))
                ),
                column(12,
                       DTOutput("data1")
                ),
                column(6,
                       selectInput("selected_gene2", "Select Gene (Tech 2) :", choices = rownames(data2))
                ),
                column(12,
                       DTOutput("data2")
                )
              )
      ),
      tabItem("tech_analysis_info",
              fluidRow(
                  div(
                    style = "margin: 10px;",
                    h3("Technology Analysis :"),
                    HTML("<p class='small-text'>We start by comparing their respective sequencing depths, corresponding to the number of times each nucleotide in a genome is read during the sequencing experiment. 
                         These results were obtained by summing the reads of all the genes (rows) for each sample (column) and computing the average of these sums. 
                         This process gives an average of the number of reads across all the samples. We observe that tech1 has a much higher sequencing depth compared to tech2, indicating that the samples were sequenced deeper using the technology 1.</p>"),
                    HTML("<p class='small-text'>We then look at the number of detected genes (NODG), corresponding to the number of unique genes for which expression was detected.
                         For tech1, the number of detected genes is 17550, which corresponds to the total number of genes present in the count matrix. 
                         For tech2 however, this number is smaller, indicating that technology 1 detected more gene expression.</p>")
                  ),
                  box(textOutput("tech_analysis_info1")),
                  box(textOutput("tech_analysis_info2")),
                  div(
                    style = "margin: 10px;",
                  HTML("<p class='small-text'>Before drawing any conclusions, let's look at the correlation between replicates in each matrix. This step will help us evaluate the relationship between the different replicates, ensuring the reliability and consistency of the data.</p>")
              ),
              box("Click ",
                  tags$a(id = "link-to-tab", "here"),
                  "to go to Correlation Heatmaps.")
              )
      ),
      tabItem("correlation_heatmaps",
              fluidRow(
                div(
                  style = "margin: 10px;",
                  HTML("<p class='small-text'>We first note that in both matrices, the replicates within a specific condition (A or B) display a strong correlation with all the other replicates within the same condition but exhibit lesser correlation with replicates from the opposing condition. 
                  In particular, we also observe that overall, the replicates exhibit a higher correlation in tech1 than in tech2. 
                  Hence, technology 1 seems to have provided more precise and accurate results, since replicates from a same condition displayed a stronger correlation.</p>"),
                  HTML("<p class='small-text'>Based on the information provided, we know that BRB-seq typically yields higher RNA content because of its RNA purification step, whereas DRUG-seq can directly work with cell lysates. 
                       Therefore, BRB-seq RNA purification step might lead to cleaner samples, resulting in a more precise and better performance. 
                       This hypothesis aligns with the higher replicates’ correlation observed in tech1, as well as the higher number of detected genes (17550) 
                       as well as the need for deeper sequencing due to the higher RNA content. Hence, technology 1 probably corresponds to BRB-seq. 
                       Similarly, DRUG-seq might lead to potential variability and lower RNA content in cell lysates since it can sequence RNA from damaged cells after drug treatment, which seems to align with the results from the analysis of tech2.</p>")
                ),
                box(title = "Tech 1 Correlation Heatmap", plotOutput("correlation_heatmap1")),
                box(title = "Tech 2 Correlation Heatmap", plotOutput("correlation_heatmap2"))
              )
      ),
      tabItem("DE_analysis",
              fluidRow(
                div(
                  style = "margin: 10px;",
                  HTML("<p class='small-text'>Based on the information provided, we know that each matrix tech1 and tech2 has numerous replicates in two different Conditions A and B. 
                       To guess which of Condition A or B was treated with a TGFb activator, let’s analyze the heatmaps obtained after performing DE Analysis between A and B.</p>"),
                  HTML("<p class='small-text'>Because each cell in the heatmap corresponds to the expression level of a gene in the specific replicate and condition, we are able to visually assess the differences in gene expression levels between Conditions A and B.
                       Here, brighter colors (red-orange) indicate higher expression and darker colors (blue) indicate lower expression. 
                       Including TGFB1 in the heatmap allowed to compare its expression levels in replicates from both conditions, which suggests that Condition B probably corresponds to the TGFB1-treated Huh7 cells. 
                       Precisely, TGFB1 expression was higher (brighter colors) in replicates from Condition B in both tech1 and tech2 heatmaps, whereas replicates from Condition A displayed darker blue colors, suggesting lower expression.</p>")),
                box(title = "Tech 1 DE Heatmap", plotOutput("DE_heatmap1")),
                box(title = "Tech 2 DE Heatmap", plotOutput("DE_heatmap2")),
                div(
                  style = "margin: 10px;",
                  HTML("<p class='small-text'>To further explore and compare differential expression analysis results obtained from both technologies, let’s look at the genes displaying the highest expression levels in tech1 and tech2. 
                  Assuming condition B corresponds to the TGFB1-treated Huh7 cells, it makes sense to focus on the genes under this condition (B), using condition A as a control. 
                  Analyzing genes that display higher expression in both technologies under this condition can provide valuable insights into genes being responsive to TGFB1 treatment. 
                  For tech1, we observe that HSPA5 is very highly expressed in condition B. This Heat shock protein is an endoplasmic reticulum chaperone that regulates cell metabolism, in particular lipid metabolism. 
                  Interestingly, TGFB1 is one of the top HSPA5 binding genes (<a href='https://www.nature.com/articles/s41435-023-00205-y'>reference</a>), supporting the hypothesis that condition B is likely to correspond to the TGFB1-treated cells.</p>"),
                  HTML("<p class='small-text'>I then observed the common differentially expressed genes between tech1 and tech2. One of the obtained genes is VIM, which displays higher expression levels in condition B compared to A, especially in tech2. 
                  It is known that this protein involved in maintaining cell structure and integrity is implicated in numerous cellular processes, and TGFB1 has been shown to induce its expression in the processing of epithelial-mesenchymal transition (EMT) (<a href='https://pubmed.ncbi.nlm.nih.gov/27466403/'>reference</a>). 
                  This is yet another example in support of the hypothesis that cells are TGFB1-treated in condition B.</p>"),
                  HTML("<p class='small-text'>To further analyze this, I then performed a Principal Components Analysis (PCA). If the samples cluster well by conditions, it suggest the presence of systematic differences in gene expression between conditions A and B :</p>")
                  ),
                
                box(title = "Tech 1 PCA Plot", plotOutput("PCA_plot1")),
                box(title = "Tech 2 PCA Plot", plotOutput("PCA_plot2")),
                div(
                  style = "margin: 10px;",
                  HTML("<p class='small-text'>Supplementary PCA plot using data from both Tech1 and Tech2 :</p>")
                ),
                box(title = "Tech 1 and Tech 2 PCA Plot", plotOutput("PCA_plot3")),
                div(
                  style = "margin: 10px;",
                  HTML("<p class='small-text'>We observe that there is a clear separation between sample clusters from condition A and B for all three of the plots, indicating that the samples cluster well by conditions. 
                  In PCA plot of replicates by condition using only tech1, the conditions are separated along PC1, which is also the case for the PCA plot for tech2. 
                  On the PCA of both tech1 and tech2 data combined, the data obtained using BRB-seq and DRUG-seq is clustered and separated along PC1, whereas the conditions A and B clusters are separated along PC2, 
                  suggesting that PC2 captures a significant source of variation that distinguishes between the two conditions. This could indicate genes or pathways that differently respond to the TGFB1 treatment.</p>"),
                  HTML("<p class='small-text'>To identify the top genes contributing to PC2 on this plot, we can also look at the components loadings (coefficients of the variables in the linear combination) : </p>")
                ),
                box(textOutput("pc_loadings"))
              )
      )
    )
  )
)

server <- function(input, output,session) {

  output$data1 <- renderDT({
    updateSelectInput(session, "selected_gene1", choices = rownames(data1))
    req(input$selected_gene1)
    selected_row <- as.data.frame(t(data1[input$selected_gene1, , drop = FALSE]))
    datatable(selected_row, rownames = TRUE)
  })
  output$data2 <- renderDT({
    updateSelectInput(session, "selected_gene2", choices = rownames(data2))
    req(input$selected_gene2)
    selected_row <- as.data.frame(t(data2[input$selected_gene2, , drop = FALSE]))
    datatable(selected_row, rownames = TRUE)
  })
  
  output$tech_analysis_info1 <- renderText({
    paste("Sequencing Depth Tech 1: ",seq_depth(data1),"    ,       Sequencing Depth Tech 2: ",seq_depth(data2))
  })
  output$tech_analysis_info2 <- renderText({
    paste("NODG Tech 1: ",NODG(data1),"    ,  NODG Tech 2: ",NODG(data2))
  })
  
  
  output$correlation_heatmap1 <- renderPlot({
    corr_heatmap(data1, "Tech1")
  })
  output$correlation_heatmap2 <- renderPlot({
    corr_heatmap(data2, "Tech2")
  })
  
  output$DE_heatmap1 <- renderPlot({
    DE_analysis(data1, "Tech1")
  })
  output$DE_heatmap2 <- renderPlot({
    DE_analysis(data2, "Tech2")
  })
  
  output$PCA_plot1 <- renderPlot({
    pca(data1,data2, "Tech 1")
  })
  output$PCA_plot2 <- renderPlot({
    pca(data1,data2, "Tech 2")
  })
  output$PCA_plot3 <- renderPlot({
    pca(data1,data2, "Tech 1 and Tech 2")
  })
  output$pc_loadings <- renderText({
    pc(data1,data2)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
