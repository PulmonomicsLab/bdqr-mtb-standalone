args=commandArgs(trailingOnly=TRUE)
# print(length(args))
if (length(args) != 3 )
  stop("Invalid number of arguments to Rscript.")

script_path <- args[1]
input_file <- args[2]
output_path <- args[3]

# print(script_path)
# print(input_file)
# print(output_path)

set.seed(1)

library(reshape2, quietly=TRUE)
library(caret, quietly=TRUE)

Input <- read.table(input_file, sep = "\t", comment.char = "", header = TRUE, stringsAsFactors=FALSE)
Ind_batch <- Input[,c("POS", "REF", "ALT", "QUAL", "SAMPLE", "AO.1", "DP.1")]
Ind_batch$AO.1[Ind_batch$AO.1=="."] <- 0
Ind_batch$DP.1[Ind_batch$DP.1=="."] <- 0
Ind_batch <- transform(Ind_batch, DP.1 = as.numeric(DP.1))
Ind_batch <- transform(Ind_batch, AO.1 = as.numeric(AO.1))
Ind_batch$MUT <- paste(Ind_batch$REF, Ind_batch$POS, Ind_batch$ALT)
Ind_batch$VAL <- Ind_batch$AO.1/Ind_batch$DP.1
Ind_batch$VAL[Ind_batch$VAL == "NaN"] <- 0.0000
Ind_batch$VAL <- round(Ind_batch$VAL, digits = 4)
Ind_batch_l <- dcast(Ind_batch[,c("SAMPLE", "MUT", "VAL")], SAMPLE~MUT, value.var = "VAL")

A_col <- unlist(read.table(paste0(script_path, "A_col.txt"), sep="\t"))

data_2_bat <- data.frame(matrix(ncol = length(A_col), nrow = nrow(Ind_batch_l)))
colnames(data_2_bat) <- A_col

for(j in 1:nrow(data_2_bat)){
  for(i in colnames(data_2_bat)){
    if(i %in% colnames(Ind_batch_l)){
      data_2_bat[j,i] <- Ind_batch_l[j,i]
    }else{
      data_2_bat[j,i] <- 0.00
    }
  }
}

mlp <- readRDS(paste0(script_path, "mlp_all_features.rds"))

pred_prob_bat <- predict(mlp, data_2_bat, type="prob")
pred_bat <- predict(mlp, data_2_bat)
# print(colnames(pred_prob_bat))
# print(pred_prob_bat)

MLP_out <- pred_prob_bat
for(i in 1:nrow(pred_prob_bat)){
  row_sum <- rowSums(pred_prob_bat[i, c("R", "S")])
  MLP_out[i, "R"] <- MLP_out[i, "R"] / row_sum
  MLP_out[i, "S"] <- MLP_out[i, "S"] / row_sum
}
# print(MLP_out)

# pred_bat <- unlist(pred_prob_svm_R_ind_bat)
output <- data.frame(matrix(ncol = 4, nrow = nrow(pred_prob_bat)))
colnames(output) <- c("Sample", "Resistance", "Susceptible", "Class")
for(i in 1:nrow(pred_prob_bat)){
  n <- which.max(pred_prob_bat[i,])
  output[i,"Sample"] <- Ind_batch_l$SAMPLE[i]
  output[i, "Resistance"] <- MLP_out[i, "R"]
  output[i, "Susceptible"] <- MLP_out[i, "S"]
  output[i, 2:3] <- round(output[i, 2:3], digits = 4)
#   output[i,"CLASS"] <- colnames(pred_prob_bat[i,])[n]
  if(colnames(MLP_out[i,])[n] == "R"){
    output[i,"Class"] <- "Resistance"
  } else{
    output[i,"Class"] <- "Susceptible"
  }
}

cat("\n")
cat("Prediction result with full model\n")
cat("\n")
print(output)
cat("\n")

cat(paste0("Saving the prediction result to \"", output_path, "prediction.tsv\" ... "))
write.table(output, file = paste0(output_path, "prediction.tsv"), row.names = FALSE)
cat("Done\n")


#### SHAP execution with 50 feature model ####

gene_mutation_map <- read.csv(paste0(script_path, "mutation_gene_map.tsv"), sep="\t", row.names = "Mutation")

suppressPackageStartupMessages(library(DALEX, quietly = TRUE))

explainer_f <- readRDS(paste0(script_path, "explainer_50F.rds"))
for(i in 1:length(pred_bat)){
  shap_MLP_f <- predict_parts_shap(explainer_f, new_observation = data_2_bat[i,], B = 4)
#   print(colnames(shap_MLP_f))
  for (j in 1:nrow(shap_MLP_f)) {
    shap_MLP_f[j, "variable_name"] <- paste0(shap_MLP_f[j, "variable_name"], " (", gene_mutation_map[shap_MLP_f$variable_name[j], "Gene"], ")")
    shap_MLP_f[j, "variable"] <- paste0(shap_MLP_f[j, "variable_name"], " = ", shap_MLP_f[j, "variable_value"])
  }

#   pred_label <- paste0("MLP.", pred_bat[i])
#   shap_MLP_class_f <- shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY", ]
  shap_MLP_f <- shap_MLP_f[shap_MLP_f$variable_name != "CATEGORY (NA)", ]
  shap_MLP_class_f_agg <- aggregate(shap_MLP_f$contribution, by=list(shap_MLP_f$variable), FUN=mean)
  colnames(shap_MLP_class_f_agg) <- c("variable", "contribution")
  shap_MLP_class_f_agg <- shap_MLP_class_f_agg[order(abs(shap_MLP_class_f_agg$contribution), decreasing = TRUE),]

  shap_MLP_f_plot <- plot(shap_MLP_f[shap_MLP_f$variable_name != "CATEGORY (NA)",], max_features = 51)

  cat("\n\n")
  cat("SHAP result with 50-feature model - ")
  cat(Ind_batch_l$SAMPLE[i])
  cat("\n\n")
  print(shap_MLP_class_f_agg)
  cat("\n")

  cat(paste0("Saving SHAP result to \"", output_path, "shap_result_50_features_", Ind_batch_l$SAMPLE[i], ".tsv\" ... "))
  write.table(shap_MLP_class_f_agg, file = paste0(output_path, "shap_result_50_features_", Ind_batch_l$SAMPLE[i], ".tsv"), row.names = FALSE)
  cat("Done\n")
  cat(paste0("Saving SHAP plot to \"", output_path, "shap_plot_50_features_", Ind_batch_l$SAMPLE[i], ".svg\" ... "))
  suppressMessages(ggsave(file = paste0(output_path, "shap_plot_50_features_", Ind_batch_l$SAMPLE[i], ".svg"), plot = shap_MLP_f_plot))
  cat("Done\n\n")
}

