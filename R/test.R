
# ============================================================
# Test the implementation of collaborative filtering.
# ============================================================

library(caret)
library(recommenderlab)

# Set environment and load algorithm implementation.
# setwd("...")
source("cf_algorithm.R")
source("similarity_measures.R")

# Calculates root-mean-square-error.
rmse_function <- function(predicted_values, actual_values){
  differences <- predicted_values - actual_values
  sqrt ( (1/length(differences)) * sum(differences^2) ) 
}

# Evaluates algorithm using cross validation with given number of folds.
evaluate_cf <- function(ratings_matrix, number_of_folds, ...) {
  
  positions_non_zero_values <- as.data.frame(which(ratings_matrix != 0, arr.ind = T))
  names(positions_non_zero_values) <- c("row", "column")
  size <- nrow(positions_non_zero_values)
  fold_train_indices <- createFolds(y = 1:size, k = number_of_folds, list = TRUE, returnTrain = TRUE)
  
  # Create an empty matrix that will contain all predictions at the end.
  all_predictions <- as(sparseMatrix(i = c(), j = c(), dims = ratings_matrix@Dim, dimnames = ratings_matrix@Dimnames), "dgCMatrix")
  
  for(i in 1:number_of_folds){
    print(paste("iteration", i))
    # Create training and test matrix
    train_index <- fold_train_indices[[i]]
    test_index <- setdiff(1:size, train_index)
    train_subset <- as.matrix(positions_non_zero_values[train_index, ])
    test_subset <- as.matrix(positions_non_zero_values[test_index, ])
    
    train_matrix <- ratings_matrix
    train_matrix[test_subset] <- 0
    train_matrix <- drop0(train_matrix)
    
    predictions_matrix <- predict_cf(train_matrix, test_subset, ...)
    all_predictions[as.matrix(test_subset)] <- predictions_matrix[as.matrix(test_subset)]
  }
  # Predictions may fall out of the required interval (for example [1, 5] if ratings are from that interval).
  all_predictions@x[all_predictions@x > 5] <- 5
  all_predictions@x[all_predictions@x < 1] <- 1
  
  prediction_indices <- as.matrix(which(all_predictions != 0, arr.ind = T))
  actual_values <- ratings_matrix[prediction_indices]
  predicted_values <- all_predictions[prediction_indices]
  rmse <- rmse_function(predicted_values, actual_values)
  
  list(all_predictions, rmse)
}


# =========== Examples (with MovieLense dataset from recommenderlab library) =================================


# =========== Example 1: cross validation, using "ibcf" 
data(MovieLense)
ratings_matrix <- MovieLense@data
start <- Sys.time()
res <- evaluate_cf(ratings_matrix, number_of_folds = 10, alg_method = "ibcf", normalization = TRUE, 
                   similarity_metric = cal_cos, k = 300, make_positive_similarities = FALSE, 
                   rowchunk_size = 1000, columnchunk_size = 2000)
end <- Sys.time()
print(end - start)
print(paste("RMSE: ", res[[2]]))

# ============ Example 2: find predictions for 3 users, using "ubcf"
data(MovieLense)
ratings_matrix <- t(MovieLense@data)
items_to_predict <- 1:nrow(ratings_matrix)
users <- c(5, 10, 30)
prediction_indices <- as.matrix(expand.grid(items_to_predict, users))
# note that predictions may fall out of the required interval (for example [1, 5])
Sys.time()
res <- predict_cf(ratings_matrix, prediction_indices, "ubcf", TRUE, cal_cos, 300, FALSE, 2000, 1000)
Sys.time()


# ==============================================================================================================



