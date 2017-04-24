
# ============================================================
# Functions used in implementation of collaborative filtering.
# ============================================================

library(Matrix)
library(recommenderlab)
library(slam)
library(data.table)
library(foreach)
source("partial_computations.R")


#' Calculates rating predictions according to CF formula.
#'
#' @param ratings_matrix (dgCMatrix)
#' @param similarity_matrix (dgCMatrix)
#' @returns Matrix of predictions (dgCMatrix)
calculate_predictions <- function(ratings_matrix, similarity_matrix){ 
  
  predictions <- ratings_matrix %*% similarity_matrix
  
  # Calculate normalization factor (sum of abs similarities).
  # Put 1 in places where ratings exist.
  ratings_matrix@x <- rep(1, length(ratings_matrix@x))
  similarity_matrix@x <- abs(similarity_matrix@x)
  sum_abs_similarities <- ratings_matrix %*% similarity_matrix
  
  # TODO: Check if sum_abs_similarities@x can contain zeros
  predictions@x <- predictions@x / sum_abs_similarities@x
  predictions
}

#' Calculates similarities between columns \code{columns_to_consider} from \code{matrix} vs all columns from \code{matrix}
#'
#' @param matrix (dgCMatrix)
#' @param columns_to_consider (vector of integers) vector of indices of columns from \code{matrix}.
#' @param similarity_metric Function used to calculate similarities. It has to accept two matrices (dgcMatrix) and  calculate similarities between columns.
#' @param make_positive_similarities (logical) Whether all similarities should be modified by a factor in order to have only positive similarities.
#' @param k (integer) Number of largest similarites to keep, per column (k nearest neighbours approach).
#' @returns similarities matrix (dgCMatrix)
#' @note  required Matrix, slam, data.table
find_similarities <- function(matrix, columns_to_consider, similarity_metric, make_positive_similarities, k){
  
  selected_columns <- matrix[, columns_to_consider, drop=FALSE]
  # similarities should be dgCMatrix with explicit zeros in places where similarity is zero.
  similarities <- similarity_metric(matrix, selected_columns)
  
  # In order to keep explicit zeros we will change them to some value close to zero.
  # Then we will set similarities of users/item to themselves to zero and drop those values.
  similarities@x [similarities@x == 0] <- 0.000001
  ind <- cbind(columns_to_consider, 1:length(columns_to_consider))
  similarities[ind] <- 0
  similarities <- drop0(similarities)
  
  # Make all similarities positive, if requested.
  if(make_positive_similarities) {
    if(min(similarities@x) < 0) similarities@x <- similarities@x + abs(min(similarities@x))
  }
  
  if(!is.null(k) && k < nrow(similarities) - 1){ # if ALL - 1 that means we need all neigbours except that user/item
    
    # We will find k nearest neighbours by using slam::simple_triplet_matrix and data.table.
    # Save dim and dimnames, in order to later reconstruct from simple_triplet_matrix. 
    dims_old <- similarities@Dim
    dimnames_old <- similarities@Dimnames
    similarities <- as.simple_triplet_matrix(similarities)
    datatable <- data.table(similarities$i, similarities$j, similarities$v)
    names(datatable) <- c("row", "column", "rating")

    # Function that finds k-th largest value in a vector.
    kthMax <- function(vector, k){
      if(length(vector) <= k) min(vector)
      else{
        sort(vector, partial = length(vector) - (k-1))[length(vector) - (k-1)]
      }
    }

    kthMaxes <- datatable[, kthMax(rating, k), by = column]
    names(kthMaxes) <- c("column", "kthMax")
    datatable <- merge(datatable, kthMaxes, by="column")
    datatable <- datatable[datatable$rating >= datatable$kthMax, ]

    similarities <- as(sparseMatrix(i = datatable$row, j = datatable$column, x = datatable$rating, dims = dims_old, dimnames = dimnames_old), "dgCMatrix")
  }
  
  similarities
}

#' This function implements memory-based collaborative filtering (neighbors method - nm) and calculates rating predictions. 
#' It divides matrix into parts and calcualtes predictions for each part iteratively.
#' This can be useful in case matrices are large and can not fit into memory.
#'
#' @param ratings_matrix (dgCMatrix) Matrix of known ratings. In case alg_method=="ubcf" it should be IU matrix (items are rows, users are columns).
#' In case In case alg_method=="ibcf" it should be UI matrix.
#' @param predictions_indices Indices of cells in ratings_matrix for which we should calculate predictions.
#' @param alg_method (string) "ubcf" or "ibcf"
#' @param normalization (logical) Whether to perform normalization. Currenlty only "center" normalization is supported (subtracting user's mean from ratings).
#' This step currenlty uses {recommendlab} implementation for normalization.
#' @param similarity_metric Function used to calculate similarities. It has to accept two matrices (dgcMatrix) and  calculate similarities between columns.
#' @param k (integer) Number of largest similarites to keep, per column (k nearest neighbours approach).
#' @param make_positive_similarities (logical) Whether all similarities should be modified by a factor in order to have only positive similarities.
#' @param rowchunk_size (integer) How many rows of rating matrix to consider in one iteration. This can be uesful if matrices are large and we want to perform calculations partially.
#' In case we want to cover all rows at once, set this parameter to be >= total number of rows in \code{ratings_matrix}.
#' @param columnchunk_size (integer) How many columns of similarity matrix to consider in one iteration. This can be uesful if matrices are large and we want to perform calculations partially.
#' In case we want to cover all columns at once, set this parameter to be >= total number of columns in \code{ratings_matrix}.
#' @returns Predictions matrix.
#' @note Returned predictions matrix may not contain predictions for all \code{predictions_indices}. This is because of CF algorithm itself 
#' (in case there are no similar users/items which can be used to find a prediction, for example)
#' @note required Matrix, recommenderlab, slam, data.table
predict_nm <- function(ratings_matrix, predictions_indices, alg_method, normalization, similarity_metric, k, make_positive_similarities, rowchunk_size, columnchunk_size){
  
  if(normalization){
    # Currently, we always use center normalization and apply it per users (subtracting user averages).
    if(alg_method == "ubcf") ratings_matrix <- normalize(as(ratings_matrix, "realRatingMatrix"), method = "center", row = FALSE)
    if(alg_method == "ibcf") ratings_matrix <- normalize(as(ratings_matrix, "realRatingMatrix"), method = "center", row = TRUE)
    ratings_matrix@data@x[ratings_matrix@data@x == 0] <- 0.000001 # Prevent droping zeros obtained after applying normalization.
    normalization_info <- ratings_matrix@normalize
    ratings_matrix <- as(ratings_matrix, "dgCMatrix")
  }
  
  # Create initial empty predictions matrix.
  predictions_matrix <- as(sparseMatrix(i = c(), j = c(), dims = ratings_matrix@Dim, dimnames = ratings_matrix@Dimnames), "dgCMatrix")
  
  # Number of splits per rows and columns. 
  num_row_splits <- ceiling(nrow(ratings_matrix)/rowchunk_size)
  num_column_splits <- ceiling(ncol(ratings_matrix)/columnchunk_size) 
  
  # Iterate over columns first, so that each chunk of similarities is calcualated only once.
  for(i in 1:num_column_splits){
    
    start_column <- columnchunk_size * (i-1) + 1 # Start column for the current chunk.
    end_column <- columnchunk_size * i # End column for the current chunk.
    if(ncol(ratings_matrix) < end_column){
      end_column <- ncol(ratings_matrix)
    }
    
    columns_to_consider <- intersect(start_column:end_column, predictions_indices[, 2])
    if(length(columns_to_consider) == 0) next
    
    # Set names of rows and columns to be numbers (indices). 
    # This way similarities and part_predictions, calculated in next steps, will use these names.
    ratings_matrix@Dimnames[[1]] <- as.character(1:nrow(ratings_matrix))
    ratings_matrix@Dimnames[[2]] <- as.character(1:ncol(ratings_matrix))
    
    similarities <- find_similarities(ratings_matrix, columns_to_consider, similarity_metric, make_positive_similarities, k)
    
    for(j in 1:num_row_splits){
      
      start_row <- rowchunk_size * (j-1) + 1 # Start row for the current chunk.
      end_row <- rowchunk_size * j # End row for the current chunk.
      if(nrow(ratings_matrix) < end_row){
        end_row <- nrow(ratings_matrix)
      }
      
      rows_to_consider <- intersect(start_row:end_row, predictions_indices[, 1])
      if(length(rows_to_consider) == 0) next
      
      # print(paste("Current chunk: ", start_row, end_row, start_column, end_column, sep = ","))
      part_predictions <- calculate_predictions(ratings_matrix[rows_to_consider, , drop = FALSE], similarities) # drop = FALSE because of the case when we have only one row, make it dgCMatrix.
      
      # Fill predictions matrix with predictions calculated in this iteration.
      predictions_indices_to_consider <- subset(predictions_indices, predictions_indices[, 1] %in% rows_to_consider & predictions_indices[, 2] %in% columns_to_consider)
      predictions_matrix <- add_predictions_to_prediction_matrix(predictions_matrix, part_predictions, predictions_indices_to_consider)
    }
    
  }
  
  if(normalization){
    temp <- as(predictions_matrix, "realRatingMatrix")
    temp@normalize <- normalization_info
    predictions_matrix <- denormalize(temp)
    predictions_matrix <- as(predictions_matrix, "dgCMatrix")
  }
  
  predictions_matrix
}


#' This function calculates predictions for given parameters, i.e. latent factors. 
#' It calculates predictions part-by-part in parallel. This can be useful in case matrices are large and can not fit memory at once.
#'
#' @param ratings_matrix (dgCMatrix) Matrix of known ratings. It is IU matrix (items should be rows).
#' @param Theta Matrix of latent factors for items.
#' @param X Matrix of latent factors for users.
#' @param rowchunk_size (integer) When calculating predictions, how many rows of ratings matrix to cover in one part.
#' This can be useful if rating matrix is large and we want to calculate predictions partially.
#' In case we want to cover all rows at once, set this parameter to be >= total number of rows in \code{ratings_matrix}.
#' @param columnchunk_size (integer) HWhen calculating predictions, how many columns of ratings matrix to cover in one part.
#' In case we want to cover all columns at once, set this parameter to be >= total number of columns in \code{ratings_matrix}.
#' @returns Predictions matrix.
#' @note required Matrix, recommenderlab, foreach
predict_mf <- function(ratings_matrix, Theta, X, rowchunk_size, columnchunk_size){
  
  # Test if dimensions of ratings_matrix, Theta, X fit ...
  
  # Take indices of ratings_matrix cells that are not empty (0 is a code for missing values in ratings_matrix).
  predictions_indices <- as.data.frame(which(ratings_matrix != 0, arr.ind = T))
  
  X = as(X, "dgCMatrix")
  X@Dimnames[[1]] <- as.character(1:nrow(X))
  X@Dimnames[[2]] <- as.character(1:ncol(X))
  Theta = as(Theta, "dgCMatrix")
  Theta@Dimnames[[1]] <- as.character(1:nrow(Theta))
  Theta@Dimnames[[2]] <- as.character(1:ncol(Theta))
  
  # Create initial empty predictions matrix.
  predictions_matrix <- as(sparseMatrix(i = c(), j = c(), dims = ratings_matrix@Dim, dimnames = ratings_matrix@Dimnames), "dgCMatrix")
  # Calculate chunk indices (for each chunk find (start_row, end_row, start_column, end_column)).
  all_quarters <- generate_chunk_quarters(nrow(ratings_matrix), ncol(ratings_matrix), rowchunk_size, columnchunk_size)
  
  foreach(chunk_indices = all_quarters) %dopar% {
    start_row <- chunk_indices[1]
    end_row <- chunk_indices[2]
    start_column <- chunk_indices[3]
    end_column <- chunk_indices[4]
    
    # print(paste("Current chunk: ", start_row, end_row, start_column, end_column, sep = ","))
    # Calculate predictions for this chunk.
    part_predictions <- X[start_row:end_row, ] %*% t(Theta)[, start_column:end_column]
    # Fill predictions matrix with predictions calculated in this iteration.
    predictions_indices_to_consider <- subset(predictions_indices, predictions_indices[, 1] %in% c(start_row:end_row) & predictions_indices[, 2] %in% c(start_column:end_column))
    predictions_matrix <- add_predictions_to_prediction_matrix(predictions_matrix, part_predictions, predictions_indices_to_consider)
  }
  
  predictions_matrix
}


