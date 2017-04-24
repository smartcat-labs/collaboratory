#' For given matrix dimensions, calculate chunk indices (for each chunk find (start_row, end_row, start_column, end_column))
#'
#' @param nrows Number of rows in the Matrix
#' @param ncols Number of columns in the Matrix
#' @param rowchunk_size Number of rows per chunk
#' @param rowchunk_size Number of columns per chunk
generate_chunk_quarters <- function(nrows, ncols, rowchunk_size, columnchunk_size){
  
  # Number of splits per rows and columns.
  num_row_splits <- ceiling(nrows/rowchunk_size)
  num_column_splits <- ceiling(ncols/columnchunk_size)
  
  all_quarters <- list()
  counter <- 1
  
  for(i in 1:num_column_splits){
    start_column <- columnchunk_size * (i-1) + 1 # Start column for the current chunk.
    end_column <- columnchunk_size * i # End column for the current chunk.
    if(ncols < end_column){
      end_column <- ncols
    }
    for(j in 1:num_row_splits){
      start_row <- rowchunk_size * (j-1) + 1 # Start row for the current chunk.
      end_row <- rowchunk_size * j # End row for the current chunk.
      if(nrows < end_row){
        end_row <- nrows
      }
      all_quarters[[counter]] <- c(start_row, end_row, start_column, end_column)
      counter <- counter + 1
    }
  }
  
  all_quarters
}


#' @param predictions_matrix (dgCMatrix) Matrix where all predictions are stored
#' @param part_predictions (dgCMatrix) Part of all predictions. It has names of rows and columns (Dimnames) corresponding to real indices in predictions matrix.
#' It contains predictions only for those rows and columns that exist in \code{predictions_matrix_indices}.
#' @param predictions_matrix_indices Indices of predictions matrix where predictions should be stored.
#' @returns Predictions matrix with added predictions (dgCMatrix)
add_predictions_to_prediction_matrix <- function(predictions_matrix, part_predictions, predictions_matrix_indices){
  
  require(Matrix)
  # require(recommenderlab)
  
  row_names <- as.integer(unlist(part_predictions@Dimnames[1])) # Real row indices from predictions matrix.
  columns_names <- as.integer(unlist(part_predictions@Dimnames[2]))
  row_info <- cbind(row_name = row_names, row_index = 1:length(row_names)) # row_index = row indices from part_predictions.
  column_info <- cbind(column_name = columns_names, column_index = 1:length(columns_names))
  
  all_indices <- predictions_matrix_indices
  colnames(all_indices) <- c("row_name", "column_name")
  all_indices <- merge(all_indices, row_info)
  all_indices <- merge(all_indices, column_info)
  
  predictions_matrix_indices <- all_indices[, c("row_name", "column_name")]
  part_matrix_indices <- all_indices[, c("row_index", "column_index")]
  if(nrow(predictions_matrix_indices) > 0){
    predictions_matrix[as.matrix(predictions_matrix_indices)] <- part_predictions[as.matrix(part_matrix_indices)]
  }
  
  predictions_matrix
}
