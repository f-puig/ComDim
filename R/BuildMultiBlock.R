BuildMultiBlock <- function(newBlock = newBlock, newVars = NULL, newSamples = NULL, growingMB = NULL, batches = NULL, equalSampleNumber = T, ...) {
  

#' BuildMultiBlock
#'
#' Creates a multi-block structure compatible with ComDim_PCA_2020.R
#' @param newBlock The data to be incorporated into the multi-block structure.
#' @param newVars The variable names (optional)
#' @param newSamples The sample names. The sample replicates should have the same name (optional).
#' @param growingMB The name of the multi-block structure it this already exists (optional).
#' @param batches A vector indicating the batch number for every sample (optional).
#' @param equalSampleNumber Logical to indicate wether the number of samples is the same (T) or not (F) in all the blocks.
#' @param ... Additional objects descriptive of newBlock can be introduced, as long as they have the same length as the sample size.
#' @return The multi-block
#' @examples 
#'  b1 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
#'  b2 = matrix(rnorm(800),10,80) # 10 rows and 80 columns
#'  b3 = matrix(rnorm(700),10,70) # 10 rows and 70 columns
#'  # Build multi-block by adding one data block at a time:
#'  mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:10))
#'  mb <- BuildMultiBlock(b2, growingMB = mb)
#'  mb <- BuildMultiBlock(b3, growingMB = mb)
#' @export  

  formattedMB <- list()
  
  if(!(is.matrix(newBlock))){
    stop('newBlock is not a matrix.')
  }
  
  if(length(as.character(substitute(newBlock))) > 1){
    print('The provided newBlock is contained in a list, while it must be a stand-alone matrix object.')
    stop('Newblock must be a matrix object.')
  }
  
  if(length(list(...)) != 0){
    add_param <- list(...) # The additional parameters (ex: Classes, Levels)
  } else {
    add_param <- NULL
  }
  
  if (is.null(growingMB)) {
    growingMB <- list()
    block <- 1
  } else {
    block <- length(growingMB)+1
    if(as.character(substitute(newBlock)) %in% names(growingMB)){ # Don't execute the code if the newBlock name already exists in growingMB
      stop(sprintf('There already exists a block named %s in the multi-block.', as.character(substitute(newBlock))))
    }
  }
  
  # Check newBlock
  if(equalSampleNumber){
    if ((block != 1) && (nrow(growingMB[[1]]$Data) != nrow(newBlock))) {
      stop('The number of samples in newBlock is incorrect.')
    }
  }
  #formattedMB$Name <- as.character(substitute(newBlock)) # The name of the block
  formattedMB$Data <- newBlock
  

  # Check newVars
  if (is.null(newVars)) {
    formattedMB$Variables <- seq(1,ncol(newBlock),1)
  } else {
    if (ncol(newBlock) == length(newVars)) {
      formattedMB$Variables <- newVars
    } else {
      print(sprintf('Please use a newVars object with %s items.', ncol(newBlock)))
      stop('The length of newVars is incorrect.')
    }
  }
  
  # Check newSamples
  if (is.null(newSamples)) {
    formattedMB$Samples <- seq(1,nrow(newBlock),1)
  } else {
    if (nrow(newBlock) == length(newSamples)) {
      formattedMB$Samples <- newSamples
    } else {
      print(sprintf('Please use a newSamples object with %s items.', nrow(newBlock)))
      stop('The length of newSamples is incorrect.')
    }
  }
  
  # Check Batches
  if (!(is.null(batches))) {
    if (!(is.numeric(batches))){
      print('batches should be a numeric vector.')
      stop('The batches object is incorrect.')
    } else {
      if (nrow(newBlock) == length(batches)) {
        formattedMB$Batch <- batches
      } else {
        print(sprintf('Please use a batches object with %s items.', nrow(newBlock)))
        stop('The length of batches is incorrect.')
      }
    }
  } else {
    formattedMB$Batch <- NULL
  }
  
  
  if (!(is.null(add_param))){
    for(i in 1:length(add_param)){
      if(length(add_param[[i]]) == nrow(newBlock)){
        formattedMB[[names(add_param[i])]] <- add_param[[i]]
      } else {
        warning(sprintf('The additional object named %s has an incorrect size.', names(add_param[i])))
        print(sprintf('The parameter %s was ignored.', names(add_param[i])))
      }
    }
  }
  
  # Add the new formattedMB to growingMB
  #growingMB[[block]]<- formattedMB
  growingMB[[ as.character(substitute(newBlock)) ]]<- formattedMB
  
  return(growingMB)
}