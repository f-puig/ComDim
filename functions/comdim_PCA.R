if (!requireNamespace("pracma", quietly = TRUE))
  install.packages("pracma")

comdim_PCA <- function(data = data, ndim = NULL, normalise = TRUE, threshold = 1e-10, loquace = FALSE, output = 'TPL', CompMethod = 'Normal', Partitions = 1) {

#comdim	- Finding common dimensions in multitable data (saisiR format)
#
# With data compression methods 
#
#
# USAGE :
#---------------
# results <- comdim_PCA(data = data, ndim = NULL)
#
# INPUT :
#---------------
# data : vector of saisiR files (the numbers of rows in each table must be equal).
#        Multitable data can be converted to saisiR format with the BuildMultiBlock function.
#        If the data contains replicate samples, the resulting multitable from BuildMultiBlock
#        can be split into replicate blocks with the splitRB function.
#
# Optional inputs:
# ndim : Number of common dimensions
# normalise : FALSE == no; TRUE == yes (default)
# threshold : If the "difference of fit" < threshold (1e-10 as default)
# loquace : TRUE == displays calculation times;  FALSE == no (default)
# output: 
# if Output== 'TPL' , all extra parameters are calculated (default)
# if Output=='[]', no extra parameters are calculated
# 'T' : Local Scores
# 'P' : Scaled Loadings
# 'L' : Unscaled Loadings

# # # # To accelerate the ComDim calculations with big matrices
# CompMethod = 'Normal' (default), 'Kernel', 'PCT', 'Tall' or 'Wide'
# Partitions = number of partition (1 as default) if CompMethod is 'Tall' or 'Wide'

# OUTPUT :
#-----------------
# res_calib with fields:
# Q : Global scores_calib (nrow x ndim)
# P : Scaled Global loadings calculated from Q and tables (nvars x ndim)
# P_Loc : Structure with ntable of Scaled Local loadings in matrices of (local vars x ndim)
# T_Loc : Structure with ntable of Local scores_calib calculated from local Loadings and tables (nrow x ndim)
# Lx : Unscaled Global loadings for X
# Lx_Loc : Unscaled Local loadings for X
# saliences : weight of the original tables in each dimension (ntable x ndim)
# Sum_saliences_Tab - Sum of Saliences for each Table in a Dimension
# Sum_saliences_Dim - Sum of Saliences for each Dimension for a Table
# b : regression coefficients between Local and Global Scores
# explained : percentage explanation given by each dimension (1 x ndim)
#
# NormMB : ntable norms of the tables (1 x 1)
# MeanMB : ntable means of tables (1 x nvars)
#
# SingVal : vector of singular values (1 x ndim)
# 
#
#
# REFERENCES :
#-----------------
# Method published by:
# E.M. Qannari, I. Wakeling, P. Courcoux and H. J. H. MacFie
# in Food quality and Preference 11 (2000) 151-154
#
# Example application
# Chemometric Tools to Highlight Possible Migration of Compounds
# from Packaging to Sunflower Oils
# J. Maalouly, N. Hayeck, A. Kassouf,D. N. Rutledge and V. Ducruet
# J. Agric. Food Chem. 2013, 61, 10565?10573
#
# Code translated from the following MATLAB function: comdim_PCA_2020.m
#
# INITIALISATION
  
  library(pracma)
  
  progress_bar = txtProgressBar(min=0, max=80, style = 3, char = "=")
  
  start_ini_time <- Sys.time()  # To start counting calculation time for the initialization
  
  ntable = length(data) # number of blocks
  
  # Check that data is a suitable multi-block structure
  properties <- c('Data', 'Variables', 'Samples')
  
  print('Checking each block individually.')
  give_error <- 0
  
  sample_number <- as.vector(NULL)
  variable_number <- as.vector(NULL)
  for (i in 1:ntable){
    
    sample_number[i] <- length(data[[1]]$Samples)
    variable_number[i] <- length(data[[i]]$Variables)
    
    if(sum(properties %in% names(data[[i]])) != 3){ # Checking if a field is missing
      print(sprintf('Block %s has missing fields.', i))
      give_error <- 1
    }
    
    if (nrow(data[[i]]$Data) != length(data[[i]]$Samples)) {
      print(sprintf('Block %s has incorrect number of samples.', i))
      give_error <- 1
    }
    
    if (ncol(data[[i]]$Data) != length(data[[i]]$Variables)) {
      print(sprintf('Block %s has incorrect number of variables.', i))
      give_error <- 1
    }  
  }

  print('Checking sample correspondence across blocks.')

  if(length(unique(sample_number)) != 1){
    print('Not all the blocks have the same number of samples.')
    give_eror <- 1
  }
  nrowMB <- sample_number[1] # I cannot use nrow because this is already a function in R.
  
  if (give_error) {
    stop('The data is not ready for ComDim.')
  } else {
    print('The data can be used for ComDim.')
  }
  
  
  if(is.null(ndim)){
    ndim <- ntable # If the number of components is not defined, 
                   # the number of components to extract is equal to the number of blocks. 
  }
  
  pieceBar <- 4+2*ntable+ndim # Number of updates in the progress bar.
  pieceBar <- 80/pieceBar
  total_progress <- pieceBar
  
  DimLabels <- paste0('CC',1:ndim)   # One label per component.
  TableLabels<- names(data) # One label per block.
  
  isT <- grepl(output, 'T')
  isL <- grepl(output, 'L')
  isP <- grepl(output, 'P')
  isLP <- c(isP,isL)
  isTLP <- c(isT,isP,isL)
  
  if(CompMethod == 'Tall' && any(variable_number > 1000)){
    print("The method 'Tall' is not recommeded for large blocks (> 1000 columns).")
    print("Do you still want to continue with the analysis?")
    print("You can change now the CompMethod if you prefer: (yes/no)")
    ans1 <- tolower(readline("Type your choice: (y/n)"))
    if (ans1 == 'y' || ans1 == 'yes') {
      print("Which method do you want to use instead (Normal/Kernel/PCT/Wide)?")
      ans2 <- tolower(readline("Type your choice: (n/k/p/w)"))
      if (ans2 == 'n' || ans2 == 'normal'){
        CompMethod = 'Normal'
      } else if(ans2 == 'k' || ans2 == 'kernel'){
        CompMethod = 'Kernel'
      } else if(ans2 == 'p' || ans2 == 'pct'){
        CompMethod = 'PCT'
      } else if(ans2 == 'w' || ans2 == 'wide'){
        CompMethod = 'Wide'
      } else {
        stop("Wrong answer typed.")
      }
    } else if (ans1 == 'n' || ans1 == 'no'){
      print("The method 'Tall' will be used")
    } else {
      stop("Wrong answer typed.")
    }
  }
  
  end_ini_time <- Sys.time() # To end the count of the calculation time.
  
  if (loquace) {
    print(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time)*1000))
  }
 
  setTxtProgressBar(progress_bar, value = total_progress)
  
# NORMALISATION  
  
  #start_norm_time <- Sys.time()  # To start counting calculation time for the initialization
  
  X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
 
  res_calib <- list()
  temp_tabCalib <- list()
  s_n <- list()
  res_calib$SingVal$Data <-as.vector(NULL)
  res_calib$NormMB$Data <-as.vector(NULL)
  res_calib$MeanMB$Data <-list()
  
  for (i in 1:ntable) {
      
    if (normalise){
      # Normalise original blocks
      res_calib$MeanMB$Data[[TableLabels[i]]] <- colMeans(data[[i]]$Data)
      res_calib$MeanMB$Variables[[TableLabels[i]]] <- data[[i]]$Variables
      #res_calib$MeanMB$i <- TableLabels 
      
      X_mean <- data[[i]]$Data-as.matrix(x = rep(1,nrowMB), ncol = 1, nrow = nrowMB) %*% res_calib$MeanMB$Data[[i]]
      XX <- X_mean*X_mean
      
      Norm_X <- sqrt(sum(XX))
      X_Normed <- X_mean/Norm_X
      
      res_calib$NormMB$Data[i] <- Norm_X
      
      temp_tabCalib[[i]] <- X_Normed
      s_n [[i]] <- X_Normed
      
    } else {
      res_calib$NormMB$Data[[TableLabels[i]]] <- rep(1,length(data[[i]]$Variables))
      res_calib$MeanMB$Variables <- data[[i]]$Variables
      #res_calib$MeanMB$i <- TableLabels 
      
      temp_tabCalib[[i]] <- data[[i]]$Data
      s_n [[i]] <- data[[i]]$Data
    }
    
    if (i==1){
      X_mat[,1:variable_number[1]] = data[[i]]$Data
      Xnorm_mat[,1:variable_number[1]] = temp_tabCalib[[i]]
    } else {
      beg <- sum(variable_number[1:(i-1)])+1
      ending <- sum(variable_number[1:i])
      X_mat[,beg:ending] = data[[i]]$Data
      Xnorm_mat[,beg:ending] = temp_tabCalib[[i]]
    }
    
    norm_comdim <- Sys.time()
    if(loquace){
      
      print(sprintf("Normalization of block %s finished after : %s millisecs", i, (norm_comdim - start_ini_time)*1000))
      
    }
    total_progress <- total_progress + pieceBar
    setTxtProgressBar(progress_bar, value = total_progress)

  }
  
  #res_calib$NormMB$i <- 'Norm'
  names(res_calib$NormMB$Data) <- TableLabels
  
  nR <- nrow(Xnorm_mat)
  nC <- ncol(Xnorm_mat)

# COMPRESSION
  
  s_r <- Compress_Data_2020(s_n, CompMethod, Partitions)
  
  end_comp_time <- Sys.time()
  
  if(loquace){
    
    print(sprintf("Compression finished after : %s millisecs", (end_comp_time - start_ini_time)*1000))
  
  }
  
  total_progress <- total_progress + pieceBar*ntable
  setTxtProgressBar(progress_bar, value = total_progress)
  
  ## DO ComDim
  #ini_comdim <- Sys.time()
  
  saliences <- list()
  saliences$Data <- matrix(, ncol = ndim, nrow  = ntable)
  Q <- list()
  Q$Data <- matrix(, ncol = ndim, nrow  = nrowMB)
  unexplained <- ntable; # Warning!: true only if the tables were set to norm=1
  varexp <- as.vector(NULL)
  
  for(dims in 1:ndim){
    previousfit <- 10000
    deltafit <- 1000000
    lambda <- matrix(rep(1,ntable),nrow=ntable,ncol=1)
  
    qini <- as.vector(s_r[[1]][,1])
    qini <- qini/sqrt(as.vector(qini %*% qini)) #random initialization of t + set to unit length
    qold <- 100
  
    iters <- 0
    while(norm(qold-qini, "2") > threshold && iters < 100){
      iters <- iters + 1
      qold <- qini
      q <- 0
    
      W <- matrix(, nrow = nrowMB, ncol = 0)
      for(j in 1:ntable){
        W <- cbind(W,sqrt(lambda[j])*s_r[[j]])
      }
    
      usv <- svd(W)
      qtmp <- usv$u[,1]
    
      for(j in 1:ntable){
        #takes each table into account for lambda after PCA
        lambda[j] <- t(qtmp)%*%(s_r[[j]]%*%t(s_r[[j]]))%*%qtmp
        q <- q + lambda[j]*qtmp
      }
    
      q <- q/sqrt(as.vector(t(q)%*%q)) #standardizes t
      if (abs(min(q)) > abs(max(q))){
        q <- -q
      }
      qini <- q
    } #deltafit>threshold
  
    saliences$Data[,dims] <- lambda
    Q$Data[,dims] <- q

  
    #'SingVal' is square of the SVD singular value
    # because here it is based on X, not X.X'
    # Does NOT work for Kernel & Wide & Tall
    if(CompMethod == 'Normal' || CompMethod == 'PCT'){
      res_calib$SingVal$Data[dims] <- usv$d[1]*usv$d[1] #from SVD on W
      varexp[dims] <- res_calib$SingVal$Data[dims]^2
    } else {
      varexp[dims] <- 0 #This will be overwritten at the end
      res_calib$SingVal$Data[dims] <- 0 #This will be overwritted at the end
    }
    

    # Deflate blocks
    aux <- diag(nrowMB)-as.matrix(q)%*%t(as.matrix(q))
    for(j in 1:ntable){
      s_r[[j]] <- aux%*%s_r[[j]]
    }
    
    iter_comdim <- Sys.time()
    if(loquace){
      
      print(sprintf("Component %s determined after : %s millisecs", dims, (iter_comdim - start_ini_time)*1000))
      
    }
    total_progress <- total_progress + pieceBar
    setTxtProgressBar(progress_bar, value = total_progress)
  
    }
  
  res_calib$SingVal$i <- 'Singular value'
  names(res_calib$SingVal$Data) <- DimLabels #Dimensions
  
  Q$i <- data[[1]]$Samples
  Q$v <- DimLabels
  colnames(Q$Data) <- DimLabels
  rownames(Q$Data) <- data[[1]]$Samples
  res_calib$Q <- Q
  
  ## Adding metadata to Q --> Metadata extracted from the first block
  if(length(names(data[[1]])) > 4){ # In case there are additional fields
    newFields <- setdiff(names(data[[1]]),c(properties,'Batch'))
    print('Adding metadata to Q. Metadata copied from the first block.')  
    if (length(newFields) != 0){
      for(fields in 1:length(newFields)) {
        if(length(data[[1]][[newFields[fields]]]) == nrowMB) {
          res_calib$Q[[newFields[fields]]] <- data[[1]][[newFields[fields]]]
        } else {
          warning(sprintf('The list %s has an incorrect size.', newFields[fields]))
          print(sprintf('The list %s was ignored.', newFields[fields]))
        }
      }
    }
  }
  
  end_comdim <- Sys.time()
  if (loquace) {
    print(sprintf("Scores finished after : %s millisecs", (end_comdim - start_ini_time)*1000))
  }
  total_progress <- total_progress + pieceBar
  setTxtProgressBar(progress_bar, value = total_progress)
  
  meanW <- as.vector(NULL)
  for (j in 1:ntable){
    meanW <- c(meanW,mean(s_r[[j ]]))
  }
  #res_calib$s_n <- s_n
  #res_calib$s_r <- s_r
  rm(s_n)
  rm(s_r)
  
  ##
  explained <- list()
  explained$Data <- varexp/sum(varexp)*100
  names(explained$Data) <- DimLabels
  explained$i <-'% explained'
  #explained$v <- DimLabels #Dimensions
  res_calib$explained <- explained
  rm(varexp)
  rm(explained)
  
  ## Calculate Sums of saliences for each Dim
  Sum_saliences_Dim <- list()
  Sum_saliences_Dim$Data <- colSums(saliences$Data)
  names(Sum_saliences_Dim$Data) <- DimLabels
  Sum_saliences_Dim$i <-'Sum Dim Saliences'
  #Sum_saliences_Dim$v <- DimLabels #Dimensions
  
  res_calib$Sum_saliences_Dim <- Sum_saliences_Dim
  rm(Sum_saliences_Dim)
  
  #Calculate Sums of saliences for each Table
  Sum_saliences_Tab <- list()
  Sum_saliences_Tab$Data <- rowSums(saliences$Data)
  Sum_saliences_Tab$i <- 'Sum Tab Saliences'
  names(Sum_saliences_Tab$Data) <- TableLabels
  #Sum_saliences_Tab$v <- TableLabels #tables
  
  res_calib$Sum_saliences_Tab <- Sum_saliences_Tab
  rm(Sum_saliences_Tab)
  
  saliences$i <- TableLabels #tables
  saliences$v <- DimLabels #Dimensions
  res_calib$saliences <- saliences
  colnames(res_calib$saliences$Data) <- DimLabels
  rownames(res_calib$saliences$Data) <- TableLabels
  
  ##Calculate Normalised concatenated Xs ('Calib') from col
  # Calculate concatenated CD loadings
  # Reorganise Loadings - 1 matrix / LV
  
  ## If Output = NULL, nothing else is calculated
  if(!(is.null(output))){
      L_CD_Vec <- NULL
      L_X_Vec <- NULL
      Calib <- NULL
   
      nCalib <-  nrowMB
  }

  ## Already defined in during ComDim initiallization
  #isT <- grepl(output, 'T')
  #isL <- grepl(output, 'L')
  #isP <- grepl(output, 'P')
  #isLP <- c(isP,isL)
  #isTLP <- c(isT,isP,isL)
  
  # tic
  # Calculate concatenated CD loadings
  # Reorganise Loadings - 1 matrix / LV
  L_CD <- list()
  L_X <- list()
  T_Loc <- list()
  for(i in 1:ntable){ # Prepare lists
    L_CD[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    L_X[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    T_Loc[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = nrowMB)
  }
  b <- matrix(,ncol = ndim, nrow = ntable)
  
  for(j in 1:ndim){
    T_mat <- matrix(, nrow = nrowMB, ncol = 0)
    
    for(i in 1:ntable){
    # Q.d are orthonormal in ICA & PCA

      temp <- t(temp_tabCalib[[i]])%*%Q$Data[,j] #Scaled CD Loadings "locaux" calculés
      
      if  (!(is.null(isP))){
        L_CD[[TableLabels[i]]][,j] <- temp
      }
      
      if  (!(is.null(isL))){
        #Q.d are orthonormal in ICA & PCA
        L_X[[TableLabels[i]]][,j] <- t(data[[i]]$Data)%*%Q$Data[,j] #Unscaled X Loadings "locaux" calculés;
      }
      
      T_mat <- cbind(T_mat,temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp))) #Scores "locaux"];

      if(!(is.null(isT))){
        T_Loc[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp)) #Scores "locaux"
      }
      
      
      # Deflate each temp_tabCalib
      temp_tabCalib[[i]] <- temp_tabCalib[[i]]-Q$Data[,j]%*%t(temp)
    }
    

    # For each CC
    # MLR b-coefficients between Local and Global Scores
    #[b0,b[,j]] <- mlr_DB(T_mat,Q$d[,j],0)
    # b=pinv(X'*X)*X'*Y;
    
    b[,j] <- pinv(t(T_mat)%*%T_mat)%*%t(T_mat)%*%Q$Data[,j]
    #pinv(t(aa$T_mat)%*%aa$T_mat)%*%t(aa$T_mat)%*%aa$Q$Data[,2]
    b0 <- 0
  }
  
  # If Output==[], nothing else is calculated
  if(!(is.null(output))){
    # Calculate Global Loadings
    if(!(is.null(isP))){
      L_CD_Vec <- t(Xnorm_mat)%*%Q$Data # Scaled CD Loadings "globaux"
    }
           
    if(!(is.null(isL))){
      L_X_Vec <- t(X_mat)%*%Q$Data # Unscaled X Loadings "globaux"
    }
  }
  
  #### If Output==[], nothing else is calculated
  load_comdim <- Sys.time()
  if(loquace){
    print(sprintf("Loadings finished after : %s millisecs", (load_comdim - start_ini_time)*1000))
  }
  
  total_progress <- total_progress + pieceBar
  setTxtProgressBar(progress_bar, value = total_progress)
  
  rm(X_mat)
  rm(temp_tabCalib)
  rm(temp)

  ## Output
  end_function <- Sys.time()

  res_calib$b$Data <- b
  colnames(res_calib$b$Data) <- DimLabels
  rownames(res_calib$b$Data) <- TableLabels
  res_calib$b$i <- TableLabels # tables
  res_calib$b$v <- DimLabels  # Dimensions
    
  if(!(is.null(isT))){
    res_calib$T_Loc$i <- res_calib$Q$i # samples
    res_calib$T_Loc$v <- DimLabels # Dimensions
    res_calib$T_Loc$Data <- T_Loc
    
    ## Add metadata to T_Loc
    for(i in 1:ntable){
      if(length(names(data[[i]])) > 3){ # In case there are additional fields
        newFields <- setdiff(names(data[[i]]), properties)
        if (length(newFields) != 0){
          for(fields in 1:length(newFields)) {
            if(length(data[[i]][[newFields[fields]]]) == nrowMB) {
              res_calib$T_Loc[[newFields[fields]]][[TableLabels[i]]] <- data[[i]][[newFields[fields]]]
            } else {
              warning(sprintf('The list %s in block %s has an incorrect size.', newFields[fields], as.character(i)))
              print(sprintf('The list %s from block %s was ignored.', newFields[fields],as.character(i)))
            }
          }
        }
      }
    }
  }
 
  # T_Loc is no longer needed 
  rm(T_Loc)
  rm(T_mat)
    
  if(!(is.null(isP))){
    res_calib$P$Data <- L_CD_Vec
    res_calib$P$i <- t(1:nC) # numbers of all variables;
    res_calib$P$v <- DimLabels # Dimensions
    colnames(res_calib$P$Data) <- DimLabels 
    rownames(res_calib$P$Data) <- t(1:nC)

    res_calib$P_Loc$Data <- L_CD 
    #res_calib$P_Loc$i <- TableLabels # Tables
    res_calib$P_Loc$v <- DimLabels # Dimensions
    for(i in 1:ntable){
      rownames(res_calib$P_Loc$Data[[TableLabels[i]]]) <- data[[i]]$Variables
      colnames(res_calib$P_Loc$Data[[TableLabels[i]]]) <- DimLabels
    }
  }
  
  rm(L_CD_Vec)
  rm(L_CD)
  
  if(!(is.null(isL))){
    res_calib$Lx$i <- t(c(1:nC)) # numbers of all variables;
    res_calib$Lx$v <- DimLabels # Dimensions
    res_calib$Lx$Data <- L_X_Vec
    colnames(res_calib$Lx$Data) <- DimLabels
    rownames(res_calib$Lx$Data) <- t(c(1:nC))
    #res_calib$Lx_Loc$i <- TableLabels # Tables
    res_calib$Lx_Loc$v <- DimLabels # Dimensions
    res_calib$Lx_Loc$Data <- L_X
    for (i in 1:ntable){
      colnames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- DimLabels
      rownames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- data[[i]]$Variables
    }
  }
  rm(L_X_Vec)
  rm(L_X)
    
  ##
  # Singular value calculated from b and saliences 
  # Since :
  # b_T_Q(:,j)=saliences.d(:,j).*saliences.d(:,j)/SingVal.d(1,j);
  if(CompMethod == 'Kernel' || CompMethod == 'Tall' || CompMethod == 'Wide'){
    for(dims in 1:ndim){
      res_calib$SingVal$Data[dims] <- t(b[,dims])%*%((saliences$Data[,dims]^2)/as.vector(t(b[,dims])%*%b[,dims]))
    }
    names(res_calib$SingVal$Data) <- DimLabels
    varexp <- res_calib$SingVal$Data^2
    res_calib$explained$Data <- varexp/sum(varexp)*100
    names(res_calib$explained$Data) <- DimLabels
  }
  
  rm(saliences)
  
  #### If Output==[], nothing else is calculated
  #if(normalise == 1){
  #  res_calib$Xnorm_mat$i <- res_calib$Q$i # samples
  #  res_calib$Xnorm_mat$v <- t(c(1:nC)) # numbers of all variables;
  #  res_calib$Xnorm_mat$Data <- Xnorm_mat
  #}
  #rm(Xnorm_mat)
  #### If Output==[], nothing else is calculated
  
  rm(Q)
  
  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time)
  if(loquace){
    print(sprintf("Analysis finished after : %s millisecs", running_time*1000))
  }
  res_calib$runtime <- running_time  # Total time of analysis
  
  progress_bar = txtProgressBar(min=0, max=80, style = 3, char = "=")
  setTxtProgressBar(progress_bar, value = 80)
  
  close(progress_bar)
  
  return(res_calib) 
}


Compress_Data_2020 <- function(s_n = s_n, CompMethod = CompMethod, Partitions = Partitions){
  # Do Compression
  s_r <- list()
  ntable <- length(s_n)
  nR <- nrow(s_n[[1]])
  
  if(!is.null(CompMethod)){
    if (CompMethod == 'Tall'){
      for(i in 1:ntable) {
        MaxRank <- min(c(nrow(s_n[[i]]),ncol(s_n[[i]])))
        # Tall segmented PCT
        s_r[[i]] <- PCA_Tall_PCT_DNR(Xn = s_n[[i]], PCs = MaxRank, Partitions = Partitions)
      }
    } else if (CompMethod == 'Wide'){
      for(i in 1:ntable){
        Xi <- ColumnsPartition(Xn = s_n[[i]], Partitions = Partitions)
        T_all <-   matrix(, nrow = nrow(s_n[[i]]), ncol = 0)
        for(p in 1:Partitions){
          usv <- svd(Xi[[p]])
          T_in <- usv$u%*%diag(usv$d)
          T_all <- cbind(T_all, T_in)
        }
        s_r[[i]] <- T_all
        #print(s_r[[i]][1:3,1:4])
      }
    } else if (CompMethod == 'Kernel'){
      for(i in 1:ntable){
        Kern <- (s_n[[i]]%*%t(s_n[[i]]))/nR
        usv <- svd(Kern)
        s_r[[i]] <- usv$u%*%diag(usv$d)
      }
    } else if (CompMethod == 'PCT'){
      for(i in 1:ntable){
        usv <- svd(s_n[[i]])
        s_r[[i]] <- usv$u%*%diag(usv$d)
      }
    } else if (CompMethod == 'Normal'){
      for(i in 1:ntable){
        s_r[[i]] <- s_n[[i]]
      }
    }
  }
  return(s_r)
}


PCA_Tall_PCT_DNR <-function(Xn = Xn, PCs = PCs, Partitions = Partitions){
  W <- RowsPartition(Xn, Partitions)
  cols <- ncol(Xn)
  
  D <- t(W[[1]])%*%W[[1]]
  
  if(Partitions > 1) {
    for(par in 2:Partitions){
      D <- D + (t(W[[par]])%*%W[[par]])
    }
  }

  vs <- eigen(D)
  vs_flipped <- vs$vectors[seq(length(vs$values),1,-1),]

  Tm <- matrix(, nrow=nrow(Xn), ncol=0)
  Taux <- matrix(, nrow=0, ncol = PCs)
  
  for(par in 1:Partitions){
    to <- W[[par]]%*%vs_flipped
    
    Taux <- rbind(Taux, to[,1:PCs])
  }

  Tm <- cbind(Tm, Taux)
  return(Tm)
}


ColumnsPartition <- function(Xn = Xn, Partitions = Partitions){
  cols <- ncol(Xn)
  stride <- floor(cols/Partitions)
  remaining <- cols %% Partitions
  
  count <- 1
  W <-list()
  
  for(i in 1:Partitions){
    step <- count + stride - 1
    if(i == Partitions){
      step <- step + remaining
    }
    W[[i]] <- Xn[,count:step]
    count <- count + stride
  }
  return(W)
}

RowsPartition <- function(Xn = Xn, Partitions = Partitions){
  rows <- nrow(Xn)
  stride <- floor(rows/Partitions)
  remaining <- rows %% Partitions
  
  count <- 1
  W <- list()
  
  for(i in 1:Partitions){
    step <- count + stride - 1
    if(i == Partitions){
      step <- step + remaining
    }
    W[[i]] <- Xn[count:step,]
    count <- count + stride
  }
  return(W)
}