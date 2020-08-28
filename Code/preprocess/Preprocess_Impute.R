############################################################
# This script is for data transformation and missing value imputation


setwd("~/Desktop/Projects/WQE/Data/Intensity/")

# load data
strains <- c('c57', 'balbc', 'dba', 'fvb', 'aj', 'cej')
treat <- c('ctrl', 'iso')
timepoints <- c(0,1,  3,  5,  7, 10, 14)
pro_lists <- list()
counts <- 0
for(strain in strains){
  for(tre in treat){
    counts <- counts + 1
    dat <- read.csv(paste0(strain,'_', tre, '_intensity.txt'), 
                    sep = '\t', header = F, row.names = 1)
    print(max(dat))
    dat[dat==0] <- NaN
    
    count <- apply(dat,1,function(x){sum(is.na(x))})
    pro_list <- rownames(dat)[which(count<=3)]
    pro_lists[[counts]] <- pro_list
    assign(paste0(strain,'_',tre,'.pro'), as.character(pro_list))
}}

# unique proteins that has more than 4 values
print(length(Reduce(unique, pro_lists)))
# proteins has more than 4 values in all conditions
pro_list <- Reduce(intersect, pro_lists)

# spline imputation
preprocess_spline <- function(dat,
                              timepoints=NULL,
                              dof=c("cv","cv.global"),
                              center.dat = TRUE,
                              scale.dat = FALSE,
                              verbose = FALSE,
                              seed = NULL,
                              ...) {
  
  if (is.null(seed)) set.seed(seed)
  m <- nrow(dat)
  n <- ncol(dat)
  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))
  
  # sanity check
  if(n != length(timepoints)) stop("denoise_spline: The number of time points must match the number of columns in the input data.")
  
  dat.denoise = matrix(0, nrow=m, ncol=n)
  df <- vector("numeric",m)
  if(dof == "cv") {
    if(verbose) message("\n\r Individual degrees of freedom via cross validation:")
    for(i in 1:m) {
      if(verbose) cat(i)
      if(sum(is.na(dat[i,]))>3) {dat.denoise[i,] <- dat[i,]}
      else{
        #fitting smoothing splines using smooth.spline(X,Y,df=...)
        tpi <- data.frame(x=timepoints, y=dat[i,])
        spline.fit <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,cv = TRUE))
        df[i] <- spline.fit$df
        dat.denoise[i,] <- with(tpi, predict(spline.fit, x)$y)}
    }
  } else if(dof == "cv.global") {
    if(verbose) message("\n\r Cross Validation:")
    for(i in 1:m) {
      if(verbose) cat(i)
      tpi <- data.frame(x=timepoints, y=dat[i,])
      df[i] <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,cv = TRUE)$df)
    }
    if(verbose)  message("\n\r Global degree of freedom:")
    for(i in 1:m) {
      if(verbose) cat(i)
      tpi <- data.frame(x=timepoints, y=dat[i,])
      globaldf = median(df, na.rm = TRUE)
      spline.fit <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,df=globaldf))
      dat.denoise[i,] <- with(tpi, predict(spline.fit, x)$y)
    }
    if(verbose) {
      message(paste("The median DoF is used for all variables:",globaldf))
      hist(df,20,col="grey",main="Degrees of Freedom from Cross Validation")
      abline(v=globaldf,col="red",lty="dashed")
    }
  } else if(is.numeric(dof)){
    df <- dof
    if(verbose) message("\n\r The degree of freedom (dof) is set to the user input:")
    for(i in 1:m) {
      if(verbose) cat(i)
      
      #fitting smoothing splines using smooth.spline(X,Y,df=...)
      spline.fit <- smooth.spline(x=timepoints,y=dat[i,],df=dof)
      dat.denoise[i,] <- predict(spline.fit,newdata = timepoints)$y
    }
  }
  rownames(dat.denoise) <- rownames(dat)
  colnames(dat.denoise) <- colnames(dat)
  
  dat.impute <- dat
  dat.impute[is.na(dat.impute)]  <- dat.denoise[is.na(dat.impute)]
  
  return(list(dat.impute=dat.impute,
              dat.denoise=dat.denoise,
              imputed=is.na(dat),
              df=df))
}

#' @rdname denoise_pca
#' @export
denoise_spline <- preprocess_spline




### Missing Value Imputation
pro_count = matrix(0, nrow = 12, ncol = 7)
counter=0
for(strain in strains){
  for(tre in treat){
    counter = counter+1
    dat <- read.csv(paste0(strain,'_', tre, '_intensity.txt'), 
                    sep = '\t', header = F, row.names = 1)
    dat[dat==0] <- NaN
    dat <- as.matrix(dat)
    colnames(dat) <- timepoints
    dat <- dat[match(pro_list,rownames(dat)),]
    dat <- log(dat)
    for(small_i in 1:7){
    pro_count[counter, small_i] = dim(dat)[1] - sum(is.na(dat[,small_i]))
    }
    result = denoise_spline(dat,timepoints=timepoints,
                            center.dat = FALSE,
                            scale.dat = FALSE,
                            verbose = FALSE,
                            seed = 123)
    assign(paste0(strain,'_',tre,'.dat'), result[[1]])
    write.csv(result[[1]], paste0(strain,'_',tre,'_imputed.csv'), 
              row.names = T, col.names = T)
  }}

a = rep('h',12)
counter = 0
for(strain in strains){
  for(tre in treat){
    counter = counter+1
    a[counter] = paste(strain , tre, sep = '_')
    }}
row.names(pro_count) <- a
colnames(pro_count) = timepoints



