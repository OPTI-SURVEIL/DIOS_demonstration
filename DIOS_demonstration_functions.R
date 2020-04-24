


#' Calcuate the variance-covariance matrix with exponential covariance function
#'
#' @param s.dist   distance matrix
#' @param sigma_s  sigma_s
#' @param scale_s  scale parameter rho
#' @param sigma_d  sigma_d
#' @param nugget   Whether or not to include nugget effect (sigma_d)
#'
#' @return The variance-covariance matrix
s.cov <- function(s.dist, sigma_s, scale_s, sigma_d = 0, nugget = TRUE)
{
  if(nugget)
  {
    sigma_s^2*exp(-s.dist/scale_s) + sigma_d^2*diag(nrow(s.dist))
  } else {
    sigma_s^2*exp(-s.dist/scale_s)
  }
}



#' Proposing new sites in the neighbourhood
#'
#' @param n.all         total number of sites to choose from
#' @param current.ID    IDs of existing sites (row number in alt.sites)
#'
#' @return an array of the new IDs with the same number of sites as the length of current.ID
#' @export
#'
#' @examples
neigh.fun <- function(n.all, current.ID)
{
  all.ID <- n.all    # all IDs
  choose.ID <- all.ID[!all.ID %in% current.ID]   # IDs can be chosen
  
  # remove one from current.ID then add one from the choose.ID
  c(current.ID[-sample(1:length(current.ID), 1)], choose.ID[sample(1:length(choose.ID),1)])
}



















#' Evaluate Objective function 1
#'
#' @param observed    data frame for the observed sites
#' @param candidate   data frame for the unobserved sites
#' @param disease.d   realizations of disease data
#' @param initial.ml  estimated parameters in the disease model
OBJ1 <- function(observed, candidate, disease.d, initial.ml){
  coordinates(observed) <- ~ x+y
  coordinates(candidate) <- ~x+y
  
  krige.current.sim <- krige(YValue ~ Xvalue, observed, candidate, model = vgm(psill = (initial.ml$sigmasq)^2, model = "Exp", range = initial.ml$phi, nugget = (initial.ml$tausq)^2), beta = c(initial.ml$beta[1],initial.ml$beta[2]), debug.level = 0)
  
  mean((krige.current.sim$var1.pred - disease.d)^2)
}

#' Evaluate Objective function 2
#'
#' @param observed    data frame for the observed sites
#' @param candidate   data frame for the unobserved sites
#' @param disease.d   realizations of disease data
#' @param initial.ml  estimated parameters in the disease model
# the 2-4 parameters are useless in the function, just to make it have the same parameters as OBJ1
OBJ2 <- function(observed, candidate, disease.d, initial.ml){
  observed.geodata <- as.geodata(observed, coords.col = 1:2, data.col = 5, covar.col = 4, covar.names = "Xvalue")
  coordinates(observed) <- ~ x+y
  
  current.ml <- likfit(coords = coordinates(observed), data = observed$YValue, trend = trend.spatial(~Xvalue, observed.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML", messages = FALSE)
  
  current.ml$beta.var[2,2]
}




#' simulated annealing
#'
#' @param existing.sites a data frame for initial design
#' @param alt.sites    a data frame for the alternative sites
#' @param n.choose    number of sites to choose
#' @param rho     which colum to use as the dependent variable in exisitng sites
#' @param disease.d.set   realizations of disease data
#' @param obj.fun      objective function
#' @param initial.ml   estimated parameters in the disease model
#' @param neigh.fun    neighborhood function to propose a new design
#' @param T0    initial temperature
#' @param alpha   decreasing rate of the temperature
#' @param T.thres   stopping temperature
SimAnneal <- function(existing.sites,  alt.sites, n.choose, rho = "logYValue0.1", disease.d.set, obj.fun, initial.ml, neigh.fun, T0, alpha, T.thres = 1e-6)
{
  # ncores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
  ncores = detectCores()-7
  
  criterion.result <- array(NA, ceiling(log(T.thres/T0, base = alpha)))
  
  criterion.i = 1
  
  n.all <- 1:nrow(alt.sites)
  
  # initialize
  T <- T0
  # current.count = best.count
  
  # random initial sites
  curID <- sample(n.all, n.choose)
  BestID <- curID
  
  # initial objective function
  selected.sites <- alt.sites[curID,]
  
  existing.sites <- existing.sites[,c("x","y","sampled", "Xvalue", rho)]
  colnames(existing.sites)[5] <- "YValue"
  
  obj.array <- array(NA, length(disease.d.set))
  for(i in 1:length(disease.d.set))
  {
    selected.sites <- rbind(existing.sites, cbind(alt.sites[curID,c("x","y","sampled", "Xvalue")], YValue = disease.d.set[curID, i]))
    
    candidate.sites <- cbind(alt.sites[-curID,], YValue = disease.d.set[-curID, i]) 
    
    obj.array[i] <- obj.fun(selected.sites, candidate.sites, disease.d.set[-curID, i ], initial.ml)
  }
  
  curCost <- mean(obj.array)
  
  BestCost <- Inf
  
  
  cl <- makeCluster(ncores, outfile = "log3.txt")
  registerDoParallel(cl)
  while(T > T.thres){
    # propose a new design
    NewID <- neigh.fun(n.all, curID)
    
    # evaluate the objective function across realizations
    obj.array <- foreach(i = 1:length(disease.d.set), .combine = "c", .packages = c("geoR","gstat","sp"), .export = c("OBJ1","OBJ2")) %dopar%
    {
      selected.sites <- rbind(existing.sites, cbind(alt.sites[NewID,c("x","y","sampled", "Xvalue")], YValue = disease.d.set[NewID,i]))
      
      candidate.sites <- cbind(alt.sites[-NewID,], YValue = disease.d.set[-NewID, i]) 
      
      obj.fun(selected.sites, candidate.sites, disease.d.set[-NewID, i], initial.ml)
    }
    
    
    NewCost <- mean(obj.array)
    det.cost <- NewCost - curCost
    
    
    if(det.cost < 0){  # accept the new design if it performs better
      curID <- NewID
      curCost <- NewCost
      if(NewCost < BestCost)
      {
        BestID <- NewID
        BestCost <- NewCost
      }
    } else    # accept the new design with a probability if it performs worse
    {
      if(runif(1) < exp(-(det.cost/T))){
        curID <- NewID
        curCost <- NewCost
      }
    }

    T <- T*alpha      # update T
    
    print(paste("T = ", round(T, 6), "; BestID:", paste(BestID, collapse = " "), "; BestCost:", round(BestCost, 6), "; CurrentCost:", round(curCost, 6), sep = ""))
    
    criterion.result[criterion.i] = curCost
    criterion.i  = criterion.i  + 1
  }
  
  stopCluster(cl)
  
  c(BestID = paste(BestID,collapse = ","), BestVar = BestCost, criterion.result = paste(criterion.result, collapse = ","))
}