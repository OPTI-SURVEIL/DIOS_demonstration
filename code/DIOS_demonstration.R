#################################################################################
#
#   This R script is for the demonstration in 
#   "The DIOS framework for optimizing infectious disease surveillance: 
#   numerical methods for simulation and multi-objective optimization of 
#   surveillance network architectures"
#
#                      Qu Cheng (qcheng@berkeley.edu)
#
#################################################################################  

library(ggplot2)
library(raster)
library(mvnfast)
library(gstat)
library(cowplot)
library(geoR)
library(foreach)
library(doParallel)
library(tidyverse)
library(rPref)
library(RColorBrewer)


source("DIOS_demonstration_functions.R")  
set.seed(123)

#==========================================================
#                Demonstration setting
#==========================================================
#================== Parameters ======================
# parameters for points
n.total <- 100
n.initial <- 30

# parameters for linear part
beta0 <- 3 - log(100000)
beta1 <- 1

# parameters for spatial correlation
Y.sigmas <- 1             # sigma_s
Y.range <- c(0.1, 0.3)   # rho
Y.sigmad <- 0.1          # sigma_d

# parameters for the risk factor X
sigma <- 0.2

# number of realizations of the disease model to generate
n.real <- 1000


#=========== generate point locations ================
coord <- data.frame(x = runif(n.total), y = runif(n.total), sampled = 0) # generate n.total random points in the unit grid
sampled.id <- sample(1:n.total, n.initial)   # sample n.initial points from them as the intial design
coord$sampled[sampled.id] <- 1
coord.sampled <- coord[coord$sampled == 1,]   # initial sites
coord.unsampled <- coord[coord$sampled == 0,] # alternative sites


#=========== generate risk factor X in the unit grid ==================
X.coord <- expand.grid(x = seq(0, 1, 0.01), y = seq(0, 1, 0.01)) # gernerate grid

# location of the two point sources. They are set to be far apart from each other.
center1 <- data.frame(x = 0.82, y = 0.55)
center2 <- data.frame(x = 0.24, y = 0.24)

# generate risk factor level
X.coord$Xvalue <- dnorm(sqrt((X.coord$x - center1$x)^2 + (X.coord$y - center1$y)^2), mean = 0, sd = sigma) + dnorm(sqrt((X.coord$x - center2$x)^2 + (X.coord$y - center2$y)^2), mean = 0, sd = sigma)


# plot risk factor
ggplot() + geom_tile(data = X.coord, aes(x = x, y = y, fill = Xvalue)) + coord_equal() + theme_bw() + scale_fill_distiller(palette = "RdYlGn") + geom_point(data = coord[coord$sampled == 0,], aes(x, y), col = "gray60") + geom_point(data = coord[coord$sampled == 1,], aes(x, y), shape = 17, size = 2) + scale_color_manual(values = c("grey","black")) + labs(color = "", fill = "X") + guides(color = FALSE) + xlab("") + ylab("") + geom_point(data = center1, aes(x = x, y = y), col = "blue") + geom_point(data = center2, aes(x = x, y = y), col = "blue") + scale_x_continuous(limits = c(0,1), expand = c(0,0)) + scale_y_continuous(limits = c(0,1), expand = c(0,0))

#============= extract X for sampled and alternative sites ==========
X.raster.raster <- rasterFromXYZ(X.coord)
coord.sampled$Xvalue <- raster::extract(X.raster.raster$Xvalue, coord.sampled[,1:2])
coord.unsampled$Xvalue <- raster::extract(X.raster.raster$Xvalue, coord.unsampled[,1:2])

#========= generate log(Y) when rho = 0.1 or 0.3 ================
sampled.dist <- as.matrix(dist(coord.sampled[,c("x","y")]))

for(i in 1:length(Y.range))
{
  set.seed(123)
  current.cov.sampled <- s.cov(sampled.dist, Y.sigmas, Y.range[i], Y.sigmad)
  coord.sampled[paste("logYValue",Y.range[i], sep = "")] <- beta0 + beta1*coord.sampled$Xvalue + rmvn(n = 1, mu = rep(0, n.initial), sigma = current.cov.sampled, ncores = 7)[1,]
}

#========================= Figure 2 ==============================
Figure2A <- ggplot() + 
  geom_tile(data = X.coord, aes(x = x, y = y, fill = Xvalue)) + 
  scale_fill_distiller(palette = "Blues", direction = 1) +   # background color
  geom_point(data = coord[coord$sampled == 0,], aes(x, y), col = "gray50", shape = 1) +  # alternative sites
  geom_point(data = coord.sampled, aes(x, y),shape = 2, size = 3) +   # initial design
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "white") +   # point source 1
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "white") +   # point source 2
  labs(fill = "Risk factor", x = "", y = "") + 
  ggtitle("(A) Risk factor X") +
  coord_equal() + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

# # kriging interpolation
coord.sampled.sp <- coord.sampled
coordinates(coord.sampled.sp) = ~ x+y

coord.grid <- X.coord
gridded(coord.grid) = ~x+y


krige0.1 <- krige(logYValue0.1 ~ Xvalue, coord.sampled.sp, coord.grid, model = vgm(psill = Y.sigmas^2, model = "Exp", range = Y.range[1], nugget = Y.sigmad^2), nsim = 1, nmax = 60, beta = c(beta0, beta1))
krige0.3 <- krige(logYValue0.3 ~ Xvalue, coord.sampled.sp, coord.grid, model = vgm(psill = Y.sigmas^2, model = "Exp", range = Y.range[2], nugget = Y.sigmad^2), nsim = 1, nmax = 60, beta = c(beta0, beta1))

X.coord$Y0.1 <- krige0.1$sim1
X.coord$Y0.3 <- krige0.3$sim1

Figure2B <- ggplot() + 
  geom_tile(data = X.coord, aes(x = x, y = y, fill = Y0.1)) +   # background color
  geom_point(data = coord[coord$sampled == 0,], aes(x, y), col = "gray50", shape = 1) +  # alternative sites
  geom_point(data = coord.sampled, aes(x, y, fill = logYValue0.1),size = 3, shape = 2) + # initial design
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "black") +  # point source 1
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "black")  +  # point source 2
  geom_contour(data = X.coord, aes(x, y, z = Xvalue), col = "gray") +  # X
  scale_fill_distiller(palette = "RdYlGn", limits = c(min(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))-0.1, max(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))+0.1)) +
  labs(x = "", y = "", fill = "Log prevalence rate", title = expression("(B) Log prevalence rate when "~rho~" = 0.1")) + 
  coord_equal() + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))


Figure2C <- ggplot() + 
  geom_tile(data = X.coord, aes(x = x, y = y, fill = Y0.3)) + 
  geom_point(data = coord[coord$sampled == 0,], aes(x, y), col = "gray50", shape = 1) + 
  geom_point(data = coord.sampled, aes(x, y),shape = 2, size = 3) + 
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "black") + 
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "black") + 
  geom_contour(data = X.coord, aes(x, y, z = Xvalue), col = "gray") + 
  scale_fill_distiller(palette = "RdYlGn", limits = c(min(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))-0.1, max(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))+0.1)) +
  labs(x = "", y = "", fill = "Log prevalence rate", title = expression("(C) Log prevalence rate when "~rho~" = 0.3")) + 
  coord_equal() + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

plot_grid(Figure2A, Figure2B, Figure2C, nrow = 1)
dev.off()





#==========================================================
#               Fitting disease model
#==========================================================
# estimate parameters when rho = 0.1
coord.sampled.geodata <- as.geodata(coord.sampled, coords.col = 1:2, data.col = 5, covar.col = 4, covar.names = "Xvalue")

initial.ml.1 <- likfit(coords = coord.sampled[, c("x","y")], data = coord.sampled$logYValue0.1, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 1), fix.nugget = FALSE, nugget = 1, cov.model = "exponential", lik.method = "ML")

# estimate parameters when rho = 0.3
coord.sampled.geodata <- as.geodata(coord.sampled, coords.col = 1:2, data.col = 5, covar.col = 4, covar.names = "Xvalue")

initial.ml.3 <- likfit(coords = coord.sampled[, c("x","y")], data = coord.sampled$logYValue0.3, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML")






#==========================================================
#    Generate n.real realizations of the disease data
#==========================================================
#============ A. generate n.real sets of the unmeasured 70 sites ==========
# rho = 0.1
coord.sampled.sp <- coord.sampled
coordinates(coord.sampled.sp) = ~ x+y

coord.unsampled.sp <- coord.unsampled
coordinates(coord.unsampled.sp) = ~ x+y

krige0.1.sim <- krige(logYValue0.1 ~ Xvalue, coord.sampled.sp, coord.unsampled.sp, model = vgm(psill = (initial.ml.1$sigmasq)^2, model = "Exp", range = initial.ml.1$phi, nugget = (initial.ml.1$tausq)^2), nsim = n.real, nmax = 60, beta = c(initial.ml.1$beta[1],initial.ml.1$beta[2]))

# rho = 0.3
krige0.3.sim <- krige(logYValue0.3 ~ Xvalue, coord.sampled.sp, coord.unsampled.sp, model = vgm(psill = (initial.ml.3$sigmasq)^2, model = "Exp", range = initial.ml.3$phi, nugget = (initial.ml.3$tausq)^2), nsim = n.real, nmax = 60, beta = c(initial.ml.3$beta[1],initial.ml.3$beta[2]))


#=========== B. generate n.real realizations at the grid for Figure 3=========
# this is not part of DIOS, just for better visualization
# Generate the 41*41 grid 
X.coord <- expand.grid(x = seq(0, 1, 0.025), y = seq(0, 1, 0.025))
X.coord$Xvalue <- dnorm(sqrt((X.coord$x - center1$x)^2 + (X.coord$y - center1$y)^2), mean = 0, sd = sigma) + dnorm(sqrt((X.coord$x - center2$x)^2 + (X.coord$y - center2$y)^2), mean = 0, sd = sigma)

X.coord.sp <- X.coord
coordinates(X.coord.sp) = ~ x + y

# rho = 0.1
krige0.1.sim.grid <- krige(logYValue0.1 ~ Xvalue, coord.sampled.sp, X.coord.sp, model = vgm(psill = (initial.ml.1$sigmasq)^2, model = "Exp", range = initial.ml.1$phi, nugget = (initial.ml.1$tausq)^2), nsim = n.real, nmax = 60, beta = c(initial.ml.1$beta[1],initial.ml.1$beta[2]))

# extract values for the unobserved sites
krige0.1.grid.extract <- list()
for(i in 1:n.real)
{
  current.data <- data.frame(x = X.coord$x, y = X.coord$y, YValue = krige0.1.sim.grid[[i]])
  current.raster <- rasterFromXYZ(current.data)
  krige0.1.grid.extract[[i]] <- raster::extract(current.raster$YValue, coord.unsampled[,1:2])
}

# rho = 0.3
krige0.3.sim.grid <- krige(logYValue0.3 ~ Xvalue, coord.sampled.sp, X.coord.sp, model = vgm(psill = (initial.ml.3$sigmasq)^2, model = "Exp", range = initial.ml.3$phi, nugget = (initial.ml.3$tausq)^2), nsim = n.real, nmax = 60, beta = c(initial.ml.3$beta[1],initial.ml.3$beta[2]))

# extract values for the unobserved sites
krige0.3.grid.extract <- list()
for(i in 1:n.real)
{
  current.data <- data.frame(x = X.coord$x, y = X.coord$y, YValue = krige0.3.sim.grid[[i]])
  current.raster <- rasterFromXYZ(current.data)
  krige0.3.grid.extract[[i]] <- raster::extract(current.raster$YValue, coord.unsampled[,1:2])
}


#==========================================================
#           choose one optimal set
#==========================================================
# Since the design space is small, exhaustive search is done here

#=============== A. Evaluate objective function at alternative sites ========
#======== rho = 0.1 =========
# Took ~15 mins using 7 threads of Intel Core i7-7700 CPU
cl <- makeCluster(detectCores()-1, outfile = "log.txt")
registerDoParallel(cl)
# loop over all realizations
OFV0.1 <- foreach(i = 1:n.real, .combine = "rbind", .packages = c("sp","gstat","geoR")) %dopar%
{
  cat(i, format(Sys.time(),usetz = TRUE), "\n")
  
  # loop over all 70 possible designs
  result <- data.frame(dataID = i, siteID = 1:70, OFV1 = NA, OFV2 = NA)
  for(j in 1:70)
  {
    # combine intial design and the "selected" site
    coord.sampled.current.sp <- rbind(coord.sampled[,c("x","y","Xvalue","logYValue0.1")], data.frame(x = coord.unsampled[j,"x"], y = coord.unsampled[j,"y"], Xvalue = coord.unsampled[j,"Xvalue"], logYValue0.1 = krige0.1.sim[[i]][j]))
    coord.sampled.geodata <- as.geodata(coord.sampled.current.sp, coords.col = 1:2, data.col = 4, covar.col = 3, covar.names = "Xvalue")
    coordinates(coord.sampled.current.sp) <- ~x+y
    
    # remove the "selected" site from the alternative site
    coord.unsampled.current.sp <- coord.unsampled[-j,]
    coordinates(coord.unsampled.current.sp) = ~x+y
    
    # predict the 69 unobserved sites
    krige.current.sim <- krige(logYValue0.1 ~ Xvalue, coord.sampled.current.sp, coord.unsampled.current.sp, model = vgm(psill = (initial.ml.1$sigmasq)^2, model = "Exp", range = initial.ml.1$phi, nugget = (initial.ml.1$tausq)^2), beta = c(initial.ml.1$beta[1],initial.ml.1$beta[2]), debug.level = 0)
    
    # calculate MSE of the spatial prediction - objective 1
    result$OFV1[j] <- mean((krige.current.sim$var1.pred - krige0.1.sim[[i]][-j])^2)
    
    # calculate var(beta1) - Objective function 2
    current.ml <- likfit(coords = coordinates(coord.sampled.current.sp), data = coord.sampled.current.sp$logYValue0.1, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML", messages = FALSE)
    
    result$OFV2[j] <- current.ml$beta.var[2,2]
  }
  result
}
stopCluster(cl)


#======== rho = 0.3 =============
# Took ~15 mins using 7 threads of Intel Core i7-7700 CPU
cl <- makeCluster(detectCores()-1, outfile = "log.txt")
registerDoParallel(cl)
OFV0.3 <- foreach(i = 1:n.real, .combine = "rbind", .packages = c("sp","gstat","geoR")) %dopar%
{
  cat(i, format(Sys.time(),usetz = TRUE), "\n")
  
  result <- data.frame(dataID = i, siteID = 1:70, OFV1 = NA, OFV2 = NA)
  for(j in 1:70)
  {
    # MSE
    coord.sampled.current.sp <- rbind(coord.sampled[,c("x","y","Xvalue","logYValue0.3")], data.frame(x = coord.unsampled[j,"x"], y = coord.unsampled[j,"y"], Xvalue = coord.unsampled[j,"Xvalue"], logYValue0.3 = krige0.3.sim[[i]][j]))
    coord.sampled.geodata <- as.geodata(coord.sampled.current.sp, coords.col = 1:2, data.col = 4, covar.col = 3, covar.names = "Xvalue")
    coordinates(coord.sampled.current.sp) <- ~x+y
    
    coord.unsampled.current.sp <- coord.unsampled[-j,]
    coordinates(coord.unsampled.current.sp) = ~x+y
    
    krige.current.sim <- krige(logYValue0.3 ~ Xvalue, coord.sampled.current.sp, coord.unsampled.current.sp, model = vgm(psill = (initial.ml.3$sigmasq)^2, model = "Exp", range = initial.ml.3$phi, nugget = (initial.ml.3$tausq)^2), beta = c(initial.ml.3$beta[1],initial.ml.3$beta[2]), debug.level = 0)
    
    result$OFV1[j] <- mean((krige.current.sim$var1.pred - krige0.3.sim[[i]][-j])^2)
    # result$OFV1[j] <- mean((krige.current.sim$sim1 - krige0.3.mean[-j])^2)
    
    # var beta1
    current.ml <- likfit(coords = coordinates(coord.sampled.current.sp), data = coord.sampled.current.sp$logYValue0.3, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML", messages = FALSE)
    
    result$OFV2[j] <- current.ml$beta.var[2,2]
  }
  
  result
}
stopCluster(cl)




#========= B. evaluate objective function at 41*41 grids as in Figure 3 ========
#======== rho = 0.1 =============
# Took ~9.5 hrs mins using 7 threads of Intel Core i7-7700 CPU
cl <- makeCluster(detectCores()-1, outfile = "log.txt")
registerDoParallel(cl)
OFV0.1.grid <- foreach(i = 1:n.real, .combine = "rbind", .packages = c("sp","gstat","geoR")) %dopar%   # loop over all realizations
{
  cat(i, format(Sys.time(),usetz = TRUE), "\n")
  
  result <- data.frame(dataID = i, siteID = 1:nrow(X.coord), OFV1 = NA, OFV2 = NA)
  for(j in 1:nrow(X.coord))   # loop over all pixels
  {
    # MSE
    coord.sampled.current.sp <- rbind(coord.sampled[,c("x","y","Xvalue","logYValue0.1")], data.frame(x = X.coord[j,"x"], y = X.coord[j,"y"], Xvalue = X.coord[j,"Xvalue"], logYValue0.1 = krige0.1.sim.grid[[i]][j]))
    coord.sampled.geodata <- as.geodata(coord.sampled.current.sp, coords.col = 1:2, data.col = 4, covar.col = 3, covar.names = "Xvalue")
    coordinates(coord.sampled.current.sp) <- ~x+y
    
    coord.unsampled.current.sp <- coord.unsampled
    coordinates(coord.unsampled.current.sp) = ~x+y
    
    krige.current.sim <- krige(logYValue0.1 ~ Xvalue, coord.sampled.current.sp, coord.unsampled.current.sp, model = vgm(psill = (initial.ml.1$sigmasq)^2 , model = "Exp", range = initial.ml.1$phi, nugget = (initial.ml.1$tausq)^2), beta = c(initial.ml.1$beta[1],initial.ml.1$beta[2]), debug.level = 0)
    
    result$OFV1[j] <- mean((krige.current.sim$var1.pred - krige0.1.grid.extract[[i]])^2)
    
    # var beta1
    current.ml <- likfit(coords = coordinates(coord.sampled.current.sp), data = coord.sampled.current.sp$logYValue0.1, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML", messages = FALSE)
    
    result$OFV2[j] <- current.ml$beta.var[2,2]
  }
  
  result
}
stopCluster(cl)


#======== rho = 0.3 =============
# Took ~9.5 hrs using 7 threads of Intel Core i7-7700 CPU
cl <- makeCluster(detectCores()-1, outfile = "log.txt")
registerDoParallel(cl)
OFV0.3.grid <- foreach(i = 1:n.real, .combine = "rbind", .packages = c("sp","gstat","geoR")) %dopar%   # loop over all realizations
{
  cat(i, format(Sys.time(),usetz = TRUE), "\n")
  
  result <- data.frame(dataID = i, siteID = 1:nrow(X.coord), OFV1 = NA, OFV2 = NA)
  for(j in 1:nrow(X.coord))   # loop over all pixels
  {
    # MSE
    coord.sampled.current.sp <- rbind(coord.sampled[,c("x","y","Xvalue","logYValue0.3")], data.frame(x = X.coord[j,"x"], y = X.coord[j,"y"], Xvalue = X.coord[j,"Xvalue"], logYValue0.3 = krige0.3.sim.grid[[i]][j]))
    coord.sampled.geodata <- as.geodata(coord.sampled.current.sp, coords.col = 1:2, data.col = 4, covar.col = 3, covar.names = "Xvalue")
    coordinates(coord.sampled.current.sp) <- ~x+y
    
    coord.unsampled.current.sp <- coord.unsampled
    coordinates(coord.unsampled.current.sp) = ~x+y
    
    krige.current.sim <- krige(logYValue0.3 ~ Xvalue, coord.sampled.current.sp, coord.unsampled.current.sp, model = vgm(psill = (initial.ml.3$sigmasq)^2, model = "Exp", range = initial.ml.3$phi, nugget = (initial.ml.3$tausq)^2), beta = c(initial.ml.3$beta[1],initial.ml.3$beta[2]), debug.level = 0)
    
    result$OFV1[j] <- mean((krige.current.sim$var1.pred - krige0.3.grid.extract[[i]])^2)
    
    # var beta1
    current.ml <- likfit(coords = coordinates(coord.sampled.current.sp), data = coord.sampled.current.sp$logYValue0.3, trend = trend.spatial(~Xvalue, coord.sampled.geodata), ini.cov.pars = c(1, 0.3), fix.nugget = FALSE, nugget = 0.1, cov.model = "exponential", lik.method = "ML", messages = FALSE)
    
    result$OFV2[j] <- current.ml$beta.var[2,2]
  }
  
  result
}
stopCluster(cl)



#=========== plot results for selecting one additional site for each objective ====
#===========Figure 3 ===========
# coord.sampled <- read.csv("coord.sampled.csv")
# coord.unsampled <- read.csv("coord.unsampled.csv")
# 
# OFV0.1.grid <- read.csv("OFV 0.1 grid.csv")
# OFV0.3.grid <- read.csv("OFV 0.3 grid.csv")

OFV01.agg <- OFV0.1.grid %>% group_by(siteID) %>% summarise(OFV1 = mean(OFV1), OFV2 = mean(OFV2))
OFV03.agg <- OFV0.3.grid %>% group_by(siteID) %>% summarise(OFV1 = mean(OFV1), OFV2 = mean(OFV2))

OFV01.agg$x <- X.coord$x
OFV01.agg$y <- X.coord$y

OFV03.agg$x <- X.coord$x
OFV03.agg$y <- X.coord$y

OFV1.01.raster <- rasterFromXYZ(OFV01.agg[,c("x","y","OFV1","OFV2")])
coord.unsampled$OFV1.01 <- raster::extract(OFV1.01.raster$OFV1, coord.unsampled[,c("x","y")])
coord.unsampled$OFV2.01 <- raster::extract(OFV1.01.raster$OFV2, coord.unsampled[,c("x","y")])

OFV3.01.raster <- rasterFromXYZ(OFV03.agg[,c("x","y","OFV1","OFV2")])
coord.unsampled$OFV1.03 <- raster::extract(OFV3.01.raster$OFV1, coord.unsampled[,c("x","y")])
coord.unsampled$OFV2.03 <- raster::extract(OFV3.01.raster$OFV2, coord.unsampled[,c("x","y")])


Figure3A <- ggplot() + 
  geom_tile(data = OFV01.agg, aes(x = x, y = y, fill =  OFV1)) +  # OFV
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + # point source 1
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + # point source 2
  geom_point(data = coord.sampled, aes(x, y), shape = 2, size = 2) +  # initial design
  geom_point(data = coord.unsampled[-which.min(coord.unsampled$OFV1.01),], aes(x, y), col = "gray40", shape = 1, size = 2) +  # unselected alternative sites
  geom_point(data = coord.unsampled[which.min(coord.unsampled$OFV1.01),], aes(x = x, y = y), col = "cyan", shape = 1, stroke = 2, size = 2) +   # optimal site to add
  scale_fill_distiller(palette = "OrRd", direction = -1) + 
  scale_color_manual(values = c("grey","black")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  labs(fill = "OFV 1", x = "", y = "") + 
  ggtitle("(A) Spatial interpolation, ρ = 0.1") +
  coord_equal() + 
  theme_bw()+ 
  theme(legend.position="bottom", legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.4,"cm"),legend.margin=margin(-20,0,0,0)) + 
  guides(color = FALSE, fill = guide_colourbar(label.theme = element_text(angle = 30, hjust = 1, size = 8), reverse = FALSE))

Figure3B <- ggplot() + 
  geom_tile(data = OFV03.agg, aes(x = x, y = y, fill =  OFV1)) + 
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = coord.sampled, aes(x, y), shape = 2, size = 2) + 
  geom_point(data = coord.unsampled[-which.min(coord.unsampled$OFV1.03),], aes(x, y), col = "gray40", shape = 1, size = 2) + 
  geom_point(data = coord.unsampled[which.min(coord.unsampled$OFV1.03),], aes(x = x, y = y), size = 2, col = "cyan", shape = 1, stroke = 2) + 
  scale_fill_distiller(palette = "OrRd", direction = -1) +
  scale_color_manual(values = c("grey","black")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+ 
  labs(fill = "OFV 1", x = "", y = "") + 
  ggtitle("(B) Spatial interpolation, ρ = 0.3") + 
  coord_equal() + 
  theme_bw()+ 
  theme(legend.position="bottom", legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.4,"cm"),legend.margin=margin(-20,0,0,0)) + 
  guides(color = FALSE, fill = guide_colourbar(label.theme = element_text(angle = 30, hjust = 1, size = 8), reverse = FALSE))

Figure3C <- ggplot() + 
  geom_tile(data = OFV01.agg, aes(x = x, y = y, fill =  OFV2)) + 
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = coord.sampled, aes(x, y), shape = 2, size = 2) + 
  geom_point(data = coord.unsampled[-which.min(coord.unsampled$OFV2.01),], aes(x, y), size = 2, col = "gray40", shape = 1) + 
  geom_point(data = coord.unsampled[which.min(coord.unsampled$OFV2.01),], aes(x = x, y = y), size = 2, col = "cyan", shape = 1, stroke = 2) + 
  scale_fill_distiller(palette = "OrRd", direction = -1) + 
  scale_color_manual(values = c("grey","black")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+ 
  labs(fill = "OFV 2", x = "", y = "") + 
  ggtitle("(C) Effect estimation, ρ = 0.1") + 
  coord_equal() + 
  theme_bw()+ 
  theme(legend.position="bottom", legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.4,"cm"), legend.margin=margin(-20,0,0,0)) + 
  guides(color = FALSE, size = guide_legend(nrow = 2, byrow = TRUE), fill = guide_colourbar(label.theme = element_text(angle = 30, hjust = 1, size = 8), reverse = FALSE))

Figure3D <- ggplot() + 
  geom_tile(data = OFV03.agg, aes(x = x, y = y, fill =  OFV2)) + 
  geom_point(data = center1, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = center2, aes(x = x, y = y), shape = 3, size = 2, col = "white", stroke = 2) + 
  geom_point(data = coord.sampled, aes(x, y), shape = 2, size = 2) + 
  geom_point(data = coord.unsampled[-which.min(coord.unsampled$OFV2.03),], aes(x, y), size = 2, col = "gray40", shape = 1) + 
  geom_point(data = coord.unsampled[which.min(coord.unsampled$OFV2.03),], aes(x = x, y = y), size = 2, col = "cyan", shape = 1, stroke = 2) + 
  scale_fill_distiller(palette = "OrRd", direction = -1) + 
  scale_color_manual(values = c("grey","black")) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+ 
  labs(fill = "OFV 2", x = "", y = "") + 
  ggtitle("(D) Effect estimation, ρ = 0.3") +
  coord_equal() + 
  theme_bw()+ 
  theme(legend.position="bottom", legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.4,"cm"), legend.margin=margin(-20,0,0,0)) + 
  guides(size = guide_legend(nrow = 2, byrow = TRUE), fill = guide_colourbar(label.theme = element_text(angle = 30, hjust = 1, size = 8), reverse = FALSE))

plot_grid(Figure3A, Figure3B, Figure3C, Figure3D)




#====== plot results for selecting one additional site for multiple objectives ====
#======= Figure 4 =======
# OFV0.3 <- read.csv("OFV 0.3.csv")
site.agg <- OFV0.3 %>% group_by(siteID) %>% summarize(OFV1 = mean(OFV1), OFV2 = mean(OFV2))

# OFV0.1 <- read.csv("OFV 0.1.csv")
# site.agg <- OFV0.1 %>% group_by(siteID) %>% summarize(OFV1 = mean(OFV1), OFV2 = mean(OFV2))

site.agg$x <- coord.unsampled$x
site.agg$y <- coord.unsampled$y

sky1 <- psel(site.agg, low(site.agg$OFV1)*low(site.agg$OFV2))
point.col <- brewer.pal(n = nrow(sky1), name = "Set1")

m.obj.p1 <- ggplot(site.agg, aes(x = OFV1, y = OFV2)) + 
  geom_point(shape = 21) + 
  geom_point(data = sky1, size = 3, aes(col = as.factor(siteID))) + 
  theme_bw() + 
  xlab("OFV 1 (spatial interpolation)") + 
  ylab("OFV 2 (effect estimation)") + 
  ggtitle("(A) Pareto set") + 
  scale_color_manual(values = point.col) + 
  guides(col = FALSE)

X.coord <- expand.grid(x = seq(0, 1, 0.01), y = seq(0, 1, 0.01)) # gernerate grid

# generate risk factor level
X.coord$Xvalue <- dnorm(sqrt((X.coord$x - center1$x)^2 + (X.coord$y - center1$y)^2), mean = 0, sd = sigma) + dnorm(sqrt((X.coord$x - center2$x)^2 + (X.coord$y - center2$y)^2), mean = 0, sd = sigma)

coord.grid <- X.coord
gridded(coord.grid) = ~x+y


krige0.1 <- krige(logYValue0.1 ~ Xvalue, coord.sampled.sp, coord.grid, model = vgm(psill = Y.sigmas^2, model = "Exp", range = Y.range[1], nugget = Y.sigmad^2), nsim = 1, nmax = 60, beta = c(beta0, beta1))
krige0.3 <- krige(logYValue0.3 ~ Xvalue, coord.sampled.sp, coord.grid, model = vgm(psill = Y.sigmas^2, model = "Exp", range = Y.range[2], nugget = Y.sigmad^2), nsim = 1, nmax = 60, beta = c(beta0, beta1))

X.coord$Y0.1 <- krige0.1$sim1
X.coord$Y0.3 <- krige0.3$sim1

m.obj.p2 <- ggplot() + 
  geom_tile(data = X.coord, aes(x = x, y = y, fill = Y0.3)) + 
  scale_fill_distiller(palette = "RdYlGn", limits = c(min(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))-0.1, max(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))+0.1)) + 
  coord_equal() + 
  theme_bw() + 
  geom_contour(data = X.coord, aes(x, y, z = Xvalue), col = "gray") + 
  geom_point(data = site.agg, aes(x, y), col = "gray60") + 
  geom_point(data = coord.sampled, aes(x, y), shape = 17, size = 2) + 
  labs(color = "", fill = "log(prev)") + guides(color = FALSE, fill = FALSE) + 
  xlab("") + 
  ylab("") + 
  geom_point(data = sky1, aes(x = x, y = y, col = as.factor(siteID)), size = 4, shape = 17) + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  ggtitle("(B) Locations of the Pareto set") + scale_color_manual(values = point.col) + 
  guides(col = FALSE)

plot_grid(m.obj.p1, m.obj.p2)


#==========================================================
#       Optimization with simulated annealing
#==========================================================
# ================ SA ===============
# estimate initial temperature according to https://www.barnesandnoble.com/p/iterative-computer-algorithms-with-applications-in-engineering-sadiq-m-sait/1101198004/2660895534523?st=PLA&sid=BNB_ADL+Marketplace+Good+Used+Textbooks+-+Desktop+Low&sourceId=PLAGoNA&dpid=tdtve346c&2sid=Google_c&gclid=EAIaIQobChMIzrzJj8P-4wIVCMJkCh1i1ArmEAYYASABEgIObfD_BwE
# ~9.2 hours on two Intel Xeon 12-core Haswell processors (24 cores in total) 
chain1 <- SimAnneal(existing.sites = coord.sampled, alt.sites = coord.unsampled, n.choose = 3, rho = "logYValue0.3", disease.d.set = krige0.3.sim@data, obj.fun = OBJ1, initial.ml = initial.ml.3, neigh.fun = neigh.fun, T0 = 0.075, alpha = 0.9987, T.thres = 1e-6)
chain2 <- SimAnneal(existing.sites = coord.sampled, alt.sites = coord.unsampled, n.choose = 3, rho = "logYValue0.3", disease.d.set = krige0.3.sim@data, obj.fun = OBJ1, initial.ml = initial.ml.3, neigh.fun = neigh.fun, T0 = 0.075, alpha = 0.9987, T.thres = 1e-6)
chain3 <- SimAnneal(existing.sites = coord.sampled, alt.sites = coord.unsampled, n.choose = 3, rho = "logYValue0.3", disease.d.set = krige0.3.sim@data, obj.fun = OBJ1, initial.ml = initial.ml.3, neigh.fun = neigh.fun, T0 = 0.075, alpha = 0.9987, T.thres = 1e-6)
SA1 <- rbind(chain1, chain2, chain3)
SA1 <- data.frame(SA1)



#==== plot results for selecting three additional sites for spatial prediction ====
#======= Figure 5 =======
# rho = 0.3, OFV 1
# SA1 <-  read.csv("SA.OFV1.03.csv")
SA1$BestID <- as.character(SA1$BestID)
SA1$criterion.result <- as.character(SA1$criterion.result)

# plot the chain
result.SA <- NULL
for(i in 1:nrow(SA1))
{
  result.SA <- rbind(result.SA, data.frame(iter = rep(i, each = length(as.numeric(unlist(strsplit(SA1[1,3],","))))), step = 1:length(as.numeric(unlist(strsplit(SA1[1,3],",")))), OFV = as.numeric(unlist(strsplit(SA1[i,3],",")))))
}

Figure5A <- ggplot(result.SA, aes(x = step, y = OFV, col = as.factor(iter), group = iter)) + 
  geom_line() + 
  theme_bw() + 
  xlab("Number of iterations") + 
  ylab("OFV1") + 
  guides(col = FALSE) + 
  ggtitle("(A) Examples of SA runs")

BestID <- as.numeric(unlist(strsplit(SA1$BestID[1], ",")))
Figure5B <- ggplot() + 
  geom_tile(data = X.coord, aes(x = x, y = y, fill = Y0.3)) + 
  geom_contour(data = X.coord, aes(x, y, z = Xvalue), col = "gray")+ 
  geom_point(data = coord.unsampled, aes(x, y), col = "gray60") + 
  geom_point(data = coord.sampled, aes(x, y), shape = 17, size = 2) + 
  geom_point(data = coord.unsampled[BestID,], aes(x = x, y = y), size = 4, shape = 17, col = "blue") + 
  scale_fill_distiller(palette = "RdYlGn", limits = c(min(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))-0.1, max(c(coord.sampled$YValue0.1, coord.sampled$YValue0.3, X.coord$Y0.1, X.coord$Y0.3))+0.1))  + 
  scale_color_manual(values = c("grey","black")) + labs(color = "", fill = "var.both") + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  ggtitle("(B) Locations of 3 optimal sites") +
  guides(color = FALSE, fill = FALSE) + 
  xlab("") + 
  ylab("") +
  coord_equal() + 
  theme_bw()

plot_grid(Figure5A, Figure5B)
