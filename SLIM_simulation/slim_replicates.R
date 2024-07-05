# STEP 1
#Runs SLIM scripts replicates in parallel

writeLines("","log.txt")

selection <- c("NFS")
theta_slope <- c("constant","2x")
replicates <- 1:20

library(foreach)
library(doParallel)
no_cores <- 14 #number of cores to use
cl <- makeCluster(no_cores)
registerDoParallel(cl)

foreach(s = selection) %:%
  foreach(t = theta_slope) %:%
  foreach(rep  = replicates) %dopar% {
    
    sink("log.txt",append = T)
    print(paste("selection",s,"theta",t,"replicate",rep,collapse = "\t"))
    sink()
    
    my_dir <- paste("selection",s,"theta",t,"replicate",rep,sep = "_")
    dir.create(my_dir)
    dir.create(paste0(my_dir,"/original"))
    setwd(paste0(my_dir,"/original"))
    system(paste("slim ","../../Main_recipe_",s,"_",t,"_reducedTheta_10x.slim",sep = ""))
    setwd("../../")
  }
stopCluster(cl)
