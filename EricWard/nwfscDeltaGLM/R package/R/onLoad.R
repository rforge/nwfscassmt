.onLoad = function(libname, pkgname){
# load default data into workspace
library(R2jags)
runif(1)
load.module("glm")
options(stringsAsFactors=TRUE)
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))

# assign to this environment to keep from overwriting user's workspace
assign("Letters", Letters, envir=.GlobalEnv)

}