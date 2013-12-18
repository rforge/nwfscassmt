.onLoad = function(libname, pkgname){
# load default data into workspace
require(rjags)
require(R2jags)
require(runjags)
require(superdiag)
require(pscl)
require(statmod)
require(stats)

options(stringsAsFactors=TRUE)
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))

# assign to this environment to keep from overwriting user's workspace
assign("Letters", Letters, envir=.GlobalEnv)

}