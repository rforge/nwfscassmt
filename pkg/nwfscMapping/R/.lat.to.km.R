.lat.to.km <-
function(lat){
# lat in degrees
    lat.rad <- (lat * pi)/180
    return(111.14 - 0.56 * cos(2 * lat.rad))
}
