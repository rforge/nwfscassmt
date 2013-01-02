.long.to.km <-
function(lat){
# lat in degrees
    lat.rad <- (lat * pi)/180
    return(111.41 * cos(lat.rad) - 0.1 * cos(3 * lat.rad))
}
