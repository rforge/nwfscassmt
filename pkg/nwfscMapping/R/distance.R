distance <-
function(long1, lat1, long2, lat2){
    sqrt(((long2 - long1) * .long.to.km((lat1 + lat2)/2))^2 + ((lat2 - lat1) * .lat.to.km((lat1 + lat2)/2))^2)
}
