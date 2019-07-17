#data_url <- "https://files.datapress.com/london/dataset/2011-boundary-files/2011_london_boundaries.zip"
#zip_fn <- "london.zip"
#download.file(data_url, zip_fn)
#unzip(zip_fn)

library(maptools)
library(rgdal)
library(data.table)
library(ggplot2)
library(plyr)
library(gganimate)


set.seed(33333334)

crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

lon <- maptools::readShapePoly("LSOA_2011_BFE_City_of_London.shp")
cam <- readShapePoly("LSOA_2011_BFE_Camden.shp")
lon <- rbind(lon, cam)
plot(lon)

lon_coords <- sp::coordinates(lon)  # Coordinates of polygon centers
lon_nbs <- spdep::poly2nb(lon)  # Neighbour relationships based on a triangulation

rownames(lon_coords) <- lon$LSOA11NM
plot(lon)
points(lon_coords)
points(lon_coords[2,1], lon_coords[2,2], col="red")

lon_w_mat <- spdep::nb2mat(lon_nbs)  # Weight matrix between the nodes.
colnames(lon_w_mat) <- rownames(lon_w_mat)
hist(Filter(f = function(x) x > 0, x = as.vector(lon_w_mat)))

num_nbs <- length(lon_nbs)
threshold <- 0.16
plot(lon)
points(lon_coords)

adj_matrix <- lon_w_mat>threshold
dimnames(adj_matrix) <- NULL

edges <- which(adj_matrix, arr.ind = TRUE)
apply(edges, 1, function(r){
  lines(lon_coords[r,1], lon_coords[r,2], col="blue")
})

#transmission rate
beta <- 3

#recovery rate
gamma <- 1

n.nodes <- ncol(adj_matrix)
status.vector <- rep("S", n.nodes)

status.vector[sample(length(status.vector), 1)] <- "I"
names(status.vector) <- 1:n.nodes
colnames(adj_matrix) <- rownames(adj_matrix) <- 1:n.nodes

#matrix and vector to store results at each time step
time.steps <- rep(NA, 2*n.nodes)
status.matrix <- matrix(NA, nrow = 2*n.nodes, ncol = n.nodes)

time.steps[[1]] <- 0
status.matrix[1,] <- status.vector

#run until
stopping.time <- 10

#initialise time
current.time <- 0
time.step <- 1

while(current.time < stopping.time){
  # message for debugging
  # message(paste(c("time step:", time.step, "current time:", current.time), collapse = " "))
  
  #calculate vector recovery rates for each infected person
  recovery.rates <- rep(gamma, sum(status.vector=="I"))
  names(recovery.rates) <- c(1:n.nodes)[status.vector=="I"]
  
  #infection rates for each suceptiable person
  infection.rates <- rowSums(adj_matrix[,status.vector=="I",drop=FALSE])[status.vector=="S"]
  infection.rates <- infection.rates*beta
  names(infection.rates) <- c(1:n.nodes)[status.vector=="S"]
  
  #sample event
  event.rate <- sum(recovery.rates, infection.rates)
  if(event.rate<=0){
    #we've reached the end of the simulation
    break
  }
  
  p <- c(recovery.rates, infection.rates)/event.rate
  event.num <- sample(names(p), 1, prob = p)
  
  #update status
  if(status.vector[event.num]=="S"){
    status.vector[event.num] <- "I"
  } else if (status.vector[event.num]=="I") {
    status.vector[event.num] <- "R"
  } else {
    stop("bad index!")
  }
  
  #sample time to event
  current.time <- current.time + rexp(1, rate = event.rate)
  time.step <- time.step + 1
  
  #update results
  time.steps[[time.step]] <- current.time
  status.matrix[time.step, ] <- status.vector
}

time.steps <- time.steps[1:time.step]
status.matrix <- status.matrix[1:time.step,]

results.df <- data.frame(t(apply(status.matrix, 1, function(r){
  table(factor(r, levels = list("S", "I", "R")))
  })), stringsAsFactors = FALSE)

results.df$time <- time.steps
results.df <- melt(results.df, id.vars="time")

# Plot
ggplot(results.df, aes(x=time, y=value, col=variable)) + geom_line()

results.df <- melt(status.matrix)
colnames(results.df) <- c("time step", "node", "status")
results.df <- cbind(results.df, lon_coords[results.df$node,], rownames(lon_coords))
results.df$node <- rownames(lon_coords)[results.df$node]
colnames(results.df) <- c("time step", "node", "status", "lat", "long", "area")

lon@data$id = rownames(lon@data)
london.points = fortify(lon, region="id")
london.df = join(london.points, lon@data, by="id")

plot.df <- do.call(rbind, lapply(unique(results.df$`time step`), function(t){
  temp.df <- results.df[results.df$`time step`==t,]
  london.df$status <- temp.df$status[match(as.character(london.df$LSOA11NM), temp.df$node)]
  london.df$time <- t
  return(london.df)
}))

gg <- ggplot(plot.df) + 
  aes(long,lat, group=group, fill=status) + 
  geom_polygon() +
  geom_path(color="white") +
  coord_equal() + 
  # Here comes the gganimate code
  transition_manual(time)

anim_save("london.gif",animate(gg))

contents <- base64enc::base64encode("london.gif")
tag <- '<img src="data:image/gif;base64,%s">'
IRdisplay::display_html(sprintf(tag, contents))
