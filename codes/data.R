#import prevalence data

data <- read.csv(".\\raw_data\\data.csv")
dim(data)
head(data)
points(data$long,data$lat,col="red")

dim(data)
data <- data[complete.cases(data),]
dim(data)

data$obs.prev <- data$y/data$n
head(data)
hist(data$obs.prev,nclass = 20)


#coordinate transfer
coords.ll <- st_as_sf(data,coords = c("long", "lat"), crs = st_crs(map)) #coordinates in long/lat (degrees)
coords.web <-st_transform(coords.ll,3857) #coordinates in web (in meters)
st_coordinates(coords.web)

data$long <- st_coordinates(coords.web)[,1]
data$lat  <- st_coordinates(coords.web)[,2]

table(duplicated(data[,c("long","lat")]))
coords<- jitter2d(data[,c("long","lat")],max=5000)


#import map
library(sf)
#map <-st_read(file.choose())
#map<- map[map$ADM0_NAME %in% c("Rwanda","Uganda","United Republic of Tanzania","Burundi","Democratic Republic of the Congo"), ]
#st_write(map, ".\\raw_data\\map.shp")

map<-st_read(".\\raw_data\\map.shp")
sort(unique(map$ADM0_NAME))
plot(st_geometry(map))
st_coordinates(map)
#map<- map[map$ADM0_NAME=="Rwanda",]

#coordinate transfer
map<- st_transform(map,3857)
plot(st_geometry(map))


#prediction grid
grid.pred <- st_make_grid(map,n=c(50,50),crs=3857,what="centers")
grid.pred

plot(grid.pred)
plot(st_geometry(map), add=T)
plot(grid.pred[map], col = '#ff000088', add = TRUE)

grid.pred <- grid.pred[map]
grid.pred
plot(grid.pred)
plot(st_geometry(map), add=T)

dim(st_coordinates(grid.pred)) #see these are polygons and have more than grid points
new.coords <- st_coordinates(grid.pred)
head(new.coords)

plot(new.coords,asp=1)
plot(st_geometry(map), add=T)


#Jittering codes
baseline.jits <- data[,3:4]
base.dups <- which(duplicated(data[,c("long","lat")]))
baseline.jits[base.dups,"X"] <- data[base.dups,"X"] + rnorm(length(base.dups), mean = 0, sd=1)
baseline.jits[base.dups,"Y"] <- data[base.dups,"Y"] + rnorm(length(base.dups), mean = 0, sd=1)

cbind(data[,3:4],baseline.jits)
