#import map 
library(sf)
#map <-st_read(file.choose())
#map<- map[map$ADM0_NAME %in% c("Rwanda","Uganda","United Republic of Tanzania","Burundi","Democratic Republic of the Congo"), ]
#st_write(map, ".\\raw_data\\map.shp")

map<-st_read(".\\raw_data\\map.shp")
sort(unique(map$ADM0_NAME))
plot(st_geometry(map))
st_coordinates(map)

#coordinate transfer
map<- st_transform(map,3857)
st_coordinates(map)
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

dim(st_coordinates(grid.pred)) 




#import prevalence data

data <- read.csv(".\\raw_data\\data.csv")
dim(data)
head(data)
points(data$long,data$lat,col="red")

dim(data)
data <- data[complete.cases(data),]
dim(data)


table(duplicated(data[,c("long","lat")]))
#coords<- jitter2d(data[,c("long","lat")],max=5000)


#Jittering codes
names(data)
baseline.jits <- data[,c("long","lat")]
base.dups <- which(duplicated(data[,c("long","lat")]))

baseline.jits[base.dups,"long"] <- data[base.dups,"long"] + rnorm(length(base.dups), mean = 0, sd=1)
baseline.jits[base.dups,"lat"] <- data[base.dups,"lat"] + rnorm(length(base.dups), mean = 0, sd=1)

baseline.jits

cbind(data[,3:4],baseline.jits)
#write.csv(baseline.jits,".\\outputs\\baseline.jits.csv")



