argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

#plots peaks for each image
plotter.peaks <- function(x,y,w, span,lane,gel) { #w is the window;span smooths data
  peaks <- argmax(x=x, y=y, w=w, span=span)
  
  wd <- getwd()
  name <- paste("./figures/",lane=lane,".png",sep="")
  pdf(paste("./figures/",lane=lane,".pdf",sep=""), width=8, height=7.5)
      #, units="in", res=1200)
  fig <- plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""),las = 1)
  lines(x, peaks$y.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]),
                                    col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
  dev.off()
  
  results <- data.frame(x[peaks$i], peaks$y.hat[peaks$i])
  colnames(results) <- c("band.point","fluorescence")
  #str(results)
  write.csv(file=paste("./output/",gel=gel,lane=lane,".csv",sep=""),results,row.names =FALSE)
  return(results)
}


# function to peak peaks for each lane in gel image - data needs to be incorporated into a single 
peak.finder <- function(input = " ", lanes = "",gel ="",w ,span, output = ""){
  data.in <- read.table(input, sep ="", header = F)
  lane.ids <- read.table(lanes,sep="", header = F)
  #std.lin <- read.table(std, sep = "", header = T)
  
  
  #determines the line distance for each sample
  num.points <- which(data.in[,1] == max(data.in[,1]))[1] 
  print(num.points)
  
  #set image threshold for peak detection
  
  threshold <- mean(data.in[1:num.points*3,2])
  print(threshold)
  
  
  #create data indexing
  start.row <- num.points*3 +1 #rows where sample starts
  
  par(mfrow = c(3,4))
  
  #initate data collection by image
  strain.list <- c()
  band.list <- c()
  peak.list <- c()
  
  for(j in 1:dim(lane.ids)[1]){
      temp <- data.in[start.row:(start.row + num.points-1),]
      head(temp);tail(temp)
      lane <- lane.ids[j,2]
      print(lane)
      
      #x,y coords for peak finding function
      #x.unscaled <- 
      x <- temp[,1]
      y <- temp[,2]
      
      print(head(x));print(tail(x))
      
      plotter.peaks(x=x, y=y, w=w,span=span,gel=gel,lane=lane)
      
      
      peaks <- argmax(x=x, y=y, w=w, span=span)
      plot(x, y, cex=0.75, col="Gray", main=paste("lane = ",lane,", w = ", w, ", span = ", span, sep=""),las = 1)
      lines(x, peaks$y.hat,  lwd=2) #$
      y.min <- min(y)
      sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]),
                                        col="Red", lty=2))
      
      points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
      text(seq(1,length(x[peaks$i]),1),x = x[peaks$i], y = peaks$y.hat[peaks$i]+8)
      results <- data.frame(x[peaks$i], peaks$y.hat[peaks$i])
      #false.pos <- which(results[,2] < threshold)
      true.pts <- results[which(results[,2] > threshold),]
      points(true.pts[,1], true.pts[,2], col="blue", pch=19, cex=1.25)
      
      #append data to list
      strain.list <- append(strain.list,rep(as.character(lane),length(true.pts[,1])))
      band.list <- append(band.list, true.pts[,1])
      peak.list <- append(peak.list, true.pts[,2])
      
      #output.strain <- data.frame(rep(lane,length(results)),results)
      
      start.row <- start.row + num.points 
      print(start.row)
      
      #detect peaks and return value
      #standardize band intensity by standard
      #sum of band intensity/total concentration multiplied by each band 
  }

  results.image <- data.frame(strain.list,band.list,peak.list)
  write.csv(file = output, results.image, row.names = F)
  return(results.image)
}