carestream.convert = function(gel.input = ""){
  
  # import data
  dat = read.delim(paste0("./data/",gel.input,".txt"),header = FALSE,sep = "\t")
  dat = as.matrix(dat)
  dat = dat[,2:ncol(dat)]
  
  # parameters
  bands = nrow(dat) - 2
  lanes = ncol(dat)/3
  
  # Initate data collection
  output = matrix(nrow = lanes * bands, ncol = 7,NA)
  
  colstart = 1
  rowstart = 1
  
  
  for(i in 1:lanes){
    
    # Make identifiers
    laneID = dat[1,colstart]
    
    output[rowstart:(rowstart+bands-1),2]= strsplit(laneID, ":")[[1]][1]
    
    # Add data
    output[rowstart:(rowstart+bands-1),1] = laneID
    output[rowstart:(rowstart+bands-1),3] = seq(1,bands,1)
    
    if(dat[2,1] == "MW (d)"){
      output[rowstart:(rowstart+bands-1),4] = as.numeric(dat[3:nrow(dat),colstart])/1000
      output[rowstart:(rowstart+bands-1),5:6] = dat[3:nrow(dat),(colstart+1):(colstart+2)]
      
    } else {
      output[rowstart:(rowstart+bands-1),4:6] = dat[3:nrow(dat),colstart:(colstart+2)]
    }
    
    colstart = colstart+3
    rowstart = rowstart + bands  
  }
  
  
  output[,7] = rep(strsplit(gel.input,"_")[[1]][3],nrow(output))
  
  
  
  # trim out bands that do not exist
  output = output[which(output[,4] != ""),]
  
  colnames(output) = c("lane.id","lane","Band","MW","intensity","rel.int","gel.id")
  
  write.csv(output,paste0("./output/carestream-output_",gel.input,".csv"), row.names = FALSE)
  print(head(output))
  return(output)
}