


#'Functions to get the part of string after string split
#'names: a vector contains strings 
#'i the part of the strings you want to extract
StrExtract=function(strings,split="_",i){
  y=strsplit(strings,split)
  len=unlist(lapply(y,length))
  if(length(which(len<i))>0) {
    warning("Strings cannot be extracted in some items")
  }else{}
  
  y=lapply(y,function(x) {
    if(length(x)<i){
      ori=paste(x,collapse=split)
      return(ori)
    }else{return(x[i])}
  })
  y=unlist(y)
  return(y)
}

#Margin removal function

removeMargins<- function(f,chans,sens=1, debris=FALSE,neg=500, verbose = T,return.gate=F)
{
  neg <-rep(neg,length(chans))
  data <- exprs(f)
  margins <- c()
  marg.gates<-c()
  if(!debris)
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.max <-max(data[,chan],na.rm = T)
      margins <- which ( data[, chan] >= stain.max*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
      marg.gates <- append(marg.gates,stain.max*sens-1 )
    }
    if(return.gate)
      return(marg.gates)
  }else
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.min <-min(data[,chan],na.rm=T)
      margins <- which ( data[, chan] <= stain.min*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
  }
  for(i in 1:length(chans))
  {
    if (neg[i]<500)
    {
      ch<-ifelse (is.character(chans[i]),yes = which(colnames(f)==chans[i]),no = chans[i])
      negs <- which ( data[, ch] < neg[i])
      margins <- negs
      if (length(margins)!=0){
        data <- data[ -margins, ]
      }
      if(verbose == T){print(paste(length(margins), "negative events in",colnames(f)[ch], "will be removed.",sep =" "))}
    }
  } 
  exprs(f) <- data
  
  return(f)
  
  
    
}




#Creating estimateLogicle transformation filter for flowFrame, flowSet or GatingSet
#updated to fix the problem of trans.chans 

transform.flow_pre_cal <- function(obj,trans.chans=NULL, p_preC)
{
  #obj: flowFrame, flowSet, or a GatignSet
  #remove.outlier: Default to T, removing far away cells.
  #sd.coeff: Default to 4.5 for removing cells that are below or abour 4.5*sd
  #trans.chans: Channels to be transformed. if it's null, then it tries to find LOG channel in the frame
  if (class(obj)=="flowFrame")
  {
    f <- obj
    f.t <- f
  }else{
    temp <-  obj
    if (class(obj)=="GatingSet")
      temp <- getData(obj, tail(getNodes(obj),1))
    f.t <- temp[[1]]
    f<-getGlobalFrame(temp)
    
  }
  if (is.null(trans.chans))
  {
    # This section is changed to fix the column selection problem introduced in keep.col
    #-----------------------------------------------
    #log.channels <- paste("P",1:ncol(f.t),"DISPLAY",sep="")
    log.channels <- names(f.t@description)[grepl("DISPLAY", names(f.t@description) )]
    
    temp_order <- gsub("([0-9]+).*$", "\\1", log.channels)
    temp_order <- as.numeric(gsub("P", "" ,temp_order))
    log.channels <- log.channels[order(temp_order)]
    
    
    trans.chans <- which(f.t@description[log.channels]=="LOG")
    trans.chans <- gsub("DISPLAY", "N", names(trans.chans))   
    trans.chans <- paste("$",trans.chans,sep="")
    trans.chans <- match(unlist(f@description[trans.chans]), colnames(f.t))
    # make sure the index don't mess up with 
    #----
    if (length(trans.chans)==0)
    {
      warnings("Couldn't find Log channels, all channels except FSC, SSC, and time will be transformed.")
      trans.chans <- 1:ncol(f.t)  
      trans.chans <- trans.chans[-c(grep(colnames(f.t),pattern = "FSC*"),grep(colnames(f.t),pattern = "SSC*"),grep(colnames(f.t),pattern = "Time*"))] 
    }
  }
  
  
  stopifnot(identical(rownames(p_preC$parameters), colnames(f)[trans.chans]))
  
  trans <- p_preC$trans
  #trans <- flowJo_biexp_trans()
  #trans <- transformerList(colnames(f)[trans.chans], trans)
  
  return(transform(obj,trans))
}






#sampling from all flowFrames


getGlobalFrame <- function(fs, sample.length=NA, all.cells=F){
  set.seed(123)
  if (is(fs, 'flowFrame')){
    return (frame)
  }
  n <- length(fs)
  sample.n <- ifelse(is.na(sample.length),n,sample.length)
  global.frame <- fsApply(fs[sample(n, sample.n)],
                          function(frame){
                            m <- nrow(frame)
                            sample.size <- ifelse (all.cells,yes = m,no =min(m,ceiling(m/sample.n )))
                            
                            exprs(frame )<- exprs(frame)[sample(m, sample.size),,drop =F]
                            return (frame)
                          })
  global.frame <- as(global.frame, 'flowFrame')
  return (global.frame)
}





##Finds markers in the FCS file
Find.markers <- function(frame,marker.list)
{
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  if(class(frame)=="cytoframe") {
    frame <- cytoframe_to_flowFrame(frame)
  } # flowWorkspace 4 uses cytoframe instead of flwoframe, so need to convert
  
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}