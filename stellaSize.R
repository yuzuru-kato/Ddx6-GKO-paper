# stellaSize.R: computing the oocyte volume from Stella-EGFP images
# x: A Masked image of stella-EGFP created by stellaMask.R

kern <- makeBrush(3, shape="diamond")
analyzed_mask <- NULL
psize <- NULL
f1<- NULL
slicenum <- 0
dist2 <- NULL
mask <- NULL

stellaSize <- function(x){
  x <- x > 0
  x <- bwlabel(x)
  dm <- dim(x)[3]-1
  for(var1 in 1:dm){
    slice1 <- x[, , var1]
    print(var1)
    slice1b <- fillHull(slice1)
    cccxy <- computeFeatures.moment(slice1b)[, 1:2] 
    for(var2 in 1:max(slice1)){
      print(c(var1, var2))
      if(max(slice1) == 0){
        break
      }
      if(max(slice1) == 1){
        ccc <- cccxy
      } 
      else {
        ccc <- cccxy[var2, ]　
      }
      getccx1 <- ccc[1] 
      getccy1 <- ccc[2] 
      area <- computeFeatures.shape(slice1)[, 1]
      area2 <- area[var2]
      psize <- append(psize, area2)
      slicenum <- slicenum + 1  
      coord <- which(slice1 == var2, arr.ind=TRUE) 
      xy2 <- 1	
      for(var3 in var1:dm){
        if(xy2 == 0){
          break
        }
        if(var3 == var1){
          getccx2 <- getccx1
          getccy2 <- getccy1
        } 
        else{
          getccx2 <- getccx3
          getccy2 <- getccy3
        }
        nextslice <- x[, , var3 + 1]
        for(var4 in 1:nrow(coord)){
          xy1 <- coord[var4, ] 
          xy2 <- nextslice[xy1[1], xy1[2]] 
          xy2 <- as.numeric(xy2)
          if(xy2 > 0){
            if(var3 + 1 == dim(x)[3]){
              psize <- NULL
              dist2 <- NULL
              slicenum <- 0
            }
            nextsliceb <- fillHull(nextslice) 
            nccxy <- computeFeatures.moment(nextsliceb)[, 1:2] 
            if(max(nextslice) == 1){
              ncc <- nccxy
            } 
            else{
              ncc <- nccxy[xy2, ] 
            }
            getccx3 <- ncc[1] 
            getccy3 <- ncc[2]  
            dist <- sqrt(abs(getccx2 - getccx3)^2 + abs(getccy2 - getccy3)^2) 
            dist2 <- append(dist2, dist)
            nextarea <- computeFeatures.shape(nextslice)[, 1]
            nar <- nextarea[xy2]
            psize <- append(psize, nar)
            slicenum <- slicenum + 1
            coord <- which(nextslice == xy2, arr.ind=TRUE) 
            remove1 <- rmObjects(x[, , var3 + 1], xy2) 
            x[, , var3 + 1] <- remove1
            break
            }
          if((var4 == nrow(coord)) && (xy2 == 0)){
            if(var1 == 1){
              psize <- NULL
              dist2 <- NULL
              slicenum <- 0
            } 
            dist3 <- max(dist2)
            volume_per_pixel <- 0.41400391*0.41400391*1.16
            obsize <- sum(psize)
            oocyte_volume <- obsize*volume_per_pixel
            meanpsize <- obsize/slicenum
            psize2 <- max(psize)
            psize_length <- length(psize)
            psize_last <- psize[psize_length]
            if((3000 < obsize) && (obsize < 14000) && (meanpsize < 1000) && (7 < slicenum) && (slicenum < 20) && (0 < dist3) && (dist3 < 12) && (400 < psize2) && (psize2 < 1400) && (area2 < 600) && (psize_last < 600)){                
              v1 <- c(var1, var2, oocyte_volume, slicenum)
              myname <- c("Slice_number", "Mask_number", "Oocyte_volume", "Number_of_slices")
              names(v1) <- myname
              f1 <- rbind(f1, v1)
              mask <- append(mask, var2)
            } 
            dist2 <- NULL
            psize <- NULL
            slicenum <- 0
          }
        }
      }
    }
    remove2 <- rmObjects(slice1, mask, reenumerate = FALSE)
    remove3 <- slice1 - remove2
    analyzed_mask <- array(append(analyzed_mask, remove3), c(nrow(x), ncol(x), var1))
    mask <- NULL
  }
  writeImage(analyzed_mask, file="lead_analyzed_mask.tif")
  write.csv(f1, file="result.csv")
  return(f1)
  }

