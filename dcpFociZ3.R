# dcpFociZ3.R: computing the volume of oocytes and of small DCP1A-foci in oocytes from z-stack images of DCP1A-EGFP.
# Setting: objective lens (100x), Zoom (1x), pixel size (1024 x 1024), Z-stack images (0.42 um) 
# Channel1: DCP1A-EGFP


kern <- makeBrush(3, shape="diamond")
rmask1 <- NULL
rmask2 <- NULL
rmask3 <- NULL
dcp1amask7 <- NULL
pbody_in_oocyte9 <- NULL
osize <- NULL
gsize <- NULL
slicenum <- 0
dist2 <- NULL
mask <- NULL
f1<- NULL
disp <- NULL
o_disp <- NULL
g_disp <- NULL

dcpFociZ3 <- function(x){
  #Creating oocyte masks from DCP1A-EGFP images
  for(var1_1 in 1:dim(x)[3]){
    print(var1_1)
    img <- x[, , var1_1]
    dgauss <- gblur(img, 99, radius=2*ceiling(3*99)+1)
    dblur <- img - dgauss
    dmask1 <- dblur > 0
    dmint <- computeFeatures.basic(dmask1, img)[, 1]
    dmask2 <- dblur > dmint/8
    dmask3 <- opening(closing(dmask2, kern),kern)
    dlabel1 <- bwlabel(dmask3)
    darea1 <- computeFeatures.shape(dlabel1)[,1]
    d1000 <- which(darea1 < 1000)
    dremove1 <- rmObjects(dlabel1, d1000)
    dcp1amask1 <- dremove1 > 0
    dcp1amask2 <- bwlabel(dcp1amask1)
    for(var1_2 in 1:max(dcp1amask2)){
      dcoord <- which(dcp1amask2 == var1_2, arr.ind=TRUE)
      for(var1_3 in 1: nrow(dcoord)){
        xy <- dcoord[var1_3,  ]
        if((xy[1] == 1) || (xy[1] == 1024) || (xy[2] == 1) || (xy[2] == 1024)){
          rmask1 <- append(rmask1, var1_2)
          break
        } 
      }
    }
    dcp1amask3 <- rmObjects(dcp1amask2, rmask1)
    dcp1amask4 <- dcp1amask3 > 0
    dcp1amask4r <- 1 - dcp1amask4
    dcp1amask4r2 <- bwlabel(dcp1amask4r)
    rarea1 <- computeFeatures.shape(dcp1amask4r2)[,1]
    rd500 <- which(rarea1 > 500)
    rremove <- rmObjects(dcp1amask4r2, rd500)
    dcp1amask5 <- dcp1amask4 + rremove
    dcp1amask6 <- bwlabel(dcp1amask5)
    dcp1amask7 <- array(append(dcp1amask7, dcp1amask6), c(1024, 1024, var1_1))
    
    #Creating masks of small DCP1A-foci in oocytes
    gauss <- gblur(img, 50, radius=2*ceiling(3*50)+1)
    blur <- img - gauss
    blur2 <- blur > 0
    meanint <- computeFeatures.basic(blur2, blur)[, 1]
    localth <- thresh(blur, 50, 50, meanint/10)
    gauss2 <- gblur(img, 3, radius=2*ceiling(3*3)+1)
    blur3 <- img - gauss2
    meanint2 <- computeFeatures.basic(localth, blur3)[, 1]
    pbodythresh <- 8*meanint2
    pbodymask <- thresh(blur3, 50, 50, pbodythresh)
    pbodymask2 <- fillHull(pbodymask)
    pbody <- pbodymask2 > 0
    pbody2 <- bwlabel(pbody)
    pbodysize <- computeFeatures.shape(pbody2)[, 1]
    if(is.null(pbodysize) == TRUE){
      rmask2 <- NULL
    }
    if(is.null(pbodysize) == FALSE){
      for(var1_4 in 1:length(pbodysize)){
        foci <- pbodysize[var1_4]
        if(foci > 80){
          rmask2 <- append(rmask2, var1_4)
        }
      }
    }
    pbody3 <- rmObjects(pbody2, rmask2)
    pbody4 <- pbody3 > 0
    largegranule <- pbody - pbody4
    blur5 <- blur2 - largegranule
    meanint4 <- computeFeatures.basic(blur5, img)[, 1]
    localth2 <- thresh(blur, 50, 50, meanint4/10)
    gauss4 <- gblur(img, 3, radius=2*ceiling(3*3)+1)
    blur6 <- img - gauss4
    meanint5 <- computeFeatures.basic(localth2, blur6)[, 1]
    pbodythresh2 <- 7*meanint5
    pbodymask3 <- thresh(blur6, 50, 50, pbodythresh2)
    pbodymask4 <- fillHull(pbodymask3)
    pbody_in_oocyte <- pbodymask4*dcp1amask6
    pbody_in_oocyte2 <- pbody_in_oocyte > 0
    pbody_in_oocyte3 <- bwlabel(pbody_in_oocyte2)
    area2 <- computeFeatures.shape(pbody_in_oocyte3)[, 1]
    d1 <- which(area2 < 3)
    pbody_in_oocyte4 <- rmObjects(pbody_in_oocyte3, d1)
    pbody_in_oocyte5 <- pbody_in_oocyte4 > 0
    pbody_in_oocyte6 <- bwlabel(pbody_in_oocyte5)
    pbodysize2 <- computeFeatures.shape(pbody_in_oocyte6)[, 1]
    if(is.null(pbodysize2) == TRUE){
      rmask3 <- NULL
    }
    if(is.null(pbodysize2) == FALSE){
      for(var1_5 in 1:length(pbodysize2)){
        foci <- pbodysize2[var1_5]
        if(foci > 60){
          rmask3 <- append(rmask3, var1_5)
        }
      }
    }
    pbody_in_oocyte7 <- rmObjects(pbody_in_oocyte6, rmask3)
    pbody_in_oocyte8 <- pbody_in_oocyte7 > 0
    rmask1 <- NULL
    rmask2 <- NULL
    rmask3 <- NULL
    pbody_in_oocyte9 <- array(append(pbody_in_oocyte9, pbody_in_oocyte8), c(1024, 1024, var1_1))
  }
  writeImage(dcp1amask7, file="oocyte.tif")
  writeImage(pbody_in_oocyte9, file="granule.tif")
  
  #Calculating oocyte and granule volumes from stack images
  dm <- dim(dcp1amask7)[3]-1
  for(var2 in 1:dm){
    print(var2)
    oocyte <- dcp1amask7[, , var2]
    granule <- pbody_in_oocyte9[, , var2]
    foocyte <- fillHull(oocyte)
    cccxy <- computeFeatures.moment(foocyte)[, 1:2]
    for(var3 in 1:max(oocyte)){
      if(max(oocyte) == 0){
        break
      }
      if(max(oocyte) == 1){
        ccc <- cccxy
      } 
      else{
        ccc <- cccxy[var3, ]
      }
      getccx1 <- ccc[1]
      getccy1 <- ccc[2]
      marea <- computeFeatures.shape(oocyte)[, 1]
      marea2 <- marea[var3]
      osize <- append(osize, marea2)
      num_of_masks <- seq(1:max(oocyte))
      other_masks <- num_of_masks[-var3]
      mask_for_analysis <- rmObjects(oocyte, other_masks)
      foci_for_analysis <- mask_for_analysis*granule
      foci_for_analysis2 <- foci_for_analysis > 0
      area_of_foci <- computeFeatures.shape(foci_for_analysis2)[1]
      if(is.null(area_of_foci) == TRUE){
        area_of_foci <- 0
      }
      gsize <- append(gsize, area_of_foci)
      slicenum <- slicenum +1 
      o_disp <- array(append(o_disp, mask_for_analysis), c(1024, 1024, slicenum))
      g_disp <- array(append(g_disp, foci_for_analysis), c(1024, 1024, slicenum))
      coord <- which(oocyte == var3, arr.ind=TRUE)
      xy2 <- 1	
      for(var4 in var2:dm){
        if(var4 == var2){
          getccx2 <- getccx1
          getccy2 <- getccy1
        } 
        else{
          getccx2 <- getccx3
          getccy2 <- getccy3
        }
        if(xy2 == 0){
          break
        }
        next_oocyte <- dcp1amask7[, , var4 + 1]
        next_granule <- pbody_in_oocyte9[, , var4 + 1]
        for(var5 in 1:nrow(coord)){
          xy1 <- coord[var5, ]
          xy2 <- next_oocyte[xy1[1], xy1[2]]
          if(xy2 > 0){
            if(var4 + 1 == dim(dcp1amask7)[3]){
              osize <- NULL
              gsize <- NULL
              dist2 <- NULL
              slicenum <- 0
            }
            next_foocyte <- fillHull(next_oocyte) 
            nccxy <- computeFeatures.moment(next_foocyte)[, 1:2]
            if(max(next_oocyte) == 1){
              ncc <- cccxy
            } 
            else {
              ncc <- nccxy[xy2, ]
            }
            getccx3 <- ncc[1]
            getccy3 <- ncc[2]
            dist <- sqrt(abs(getccx2 - getccx3)^2 + abs(getccy2 - getccy3)^2)
            dist2 <- append(dist2, dist)
            nextarea <- computeFeatures.shape(next_oocyte)[, 1]
            nar <- nextarea[xy2]
            osize <- append(osize, nar)
            next_num_of_masks <- seq(1:max(next_oocyte))
            next_other_masks <- next_num_of_masks[-xy2]
            next_mask_for_analysis <- rmObjects(next_oocyte, next_other_masks)
            next_foci_for_analysis <- next_mask_for_analysis*next_granule
            next_foci_for_analysis2 <- next_foci_for_analysis > 0
            next_area_of_foci <- computeFeatures.shape(next_foci_for_analysis2)[1]
            if(is.null(next_area_of_foci) == TRUE){
              next_area_of_foci <- 0
            }
            gsize <- append(gsize, next_area_of_foci)
            slicenum <- slicenum + 1
            o_disp <- array(append(o_disp, next_mask_for_analysis), c(1024, 1024, slicenum))
            g_disp <- array(append(g_disp, next_foci_for_analysis), c(1024, 1024, slicenum))
            coord <- which(next_oocyte == xy2, arr.ind=TRUE)
            remove1 <- rmObjects(dcp1amask7[, , var4 + 1], xy2)
            dcp1amask7[, , var4 + 1] <- remove1
            break
          }
          if((var5 == nrow(coord)) && (xy2 == 0)){
            if(var2 == 1){
              osize <- NULL
              gsize <- NULL
              dist2 <- NULL
              slicenum <- 0
            } 
            dist3 <- max(dist2)
            volume_per_pixel <- 0.12400391*0.12400391*0.42
            sum_osize <- sum(osize)
            sum_gsize <- sum(gsize)
            area_of_oocyte_nm <- sum_osize*15625
            area_of_foci_nm <- sum_gsize*15625
            area_of_oocyte_um <- area_of_oocyte_nm/1000000
            area_of_foci_um <- area_of_foci_nm/1000000
            oocyte_volume <- volume_per_pixel*sum_osize
            foci_volume <- volume_per_pixel*sum_gsize
            unit_focivolume <- foci_volume/oocyte_volume
            mean_osize <- sum_osize/slicenum
            mean_gsize <- sum_gsize/slicenum
            max_osize <- max(osize)
            osize_length <- length(osize)
            osize_last <- osize[osize_length]
            if((1500 < area_of_oocyte_um) && (25 < slicenum) && (50 > dist3)){
              v1 <- data.frame(var2, var3, oocyte_volume, foci_volume, unit_focivolume, row.names = NULL)
              myname <- c("Slice_number", "Mask_number", "Oocyte_volume", "Foci_volume", "Foci_volume_per_Oocyte_volume")
              names(v1) <- myname
              f1 <- rbind(f1, v1)
              mask <- append(mask, var3)
            }
            dist2 <- NULL
            osize <- NULL
            gsize <- NULL
            slicenum <- 0
            o_disp <- NULL
            g_disp <- NULL
          }
        }
      }
    }
    remove2 <- rmObjects(oocyte, mask, reenumerate = FALSE)
    remove3 <- oocyte - remove2
    disp <- array(append(disp, remove3), c(nrow(x), ncol(x), var2))
    mask <- NULL
  }
  writeImage(disp, file="lead_analyzed_oocyte.tif")
  write.csv(f1, file="results.csv")
  return(f1)
}