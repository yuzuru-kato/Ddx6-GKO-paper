# stellaFox.R: computing the nuclear to cytoplasmic ratio of FOXO3A signal in oocytes
# Setting: objective lens (100x), Zoom (1x), pixel size (1024 x 1024) 
# Channel1: DAPI
# Channel2: Stella-EGFP
# Channel3: FOXO3A staining

rmaskg <- NULL
rmaskb <- NULL
roocyte <- NULL
roocyte2 <- NULL
oocyte_pixel <- NULL
oocyte_um <- NULL
mint_nuc <- NULL
mint_cyto <- NULL
ratio2 <- NULL

kern <- makeBrush(3, shape="diamond")

stellaFox <- function(x){
  # Creating masks of Stella-EGFP
  imgg <- x[, , 2]
  gaussg <- gblur(imgg, 99, radius=2*ceiling(3*99)+1)
  blurg <- imgg - gaussg
  blurg2 <- blurg > 0
  meanintg <- computeFeatures.basic(blurg2, imgg)[, 1]
  localthg <- thresh(blurg, 50, 50, meanintg/10)
  filledg <- fillHull(localthg)
  maskg <- closing(opening(filledg,kern), kern)
  maskg_label <- bwlabel(maskg)
  areag <- computeFeatures.shape(maskg_label)[,1]
  d5000 <- which(areag < 5000)
  removeg <- rmObjects(maskg_label, d5000)
  removeg2 <- removeg > 0
  removeg3 <- bwlabel(removeg2)
  for(var in 1:max(removeg3)){ 
    coordg <- which(removeg3 == var, arr.ind=TRUE)
    for(var2 in 1:nrow(coordg)){
      xy <- coordg[var2,  ]
      if((xy[1] == 1) || (xy[1] == 1024) || (xy[2] == 1) || (xy[2] == 1024)){
        rmaskg <- append(rmaskg, var)
        break
      } 
    }
  }
  removeg4 <- rmObjects(removeg3, rmaskg)
  removeg5 <- removeg4 > 0
  removeg5_label <- bwlabel(removeg5)
  eccent <- computeFeatures.moment(removeg5_label)[, 4]
  d_eccent <- which(eccent >= 0.8)
  removeg6 <- rmObjects(removeg5_label, d_eccent)
  removeg7 <- removeg6 > 0
  writeImage(removeg7, file="stella_mask.tif")
  
  # Creating DAPI masks in Stella-EGFP masks
  imgb <- x[, , 1]
  imgb_stella <- imgb*removeg7
  gaussb <- gblur(imgb_stella, 35)
  blurb <- imgb_stella - gaussb
  blurb2 <- blurb > 0
  meanintb <- computeFeatures.basic(blurb2, imgb_stella)[,1]
  blurb3 <- blurb > meanintb/6
  maskb <- closing(opening(blurb3, kern), kern)
  maskb2 <- maskb > 0
  maskb3 <- bwlabel(maskb2)
  maskb4 <- computeFeatures.shape(maskb3)[,1]
  d3000 <- which(maskb4 < 3000)
  removeb <- rmObjects(maskb3, d3000)
  removeb2 <- fillHull(removeb)
  removeb2 <- removeb2 > 0
  removeb3 <- bwlabel(removeb2)
  areab <- computeFeatures.shape(removeb3)[,1]
  perimeterb <- computeFeatures.shape(removeb3)[,2]
  matrixb <- rbind(areab, perimeterb)
  for(var3 in 1:ncol(matrixb)){
    circular <- 4*pi*matrixb[1,var3]/(matrixb[2,var3]*matrixb[2,var3])
    if(circular < 1){
        rmaskb <- append(rmaskb, var3)
    }
  }
  removeb4 <- rmObjects(removeb3, rmaskb)
  removeb5 <- removeb4 > 0
  writeImage(removeb5, file="nuclear_mask.tif")
  
  # Creating masks of oocyte cytoplasm
  oocyte <- removeg7 - removeb5
  oocyte2 <- oocyte > 0
  oocyte3 <- bwlabel(oocyte2)
  oocyte4 <- fillHull(oocyte3)
  area1 <- computeFeatures.shape(oocyte3)[,1]
  area2 <- computeFeatures.shape(oocyte4)[,1]
  for(var4 in 1:length(area1)){
    s_area <- area1[var4]
    o_area <- area2[var4]
    if(s_area == o_area){
      roocyte <- append(roocyte, var4)
    }
  }
  oocyte5 <- rmObjects(oocyte4, roocyte)
  oocyte6 <- oocyte5 > 0
  oocyte7 <- oocyte6 - removeb5
  oocyte7 <- oocyte7 > 0
  oocyte8 <- bwlabel(oocyte7)
  oocyte9 <- bwlabel(oocyte6)
  area3 <- computeFeatures.shape(oocyte8)[,1]
  area4 <- computeFeatures.shape(oocyte9)[,1]
  for(var5 in 1:length(area3)){
    s_area2 <- area3[var5]
    o_area2 <- area4[var5]
    if(s_area2 == o_area2){
      roocyte2 <- append(roocyte2, var5)
    }
  }
  oocyte10 <- rmObjects(oocyte8, roocyte2)
  oocyte11 <- oocyte10 > 0
  oocyte12 <- bwlabel(oocyte11)
  writeImage(oocyte12, file="cytoplasmic_mask.tif")
  
  # Creating FOXO3A positive pixels
  imgr <- x[, , 3]
  display(normalize(imgr))
  gaussr <- gblur(imgr, 12, radius=2*ceiling(3*12)+1)
  blurr <- imgr - gaussr
  maskr1 <- blurr > 0
  
  # Analysis of FOXO3A in each oocyte
  oocyte12_area <- computeFeatures.shape(oocyte12)[,1]
  for(var6 in 1:max(oocyte12)){
    num_of_masks <- seq(1:max(oocyte12))
    other_masks <- num_of_masks[-var6]
    mask_for_analysis <- rmObjects(oocyte12, other_masks)
    area_of_oocyte <- oocyte12_area[var6]
    area_of_oocyte_nm <- area_of_oocyte*15625     # 125 nm/pixel, 15625 nm2/pixel
    area_of_oocyte_um <- area_of_oocyte_nm/1000000
    oocyte_pixel <- rbind(oocyte_pixel, area_of_oocyte)
    oocyte_um <- rbind(oocyte_um, area_of_oocyte_um)
    filled_mask <- fillHull(mask_for_analysis)
    nucleus <- filled_mask - mask_for_analysis
    foxo3a_for_analysis_nucleus <- nucleus*maskr1
    foxo3a_for_analysis_cytoplasm <- mask_for_analysis*maskr1
    mint_foxo3a_nucleus <- computeFeatures.basic(foxo3a_for_analysis_nucleus, blurr)[1]
    mint_foxo3a_cytoplasm <- computeFeatures.basic(foxo3a_for_analysis_cytoplasm, blurr)[1]
    ratio1 <- mint_foxo3a_nucleus/mint_foxo3a_cytoplasm
    mint_nuc <- rbind(mint_nuc, mint_foxo3a_nucleus)
    mint_cyto <- rbind(mint_cyto, mint_foxo3a_cytoplasm)
    ratio2 <- rbind(ratio2, ratio1)
  }
  output <- data.frame(oocyte_um, mint_nuc, mint_cyto, ratio2, row.names = NULL)
  myname <- c("oocyte_area_um2", "mean_FOXO3A_nucleus", "mean_FOXO3A_cytoplasm", "FOXO3A_NC_ratio")
  names(output) <- myname
  write.csv(output, file="result.csv")
  print(output)
  return(output)
}
