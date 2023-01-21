# stellaMask.R: Creating masks of oocyte from Stella-EGFP images
# Setting: objective lens: 20x, Zoom setting: 1.5x, Pixel size: 1024x1024
# Channel1: Stella-EGFP

kern <- makeBrush(3, shape="diamond")
arr <- NULL

stellaMask <- function(x){
	for(var in 1:dim(x)[3]){
	  print(var)
		rmask <- NULL
		imgv <- x[, , var]
		gauss <- gblur(imgv, 7, radius=2*ceiling(3*7)+1)
		blur <- imgv - gauss
		mask1 <- blur > 0
		meanint <- computeFeatures.basic(mask1, imgv)[, 1]
		mask2 <- blur > meanint/10
		mask3 <- closing(opening(mask2, kern),kern)
		label1 <- bwlabel(mask3)
		area1 <- computeFeatures.shape(label1)[,1]
		d100 <- which(area1 < 100)
		remove1 <- rmObjects(label1, d100)
		remove1 <- remove1 > 0
		label2 <- bwlabel(remove1)
		for(var2 in 1: max(label2)){ 
			coord <- which(label2 == var2, arr.ind=TRUE)
			for(var3 in 1: nrow(coord)){
				xy <- coord[var3,  ]
				if((xy[1] == 1) || (xy[1] == 1024) || (xy[2] == 1) || (xy[2] == 1024)){
					rmask <- append(rmask, var2)
					break
				} 
			}
		}
		label3 <- rmObjects(label2, rmask)
		label3 <- label3 > 0
		label3 <- bwlabel(label3)
		arr <- array(append(arr, label3), c(nrow(x), ncol(x), var))
	}
  writeImage(arr, file="stella_mask.tif")
	return(arr)
	} 