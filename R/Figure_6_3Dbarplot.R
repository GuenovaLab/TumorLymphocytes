install.packages("plot3D")
install.packages("barplot3d")

library(plot3D)
library(rgl)
library(barplot3d)


# Data description and importation

setwd("M:/DER/RECHERCHE/GUENOVA_LAB/Project_Yun-Tsan_ADCC_paper/raw data and R code/Figure_6/6B/")

data <- read.csv("TPM_IL32_all_sum_result_expression.csv", header=T, row.names =  1, check.names=FALSE)

# If drug data (from column #2) have integers, convert them to numeric
data[1:length(data)] <- lapply(data[1:length(data)], as.numeric)
str(data)



# x, y and z coordinates
x <- rownames(data) #time
y <- colnames(data) #IL4conc
z <- as.matrix(data)

#
barplot3d(rows=7,cols=5,x, y, z,scalexy=500,alpha=0.2,theta=50,phi=20,
          topcolors=c("green3", "skyblue2", "dodgerblue2","violetred2", "red2"),sidecolors = c("gray80"), linecolors = c("gray40"),
          xlabels = rownames(data), ylabels=colnames(data), gap=0.2, gridlines = F, 
          xsub="",ysub="",zsub=""
          )
 
rgl.postscript('3dplot.pdf', fmt = 'pdf')
rgl.snapshot('3dplot.png', fmt = 'png')


barplot3d(rows=3,cols=5,z=inputdata,scalexy=5,alpha=0.4,theta=30,phi=50,
          topcolors=rainbow(15),xlabels = 1:5,ylabels=LETTERS[1:3],
          xsub="Numbers",ysub="Letters",zsub="Count")


  
rgl.postscript('barplot3d.pdf', fmt = 'pdf')



# not working
hist3D (x, y , z,
        bty = "g", phi = 20,  theta = 50,
        xlab = "time", ylab = "IL4 conc. (ng/ml)", zlab = "", main = "IL4 induction",
        col = "#0072B2", border = "grey", shade = 0.05,
        ticktype = "detailed", space = 0.15, d = 2, 
        length = 10, width = 10, 
        cex.lab=1.5, cex.axis=1, cex.main=2) 



hist3D(x,y,z, zlim=c(0,50), theta=40, phi=40, axes=TRUE,label=TRUE, nticks=5, 
       ticktype="detailed", space=0.5, lighting=TRUE, light="diffuse", shade=0.5)


hist3D (x, y , z)




library(plot3D)

# X,Y and Z values
x = c(1,2,3,4,5)
y = c(1,2,3,4,5)

zval = c(20.8, 22.3, 22.7, 11.1, 20.1,
         2.2,  6.7,  14.1, 6.6,  24.7,
         15.7, 15.1, 9.9,  9.3,  14.7,
         8.0,  14.3, 5.1,  6.5,  19.7,
         21.9, 11.2, 11.6, 3.9,  14.8 )

# Convert Z values into a matrix.
z = matrix(zval, nrow=5, ncol=5, byrow=TRUE)


hist3D(x,y,z, zlim=c(0,50), theta=40, phi=40, axes=TRUE,label=TRUE, nticks=5, 
       ticktype="detailed", space=0.5, lighting=TRUE, light="diffuse", shade=0.5)















