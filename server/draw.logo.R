################################################################

library(plotrix) ## to draw arcs

set.seed(19)

################################################################

## simulate interactions, two viewpoints

xlim <- c(1,20000)

viewpoints <- round(runif(2, min=xlim[1], max=xlim[2]))

interactions <- list()
interactions[[1]] <- round(runif(5, min=xlim[1], max=xlim[2]))
interactions[[2]] <- round(runif(7, min=xlim[1], max=xlim[2]))

## add some common points

interactions[[1]] <- c(interactions[[1]], interactions[[2]][1])
interactions[[2]] <- c(interactions[[2]], interactions[[1]][4])

################################################################

## compute radii

radius <- list()
radius[[1]] <- abs(viewpoints[1]-interactions[[1]])/2
radius[[2]] <- abs(viewpoints[2]-interactions[[2]])/2

centers <- list()
centers[[1]] <- (viewpoints[1]+interactions[[1]])/2
centers[[2]] <- (viewpoints[2]+interactions[[2]])/2

################################################################

png(file="img/GOntact_uncropped.png", width=10, height=10, unit="cm", res=800)

################################################################

ymax <- max(unlist(radius))/4
ylim <- c(-ymax*1.1, ymax*1.1)

## canvas

par(mar=c(1,1,1,1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)


## symbolize genes

width <- diff(xlim)/100
height <- diff(ylim)/100

abline(h=-height*1.5, lty=1, col="gray70")

draw.arc(x=centers[[1]], y=rep(0, length(centers[[1]])), radius=radius[[1]], deg1=0, deg2=180, col="steelblue", lwd=1.5)
draw.arc(x=centers[[2]], y=rep(0, length(centers[[2]])), radius=radius[[2]], deg1=0, deg2=180, col="steelblue", lwd=1.5)

################################################################

## schema for gene ontology

xlevel1 <- mean(xlim)
xlevel2 <- seq(from=xlim[1], to=xlim[2], length=5)[-c(1,5)]
xlevel3 <- seq(from=xlim[1], to=xlim[2], length=8)[-c(1,8)]
xlevel4 <- seq(from=xlim[1], to=xlim[2], length=14)[-c(1,14)]

ydiff <- diff(ylim)/15
ylevel4 <- mean(ylim)-ydiff
ylevel3 <- ylevel4-ydiff
ylevel2 <- ylevel3-ydiff
ylevel1 <- ylevel2-ydiff*1.2

## segments that link GO categories among them

for(i in 1:length(xlevel2)){
    segments(xlevel1, ylevel1, xlevel2[i], ylevel2, lty=1, col="gray40")

    for(j in 1:(length(xlevel3)/3)){
        segments(xlevel2[i], ylevel2, xlevel3[2*(i-1)+j], ylevel3,lty=1, col="gray40")
    }
}

for(i in 1:length(xlevel3)){
    for(j in 1:(length(xlevel4)/6)){
        segments(xlevel3[i], ylevel3, xlevel4[2*(i-1)+j], ylevel4,lty=1, col="gray40")
    }
}

## segments that link genes to GO categories

connected <- list()
connected[[1]] <- xlevel4[c(1,3,4,6)]
connected[[2]] <- xlevel4[c(4,6,11,12)]

segments(viewpoints[1], -2*height, connected[[1]], ylevel4, lty=3, col="gray40")
segments(viewpoints[2], -2*height, connected[[2]], ylevel4, lty=3, col="gray40")

################################################################

## rectangles for the genes

rect(viewpoints[1]-width, -2*height, viewpoints[1]+width, -height, col="steelblue", border="black")
rect(viewpoints[2]-width, -2*height, viewpoints[2]+width, -height, col="steelblue", border="black")

################################################################


## circles for GO categories

go.radius <- ymax/3

draw.circle(x=xlevel1, y=ylevel1, radius=go.radius, col="darkgoldenrod3", border="gray40")

for(i in 1:length(xlevel2)){

    draw.circle(x=xlevel2[i], y=ylevel2, radius=go.radius*0.75, col="darkgoldenrod2", border="gray40")
}

for(i in 1:length(xlevel3)){
    draw.circle(x=xlevel3[i], y=ylevel3, radius=go.radius*0.5, col="darkgoldenrod1", border="gray40")
}

for(i in 1:length(xlevel4)){
    draw.circle(x=xlevel4[i], y=ylevel4, radius=go.radius*0.25, col="gold", border="gray40" )
}

################################################################

dev.off()

################################################################
