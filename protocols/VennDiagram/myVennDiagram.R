venn.overlap <- 
function(r, a, b, target = 0)
{
#
# calculate the overlap area for circles of radius a and b 
# with centers separated by r
# target is included for the root finding code
#
	pi = acos(-1)
	if(r >= a + b) {
		return( - target)
	}
	if(r < a - b) {
		return(pi * b * b - target)
	}
	if(r < b - a) {
		return(pi * a * a - target)
	}
	s = (a + b + r)/2
	triangle.area = sqrt(s * (s - a) * (s - b) * (s - r))
	h = (2 * triangle.area)/r
	aa = 2 * atan(sqrt(((s - r) * (s - a))/(s * (s - b))))
	ab = 2 * atan(sqrt(((s - r) * (s - b))/(s * (s - a))))
	sector.area = aa * (a * a) + ab * (b * b)
	overlap = sector.area - 2 * triangle.area
	return(overlap - target)
}


myVennDiagram<-function(x,y,digits=2)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr= c(-4,4,-4,4))
  r<- cor(x[x>0 & y>0], y[x>0 & y>0])
  
  txt<-format(r,digits=digits)
  text(0,3,paste("R=",txt),cex=1.5)
  
  tt<-vennCounts(cbind(x,y))
  overlap<-tt[4,3];
  a.uniq<-tt[2,3];
  b.uniq<-tt[3,3];
  sum<-a.uniq+b.uniq+overlap;
  a.frac<-(a.uniq+overlap)/sum;
  b.frac<-(b.uniq+overlap)/sum;
  pi=acos(-1);
  a.r<-sqrt(a.frac/pi);
  b.r<-sqrt(b.frac/pi);
  ab.frac<-overlap/sum;
  rab<-uniroot(venn.overlap,interval=c(max(a.r -b.r,b.r-a.r,0)+0.01, a.r+b.r-0.01), a=a.r,b=b.r,target=ab.frac)$root
  
  xpos <- c(-rab/2*4,rab/2*4);
  ypos <- c(0,0);
  r<-c(a.r,b.r)

  r<-r*4

  theta <- 2 * pi * (1:360)/360
  
 for (i in 1:2)
 {
  lines(xpos[i]+r[i]*cos(theta),ypos[i]+r[i]*sin(theta))
  }

  text(0,0,tt[4,3],cex=1.0)
  text(-rab/2*4-1,0,tt[2,3],cex=1.0)

  text(rab/2*4+1,0,tt[3,3],cex=1.0)
 }





myNewVennDiagram<-function(x,y,digits=2)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr= c(-4,4,-4,4))
  r<- cor(x[x>0 & y>0], y[x>0 & y>0])
  
  txt<-format(r,digits=digits)
  text(0,3,paste("R=",txt),cex=1.5)
  
  tt<-vennCounts(cbind(x,y))
  overlap<-tt[4,3];
  a.uniq<-tt[2,3];
  b.uniq<-tt[3,3];
  sum<-a.uniq+b.uniq+overlap;
  a.frac<-(a.uniq+overlap)/sum;
  b.frac<-(b.uniq+overlap)/sum;
  pi=acos(-1);
  a.r<-sqrt(a.frac/pi);
  b.r<-sqrt(b.frac/pi);
  ab.frac<-overlap/sum;
  rab<-uniroot(venn.overlap,interval=c(max(a.r -b.r,b.r-a.r,0)+0.01, a.r+b.r-0.01), a=a.r,b=b.r,target=ab.frac)$root
  
  xpos <- c(-rab/2*4,rab/2*4);
  ypos <- c(0,0);
  r<-c(a.r,b.r)

  r<-r*4

  theta <- 2 * pi * (1:720)/720
  cc<-c("#FF000030","#00FF0030","#0000FF30","#FFFF0030")
  number<-c(2735,1861,2451,2098)
  ids<-c("LIN9","LIN54","p130","E2F4")
  index<-c(0,0)
  for(i in 1:4)
  {
  	if(a.uniq+overlap == number[i]) {index[1]=i}
  	if(b.uniq+overlap == number[i]) {index[2]=i}
  }

 for (i in 1:2)
 {
  for(j in 1:360)
  {
  if(i==1)
  {	
  lines(c(xpos[i]+r[i]*cos(theta[j]),xpos[i]+r[i]*cos(theta[j])),c(ypos[i]-r[i]*sin(theta[j]),ypos[i]+r[i]*sin(theta[j])),col=cc[index[1]])
  }
  else
  {
  	lines(c(xpos[i]+r[i]*cos(theta[j]),xpos[i]+r[i]*cos(theta[j])),c(ypos[i]-r[i]*sin(theta[j]),ypos[i]+r[i]*sin(theta[j])),col=cc[index[2]])

  	}
  }
 }

  text(0,0,tt[4,3],cex=1.2)
  text(-rab/2*4-2*a.r,0,tt[2,3],cex=1.2)

  text(rab/2*4+2*b.r,0,tt[3,3],cex=1.2)

  
 }
 
 
 
 
 myTextPanel<-function(x,pos,label,...)
 {
 	usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    #text(0.5,0.2,pos,cex=1.5)  
    s=0
     cc<-c("#FF0000A0","#00FF00A0","#0000FFA0","#FFA000A0")
  number<-c(2735,1861,2451,2098)
  ids<-c("LIN9","LIN54","p130","E2F4")
   index=0
   for (i in ids)
   {
   	index=index+1
   	if(i == label) {myid=index;}
   }
        text(0.5, 0.5, label, cex=2.0,col=cc[myid]) 	
 	 text(0.5,0.3,number[myid],cex=1.5)
 	
 }
 
 
 
 
 
 mySTextPanel<-function(x,pos,label,...)
 {
 	usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    #text(0.5,0.2,pos,cex=1.5)  
    s=0
    
     cc<-c("#FF0000A0","#00FF00A0","#0000FFA0","#FFA000A0")
  number<-c(2517,1046,1773,1428)
  ids<-c("S_LIN9","S_LIN54","S_p130","S_E2F4")
   index=0
   for (i in ids)
   {
   	index=index+1
   	if(i == label) {myid=index;}
   }
        text(0.5, 0.5, label, cex=2.0,col=cc[myid]) 	
 	 text(0.5,0.3,number[myid],cex=1.5)
 	
 }
 
 
 
 mySNewVennDiagram<-function(x,y,digits=2)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr= c(-4,4,-4,4))
  r<- cor(x[x>0 & y>0], y[x>0 & y>0])
  
  txt<-format(r,digits=digits)
  text(0,3,paste("R=",txt),cex=1.5)
  
  tt<-vennCounts(cbind(x,y))
  overlap<-tt[4,3];
  a.uniq<-tt[2,3];
  b.uniq<-tt[3,3];
  sum<-a.uniq+b.uniq+overlap;
  a.frac<-(a.uniq+overlap)/sum;
  b.frac<-(b.uniq+overlap)/sum;
  pi=acos(-1);
  a.r<-sqrt(a.frac/pi);
  b.r<-sqrt(b.frac/pi);
  ab.frac<-overlap/sum;
  rab<-uniroot(venn.overlap,interval=c(max(a.r -b.r,b.r-a.r,0)+0.01, a.r+b.r-0.01), a=a.r,b=b.r,target=ab.frac)$root
  
  xpos <- c(-rab/2*4,rab/2*4);
  ypos <- c(0,0);
  r<-c(a.r,b.r)

  r<-r*4

  theta <- 2 * pi * (1:720)/720
  cc<-c("#FF000030","#00FF0030","#0000FF30","#FFFF0030")
  number<-c(2517,1046,1773,1428)
  ids<-c("LIN9","LIN54","p130","E2F4")
  index<-c(0,0)
  for(i in 1:4)
  {
  	if(a.uniq+overlap == number[i]) {index[1]=i}
  	if(b.uniq+overlap == number[i]) {index[2]=i}
  }

 for (i in 1:2)
 {
  for(j in 1:360)
  {
  if(i==1)
  {	
  lines(c(xpos[i]+r[i]*cos(theta[j]),xpos[i]+r[i]*cos(theta[j])),c(ypos[i]-r[i]*sin(theta[j]),ypos[i]+r[i]*sin(theta[j])),col=cc[index[1]])
  }
  else
  {
  	lines(c(xpos[i]+r[i]*cos(theta[j]),xpos[i]+r[i]*cos(theta[j])),c(ypos[i]-r[i]*sin(theta[j]),ypos[i]+r[i]*sin(theta[j])),col=cc[index[2]])

  	}
  }
 }

  text(2*a.r-2*b.r,0,tt[4,3],cex=1.2)
  text(-rab/2*4-2*a.r,0,tt[2,3],cex=1.2)

  text(rab/2*4+2*b.r,0,tt[3,3],cex=1.2)

  
 }
