# --------------------- initialize variables -----------------
stepsize <- 5 # [Input] 
nt <- 100 # [Input] number of steps
nl <- 100 # [Input] number of lubricant molecules
maxBook <- 1 # maximum number of adsorbed molecules per site
probAd <- 0.6 # [Input] adsorption probability when molecule is on site

# define boundary, i.e. "b"
max_b_x <- 100
max_b_y <- 100

# define adsorption site, i.e. "s"
min_s_x <- 20 # [Input] 
max_s_x <- 100 - min_s_x
min_s_y <- min_s_x
max_s_y <- max_s_x

# initialize vectors or matrices
lx <- matrix(0,nt,nl) # x-coordinate of lubricant molecules for every step
ly <- matrix(0,nt,nl) # y-coordinate of lubricant molecules for every step
lcheck <- matrix(0,1,nl) # "0" means not adsorbed, ">0" means adsorbed - record total steps to adsorption
dx <- matrix(round(runif(nt*nl,-stepsize,stepsize)),nt,nl) # molecule displacement in x-axis
dy <- matrix(round(runif(nt*nl,-stepsize,stepsize)),nt,nl)# molecule displacement in y-axis
sbook <- sbookGen(min_s_x, max_s_x) # generate coordinates of adsorption site
# ********************************************************/
# ---------------------initialize position of molecules ----------------------
for (i in 1:nl)
{
  # generate random coordinate
  tempx <- floor(runif(1,1,max_b_x))
  tempy <- floor(runif(1,1,max_b_y))
  
  # check if coordinate lies in adsorption site
  tempCheck <- 0 # set default value
  tempCheck <- inSite(tempx,tempy,min_s_x, min_s_y, max_s_x, max_s_y)
  
  # regenerate coordinate if it lies in adsorption site
  while (tempCheck==1)
  {
    # generate random coordinate
    tempx <- floor(runif(1,1,max_b_x))
    tempy <- floor(runif(1,1,max_b_y))
    
    # check if coordinate lies in adsorbtion site
    # if inSite() returns '0', means NOT on adsorption site
    # if inSite() returns '1', means IS on adsorption site
    tempCheck <- inSite(tempx,tempy,min_s_x, min_s_y, max_s_x, max_s_y)
  }
  
  # assign coordinate
  lx[1,i] <- tempx
  ly[1,i] <- tempy
  
}
# ********************************************************/
# ---------------- moving the molecules ------------------

# note that initially, all molecules are not at the adsorption site
# "ts" is assigned as time index
# "i" is assigned as lubricant molecule index
for (ts in 2:nt) # time loop
{
  for(i in 1:nl) # molecule loop
  {
  
  	# refresh adsorption flag
  	flagL <- 0
  	flagR <- 0
	
	
    if (lcheck[i]>0) # if the molecule is already adsorbed
    {
	
      # here, the molecule is at the adsorption site, i.e. adsorbed
      lx[ts,i] <- lx[ts-1,i] # no displacement
      ly[ts,i] <- ly[ts-1,i] # no displacement
	  
    } else { # if the molecule is not yet adsorbed
	
      lx[ts,i] <- lx[ts-1,i] + dx[ts,i] # x-axis displacement
      ly[ts,i] <- ly[ts-1,i] + dy[ts,i] # y-axis displacement
      
      # adjust if displacement exceed boundary "b"
      # 1. periodic boundary condition on x-axis
      if (lx[ts,i]<0) # if it falls below x=0 line
      {
        lx[ts,i] <- lx[ts,i] + max_b_x # period boundary condition applied
      }
      if (lx[ts,i]>max_b_x) # if it falls above x=max_b_x line
      {
        lx[ts,i] <- lx[ts,i] - max_b_x # period boundary condition applied
      }
      # 2. periodic boundary condition on y-axis
      if (ly[ts,i]<0) # if it falls below y=0 line
      {
        ly[ts,i] <- ly[ts,i] + max_b_y # period boundary condition applied
      }
      if (ly[ts,i]>max_b_y) # if it falls above y=max_b_y line
      {
        ly[ts,i] <- ly[ts,i] - max_b_y # period boundary condition applied
      }
      # *****************************************************
		  # ---- check for adsorption for LEFT tail -----------
		  
		  # check if tail's edge lays on a site, i.e. clay / adsorption site
		  tempx <- lx[ts,i] - 1
		  tempy <- ly[ts,i]
		  tempCheck = inSite(tempx,tempy,min_s_x, min_s_y, max_s_x, max_s_y)
		  
		  # give chance to adsorb using probAd value
		  if (runif(1,0,1) <= probAd && tempCheck==1)
		  {
		    tempCheck <- 1
		  } else{
		    tempCheck <- 0
		  }
		  
		  if (tempCheck==1) # it IS on a site
		  {
  			# see if the site has vacancy
  			sbooknum <- sbooksearch(tempx,tempy,sbook)
			
  			if (sbooknum[1]<maxBook)
  			{
  				sbooknumL <- sbooknum
  				
  				#update flagL
  				flagL <- flagL + 1
  			}
      }
		  # *****************************************************
		  
      # *****************************************************
		  # ---- check for adsorption for RIGHT tail -----------
		  
		  # check if tail's edge lays on a site, i.e. clay / adsorption site
		  tempx <- lx[ts,i] + 1
		  tempy <- ly[ts,i]
		  tempCheck = inSite(tempx,tempy,min_s_x, min_s_y, max_s_x, max_s_y)
		  
		  # give chance to adsorb using probAd value
		  if (runif(1,0,1) <= probAd && tempCheck==1)
		  {
		    tempCheck <- 1
		  } else{
		    tempCheck <- 0
		  }
		  
		  if (tempCheck==1) # it IS on a site
		  {
		  
  			# see if the site has vacancy
  			sbooknum <- sbooksearch(tempx,tempy,sbook)
			
  			if (sbooknum[1]<maxBook)
  			{
  				sbooknumR <- sbooknum
  				
  				#update flagR
  				flagR <- flagR + 1
  			}
      }
		  
		  # *****************************************************
		  # Note:	update sbook (list of clay sites) 
		  #			if the current site has adsorbed the molecule
		  if( flagL>0 && flagR>0 )
		  {
  			# tag the adsorbed molecule & record total number of steps to adsorption
  			lcheck[i] <- ts
  			
  			# update sbook, so that we reduce the vacancy for the clay's site
  			sbook[sbooknumL[2],3] <- sbook[sbooknumL[2],3] + 1
  			sbook[sbooknumR[2],3] <- sbook[sbooknumR[2],3] + 1
		  }
          
    }
    
  }

}

# ********************************************************/
# ---------------- analysis of results ------------------

par(mfrow=c(1,3))

# zoom on adsorption site
plot(lx[,1],ly[,1],type="l", main="Adsorption sites",
     xlim=c(0,max_b_x), ylim=c(0,max_b_y), xlab="X", ylab="Y", col="gray")
points(sbook[,1],sbook[,2],col="red",pch=16,cex=0.1) # clay sites
points(lx[1,1],ly[1,1],col="green",pch=16,cex=1.5) # start point
#points(lx[nt,1],ly[nt,1],col="red",pch=16,cex=1.5) # end point
points(lx[nt,1]-1,ly[nt,1],col="blue",pch=16,cex=0.8) # left tail's edge
points(lx[nt,1]+1,ly[nt,1],col="blue",pch=16,cex=0.8) # right tail's edge
lines(c(lx[nt,1]-1,lx[nt,1]+1), c(ly[nt,1], ly[nt,1]), col="blue")

# plotting paths of next molecule
for (i in 1:nl)
{
  # lines(lx[,i],ly[,i],type="l",col="grey")
  # points(sbook[,1],sbook[,2],col="grey",pch=16,cex=0.1)
  # points(lx[1,i],ly[1,i],col="green",pch=16,cex=1.5) # start point
  # points(lx[nt,i],ly[nt,i],col="red",pch=16,cex=1.5) # end point
  points(lx[nt,i]-1,ly[nt,i],col="blue",pch=16,cex=0.8) # left tail's edge
  points(lx[nt,i]+1,ly[nt,i],col="blue",pch=16,cex=0.8) # right tail's edge
  lines(c(lx[nt,i]-1,lx[nt,i]+1), c(ly[nt,i], ly[nt,i]), col="blue")
}

# plotting paths of 1st molecule (zoom)
plot(lx[,1],ly[,1],type="l", main="Molecules' paths",
     xlim=c(min_s_x-5,max_s_x+5), ylim=c(min_s_x-5,max_s_x+5), xlab="X", ylab="Y", col="white")
points(sbook[,1],sbook[,2],col="red",pch=16,cex=0.1) # clay sites
points(lx[1,1],ly[1,1],col="green",pch=16,cex=1.5) # start point
#points(lx[nt,1],ly[nt,1],col="red",pch=16,cex=1.5) # end point
points(lx[nt,1]-1,ly[nt,1],col="blue",pch=16,cex=0.8) # left tail's edge
points(lx[nt,1]+1,ly[nt,1],col="blue",pch=16,cex=0.8) # right tail's edge
lines(c(lx[nt,1]-1,lx[nt,1]+1), c(ly[nt,1], ly[nt,1]), col="blue")

# plotting paths of next molecule
for (i in 1:nl)
{
  # lines(lx[,i],ly[,i],type="l",col="grey")
  # points(sbook[,1],sbook[,2],col="grey",pch=16,cex=0.1)
  # points(lx[1,i],ly[1,i],col="green",pch=16,cex=1.5) # start point
  # points(lx[nt,i],ly[nt,i],col="red",pch=16,cex=1.5) # end point
  points(lx[nt,i]-1,ly[nt,i],col="blue",pch=16,cex=0.8) # left tail's edge
  points(lx[nt,i]+1,ly[nt,i],col="blue",pch=16,cex=0.8) # right tail's edge
  lines(c(lx[nt,i]-1,lx[nt,i]+1), c(ly[nt,i], ly[nt,i]), col="blue")
}


# plot molecule adsorption status
hist(t(lcheck[lcheck != 0]), xlab="Number of steps", ylab="Number of molecules adsorbed", main="Adsorption Time Distribution")

# compute percentage of clay sites occupied by lubricants
print("-------------------------------------")
print("% clay sites occupied by lubricants")
print(sum(ifelse(sbook[,3]>0,1,0)/(max_s_x-min_s_x)^2))
print("-------------------------------------")

# compute percentage of lubricants adsorbed
print("-------------------------------------")
print("% lubricants adsorbed")
print(sum(ifelse(lcheck>0,1,0))/nl)
print("-------------------------------------")
 

