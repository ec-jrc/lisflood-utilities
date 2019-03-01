
##purpose: estimating return levels from a Gumbel extreme value distribution functions fit to annual maxima 
##Created: Zuzanna, August 2015; 
##Modified: Lorenzo:
##	  Nov 2015  - user input is now read inline as argument to be run as batch script;
##	  Nov 2015  - bug fix for warm_up year and addition of Nyears;
##	  Nov 2015  - header for Ascii file is added on top of the output files. Output files can be directly converted to PCRaster with: asc2map --clone clone.map -S -a ASCIImap.txt ASCIImap.map
##    Jan 2018  - add option to remove values smaller than the long term average discharge
##    Jan 2018  - fixed bug in the threshold estimation in presence of NA and replaced for loops with parallel computing
##    Sep 2018  - fixed bug in the estimation of the gumbel parameters

# USAGE: /usr/bin/Rscript gfit2.r <input.file> <output.file.path> <Return.Periods> <warmup_years> <N of years> <MaxGTavgDis>
# e.g. : /usr/bin/Rscript gfit2.r "/myInputPath/disAMax.nc" "/myOutputPath" 2-5-10-20-50-100-200-500 1 


######################################################
############   Input arguments   #####################

args <- commandArgs(TRUE)
input.file<-as.character(args[1])      #File with annual maxima (from cdo yearmax command)
out.file.path<-as.character(args[2])   #Output folder
rp1<-try(as.character(args[3]))        #Vector of output return periods

# Following options are not mandatory. If not specified, default values are used
warmup_yrs<-try(as.numeric(args[4]))   #default=0         Warm-up years are excluded from the calculations of return levels. Normally this value is zero as Warmup years are removed already in the Gfit2.sh
MaxGTavgDis<-try(as.numeric(args[5]))  #default=1 (Yes)   Option (1 or 0) to use only annual maxima larger than the long term average discharge. The latter must be put in out.file.path and named <whatever>Avg.nc (e.g., disAvg.nc)
Nyears<-try(as.numeric(args[6]))       #default=all       Calculations of return levels is done on this number of years (normally 20-30 years)

ptm <- proc.time() # used for estimating processing time
dir.create(out.file.path, showWarnings = FALSE, recursive = TRUE)

####################################################################################################

rp<-as.numeric(unlist(strsplit(rp1, split='-', fixed=TRUE)))
print(rp)

library(ncdf4) ###reading and saving netcdf files 
library(lmomco) ### L-moments for fitting distributions
library(ismev) ### for fitting Gumbel distribution, when L-moments is not valid, in such case GLM method is used within - gum.fit

####define function for lmom (to optimize calculation of lmoments, lmom.ub calculates up L5 and other stuff that we don't need, we need only L1 and L2)
lmom.ub.LA <- function (x) 
{
    x2<-x[!is.na(x)]
	n <- length(x2)
    if (n == 1) 
        stop("use mean() for data with one value")
    if (n < 5) 
        stop("a minimum of 5 data values are required\n                   because 5 lmoments are to be computed")
    #if (length(unique(x2)) == 1) 
    #    stop("all values are equal--lmoments can not be computed")
    x2 <- sort(x2)
    L1 = 0
    L2 = 0
    for (i in seq(1, n)) {
        CL1 <- i - 1
        CR1 <- n - i
        L1 <- L1 + x2[i]
        L2 <- L2 + x2[i] * (CL1 - CR1)
    }
    C1 <- n
    C2 <- C1 * (n - 1)/2
    L1 <- L1/C1
    L2 <- L2/C2/2
     z <- list(L1 = L1, L2 = L2, source = "lmom.ub")
    return(z)
}


#extract data from source nc file
nc = nc_open(filename = input.file,write=FALSE,suppress_dimvals=FALSE)	#open file with annual maxima from cdo
xDim = nc$var[[1]]$dim[[1]]
yDim = nc$var[[1]]$dim[[2]]
data = ncvar_get(nc)

Units<-nc$var[[1]]$units
nCols<-nc$var[[1]]$size[1]   
nRows<-nc$var[[1]]$size[2]
xllCenter<-min(nc$dim[[1]]$vals)
yllCenter<-min(nc$dim[[2]]$vals)
cellSize<-abs(nc$dim[[1]]$vals[1]-nc$dim[[1]]$vals[2])
MV<-nc$var[[1]]$missval

Header<-list(paste0("ncols      ",nCols),paste0("nrows      ",nRows),paste0("xllcenter  ",xllCenter),paste0("yllcenter  ",yllCenter),paste0("cellsize   ",cellSize),paste0("NODATA_value  ",MV))

nc_close(nc)


# Remove warm up years from the database
if (!exists("warmup_yrs")){warmup_yrs<-0}
if(warmup_yrs>=1){
  data<-data[,,-c(1:warmup_yrs)]
}

# Limit the analysis in the time slice from year 1 to Nyears 
if (!exists("Nyears")){Nyears<-NA}
if(!is.na(Nyears)){
  data<-data[,,c(1:Nyears)]
}

# Remove values smaller than the climatological average discharge
data2use<-data
if(!exists("MaxGTavgDis")){MaxGTavgDis<-1}    #default=1 (Yes)   Option (1 or 0) to use only annual maxima larger than the long term average discharge.

if(MaxGTavgDis==0){

}else{

avgFile<-dir(path=out.file.path, pattern="Avg.nc")[1]
avgDisnc = nc_open(filename = paste0(out.file.path,"/",avgFile),write=FALSE,suppress_dimvals=FALSE)	#open file with annual maxima from cdo
avgData = ncvar_get(avgDisnc)
nc_close(nc)
Ny<-dim(data)[3]
for (k in 1:Ny){
  dx<-data2use[,,k]
  dx[dx<avgData]<-NA
  data2use[,,k]<-dx
}
}


# import grid dimensions from the source nc file
num_grid_x = xDim$len
num_grid_y = yDim$len
num__grid = num_grid_x * num_grid_y

out_return_level = array(MV,dim=c(num_grid_x,num_grid_y,length(rp)))	# define array used for storing calculated return levels


FunGumFit <- function(yrMax){
  if(sum(!is.na(yrMax))>=5)	#process only for cells within a mask (unlike observations that may have some missing data, simulation results have only !NA values)
    {
      try(xmom<-pargum(lmom.ub.LA(yrMax)))	# calculate L-moments,with modified function; only first 2 - this is all we need.
	}
}   

gumParams<-apply(data2use,c(1,2),FunGumFit)
mu<-sapply(sapply(gumParams, "[[", 2),"[[",1)
sig<-sapply(sapply(gumParams, "[[", 2),"[[",2)

mu[sapply(mu, is.null)] <- NA
mu2<-matrix(unlist(mu),nrow=num_grid_x, ncol=num_grid_y)


sig[sapply(sig, is.null)] <- NA
sig2<-matrix(unlist(sig),nrow=num_grid_x, ncol=num_grid_y)


q = 1/rp
for (i in 1:length(rp)){

z = mu2 - sig2 * log(-log(1-q[i])) # Gumbel return level
z[is.na(z)]<-MV
out_return_level[,,i] = z
		
}

mu2[is.na(mu2)]<-MV 
sig2[is.na(sig2)]<-MV

### saving results as .nc file 
# specify netCDF variables
rpVars = vector("list", length(rp))
for (k in 1:length(rp))
	{
	rpVars[[k]] = ncvar_def(paste("rl", rp[k], sep=""), Units, list(xDim, yDim), MV, longname=paste("Return Level", rp[k]), prec="double")
	}
	
# create a netCDF file with all variables
rpNcdf <- nc_create(paste0(out.file.path,"/return_levels.nc"), rpVars)

# write values to each variable on disk.
for (k in 1:length(rp))
	{
	ncvar_put( rpNcdf, rpVars[[k]], out_return_level[, , k], verbose=FALSE)
	}

# add source info metadata to file
ncatt_put( rpNcdf, 0, "source", "return levels from cdo annual maxima")

nc_close(rpNcdf)

############################################################
### option of saving results as .txt files (Comment out if not of interest)
 for (k in 1:length(rp))
	 {
	 if(rp[k]%%1!=0){
	 write.table(sapply(Header, "[[", 1),file=paste0(out.file.path,"/rl_", sprintf("%03.1f",rp[k]), ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
	 write.table(t(out_return_level[, , k]), file=paste0(out.file.path,"/rl_", sprintf("%03.1f",rp[k]), ".txt"), row.names = FALSE, col.names = FALSE, append=TRUE)
	 }else{
	 write.table(sapply(Header, "[[", 1),file=paste0(out.file.path,"/rl_", sprintf("%03.0f",rp[k]), ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
	 write.table(t(out_return_level[, , k]), file=paste0(out.file.path,"/rl_", sprintf("%03.0f",rp[k]), ".txt"), row.names = FALSE, col.names = FALSE, append=TRUE)
	 }
	 }

### save gumbel parameters as ascii files 
write.table(sapply(Header, "[[", 1),file=paste0(out.file.path,"/gumbel_xi.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(t(mu2), file=paste0(out.file.path,"/gumbel_xi.txt"), row.names = FALSE, col.names = FALSE, append=TRUE)

write.table(sapply(Header, "[[", 1),file=paste0(out.file.path,"/gumbel_alpha.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(t(sig2), file=paste0(out.file.path,"/gumbel_alpha.txt"), row.names = FALSE, col.names = FALSE, append=TRUE)
############################################################

print ("time [sec]") # display processing time
print (proc.time() - ptm)
