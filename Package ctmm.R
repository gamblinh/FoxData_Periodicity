##############################################################
# Periodograms, Variograms, and AKDE home ranges using ctmm #
##############################################################


# Reload workspace with data and associated models for Randy and Woody
load("2023FoxData")

library(move)
library(ctmm)
library(dplyr)
library(tidyverse)
library(lubridate)

#--------------------------#
# Import and prepare data  #
#--------------------------#


data<-read.csv("data/fox_data.csv")

#Time as a POSIXct object

local.tz<-"America/Winnipeg"
data$loc_datetime<- paste0(data$date, " ", data$time)
data$loc_datetime <- as.POSIXct(data$loc_datetime, tz=local.tz, format = "%Y-%m-%d %H:%M")
data$date <- as.POSIXct(data$date, tz=local.tz, format = "%Y-%m-%d")

data<- data %>%
  mutate(dm = format(as.Date(date), "%Y-%m-%d"))
data$dm <- as.POSIXct(data$dm, tz=local.tz, format = "%Y-%m-%d")


#### Create movestack and convert to list of telemetry objects####

# Create movestack
data.move<- move(x=data$easting, y=data$northing, time=data$loc_datetime, 
                 proj=CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=GRS80"), 
                 data=data, animal=data$name)

class(data.move)

#Create list of telemetry objects
fox.ctmm <- as.telemetry(data.move, drop = FALSE)
plot(fox.ctmm)

# Create individual dataframes

list2env(setNames(fox.ctmm,names(fox.ctmm)), envir = parent.frame()) 


#----------------------#
# Check each fox track #
#----------------------#

#******#
# RF Woody #
#******#

plot(Woody)

#Calculate the variogram
Woody.vg<-variogram(Woody)
plot(Woody.vg)
plot(Woody.vg, fraction = 0.25) #fraction represents the proportion of the variogram object plotted
#The oscillation in the variogram suggests Woody performs a periodic behavior.
plot(Woody.vg,fraction=0.0075) #short lag 24h to check for velocity autocorrelation

#Calculate the periodogram
Woody.pg <- periodogram(Woody)
#Plot periodogram without diagnostics
plot(Woody.pg)
#Plot periodogram with diagnostics and overlay on the data periodogram
plot(Woody.pg, diagnostic = TRUE)

#Use the sliders provided by variogram.fit to specify starting values.
variogram.fit(Woody.vg)
#save the parameters in the manipulator window to GUESS

#Automatically fit the range-resident models via maximum likelihood
Woody.fit.mods<-ctmm.select(Woody,CTMM=GUESS, verbose=TRUE)
summary(Woody.fit.mods) 

# OUF Anisotropic is top model for Woody
Woody.OUF<-Woody.fit.mods[[1]]
summary(Woody.OUF)

plot(Woody.vg, CTMM = Woody.OUF)
plot(Woody.vg, fraction=0.25, CTMM = Woody.OUF)
plot(Woody.vg, fraction=0.0075, CTMM = Woody.OUF)


#******#
# AF Randy #
#******#

# A note on Randy's data: it appears Randy was captured and collared during either an exploratory movement or a dispersal, as the next two days it settles ~25 km from the capture location. It stays in this area for ~ 1 month and then disperses northwest of the study area, where it appears to have established a home range.  

plot(Randy)

#Calculate the variogram
Randy.vg<-variogram(Randy)
plot(Randy.vg)
plot(Randy.vg, fraction = 0.95)
plot(Randy.vg,fraction=0.015) 

#Calculate the periodogram
Randy.pg <- periodogram(Randy)
#Plot periodogram without diagnostics
plot(Randy.pg)
#Plot periodogram with diagnostics and overlay on the data periodogram
plot(Randy.pg, diagnostic = TRUE)

#Use the sliders provided by variogram.fit to specify starting values.
#The default choices are usually acceptable.
variogram.fit(Randy.vg)
#save the parameters in the manipulator window to GUESS

#Automatically fit the range-resident models via maximum likelihood
Randy.fit.mods<-ctmm.select(Randy,CTMM=GUESS, verbose=TRUE)
summary(Randy.fit.mods) 

Randy.OUF<-Randy.fit.mods[[1]]

# Estimated area is crazy high (35854 km2), likely due to huge spread of points
summary(Randy.OUF)

plot(Randy.vg, CTMM = Randy.OUF)
#model appears way off

plot(Randy.vg, fraction=0.25, CTMM = Randy.OUF)
plot(Randy.vg, fraction=0.0075, CTMM = Randy.OUF)

#########################################################################################################
# Function from Appendix
#########################################################################################################
period.test=function(data,model,period=NULL,n=100,deltat=60,...) {
  # initiate output vector
  RES=NULL
  # compute periodogram of the true data
  per=periodogram(data,...)
  # set period of interest at the periodogram peak if not provided by user
  if(is.null(period)) period=1/per$f[which.max(per$LSP)]
  else period=period[1]
  # MC loop
  for (k in 1:n) {
    # simulate dataset using model provided by user and time grid from the true data
    s=simulate(object=model,t=data$t-data$t[1])
    # compute periodogram of the simulated dataset
    s.per=periodogram(s,...)
    # compute whether simulated value is higher than true value
    RES=c(RES, max(s.per$LSP[abs(1/s.per$f - period)<=deltat])<
            max(per$LSP[abs(1/per$f - period)<=deltat]) )
  }
  return(list(period=period,P.value=1-sum(RES)/sum(!is.na(RES))))
}
#########################################################################################################
### returns a list with two items: tested period in seconds and P-value of that period

period.test(Woody, Woody.OUF, n = 100) # period (days) = 14745600/(86400) = 170, p = 0
period.test(Randy, Randy.OUF, n = 100) # period (days) = 14745600/(86400) = 170, p = 0

save.image("2023FoxData")
