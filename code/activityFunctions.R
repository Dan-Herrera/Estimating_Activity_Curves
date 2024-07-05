#The functions in this script can be used to analyze and predict activity densities of animals.
#Please contact Dan Herrera (herrerawildlife@gmail.com) with any questions or issues
#########################################

format_idealPDF <- function(time.start, time.stop, observedPDF, dist, param1 = NA, param2 = NA, param3 = NA){ #begin function (used to be: time.start, time.stop, shape, by=10)
  
  # Error messages for start and stop times which are not character strings with 4 numbers split by a colon
  if(class(time.start) %in% "character" == FALSE){stop("time.start must be a character string in format: 'HH:MM'.")}
  if(class(time.stop) %in% "character" == FALSE){stop("time.stop must be a character string in format: 'HH:MM'.")}
  if(nchar(time.start) != 5){stop("time.start must be a character string in format: 'HH:MM'.")}
  if(nchar(time.stop) != 5){stop("time.stop must be a character string in format: 'HH:MM'.")}
  
  #Error messages for shape 
  shapes <- c("uniform","normal","binomial","beta","bimodal","kumaraswamy","triangular")
  if(dist %in% shapes == FALSE){stop("the following shapes are permissable (but custom shapes can be used as long as the vector length matches that of the disturbance vector and the vector sums to 1): beta, binomial, bimodal, kumaraswamy, normal, triangular, uniform.")}
  
  #Find time of day, as an integer, for start times
  time.start <- as.character(time.start) #Convert time.start to a character if it is not already
  time.start <- sub(x = time.start, #Remove all text before and including the space between date and time, leaving only time as a character
                    pattern = ".* ", #Removes anything before a space
                    replacement = "")
  start.hours <- sub(x = time.start, #Remove all text after the colon, leaving hours as a character
                     pattern = ":.*", #Removes anything after a colon
                     replacement = "")
  start.minutes <- sub(x = time.start, #Remove all text before the colon, leaving minutes as a character
                       pattern = ".*:", #Removes anything before a colon
                       replacement = "")
  start.hours <- as.numeric(start.hours)*60 #Multiply hours by 60 to convert to minutes
  start.time <- start.hours + as.numeric(start.minutes) #Add hours (in minutes) to minutes to determine number of minutes into the day
  
  
  #Find time of day, as an integer, for start times
  time.stop <- as.character(time.stop) #Convert time.start to a character if it is not already
  time.stop <- sub(x = time.stop, #Remove all text before and including the space between date and time, leaving only time as a character
                   pattern = ".* ",
                   replacement = "")
  stop.hours <- sub(x = time.stop, #Remove all text after the colon, leaving hours as a character
                    pattern = ":.*", #Removes anything after a colon
                    replacement = "")
  stop.minutes <- sub(x = time.stop, #Remove all text before the colon, leaving minutes as a character
                      pattern = ".*:", #Removes anything before a colon
                      replacement = "")
  
  #Adjust for stretches of time which include midnight
  if(as.numeric(substr(time.start,0,2)) > as.numeric(substr(time.stop,0,2))) { #start instructions specific to time stretches with midnight
    stop.hours <- as.numeric(stop.hours)+24
  } #end instructions specific to time stretches with midnight
  
  #Adjust for 24-hour stretches of time
  if(as.numeric(substr(time.start,0,2)) == as.numeric(substr(time.stop,0,2))) { #start instructions specific to time stretches with midnight
    stop.hours <- as.numeric(stop.hours)+24
  } #end instructions specific to time stretches with midnight
  
  stop.hours <- as.numeric(stop.hours)*60 #Multiply hours by 60 to convert to minutes
  stop.time <- stop.hours + as.numeric(stop.minutes) #Add hours (in minutes) to minutes to determine number of minutes into the day
  
  #Create data frame of values between start and stop times
  df <- data.frame(t = seq(from = start.time, 
                           to = stop.time, 
                           length.out = length(observedPDF)),
                   observed = observedPDF)
  
  if(dist %in% "uniform"){#start instructions for uniform distribution
    df$PDF <- 1/nrow(df)
  }
  
  if(dist %in% "normal"){#start instructions for normal distribution
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), often called the mean (as a proportion)")}
    if(is.na(param2)){stop("this distribution requires a shape parameter (param2), often called the standard deviation (as a proportion)")}
    if(param1 < 0){warning("param1 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param1 > 1){warning("param1 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param2 < 0){warning("param2 (sd) should be a proportion (0-1) of the active period")}
    if(param2 > 1){warning("param2 (sd) should be a proportion (0-1) of the active period")}
    
    
    df$PDF <- dnorm(x = as.vector(df$t),
                    mean = param1*(max(df$t)-min(df$t))+min(df$t),
                    sd = param2*nrow(df))
    
    df$PDF <- df$PDF/sum(df$PDF) #ensure that PDF sums to 1
  }
  
  if(dist %in% "binomial"){#start instructions for binomial distribution
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), often called the probability of success")}
    if(param1 < 0){stop("the shape parameter is a probability and must be between 0 and 1")}
    if(param1 > 1){stop("the shape parameter is a probability and must be between 0 and 1")}
    
    
    df$PDF <- dbinom(x = as.numeric(rownames(df)),
                     size = nrow(df),
                     prob = param1)
    
    df$PDF <- df$PDF/sum(df$PDF) #ensure that PDF sums to 1
  }
  
  if(dist %in% "beta"){#start instructions for beta distribution
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), often called alpha")}
    if(is.na(param2)){stop("this distribution requires a second shape parameter (param2), often called beta")}
    
    if(param1 < 0){stop("both shape parameters must be positive numbers")}
    if(param2 < 0 ){stop("both shape parameters must be positive numbers")}
    
    tprop <- as.vector(as.numeric(row.names(df))/nrow(df)) #create a column where t is shown as a proportion, not an integer
    
    df$PDF <- dbeta(tprop,
                    shape1 = param1,
                    shape2 = param2)
    
    df$PDF <- df$PDF/sum(df$PDF) #ensure that PDF sums to 1
    
  }
  
  if(dist %in% "bimodal"){#start instructions for bimodal distribution
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), which is the mean of the first peak (as a proportion)")}
    if(is.na(param2)){stop("this distribution requires a second shape parameter (param2), which is the mean of the second peak (as a proportion)")}
    if(is.na(param3)){stop("this distribution requires a third shape parameter (param3), which is a single standard deviation for both peaks (as a proportion")}
    
    if(param1 < 0){warning("param1 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param1 > 1){warning("param1 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param2 < 0){warning("param2 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param2 > 1){warning("param2 (mean) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param3 < 0){warning("param3 (sd) should be a proportion (0-1) of the active period")}
    if(param3 > 1){warning("param3 (sd) should be a proportion (0-1) of the active period")}
    if(param1 == param2){warning("param1 and param2 are currently identical, which will produce only one peak")}
    
    peak.1 <- dnorm(x = as.vector(df$t),
                    mean = param1*(max(df$t)-min(df$t))+min(df$t),
                    sd = param3*nrow(df))
    
    peak.2 <- dnorm(x = as.vector(df$t),
                    mean = param2*(max(df$t)-min(df$t))+min(df$t),
                    sd = param3*nrow(df))
    
    both.peaks <- peak.1 + peak.2
    
    df$PDF <- both.peaks/sum(both.peaks) #ensure that PDF sums to 1
    
  }
  
  if(dist %in% "kumaraswamy"){#start instructions for kumaraswamy distribution
    if(system.file(package = "extraDistr") %in% ""){install.packages("extraDistr")}
    library(extraDistr)
    
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), often called alpha")}
    if(is.na(param2)){stop("this distribution requires a second shape parameter (param2), often called beta")}
    
    tprop <- as.vector(as.numeric(row.names(df))/nrow(df)) #create a column where t is shown as a proportion, not an integer
    
    kumar.vect <- dkumar(x = tprop,
                         a= param1,
                         b= param2) 
    
    df$PDF <- kumar.vect/sum(kumar.vect)
    
  }
  
  if(dist %in% "triangular"){#start instructions for kumaraswamy distribution
    if(system.file(package = "extraDistr") %in% ""){install.packages("extraDistr")}
    library(extraDistr)
    
    if(is.na(param1)){stop("this distribution requires a shape parameter (param1), often called the minimum value (as a proportion)")}
    if(is.na(param1)){stop("this distribution requires a second shape parameter (param2), often called the maximum value (as a proportion)")}
    if(is.na(param1)){stop("this distribution requires a third shape parameter (param3), often called the mode of the distribution (as a proportion)")}
    
    if(param1 < 0){warning("param1 (min) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param1 > 1){warning("param1 (min) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param2 < 0){warning("param2 (max) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param2 > 1){warning("param2 (max) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param3 < 0){warning("param3 (mode) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    if(param3 > 1){warning("param3 (mode) should be a proportion (0-1) of the active period. For instance, noon would equal 0.5 for a 24-hour period")}
    
    tprop <- as.vector(as.numeric(row.names(df))/nrow(df)) #create a column where t is shown as a proportion, not an integer
    
    triang.vect <- dtriang(x = tprop,
                           a = param1,
                           b = param2,
                           c = param3) 
    
    df$PDF <- triang.vect/sum(triang.vect)
    
  }
  
  #Convert times to be readable by R
  time.start <- as.POSIXct(x = time.start,
                           format = "%H:%M")
  
  time.stop <- as.POSIXct(x = time.stop,
                          format = "%H:%M")
  
  if(time.stop < time.start) {time.stop <- time.stop + 86400} #if the end time is smaller than the start time (midnight happens during diel period), assign the day after the initial time.
  if(time.stop == time.start) {time.stop <- time.stop + 86400} #if the end time is equal to the start time (24-hour period), assign the day after the initial time.
  
  time.quarter <- time.start + ((time.stop - time.start)/4) #determine what time occurs a quarter of the way between start and stop
  time.half <- time.start + ((time.stop - time.start)/2) #determine what time occurs half way between start and stop
  time.threequarter <- time.half + ((time.stop - time.start)/4) #determine what time occurs three quarters of the way between start and stop
  
  #Convert times to characters for the purpose of plotting
  plot.start.time <- as.character(substr(time.start, 12,16))
  plot.quarter.time <- as.character(substr(time.quarter, 12,16))
  plot.half.time <- as.character(substr(time.half, 12,16))
  plot.threequarter.time <- as.character(substr(time.threequarter, 12,16))
  plot.stop.time <- as.character(substr(time.stop, 12,16))
  
  max.height <- max(as.vector(c(df$observed, df$PDF))) #identify maximum height of either line on graph
  
  #Plot distribution for user
  plot(df$PDF ~ rownames(df),
       type = "l",
       main = "Idealized Distribution",
       sub = dist,
       xlab = "Time",
       ylab = "Probability Mass",
       xaxt="n", #don't plot any ticks or values on x axis (see next argument)
       col = "black",
       ylim = c(0,(max.height + 0.1*max.height))) 
  axis(side = 1, 
       at = c(0,(0.25*nrow(df)),(0.5*nrow(df)),(0.75*nrow(df)),nrow(df)),
       #at = c(0, 21, 42, 63, 84),
       labels = c(plot.start.time,
                  plot.quarter.time,
                  plot.half.time,
                  plot.threequarter.time,
                  plot.stop.time))
  lines(df$observed ~ rownames(df),
        type = "l",
        col = "red")
  legend(x = 0, y = max.height,
         legend = c("proposed", "observed"),
         col = c("black","red"),
         lty = 1,
         cex=.5)
  
  #Leave final distribution
  PDF <- df$PDF
  
} #end function

format_observedPDF <- function(df, sp.col, sp.name, time.col, time.start, time.stop, by = 10){#begin function
  
  #Error messages for start and stop times which are not character strings with 4 numbers split by a colon
  if(class(time.start) %in% "character" == FALSE){stop("time.start must be a character string in format: 'HH:MM'.")}
  if(class(time.stop) %in% "character" == FALSE){stop("time.stop must be a character string in format: 'HH:MM'.")}
  if(nchar(time.start) != 5){stop("time.start must be a character string in format: 'HH:MM'.")}
  if(nchar(time.stop) != 5){stop("time.stop must be a character string in format: 'HH:MM'.")}
  
  #Error message for time columns which are not formatted
  posit.test <- df[1,time.col] #isolate time column to test data class
  if(inherits(posit.test, "POSIXct")){stop("time.col must be class 'POSIXlt' or 'POSIXt'. Use the 'strptime' function to convert.")}
  
  #Error messages for disturbances not listed in disturbance column
  if(!(sp.name %in% df[,sp.col])){stop(paste("'", sp.name,"'", " is not contained within column titled: ", sp.col, sep=""))}
  
  #Load necessary package 'overlap'
  if(!("overlap" %in% as.vector(installed.packages()))){stop("This function requires the use of the 'overlap' package. Please install and load the 'overlap' package to proceed")} #check to see if the package is installed. If not, stop the function.
  library(overlap)
  
  #Format time to be readable by 'overlap' package
  temp.df <- df[which(df[,sp.col] %in% sp.name), ] #isolate data set to only the species of interest
  temp.df$time.char <- strftime(temp.df[,time.col]) #convert time to character
  temp.df$hour <- as.numeric(substr(temp.df$time.char, 12,13)) #isolate hour value
  temp.df$hour <- temp.df$hour * 60 #calculate how many minutes have passed based on hour value
  temp.df$minute <- as.numeric(substr(temp.df$time.char, 15,16)) #isolate true minute value
  temp.df$time.prop <- temp.df$hour + temp.df$minute #add minutes from hours to regular minutes
  temp.df$time.prop <- temp.df$time.prop/(24*60) #calculate time as a proportion (0-1)
  time.radians <- temp.df$time.prop * 2 * pi #convert time to radians
  
  #Calculate activity density and format for manipulation
  sp.density <- overlap::densityPlot(time.radians, #calculate activity density across 24-hour period
                                     adjust = .5) #adjusts to be more accurate to the data, rather than the smoothing
  sp.density <- sp.density[which(sp.density$x >= 0), ] #remove density outside the bounds of a 24-hour day
  sp.density <- sp.density[which(sp.density$x <= 24), ] #remove density outside the bounds of a 24-hour day
  
  #save start and stop time for later use
  plot.start.time <- time.start
  plot.stop.time <- time.stop
  
  #Determine lower and upper time limits based on user input
  start.hour <- as.numeric(substr(time.start,1,2)) #isolate starting hour
  start.min <- as.numeric(substr(time.start,4,5))/60 #isolate starting minute and convert it to a fraction of an hour
  time.start <- start.hour + start.min #make start time a decimal
  
  stop.hour <- as.numeric(substr(time.stop,1,2)) #isolate starting hour
  stop.min <- as.numeric(substr(time.stop,4,5))/60 #isolate starting minute and convert it to a fraction of an hour
  time.stop <- stop.hour + stop.min #make start time a decimal
  
  #Truncate mass to the lower and upper time limits specified by user
  if(time.start < time.stop){ #If starting time is less than stopping time...
    
    sp.density <- sp.density[which(sp.density$x >= time.start), ] #remove mass that occurs prior to the start time
    sp.density <- sp.density[which(sp.density$x <= time.stop), ] #remove mass that occurs after the stop time
    
    
  }
  
  if(time.start > time.stop){ #If starting time is less than stopping time...
    
    sp.density.beforemidnight <- sp.density[which(sp.density$x >= time.start), ] #remove mass that occurs prior to the start time
    sp.density.aftermidnight <- sp.density[which(sp.density$x <= time.stop), ] #remove mass that occurs after the stop time
    sp.density.aftermidnight$x <- sp.density.aftermidnight$x + 24 #adjust for stop time being smaller than start time
    
    sp.density <- rbind(sp.density.beforemidnight,#combine data sets
                        sp.density.aftermidnight)
    
  }
  
  if(time.start > time.stop){time.stop <- time.stop+24} #adjust stop time to be greater than start time, if needed.
  
  n.timesteps <- floor(((time.stop*60)-(time.start*60))/by) #Calculate how many timesteps there will be based on how time is binned
  
  #create an empty data frame to be filled with re-binned density value
  binned.mass <- data.frame()
  
  for(i in 1:n.timesteps){ #begin loop to re-bin density data
    
    #set start and stop times for bin i
    bin.start <- (time.start * 60) + ((i-1) * by) #identify minute in which bin should start
    bin.stop <- (time.start * 60) + (i * by) #identify minute in which bin should start
    
    #convert bin times to match sp.density format
    bin.start <- bin.start/60
    bin.stop <- bin.stop/60
    
    #isolate activity density within those times
    local.df <- sp.density[which(sp.density$x >= bin.start), ] #keep only records whose time is greater than local start time
    local.df <- local.df[which(local.df$x <= bin.stop), ] #keep only records whose time is less than local stop time
    
    sum.mass <- sum(local.df$y) #sum all activity mass present in bin
    if(is.na(sum.mass)){sum.mass <- 0} #if there was no mass to sum, make value 0 instead of NA
    
    #create row of data to add to data frame
    data.row <- data.frame(timestep = i,
                           mass = sum.mass)
    
    #add row to data frame
    binned.mass <- rbind(binned.mass, data.row)
    
  } #end loop to re-bin density data
  
  for(i in 1:nrow(binned.mass)){ #begin loop to adjust for inaccurate drops in activity when calculated mass is not within a time step bin
    
    if(i > 1){ #if looking at any time step greater than 1, impute missing value based on previous and next value
      if(binned.mass[i,"mass"] == 0){binned.mass[i,"mass"] <- mean(c(binned.mass[i-1,"mass"], binned.mass[i+1,"mass"]))}
    }
    
    if(i == nrow(binned.mass)){ #if looking at the last time step bin, impute based on previous value
      if(binned.mass[i,"mass"] == 0){binned.mass[i,"mass"] <- binned.mass[i-1,"mass"]}
    }
  } #end loop to adjust for inaccurate drops in activity when calculated mass is not within a time step bin
  
  binned.mass$mass <- binned.mass$mass/sum(binned.mass$mass) #Convert to proportion of activity mass (sums to 1)
  
  
  #convert start/stop times to be readable by R
  plot.start.time <- strptime(plot.start.time,
                              format = "%H:%M")
  plot.stop.time <- strptime(plot.stop.time,
                             format = "%H:%M")
  
  if(plot.stop.time < plot.start.time) {plot.stop.time <- plot.stop.time + 86400} #if the end time is smaller than the start time (midnight happens during diel period), assign the day after the initial time.
  if(plot.stop.time == plot.start.time) {plot.stop.time <- plot.stop.time + 86400} #if the end time is equal to the start time (24-hour period), assign the day after the initial time.
  
  #rename start/stop time to be consistent with other functions
  time.start <- plot.start.time
  time.stop <- plot.stop.time
  
  time.quarter <- time.start + ((time.stop - time.start)/4) #determine what time occurs a quarter of the way between start and stop
  time.half <- time.start + ((time.stop - time.start)/2) #determine what time occurs half way between start and stop
  time.threequarter <- time.half + ((time.stop - time.start)/4) #determine what time occurs three quarters of the way between start and stop
  
  #Convert times to characters for the purpose of plotting
  plot.start.time <- as.character(substr(time.start, 12,16))
  plot.quarter.time <- as.character(substr(time.quarter, 12,16))
  plot.half.time <- as.character(substr(time.half, 12,16))
  plot.threequarter.time <- as.character(substr(time.threequarter, 12,16))
  plot.stop.time <- as.character(substr(time.stop, 12,16))
  
  #Plot distribution for user
  plot(binned.mass$mass ~ binned.mass$timestep,
       type = "l",
       main = "Observed (smoothed) Distribution",
       xlab = "Time",
       ylab = "Activity Density",
       xaxt="n") #don't plot any ticks or values on x axis (see next line)
  axis(side = 1, 
       at = c(0,(0.25*nrow(binned.mass)),(0.5*nrow(binned.mass)),(0.75*nrow(binned.mass)),nrow(binned.mass)),
       #at = c(0, 21, 42, 63, 84),
       labels = c(plot.start.time,
                  plot.quarter.time,
                  plot.half.time,
                  plot.threequarter.time,
                  plot.stop.time))
  
  final <- as.vector(binned.mass$mass) #Save observed activity mass column as a vector and leave here for R to return
  
} #end function

estimate_activity <- function(idealPDF, disturbance, observedPDF, confInt = FALSE, nIter = NA, save.boots = FALSE, AIC = FALSE){#begin function
  
  #ERROR MESSAGES AND WARNINGS
  pdf.length <- as.vector(c(length(idealPDF))) #determine length of probability mass function
  
  for(i in 1:length(disturbance)){ #check for mismatch
    
    disturbance.length <- length(disturbance[[i]]) #determine length of disturbance data
    
    #set acceptable lengths for disturbance data based on PDF data (allowing for binning errors)
    margin.of.error <- c(pdf.length-2,
                         pdf.length-1,
                         pdf.length,
                         pdf.length+1,
                         pdf.length+2)
    
    if(!(pdf.length %in% margin.of.error)){stop("idealPDF and at least one disturbance vector are of different lengths.")} #throw error if mismatach
  }
  
  observed.PDF.length <- length(observedPDF)
  if(!(observed.PDF.length %in% margin.of.error)){stop("idealPDF and observedPDF are of different lengths.")} #throw error if mismatach
  
  if(confInt == TRUE & is.na(nIter)){stop("number of iterations (nIter) must be spcecified to estimate confidence intervals (confInt = TRUE)")}
  if(confInt == TRUE){message("please note that estimation of confidence intervals via bootstrapping can be slow - progress will be indicated")}
  
  if(length(disturbance) == 2){message("note that optimization is finicky for two-predator analyses and may fail during bootstrapping")}
  if(length(disturbance) > 2){stop("while several predators are possible, this code only allows for two predators at a time (and often fails any time there is more than one predator/disturbance)")}
  
  #pack data into a list to pass to function
  data.for.function <- list()
  data.for.function[[1]] <- idealPDF
  data.for.function[[2]] <- disturbance
  data.for.function[[3]] <- observedPDF
  
  #data <- data.for.function
  
  min.LL.function <- function(data, par){
    
    #put searching param values into a list for compatability with prediction function
    par.list <- list()
    par.list[[1]] <- par[1]
    if(length(par) == 2){par.list[[2]] <- par[2]}
    
    predict_activity(idealPDF = data[[1]],
                     disturbance = data[[2]],
                     sensitivity = par.list,
                     output = "LL",
                     observedPDF = data[[3]],
                     plot = FALSE)
    
  }
  
  
  if(length(disturbance) == 1){ #begin optimization if only one predator
    results <- optim(par = runif(n = 1, min=0, max=1), #assign a random starting value
                     fn = min.LL.function,
                     data = data.for.function,
                     method = "Brent", #this method allows limits on estimated parameter values but facilitates only one paramater (but is more stable in estimation than L-BFGS-B method)
                     lower = 0,
                     upper = 1)
  } #end optimization if only one predator
  
  if(length(disturbance) == 2){ #begin optimization if two predators
    results <- optim(par = c(runif(n = 1, min=0, max=1),runif(n = 1, min=0, max=1)), #assign a random starting value
                     fn = min.LL.function,
                     data = data.for.function,
                     method = "L-BFGS-B", #this method allows limits on estimated parameter values and facilitates multiple paramaters
                     lower = c(0,0),
                     upper = c(1,1))
  } #end optimization if two predators
  
  #bootstrap to estimate confidence intervals
  if(confInt == TRUE) {#begin bootstrapping instructions
    
    all.indices <- as.vector(rownames(as.data.frame(observedPDF))) #create a vector of indices for data distributions
    bootstrap.vals <- data.frame() #create empty dataframe to fill
    
    for(i in 1:nIter){ #begin looping through iterations
      
      #determine which steps will be retained in this run
      keep.indices <- sample(x = all.indices,
                             size = floor(0.8 * length(all.indices)),
                             replace = FALSE)
      keep.indices <- sort(as.numeric(keep.indices)) #re-order steps numerically
      
      #thin data according to indices
      bsIdeal <- idealPDF[keep.indices]
      bsDisturbance <- list()
      bsDisturbance[[1]] <- disturbance[[1]][keep.indices]
      if(length(disturbance) == 2){bsDisturbance[[1]] <- disturbance[[2]][keep.indices]}
      bsObserved <- observedPDF[keep.indices]
      
      #pack data into a list to pass to function
      data.for.bs <- list()
      data.for.bs[[1]] <- bsIdeal
      data.for.bs[[2]] <- bsDisturbance
      data.for.bs[[3]] <- bsObserved
      
      if(length(disturbance) == 1){ #begin optimization if only one predator
        bs.results <- optim(par = runif(n = 1, min=0, max=1), #assign a random starting value
                            fn = min.LL.function,
                            data = data.for.bs,
                            method = "L-BFGS-B",
                            lower = 0,
                            upper = 1)
        
        local.results <- data.frame(run = i,
                                    S = bs.results$par)
      } #end optimization if only one predator
      
      if(length(disturbance) == 2){ #begin optimization if two predators
        bs.results <- optim(par = c(runif(n = 1, min=0, max=1),runif(n = 1, min=0, max=1)), #assign a random starting value
                            fn = min.LL.function,
                            data = data.for.function,
                            method = "L-BFGS-B",
                            lower = 0,
                            upper = 1)
        
        local.results <- data.frame(run = i,
                                    S1 = bs.results$par[1],
                                    S2 = bs.results$par[2])
      } #end optimization if two predators
      
      #add estimated values to the data frame
      bootstrap.vals <- rbind(bootstrap.vals, local.results)
      
      #print progress message
      print(paste(i," of ", nIter," iterations complete (", round((i/nIter*100),2),"%)", sep=""))
    } #end looping through iterations
    
    if(length(disturbance) == 1){ #begin summarizing from bootstraped iterations (single predator)
      
      S.CI <- quantile(x = bootstrap.vals$S,
                       probs = c(0.025, 0.975))
      
      report <- data.frame(S = results$par[1],
                           CI95 = paste("(",round(S.CI[1],4)," - ",round(S.CI[2],4),")", sep=""),
                           negLL = results$value)
      
    } #end summarizing from bootstraped iterations (single predator)
    
    if(length(disturbance) == 2){ #begin summarizing from bootstraped iterations (multiple predators)
      
      S1.CI <- quantile(x = bootstrap.vals$S1,
                        probs = c(0.025, 0.975))
      
      S2.CI <- quantile(x = bootstrap.vals$S2,
                        probs = c(0.025, 0.975))
      
      report <- data.frame(S1 = results$par[1],
                           CI95.1 = paste("(",round(S1.CI[1],4)," - ",round(S1.CI[2],4),")", sep=""),
                           S2 = results$par[2],
                           CI95.2 = paste("(",round(S2.CI[1],4)," - ",round(S2.CI[2],4),")", sep=""),
                           negLL = results$value)
      
    } #end summarizing from bootstraped iterations (multiple predators)
  } #end bootstrapping instructions
  
  if(confInt == FALSE){ #format results for user without confidence intervals if requested
    
    if(length(disturbance)==1){ #formatting specific to a single predator
      report <- data.frame(S = results$par[1],
                           negLL = results$value)
    }
    
    if(length(disturbance)==2){ #formatting specific to a single predator
      report <- data.frame(S1 = results$par[1],
                           S2 = results$par[2],
                           negLL = results$value)
    }
  } #end re-formatting results without confidence intervals
  
  if(AIC == TRUE){ #begin AIC calculation
    
    #note that these are specific to the present formulat. If future user add parameters then these will need to be updated
    if(length(disturbance) == 1){n.params <- 9}
    if(length(disturbance) == 2){n.params <- 11}
    
    AIC.score <- (2*n.params) - (2*(results$value)) #AIC = 2K - 2ln(L)
    
    report$AIC <- AIC.score #add AIC score to the report
    
  } #end AIC calculation
  
  if(save.boots == TRUE){
    bootstrap.vals <<- bootstrap.vals
    message("bootstrap values ('bootstrap.vals') have been saved to your environment")
  }
  #remove row number
  rownames(report) <- ""
  
  #print report in console
  print(report)
  
  estimated.vals <- report #leave here for the user
} #end function


predict_activity <- function(idealPDF, disturbance, sensitivity, time.start = NULL, time.stop = NULL, plot=TRUE, output = "complete", observedPDF = NULL){#begin function
  
  #ERROR MESSAGES
  pdf.length <- as.vector(c(length(idealPDF))) #determine length of probability mass function
  
  for(i in 1:length(disturbance)){ #check for mismatch
    
    disturbance.length <- length(disturbance[[i]]) #determine length of disturbance data
    
    #set acceptable lengths for disturbance data based on PDF data (allowing for binning errors)
    margin.of.error <- c(pdf.length-2, 
                         pdf.length-1, 
                         pdf.length,
                         pdf.length+1,
                         pdf.length+2)
    
    if(!(pdf.length %in% margin.of.error)){stop("idealPDF and at least one disturbance vector are of different lengths.")} #throw error if mismatach
  }
  
  for(i in 1:length(sensitivity)){ #check for incorrectly formatted sensitivity data
    
    sensitivity.length <- length(sensitivity[[i]])
    
    if(sensitivity.length != 1){stop("sensitivity has a length greater than one. Ensure that only one sensitivity value is provided per disturbance, and each value holds a different list slot.")} #throw error if mismatach
  }
  
  for(i in 1:length(sensitivity)){ #check for incorrectly formatted sensitivity data
    
    sensitivity.value <- sensitivity[[i]]
    
    if(sensitivity.value < 0 ){stop("sensitivity values may range between zero and one. Ensure your data are within these limits.")} #throw error if S < 0 
    if(sensitivity.value > 1 ){stop("sensitivity values may range between zero and one. Ensure your data are within these limits.")} #throw error if S > 1 
  }
  
  if(!output %in% c("complete","LL")){stop("output must be either 'complete' (default) or 'LL' which returns negative log likelihood")}
  if(output %in% "LL" & is.null(observedPDF)){stop("observedPDF is required to calculate and maximize likelihood")} #throw error if observed data is not included but LL is the objective
  if(output %in% "LL" & !(pdf.length == length(observedPDF))){stop("idealPDF and observedPDF are of different lengths.")} #throw error if mismatach
  
  #INITIAL VALUES
  df <- data.frame() #create empty data frame to fill with values
  Lt <- 0 #at the start, there is no lost activity
  bigT <- length(disturbance[[1]]) #identify how many time steps are in the period
  
  #Here note that both T and t have meanings in R, so we use bigT in place of T, and i in place of t
  #However, we still use t as a subscript in objecting naming
  
  for(i in 1:bigT){ #begin loop for each time step
    
    cumulative.pressure.t <- vector() #create an empty vector to fill with predation activity loss values
    cumulative.compensated.t <- vector()
    dDot <- idealPDF[i] #identify the expected activity density at time step i based on the ideal distribution
    
    for(j in 1:length(disturbance)){#start loop for each disturbance
      
      #LOST ACTIVITY
      
      a.t <- disturbance[[j]][i] / max(disturbance[[j]]) #convert predator probability mass to a proportion of the maximum predator activity
      pt <- sensitivity[[j]] * a.t #calculate pt (disturbance pressure)
      ltj <- pt * dDot #calculate lt from disturbance j (activity lost) 
      
      cumulative.pressure.t <- c(cumulative.pressure.t, ltj) #add local activity loss due to predator j at time i to the vector of activity lost for time step i
      
      Dtj <- 1 - a.t #determines the severity of disturbance at current time step to inform level of activity compensation
      
      ctj <- Dtj * Lt * (i/bigT) #calculate ct at time t for disturbance j based on missing activity at time t (Lt), proportion of the active period that has passed (t/T), and severity of threat
      
      cumulative.compensated.t <- c(cumulative.compensated.t, ctj) #add local activity loss due to predator j at time i to the vector of activity lost for time step i
    } #end loop for each disturbance
    
    #UPDATE lt and ct
    
    lt <- sum(cumulative.pressure.t) #add activity densities lost from each disturbance (if multiple)
    ct <- sum(cumulative.compensated.t) #add compensated activity
    
    
    #COMPENSATED ACTIVITY
    
    #determine realized activity level 
    dRealized <- dDot - lt + ct
    if(dRealized < 0){dRealized <- 0} #if dRealized is less than 0, which is biologically impossible, keep it at zero.
    
    Lt <- Lt + lt - ct #Calculate total lost activity (Lt) by adding lost activity at time t and subtracting compensated activity at time t
    #Note that this will not be used again until the next time step, avoiding a circular equation.
    
    #SAVE RESULTS
    temp.df <- data.frame(timestep = i,
                          ideal = dDot,
                          Lt = Lt,
                          lt = lt,
                          ct = ct,
                          activity.dens = dRealized)
    
    df <- rbind(df, temp.df) #Add this row to the dataframe
    
  } #end loop for time step i
  
  df <- as.data.frame(df) #save realized activity curve
  
  if(output %in% "LL"){ #begin calculating likelihood
    
    #attach observed values to prediction
    df$observed <- observedPDF
    
    #calculate residuals based on candidate sensitivity value
    residuals <- df$observed - df$activity.dens
    
    #calculate negative log likelihood assuming a normal distribution of residuals
    negLL <- -1*(sum(dnorm(x = residuals,
                           mean = 0,
                           sd = sd(residuals),
                           log = TRUE)))
    
  }#end calculating likelihood
  
  #PLOTTING
  if(plot == TRUE){ #begin plotting instructions
    df$activity.dens <- ifelse(df$activity.dens <0, 0, df$activity.dens) #if density goes negative, keep limit it to zero
    
    max.prob <- max(c(df$activity.dens, df$ideal))+(.33*max(c(df$activity.dens, df$ideal))) #calculate a y limit so the graph stops cutting off data and leaves room for a legend
    
    if(!is.null(time.start)){
      
      #Convert times to be readable by R
      time.start <- as.POSIXct(x = time.start,
                               format = "%H:%M")
      
      time.stop <- as.POSIXct(x = time.stop,
                              format = "%H:%M")
      
      if(time.stop < time.start) {time.stop <- time.stop + 86400} #if the end time is smaller than the start time (midnight happens during diel period), assign the day after the initial time.
      if(time.stop == time.start) {time.stop <- time.stop + 86400} #if the end time is equal to the start time (24-hour period), assign the day after the initial time.
      
      time.quarter <- time.start + ((time.stop - time.start)/4) #determine what time occurs a quarter of the way between start and stop
      time.half <- time.start + ((time.stop - time.start)/2) #determine what time occurs half way between start and stop
      time.threequarter <- time.half + ((time.stop - time.start)/4) #determine what time occurs three quarters of the way between start and stop
      
      #Convert times to characters for the purpose of plotting
      plot.start.time <- as.character(substr(time.start, 12,16))
      plot.quarter.time <- as.character(substr(time.quarter, 12,16))
      plot.half.time <- as.character(substr(time.half, 12,16))
      plot.threequarter.time <- as.character(substr(time.threequarter, 12,16))
      plot.stop.time <- as.character(substr(time.stop, 12,16))
      
      #Initiate plot with custom x axis, as defined by time.start and time.stop
      plot(df$ideal ~ df$timestep, #plot expected, ideal distribution
           type = "l", #line graph
           col ="#EDEDED", #color the line black
           xlab="Timestep", #x axis label
           ylab="Activity Density", #y axis label
           ylim=c(0,max.prob), #upper limit of y axis
           lwd=2,#line size 2
           xaxt="n") #don't include axes
      axis(side = 1, 
           at = c(0,(0.25*nrow(df)),(0.5*nrow(df)),(0.75*nrow(df)),nrow(df)),
           labels = c(plot.start.time,
                      plot.quarter.time,
                      plot.half.time,
                      plot.threequarter.time,
                      plot.stop.time))
    }
    
    if(is.null(time.start)){
      
      #Initiate plot with standard axis if time.start and time.stop are not provided
      plot(df$ideal ~ df$timestep, #plot expected, ideal distribution
           type = "l", #line graph
           col ="#EDEDED", #color the line black
           xlab="Time Step", #x axis label
           ylab="Activity Density", #y axis label
           ylim=c(0,max.prob), #upper limit of y axis
           lwd=2)#line size 2
    }
    
    #Plot raw realized activity
    lines(x=df$timestep, 
          y=df$activity.dens,
          type = "l", 
          col="red", 
          lwd=2)
    
    legend(x = 0, y = max.prob,
           legend = c("ideal","predicted"),
           col = c("#EDEDED","red"),
           lty = 1,
           cex=.5)
    
  } #end plotting instructions
  
  if(output %in% "complete"){final <- df} #leaving here for the user, if they so desire
  if(output %in% "LL"){final <- negLL}
  #output <- ifelse(output == "LL", negLL, output) #this seemed to fix the error from above lines
  
  predicted.vals <- final #leaving this here for the user
  
} #end function

UtahData <- function(){#begin data function
  
  sample.data <- data.frame(project = "Snapshot USA 2019",
             array = "Wasatch Wildlife Watch",
             site = c("Creekside Park.USA.2019",       "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019"   ,    "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019" ,      "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019"  ,     "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019" ,      "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019"  ,     "Creekside Park.USA.2019"  ,    
                      "Creekside Park.USA.2019"   ,    "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"    ,   "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"   ,    "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"    ,   "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"   ,    "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"    ,   "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"     ,  "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019",       "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019" ,      "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"  ,     "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"   ,    "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019"    ,   "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019",       "Crestwood Park.USA.2019"  ,    
                      "Crestwood Park.USA.2019" ,      "Crestwood Park.USA.2019"  ,    
                      "DD_01.USA.2019"           ,     "DD_01.USA.2019"           ,    
                      "DD_01.USA.2019"            ,    "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019",                "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019" ,               "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"  ,              "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"   ,             "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"    ,            "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"     ,           "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"      ,          "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"       ,         "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"        ,        "DD_02.USA.2019"           ,    
                      "DD_02.USA.2019"         ,       "DD_02.USA.2019"           ,    
                      "DD_03.USA.2019",                "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019" ,               "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"  ,              "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"   ,             "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"    ,            "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"     ,           "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"      ,          "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"       ,         "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"        ,        "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"         ,       "DD_03.USA.2019"           ,    
                      "DD_03.USA.2019"          ,      "DD_03.USA.2019"           ,    
                      "DD_04.USA.2019"           ,     "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"            ,    "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"             ,   "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019",                "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019" ,               "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"  ,              "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"   ,             "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"    ,            "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"     ,           "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"      ,          "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"       ,         "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"        ,        "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"         ,       "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"          ,      "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"           ,     "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"            ,    "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"             ,   "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019",                "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019" ,               "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"  ,              "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"   ,             "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"    ,            "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"     ,           "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"      ,          "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"       ,         "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"        ,        "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"         ,       "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"          ,      "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"           ,     "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"            ,    "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019",                "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019" ,               "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"  ,              "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"   ,             "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"    ,            "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"     ,           "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"      ,          "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"       ,         "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"        ,        "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"         ,       "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"          ,      "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"           ,     "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"            ,    "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019",                "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019" ,               "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"  ,              "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"   ,             "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"    ,            "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"     ,           "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"      ,          "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"       ,         "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"        ,        "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"         ,       "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"          ,      "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"           ,     "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"            ,    "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019",                "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019" ,               "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"  ,              "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"   ,             "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"    ,            "DD_04.USA.2019"           ,    
                      "DD_04.USA.2019"     ,           "DD_07.USA.2019"           ,    
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019" , "General Holm Park N.USA.2019", 
                      "General Holm Park N.USA.2019",  "General Holm Park S.USA.2019", 
                      "General Holm Park S.USA.2019",  "General Holm Park S.USA.2019", 
                      "General Holm Park S.USA.2019",  "General Holm Park S.USA.2019", 
                      "General Holm Park S.USA.2019",  "General Holm Park S.USA.2019", 
                      "General Holm Park S.USA.2019",  "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019",                "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019" ,               "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"  ,              "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"   ,             "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"    ,            "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"     ,           "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"      ,          "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"       ,         "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"        ,        "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"         ,       "JR_01.USA.2019"             ,  
                      "JR_01.USA.2019"          ,      "JR_02.USA.2019"             ,  
                      "JR_02.USA.2019"           ,     "JR_02.USA.2019"             ,  
                      "JR_02.USA.2019"            ,    "JR_02.USA.2019"             ,  
                      "JR_02.USA.2019"             ,   "JR_02.USA.2019"             ,  
                      "JR_02.USA.2019"              ,  "JR_02.USA.2019"             ,  
                      "JRP_02.USA.2019" ,              "JRP_02.USA.2019"            ,  
                      "JRP_02.USA.2019"  ,             "JRP_02.USA.2019"             , 
                      "JRP_02.USA.2019"   ,            "JRP_02.USA.2019"            ,  
                      "JRP_02.USA.2019"    ,           "JRP_02.USA.2019"            ,  
                      "JRP_02.USA.2019"     ,          "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"      ,         "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"       ,        "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"        ,       "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"         ,      "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"          ,     "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"           ,    "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"            ,   "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019" ,              "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"  ,             "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"   ,            "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"    ,           "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"     ,          "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"      ,         "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"       ,        "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"        ,       "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"         ,      "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"          ,     "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"           ,    "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019" ,              "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"  ,             "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"   ,            "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"    ,           "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"     ,          "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"      ,         "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"       ,        "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"        ,       "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"         ,      "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"          ,     "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"           ,    "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019" ,              "LCC_USA_01.2019"            ,  
                      "LCC_USA_01.2019"  ,             "Little Confluence.USA.2019" ,  
                      "Little Confluence.USA.2019",    "Little Confluence.USA.2019" ,  
                      "Little Confluence.USA.2019",    "Little Confluence.USA.2019" ,  
                      "Little Confluence.USA.2019",    "Little Confluence.USA.2019" ,  
                      "Little Confluence.USA.2019",    "Little Confluence.USA.2019" ,  
                      "Oxbow.USA.2019"            ,    "Parley's Nature Park.USA.2019",
                      "Parley's Nature Park.USA.2019", "Parley's Nature Park.USA.2019",
                      "Parley's Nature Park.USA.2019", "Parley's Nature Park.USA.2019",
                      "Redwood Natural Area.USA.2019", "Redwood Natural Area.USA.2019",
                      "Redwood Natural Area.USA.2019", "Redwood Natural Area.USA.2019",
                      "Redwood Natural Area.USA.2019", "Redwood Natural Area.USA.2019",
                      "Redwood Natural Area.USA.2019", "Redwood Natural Area.USA.2019",
                      "Redwood Natural Area.USA.2019", "SLC Cemetery.USA.2019"        ,
                      "SLC Cemetery.USA.2019"       ,  "WF.USA.2019"                  ,
                      "WF.USA.2019"       ,            "WF.USA.2019"                  ,
                      "WF.USA.2019"        ,           "WF.USA.2019"                  ,
                      "WF.USA.2019"         ,          "WF.USA.2019"                  ,
                      "WF.USA.2019"          ,         "WF.USA.2019"                  ,
                      "WF.USA.2019"           ,        "WF.USA.2019"                  ,
                      "WF.USA.2019"            ,       "WF.USA.2019"                  ,
                      "WF.USA.2019"             ,      "WF.USA.2019"                  ,
                      "WF.USA.2019"              ,     "WF.USA.2019"                  ,
                      "WF.USA.2019"               ,    "WF.USA.2019"                  ,
                      "WF.USA.2019"                ,   "WF.USA.2019"                  ,
                      "WF.USA.2019" ,                  "WF.USA.2019"                  ,
                      "WF.USA.2019"  ,                 "WF.USA.2019"                  ,
                      "WF.USA.2019"   ,                "WF.USA.2019"                  ,
                      "WF.USA.2019"    ,               "WF.USA.2019"                  ,
                      "WF.USA.2019"     ,              "WF.USA.2019"                  ,
                      "WF.USA.2019"      ,             "WF.USA.2019"                  ,
                      "WF.USA.2019"       ,            "WF.USA.2019"                  ,
                      "WF.USA.2019"        ,           "WF.USA.2019"                  ,
                      "WF.USA.2019"         ,          "WF.USA.2019"                  ,
                      "WF.USA.2019"          ,         "WF.USA.2019"                  ,
                      "WF.USA.2019"           ,        "WF.USA.2019"                  ,
                      "WF.USA.2019"            ,       "WF.USA.2019"                  ,
                      "WF.USA.2019"             ,      "WF.USA.2019"                  ,
                      "WF.USA.2019"              ,     "WF.USA.2019"                  ,
                      "WF.USA.2019"               ,    "WF.USA.2019"                  ,
                      "WF.USA.2019",                   "WF.USA.2019"                  ,
                      "WF.USA.2019" ,                  "WF.USA.2019"                  ,
                      "WF.USA.2019"  ,                 "WF.USA.2019"                  ,
                      "WF.USA.2019"   ,                "WF.USA.2019"                  ,
                      "WF.USA.2019"    ,               "WF.USA.2019"                  ,
                      "WF.USA.2019"     ,              "WF.USA.2019"                  ,
                      "WF.USA.2019"      ,             "WF.USA.2019"                  ,
                      "WF.USA.2019"       ,            "WF.USA.2019"                  ,
                      "WF.USA.2019"        ,           "WF.USA.2019"                  ,
                      "WF.USA.2019"         ,          "WF.USA.2019"                  ,
                      "WF.USA.2019"          ,         "WF.USA.2019"                  ,
                      "WF.USA.2019"),
             obs.time = c( "9/11/19 3:07",   "9/13/19 2:17"  , "9/16/19 0:35",  
                           "9/16/19 2:29" ,  "9/16/19 5:33" ,  "9/18/19 3:32",  
                          "9/18/19 21:22"  ,"9/19/19 3:09"  , "9/20/19 22:49", 
                           "9/20/19 22:50",  "9/26/19 0:32" ,  "9/27/19 21:13" ,
                          "9/27/19 22:42" , "9/9/19 3:11"   , "9/10/19 3:51",  
                        "9/11/19 21:14"  ,"9/12/19 0:21"   ,"9/12/19 0:31",  
                           "9/13/19 2:07",   "9/13/19 3:42" ,  "9/13/19 21:53", 
                           "9/13/19 22:29",  "9/14/19 1:16" ,  "9/14/19 4:29" , 
                           "9/14/19 6:25"  , "9/15/19 4:27" ,  "9/17/19 23:01", 
                           "9/18/19 0:51"  , "9/18/19 4:13" ,  "9/20/19 2:10" , 
                           "9/21/19 23:24"  ,"9/22/19 0:52" ,  "9/22/19 3:57" , 
                           "9/23/19 2:48"   ,"9/26/19 23:03",  "9/27/19 1:25" , 
                           "9/27/19 2:22",   "9/28/19 2:00" ,  "9/29/19 1:00" , 
                           "9/29/19 1:57" ,  "10/3/19 0:34" ,  "10/3/19 21:46", 
                           "10/9/19 22:52" , "10/9/19 4:10" ,  "10/11/19 21:27",
                           "10/11/19 23:36" ,"10/11/19 23:48", "10/13/19 2:34" ,
                           "10/13/19 2:56",  "10/15/19 21:47", "10/16/19 5:55", 
                           "10/16/19 7:16" , "10/19/19 7:47",  "10/21/19 1:33", 
                           "10/22/19 20:08", "10/27/19 20:32", "11/2/19 3:06" , 
                           "11/2/19 15:10",  "11/5/19 21:24",  "11/7/19 21:35", 
                           "11/9/19 7:56" ,  "11/13/19 3:35",  "11/16/19 20:45",
                           "11/17/19 10:05", "10/5/19 20:26",  "10/6/19 0:03" , 
                           "10/8/19 0:41" ,  "10/9/19 20:45",  "10/10/19 23:47",
                           "10/14/19 7:15",  "10/17/19 3:41",  "10/17/19 20:01",
                           "10/22/19 0:59",  "10/23/19 10:25", "10/25/19 4:59" ,
                           "10/26/19 0:34",  "10/27/19 9:19",  "10/29/19 19:56",
                           "11/4/19 7:39" ,  "11/4/19 21:12",  "11/7/19 5:42"  ,
                           "11/8/19 20:11",  "11/9/19 4:17" ,  "11/10/19 21:54",
                           "11/14/19 19:49", "11/14/19 19:51", "10/3/19 8:27"  ,
                           "10/3/19 19:20",  "10/4/19 2:58" ,  "10/4/19 19:50" ,
                           "10/5/19 5:07" ,  "10/5/19 18:30" , "10/6/19 3:46"  ,
                           "10/8/19 5:32" ,  "10/8/19 6:28" ,  "10/9/19 22:50" ,
                           "10/9/19 22:54",  "10/9/19 23:12",  "10/10/19 6:46" ,
                           "10/10/19 7:51",  "10/10/19 21:19", "10/10/19 23:13",
                           "10/11/19 1:08",  "10/11/19 1:10",  "10/11/19 2:21" ,
                           "10/11/19 11:00", "10/11/19 11:02", "10/11/19 21:48",
                           "10/11/19 22:18", "10/11/19 22:33", "10/11/19 22:38",
                           "10/12/19 2:01",  "10/12/19 2:02",  "10/12/19 20:30",
                           "10/12/19 23:45", "10/13/19 5:28",  "10/13/19 5:30" ,
                           "10/13/19 7:16"  ,"10/18/19 0:08",  "10/18/19 6:01" ,
                           "10/18/19 6:13",  "10/19/19 21:02", "10/20/19 7:23" ,
                           "10/21/19 4:03" , "10/21/19 4:20",  "10/21/19 23:57",
                           "10/22/19 6:20"  ,"10/22/19 6:21",  "10/22/19 19:26",
                           "10/23/19 21:10", "10/24/19 11:51", "10/24/19 19:15",
                           "10/25/19 3:52"  ,"10/25/19 4:10",  "10/25/19 9:16" ,
                           "10/25/19 10:05", "10/25/19 21:21", "10/26/19 2:07" ,
                          "10/26/19 4:12",  "10/26/19 7:31",  "10/27/19 5:06" ,
                           "10/27/19 22:42", "10/27/19 23:13", "10/27/19 23:36",
                           "10/28/19 4:05",  "10/29/19 5:07",  "10/30/19 0:53" ,
                           "10/30/19 3:31",  "10/31/19 0:29",  "10/31/19 12:18",
                           "11/1/19 5:28" ,  "11/1/19 17:49",  "11/2/19 3:17"  ,
                           "11/3/19 2:31" ,  "11/3/19 13:07",  "11/5/19 20:46" ,
                           "11/5/19 21:20",  "11/5/19 22:30",  "11/6/19 10:15" ,
                           "11/7/19 6:14" ,  "11/7/19 7:26" ,  "11/7/19 7:28"  ,
                           "11/7/19 22:07",  "11/7/19 22:33",  "11/8/19 0:11"  ,
                           "11/8/19 0:56" ,  "11/8/19 6:32" ,  "11/8/19 19:14" ,
                           "11/8/19 19:18",  "11/8/19 19:33",  "11/8/19 19:39" ,
                           "11/8/19 19:57",  "11/9/19 0:44" ,  "11/9/19 3:18"  ,
                           "11/9/19 7:08" ,  "11/9/19 20:10",  "11/9/19 20:16" ,
                           "11/9/19 20:23",  "11/9/19 20:37",  "11/9/19 20:42" ,
                           "11/9/19 20:43",  "11/9/19 20:46",  "11/9/19 21:06" ,
                          "11/9/19 22:36" , "11/9/19 22:54" , "11/9/19 23:48" ,
                           "11/9/19 23:51",  "11/10/19 0:01",  "11/10/19 4:12" ,
                           "11/10/19 20:30", "11/10/19 21:53", "11/10/19 22:09",
                           "11/11/19 3:22",  "11/11/19 3:29",  "11/11/19 7:36" ,
                           "11/12/19 1:06",  "11/12/19 2:03",  "11/12/19 21:40",
                           "11/12/19 21:41", "11/13/19 0:14",  "11/14/19 0:44" ,
                           "11/14/19 11:45", "11/15/19 1:07",  "11/15/19 7:10" ,
                           "11/15/19 7:11",  "11/15/19 12:07", "11/16/19 2:01" ,
                           "11/16/19 4:33" , "11/16/19 5:49",  "11/14/19 0:47" ,
                           "9/7/19 3:39"    ,"9/9/19 4:27"  ,  "9/12/19 5:40"  ,
                           "9/14/19 6:01" ,  "9/18/19 5:25" ,  "9/21/19 3:44"  ,
                           "9/21/19 6:38" ,  "9/21/19 6:45" ,  "9/21/19 6:50"  ,
                           "9/22/19 1:30" ,  "9/22/19 5:15" ,  "9/22/19 6:34"  ,
                           "9/30/19 0:38" ,  "9/30/19 4:39" ,  "9/30/19 4:59"  ,
                           "9/6/19 3:11"  ,  "9/8/19 5:15"  ,  "9/9/19 2:14"   ,
                           "9/10/19 2:36" ,  "9/11/19 20:26",  "9/14/19 1:45"  ,
                           "9/18/19 2:11" ,  "9/22/19 2:06" ,  "9/10/19 1:13"  ,
                           "9/11/19 1:20" ,  "9/11/19 4:36" ,  "9/11/19 6:12"  ,
                           "9/11/19 22:49",  "9/12/19 20:08",  "9/12/19 22:30" ,
                           "9/13/19 4:05"  , "9/13/19 20:48",  "9/14/19 7:03"  ,
                           "9/14/19 19:48",  "9/17/19 4:13" ,  "9/17/19 5:12"  ,
                           "9/19/19 4:20" ,  "9/20/19 1:08" ,  "9/23/19 6:10"  ,
                           "9/23/19 11:05",  "9/24/19 5:41" ,  "9/25/19 23:02" ,
                           "9/27/19 19:51",  "9/28/19 4:57" ,  "9/28/19 22:05" ,
                           "9/8/19 22:04" ,  "9/8/19 22:36" ,  "9/9/19 13:46"  ,
                           "9/16/19 5:24" ,  "9/21/19 3:06" ,  "9/25/19 2:30"  ,
                           "9/27/19 0:08" ,  "9/28/19 6:12" ,  "9/28/19 20:09" ,
                           "9/12/19 11:23",  "9/13/19 20:38",  "9/14/19 6:42"  ,
                           "9/14/19 20:08",  "9/16/19 20:08",  "9/19/19 6:33"  ,
                           "9/24/19 19:42",  "9/25/19 19:30",  "9/28/19 5:27"  ,
                           "10/4/19 23:49",  "10/4/19 23:58",  "10/5/19 0:49"  ,
                           "10/5/19 4:50" ,  "10/6/19 4:16" ,  "10/6/19 4:58"  ,
                           "10/6/19 20:51",  "10/6/19 20:53",  "10/6/19 21:31" ,
                           "10/6/19 23:57",  "10/7/19 4:25" ,  "10/7/19 6:12"  ,
                           "10/7/19 11:10",  "10/7/19 23:51",  "10/8/19 0:26"  ,
                           "10/8/19 5:00" ,  "10/8/19 20:09",  "10/8/19 20:46" ,
                           "10/8/19 23:39",  "10/9/19 1:27" ,  "10/9/19 1:59"  ,
                           "10/10/19 4:10",  "10/10/19 6:01",  "10/10/19 19:29",
                           "10/10/19 20:12", "10/10/19 20:49", "10/11/19 1:34" ,
                           "10/11/19 3:52",  "10/11/19 6:20",  "10/12/19 1:32" ,
                           "10/12/19 3:06",  "10/12/19 5:09",  "10/12/19 6:03" ,
                           "10/13/19 2:03",  "10/13/19 4:35",  "10/13/19 5:30" ,
                           "10/14/19 0:09",  "10/14/19 0:48",  "10/14/19 11:44",
                           "10/14/19 20:28", "10/14/19 22:46", "10/15/19 0:51" ,
                           "10/15/19 0:56",  "10/15/19 4:50",  "10/15/19 5:27" ,
                           "10/15/19 21:23", "10/15/19 23:11", "10/15/19 23:42",
                           "10/16/19 0:54",  "10/16/19 3:30",  "10/17/19 3:35" ,
                           "10/17/19 5:46",  "10/17/19 22:35", "10/18/19 1:38" ,
                           "10/18/19 4:16",  "10/18/19 19:40", "10/18/19 19:41",
                           "10/19/19 6:23",  "10/20/19 7:37",  "10/21/19 3:05" ,
                           "10/21/19 19:19", "10/23/19 19:27", "10/7/19 0:15"  ,
                           "10/11/19 1:55",  "10/11/19 1:59",  "10/23/19 23:09",
                           "10/27/19 4:32",  "10/31/19 5:30",  "11/4/19 4:43"  ,
                           "11/8/19 4:03" ,  "11/8/19 4:12" ,  "9/8/19 4:36"   ,
                           "9/22/19 0:52" ,  "9/22/19 1:57" ,  "9/24/19 1:45"  ,
                           "9/24/19 21:21",  "9/26/19 1:50" ,  "9/8/19 4:26"   ,
                           "9/10/19 21:22",  "9/11/19 3:32" ,  "9/13/19 6:20"  ,
                           "9/13/19 6:29" ,  "9/20/19 2:40" ,  "9/20/19 3:00"  ,
                           "9/29/19 20:13",  "9/30/19 3:06" ,  "9/11/19 3:46"  ,
                           "9/16/19 1:06" ,  "9/9/19 1:45"  ,  "9/9/19 6:06"   ,
                           "9/10/19 0:18" ,  "9/10/19 2:38" ,  "9/10/19 2:41"  ,
                           "9/10/19 2:43" ,  "9/10/19 6:56" ,  "9/10/19 20:47" ,
                           "9/10/19 23:48",  "9/11/19 3:06" ,  "9/11/19 21:09" ,
                           "9/12/19 0:03" ,  "9/12/19 1:37" ,  "9/12/19 21:54" ,
                           "9/12/19 22:07",  "9/13/19 6:40" ,  "9/13/19 21:35" ,
                           "9/14/19 0:24" ,  "9/14/19 4:24" ,  "9/14/19 4:40"  ,
                           "9/14/19 23:37",  "9/15/19 6:06" ,  "9/15/19 22:35" ,
                           "9/16/19 2:02" ,  "9/16/19 2:04" ,  "9/16/19 2:06"  ,
                           "9/16/19 2:13" ,  "9/16/19 2:14" ,  "9/16/19 3:50"  ,
                           "9/16/19 4:09" ,  "9/16/19 21:31",  "9/16/19 21:34" ,
                           "9/17/19 4:54" ,  "9/18/19 21:58",  "9/19/19 0:24"  ,
                           "9/19/19 1:13" ,  "9/19/19 1:57" ,  "9/19/19 4:49"  ,
                           "9/19/19 21:15",  "9/19/19 21:21",  "9/19/19 21:39" ,
                           "9/19/19 23:57",  "9/20/19 2:18" ,  "9/20/19 21:40" ,
                           "9/20/19 21:45",  "9/20/19 21:55",  "9/20/19 23:55" ,
                           "9/20/19 23:59",  "9/21/19 0:15" ,  "9/21/19 0:17"  ,
                           "9/21/19 0:58" ,  "9/21/19 1:15" ,  "9/22/19 21:22" ,
                           "9/22/19 21:58",  "9/22/19 22:24",  "9/23/19 3:52"  ,
                           "9/24/19 22:08",  "9/24/19 22:21",  "9/24/19 22:22" ,
                           "9/24/19 22:57",  "9/25/19 5:27" ,  "9/25/19 22:00" ,
                           "9/27/19 4:34" ,  "9/27/19 21:20",  "9/27/19 21:55" ,
                           "9/27/19 22:37",  "9/28/19 4:02" ,  "9/28/19 6:00"  ,
                           "9/28/19 21:48",  "9/28/19 22:52",  "9/29/19 5:25"  ,
                           "9/29/19 5:26" ,  "9/29/19 6:12" ,  "9/29/19 8:40" ),
             species = c(  "Red Fox"  ,           "Red Fox",            
                           "Red Fox"   ,          "Red Fox" ,           
                          "Red Fox"     ,        "Red Fox"   ,         
                           "Red Fox"     ,        "Red Fox"   ,         
                           "Red Fox",             "Red Fox"    ,        
                           "Red Fox"  ,           "Red Fox"   ,         
                           "Red Fox" ,            "Red Fox"   ,         
                           "Red Fox"   ,          "Red Fox"   ,         
                           "Red Fox"    ,         "Red Fox"   ,         
                           "Red Fox"     ,        "Red Fox"   ,         
                           "Red Fox"      ,       "Red Fox"   ,         
                           "Red Fox"       ,      "Red Fox"   ,         
                           "Red Fox"        ,     "Red Fox"   ,         
                           "Red Fox"         ,    "Red Fox"   ,         
                           "Red Fox"          ,   "Red Fox"   ,         
                           "Red Fox"           ,  "Red Fox"   ,         
                           "Red Fox",             "Red Fox"   ,         
                           "Red Fox" ,            "Red Fox"   ,         
                           "Red Fox"  ,           "Red Fox"   ,         
                           "Red Fox"   ,          "Red Fox"   ,         
                           "Coyote"     ,         "Coyote"    ,         
                           "Coyote"      ,        "Coyote"    ,         
                           "Coyote"       ,       "Coyote"    ,         
                           "Coyote"        ,      "Coyote"    ,         
                           "Coyote"         ,     "Coyote"    ,         
                           "Coyote"          ,    "Coyote"    ,         
                           "Coyote"           ,   "Coyote"    ,         
                           "Coyote"            ,  "Coyote"    ,         
                           "Coyote" ,             "Coyote"    ,         
                           "Coyote"  ,            "Coyote"     ,        
                           "Coyote"   ,           "Coyote"    ,         
                           "Coyote"    ,          "Coyote"    ,         
                           "Coyote"     ,         "Coyote"    ,         
                           "Coyote"      ,        "Coyote"    ,         
                           "Coyote"       ,       "Coyote"    ,         
                           "Coyote"        ,      "Coyote"    ,         
                           "Coyote"         ,     "Coyote"    ,         
                           "Coyote"          ,    "Coyote"    ,         
                           "Coyote"           ,   "Coyote"    ,         
                           "Coyote",              "Coyote"     ,        
                           "Coyote" ,             "Coyote"      ,       
                           "Mountain Cottontail", "Coyote"       ,      
                           "Coyote",              "Coyote"        ,     
                           "Coyote" ,             "Coyote"         ,    
                           "Coyote"  ,            "Coyote"          ,   
                           "Coyote"   ,           "Coyote"           ,  
                           "Mountain Cottontail", "Coyote"            , 
                           "Mountain Cottontail" ,"Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote",              "Coyote"            , 
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Coyote"             ,
                           "Coyote",              "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Coyote"             ,
                           "Coyote"             , "Coyote"             ,
                          "Coyote"              ,"Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Mountain Cottontail",
                           "Coyote"             , "Coyote"             ,
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Coyote"             ,
                           "Coyote"             , "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Coyote"             ,
                           "Mountain Cottontail", "Coyote"             ,
                           "Mountain Cottontail", "Coyote"             ,
                          "Coyote"              ,"Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Coyote"             , "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                          "Mountain Cottontail" ,"Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Mountain Cottontail",
                          "Mountain Cottontail" ,"Mountain Cottontail",
                           "Mountain Cottontail", "Coyote"             ,
                           "Mountain Cottontail", "Mountain Cottontail",
                           "Mountain Cottontail", "Coyote"             ,
                          "Mountain Cottontail" ,"Mountain Cottontail",
                           "Mountain Cottontail", "Coyote"             ,
                          "Red Fox"             ,"Red Fox"            ,
                          "Red Fox"             ,"Red Fox"            ,
                        "Red Fox"             ,"Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                        "Red Fox"             ,"Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                        "Red Fox"             ,"Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                         "Red Fox"            , "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"           , 
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Coyote"             ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                          "Red Fox"           ,  "Red Fox"            ,
                           "Red Fox"          ,   "Red Fox"            ,
                           "Red Fox" ))
  
  sample.data$obs.time <- strptime(sample.data$obs.time,
                                   format = "%m/%d/%y %H:%M")
  
  UtahData <<- sample.data
}#end data function
