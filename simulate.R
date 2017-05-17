## Functions for simulating data



#' Get instance data for a CVRP problem with one depot and homogeneous trucks.
#'
#' @param trucks Maximum number of trucks.
#' @param capacity Truck capacity
#' @param customers Number of customers.
#' @param demands Min and max demand from customers (chosen uniform).
#' @param distances Min and max travel time / distance / cost to customers (chosen uniform).
#' @param seed Seed for random number generation.
#'
#' @return A list with data. In the distance matrix the first entry is the depot. In the demand
#' vector the first entry is the first customer.
simulateCVRP<-function(trucks = 5, capacity = 10000, customers = 10, demands = c(1000,4000),
                       distances = c(1,300), seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   nbTrucks <- as.integer(trucks)
   truckCapacity <- as.integer(capacity)
   nbNodes <- as.integer(customers + 1)
   nbCustomers <- as.integer(customers)
   demands <- as.integer(  # not including the depot
      sample(demands[1]:demands[2], customers, replace = TRUE)
   )

   # distanceMatrix[i,j] = distance between customer i and j.
   distanceMatrix<-matrix(sample(distances[1]:distances[2], nbNodes^2, replace = TRUE), nbNodes)
   distanceMatrix[lower.tri(distanceMatrix)] <- t(distanceMatrix)[lower.tri(distanceMatrix)]
   diag(distanceMatrix) <- 0
   storage.mode(distanceMatrix) <- "integer"
   return(list(nbTrucks=nbTrucks, truckCapacity=truckCapacity, nbNodes=nbNodes,
               nbCustomers=nbCustomers, demands=demands, distanceMatrix=distanceMatrix))
}



#' Minutes since midnigh
toMin<-function(hour, min=0) {60*hour + min}



#' Get instance data for a VRP problem with drivers and homogeneous trucks and working times.
#'
#' @param maxRoutes Maximum number of routes per driver.
#' @param capacity Truck capacity
#' @param customers Number of customers.
#' @param demands Min and max demand from customers (chosen uniform inbetween).
#' @param distances Min and distance / cost to customers (chosen uniform inbetween).
#' @param handling Handling time at customer.
#' @param seed Seed for random number generation.
#' @param drivers Number of drivers.
#' @param travelTimes Min and max travel time to customers (chosen uniform inbetween).
#' @param workday Earliest start and latest finish for drivers.
#'
#' @return A list with data. In the distance matrix the first entry is the depot. In the demand
#' vector the first entry is the first customer.
simulateDriver<-function(maxRoutes = 5, capacity = 10000, customers = 10, demands = c(1000,4000),
                       distances = c(1,300), handling = 30, drivers = 2, travelTimes = c(5,30),
                       workday = c(toMin(6),toMin(20)), seed = NULL)
{
   dat <- simulateCVRP(trucks = maxRoutes, capacity = capacity, customers = customers,
                       demands = demands, distances = distances, seed = seed)
   dat$maxRoutes <- as.integer(maxRoutes)
   dat$nbTrucks <- NULL
   dat$handlingTime <- as.integer(handling)   # handling time at node in minutes
   dat$drivers <- as.integer(drivers)
   dat$depotEarliestStart <- as.integer(workday[1])  # earliest and latest depature time
   dat$depotLatestArrival <- as.integer(workday[2])

   travelTimeMatrix<-matrix(sample(travelTimes[1]:travelTimes[2], dat$nbNodes^2, replace = TRUE), dat$nbNodes)
   travelTimeMatrix[lower.tri(travelTimeMatrix)] <- t(travelTimeMatrix)[lower.tri(travelTimeMatrix)]
   diag(travelTimeMatrix) <- 0
   storage.mode(travelTimeMatrix) <- "integer"
   dat$travelTimeMatrix <- travelTimeMatrix
   return(dat)
}




#' Get instance data for a VRPTW problem with drivers and homogeneous trucks and working times.
#'
#' @param maxRoutes Maximum number of routes per driver.
#' @param capacity Truck capacity
#' @param customers Number of customers.
#' @param demands Min and max demand from customers (chosen uniform inbetween).
#' @param distances Min and distance / cost to customers (chosen uniform inbetween).
#' @param handling Handling time at customer.
#' @param seed Seed for random number generation.
#' @param drivers Number of drivers.
#' @param travelTimes Min and max travel time to customers (chosen uniform inbetween).
#' @param workday Earliest start and latest finish for drivers.
#' @days Number of days considered to do the routes
#' @param timeWindowSampleHours Start hours for time windows.
#' @param timeWindowSampleGapHours Gap in hours added to time window.
#'
#' @return A list with data. In the distance matrix the first entry is the depot. In the demand
#' vector the first entry is the first customer.
simulateVRPTW<-function(maxRoutes = 5, capacity = 10000, customers = 10, demands = c(1000,4000),
                         distances = c(1,300), handling = 30, drivers = 2, travelTimes = c(5,30),
                         workday = c(toMin(6),toMin(20)), seed = NULL, days = 1,
                        timeWindowSampleHours = c(5,10,15), timeWindowSampleGapHours = c(1,2))
{
   dat <- simulateDriver(maxRoutes =maxRoutes, capacity = capacity, customers = customers,
                       demands = demands, distances = distances, seed = seed, handling = handling,
                       drivers = drivers, travelTimes = travelTimes, workday = workday)

   dat$days <- as.integer(days)

   validDays<-matrix(NA, nrow = customers, ncol = days)
   for (i in 1:customers) {
      ctr<-sample(1:days, 1)
      tmp<-c(rep(0,days-ctr), rep(1,ctr), rep(0,days-ctr))
      ctr<-sample(1:(days-ctr+1), 1)
      validDays[i,]<-tmp<-tmp[ctr:(ctr+days-1)]
   }
   storage.mode(validDays) <- "integer"
   dat$validDays <- validDays

   timeWindowsLB <- sample(timeWindowSampleHours, dat$nbNodes, replace = TRUE)
   timeWindowsUB <- timeWindowsLB + sample(timeWindowSampleGapHours, dat$nbNodes, replace = TRUE)
   timeWindowsLB <- as.integer(toMin(timeWindowsLB))
   timeWindowsUB <- as.integer(toMin(timeWindowsUB))
   timeWindowsLB[1] <- dat$depotEarliestStart
   timeWindowsUB[1] <- dat$depotLatestArrival
   dat$timeWindowsLB <- timeWindowsLB
   dat$timeWindowsUB <- timeWindowsUB

   return(dat)
}
