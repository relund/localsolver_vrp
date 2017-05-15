## Functions for simulating data



#' Get instance data for a CVRP problem with one depot and homogeneous trucks.
#'
#' @param trucks Maximum number of trucks.
#' @param capacity Truck capacity
#' @param customers Number of customers.
#' @param demands Min and max demand from customers (chosen uniform).
#' @param distances Min and max travel time / distance / cost to customers (chosen uniform).
#' @param handling Handling time at customer.
#' @param seed Seed for random number generation.
#'
#' @return A list with data. In the distance matrix the first entry is the depot. In the demand
#' vector the first entry is the first customer.
simulateCVRP<-function(trucks = 5, capacity = 10000, customers = 10, demands = c(1000,4000),
                       distances = c(1,300), handling = 30, seed = NULL) {
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
   handlingTime <- as.integer(handling)   # handling time at node in minutes
   return(list(nbTrucks=nbTrucks, truckCapacity=truckCapacity, nbNodes=nbNodes,
               nbCustomers=nbCustomers, demands=demands, distanceMatrix=distanceMatrix,
               handlingTime=handlingTime))
}
