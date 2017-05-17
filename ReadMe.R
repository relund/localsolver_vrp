## ----cvrpData------------------------------------------------------------
source("simulate.R")
data <- simulateCVRP(customers = 10, seed = 578)
str(data)

## ----cvrpModel, warning=FALSE--------------------------------------------
library(localsolver)
model <- "function model(){
   routeSequences[k in 0..nbTrucks-1] <- list(nbCustomers);      // decision variable
   constraint partition[k in 0..nbTrucks-1](routeSequences[k]);  // visited by exaclty one route

   trucksUsed[k in 0..nbTrucks-1] <- count(routeSequences[k]) > 0;  
   nbTrucksUsed <- sum[k in 0..nbTrucks-1](trucksUsed[k]);

   for [k in 0..nbTrucks-1] {
      local sequence <- routeSequences[k];
      local c <- count(sequence);

      // The quantity needed in each route must not exceed the truck capacity
      routeQuantity[k] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k] <= truckCapacity;

      // Distance travelled by truck k
      routeDistances[k] <- sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]) +
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0);
   }

   // Total distance travelled
   totalDistance <- sum[k in 0..nbTrucks-1](routeDistances[k]);

   // Objective: minimize the number of trucks used, then minimize the distance travelled
   minimize nbTrucksUsed;   
   minimize totalDistance;
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(20,20), indexFromZero=TRUE, lsNbThreads = 4)

## ----cvrpSolve-----------------------------------------------------------
lsp <- set.temp.dir(lsp, path=getwd())    # add temporay files to this folder
lsp <- add.output.expr(lsp, "routeSequences", c(data$nbTrucks, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$nbTrucks)
lsp <- add.output.expr(lsp, "routeQuantity", data$nbTrucks)
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
sol<-ls.solve(lsp, data)

## ----cvrpCheck-----------------------------------------------------------
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0

## ----cvrpPostprocess-----------------------------------------------------
printSolution<-function(sol) {
   cat("Number of routes:", sol$nbTrucksUsed, "\n")
   cat("Total distance:", sol$totalDistance, "\n")
   cat("Routes:\n")
   for (r in 1:length(sol$routeDistances)) {
      if (max(sol$routeSequences[r,])>=0) {
         cat("Customers:",
             paste0(sol$routeSequences[r,which(sol$routeSequences[r,]!=-1)], collapse = "-"),
             "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r],
             "\n")
      }
   }
}
printSolution(sol)

## ----oneDriverData-------------------------------------------------------
data <- simulateDriver(customers = 10, seed = 578, drivers = 1, maxRoutes = 5)
str(data)

## ----oneDriverModel1-----------------------------------------------------
model <- "function model(){
   routeSequences[k in 0..maxRoutes-1] <- list(nbCustomers);  
   routeStops[k in 0..maxRoutes-1][0] = 0;
   routeStops[k in 0..maxRoutes-1][nbNodes] = 0;
   departureTimes[k in 0..maxRoutes-1][i in 0..nbNodes] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..maxRoutes-1][i in 0..nbNodes] <- int(0,depotLatestArrival);

   constraint partition[k in 0..maxRoutes-1](routeSequences[k]);  

   for [k in 0..maxRoutes-1] {
      local sequence <- routeSequences[k];
      local c <- count(sequence);
      routesUsed[k] <- c > 0;

      routeQuantity[k] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k] <= truckCapacity;

      routeDistances[k] <- sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]) +
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0);

      if (k==0) constraint departureTimes[k][0] >= depotEarliestStart;  // first route 
      for [i in 1..nbNodes] {
         routeStops[k][i] <- routeSequences[k][i-1] + 1;
         arrivalTimes[k][i] <- departureTimes[k][i-1] + travelTimeMatrix[routeStops[k][i-1]][routeStops[k][i]];
         departureTimes[k][i] <- (routeStops[k][i] > 0  ? arrivalTimes[k][i] + handlingTime : 0);
      }
      routeArrivalDepot[k] <- departureTimes[k][0] + (c > 0 ? (travelTimeMatrix[0][sequence[0]+1] +
         travelTimeMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => travelTimeMatrix[sequence[i-1]+1][sequence[i]+1]) + c*handlingTime;
      
      if (k>0) {
         constraint departureTimes[k][0] >= routeArrivalDepot[k-1];
         constraint routeArrivalDepot[k] <= depotLatestArrival;
      } else {
         constraint routeArrivalDepot[k] <= depotLatestArrival;
      }
   }

   totalRoutesUsed <- sum[k in 0..maxRoutes-1](routesUsed[k]);
   totalDistance <- sum[k in 0..maxRoutes-1](routeDistances[k]);
   latestArrivalTime <- max[k in 0..maxRoutes-1](routeArrivalDepot[k]);

   minimize totalDistance;
   minimize latestArrivalTime;
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=20, indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$maxRoutes)
lsp <- add.output.expr(lsp, "routeQuantity", data$maxRoutes)
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", data$maxRoutes)
sol<-ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0

## ----oneDriverPostprocess------------------------------------------------
printSolution<-function(sol) {
   cat("Number of routes:", sol$totalRoutesUsed, "\n")
   cat("Total distance:", sol$totalDistance, "\n")
   cat("Latest depot arrival (min from midnight): ", sol$latestArrivalTime, 
       " (constraint <=", data$depotLatestArrival, ")\n", sep = "")
   cat("Routes:\n")
   for (r in 1:length(sol$routeArrivalDepot)) {
      if (max(sol$routeStops[r,])>0) {
         cat("Customers:",
             paste0(sol$routeStops[r,1:(which(sol$routeStops[r,]==0)[2])], collapse = "-"),
             "departure:", sol$departureTimes[r,1], 
             "finished:", sol$arrivalTimes[r, which(sol$routeStops[r,]==0)[2]],
             "distance:", sol$routeDistances[r], 
             "quantity:", sol$routeQuantity[r], "\n")
         cat("  Details:", sol$routeStops[r,1], paste0("(d:", sol$departureTimes[r,1], ") ->"))
         for (i in 2:(which(sol$routeStops[r,]==0)[2]-1)) {
            cat("\n           ", sol$routeStops[r,i], 
                " (tt:", data$travelTimeMatrix[sol$routeStops[r,i-1] + 1, sol$routeStops[r,i]+1],
                " a:", sol$arrivalTimes[r,i], " d:", sol$departureTimes[r,i], ") ", "->", sep="")
         }
         cat("\n           ", sol$routeStops[r,which(sol$routeStops[r,]==0)[2]],
             " (tt:", data$travelTimeMatrix[sol$routeStops[r,i-1]+1,sol$routeStops[r,i]+1],
             " a:", sol$arrivalTimes[r,which(sol$routeStops[r,]==0)[2]], ")\n", sep="")
      }
   }
}
printSolution(sol)

## ----oneDriverModel2-----------------------------------------------------
data <- simulateDriver(customers = 20, seed = 578, drivers = 1, maxRoutes = 12)
lsp <- clear.output.exprs(lsp)
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$maxRoutes)
lsp <- add.output.expr(lsp, "routeQuantity", data$maxRoutes)
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$nbNodes+1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", data$maxRoutes)
sol<-ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0
printSolution(sol)

## ----driversModel--------------------------------------------------------
data$drivers <- as.integer(3)
model <- "function model(){
   routeSequences[k in 0..maxRoutes-1][d in 0..drivers-1] <- list(nbCustomers);  
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][0] = 0;
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][nbNodes] = 0;
   departureTimes[k in 0..maxRoutes-1][d in 0..drivers-1][i in 0..nbNodes] <- int(-1,depotLatestArrival);
   arrivalTimes[k in 0..maxRoutes-1][d in 0..drivers-1][i in 0..nbNodes] <- int(-1,depotLatestArrival);

   constraint partition[k in 0..maxRoutes-1][d in 0..drivers-1](routeSequences[k][d]);

   for [k in 0..maxRoutes-1][d in 0..drivers-1] {
      local sequence <- routeSequences[k][d];
      local c <- count(sequence);
      routesUsed[k][d] <- c > 0;

      routeQuantity[k][d] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d] <= truckCapacity;

      routeDistances[k][d] <- (c > 0 ? (distanceMatrix[0][sequence[0]+1] + 
         distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][0] >= depotEarliestStart;
      for [i in 1..nbNodes] {
         routeStops[k][d][i] <- sequence[i-1] + 1;
         arrivalTimes[k][d][i] <- departureTimes[k][d][i-1] + 
         travelTimeMatrix[routeStops[k][d][i-1]][routeStops[k][d][i]];
         departureTimes[k][d][i] <- (routeStops[k][d][i] > 0  ? arrivalTimes[k][d][i] + handlingTime : -1);
      }
      routeArrivalDepot[k][d] <- departureTimes[k][d][0] + (c > 0 ? (travelTimeMatrix[0][sequence[0]+1] + travelTimeMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => travelTimeMatrix[sequence[i-1]+1][sequence[i]+1]) + c*handlingTime;

      if (k>0) {
         constraint departureTimes[k][d][0] >= routeArrivalDepot[k-1][d];
         constraint routeArrivalDepot[k][d] <= depotLatestArrival;
      } else {
         constraint routeArrivalDepot[k][d] <= depotLatestArrival;
      }
   }

   totalRoutesUsed <- sum[k in 0..maxRoutes-1][d in 0..drivers-1] (routesUsed[k][d]);
   totalDistance <- sum[k in 0..maxRoutes-1][d in 0..drivers-1](routeDistances[k][d]);
   latestArrivalTime <- max[k in 0..maxRoutes-1][d in 0..drivers-1] (routeArrivalDepot[k][d]);

   minimize totalDistance;
   minimize latestArrivalTime;
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(20,20), indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$drivers, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", c(data$maxRoutes, data$drivers))
lsp <- add.output.expr(lsp, "routeQuantity", c(data$maxRoutes, data$drivers))
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$drivers, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$drivers, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$drivers, data$nbNodes+1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", c(data$maxRoutes, data$drivers))
sol<-ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0

## ----driversPostProcess--------------------------------------------------
printSolution<-function(sol) {
   cat("Number of routes:", sol$totalRoutesUsed, "\n")
   cat("Total distance:", sol$totalDistance, "\n")
   cat("Latest depot arrival (min from midnight):", sol$latestArrivalTime, "\n")
   cat("Routes:\n")
   for (d in 1:length(sol$routeArrivalDepot[1,])) {
      for (r in 1:length(sol$routeArrivalDepot[,1])) {
         if (max(sol$routeStops[r,d,])>0) {
            cat("Driver:", d, "Customers:",
                paste0(sol$routeStops[r,d,1:(which(sol$routeStops[r,d,]==0)[2])], collapse = "-"),
                "departure:", sol$departureTimes[r,d,1], "finished:", sol$arrivalTimes[r,d, which(sol$routeStops[r,d,]==0)[2]],
                "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r],
                "\n")
            cat("  Details:", sol$routeStops[r,d,1], paste0("(d:", sol$departureTimes[r,d,1], ") ->"))
            for (i in 2:(which(sol$routeStops[r,d,]==0)[2]-1)) {
               cat("\n           ", sol$routeStops[r,d,i], " (tt:", 
                   data$travelTimeMatrix[sol$routeStops[r,d,i-1]+1,sol$routeStops[r,d,i]+1],
                   " a:", sol$arrivalTimes[r,d,i], " d:", sol$departureTimes[r,d,i], ") ", "->", sep="")
            }
            cat("\n           ", sol$routeStops[r,d,which(sol$routeStops[r,d,]==0)[2]],
                " (tt:", data$travelTimeMatrix[sol$routeStops[r,d,i-1]+1,sol$routeStops[r,d,i]+1],
                " a:", sol$arrivalTimes[r,d,which(sol$routeStops[r,d,]==0)[2]], ")\n", sep="")
         }
      }
   }
}
printSolution(sol)

## ----daysData, cache=FALSE-----------------------------------------------
data<-simulateVRPTW(drivers = 2, maxRoutes = 5, customers = 30, days = 7, seed = 789)
str(data)
data$validDays

## ----daysModel-----------------------------------------------------------
model <- "
function model(){
   routeSequences[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] <- list(nbCustomers); 
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][0] = 0;
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][nbNodes] = 0;
   departureTimes[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][0] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,depotLatestArrival);

   constraint partition[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeSequences[k][d][t]);

   for [k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] {
      local sequence <- routeSequences[k][d][t];
      local c <- count(sequence);
      routesUsed[k][d][t] <- c > 0;

      routeQuantity[k][d][t] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d][t] <= truckCapacity;

      routeDistances[k][d][t] <-
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][t][0] >= depotEarliestStart;
      for [i in 1..nbNodes] {
         routeStops[k][d][t][i] <- sequence[i-1] + 1;
         arrivalTimes[k][d][t][i] <- (routeStops[k][d][t][i-1] > 0 || (i==1 && routeStops[k][d][t][i]>0) ?
            departureTimes[k][d][t][i-1] +
            travelTimeMatrix[routeStops[k][d][t][i-1]][routeStops[k][d][t][i]] : 0);
         departureTimes[k][d][t][i] <-
            (routeStops[k][d][t][i] > 0  ? arrivalTimes[k][d][t][i] + handlingTime : 0);
      }
      routeArrivalDepot[k][d][t] <- max[i in 1..nbCustomers](arrivalTimes[k][d][t][i]);

      if (k>0) {
         constraint departureTimes[k][d][t][0] >= routeArrivalDepot[k-1][d][t];
         constraint routeArrivalDepot[k][d][t] <= depotLatestArrival;
      } else {
         constraint routeArrivalDepot[k][d][t] <= depotLatestArrival;
      }
   }

   for [t in 0..days-1][i in 0..nbCustomers-1] {
      visitAtDay[t][i] <- sum[k in 0..maxRoutes-1][d in 0..drivers-1](indexOf(routeSequences[k][d][t],i) > -1);
      constraint visitAtDay[t][i] <= validDays[i][t];
   }

   totalRoutesUsed <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] (routesUsed[k][d][t]);
   totalDistance <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeDistances[k][d][t]);
   latestArrival <- max[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeArrivalDepot[k][d][t]);

   minimize totalDistance;
   minimize latestArrival;

   // modify results so can output to R (cannot output array with more than 3 index)
   for [k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] {
         local r = k + d*maxRoutes + t*maxRoutes*drivers;
         routes[r][i in 0..nbNodes] <- routeStops[k][d][t][i];
         arrivalT[r][i in 0..nbNodes] <- arrivalTimes[k][d][t][i];
         routeDepotArrival[r] <- routeArrivalDepot[k][d][t];
         departureT[r][i in 0..nbNodes] <- departureTimes[k][d][t][i];
         routeDist[r] <- routeDistances[k][d][t];
         routeQuant[r] <- routeQuantity[k][d][t];
         driver[r] <- d+1;
         day[r] <- t+1;
   }
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(10,10), indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$maxRoutes*data$days*data$drivers
lsp <- add.output.expr(lsp, "routes", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "driver", maxR)
lsp <- add.output.expr(lsp, "day", maxR)
lsp <- add.output.expr(lsp, "routeDist", maxR)
lsp <- add.output.expr(lsp, "routeDepotArrival", maxR)
lsp <- add.output.expr(lsp, "routeQuant", maxR)
sol<-ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0

## ----daysPostProcess-----------------------------------------------------
printSolution<-function(sol) {
   cat("Number of routes:", sol$totalRoutesUsed, "\n")
   cat("Total distance:", sol$totalDistance, "\n")
   cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
   cat("Routes:\n")
   for (r in 1:length(sol$driver)) {
      if (max(sol$routes[r,])>0) {
         cat("Day:", sol$day[r],"driver:", sol$driver[r], "customers:",
             paste0(sol$routes[r,1:(which(sol$routes[r,]==0)[2])], collapse = "-"),
             "departure:", sol$departureT[r,1], "finished:", sol$arrivalT[r, which(sol$routes[r,]==0)[2]],
             "distance:", sol$routeDist[r], "quantity:", sol$routeQuant[r],
             "\n")
         cat("  Details:", sol$routes[r,1], paste0("(d:", sol$departureT[r,1], ") ->"))
         for (i in 2:(which(sol$routes[r,]==0)[2]-1)) {
            cat("\n           ", sol$routes[r,i], " (tt:", 
                data$travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1], 
                " a:", sol$arrivalT[r,i], " d:", sol$departureT[r,i], ") ", "->", sep="")
         }
         cat("\n           ", sol$routes[r,which(sol$routes[r,]==0)[2]],
             " (tt:", data$travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1],
             " a:", sol$arrivalT[r,which(sol$routes[r,]==0)[2]], ")\n", sep="")
      }
   }
}
printSolution(sol)

## ----twModel-------------------------------------------------------------
model <- "
function model(){
   routeSequences[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] <- list(nbCustomers); 
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][0] = 0;
   routeStops[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][nbNodes] = 0;
   departureTimes[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][0] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,depotLatestArrival);
   waitingTimes[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,600);
   constraint partition[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeSequences[k][d][t]);

   for [k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] {
      local sequence <- routeSequences[k][d][t];
      local c <- count(sequence);
      routesUsed[k][d][t] <- c > 0;

      routeQuantity[k][d][t] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d][t] <= truckCapacity;

      routeDistances[k][d][t] <-
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][t][0] >= depotEarliestStart;
      for [i in 1..nbNodes] {
         routeStops[k][d][t][i] <- sequence[i-1] + 1;
         arrivalTimes[k][d][t][i] <- (routeStops[k][d][t][i-1] > 0 || (i==1 && routeStops[k][d][t][i]>0) ?
            departureTimes[k][d][t][i-1] +
            travelTimeMatrix[routeStops[k][d][t][i-1]][routeStops[k][d][t][i]] +
            waitingTimes[k][d][t][i] : 0);
         // If hard time window constraints
         //constraint arrivalTimes[k][d][t][i] >= (routeStops[k][d][t][i] > 0) * timeWindowsLB[routeStops[k][d][t][i]];
         //constraint arrivalTimes[k][d][t][i] <= timeWindowsUB[routeStops[k][d][t][i]];
         departureTimes[k][d][t][i] <-
            (routeStops[k][d][t][i] > 0  ? arrivalTimes[k][d][t][i] + handlingTime : 0);
      }
      routeArrivalDepot[k][d][t] <- max[i in 1..nbCustomers](arrivalTimes[k][d][t][i]);

      if (k>0) {
         constraint departureTimes[k][d][t][0] >= routeArrivalDepot[k-1][d][t];
         constraint routeArrivalDepot[k][d][t] <= depotLatestArrival;
      } else {
         constraint routeArrivalDepot[k][d][t] <= depotLatestArrival;
      }
   }

   for [t in 0..days-1][i in 0..nbCustomers-1] {
      visitAtDay[t][i] <- sum[k in 0..maxRoutes-1][d in 0..drivers-1](indexOf(routeSequences[k][d][t],i) > -1);
      constraint visitAtDay[t][i] <= validDays[i][t];
   }

   totalRoutesUsed <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] (routesUsed[k][d][t]);
   totalDistance <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeDistances[k][d][t]);
   totalWaiting <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes]
      (waitingTimes[k][d][t][i]);
   totalViolateTWCtr <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 1..nbCustomers](
      routeStops[k][d][t][i] > 0 && (arrivalTimes[k][d][t][i] < timeWindowsLB[routeStops[k][d][t][i]]  ||
      arrivalTimes[k][d][t][i] > timeWindowsUB[routeStops[k][d][t][i]])
   );
   totalViolateTWMin <- sum[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1][i in 1..nbCustomers](
      (routeStops[k][d][t][i] > 0 && arrivalTimes[k][d][t][i] < timeWindowsLB[routeStops[k][d][t][i]]) *
      (timeWindowsLB[routeStops[k][d][t][i]] - arrivalTimes[k][d][t][i]) +
      (routeStops[k][d][t][i] > 0 && arrivalTimes[k][d][t][i] > timeWindowsUB[routeStops[k][d][t][i]]) *
      (arrivalTimes[k][d][t][i] - timeWindowsUB[routeStops[k][d][t][i]])
   );
   latestArrival <- max[k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1](routeArrivalDepot[k][d][t]);

   minimize totalViolateTWMin;   // total waiting time in minutes
   minimize totalViolateTWCtr;   // times we wait
   minimize totalRoutesUsed;
   minimize totalDistance;
   minimize totalWaiting;
   minimize latestArrival;

   // modify results so can output to R (cannot output array with more than 3 index)
   for [k in 0..maxRoutes-1][d in 0..drivers-1][t in 0..days-1] {
         local r = k + d*maxRoutes + t*maxRoutes*drivers;
         routes[r][i in 0..nbNodes] <- routeStops[k][d][t][i];
         arrivalT[r][i in 0..nbNodes] <- arrivalTimes[k][d][t][i];
         routeDepotArrival[r] <- routeArrivalDepot[k][d][t];
         departureT[r][i in 0..nbNodes] <- departureTimes[k][d][t][i];
         waitingT[r][i in 0..nbNodes] <- waitingTimes[k][d][t][i];
         routeDist[r] <- routeDistances[k][d][t];
         routeQuant[r] <- routeQuantity[k][d][t];
         driver[r] <- d+1;
         day[r] <- t+1;
   }
}"
rm(lsp)
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(20,20,20,20,20,20), indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "totalWaiting")
lsp <- add.output.expr(lsp, "totalViolateTWCtr")
lsp <- add.output.expr(lsp, "totalViolateTWMin")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$maxRoutes*data$days*data$drivers
lsp <- add.output.expr(lsp, "routes", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "waitingT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "driver", maxR)
lsp <- add.output.expr(lsp, "day", maxR)
lsp <- add.output.expr(lsp, "routeDist", maxR)
lsp <- add.output.expr(lsp, "routeDepotArrival", maxR)
lsp <- add.output.expr(lsp, "routeQuant", maxR)
sol<-ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE))>0

## ----twPostProcess-------------------------------------------------------
printSolution<-function(sol) {
   cat("Number of routes:", sol$totalRoutesUsed, "\n")
   cat("Total distance:", sol$totalDistance, "\n")
   cat("Total waiting time:", sol$totalWaiting, "\n")
   cat("Total time window violations:", sol$totalViolateTWCtr, "\n")
   cat("Total time window violations (min):", sol$totalViolateTWMin, "\n")
   cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
   cat("Routes:\n")
   for (r in 1:length(sol$driver)) {
      if (max(sol$routes[r,])>0) {
         cat("Day:", sol$day[r],"driver:", sol$driver[r], "customers:",
             paste0(sol$routes[r,1:(which(sol$routes[r,]==0)[2])], collapse = "-"),
             "departure:", sol$departureT[r,1], "finished:", sol$arrivalT[r, which(sol$routes[r,]==0)[2]],
             "distance:", sol$routeDist[r], "quantity:", sol$routeQuant[r], "\n")
         cat("  Details:", sol$routes[r,1], paste0("(d:", sol$departureT[r,1], ") ->"))
         for (i in 2:(which(sol$routes[r,]==0)[2]-1)) {
            cat("\n           ", sol$routes[r,i], " (tt:", data$travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1],
                " w:", sol$waitingT[r,i], " a:", sol$arrivalT[r,i], " d:", sol$departureT[r,i], ") ",
                " tw:[", data$timeWindowsLB[sol$routes[r,i]+1], ",", data$timeWindowsUB[sol$routes[r,i]+1], "] ->", sep="")
         }
         cat("\n           ", sol$routes[r,which(sol$routes[r,]==0)[2]],
             " (tt:", data$travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1],
             " a:", sol$arrivalT[r,which(sol$routes[r,]==0)[2]], ")\n", sep="")
      }
   }
}
printSolution(sol)

