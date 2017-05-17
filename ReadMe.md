
Using LocalSolver to solve VRP problems
=======================================

The last days I have been playing a bit with [LocalSolver](http://www.localsolver.com/) and tried to model various vehicle routing problems (VRP). LocalSolver is a heuristic solver based on a heuristic search approach combining different optimization techniques. It includes a math modeling language (LSP) and there is a package so that it can be used from [R](https://www.r-project.org/). You need to download and install LocalSolver (I used v7.0) and get an academic license from their webpage.

Below I will formulate models for different VRPs and solve them using LocalSolver. The purpose of these tests is not to compare them with an optimal solution to the problems, but to examine the modelling capabilities of LSP. All the code used can be found in files `simulate.R` and `ReadMe.R`.

Capacitated VRP (CVRP)
----------------------

Let us start by consider the classical CVRP. First, we simulate some data:

``` r
source("simulate.R")
data <- simulateCVRP(customers = 10, seed = 578)
str(data)
#> List of 6
#>  $ nbTrucks      : int 5
#>  $ truckCapacity : int 10000
#>  $ nbNodes       : int 11
#>  $ nbCustomers   : int 10
#>  $ demands       : int [1:10] 2514 2585 2938 3935 2649 1429 3766 3706 2778 2547
#>  $ distanceMatrix: int [1:11, 1:11] 0 144 5 31 14 97 167 161 68 289 ...
```

Note that all data are integers since I have chosen to use integers in LocalSolver (input data must be the same data type). Next, we formulate the model using [LSP](http://www.localsolver.com/documentation/lspreference/index.html). Here I used the example on LocalSolver's webpage as a starting point. The learning curve for LSP is quite steep; however, the [examples](http://www.localsolver.com/documentation/exampletour/index.html) help. I was missing a forum where to ask question though. An important issue is that you must choose where the index in arrays start. I choose C style and let all index start from zero. We need decision variables `routeSequences[k]` containing a list of numbers in the range zero to `customers-1`. That is, if `routeSequences[k]` equals {0,4,6}, then truck mumber `k` visit customer 1, 5 and 7 (remember index start from zero).

``` r
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
lsp <- set.params(lsp, lsTimeLimit = c(20, 20), indexFromZero = TRUE, lsNbThreads = 4)
```

The model is added to LocalSolver using the [localsolver package](https://cran.r-project.org/web/packages/localsolver/index.html). I set the time limit to 20 seconds for each objective. Output returned to R can be added using function `add.output.expr`.

``` r
lsp <- set.temp.dir(lsp, path = getwd())  # add temporay files to this folder
lsp <- add.output.expr(lsp, "routeSequences", c(data$nbTrucks, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$nbTrucks)
lsp <- add.output.expr(lsp, "routeQuantity", data$nbTrucks)
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
sol <- ls.solve(lsp, data)
```

The results are now in the list `sol`. Three temporay files are used by the solver. `input.lsp` contains the model (with input and output) passed to the solver, `output.txt` the solver log and `data.txt` the input data. Unfortuantly, the solver status is not returned! So you may have to check the log to see if a feasible solution was found:

``` r
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] TRUE
```

We postprocess the solution to get a result in a more readable format:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$routeDistances)) {
        if (max(sol$routeSequences[r, ]) >= 0) {
            cat("Customers:", paste0(sol$routeSequences[r, which(sol$routeSequences[r, ] != -1)], collapse = "-"), 
                "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r], "\n")
        }
    }
}
printSolution(sol)
#> Number of routes: 3 
#> Total distance: 1147 
#> Routes:
#> Customers: 4-1-3 distance: 222 quantity: 9169 
#> Customers: 7-6-0 distance: 559 quantity: 9986 
#> Customers: 2-5-8-9 distance: 366 quantity: 9692
```

Extension - Workday constraints
-------------------------------

Assume that only a single driver is avaliabele and that the maximum possible workday is specified. We need workday length, travel times and handling time. Note that we now use `maxRoutes` instead of `trucks` which denote the maximum number of routes per driver.

``` r
data <- simulateDriver(customers = 10, seed = 578, drivers = 1, maxRoutes = 5)
str(data)
#> List of 11
#>  $ truckCapacity     : int 10000
#>  $ nbNodes           : int 11
#>  $ nbCustomers       : int 10
#>  $ demands           : int [1:10] 2514 2585 2938 3935 2649 1429 3766 3706 2778 2547
#>  $ distanceMatrix    : int [1:11, 1:11] 0 144 5 31 14 97 167 161 68 289 ...
#>  $ maxRoutes         : int 5
#>  $ handlingTime      : int 30
#>  $ drivers           : int 1
#>  $ depotEarliestStart: int 360
#>  $ depotLatestArrival: int 1200
#>  $ travelTimeMatrix  : int [1:11, 1:11] 0 14 28 29 6 26 24 29 13 11 ...
```

To handle travel time constraints we add variables `routeStops[k][i]` which is the node visited at stop i (depot = node 0) and similar `arrivalTimes[k][i]` and `departureTimes[k][i]`.

``` r
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
lsp <- set.params(lsp, lsTimeLimit = 20, indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$maxRoutes)
lsp <- add.output.expr(lsp, "routeQuantity", data$maxRoutes)
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", data$maxRoutes)
sol <- ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] TRUE
```

We have minimized the total distance and next the lastest arrival time at the depot and postprocess the solution to get the result in a more readable format:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$totalRoutesUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight): ", sol$latestArrivalTime, " (constraint <=", data$depotLatestArrival, 
        ")\n", sep = "")
    cat("Routes:\n")
    for (r in 1:length(sol$routeArrivalDepot)) {
        if (max(sol$routeStops[r, ]) > 0) {
            cat("Customers:", paste0(sol$routeStops[r, 1:(which(sol$routeStops[r, ] == 0)[2])], collapse = "-"), 
                "departure:", sol$departureTimes[r, 1], "finished:", sol$arrivalTimes[r, which(sol$routeStops[r, 
                  ] == 0)[2]], "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r], 
                "\n")
            cat("  Details:", sol$routeStops[r, 1], paste0("(d:", sol$departureTimes[r, 1], ") ->"))
            for (i in 2:(which(sol$routeStops[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routeStops[r, i], " (tt:", data$travelTimeMatrix[sol$routeStops[r, 
                  i - 1] + 1, sol$routeStops[r, i] + 1], " a:", sol$arrivalTimes[r, i], " d:", sol$departureTimes[r, 
                  i], ") ", "->", sep = "")
            }
            cat("\n           ", sol$routeStops[r, which(sol$routeStops[r, ] == 0)[2]], " (tt:", data$travelTimeMatrix[sol$routeStops[r, 
                i - 1] + 1, sol$routeStops[r, i] + 1], " a:", sol$arrivalTimes[r, which(sol$routeStops[r, 
                ] == 0)[2]], ")\n", sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 4 
#> Total distance: 795 
#> Latest depot arrival (min from midnight): 908 (constraint <=1200)
#> Routes:
#> Customers: 0-8-0 departure: 360 finished: 416 distance: 136 quantity: 3706 
#>   Details: 0 (d:360) ->
#>            8 (tt:13 a:373 d:403) ->
#>            0 (tt:13 a:416)
#> Customers: 0-1-10-7-0 departure: 416 finished: 572 distance: 377 quantity: 8827 
#>   Details: 0 (d:416) ->
#>            1 (tt:14 a:430 d:460) ->
#>            10 (tt:6 a:466 d:496) ->
#>            7 (tt:17 a:513 d:543) ->
#>            0 (tt:17 a:572)
#> Customers: 0-4-2-0 departure: 572 finished: 675 distance: 47 quantity: 6520 
#>   Details: 0 (d:572) ->
#>            4 (tt:6 a:578 d:608) ->
#>            2 (tt:9 a:617 d:647) ->
#>            0 (tt:9 a:675)
#> Customers: 0-5-6-9-3-0 departure: 675 finished: 908 distance: 235 quantity: 9794 
#>   Details: 0 (d:675) ->
#>            5 (tt:26 a:701 d:731) ->
#>            6 (tt:9 a:740 d:770) ->
#>            9 (tt:22 a:792 d:822) ->
#>            3 (tt:27 a:849 d:879) ->
#>            0 (tt:27 a:908)
```

Note currently there is no time added for filling the truck at the depot. Let us try another instance.

``` r
data <- simulateDriver(customers = 20, seed = 578, drivers = 1, maxRoutes = 12)
lsp <- clear.output.exprs(lsp)
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$maxRoutes)
lsp <- add.output.expr(lsp, "routeQuantity", data$maxRoutes)
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", data$maxRoutes)
sol <- ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] FALSE
printSolution(sol)
#> Number of routes: 7 
#> Total distance: 3682 
#> Latest depot arrival (min from midnight): 1200 (constraint <=1200)
#> Routes:
#> Customers: 0-6-14-16-5-0 departure: 360 finished: 514 distance: 708 quantity: 9133 
#>   Details: 0 (d:360) ->
#>            6 (tt:5 a:365 d:395) ->
#>            14 (tt:5 a:400 d:430) ->
#>            16 (tt:11 a:441 d:471) ->
#>            5 (tt:5 a:476 d:506) ->
#>            0 (tt:5 a:514)
#> Customers: 0-3-19-9-0 departure: 514 finished: 644 distance: 463 quantity: 7914 
#>   Details: 0 (d:514) ->
#>            3 (tt:6 a:520 d:550) ->
#>            19 (tt:7 a:557 d:587) ->
#>            9 (tt:10 a:597 d:627) ->
#>            0 (tt:10 a:644)
#> Customers: 0-1-20-11-0 departure: 644 finished: 780 distance: 545 quantity: 9844 
#>   Details: 0 (d:644) ->
#>            1 (tt:11 a:655 d:685) ->
#>            20 (tt:10 a:695 d:725) ->
#>            11 (tt:12 a:737 d:767) ->
#>            0 (tt:12 a:780)
#> Customers: 0-7-17-0 departure: 780 finished: 869 distance: 437 quantity: 7642 
#>   Details: 0 (d:780) ->
#>            7 (tt:15 a:795 d:825) ->
#>            17 (tt:9 a:834 d:864) ->
#>            0 (tt:9 a:869)
#> Customers: 0-4-12-8-0 departure: 869 finished: 987 distance: 542 quantity: 8861 
#>   Details: 0 (d:869) ->
#>            4 (tt:5 a:874 d:904) ->
#>            12 (tt:11 a:915 d:945) ->
#>            8 (tt:6 a:951 d:981) ->
#>            0 (tt:6 a:987)
#> Customers: 0-2-10-18-0 departure: 973 finished: 1123 distance: 536 quantity: 9062 
#>   Details: 0 (d:973) ->
#>            2 (tt:20 a:993 d:1023) ->
#>            10 (tt:9 a:1032 d:1062) ->
#>            18 (tt:7 a:1069 d:1099) ->
#>            0 (tt:7 a:1123)
#> Customers: 0-13-15-0 departure: 1123 finished: 1200 distance: 451 quantity: 5671 
#>   Details: 0 (d:1123) ->
#>            13 (tt:6 a:1129 d:1159) ->
#>            15 (tt:6 a:1165 d:1195) ->
#>            0 (tt:6 a:1200)
```

Note no feasible solution can be found, since the driver can not keep the latest arrival at depot limit. That is, the printed solution is the unfeasible soltion found when the solver stopped. This is not the same as a proof of that no feasible solution exists since it is a heuristic solver we use.

Extension - Multiple drivers
----------------------------

Assume that multiple drivers are allowed and for simplicity assume that they have the same working time. To handle multiple drivers we need to add a driver index to the model.

``` r
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
lsp <- set.params(lsp, lsTimeLimit = c(20, 20), indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$maxRoutes, data$drivers, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", c(data$maxRoutes, data$drivers))
lsp <- add.output.expr(lsp, "routeQuantity", c(data$maxRoutes, data$drivers))
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$maxRoutes, data$drivers, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$maxRoutes, data$drivers, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$maxRoutes, data$drivers, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", c(data$maxRoutes, data$drivers))
sol <- ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] TRUE
```

We consided the same instance as before but now with 3 drivers, minimized the total distance and next the latest arrival time at the depot. The solution are:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$totalRoutesUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrivalTime, "\n")
    cat("Routes:\n")
    for (d in 1:length(sol$routeArrivalDepot[1, ])) {
        for (r in 1:length(sol$routeArrivalDepot[, 1])) {
            if (max(sol$routeStops[r, d, ]) > 0) {
                cat("Driver:", d, "Customers:", paste0(sol$routeStops[r, d, 1:(which(sol$routeStops[r, 
                  d, ] == 0)[2])], collapse = "-"), "departure:", sol$departureTimes[r, d, 1], "finished:", 
                  sol$arrivalTimes[r, d, which(sol$routeStops[r, d, ] == 0)[2]], "distance:", sol$routeDistances[r], 
                  "quantity:", sol$routeQuantity[r], "\n")
                cat("  Details:", sol$routeStops[r, d, 1], paste0("(d:", sol$departureTimes[r, d, 1], 
                  ") ->"))
                for (i in 2:(which(sol$routeStops[r, d, ] == 0)[2] - 1)) {
                  cat("\n           ", sol$routeStops[r, d, i], " (tt:", data$travelTimeMatrix[sol$routeStops[r, 
                    d, i - 1] + 1, sol$routeStops[r, d, i] + 1], " a:", sol$arrivalTimes[r, d, i], " d:", 
                    sol$departureTimes[r, d, i], ") ", "->", sep = "")
                }
                cat("\n           ", sol$routeStops[r, d, which(sol$routeStops[r, d, ] == 0)[2]], " (tt:", 
                  data$travelTimeMatrix[sol$routeStops[r, d, i - 1] + 1, sol$routeStops[r, d, i] + 1], 
                  " a:", sol$arrivalTimes[r, d, which(sol$routeStops[r, d, ] == 0)[2]], ")\n", sep = "")
            }
        }
    }
}
printSolution(sol)
#> Number of routes: 7 
#> Total distance: 1476 
#> Latest depot arrival (min from midnight): 728 
#> Routes:
#> Driver: 1 Customers: 0-8-18-6-0 departure: 364 finished: 501 distance: 153 quantity: 9065 
#>   Details: 0 (d:364) ->
#>            8 (tt:6 a:370 d:400) ->
#>            18 (tt:15 a:415 d:445) ->
#>            6 (tt:21 a:466 d:496) ->
#>            0 (tt:21 a:501)
#> Driver: 1 Customers: 0-12-9-14-13-0 departure: 519 finished: 719 distance: 319 quantity: 8426 
#>   Details: 0 (d:519) ->
#>            12 (tt:15 a:534 d:564) ->
#>            9 (tt:30 a:594 d:624) ->
#>            14 (tt:19 a:643 d:673) ->
#>            13 (tt:10 a:683 d:713) ->
#>            0 (tt:10 a:719)
#> Driver: 2 Customers: 0-3-15-10-0 departure: 379 finished: 549 distance: 0 quantity: 0 
#>   Details: 0 (d:379) ->
#>            3 (tt:6 a:385 d:415) ->
#>            15 (tt:20 a:435 d:465) ->
#>            10 (tt:28 a:493 d:523) ->
#>            0 (tt:28 a:549)
#> Driver: 2 Customers: 0-2-5-4-0 departure: 575 finished: 711 distance: 0 quantity: 0 
#>   Details: 0 (d:575) ->
#>            2 (tt:20 a:595 d:625) ->
#>            5 (tt:5 a:630 d:660) ->
#>            4 (tt:16 a:676 d:706) ->
#>            0 (tt:16 a:711)
#> Driver: 3 Customers: 0-7-0 departure: 360 finished: 420 distance: 0 quantity: 0 
#>   Details: 0 (d:360) ->
#>            7 (tt:15 a:375 d:405) ->
#>            0 (tt:15 a:420)
#> Driver: 3 Customers: 0-19-11-20-0 departure: 420 finished: 577 distance: 0 quantity: 0 
#>   Details: 0 (d:420) ->
#>            19 (tt:20 a:440 d:470) ->
#>            11 (tt:7 a:477 d:507) ->
#>            20 (tt:12 a:519 d:549) ->
#>            0 (tt:12 a:577)
#> Driver: 3 Customers: 0-1-17-16-0 departure: 577 finished: 728 distance: 0 quantity: 0 
#>   Details: 0 (d:577) ->
#>            1 (tt:11 a:588 d:618) ->
#>            17 (tt:14 a:632 d:662) ->
#>            16 (tt:17 a:679 d:709) ->
#>            0 (tt:17 a:728)
```

Note with 3 drivers we can keep the latest arrival at depot limit.

Extension - Delivery over multiple days
---------------------------------------

Assume that the orders can be delivered at specific days (we ignore the time windows for the moment):

``` r
data <- simulateVRPTW(drivers = 2, maxRoutes = 5, customers = 30, days = 7, seed = 789)
str(data)
#> List of 15
#>  $ truckCapacity     : int 10000
#>  $ nbNodes           : int 31
#>  $ nbCustomers       : int 30
#>  $ demands           : int [1:30] 3100 1280 1035 2775 2476 1060 2718 1497 2077 2039 ...
#>  $ distanceMatrix    : int [1:31, 1:31] 0 212 1 214 236 78 226 240 300 86 ...
#>  $ maxRoutes         : int 5
#>  $ handlingTime      : int 30
#>  $ drivers           : int 2
#>  $ depotEarliestStart: int 360
#>  $ depotLatestArrival: int 1200
#>  $ travelTimeMatrix  : int [1:31, 1:31] 0 10 27 19 20 11 9 29 29 11 ...
#>  $ days              : int 7
#>  $ validDays         : int [1:30, 1:7] 0 1 1 1 0 1 1 0 0 0 ...
#>  $ timeWindowsLB     : int [1:31] 360 900 600 900 300 300 900 600 900 600 ...
#>  $ timeWindowsUB     : int [1:31] 1200 1020 660 1020 420 360 1020 720 1020 720 ...
data$validDays
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#>  [1,]    0    1    1    1    1    1    1
#>  [2,]    1    1    1    1    1    1    0
#>  [3,]    1    1    1    1    1    1    1
#>  [4,]    1    1    1    1    1    0    0
#>  [5,]    0    0    0    0    0    1    0
#>  [6,]    1    1    1    1    1    1    1
#>  [7,]    1    1    1    1    1    1    1
#>  [8,]    0    1    1    1    1    1    0
#>  [9,]    0    0    1    1    0    0    0
#> [10,]    0    0    0    0    1    1    0
#> [11,]    1    1    1    1    1    1    1
#> [12,]    0    1    1    1    1    1    1
#> [13,]    1    1    1    1    1    1    0
#> [14,]    1    1    1    1    1    1    1
#> [15,]    0    0    1    1    1    0    0
#> [16,]    0    0    0    1    1    1    1
#> [17,]    0    0    0    0    1    1    0
#> [18,]    0    1    1    1    0    0    0
#> [19,]    0    1    1    0    0    0    0
#> [20,]    0    0    1    1    0    0    0
#> [21,]    1    1    0    0    0    0    0
#> [22,]    0    0    0    1    0    0    0
#> [23,]    1    0    0    0    0    0    0
#> [24,]    1    1    1    1    1    1    1
#> [25,]    1    1    1    1    1    1    1
#> [26,]    0    0    1    1    1    1    1
#> [27,]    1    1    1    1    1    1    1
#> [28,]    0    1    1    1    1    1    0
#> [29,]    0    0    0    1    1    1    1
#> [30,]    0    1    1    1    0    0    0
```

Valid days for an order corresponds to a one in the matrix `validDays`. For instance, customer 29 can be visited in days 4, 5, 6, 7. We add a constraint that customer i must be visited at valid days.

``` r
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
lsp <- set.params(lsp, lsTimeLimit = c(10, 10), indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$maxRoutes * data$days * data$drivers
lsp <- add.output.expr(lsp, "routes", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalT", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureT", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "driver", maxR)
lsp <- add.output.expr(lsp, "day", maxR)
lsp <- add.output.expr(lsp, "routeDist", maxR)
lsp <- add.output.expr(lsp, "routeDepotArrival", maxR)
lsp <- add.output.expr(lsp, "routeQuant", maxR)
sol <- ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] TRUE
```

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$totalRoutesUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$driver)) {
        if (max(sol$routes[r, ]) > 0) {
            cat("Day:", sol$day[r], "driver:", sol$driver[r], "customers:", paste0(sol$routes[r, 1:(which(sol$routes[r, 
                ] == 0)[2])], collapse = "-"), "departure:", sol$departureT[r, 1], "finished:", sol$arrivalT[r, 
                which(sol$routes[r, ] == 0)[2]], "distance:", sol$routeDist[r], "quantity:", sol$routeQuant[r], 
                "\n")
            cat("  Details:", sol$routes[r, 1], paste0("(d:", sol$departureT[r, 1], ") ->"))
            for (i in 2:(which(sol$routes[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routes[r, i], " (tt:", data$travelTimeMatrix[sol$routes[r, i - 
                  1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, i], " d:", sol$departureT[r, 
                  i], ") ", "->", sep = "")
            }
            cat("\n           ", sol$routes[r, which(sol$routes[r, ] == 0)[2]], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                i - 1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, which(sol$routes[r, ] == 0)[2]], 
                ")\n", sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 7 
#> Total distance: 1778 
#> Latest depot arrival (min from midnight): 290 
#> Routes:
#> Day: 1 driver: 1 customers: 0-25-4-23-6-21-0 departure: 0 finished: 290 distance: 246 quantity: 9628 
#>   Details: 0 (d:0) ->
#>            25 (tt:19 a:19 d:49) ->
#>            4 (tt:30 a:79 d:109) ->
#>            23 (tt:29 a:138 d:168) ->
#>            6 (tt:28 a:196 d:226) ->
#>            21 (tt:20 a:246 d:276) ->
#>            0 (tt:20 a:290)
#> Day: 3 driver: 1 customers: 0-2-1-28-26-0 departure: 20 finished: 271 distance: 68 quantity: 9421 
#>   Details: 0 (d:20) ->
#>            2 (tt:27 a:47 d:77) ->
#>            1 (tt:26 a:103 d:133) ->
#>            28 (tt:26 a:159 d:189) ->
#>            26 (tt:25 a:214 d:244) ->
#>            0 (tt:25 a:271)
#> Day: 3 driver: 2 customers: 0-9-19-20-18-0 departure: 13 finished: 222 distance: 319 quantity: 9657 
#>   Details: 0 (d:13) ->
#>            9 (tt:11 a:24 d:54) ->
#>            19 (tt:28 a:82 d:112) ->
#>            20 (tt:21 a:133 d:163) ->
#>            18 (tt:13 a:176 d:206) ->
#>            0 (tt:13 a:222)
#> Day: 4 driver: 1 customers: 0-12-8-30-16-0 departure: 89 finished: 290 distance: 237 quantity: 7582 
#>   Details: 0 (d:89) ->
#>            12 (tt:17 a:106 d:136) ->
#>            8 (tt:23 a:159 d:189) ->
#>            30 (tt:30 a:219 d:249) ->
#>            16 (tt:5 a:254 d:284) ->
#>            0 (tt:5 a:290)
#> Day: 4 driver: 2 customers: 0-14-15-13-22-0 departure: 71 finished: 286 distance: 135 quantity: 8790 
#>   Details: 0 (d:71) ->
#>            14 (tt:18 a:89 d:119) ->
#>            15 (tt:24 a:143 d:173) ->
#>            13 (tt:20 a:193 d:223) ->
#>            22 (tt:8 a:231 d:261) ->
#>            0 (tt:8 a:286)
#> Day: 5 driver: 2 customers: 0-27-24-11-7-17-0 departure: 10 finished: 290 distance: 446 quantity: 9975 
#>   Details: 0 (d:10) ->
#>            27 (tt:20 a:30 d:60) ->
#>            24 (tt:27 a:87 d:117) ->
#>            11 (tt:21 a:138 d:168) ->
#>            7 (tt:24 a:192 d:222) ->
#>            17 (tt:30 a:252 d:282) ->
#>            0 (tt:30 a:290)
#> Day: 6 driver: 2 customers: 0-5-3-29-10-0 departure: 1 finished: 213 distance: 327 quantity: 9020 
#>   Details: 0 (d:1) ->
#>            5 (tt:11 a:12 d:42) ->
#>            3 (tt:20 a:62 d:92) ->
#>            29 (tt:9 a:101 d:131) ->
#>            10 (tt:27 a:158 d:188) ->
#>            0 (tt:27 a:213)
```

Note there may not always be deliveries on a given day.

Time windows
------------

Finally, let us try to add time windows to the model above. I here formulate the time windows as soft constraints and try to minimize the total time not satisfing the time window.

``` r
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
lsp <- set.params(lsp, lsTimeLimit = c(20, 20, 20, 20, 20, 20), indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "totalRoutesUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "totalWaiting")
lsp <- add.output.expr(lsp, "totalViolateTWCtr")
lsp <- add.output.expr(lsp, "totalViolateTWMin")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$maxRoutes * data$days * data$drivers
lsp <- add.output.expr(lsp, "routes", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalT", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureT", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "waitingT", c(maxR, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "driver", maxR)
lsp <- add.output.expr(lsp, "day", maxR)
lsp <- add.output.expr(lsp, "routeDist", maxR)
lsp <- add.output.expr(lsp, "routeDepotArrival", maxR)
lsp <- add.output.expr(lsp, "routeQuant", maxR)
sol <- ls.solve(lsp, data)
length(grep("Feasible solution", readLines("output.txt"), value = TRUE, ignore.case = FALSE)) > 0
#> [1] TRUE
```

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$totalRoutesUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Total waiting time:", sol$totalWaiting, "\n")
    cat("Total time window violations:", sol$totalViolateTWCtr, "\n")
    cat("Total time window violations (min):", sol$totalViolateTWMin, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$driver)) {
        if (max(sol$routes[r, ]) > 0) {
            cat("Day:", sol$day[r], "driver:", sol$driver[r], "customers:", paste0(sol$routes[r, 1:(which(sol$routes[r, 
                ] == 0)[2])], collapse = "-"), "departure:", sol$departureT[r, 1], "finished:", sol$arrivalT[r, 
                which(sol$routes[r, ] == 0)[2]], "distance:", sol$routeDist[r], "quantity:", sol$routeQuant[r], 
                "\n")
            cat("  Details:", sol$routes[r, 1], paste0("(d:", sol$departureT[r, 1], ") ->"))
            for (i in 2:(which(sol$routes[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routes[r, i], " (tt:", data$travelTimeMatrix[sol$routes[r, i - 
                  1] + 1, sol$routes[r, i] + 1], " w:", sol$waitingT[r, i], " a:", sol$arrivalT[r, i], 
                  " d:", sol$departureT[r, i], ") ", " tw:[", data$timeWindowsLB[sol$routes[r, i] + 1], 
                  ",", data$timeWindowsUB[sol$routes[r, i] + 1], "] ->", sep = "")
            }
            cat("\n           ", sol$routes[r, which(sol$routes[r, ] == 0)[2]], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                i - 1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, which(sol$routes[r, ] == 0)[2]], 
                ")\n", sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 9 
#> Total distance: 2818 
#> Total waiting time: 3249 
#> Total time window violations: 0 
#> Total time window violations (min): 0 
#> Latest depot arrival (min from midnight): 1029 
#> Routes:
#> Day: 1 driver: 2 customers: 0-23-24-6-21-0 departure: 616 finished: 994 distance: 304 quantity: 6474 
#>   Details: 0 (d:616) ->
#>            23 (tt:22 w:0 a:638 d:668)  tw:[600,660] ->
#>            24 (tt:8 w:44 a:720 d:750)  tw:[600,720] ->
#>            6 (tt:16 w:134 a:900 d:930)  tw:[900,1020] ->
#>            21 (tt:20 w:0 a:950 d:980)  tw:[900,1020] ->
#>            0 (tt:20 a:994)
#> Day: 3 driver: 1 customers: 0-11-19-9-14-0 departure: 122 finished: 948 distance: 447 quantity: 9728 
#>   Details: 0 (d:122) ->
#>            11 (tt:15 w:193 a:330 d:360)  tw:[300,360] ->
#>            19 (tt:16 w:255 a:631 d:661)  tw:[600,720] ->
#>            9 (tt:28 w:31 a:720 d:750)  tw:[600,720] ->
#>            14 (tt:10 w:140 a:900 d:930)  tw:[900,960] ->
#>            0 (tt:10 a:948)
#> Day: 3 driver: 1 customers: 0-20-18-12-0 departure: 384 finished: 1000 distance: 299 quantity: 6791 
#>   Details: 0 (d:384) ->
#>            20 (tt:28 w:8 a:420 d:450)  tw:[300,420] ->
#>            18 (tt:13 w:437 a:900 d:930)  tw:[900,1020] ->
#>            12 (tt:23 w:0 a:953 d:983)  tw:[900,960] ->
#>            0 (tt:23 a:1000)
#> Day: 4 driver: 1 customers: 0-2-1-8-29-0 departure: 392 finished: 1029 distance: 246 quantity: 9347 
#>   Details: 0 (d:392) ->
#>            2 (tt:27 w:241 a:660 d:690)  tw:[600,660] ->
#>            1 (tt:26 w:184 a:900 d:930)  tw:[900,1020] ->
#>            8 (tt:26 w:0 a:956 d:986)  tw:[900,1020] ->
#>            29 (tt:5 w:0 a:991 d:1021)  tw:[900,1020] ->
#>            0 (tt:5 a:1029)
#> Day: 4 driver: 2 customers: 0-7-30-22-0 departure: 582 finished: 955 distance: 333 quantity: 7568 
#>   Details: 0 (d:582) ->
#>            7 (tt:29 w:0 a:611 d:641)  tw:[600,720] ->
#>            30 (tt:19 w:0 a:660 d:690)  tw:[600,660] ->
#>            22 (tt:17 w:193 a:900 d:930)  tw:[900,960] ->
#>            0 (tt:17 a:955)
#> Day: 5 driver: 1 customers: 0-26-28-27-0 departure: 73 finished: 950 distance: 351 quantity: 6979 
#>   Details: 0 (d:73) ->
#>            26 (tt:27 w:560 a:660 d:690)  tw:[600,660] ->
#>            28 (tt:25 w:5 a:720 d:750)  tw:[600,720] ->
#>            27 (tt:25 w:125 a:900 d:930)  tw:[900,1020] ->
#>            0 (tt:25 a:950)
#> Day: 5 driver: 2 customers: 0-17-15-13-0 departure: 283 finished: 454 distance: 265 quantity: 5518 
#>   Details: 0 (d:283) ->
#>            17 (tt:8 w:9 a:300 d:330)  tw:[300,420] ->
#>            15 (tt:16 w:0 a:346 d:376)  tw:[300,360] ->
#>            13 (tt:20 w:0 a:396 d:426)  tw:[300,420] ->
#>            0 (tt:20 a:454)
#> Day: 5 driver: 2 customers: 0-16-4-25-0 departure: 309 finished: 462 distance: 415 quantity: 6118 
#>   Details: 0 (d:309) ->
#>            16 (tt:6 w:0 a:315 d:345)  tw:[300,360] ->
#>            4 (tt:8 w:0 a:353 d:383)  tw:[300,420] ->
#>            25 (tt:30 w:0 a:413 d:443)  tw:[300,420] ->
#>            0 (tt:30 a:462)
#> Day: 6 driver: 2 customers: 0-5-3-10-0 departure: 149 finished: 995 distance: 158 quantity: 5550 
#>   Details: 0 (d:149) ->
#>            5 (tt:11 w:200 a:360 d:390)  tw:[300,360] ->
#>            3 (tt:20 w:490 a:900 d:930)  tw:[900,1020] ->
#>            10 (tt:10 w:0 a:940 d:970)  tw:[900,1020] ->
#>            0 (tt:10 a:995)
```

All time windows has been satisfied.

Pitfalls
========

To summarize some of the pitfalls I experienced during coding:

-   All input data must have be of the same data type as specified in the model (e.g. integers).
-   Remember when index from zero, this holds for all variables/data.
-   In `set.params` the parameter `lsTimeLimit=c(10,10)` must have the same length as the number of objectives. If e.g. `lsTimeLimit=10` it corresponds to `lsTimeLimit=c(0,10)` (only the last objective is optimized) and not `lsTimeLimit=c(10,10)` which seems more obvious.
-   You must specify the dimensions of output. Hence, when you run a new instance, you have to add output expressions again.
-   The dimension of output arrays cannot be more than 3 for some strange reason. So you may have to transform your output.
