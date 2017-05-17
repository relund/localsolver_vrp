
Using LocalSolver to solve VRP problems
=======================================

The last days I have been playing a bit with [LocalSolver](http://www.localsolver.com/) and tried to model various vehicle routing problems (VRP). LocalSolver is a heuristic solver based on a heuristic search approach combining different optimization techniques. It includes a math modeling language and there is a package so that it can be used from [R](https://www.r-project.org/). You need to download and install LocalSolver (I used v7.0) and get an academic license from their webpage.

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

Note that all data are integers since I have chosen to use integers in LocalSolver (input data must be the same data type). Next, we formulate the model using [LSP](http://www.localsolver.com/documentation/lspreference/index.html). Here I used the example on LocalSolver's webpage as a starting point. The learning curve for LSP is quite steep; however, the [examples](http://www.localsolver.com/documentation/exampletour/index.html) help. I was missing a forum where to ask question though. An important issue is that you must choose from where the index in arrays start. I choose C style and let all index start from zero. We need decision variables `routeSequences[k]` containing a list of numbers in the range zero to `customers-1`. That is, if `routeSequences[k]` equals {0,4,6}, then truck mumber `k` visit customer 1, 5 and 7 (remember index start from zero).

``` r
library(localsolver)
model <- "function model(){
   routeSequences[k in 0..nbTrucks-1] <- list(nbCustomers);      // decision variable
   constraint partition[k in 0..nbTrucks-1](routeSequences[k]);  // each customer visited by exaclty one route

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
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE)) > 
    0
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
            cat("Customers:", paste0(sol$routeSequences[r, which(sol$routeSequences[r, ] != 
                -1)], collapse = "-"), "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r], 
                "\n")
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

Assume that only a single driver is avaliabele and that the maximum possible workday is specified. We need workday length, travel times and handling time.

``` r
library(localsolver)
source("simulate.R")
data <- simulateDriver(customers = 10, seed = 578, drivers = 1)
str(data)
#> List of 11
#>  $ nbTrucks          : int 5
#>  $ truckCapacity     : int 10000
#>  $ nbNodes           : int 11
#>  $ nbCustomers       : int 10
#>  $ demands           : int [1:10] 2514 2585 2938 3935 2649 1429 3766 3706 2778 2547
#>  $ distanceMatrix    : int [1:11, 1:11] 0 144 5 31 14 97 167 161 68 289 ...
#>  $ handlingTime      : int 30
#>  $ drivers           : int 1
#>  $ depotEarliestStart: int 360
#>  $ depotLatestArrival: int 3000
#>  $ travelTimeMatrix  : int [1:11, 1:11] 0 14 28 29 6 26 24 29 13 11 ...
```

To handle travel time constraints we add variables `routeStops[k][i]` which is the node visited at stop i (depot = node 0) and similar `arrivalTimes[k][i]` and `departureTimes[k][i]`.

``` r
model <- "function model(){
   routeSequences[k in 0..nbTrucks-1] <- list(nbCustomers);  
   routeStops[k in 0..nbTrucks-1][0] = 0;
   routeStops[k in 0..nbTrucks-1][nbNodes] = 0;
   departureTimes[k in 0..nbTrucks-1][i in 0..nbNodes] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..nbTrucks-1][i in 0..nbNodes] <- int(0,depotLatestArrival);

   constraint partition[k in 0..nbTrucks-1](routeSequences[k]);  // each customer visited by exaclty one route

   for [k in 0..nbTrucks-1] {
      local sequence <- routeSequences[k];
      local c <- count(sequence);
      trucksUsed[k] <- count(routeSequences[k]) > 0;
      // The quantity needed in each route must not exceed the truck capacity
      routeQuantity[k] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k] <= truckCapacity;

      // Distance travelled by truck k
      routeDistances[k] <- sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]) +
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0);

      if (k==0) constraint departureTimes[k][0] >= depotEarliestStart;  // first route 

      for [i in 0..nbNodes] {
         if (i!=0) {
            routeStops[k][i] <- routeSequences[k][i-1] + 1;
            arrivalTimes[k][i] <- departureTimes[k][i-1] + travelTimeMatrix[routeStops[k][i-1]][routeStops[k][i]];
            departureTimes[k][i] <- (routeStops[k][i] > 0  ? arrivalTimes[k][i] + handlingTime : 0);
         }
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

   nbTrucksUsed <- sum[k in 0..nbTrucks-1](trucksUsed[k]);
   totalDistance <- sum[k in 0..nbTrucks-1](routeDistances[k]);
   latestArrivalTime <- max[k in 0..nbTrucks-1](routeArrivalDepot[k]);

   minimize totalDistance;
   minimize latestArrivalTime;
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit = 20, indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$nbTrucks, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", data$nbTrucks)
lsp <- add.output.expr(lsp, "routeQuantity", data$nbTrucks)
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$nbTrucks, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$nbTrucks, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$nbTrucks, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", data$nbTrucks)
sol <- ls.solve(lsp, data)
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE)) > 
    0
#> [1] TRUE
```

We have minimized the total distance and next the lastest arrival time at the depot and postprocess the solution to get a result in a more readable format:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrivalTime, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$routeArrivalDepot)) {
        if (max(sol$routeStops[r, ]) > 0) {
            cat("Customers:", paste0(sol$routeStops[r, 1:(which(sol$routeStops[r, ] == 0)[2])], 
                collapse = "-"), "departure:", sol$departureTimes[r, 1], "finished:", sol$arrivalTimes[r, 
                which(sol$routeStops[r, ] == 0)[2]], "distance:", sol$routeDistances[r], "quantity:", 
                sol$routeQuantity[r], "\n")
            cat("  Details:", sol$routeStops[r, 1], paste0("(d:", sol$departureTimes[r, 1], 
                ") ->"))
            for (i in 2:(which(sol$routeStops[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routeStops[r, i], " (tt:", data$travelTimeMatrix[sol$routeStops[r, 
                  i - 1] + 1, sol$routeStops[r, i] + 1], " a:", sol$arrivalTimes[r, i], " d:", 
                  sol$departureTimes[r, i], ") ", "->", sep = "")
            }
            cat("\n           ", sol$routeStops[r, which(sol$routeStops[r, ] == 0)[2]], " (tt:", 
                data$travelTimeMatrix[sol$routeStops[r, i - 1] + 1, sol$routeStops[r, i] + 
                  1], " a:", sol$arrivalTimes[r, which(sol$routeStops[r, ] == 0)[2]], ")\n", 
                sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 4 
#> Total distance: 793 
#> Latest depot arrival (min from midnight): 901 
#> Routes:
#> Customers: 0-7-10-1-0 departure: 360 finished: 516 distance: 377 quantity: 8827 
#>   Details: 0 (d:360) ->
#>            7 (tt:29 a:389 d:419) ->
#>            10 (tt:17 a:436 d:466) ->
#>            1 (tt:6 a:472 d:502) ->
#>            0 (tt:6 a:516)
#> Customers: 0-5-9-3-0 departure: 516 finished: 704 distance: 158 quantity: 8365 
#>   Details: 0 (d:516) ->
#>            5 (tt:26 a:542 d:572) ->
#>            9 (tt:16 a:588 d:618) ->
#>            3 (tt:27 a:645 d:675) ->
#>            0 (tt:27 a:704)
#> Customers: 0-4-0 departure: 704 finished: 746 distance: 28 quantity: 3935 
#>   Details: 0 (d:704) ->
#>            4 (tt:6 a:710 d:740) ->
#>            0 (tt:6 a:746)
#> Customers: 0-8-6-2-0 departure: 746 finished: 901 distance: 230 quantity: 7720 
#>   Details: 0 (d:746) ->
#>            8 (tt:13 a:759 d:789) ->
#>            6 (tt:12 a:801 d:831) ->
#>            2 (tt:12 a:843 d:873) ->
#>            0 (tt:12 a:901)
```

Note currently there is no time added for filling the truck at the depot. Here the number of trucks is the maximum number of routes that can be used. Let us try another instance.

``` r
data <- simulateDriver(customers = 30, seed = 578, drivers = 1, trucks = 10)
sol <- ls.solve(lsp, data)
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE)) > 
    0
#> [1] TRUE
printSolution(sol)
#> Number of routes: 10 
#> Total distance: 3533 
#> Latest depot arrival (min from midnight): 1900 
#> Routes:
#> Customers: 0-28-29-13-27-0 departure: 360 finished: 562 distance: 259 quantity: 9701 
#>   Details: 0 (d:360) ->
#>            28 (tt:19 a:379 d:409) ->
#>            29 (tt:20 a:429 d:459) ->
#>            13 (tt:24 a:483 d:513) ->
#>            27 (tt:13 a:526 d:556) ->
#>            0 (tt:13 a:562)
#> Customers: 0-14-11-26-6-0 departure: 562 finished: 766 distance: 520 quantity: 8073 
#>   Details: 0 (d:562) ->
#>            14 (tt:23 a:585 d:615) ->
#>            11 (tt:7 a:622 d:652) ->
#>            26 (tt:16 a:668 d:698) ->
#>            6 (tt:14 a:712 d:742) ->
#>            0 (tt:14 a:766)
#> Customers: 0-8-3-16-0 departure: 766 finished: 939 distance: 296 quantity: 9980 
#>   Details: 0 (d:766) ->
#>            8 (tt:25 a:791 d:821) ->
#>            3 (tt:29 a:850 d:880) ->
#>            16 (tt:8 a:888 d:918) ->
#>            0 (tt:8 a:939)
#> Customers: 0-19-9-25-0 departure: 939 finished: 1070 distance: 295 quantity: 7252 
#>   Details: 0 (d:939) ->
#>            19 (tt:14 a:953 d:983) ->
#>            9 (tt:10 a:993 d:1023) ->
#>            25 (tt:11 a:1034 d:1064) ->
#>            0 (tt:11 a:1070)
#> Customers: 0-4-21-0 departure: 1070 finished: 1159 distance: 280 quantity: 6511 
#>   Details: 0 (d:1070) ->
#>            4 (tt:8 a:1078 d:1108) ->
#>            21 (tt:8 a:1116 d:1146) ->
#>            0 (tt:8 a:1159)
```

Finally, an instance which is not feasible since the travel time is to long:

``` r
data <- simulateDriver(customers = 10, seed = 578, drivers = 1, trucks = 10, travelTimes = c(30, 
    60))
sol <- ls.solve(lsp, data)
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE)) > 
    0
#> [1] TRUE
printSolution(sol)
#> Number of routes: 5 
#> Total distance: 786 
#> Latest depot arrival (min from midnight): 1347 
#> Routes:
#> Customers: 0-8-0 departure: 360 finished: 468 distance: 136 quantity: 3706 
#>   Details: 0 (d:360) ->
#>            8 (tt:39 a:399 d:429) ->
#>            0 (tt:39 a:468)
#> Customers: 0-2-0 departure: 468 finished: 614 distance: 10 quantity: 2585 
#>   Details: 0 (d:468) ->
#>            2 (tt:58 a:526 d:556) ->
#>            0 (tt:58 a:614)
```

There seems to be an error here in LocalSolver which writes out that there is a feasible solution, but the partition constraint does not hold:

``` r
sol$routeSequences
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]   -1   -1   -1   -1   -1   -1   -1   -1   -1    -1
#> [2,]    7   -1   -1   -1   -1   -1   -1   -1   -1    -1
#> [3,]    1   -1   -1   -1   -1   -1   -1   -1   -1    -1
#> [4,]   -1   -1   -1   -1   -1   -1   -1   -1   -1    -1
#> [5,]   -1   -1   -1   -1   -1   -1   -1   -1   -1    -1
```

Extension - Multiple drivers
----------------------------

Assume that multiple drivers are allowed and for simplicity assume that they have the same working time. To handle multiple drivers we need to add an driver index to the model. Note `nbTrucks` now becomes maximum number of trucks that can be used per driver.

``` r
data$drivers <- as.integer(2)
model <- "function model(){
   routeSequences[k in 0..nbTrucks-1][d in 0..drivers-1] <- list(nbCustomers);  
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][0] = 0;
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][nbNodes] = 0;
   departureTimes[k in 0..nbTrucks-1][d in 0..drivers-1][i in 0..nbNodes] <- int(-1,depotLatestArrival);
   arrivalTimes[k in 0..nbTrucks-1][d in 0..drivers-1][i in 0..nbNodes] <- int(-1,depotLatestArrival);

   constraint partition[k in 0..nbTrucks-1][d in 0..drivers-1](routeSequences[k][d]);

   for [k in 0..nbTrucks-1][d in 0..drivers-1] {
      local sequence <- routeSequences[k][d];
      local c <- count(sequence);
      trucksUsed[k][d] <- count(routeSequences[k][d]) > 0;
      // The quantity needed in each route must not exceed the truck capacity
      routeQuantity[k][d] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d] <= truckCapacity;

      // Distance travelled by truck k
      routeDistances[k][d] <- (c > 0 ? (distanceMatrix[0][sequence[0]+1] + 
         distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][0] >= depotEarliestStart;

      for [i in 0..nbNodes] {
         if (i!=0) {
            routeStops[k][d][i] <- sequence[i-1] + 1;
            arrivalTimes[k][d][i] <- departureTimes[k][d][i-1] + 
            travelTimeMatrix[routeStops[k][d][i-1]][routeStops[k][d][i]];
            departureTimes[k][d][i] <- (routeStops[k][d][i] > 0  ? arrivalTimes[k][d][i] + handlingTime : -1);
         }
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

   nbTrucksUsed <- sum[k in 0..nbTrucks-1][d in 0..drivers-1] (trucksUsed[k][d]);
   totalDistance <- sum[k in 0..nbTrucks-1][d in 0..drivers-1](routeDistances[k][d]);
   latestArrivalTime <- max[k in 0..nbTrucks-1][d in 0..drivers-1] (routeArrivalDepot[k][d]);

   minimize totalDistance;
   minimize latestArrivalTime;
}"
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit = c(20, 20), indexFromZero = TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path = getwd())
lsp <- add.output.expr(lsp, "routeSequences", c(data$nbTrucks, data$drivers, data$nbCustomers))
lsp <- add.output.expr(lsp, "routeDistances", c(data$nbTrucks, data$drivers))
lsp <- add.output.expr(lsp, "routeQuantity", c(data$nbTrucks, data$drivers))
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrivalTime")
lsp <- add.output.expr(lsp, "routeStops", c(data$nbTrucks, data$drivers, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "arrivalTimes", c(data$nbTrucks, data$drivers, data$nbNodes + 1))
lsp <- add.output.expr(lsp, "departureTimes", c(data$nbTrucks, data$drivers, data$nbNodes + 
    1))
lsp <- add.output.expr(lsp, "routeArrivalDepot", c(data$nbTrucks, data$drivers))
sol <- ls.solve(lsp, data)
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE)) > 
    0
#> [1] TRUE
```

We consided the same instance as before by now with 2 drivers, minimized the total distances and next the latest arrival time at the depot. The solution are:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrivalTime, "\n")
    cat("Routes:\n")
    for (d in 1:length(sol$routeArrivalDepot[1, ])) {
        for (r in 1:length(sol$routeArrivalDepot[, 1])) {
            if (max(sol$routeStops[r, d, ]) > 0) {
                cat("Driver:", d, "Customers:", paste0(sol$routeStops[r, d, 1:(which(sol$routeStops[r, 
                  d, ] == 0)[2])], collapse = "-"), "departure:", sol$departureTimes[r, d, 
                  1], "finished:", sol$arrivalTimes[r, d, which(sol$routeStops[r, d, ] == 0)[2]], 
                  "distance:", sol$routeDistances[r], "quantity:", sol$routeQuantity[r], "\n")
                cat("  Details:", sol$routeStops[r, d, 1], paste0("(d:", sol$departureTimes[r, 
                  d, 1], ") ->"))
                for (i in 2:(which(sol$routeStops[r, d, ] == 0)[2] - 1)) {
                  cat("\n           ", sol$routeStops[r, d, i], " (tt:", data$travelTimeMatrix[sol$routeStops[r, 
                    d, i - 1] + 1, sol$routeStops[r, d, i] + 1], " a:", sol$arrivalTimes[r, 
                    d, i], " d:", sol$departureTimes[r, d, i], ") ", "->", sep = "")
                }
                cat("\n           ", sol$routeStops[r, d, which(sol$routeStops[r, d, ] == 0)[2]], 
                  " (tt:", data$travelTimeMatrix[sol$routeStops[r, d, i - 1] + 1, sol$routeStops[r, 
                    d, i] + 1], " a:", sol$arrivalTimes[r, d, which(sol$routeStops[r, d, ] == 
                    0)[2]], ")\n", sep = "")
            }
        }
    }
}
printSolution(sol)
#> Number of routes: 4 
#> Total distance: 732 
#> Latest depot arrival (min from midnight): 925 
#> Routes:
#> Driver: 1 Customers: 0-4-0 departure: 408 finished: 500 distance: 28 quantity: 3935 
#>   Details: 0 (d:408) ->
#>            4 (tt:31 a:439 d:469) ->
#>            0 (tt:31 a:500)
#> Driver: 1 Customers: 0-5-9-3-0 departure: 609 finished: 914 distance: 158 quantity: 8365 
#>   Details: 0 (d:609) ->
#>            5 (tt:55 a:664 d:694) ->
#>            9 (tt:44 a:738 d:768) ->
#>            3 (tt:57 a:825 d:855) ->
#>            0 (tt:57 a:914)
#> Driver: 2 Customers: 0-1-6-8-0 departure: 360 finished: 629 distance: 0 quantity: 0 
#>   Details: 0 (d:360) ->
#>            1 (tt:41 a:401 d:431) ->
#>            6 (tt:60 a:491 d:521) ->
#>            8 (tt:39 a:560 d:590) ->
#>            0 (tt:39 a:629)
#> Driver: 2 Customers: 0-10-7-2-0 departure: 629 finished: 925 distance: 0 quantity: 0 
#>   Details: 0 (d:629) ->
#>            10 (tt:58 a:687 d:717) ->
#>            7 (tt:44 a:761 d:791) ->
#>            2 (tt:46 a:837 d:867) ->
#>            0 (tt:46 a:925)
```

Extension - Delivery over multiple days
---------------------------------------

Assume that the orders can be delivered at specific days (we ignore the time windows for the moment):

``` r
data <- simulateVRPTW(drivers = 2, trucks = 5, customers = 10, days = 4)
str(data)
#> List of 15
#>  $ nbTrucks          : int 5
#>  $ truckCapacity     : int 10000
#>  $ nbNodes           : int 11
#>  $ nbCustomers       : int 10
#>  $ demands           : int [1:10] 3694 3489 2151 1034 1465 2859 2354 3691 2471 2368
#>  $ distanceMatrix    : int [1:11, 1:11] 0 206 88 6 94 106 96 6 84 210 ...
#>  $ handlingTime      : int 30
#>  $ drivers           : int 2
#>  $ depotEarliestStart: int 360
#>  $ depotLatestArrival: int 3000
#>  $ travelTimeMatrix  : int [1:11, 1:11] 0 18 8 21 17 12 8 17 18 29 ...
#>  $ days              : int 4
#>  $ validDays         : int [1:10, 1:4] 0 0 1 0 1 0 0 0 1 0 ...
#>  $ timeWindowsLB     : int [1:11] 360 600 300 300 900 600 900 900 300 600 ...
#>  $ timeWindowsUB     : int [1:11] 3000 660 360 420 960 660 1020 1020 360 720 ...
data$validDays
#>       [,1] [,2] [,3] [,4]
#>  [1,]    0    0    1    1
#>  [2,]    0    1    1    1
#>  [3,]    1    1    1    0
#>  [4,]    0    0    1    1
#>  [5,]    1    1    1    1
#>  [6,]    0    0    1    1
#>  [7,]    0    1    0    0
#>  [8,]    0    1    1    1
#>  [9,]    1    1    1    1
#> [10,]    0    0    1    0
```

Valid days for an order corresponds to a one in the matrix `validDays`.

``` r
model <- '
function model(){
   routeSequences[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] <- list(nbCustomers);  
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] = 0;
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][nbNodes] = 0;
   departureTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,depotLatestArrival);

   constraint partition[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeSequences[k][d][t]);

   for [k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] {
      local sequence <- routeSequences[k][d][t];
      local c <- count(sequence);

      trucksUsed[k][d][t] <- count(sequence) > 0;

      routeQuantity[k][d][t] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d][t] <= truckCapacity;

      routeDistances[k][d][t] <-
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][t][0] >= depotEarliestStart;

      for [i in 0..nbNodes] {
         if (i!=0) {
            routeStops[k][d][t][i] <- sequence[i-1] + 1;
            arrivalTimes[k][d][t][i] <- (routeStops[k][d][t][i-1] > 0 || (i==1 && routeStops[k][d][t][i]>0) ?
               departureTimes[k][d][t][i-1] +
               travelTimeMatrix[routeStops[k][d][t][i-1]][routeStops[k][d][t][i]] : 0);
            departureTimes[k][d][t][i] <-
               (routeStops[k][d][t][i] > 0  ? arrivalTimes[k][d][t][i] + handlingTime : 0);
         }
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
      visitAtDay[t][i] <- sum[k in 0..nbTrucks-1][d in 0..drivers-1](indexOf(routeSequences[k][d][t],i) > -1);
      constraint visitAtDay[t][i] <= validDays[i][t];
   }

   nbTrucksUsed <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] (trucksUsed[k][d][t]);
   totalDistance <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeDistances[k][d][t]);
   latestArrival <- max[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeArrivalDepot[k][d][t]);

   // Objective: minimize the number of trucks used, then minimize the distance travelled
   minimize totalDistance;
   minimize latestArrival;

   // modify results so can output to R (cannot output array with more than 3 index)
   for [k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] {
         local r = k + d*nbTrucks + t*nbTrucks*drivers;
         routes[r][i in 0..nbNodes] <- routeStops[k][d][t][i];
         arrivalT[r][i in 0..nbNodes] <- arrivalTimes[k][d][t][i];
         routeDepotArrival[r] <- routeArrivalDepot[k][d][t];
         departureT[r][i in 0..nbNodes] <- departureTimes[k][d][t][i];
         routeDist[r] <- routeDistances[k][d][t];
         routeQuant[r] <- routeQuantity[k][d][t];
         driver[r] <- d+1;
         day[r] <- t+1;
   }
}
'
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(10,10), indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$nbTrucks*data$days*data$drivers
lsp <- add.output.expr(lsp, "routes", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "departureT", c(maxR, data$nbNodes+1))
lsp <- add.output.expr(lsp, "driver", maxR)
lsp <- add.output.expr(lsp, "day", maxR)
lsp <- add.output.expr(lsp, "routeDist", maxR)
lsp <- add.output.expr(lsp, "routeDepotArrival", maxR)
lsp <- add.output.expr(lsp, "routeQuant", maxR)
sol<-ls.solve(lsp, data)
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE))>0
#> [1] TRUE
```

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$driver)) {
        if (max(sol$routes[r, ]) > 0) {
            cat("Day:", sol$day[r], "driver:", sol$driver[r], "customers:", paste0(sol$routes[r, 
                1:(which(sol$routes[r, ] == 0)[2])], collapse = "-"), "departure:", sol$departureT[r, 
                1], "finished:", sol$arrivalT[r, which(sol$routes[r, ] == 0)[2]], "distance:", 
                sol$routeDist[r], "quantity:", sol$routeQuant[r], "\n")
            cat("  Details:", sol$routes[r, 1], paste0("(d:", sol$departureT[r, 1], ") ->"))
            for (i in 2:(which(sol$routes[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routes[r, i], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                  i - 1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, i], " d:", sol$departureT[r, 
                  i], ") ", "->", sep = "")
            }
            cat("\n           ", sol$routes[r, which(sol$routes[r, ] == 0)[2]], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                i - 1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, which(sol$routes[r, 
                ] == 0)[2]], ")\n", sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 4 
#> Total distance: 894 
#> Latest depot arrival (min from midnight): 227 
#> Routes:
#> Day: 2 driver: 2 customers: 0-7-0 departure: 96 finished: 160 distance: 12 quantity: 2354 
#>   Details: 0 (d:96) ->
#>            7 (tt:17 a:113 d:143) ->
#>            0 (tt:17 a:160)
#> Day: 3 driver: 1 customers: 0-8-4-10-6-0 departure: 0 finished: 227 distance: 374 quantity: 9952 
#>   Details: 0 (d:0) ->
#>            8 (tt:18 a:18 d:48) ->
#>            4 (tt:27 a:75 d:105) ->
#>            10 (tt:29 a:134 d:164) ->
#>            6 (tt:25 a:189 d:219) ->
#>            0 (tt:25 a:227)
#> Day: 3 driver: 1 customers: 0-2-0 departure: 139 finished: 185 distance: 176 quantity: 3489 
#>   Details: 0 (d:139) ->
#>            2 (tt:8 a:147 d:177) ->
#>            0 (tt:8 a:185)
#> Day: 3 driver: 2 customers: 0-3-9-1-5-0 departure: 9 finished: 212 distance: 332 quantity: 9781 
#>   Details: 0 (d:9) ->
#>            3 (tt:21 a:30 d:60) ->
#>            9 (tt:19 a:79 d:109) ->
#>            1 (tt:22 a:131 d:161) ->
#>            5 (tt:9 a:170 d:200) ->
#>            0 (tt:9 a:212)
```

Note there are no deliveries at day 1.

Time windows
------------

Finally, let us try to add time windows to the model above

``` r
model <- '
function model(){

   routeSequences[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] <- list(nbCustomers);  
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] = 0;
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][nbNodes] = 0;
   departureTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,depotLatestArrival);
   waitingTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,600);

   constraint partition[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeSequences[k][d][t]);

   for [k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] {
      local sequence <- routeSequences[k][d][t];
      local c <- count(sequence);

      trucksUsed[k][d][t] <- count(sequence) > 0;

      // The quantity needed in each route must not exceed the truck capacity
      routeQuantity[k][d][t] <- sum(0..c-1, i => demands[sequence[i]]);
      constraint routeQuantity[k][d][t] <= truckCapacity;

      // Distance travelled by truck k
      routeDistances[k][d][t] <-
         (c > 0 ? (distanceMatrix[0][sequence[0]+1] + distanceMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => distanceMatrix[sequence[i-1]+1][sequence[i]+1]);

      if (k==0) constraint departureTimes[k][d][t][0] >= depotEarliestStart;

      for [i in 0..nbNodes] {
         if (i!=0) {
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
      visitAtDay[t][i] <- sum[k in 0..nbTrucks-1][d in 0..drivers-1](indexOf(routeSequences[k][d][t],i) > -1);
      constraint visitAtDay[t][i] <= validDays[i][t];
   }

   nbTrucksUsed <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] (trucksUsed[k][d][t]);

   totalDistance <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeDistances[k][d][t]);

   totalWaiting <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes]
      (waitingTimes[k][d][t][i]);

   totalViolateTWCtr <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 1..nbCustomers](
      routeStops[k][d][t][i] > 0 && (arrivalTimes[k][d][t][i] < timeWindowsLB[routeStops[k][d][t][i]]  ||
      arrivalTimes[k][d][t][i] > timeWindowsUB[routeStops[k][d][t][i]])
   );
   totalViolateTWMin <- sum[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 1..nbCustomers](
      (routeStops[k][d][t][i] > 0 && arrivalTimes[k][d][t][i] < timeWindowsLB[routeStops[k][d][t][i]]) *
      (timeWindowsLB[routeStops[k][d][t][i]] - arrivalTimes[k][d][t][i]) +
      (routeStops[k][d][t][i] > 0 && arrivalTimes[k][d][t][i] > timeWindowsUB[routeStops[k][d][t][i]]) *
      (arrivalTimes[k][d][t][i] - timeWindowsUB[routeStops[k][d][t][i]])
   );

   latestArrival <- max[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1](routeArrivalDepot[k][d][t]);

   minimize totalViolateTWMin;
   minimize totalViolateTWCtr;
   minimize nbTrucksUsed;
   minimize totalDistance;
   minimize totalWaiting;
   minimize latestArrival;

   // modify results so can output to R (cannot output array with more than 3 index)
   for [k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] {
         local r = k + d*nbTrucks + t*nbTrucks*drivers;
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
}
'
rm(lsp)
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(10,10,10,10,10,10), indexFromZero=TRUE, lsNbThreads = 4)
lsp <- set.temp.dir(lsp, path=getwd())

lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "totalWaiting")
lsp <- add.output.expr(lsp, "totalViolateTWCtr")
lsp <- add.output.expr(lsp, "totalViolateTWMin")
lsp <- add.output.expr(lsp, "latestArrival")
maxR <- data$nbTrucks*data$days*data$drivers
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
length(grep("feasible solution", readLines("output.txt"), value = TRUE, ignore.case = TRUE))>0
#> [1] TRUE
```

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Total waiting time:", sol$totalWaiting, "\n")
    cat("Total time window violations:", sol$totalViolateTWCtr, "\n")
    cat("Total time window violations (min):", sol$totalViolateTWMin, "\n")
    cat("Latest depot arrival (min from midnight):", sol$latestArrival, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$driver)) {
        if (max(sol$routes[r, ]) > 0) {
            cat("Day:", sol$day[r], "driver:", sol$driver[r], "customers:", paste0(sol$routes[r, 
                1:(which(sol$routes[r, ] == 0)[2])], collapse = "-"), "departure:", sol$departureT[r, 
                1], "finished:", sol$arrivalT[r, which(sol$routes[r, ] == 0)[2]], "distance:", 
                sol$routeDist[r], "quantity:", sol$routeQuant[r], "\n")
            cat("  Details:", sol$routes[r, 1], paste0("(d:", sol$departureT[r, 1], ") ->"))
            for (i in 2:(which(sol$routes[r, ] == 0)[2] - 1)) {
                cat("\n           ", sol$routes[r, i], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                  i - 1] + 1, sol$routes[r, i] + 1], " w:", sol$waitingT[r, i], " a:", sol$arrivalT[r, 
                  i], " d:", sol$departureT[r, i], ") ", " tw:[", data$timeWindowsLB[sol$routes[r, 
                  i] + 1], ",", data$timeWindowsUB[sol$routes[r, i] + 1], "] ->", sep = "")
            }
            cat("\n           ", sol$routes[r, which(sol$routes[r, ] == 0)[2]], " (tt:", data$travelTimeMatrix[sol$routes[r, 
                i - 1] + 1, sol$routes[r, i] + 1], " a:", sol$arrivalT[r, which(sol$routes[r, 
                ] == 0)[2]], ")\n", sep = "")
        }
    }
}
printSolution(sol)
#> Number of routes: 3 
#> Total distance: 1548 
#> Total waiting time: 1146 
#> Total time window violations: 0 
#> Total time window violations (min): 0 
#> Latest depot arrival (min from midnight): 1008 
#> Routes:
#> Day: 2 driver: 2 customers: 0-3-2-7-0 departure: 295 finished: 947 distance: 440 quantity: 7994 
#>   Details: 0 (d:295) ->
#>            3 (tt:21 w:2 a:318 d:348)  tw:[300,420] ->
#>            2 (tt:12 w:0 a:360 d:390)  tw:[300,360] ->
#>            7 (tt:14 w:496 a:900 d:930)  tw:[900,1020] ->
#>            0 (tt:14 a:947)
#> Day: 3 driver: 1 customers: 0-8-5-4-10-0 departure: 342 finished: 1008 distance: 485 quantity: 8558 
#>   Details: 0 (d:342) ->
#>            8 (tt:18 w:0 a:360 d:390)  tw:[300,360] ->
#>            5 (tt:11 w:220 a:621 d:651)  tw:[600,660] ->
#>            4 (tt:10 w:239 a:900 d:930)  tw:[900,960] ->
#>            10 (tt:29 w:0 a:959 d:989)  tw:[900,1020] ->
#>            0 (tt:29 a:1008)
#> Day: 3 driver: 2 customers: 0-9-1-6-0 departure: 579 finished: 938 distance: 623 quantity: 9024 
#>   Details: 0 (d:579) ->
#>            9 (tt:29 w:0 a:608 d:638)  tw:[600,720] ->
#>            1 (tt:22 w:0 a:660 d:690)  tw:[600,660] ->
#>            6 (tt:21 w:189 a:900 d:930)  tw:[900,1020] ->
#>            0 (tt:21 a:938)
```

Note there are no deliveries at day 1.
