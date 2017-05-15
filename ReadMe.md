
Using LocalSolver to solve VRP problems
=======================================

The last days I have been playing a bit with [LocalSolver](http://www.localsolver.com/) and tried to model various vehicle routing problems (VRP). LocalSolver is a heuristic solver based on a heuristic search approach combining different optimization techniques. It includes a math modeling language and there is a package so that it can be used from \[R\](<https://www.r-project.org/>. You need to download and install LocalSolver (I used v7.0) and get an academic license.

Capacitated VRP (CVRP)
----------------------

First, we simulate some data:

``` r
source("simulate.R")
data <- simulateCVRP(customers = 10, seed = 578)
str(data)
```

    ## List of 7
    ##  $ nbTrucks      : int 5
    ##  $ truckCapacity : int 10000
    ##  $ nbNodes       : int 11
    ##  $ nbCustomers   : int 10
    ##  $ demands       : int [1:10] 2514 2585 2938 3935 2649 1429 3766 3706 2778 2547
    ##  $ distanceMatrix: int [1:11, 1:11] 0 144 5 31 14 97 167 161 68 289 ...
    ##  $ handlingTime  : int 30

Note that all data are integers since I have chosen to use integers in LocalSolver (input data must be the same data type). Next, we formulate the model using [LSP](http://www.localsolver.com/documentation/lspreference/index.html). The learning curve for LSP is quite steep; however, the [examples](http://www.localsolver.com/documentation/exampletour/index.html) help. I was missing a forum where to ask question though. An important issue is that you must choose from where the index in arrays start. I choose C style and let all index start from zero. We need decision variables `customersSequences[k]` containing a list of numbers in the range zero to `customers-1`. That is, if `customersSequences[k]` equals {0,4,6}, then truck mumber `k` visit customer 1, 5 and 7 (remember index start from zero).

``` r
library(localsolver)

model <- "function model(){
   customersSequences[k in 0..nbTrucks-1] <- list(nbCustomers);      // decision variable
   constraint partition[k in 0..nbTrucks-1](customersSequences[k]);  // all customers must be visited by exaclty one truck

   trucksUsed[k in 0..nbTrucks-1] <- count(customersSequences[k]) > 0;  
   nbTrucksUsed <- sum[k in 0..nbTrucks-1](trucksUsed[k]);

   for [k in 0..nbTrucks-1] {
      local sequence <- customersSequences[k];
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
lsp <- add.output.expr(lsp, "customersSequences", c(data$nbTrucks, data$nbCustomers))
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
```

    ## [1] TRUE

We postprocess the solution to get a result in a readable format:

``` r
printSolution <- function(sol) {
    cat("Number of routes:", sol$nbTrucksUsed, "\n")
    cat("Total distance:", sol$totalDistance, "\n")
    cat("Routes:\n")
    for (r in 1:length(sol$routeDistances)) {
        if (max(sol$customersSequences[r, ]) >= 0) {
            cat("Customers:", paste0(sol$customersSequences[r, which(sol$customersSequences[r, 
                ] != -1)], collapse = "-"), "distance:", sol$routeDistances[r], 
                "quantity:", sol$routeQuantity[r], "\n")
        }
    }
}
printSolution(sol)
```

    ## Number of routes: 3 
    ## Total distance: 1147 
    ## Routes:
    ## Customers: 4-1-3 distance: 222 quantity: 9169 
    ## Customers: 7-6-0 distance: 559 quantity: 9986 
    ## Customers: 2-5-8-9 distance: 366 quantity: 9692