library(localsolver)
## Example in thesis
# NOTE index starts from zero when passed to localsolver (must be consistent with the model)
# comments below are index after passed to localsolver (think the depot as node -1)
nbTrucks <- as.integer(4)
truckCapacity <- as.integer(37000)
nbNodes <- as.integer(12)
nbCustomers <- as.integer(nbNodes - 1)
demands <- as.integer(  # not including the depot
   c(15000, 2000, 10000, 2000, 27000, 10000, 2000, 7000, 2000, 7000, 10000)
)
# distanceMatrix[i,j] = distance between customer i and j.
con <- textConnection(
"0	97	130	59	33	23	113	36	209	73	135	160
97	0	54	49	125	85	185	71	134	145	124	62
130	54	0	82	158	118	218	104	80	178	98	53
59	49	82	0	87	47	148	34	161	107	113	111
33	125	158	87	0	50	141	63	237	100	162	188
23	85	118	47	50	0	100	20	197	60	122	148
113	185	218	148	141	100	0	125	298	42	224	249
36	71	104	34	63	20	125	0	184	83	100	134
209	134	80	161	237	197	298	184	0	257	134	135
73	145	178	107	100	60	42	83	257	0	182	208
135	124	98	113	162	122	224	100	134	182	0	177
160	62	53	111	188	148	249	134	135	208	177	0")
distanceMatrix <- read.table(con, sep="\t")
close(con)
distanceMatrix <- as.matrix(distanceMatrix)
diag(distanceMatrix) <- 0
storage.mode(distanceMatrix) <- "integer"

# travelTimeMatrix[i,j] = travel time between customer i and j.
con <- textConnection(
   "0	110	130	89	74	66	115	77	172	93	134	150
110	0	89	84	128	104	158	97	130	135	130	103
130	89	0	106	150	126	179	119	100	157	139	93
89	84	106	0	107	84	137	77	146	115	128	123
74	128	150	107	0	83	133	94	189	111	151	167
66	104	126	84	83	0	107	72	167	85	129	145
115	158	179	137	133	107	0	125	221	75	183	198
77	97	119	77	94	72	125	0	161	102	121	139
172	130	100	146	189	167	221	161	0	200	158	139
93	135	157	115	111	85	75	102	200	0	160	175
134	130	139	128	151	129	183	121	158	160	0	165
150	103	93	123	167	145	198	139	139	175	165	0")
travelTimeMatrix <- read.table(con, sep="\t")
close(con)
travelTimeMatrix <- as.matrix(travelTimeMatrix)
diag(travelTimeMatrix) <- 0
# travelTimeWarehouse[i] = travel time from warehouse to customer i
#travelTimeWarehouse <- travelTimeMatrix[1,2:nbNodes]
# remove warehouse travel times
#travelTimeMatrix <- travelTimeMatrix[2:nbNodes, 2:nbNodes]
storage.mode(travelTimeMatrix) <- "integer"

handlingTime <- as.integer(30)   # handling time at node in minutes

toMin<-function(hour, min=0) {60*hour + min}  # minutes since midnight
depotEarliestStart <- as.integer(toMin(5))  # earliest and latest depature time
depotLatestArrival <- as.integer(toMin(50))

drivers <- as.integer(2)
# note nbTrucks = max number of routes done in total by EACH driver

days = as.integer(4)
validDays<-matrix(c(
1, 1, 1, 0,
1, 1, 1, 0,
0, 1, 0, 0,
0, 1, 1, 1,
0, 0, 1, 1,
0, 1, 0, 0,
0, 1, 1, 1,
0, 0, 1, 1,
0, 1, 1, 1,
0, 0, 1, 1,
0, 1, 0, 0), nrow = nbCustomers, ncol = days)
storage.mode(validDays) <- "integer"

set.seed(1234)
timeWindowsLB <- sample(c(5,10,15), nbNodes, replace = TRUE)
timeWindowsUB <- timeWindowsLB + 2
timeWindowsLB <- as.integer(toMin(timeWindowsLB))
timeWindowsUB <- as.integer(toMin(timeWindowsUB))
timeWindowsLB[1] <- depotEarliestStart
timeWindowsUB[1] <- depotLatestArrival

## CVRP model with multiple drivers and days and time windows as soft constraints

model <- '
function model(){
   // Sequence of customers visited by each truck
   routeSequences[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1] <- list(nbCustomers);  // note i in the list is node i+1
   // routeStops[k][i] = node visited at stop i (0 = depot)
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] = 0;
   routeStops[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][nbNodes] = 0;
   //arrivalTimes[k][i] = arrival time truck k at stop i (0 = depot)
   departureTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][0] <- int(0,depotLatestArrival);
   arrivalTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,depotLatestArrival);
   waitingTimes[k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1][i in 0..nbNodes] <- int(0,600);

   // All customers must be visited by the trucks
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
         /*departureTimes[k][d][t][0] +
         (c > 0 ? (travelTimeMatrix[0][sequence[0]+1] + travelTimeMatrix[sequence[c-1]+1][0]) : 0) +
         sum(1..c-1, i => travelTimeMatrix[sequence[i-1]+1][sequence[i]+1]) + c*handlingTime;*/
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

   // Objective: minimize the number of trucks used, then minimize the distance travelled
   minimize totalViolateTWMin;
   minimize totalViolateTWCtr;
   minimize nbTrucksUsed;
   minimize totalDistance;
   minimize totalWaiting;
   minimize latestArrival;

   // modify results so can output to R (cannot output array with more than 3 index)
   //routes[r in 0..days*drivers*nbTrucks-1][i in 0..nbNodes] <- -1;//int(0,nbNodes);
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

function param(){
   cust = 0;
   for [k in 0..nbTrucks-1][d in 0..drivers-1][t in 0..days-1]{
      if (validDays[cust][t]) {
         routeSequences[k][d][t].value.add(cust);
         departureTimes[k][d][t][0].value = timeWindowsLB[cust+1] - travelTimeMatrix[cust+1][0] >= depotEarliestStart ?
            timeWindowsLB[cust+1] - travelTimeMatrix[cust+1][0] : depotEarliestStart;
         waitingTimes[k][d][t][1].value = 0;
         //arrivalTimes[k][d][t][1].value = departureTimes[k][d][t][0].value + travelTimeMatrix[cust+1][0];
         //departureTimes[k][d][t][1].value = arrivalTimes[k][d][t][1].value + handlingTime;
         //arrivalTimes[k][d][t][2].value = departureTimes[k][d][t][1].value + travelTimeMatrix[cust+1][0];
         cust = cust + 1;
         if (cust == nbCustomers) break;
      }
   }


}
'
rm(lsp)
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=c(10,10,10,10,10,10), indexFromZero=TRUE, lsNbThreads = 4,
                  lsVerbosity = 1, lsAnnealingLevel = 1)
lsp <- set.temp.dir(lsp, path=getwd())
lsp <- add.output.expr(lsp, "nbTrucksUsed")
lsp <- add.output.expr(lsp, "totalDistance")
lsp <- add.output.expr(lsp, "totalWaiting")
lsp <- add.output.expr(lsp, "totalViolateTWCtr")
lsp <- add.output.expr(lsp, "totalViolateTWMin")
lsp <- add.output.expr(lsp, "latestArrival")
lsp <- add.output.expr(lsp, "routes", c(nbTrucks*days*drivers, nbNodes+1))
lsp <- add.output.expr(lsp, "arrivalT", c(nbTrucks*days*drivers, nbNodes+1))
lsp <- add.output.expr(lsp, "departureT", c(nbTrucks*days*drivers, nbNodes+1))
lsp <- add.output.expr(lsp, "waitingT", c(nbTrucks*days*drivers, nbNodes+1))
lsp <- add.output.expr(lsp, "driver", c(nbTrucks*days*drivers))
lsp <- add.output.expr(lsp, "day", c(nbTrucks*days*drivers))
lsp <- add.output.expr(lsp, "routeDist", c(nbTrucks*days*drivers))
lsp <- add.output.expr(lsp, "routeDepotArrival", c(nbTrucks*days*drivers))
lsp <- add.output.expr(lsp, "routeQuant", c(nbTrucks*days*drivers))
lsp <- add.output.expr(lsp, "timeWindowsLB", c(nbCustomers))

data <- list(nbTrucks=nbTrucks, truckCapacity=truckCapacity, nbNodes=nbNodes,
             nbCustomers=nbCustomers, demands=demands, distanceMatrix=distanceMatrix,
             depotEarliestStart=depotEarliestStart, validDays=validDays,
             depotLatestArrival=depotLatestArrival, days=days,
             travelTimeMatrix=travelTimeMatrix, timeWindowsLB=timeWindowsLB,
             timeWindowsUB=timeWindowsUB,
             handlingTime=handlingTime, drivers=drivers)
unlink("input.lsp")
unlink("data.txt")
unlink("error.txt")
unlink("output.txt")
sol<-ls.solve(lsp, data)
printSolution<-function(sol) {
   cat("Number of routes:", sol$nbTrucksUsed, "\n")
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
             "distance:", sol$routeDist[r], "quantity:", sol$routeQuant[r],
             "\n")
         cat("  Details:", sol$routes[r,1], paste0("(d:", sol$departureT[r,1], ") ->"))
         for (i in 2:(which(sol$routes[r,]==0)[2]-1)) {
            cat("\n           ", sol$routes[r,i], " (tt:", travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1],
                " w:", sol$waitingT[r,i], " a:", sol$arrivalT[r,i], " d:", sol$departureT[r,i], ") ",
                " tw:[", timeWindowsLB[sol$routes[r,i]+1], ",", timeWindowsUB[sol$routes[r,i]+1], "] ->", sep="")
         }
         cat("\n           ", sol$routes[r,which(sol$routes[r,]==0)[2]],
             " (tt:", travelTimeMatrix[sol$routes[r,i-1]+1,sol$routes[r,i]+1],
             " a:", sol$arrivalT[r,which(sol$routes[r,]==0)[2]], ")\n", sep="")
      }
   }
}
printSolution(sol)

