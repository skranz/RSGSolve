example.solve.pd = function() {
  library(RSGSolve)
  setwd("D:/libraries/dyngame")
  rsg = loadJsonSG("pd_abs.json")
  #rsg = loadJsonSG("cournot.json")

  delta = rsg$delta
  states = rsg$states
  sol = solveSG(rsg$delta, states, all.iter = TRUE)
  sol

  plot(sol$ipoints[[1]], type="p", lwd=2, col="grey")
  lines(sol$points[[1]], type="o", lwd=2, col="blue")

}

#' Solve a two-player stochastic game of perfect monitoring
#'
#' Uses Ben Brooks implementation SGSolve of the pencil sharpening
#' algorithm of Abreu, Brooks and Sanikov (2016)
#' See http://www.benjaminbrooks.net/sgsolvedoc/html/index.html
#'
#' @param delta the discount factor between 0 and 1
#' @param stages a list of stages. Each stage is a list with the following elements:
#'  - numActions: a size 2 vector that contains the number of actions for player 1 and 2. Its product is numActionProfiles.
#'  - payoffs: a matrix of dimension numActionProfiles x 2. The payoffs as function of the action profile. The actions are mapped to action profiles as described in the documentation of Brooks C++ library: A pair (a1,a2) is mapped into an action profile index using the formula a=a1+a2*numActions[s][a1].
#'  - transition: a matrix of dimension numActionProfiles x numStates. The transition probabilities. Each row has to sum up to 1.
#' @return a list where points contains the extreme points
#'         of the final payoff set approximation for every state.
#'         Using the terminology of ABS (2016): the pivots from the last revolution.
#'         While SGSolve computes much more this R Interface only returns these values so far.
#' @export
solveSG = function(delta = rsg$delta, states=rsg$states, rsg=NULL,duplicate.first.point = TRUE, all.iter=TRUE, tol=1e-12, normtol=tol,  directiontol=tol, leveltol=tol, improvetol=tol) {
  restore.point("solve_rsg")

  all.iter = as.integer(all.iter)
  sol = solve_rsg(delta, length(states), states, all.iter,normtol, directiontol, leveltol,improvetol)

  if (sol$solved==1) {
    if (duplicate.first.point) {
      for (i in seq_along(states)) {
        points = sol$points[[i]]
        sol$points[[i]] = points[c(seq.int(NROW(points)),1),]
      }
    }
  }
  return(sol)
}

#' Loads a json file rsg game specification
#' and returns the rsg object as R list
#' You then can call rsg_solve on the result
#' to get the approximations to the equilibrium payoff sets
#' @param file the file name of the json file
#' @export
loadJsonSG = function(file) {
  json = readLines(file, warn=FALSE)
  rsg = jsonlite::fromJSON(json,flatten = FALSE,simplifyDataFrame = FALSE)
  rsg
}