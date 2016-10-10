example.solve.pd = function() {
  library(RSGSolve)
  setwd("D:/libraries/dyngame")
  rsg = loadJsonSG("pd_abs.json")
  #rsg = loadJsonSG("cournot.json")

  delta = rsg$delta
  states = rsg$states
  sol = solveSG(rsg$delta, states, all.iter = TRUE)
  sol

  solveSG(rsg$delta,states, noreturn=TRUE)
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
#'  @param duplicate.first.point if TRUE (default) the first row o the point matrices will be added again to the end. This facilitates plotting of the payoff set using the lines command.
#' @param all.iter shall return value contain the field ipoints that contains the pivots of all iterations?
#' @param noreturn if set to TRUE solves the game but does not return any values. Useful for bechmarking and not taking into account the time it takes to copy results into R format. (My C++ code there may not be super efficient.)
#' @return a list with the following elements
#'    solved: TRUE if the game could be solved
#'
#'    points: a list with a matrix for every state that
#'      contains the the extreme points of the final
#'      payoff set approximation.
#'      Using the terminology of ABS (2016): the pivots
#'      from the last revolution.
#'
#'    ipoints: a list with a matrix for every state, containing
#'      the pivots (extreme points) from every iteration.
#'
#'    revolution: a vector that denotes the revolution of each
#'      point in the ipoints matrices.
#'
#' @export
solveSG = function(delta = rsg$delta, states=rsg$states, rsg=NULL,duplicate.first.point = TRUE, all.iter=TRUE, tol=1e-12, normtol=tol,  directiontol=tol, leveltol=tol, improvetol=tol, verbose=1L, noreturn=FALSE) {
  restore.point("solve_rsg")

  all.iter = as.integer(all.iter)
  verbose = as.integer(verbose)
  noreturn = as.integer(noreturn)
  sol = solve_rsg(delta, length(states), states, all.iter,verbose,noreturn, normtol, directiontol, leveltol,improvetol)
  sol$solved = (sol$solved == 1)
  if (noreturn==1) return(sol)

  if (sol$solved) {
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
#' You then can call solveSG on the result
#' to get the approximations to the equilibrium payoff sets
#' @param file the file name of the json file
#' @export
loadJsonSG = function(file) {
  json = readLines(file, warn=FALSE)
  rsg = jsonlite::fromJSON(json,flatten = FALSE,simplifyDataFrame = FALSE)
  rsg
}
