// This file is part of the SGSolve library for stochastic games
// Copyright (C) 2016 Benjamin A. Brooks
//
// SGSolve free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SGSolve is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//
// Benjamin A. Brooks
// ben@benjaminbrooks.net
// Chicago, IL

//! One state prisoner's dilemma
//! @example

#include "sg.hpp"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP solve_rsg(double delta, int numStates, Rcpp::List states, int allIterations,
  double normtol,  double directiontol, double leveltol, double improvetol
)
{

  int action, state, player;
  int numPlayers = 2;

  int set_tol = 0;
  if (normtol >0 & directiontol > 0 & leveltol >0 & improvetol >0)
    set_tol = 1;

  vector<bool> unconstrained(2,false);

  //vector< vector< int > > numActions(numStates,vector<int>(numPlayers,2));
  //vector<int> numActions_total(numStates,4);

  vector< vector< int > > numActions(numStates);
  vector<int> numActions_total(numStates);

  // (numStates, numActions_total, numPlayers)
  vector< vector< vector<double> > >
    payoffs(numStates);

  vector< vector< vector<double> > >
    probabilities(numStates);

  cout << "\nCopy payoffs and transitions" << endl;
  for (state=0; state < numStates; state++) {
    cout << "Constructing state" << state << endl;

    Rcpp::List st = states[state];
    //Rcpp::IntegerVector na = st["numActions"];
    std::vector<int> nac = st["numActions"];
    numActions[state] = nac;
    int nat = std::accumulate(nac.begin(),nac.end(), 1, std::multiplies<int>());
    numActions_total[state] = nat;

    cout << "  payoffs... ";

    NumericMatrix rpm = st["payoffs"];
    vector<vector<double> > pm(nat, vector<double>(2,0.0));
    for (action = 0; action < nat; action++) {
      for (player = 0; player < numPlayers; player++) {
        pm[action][player] = rpm(action, player);
      }
    }
    payoffs[state] = pm;
    cout << "  done!" << endl;


    cout << "  transition matrix... ";
    NumericMatrix rtr = st["transition"];
    vector<vector<double> > tr(nat, vector<double>(numStates,0.0));
    for (action = 0; action < nat; action++) {
      for (int tostate = 0; tostate < numStates; tostate++) {
        tr[action][tostate] = rtr(action, tostate);
      }
    }
    probabilities[state] = tr;
    cout << "  done!" << endl;

  }

  //cout << "payoffs" << payoffs[0][0][0];
  /*
  Rcpp::List res = Rcpp::List::create(
    Rcpp::Named("numActions") = numActions,
    Rcpp::Named("numActions_total") = numActions_total,
    Rcpp::Named("payoffs") = payoffs,
    Rcpp::Named("probabilities") = probabilities
  );
  */

  try
    {
      cout << "Constructing game object" << endl;
      SGGame game(delta,
		  numStates,
		  numActions,
		  payoffs,
		  probabilities,
		  unconstrained);

      SGEnv env;

      if (set_tol == 1) {
        env.setParam(SG::DIRECTIONTOL,directiontol);
        env.setParam(SG::NORMTOL,normtol);
        env.setParam(SG::LEVELTOL,leveltol);
        env.setParam(SG::IMPROVETOL,improvetol);
      }

      cout << "Building solver" << endl;
      SGSolver solver(env,game);

      cout << "Starting solve routine" << endl;
      solver.solve();
      cout << "Done!" << endl;

      //cout << "\nDelta = " << game.getDelta()
      //     << "\nnumStates" << numStates;

      SGSolution sol = solver.getSolution();
      int numRevolutions = solver.getSolution().getIterations().back().getRevolution();
      // compute the number of iterations in the last revolution
      int numPoints = 0;
      std::list<SGIteration>::const_reverse_iterator it;
      it = solver.getSolution().getIterations().rbegin();
      while (it->getRevolution() == numRevolutions) {
        ++numPoints;
        ++it;
      }

      //vector<NumericMatrix>
      //  mpoints (numStates,NumericMatrix(numPoints, numPlayers));
      vector<NumericMatrix> mpoints (numStates);
      cout << "\nCopy results to R... ";

      for (int state = 0; state < numStates; state++) {
        //cout << "\n  State " << state << ": ";
        NumericMatrix spoints = NumericMatrix(numPoints, numPlayers);
        it = solver.getSolution().getIterations().rbegin();
        int pointInd = 0;
        while (it->getRevolution() == numRevolutions) {
          //cout << "\nIteration " << it->getIteration()
          //     << " Point " << pointInd << ":";
          for (int player = 0; player < numPlayers; player++) {
            double val = it->getPivot()[state][player];
            //cout << val << " ";
            //points[state][pointInd][player] = val;
            spoints(pointInd,player) = val;
          }
          ++pointInd;
          ++it;
          //cout << endl;
        }
        mpoints[state] = spoints;
        //cout << endl;
      }

      if (allIterations ==0) {
        cout << "... done!\n";
        //cout << "payoffs" << payoffs[0][0][0];
        Rcpp::List res = Rcpp::List::create(
          Rcpp::Named("solved") = 1,
          Rcpp::Named("points") = mpoints
        );

        return res;
      }

      // Return all iterations
      int numIterations = solver.getSolution().getIterations().back().getIteration();
      // compute the number of iterations in the last revolution
      numPoints = 0;
      std::list<SGIteration>::const_iterator fit;
      fit = solver.getSolution().getIterations().begin();
      while (fit->getIteration() <= numIterations) {
        ++numPoints;
        ++fit;
      }

      //vector<NumericMatrix>
      //  mpoints (numStates,NumericMatrix(numPoints, numPlayers));
      vector<NumericMatrix> ipoints (numStates);
      IntegerVector rev_vec(numPoints);


      for (int state = 0; state < numStates; state++) {
        //cout << "\n  State " << state << ": ";
        NumericMatrix spoints = NumericMatrix(numPoints, numPlayers);
        fit = solver.getSolution().getIterations().begin();
        int pointInd = 0;
        while (fit->getIteration() <= numIterations) {
          //cout << "\nIteration " << fit->getIteration()
          //     << " Point " << pointInd << ":";
          for (int player = 0; player < numPlayers; player++) {
            double val = fit->getPivot()[state][player];
            //cout << val << " ";
            //points[state][pointInd][player] = val;
            spoints(pointInd,player) = val;
          }
          if (state == 0) {
            rev_vec(pointInd) = fit->getRevolution();
          }
          ++pointInd;
          ++fit;
          //cout << endl;
        }
        ipoints[state] = spoints;
        //cout << endl;
      }
      //cout << "payoffs" << payoffs[0][0][0];
      Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("solved") = 1,
        Rcpp::Named("points") = mpoints,
        Rcpp::Named("ipoints") = ipoints,
        Rcpp::Named("revolution") = rev_vec
      );
      cout << "... done!" << endl;
      return(res);


    }
  catch (SGException e)
    {
      cout << "Caught the following exception:" << endl
	   << e.what() << endl;
      Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("payoffs") = payoffs,
        Rcpp::Named("probabilities") = probabilities,
        Rcpp::Named("solved") = 0,
        Rcpp::Named("msg") = e.what()
      );
      return(res);
    }
  return 0;
}
