/* Compile it with (assuming COIN lives in $HOME/COIN):
  g++ -g -o test2mir CglTwomirTest.cpp -I$HOME/COIN/include -L$HOME/COIN/lib -Wl,-rpath,$HOME/COIN/lib -lSbb -lCgl -lOsiClp -lClp -lOsi -lCoin -lm

  Run it with:
  test2mir <mps_file_name>
*/

#include "CglTwomir.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "SbbModel.hpp"
#include "SbbHeuristic.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
  OsiClpSolverInterface *osip = new OsiClpSolverInterface ();
  int numMpsReadErrors = osip->readMps(argv[1]);

  double  value = 0.0;
  SbbModel sbb(*osip);
  sbb.initialSolve();
  sbb.setMaximumNodes(100);

  CglTwomir mir_tab;
  mir_tab.setMirScale (1,1); mir_tab.setTwomirScale (1,5); mir_tab.setAMax (5);
  mir_tab.setCutTypes (true, false, true, false);

  CglTwomir mir_form;
  mir_form.setMirScale (1,5); mir_form.setTwomirScale (1,5); mir_form.setAMax (5);
  mir_form.setCutTypes (true, false, false, true);

  CglTwomir twomir_tab;
  twomir_tab.setMirScale (1,5); twomir_tab.setTwomirScale (1,5); twomir_tab.setAMax (5);
  twomir_tab.setCutTypes (false, true, true, false);

  CglTwomir twomir_form;
  twomir_form.setMirScale (1,5); twomir_form.setTwomirScale (1,5); twomir_form.setAMax (5);
  twomir_form.setCutTypes (false, true, false, true);
  // void setCutTypes (bool mir, bool twomir, bool tab, bool form)
  
  sbb.addCutGenerator(&mir_tab,1,"Mir_tab");
  sbb.addCutGenerator(&mir_form,1,"Mir_form");
  sbb.addCutGenerator(&twomir_tab,1,"Twomir_tab");
  sbb.addCutGenerator(&twomir_form,1,"Twomir_form");

  sbb.branchAndBound();
}  


  /*
  If you want to test interaction with other cut generators:

  SbbRounding sbround(sbb);
  int res = !sbround.solution (value, result);
  CglProbing generator1;  generator1.setUsingObjective(true);  generator1.setMaxPass(3);  
  generator1.setMaxProbe(100);    generator1.setMaxLook(50);  generator1.setRowCuts(3);
  CglKnapsackCover generator3;  CglOddHole generator4;  generator4.setMinimumViolation(0.005);
  generator4.setMinimumViolationPer(0.00002);  generator4.setMaximumEntries(200);
  // Add in generators
  sbb.addCutGenerator(&generator1,-1,"Probing");
  sbb.addCutGenerator(&generator2,-99,"Gomory");
  sbb.addCutGenerator(&generator3,-1,"Knapsack");
  sbb.addCutGenerator(&generator4,-1,"OddHole");
  // Allow rounding heuristic
  SbbRounding heuristic1(sbb);
  sbb.addHeuristic(&heuristic1);
  */
