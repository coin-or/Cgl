/*
 *  CglLandPTabRow.cpp
 *  LandP
 *
 *  Created by Pierre Bonami on 2/23/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "CglLandPTabRow.h"
#include "CglLandPSimplex.hpp"
namespace LAP {
void
TabRow::print(std::ostream & os, int width, const int * nonBasics,
              int m) {
    os.width(3);
    os.precision(4);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<"idx: ";
    const double * dense = denseVector();
    for (int j = 0 ; j < m ; j++) {
        os.width(width);
        os.setf(std::ios_base::right, std::ios_base::adjustfield);
        os<<nonBasics[j]<<" ";
    }

    os<<std::endl;
    os.width(3);
    os.precision(4);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<< num <<": ";
    for (int j = 0 ; j < m ; j++) {
        os.width(width);
        os.precision(3);
        //      os.setf(std::ios_base::fixed, std::ios_base::floatfield);
        os.setf(std::ios_base::right, std::ios_base::adjustfield);
        os<<dense[nonBasics[j]]<<" ";
    }

    os.width(width);
    os.precision(4);
    //    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<rhs;

    os<<std::endl;

}


}/* Ends namespace LAP.*/
