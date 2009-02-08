/*
 *  CglLandPTabRow.hpp
 *
 *  Created by Pierre Bonami on 2/23/08.
 *  Copyright 2008 CNRS and others. All rights reserved.
 *
 */
#ifndef CglLandPTabRow_H
#define CglLandPTabRow_H

#include "CoinIndexedVector.hpp"
#include <iostream>

namespace LAP{
class CglLandPSimplex;
struct TabRow: public CoinIndexedVector {
  /** Row number.*/
  int num;
  /** Row right-hand-side.*/
  double rhs;
  /** Row of what?*/
  const CglLandPSimplex * si_;
  
  
  TabRow(const CglLandPSimplex *si):
    CoinIndexedVector(), num(-1), rhs(0), si_(si) {}
  TabRow(const TabRow & source):CoinIndexedVector(source),
  num(source.num), rhs(source.rhs) {
  }
  ~TabRow() {
  }
  
  void print(std::ostream & os, int width = 9, const int * nonBasics = NULL,
             int m = 0);
  inline
    const double& operator[](const int &index) const {
      return denseVector()[index];
    }
  
  inline
    double& operator[](const int &index) {
      return denseVector()[index];
    }
};
}/* Ends LAP Namespace.*/

#endif

