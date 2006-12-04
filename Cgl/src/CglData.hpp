// Name:     CglData.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     12/2/06
//-----------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#ifndef CglData_H
#define CglData_H

/** Class collecting pointers on data for all cut generators. Each generator
    may have a derived class to add additional pointers on data. 
    If a data member is not used by
    a generator, the data member need not be defined (or may be NULL).
    If a data member is NULL and the generator needs the data, the call to
    generateCuts() is aborted. If a data member is not NULL, the calling method
    is responsible for the accuracy of the data. 
    Ownership of the data remains with the calling method.
*/

#include "CoinPackedMatrix.hpp"
#include "CglTreeInfo.hpp"

class CglData {

public:

  /**@name Public Set/get methods */
  //@{

  /** Set nrow to the number of rows */
  virtual void setNrow(const int givenNrow);
  /** Get nrow */
  inline int getNrow() const {return nrow;};

  /** Set ncol to the number of variables */
  virtual void setNcol(const int givenNcol);
  /** Get ncol */
  inline int getNcol() const {return ncol;};

  /** Set  matrixByCol to point on the coefficient matrix ordered by 
      columns */
  virtual void setMatrixByCol(const CoinPackedMatrix *givenMatrixByCol);
  /** Get matrixByCol */
  inline const CoinPackedMatrix *getMatrixByCol() const {return matrixByCol;};

  /** Set  matrixByRow to point on the coefficient matrix ordered by 
      rows */
  virtual void setMatrixByRow(const CoinPackedMatrix *givenMatrixByRow);
  /** Get matrixByRow */
  inline const CoinPackedMatrix *getMatrixByRow() const {return matrixByRow;};

  /** Set obj to point on a vector holding the objective coefficient values */
  virtual void setObj(const double *givenObj);
  /** Get obj */
  inline const double *getObj() const {return obj;};

  /** Set  colLower to point on a vector holding the lower bounds on the 
      variables */
  virtual void setColLower(const double *givenColLower);
  /** Get colLower */
  inline const double *getColLower() const {return colLower;};

  /** Set  colUpper to point on a vector holding the upper bounds on the 
      variables */
  virtual void setColUpper(const double *givenColUpper);
  /** Get colUpper */
  inline const double *getColUpper() const {return colUpper;};

  /** Set  rowLower to point on a vector holding the lower bounds on the 
      constraints */
  virtual void setRowLower(const double *givenRowLower);
  /** Get rowLower */
  inline const double *getRowLower() const {return rowLower;};

  /** Set  rowUpper to point on a vector holding the upper bounds on the 
      constraints */
  virtual void setRowUpper(const double *givenRowUpper);
  /** Get rowUpper */
  inline const double *getRowUpper() const {return rowUpper;};

  /** Set  rowRhs to point on a vector holding the right hand side of the 
      constraints (for a ranged constraint, it contains the upper bound). */
  virtual void setRowRhs(const double *givenRowRhs);
  /** Get rowRhs */
  inline const double *getRowRhs() const {return rowRhs;};

  /** Set  rowActivity to point on a vector holding the activity of the 
      constraints (i.e. coefficient matrix times separateThis). */
  virtual void setRowActivity(const double *givenRowActivity);
  /** Get rowActivity */
  inline const double *getRowActivity() const {return rowActivity;};

  /** Set colType to point on a vector holding the type of the
      variables ('B', 'I', or 'C' for Binary, Integer and Continuous) */ 
  virtual void setColType(const char *givenColType);
  /** Get colType */
  inline const char *getColType() const {return colType;};

  /** Set separateThis to point on a vector holding the point to separate */
  virtual void setSeparateThis(const double *givenSeparateThis);
  /** Get separateThis */
  inline const double *getSeparateThis() const {return separateThis;};

  /** Set doNotSeparateThis to point on a vector holding a point that
      should not be cut; only for debug */
  virtual void setDoNotSeparateThis(
                        const double *givenDoNotSeparateThis);
  /** Get doNotSeparateThis */
  inline const double *getDoNotSeparateThis() const 
                                     {return doNotSeparateThis;};

  /** Set treeInfo to point on a CglTreeInfo object */
  virtual void setTreeInfo(const CglTreeInfo *givenTreeInfo);
  /** Get treeInfo*/
  inline const CglTreeInfo *getTreeInfo() const {return treeInfo;};
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglData(const int &givenNrow = 0, const int &givenNcol = 0,
	  const CoinPackedMatrix *givenMatrixByCol = NULL,
	  const CoinPackedMatrix *givenMatrixByRow = NULL,
	  const double *givenObj = NULL,
	  const double *givenColLower = NULL, 
	  const double *givenColUpper = NULL,   
	  const double *givenRowLower = NULL, 
	  const double *givenRowUpper = NULL,
	  const double *givenRowRhs = NULL,
	  const double *givenRowActivity = NULL,
	  const char *givenColType = NULL,
	  const double *givenSeparateThis = NULL,
	  const CglTreeInfo *givenTreeInfo = NULL,
	  const double *givenDoNotSeparateThis = NULL);

  /// Copy constructor
  CglData(const CglData&);

  /// Clone
  virtual CglData* clone() const;

  /// Assignment operator 
  virtual CglData& operator=(const CglData &rhs);

  /// Destructor 
  virtual ~CglData();
  //@}

protected:

  // Private member data

  /**@name Private member data */

  //@{

  // Number of constraints
  int nrow;

  // Number of variables.
  int ncol;

  // Pointer on matrix of coefficients (ordered by columns).
  CoinPackedMatrix const *matrixByCol;

  // Pointer on matrix of coefficients (ordered by rows).
  CoinPackedMatrix const *matrixByRow;

  // Pointer on vector of objective coefficients. 
  const double *obj;

  // Pointer on vector of lower bounds on variables.
  const double *colLower; 

  // Pointer on vector of upper bounds for variables.
  const double *colUpper;   

  // Pointer on vector of lower bounds for constraints.
  const double *rowLower; 

  // Pointer on vector of upper bounds for constraints.
  const double *rowUpper;

  // Pointer on vector of upper bounds for constraints.
  const double *rowRhs;

  // Pointer on vector of activity of constraints (i.e. coefficient matrix 
  // times separateThis)..
  const double *rowActivity;

  /** Pointer on vector of characters for columns types.
      colType[i] can have values
      <UL>
      <LI> 'C' : continuous
      <LI> 'B' : binary
      <LI> 'I' : integer
      </UL>
  */
  const char *colType;

  // Pointer on vector for point to separate.  
  const double *separateThis;

  // Pointer on tree information.
  const CglTreeInfo *treeInfo;

  /// Pointer on vector for point that should not be cut; only for debug. 
  const double *doNotSeparateThis;

  //@}
};

#endif
