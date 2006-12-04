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
    If a data member is NULL and the generator needs the data, an
    error is raised. If a data member is not NULL, the calling method
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

  /** Set  matrixByColPtr to point on the coefficient matrix ordered by 
      columns */
  virtual void setMatrixByColPtr(const CoinPackedMatrix **givenMatrixByColPtr);
  /** Get matrixByColPtr */
  inline const CoinPackedMatrix ** getMatrixByColPtr() const {return matrixByColPtr;};

  /** Set  matrixByRowPtr to point on the coefficient matrix ordered by 
      rows */
  virtual void setMatrixByRowPtr(const CoinPackedMatrix **givenMatrixByRowPtr);
  /** Get matrixByRowPtr */
  inline const CoinPackedMatrix ** getMatrixByRowPtr() const {return matrixByRowPtr;};

  /** Set objPtr to point on a vector holding the objective coefficient values */
  virtual void setObjPtr(const double **givenObjPtr);
  /** Get objPtr */
  inline const double ** getObjPtr() const {return objPtr;};

  /** Set  colLowerPtr to point on a vector holding the lower bounds on the 
      variables */
  virtual void setColLowerPtr(const double **givenColLowerPtr);
  /** Get colLowerPtr */
  inline const double ** getColLowerPtr() const {return colLowerPtr;};

  /** Set  colUpperPtr to point on a vector holding the upper bounds on the 
      variables */
  virtual void setColUpperPtr(const double **givenColUpperPtr);
  /** Get colUpperPtr */
  inline const double ** getColUpperPtr() const {return colUpperPtr;};

  /** Set  rowLowerPtr to point on a vector holding the lower bounds on the 
      constraints */
  virtual void setRowLowerPtr(const double **givenRowLowerPtr);
  /** Get rowLowerPtr */
  inline const double ** getRowLowerPtr() const {return rowLowerPtr;};

  /** Set  rowUpperPtr to point on a vector holding the upper bounds on the 
      constraints */
  virtual void setRowUpperPtr(const double **givenRowUpperPtr);
  /** Get rowUpperPtr */
  inline const double ** getRowUpperPtr() const {return rowUpperPtr;};

  /** Set  rowRhsPtr to point on a vector holding the right hand side of the 
      constraints (for a ranged constraint, it contains the upper bound). */
  virtual void setRowRhsPtr(const double **givenRowRhsPtr);
  /** Get rowRhsPtr */
  inline const double ** getRowRhsPtr() const {return rowRhsPtr;};

  /** Set  rowActivityPtr to point on a vector holding the activity of the 
      constraints (i.e. coefficient matrix times separateThis). */
  virtual void setRowActivityPtr(const double **givenRowActivityPtr);
  /** Get rowActivityPtr */
  inline const double ** getRowActivityPtr() const {return rowActivityPtr;};

  /** Set colTypePtr to point on a vector holding the type of the
      variables ('B', 'I', or 'C' for Binary, Integer and Continuous) */ 
  virtual void setColTypePtr(char * const *givenColTypePtr);
  /** Get colTypePtr*/
  inline char *const *getColTypePtr() const {return colTypePtr;};

  /** Set separateThisPtr to point on a vector holding the point to separate */
  virtual void setSeparateThisPtr(const double **givenSeparateThisPtr);
  /** Get separateThisPtr */
  inline const double ** getSeparateThisPtr() const {return separateThisPtr;};

  /** Set doNotSeparateThisPtr to point on a vector holding a point that
      should not be cut; only for debug */
  virtual void setDoNotSeparateThisPtr(
                        const double **givenDoNotSeparateThisPtr);
  /** Get doNotSeparateThisPtr */
  inline const double ** getDoNotSeparateThisPtr() const 
                                     {return doNotSeparateThisPtr;};

  /** Set treeInfoPtr to point on a CglTreeInfo object */
  virtual void setTreeInfoPtr(const CglTreeInfo **givenTreeInfoPtr);
  /** Get treeInfoPtr*/
  inline const CglTreeInfo ** getTreeInfoPtr() const {return treeInfoPtr;};
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglData(const int &givenNrow = 0, const int &givenNcol = 0,
	  const CoinPackedMatrix **givenMatrixByColPtr = NULL,
	  const CoinPackedMatrix **givenMatrixByRowPtr = NULL,
	  const double **givenObjPtr = NULL,
	  const double **givenColLowerPtr = NULL, 
	  const double **givenColUpperPtr = NULL,   
	  const double **givenRowLowerPtr = NULL, 
	  const double **givenRowUpperPtr = NULL,
	  const double **givenRowRhsPtr = NULL,
	  const double **givenRowActivityPtr = NULL,
	  char * const *givenColTypePtr = NULL,
	  const double **givenSeparateThisPtr = NULL,
	  const CglTreeInfo **givenTreeInfoPtr = NULL,
	  const double **givendoNotSeparateThisPtr = NULL);

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
  CoinPackedMatrix const **matrixByColPtr;

  // Pointer on matrix of coefficients (ordered by rows).
  CoinPackedMatrix const **matrixByRowPtr;

  // Pointer on vector of objective coefficients. 
  const double **objPtr;

  // Pointer on vector of lower bounds on variables.
  const double **colLowerPtr; 

  // Pointer on vector of upper bounds for variables.
  const double **colUpperPtr;   

  // Pointer on vector of lower bounds for constraints.
  const double **rowLowerPtr; 

  // Pointer on vector of upper bounds for constraints.
  const double **rowUpperPtr;

  // Pointer on vector of upper bounds for constraints.
  const double **rowRhsPtr;

  // Pointer on vector of activity of constraints (i.e. coefficient matrix 
  // times separateThis)..
  const double **rowActivityPtr;

  /** Pointer on vector of characters for columns types.
      colType[i] can have values
      <UL>
      <LI> 'C' : continuous
      <LI> 'B' : binary
      <LI> 'I' : integer
      </UL>
  */
  char * const *colTypePtr;

  // Pointer on vector for point to separate.  
  const double **separateThisPtr;

  // Pointer on tree information.
  const CglTreeInfo **treeInfoPtr;

  /// Pointer on vector for point that should not be cut; only for debug. 
  const double **doNotSeparateThisPtr;

  //@}
};

#endif
