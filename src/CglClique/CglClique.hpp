// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef _CglClique_h_
#define _CglClique_h_

#include "CglCutGenerator.hpp"

/*! \brief Classical cliques

  We're looking for violated classical clique constraints on binary variables:
  \verbatim
  SUM{j} x(j) = 1 or SUM{j} x(j) <= 1.
  \endverbatim
  The general approach is to identify binary variables with fractional values
  in the current solution, form an adjacency graph on these variables, and
  look for cliques in the adjacency graph. Two variables are adjacent if both
  occur together in at least one suitable constraint. A suitable constraint is
  defined as
  - The constraint is an `=' or `<=' constraint with row upper bound of 1.0.
  - All fractional binary variables have coefficients of 1.0.
  - All other coefficients in the row are positive.

  For an explanation of the algorithms, see, for example, Hoffman, K. and
  Padberg, M., "Solving Airline Crew Scheduling Problems by Branch-and-Cut",
  Management Science 39(6), June, 1993.
*/
class CglClique : public CglCutGenerator {

  friend void CglCliqueUnitTest(const OsiSolverInterface *siP,
				const std::string mpdDir) ;

public:

  /*! \name Cut generation */
  //@{
  /*! \brief Generate cuts from cliques

    Generate cliques as requested and process them into cuts.
  */
  virtual void
  generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
	       const CglTreeInfo info = CglTreeInfo()) const ;
  //@}

  /*! \brief Possible choices for selecting the next node in the star
	     clique search.

    Minimum degree is the classic choice. The default is maximum value, then
    maximum degree.
  */
  enum scl_next_node_method {
      SCL_MIN_DEGREE,		///< minimum degree
      SCL_MAX_DEGREE,		///< maximum degree
      SCL_MAX_XJ_MAX_DEG	///< maximum value, then maximum degree
  } ;

  /*! \name Cut generation control */
  //@{

  /// Set method for selecting next node in star clique search
  inline void setStarCliqueNextNodeMethod(scl_next_node_method method)
  { scl_next_node_rule = method ; }

  /*! \brief The maximal length of the candidate list to extend a row clique.
  
    See #scl_candidate_length_threshold.
  */
  inline void setStarCliqueCandidateLengthThreshold(int maxlen)
  { scl_candidate_length_threshold = maxlen ; }

  /*! \brief The maximal length of the candidate list to extend a star clique.
  
    See #rcl_candidate_length_threshold.
  */
  inline void setRowCliqueCandidateLengthThreshold(int maxlen)
  { rcl_candidate_length_threshold = maxlen ; }

  /// Control star clique formation
  inline void setDoStarClique(bool yesno = true) { do_star_clique = yesno ; }
  /// Control row clique formation
  inline void setDoRowClique(bool yesno = true) { do_row_clique = yesno ; }

  /*! \brief Set minimum acceptable violation of clique cuts

    The default value is -1.0. If #petol = -1.0, whether by default or because
    it's been set to -1.0, the actual value used will be the solver's primal
    tolerance.
  */
  inline void setMinViolation(double minviol) { petol = minviol; }
  /// Get minimum acceptable violation of clique cuts
  inline double getMinViolation() const { return petol; }
  //@}

  /*! \name Constructors and destructors */
  //@{
  /*! \brief Default constructor.

    If \p setPacking is set to true then CglClique will assume that the
    problem describes a set packing problem, i.e.,
      - all variables are binary
      - the matrix is a 0-1 matrix
      - all constraints are '= 1' or '<= 1'
     
    Otherwise, CglClique will start the #generateCuts method by scanning
    the matrix for suitable rows.

    \p justOriginalRows can be set to true to limit the constraints that are
    considered when in the search tree; see #justOriginalRows_.
  */
  CglClique(bool setPacking = false, bool justOriginalRows = false) ;
  /// Copy constructor
  CglClique(const CglClique& rhs) ;
  /// Clone
  virtual CglCutGenerator * clone() const ;
  /// Assignment operator
  CglClique& operator=(const CglClique& rhs) ;
 
  /// Destructor
  virtual ~CglClique() {}
  //@}

  /*! \name Utility methods */
  //@{
  /// Control printing of detailed statistics on the star clique method
  inline void setStarCliqueReport(bool yesno = true)
  { scl_report_result = yesno ; }
  /// Control printing of detailed statistics on the row clique method
  inline void setRowCliqueReport(bool yesno = true)
  { rcl_report_result = yesno ; }

  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp) ;
  //@}

private:

    struct frac_graph ;
    friend struct frac_graph ;

    /*! \brief A node of the fractional graph.
    
      There is a node for every variable at fractional level.
    */
    struct fnode {
      /// Start of neighbours of this node in frac_graph::all_nbr
      int *nbrs ;
      /*! \brief 1-x_i-x_j, needed for odd holes
      
        Currently unused.
      */
      double *edgecosts ;
      /// degree of the node
      int degree ;
      /// the fractional value of the variable corresponding to this node
      double val ;
    } ;

    /*! \brief A graph corresponding to a fractional solution of an LP.
    
      Two nodes are adjacent iff their columns are non-orthogonal.
    */
    struct frac_graph {
      /*! \brief The number of nodes in the graph
      
	The number of fractional values in the LP solution
      */
      int nodenum ;
      /// The number of edges in the graph
      int edgenum ;
      /*! \brief density
      
        (actual edges)/(potential edges) = edgenum/(nodenum choose 2)
      */
      double density ;
      /// Node with minimum degree
      int min_deg_node ;
      /// Minimum degree
      int min_degree ;
      /// Node with maximum degree
      int max_deg_node ;
      /// Maximum degree
      int max_degree ;
      /*! \brief The array of the nodes in the graph

        nodes[i] points to the start of the neighbours of i in #all_nbr
      */
      fnode *nodes ;
      /*! \brief The array of all the neighbours.
      
        First the indices of the nodes adjacent to node 0 are listed, then
	those adjacent to node 1, etc.
      */
      int *all_nbr ;
      /*! \brief The array of the costs of the edges going to the neighbors
      
        Currently unused.
      */
      double *all_edgecost ;
      /// Constructor
      frac_graph() :
	  nodenum(0), edgenum(0), density(0),
	  min_deg_node(0), min_degree(0), max_deg_node(0), max_degree(0),
	  nodes(0), all_nbr(0), all_edgecost(0) {}
    } ;

protected:

    /*! \brief True if the whole matrix is a set packing problem. */
    bool setPacking_;
    /*! \brief True to consider only original rows.

      Note that this is ineffective at the root, where all rows will always be
      considered.
    */
    bool justOriginalRows_;

    /// Number of clique rows under consideration
    mutable int sp_numrows;
    /// Original row index for clique rows
    mutable int* sp_orig_row_ind;
    /// Number of fractional binary variables
    mutable int sp_numcols;
    /// Original column index for fractional binary variables
    mutable int* sp_orig_col_ind;
    /*! \brief Solution value for fractional binary variables

      Correlated with #sp_orig_col_ind.
    */
    mutable double* sp_colsol;
    /// Column starts for the set packing submatrix
    mutable int* sp_col_start;
    /// Row indices for columns of the set packing submatrix
    mutable int* sp_col_ind;
    /// Row starts for the set packing submatrix
    mutable int* sp_row_start;
    /// Column indices for the rows of the set packing submatrix
    mutable int* sp_row_ind;

    /// The intersection graph corresponding to the set packing problem
    mutable frac_graph fgraph;
    /*! \brief The node-node incidence matrix of the intersection graph.
    
      This is stored as a full sp_numcols x sp_numcols matrix.
    */
    mutable bool* node_node;

    /*! \brief Minimum acceptable clique violation
    
      The default set by the constructor is -1.0. If petol is -1.0 when
      #generateCuts is called (whether set by default or by the user), the
      actual value used will be the solver's primal tolerance.
    */
    mutable double petol;

    /*! \name Clique control parameters */
    /**@{*/
    /** True to generate row cliques. */
    bool do_row_clique;
    /** True to generate star cliques. */
    bool do_star_clique;

    /** How the next node to be added to the star clique should be selected */
    scl_next_node_method scl_next_node_rule;
    /*! \brief Maximum number of candidates for complete enumeration of
    	       star cliques.

      In the star clique method, the maximal length of the candidate list
      (those nodes that are in a star, i.e., connected to the center of the
      star) to allow complete enumeration of all maximal cliques. If the
      candidate list is longer, a greedy algorithm is used.
    */
    int scl_candidate_length_threshold;

    /** True to report detailed statistics on the star clique method */
    bool scl_report_result;

    /*! \brief Maximum number of candidates for complete enumeration of
	       row cliques.
     
      In the row clique method, the maximal length of the candidate list
      (those nodes that can extend the row clique, i.e., connected to all
      nodes in the row clique) to allow complete enumeration of all maximal
      cliques. If the candidate list is longer, a greedy algorithm is used.
    */
    int rcl_candidate_length_threshold;

    /** True to report detailed statistics on the row clique method */
    bool rcl_report_result;
    /**@}*/

    /*! \name Variables/arrays that are used across many methods */
    /**@{*/
    /** List of indices that must be in the to-be-created clique. This is just
	a pointer, it is never new'd and therefore does not need to be
	delete[]'d either. */
    mutable const int* cl_perm_indices;
    /** The length of cl_perm_indices */
    mutable int cl_perm_length;

    /** List of indices that should be considered for extending the ones listed
	in cl_perm_indices. */
    mutable int* cl_indices;
    /** The length of cl_indices */
    mutable int cl_length;

    /** An array of nodes discarded from the candidate list. These are
	rechecked when a maximal clique is found just to make sure that the
	clique is really maximal. */
    mutable int* cl_del_indices;
    /** The length of cl_del_indices */
    mutable int cl_del_length;

    /**@}*/

private:
    /*! \brief Scan through the variables and select those that are binary
    	       and have a fractional value.

      The definition of fractional is asymmetric: x* > tol1 and x* < tol2,
      where tol1 is the solver's primal tolerance and tol2 is the value set
      for #petol. It's possible this is an error; there's some special-case
      code in this method that does not make sense in context.
    */
    void selectFractionalBinaries(const OsiSolverInterface &si) const;

    /*! \brief Scan through the variables and select those that have a
    	       fractional value.

      This method assumes that all the variables are binary, hence it does not
      check this. The definition of fractional is symmetric, x* > tol1 and
      x* < tol2, tol1 = tol2 = solver's primal tolerance.
    */
    void selectFractionals(const OsiSolverInterface& si) const;
    /**  */
    void selectRowCliques(const OsiSolverInterface& si,int numOriginalRows) const;
    /**  */
    void createSetPackingSubMatrix(const OsiSolverInterface& si) const;
    /**  */
    void createFractionalGraph() const;
    /**  */
    int createNodeNode() const;
    /**  */
    void deleteSetPackingSubMatrix() const;
    /**  */
    void deleteFractionalGraph() const;
    /**  */
    void find_scl(OsiCuts& cs) const;
    /**  */
    void find_rcl(OsiCuts& cs) const;
    /**  */
    int scl_choose_next_node(const int current_nodenum,
			     const int *current_indices,
			     const int *current_degrees,
			     const double *current_values) const;
    /**  */
    void scl_delete_node(const int del_ind, int& current_nodenum,
			 int *current_indices, int *current_degrees,
			 double *current_values) const;
    /**  */
    int enumerate_maximal_cliques(int& pos, bool* scl_label, OsiCuts& cs) const;
    /**  */
    int greedy_maximal_clique(OsiCuts& cs) const;
    /**  */
    void recordClique(const int len, int* indices, OsiCuts& cs) const;
};
//#############################################################################
/** A function that tests the methods in the CglClique class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CglCliqueUnitTest(const OsiSolverInterface * siP,
		       const std::string mpdDir);

class CglProbing;
/*! This works on a fake solver i.e. invented rows

  In more words, you can load in an Osi (#fakeSolver_) of your choosing.
  When the #generateCuts method is called, the loaded Osi is used, if
  present. If it's not present, then the si supplied as a parameter to
  #generatCuts will be used. If it, too, is absent, things will end badly.

  CglFakeClique also conceals a CglProbing object, which is invoked iff the
  loaded Osi is used.
*/
class CglFakeClique : public CglClique {
  
public:
  /// Copy constructor
  CglFakeClique(const CglFakeClique& rhs);
  /// Clone
  virtual CglCutGenerator * clone() const;
  
  /// Assignment operator
  CglFakeClique& operator=(const CglFakeClique& rhs);
  
  virtual void
  generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
	       const CglTreeInfo info = CglTreeInfo()) const;
  
  /**@name Constructors and destructors */
  //@{
  /** Default constructor.
      If the setPacking argument is set to true then CglFakeClique will
      assume that the
      problem in the solverinterface passed to the generateCuts() method
      describes a set packing problem, i.e.,
      - all variables are binary
      - the matrix is a 0-1 matrix
      - all constraints are '= 1' or '<= 1'
      
      Otherwise, CglFakeClique will
      start the generateCuts() methods by scanning the matrix for them.
  */
  CglFakeClique(OsiSolverInterface * solver=NULL,bool setPacking = false);
  /// Destructor
  virtual ~CglFakeClique();
  /// Assign solver (generator takes over ownership)
  void assignSolver(OsiSolverInterface * fakeSolver);
protected:
  /// fake solver to use
  mutable OsiSolverInterface * fakeSolver_;
  /// Probing object
  mutable CglProbing * probing_;
};

#endif
