#ifndef ODD_HOLE_SEP_H
#define ODD_HOLE_SEP_H

#include "cgraph.h"
#include <vector>

typedef struct _OddHoleSep OddHoleSep;
typedef OddHoleSep *OddHoleSepPtr;

/**
 * Initialize Odd-Holes
 * search object.
 **/
OddHoleSep *oddhs_create(const CGraph *cgraph);

/**
 * finds odd holes that correspond to violated cuts.
 * now we are ignoring the ones that do not generate
 * violated cuts. Odd holes with size smaller or equal than 3
 * also are ignored.
 **/
size_t oddhs_search_odd_holes(OddHoleSep *oddhs, const double x[], const double rc[]);

/**
 * Returns the indexes of the idx-th discovered odd hole
 * in the last search.
 * Indexes are related to the original indexes of variables.
 **/
const std::vector<size_t>&oddhs_get_odd_hole(const OddHoleSep *oddhs, const size_t idx);
const size_t oddhs_get_odd_hole_size(const OddHoleSep *oddhs, const size_t idx);
/**
 * computes the violation of a odd hole
 **/
const double oddhs_get_odd_hole_violation(const OddHoleSep *oddhs, const size_t idx);
/**
 * returns the rhs of a oddhs constraint of size s
 * (s/2)
 **/
const double oddhs_get_odd_hole_rhs(const OddHoleSep *oddhs, const size_t idx);

/**
 * Returns the total number of discovered odd holes in
 * the last search.
 **/
const size_t oddhs_get_odd_hole_count(const OddHoleSep *oddhs);

/**
 * the inequality for a discovered odd hole
 * may be extended with the addition of
 * wheel centers - this function returns the
 * number of computed wheel centers for a
 * discovered odd hole
 */
size_t oddhs_get_nwc_doh(OddHoleSep *oddhs, const size_t doh);


/**
 * the inequality for a discovered odd hole
 * may be extended with the addition of
 * wheel centers - this function returns the
 * computed wheel centers for a discovered odd hole
 */
const std::vector<size_t>& oddhs_get_wc_doh(OddHoleSep *oddhs, const size_t doh);

/**
 * Frees memory related to Odd-Hole
 * search object.
 **/
void oddhs_free(OddHoleSepPtr *oddhs);

#endif
