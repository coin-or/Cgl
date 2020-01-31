#ifndef SPATHS_H
#define SPATHS_H

typedef struct _ShortestPathsFinder ShortestPathsFinder;
typedef ShortestPathsFinder *ShortestPathsFinderPtr;

/*
 * Creates Shortest Path Finder.
 * arcStart[i] indicates the position in vector
 * toNode and dist where arcs of node i start
*/
ShortestPathsFinderPtr spf_create(const size_t nodes, const size_t arcs, const size_t *arcStart, const size_t *toNode, const size_t *dist);

/* updates just one arc
 **/
void spf_update_arc(ShortestPathsFinder *spf, const size_t tail, const size_t head, const size_t cost);

/* returns the distance of one arc
 **/
size_t spf_get_arc(ShortestPathsFinder *spf, const size_t tail, const size_t head);

/*
 * queries number of nodes
 */
size_t spf_nodes(ShortestPathsFinder *spf);

/*
 * queries number of arcs
 */
size_t spf_arcs(ShortestPathsFinder *spf);

/*
 * executes the shortest path finder
 * using the Dijkstra algorithm
 */
void spf_find(ShortestPathsFinder *spf, const size_t origin);
void spf_find(ShortestPathsFinder *spf, const size_t origin, const size_t destination);

/*
 * solution query: returns distance to a node after executing spf_find
 */
size_t spf_get_dist(const ShortestPathsFinderPtr spf, const size_t node);

/*
 * solution query: returns previous node and allows one to build a path after executing spf_find
 */
size_t spf_get_previous(const ShortestPathsFinderPtr spf, const size_t node);

size_t *spf_previous(const ShortestPathsFinder *spf);

/* returns all previous nodes
 * which should be steped
 * to arrive at a given node (this node is not included)
 * returns how many nodes were filles in indexes
 */
size_t spf_get_path(const ShortestPathsFinder *spf, const size_t toNode, size_t indexes[]);

/*
 * releases Shortest Path Finder object
 */
void spf_free(ShortestPathsFinderPtr *spf);

#endif /* ifndef SPATHS_H */
