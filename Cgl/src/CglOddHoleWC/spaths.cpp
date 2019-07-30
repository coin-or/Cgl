#include <cstdio>
#include <cstring>
#include <cassert>
#include <limits>
#include <algorithm>
#include "spaths.h"
#include "node_heap.h"

#define SP_INFTY_DIST std::numeric_limits<size_t>::max()
#define NULL_NODE     std::numeric_limits<size_t>::max()

typedef struct {
    size_t node;
    size_t distance;
} Neighbor;

struct _ShortestPathsFinder {
    size_t nodes;
    size_t arcs;

    Neighbor *neighs; // all neighbors
    Neighbor **startn; // Start of neighbors for node i. The neighbor ends at startn[i+1]

    // solution:
    size_t *dist;
    size_t *previous;
    size_t *path; // temporary storage for path

    NodeHeap *nh;
};

int compNeighs(const void *n1, const void *n2) {
    const Neighbor *pn1 = (const Neighbor *) n1;
    const Neighbor *pn2 = (const Neighbor *) n2;

    return pn1->node - pn2->node;
}

bool compare_neighbors(const Neighbor &n1, const Neighbor &n2) {
    return n1.node < n2.node;
}

ShortestPathsFinderPtr spf_create(const size_t nodes, const size_t arcs, const size_t *arcStart, const size_t *toNode, const size_t *dist) {
    ShortestPathsFinder *spf = new ShortestPathsFinder;

    spf->nodes = nodes;
    spf->arcs = arcs;
    spf->startn = new Neighbor*[nodes + 1];
    spf->nh = nh_create(nodes, SP_INFTY_DIST);
    spf->previous = new size_t[nodes];
    spf->dist = new size_t[nodes];
    spf->path = new size_t[nodes];
    spf->neighs = new Neighbor[arcs];

    for (size_t idx = 0; idx < arcs; idx++) {
        spf->neighs[idx].node = toNode[idx];
        spf->neighs[idx].distance = dist[idx];
#ifdef DEBUG
        assert(spf->neighs[idx].node < spf->nodes);
#endif
    }

    for (size_t n = 0; n <= spf->nodes; n++) {
        spf->startn[n] = spf->neighs + arcStart[n];
    }

    return spf;
}

void spf_find(ShortestPathsFinder *spf, const size_t origin) {
    NodeHeap *nh = spf->nh;

    nh_reset(nh);

    for (size_t i = 0; i < spf->nodes; i++) {
        spf->dist[i] = SP_INFTY_DIST;
        spf->previous[i] = NULL_NODE;
    }

    spf->dist[origin] = 0;
    nh_update(nh, origin, 0);

    size_t topCost, topNode;
    while ((topCost = nh_remove_first(nh, &topNode)) < SP_INFTY_DIST) {
        // updating neighbors distances by iterating in all neighbors
        for (Neighbor *n = spf->startn[topNode]; n < spf->startn[topNode+1]; n++) {
            const size_t toNode = n->node;
            const size_t dist = n->distance;
            const size_t newDist = topCost + dist;

            if (spf->dist[toNode] > newDist) {
                spf->previous[toNode] = topNode;
                spf->dist[toNode] = newDist;
                nh_update(spf->nh, toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

void spf_find(ShortestPathsFinder *spf, const size_t origin, const size_t destination) {
    NodeHeap *nh = spf->nh;

    nh_reset(nh);

    for (size_t i = 0; i < spf->nodes; i++) {
        spf->dist[i] = SP_INFTY_DIST;
        spf->previous[i] = NULL_NODE;
    }

    spf->dist[origin] = 0;
    nh_update(nh, origin, 0);

    size_t topCost, topNode;
    while ((topCost = nh_remove_first(nh, &topNode)) < SP_INFTY_DIST) {
        if(topNode == destination) {
            break;
        }
        // updating neighbors distances by iterating in all neighbors
        for (Neighbor *n = spf->startn[topNode]; n < spf->startn[topNode+1]; n++) {
            const size_t toNode = n->node;
            const size_t dist = n->distance;
            const size_t newDist = topCost + dist;

            if (spf->dist[toNode] > newDist) {
                spf->previous[toNode] = topNode;
                spf->dist[toNode] = newDist;
                nh_update(spf->nh, toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

size_t spf_nodes(ShortestPathsFinder *spf) {
    return spf->nodes;
}

size_t spf_arcs(ShortestPathsFinder *spf) {
    return spf->arcs;
}

size_t spf_get_dist(const ShortestPathsFinderPtr spf, const size_t node) {
#ifdef DEBUG
    assert(node < spf->nodes);
#endif
    return spf->dist[node];
}

size_t spf_get_previous(const ShortestPathsFinderPtr spf, const size_t node) {
#ifdef DEBUG
    assert(node < spf->nodes);
#endif
    return spf->previous[node];
}

size_t *spf_previous(const ShortestPathsFinder *spf) {
    return spf->previous;
}

size_t spf_get_path(const ShortestPathsFinder *spf, const size_t toNode, size_t indexes[]) {
    size_t n = 0;
    size_t currNode = toNode;

    spf->path[n++] = currNode;

    while((currNode = spf->previous[currNode]) != NULL_NODE) {
        spf->path[n++] = currNode;
    }

    for (size_t i = 0; i < n; i++) {
        indexes[i] = spf->path[n-i-1];
    }

    return n;
}

void spf_free(ShortestPathsFinderPtr *spf) {
    delete[] (*spf)->neighs;
    delete[] (*spf)->startn;
    nh_free(&((*spf)->nh));
    delete[] (*spf)->previous;
    delete[] (*spf)->dist;
    delete[] (*spf)->path;

    delete (*spf);
    (*spf) = nullptr;
}

void spf_update_arc(ShortestPathsFinder *spf, const size_t head, const size_t tail, const size_t cost) {
    const Neighbor *start = spf->startn[head];
    const Neighbor *end = spf->startn[head + 1];
    const Neighbor key = {tail, 0};
#ifdef DEBUG
    assert(end - start > 0);
#endif
    Neighbor *result = (Neighbor *) bsearch(&key, start, end - start, sizeof(Neighbor), compNeighs);

#ifdef DEBUG
    if (!result) {
        fprintf(stderr, "Could not find arc (%ld->%ld) in graph.\n", head, tail);
        abort();
    }
    assert(result->node == tail);
#endif
    result->distance = cost;
}

size_t spf_get_arc(ShortestPathsFinder *spf, const size_t head, const size_t tail) {
    const Neighbor *start = spf->startn[head];
    const Neighbor *end = spf->startn[head + 1];
    const Neighbor key = {tail, 0};
#ifdef DEBUG
    assert(end - start > 0);
#endif
    Neighbor *result = (Neighbor *) bsearch(&key, start, end - start, sizeof(Neighbor), compNeighs);

#ifdef DEBUG
    if (!result) {
        fprintf(stderr, "Could not find arc (%ld->%ld) in graph.\n", head, tail);
        abort();
    }
    assert(result->node == tail);
#endif
    return result->distance;
}
