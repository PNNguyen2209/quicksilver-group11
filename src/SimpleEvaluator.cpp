#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include "SimpleGraph.h"
#include <queue>
#include <functional>

/**
 * Build a reverse index from a graph's forward adjacency lists.
 * reverseIndex[v] will contain (label, u) pairs for edges u -> v with that label.
 */
std::vector<std::vector<std::pair<uint32_t, uint32_t>>>
buildReverseIndex(const std::shared_ptr<SimpleGraph> &g)
{
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> reverseIndex(g->getNoVertices());
    for (uint32_t s = 0; s < g->getNoVertices(); ++s) {
        for (auto &lblTgt : g->adj[s]) {
            // lblTgt.first = label, lblTgt.second = target
            reverseIndex[lblTgt.second].push_back({lblTgt.first, s});
        }
    }
    return reverseIndex;
}

SimpleEvaluator::SimpleEvaluator(std::shared_ptr<SimpleGraph> &g) {

    // works only with SimpleGraph
    graph = g;
    est = nullptr; // estimator not attached by default
}

void SimpleEvaluator::attachEstimator(std::shared_ptr<SimpleEstimator> &e) {
    est = e;
}

void SimpleEvaluator::prepare() {

    // if attached, prepare the estimator
    if(est != nullptr) est->prepare();

    // precompute all labels
    for(uint32_t label = 0; label < graph->getNoLabels(); label++) {
        auto labelGraph = selectLabel(label, 0, false, graph);
        LabelDir labelDir {label, false};
        std::string keyFull = MyQueryExtensions::labelDirToString(labelDir);
        cache[keyFull] = labelGraph;
    }

    // prepare other things here.., if necessary

}

cardStat SimpleEvaluator::computeStats(std::shared_ptr<SimpleGraph> &g) {

    cardStat stats {};

    for(int source = 0; source < g->getNoVertices(); source++) {
        if(!g->adj[source].empty()) stats.noOut++;
    }

    stats.noPaths = g->getNoDistinctEdges();

    auto rIndex = buildReverseIndex(g);
    for (int v = 0; v < g->getNoVertices(); v++) {
        if (!rIndex[v].empty()) {
            stats.noIn++;
        }
    }

    return stats;
}

/**
 * Select all edges from a graph with a given edge label.
 * @param projectLabel Label to select.
 * @param outLabel Label to rename the selected edge labels to (used in the TC).
 * @param inverse Follow the edges in inverse direction.
 * @param in The graph to select from.
 * @return The graph which only has edges with the specified edge label.
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::selectLabel(uint32_t projectLabel, uint32_t outLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if(!inverse) {
        // going forward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, outLabel);
            }
        }
    } else {
        // going backward
        auto rIndex = buildReverseIndex(in);
        for (uint32_t node = 0; node < in->getNoVertices(); node++) {
            for (auto &lblSrc : rIndex[node]) {
                // lblSrc.first = label, lblSrc.second = "true source"
                if (lblSrc.first == projectLabel) {
                    out->addEdge(node, lblSrc.second, outLabel);
                }
            }
        }
    }

    return out;
}


std::shared_ptr<SimpleGraph> layeredTC(std::shared_ptr<SimpleGraph> &base) {
    auto out = std::make_shared<SimpleGraph>(base->getNoVertices());
    out->setNoLabels(1);

    uint32_t V = base->getNoVertices();
    bool updated = true;
    int maxIterations = V; // Maximum layers for safety

    // Copy direct edges from base to out
    for (uint32_t u = 0; u < V; ++u) {
        for (auto &edge : base->adj[u]) {
            out->addEdge(u, edge.second, 0);
        }
    }

    int iteration = 0;
    while (updated && iteration++ < maxIterations) {
        updated = false;

        for (uint32_t u = 0; u < V; ++u) {
            std::vector<uint32_t> newEdges;

            for (auto &edge : out->adj[u]) {
                uint32_t mid = edge.second;

                for (auto &next : base->adj[mid]) {
                    uint32_t w = next.second;

                    if (!out->edgeExists(u, w, 0)) {
                        newEdges.push_back(w);
                        updated = true;
                    }
                }
            }

            for (uint32_t w : newEdges) {
                out->addEdge(u, w, 0);
            }
        }
    }
    return out;
}

std::shared_ptr<SimpleGraph> optimizedBFS(std::shared_ptr<SimpleGraph> &base) {
    auto out = std::make_shared<SimpleGraph>(base->getNoVertices());
    out->setNoLabels(1);

    uint32_t V = base->getNoVertices();
    std::vector<uint32_t> visited(V, 0);
    uint32_t visitId = 0;

    for (uint32_t s = 0; s < V; ++s) {
        ++visitId;  // Use a unique ID to mark visited nodes
        std::queue<uint32_t> q;
        q.push(s);

        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();

            if (visited[u] != visitId) {
                visited[u] = visitId;
                if (s != u) {
                    out->addEdge(s, u, 0);
                }
                for (auto &e : base->adj[u]) {
                    uint32_t w = e.second;
                    if (visited[w] != visitId) {
                        q.push(w);
                    }
                }
            }
        }
    }

    return out;
}

/**
 * A naive transitive closure (TC) computation.
 * @param in Input graph
 * @return
 */
std::shared_ptr<SimpleGraph>
SimpleEvaluator::transitiveClosure(std::shared_ptr<SimpleGraph> &base, bool hasReverseDirection, cardStat stats) {

    // Case 1: Sparse Graph
    if (stats.noPaths < 0.1 * base->getNoVertices() * base->getNoVertices()) {
        return optimizedBFS(base);
    }

    // Case 2: High Out-Degree Nodes
    if (stats.noOut > 0.8 * base->getNoVertices()) {
        return optimizedBFS(base);
    }

    // Case 3: High-Diameter Graph
    if (stats.noIn < 0.5 * base->getNoVertices()) {
        return optimizedBFS(base);
    }

    if (hasReverseDirection) {
        return optimizedBFS(base);
    }

    // Default to Layered TC if no specific case matches
    return layeredTC(base);
}


/**
 * Merges a graph into another graph.
 * @param left A graph to be merged into.
 * @param right A graph to be merged from.
 * @return A number of distinct new edges added from the "right" graph into the "left" graph.
 */
uint32_t SimpleEvaluator::unionDistinct(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    uint32_t numNewAdded = 0;

    for(uint32_t source = 0; source < right->getNoVertices(); source++) {
        for (auto labelTarget : right->adj[source]) {

            auto label = labelTarget.first;
            auto target = labelTarget.second;

            if(!left->edgeExists(source, target, label)) {
                left->addEdge(source, target, label);
                numNewAdded++;
            }
        }
    }

    return numNewAdded;
}

/**
 * Simple implementation of a join of two graphs.
 * @param left A graph to be joined.
 * @param right Another graph to join with.
 * @return Answer graph for a join. Note that all labels in the answer graph are "0".
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {
    auto out = std::make_shared<SimpleGraph>(left->getNoVertices());
    out->setNoLabels(1);

    for (uint32_t leftSource = 0; leftSource < left->getNoVertices(); ++leftSource) {
        for (const auto &leftEdge : left->adj[leftSource]) {
            uint32_t intermediate = leftEdge.second; // Intermediate node
            for (const auto &rightEdge : right->adj[intermediate]) {
                uint32_t rightTarget = rightEdge.second;
                out->addEdge(leftSource, rightTarget, 0); // Add combined path
            }
        }
    }

    return out;
}

bool hasReverseDirection(const PathEntry &pe) {
    for (const auto& labelDir : pe.labels) {
        if (labelDir.reverse) {
            return true;
        }
    }
    return false;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluateUnionKleene(PathEntry &pe) {
    if(pe.kleene) {
        std::string keyFull = MyQueryExtensions::pathEntryToString(pe);

        if (cache.find(keyFull) != cache.end()) {
            return cache[keyFull];
        }

        pe.kleene = false;
        auto base = evaluateUnionKleene(pe);
        Identifier src = NO_IDENTIFIER;
        Identifier trg = NO_IDENTIFIER;
        auto estimate = est->estimateKleene(src, pe, trg);
        auto transitiveClosureGraph = transitiveClosure(base, hasReverseDirection(pe), estimate);
        cache[keyFull] = transitiveClosureGraph;

        return transitiveClosureGraph;
    } else {
        // not Kleene
        std::string keyNoKleene = MyQueryExtensions::pathEntryToString(pe);
        if (cache.find(keyNoKleene) != cache.end()) {
            return cache[keyNoKleene];
        }

        std::shared_ptr<SimpleGraph> out;

        if (pe.labels.size() == 1) {
            // Single labelDir selection
            std::string keySingleLabel = MyQueryExtensions::labelDirToString(pe.labels[0]);

            if (cache.find(keySingleLabel) != cache.end()) {
                // **Clone the cached graph instead of using it directly**
                out = std::make_shared<SimpleGraph>(*cache[keySingleLabel]);
            } else {
                out = selectLabel(pe.labels[0].label, 0, pe.labels[0].reverse, graph);
                cache[keySingleLabel] = out;
            }
        } else {
            // (left-deep) union
            std::string keyFirstLabel = MyQueryExtensions::labelDirToString(pe.labels[0]);

            if (cache.find(keyFirstLabel) != cache.end()) {
                // **Clone the cached graph instead of using it directly**
                out = std::make_shared<SimpleGraph>(*cache[keyFirstLabel]);
            } else {
                out = selectLabel(pe.labels[0].label, 0, pe.labels[0].reverse, graph);
                cache[keyFirstLabel] = out;
            }

            // Union with subsequent labels
            for (size_t i = 1; i < pe.labels.size(); ++i) {
                std::shared_ptr<SimpleGraph> rg;
                std::string keyIthLabel = MyQueryExtensions::labelDirToString(pe.labels[i]);

                // Retrieve or compute the i-th label
                if (cache.find(keyIthLabel) != cache.end()) {
                    rg = cache[keyIthLabel];
                } else {
                    rg = selectLabel(pe.labels[i].label, 0, pe.labels[i].reverse, graph);
                    cache[keyIthLabel] = rg;
                }

                // **Ensure that `rg` is not modified if it's from the cache**
                // Since `unionDistinct` only modifies `out`, and `rg` is not modified, it's safe.
                unionDistinct(out, rg);
            }
        }

        cache[keyNoKleene] = out;
        return out;
    }
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluateConcat(std::vector<PathEntry> &path) {
    if (path.empty()) {
        // No path entries, return an empty graph
        return std::make_shared<SimpleGraph>(graph->getNoVertices());
    }

    size_t n = path.size();
    // DP tables
    std::vector<std::vector<size_t>> plan(n, std::vector<size_t>(n, 0)); // Best split
    std::vector<std::vector<cardStat>> dpEstimates(n, std::vector<cardStat>(n)); // Best estimates

    Identifier src = NO_IDENTIFIER;
    Identifier trg = NO_IDENTIFIER;
    // Step 1: Initialize DP table for single segments
    for (size_t i = 0; i < n; ++i) {
        dpEstimates[i][i] = est->estimateKleene(src, path[i], trg);
        plan[i][i] = i; // Single segment, no split
    }

    // Step 2: Fill DP table for subpaths of increasing length
    for (size_t len = 2; len <= n; ++len) {
        for (size_t i = 0; i <= n - len; ++i) {
            size_t j = i + len - 1;

            // Initialize with a large cost
            cardStat bestEstimate = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
            size_t bestSplit = i;

            // Consider all possible splits
            for (size_t k = i; k < j; ++k) {
                auto leftEst = dpEstimates[i][k];
                auto rightEst = dpEstimates[k + 1][j];

                // Estimate cost for forward join
                uint32_t estimateCost = leftEst.noPaths * rightEst.noIn;

                if (estimateCost < bestEstimate.noPaths) {
                    bestEstimate = {leftEst.noOut, rightEst.noIn, estimateCost};
                    bestSplit = k;
                }
            }

            // Store the best split and estimates
            dpEstimates[i][j] = bestEstimate;
            plan[i][j] = bestSplit;
        }
    }

    // Helper function to print the plan
    auto printPlan = [&](size_t i, size_t j, const auto &rec) -> void {
        if (i == j) {
            std::cout << "Path segment [" << i << "]\n";
            return;
        }

        size_t split = plan[i][j];
        std::cout << "Join segments [" << i << "-" << split << "] and [" << (split + 1) << "-" << j << "] FORWARD\n";
        rec(i, split, rec);
        rec(split + 1, j, rec);
    };

    // Print the best plan for debugging
    // std::cout << "Best join plan:\n";
    // printPlan(0, n - 1, printPlan);

    // Step 3: Perform the joins for the best plan
    std::function<std::shared_ptr<SimpleGraph>(size_t, size_t)> evaluateBestPlan =
        [&](size_t i, size_t j) -> std::shared_ptr<SimpleGraph> {
            if (i == j) {
                return evaluateUnionKleene(path[i]);
            }

            size_t split = plan[i][j];
            auto leftGraph = evaluateBestPlan(i, split);
            auto rightGraph = evaluateBestPlan(split + 1, j);

            // Perform forward join only
            return join(leftGraph, rightGraph);
        };

    // Evaluate the best plan for the entire path
    return evaluateBestPlan(0, n - 1);
}



uint32_t SimpleEstimator::estimate(PathEntry &entry) {
    // Use label frequencies or other statistics for estimation
    uint32_t estimate = 0;
    for (auto &labelDir : entry.labels) {
        estimate += labelFrequencies[labelDir.label]; // Sum up frequencies
    }
    return estimate;
}

/**
 * Perform a selection on a source constant.
 * @param s A source constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<SimpleGraph> selectSource(Identifier &s, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    for (auto labelTarget : in->adj[s]) {

        auto label = labelTarget.first;
        auto target = labelTarget.second;

        out->addEdge(s, target, label);
    }

    return out;
}

/**
 * Perform a selection on a target constant.
 * @param s A target constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<SimpleGraph> selectTarget(Identifier &t, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    auto rIndex = buildReverseIndex(in);
    for (auto &lblSrc : rIndex[t]) {
        auto label  = lblSrc.first;
        auto source = lblSrc.second;
        out->addEdge(source, t, label);
    }
    return out;
}

/**
 * Evaluate a path query. Produce a cardinality of the answer graph.
 * @param query Query to evaluate.
 * @return A cardinality statistics of the answer graph.
 */
cardStat SimpleEvaluator::evaluate(Triple &query) {
    // Try retrieving full query from cache
    // std::string keyFull = MyQueryExtensions::tripleToString(query);
    // if (cache.find(keyFull) != cache.end()) {
    //     auto cachedResult = cache[keyFull];
    //     return computeStats(cachedResult);
    // }

    // Try retrieving query without endpoints from cache
    std::shared_ptr<SimpleGraph> res;
    std::string keyNoEndpoints = MyQueryExtensions::tripleToString(query, true);
    if (cache.find(keyNoEndpoints) != cache.end()) {
        res = cache[keyNoEndpoints];
    } else {
        res = evaluateConcat(query.path);
        cache[keyNoEndpoints] = res;
    }

    if(query.src != NO_IDENTIFIER) res = selectSource(query.src, res);
    else if(query.trg != NO_IDENTIFIER) res = selectTarget(query.trg, res);
    // cache[keyFull] = res;

    return computeStats(res);

}