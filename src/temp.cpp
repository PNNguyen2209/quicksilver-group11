std::shared_ptr<SimpleGraph> layeredTC(std::shared_ptr<SimpleGraph> &base) {
    auto out = std::make_shared<SimpleGraph>(base->getNoVertices());
    out->setNoLabels(1);

    uint32_t V = base->getNoVertices();
    bool updated = true;
    int maxIterations = V; // Safety limit on layers

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
    std::unordered_set<uint32_t> updatedNodes;

    // Step 1: Copy direct edges from base to out
    for (uint32_t u = 0; u < V; ++u) {
        for (auto &edge : base->adj[u]) {
            edges.emplace_back(u, edge.second, 0);
            updatedNodes.insert(u);
        }
    }
    out->addEdgesInBatch(edges);

    // Clear edges for reuse
    edges.clear();
    int iteration = 0;

    // Step 2: Iteratively compute transitive closure
    while (updated && iteration++ < maxIterations) {
        updated = false;
        std::unordered_set<uint32_t> newUpdatedNodes;

        for (uint32_t u : updatedNodes) {
            std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> newEdges;

            // Traverse edges from current node in `out`
            for (auto &edge : out->adj[u]) {
                uint32_t mid = edge.second;

                // Traverse edges from `mid` in `base`
                for (auto &next : base->adj[mid]) {
                    uint32_t w = next.second;

                    // Use `edgeExists(source, target, label)` to check edge existence
                    if (!out->edgeExists(u, w, 0)) {
                        newEdges.emplace_back(u, w, 0);
                        updated = true;
                        newUpdatedNodes.insert(u);
                    }
                }
            }

            // Add all new edges for this node in one batch
            edges.insert(edges.end(), newEdges.begin(), newEdges.end());
        }

        // Batch add new edges to maintain sorted adjacency lists
        out->addEdgesInBatch(edges);
        edges.clear();

        // Update the list of nodes to process in the next iteration
        updatedNodes = std::move(newUpdatedNodes);
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
SimpleEvaluator::transitiveClosure(std::shared_ptr<SimpleGraph> &base, bool hasReverseDirection, cardStat stats, size_t pathSize) {

    // Case 1: Sparse Graph
    auto condition = 0.1 * base->getNoVertices();
    if (stats.noPaths < condition) {
        return optimizedBFS(base);
    }
    // if (pathSize > 2) {
    //     return optimizedBFS(base);
    // }

    if (hasReverseDirection) {
        return optimizedBFS(base);
    }

    // std::cout << "Using Layered TC" << std::endl;
    // Default to Layered TC if no specific case matches
    return layeredTC(base);
}