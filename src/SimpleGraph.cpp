#include "SimpleGraph.h"
#include <iostream>

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

SimpleGraph::SimpleGraph(std::shared_ptr<SimpleGraph> &other)   {
    adj = other->adj;
    V = other->getNoVertices();
    L = other->getNoLabels();
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
    adj.resize(V);
}

bool SimpleGraph::anyEdgeExists(uint32_t source, uint32_t target) const {
    if (source >= adj.size()) return false;
    for (auto &edge : adj[source]) {
        // edge.first = label, edge.second = targetNode
        if (edge.second == target) {
            return true; // Found (source, ?label, target)
        }
    }
    return false;
}


uint32_t SimpleGraph::getNoEdges() const {
    uint32_t sum = 0;
    for (const auto & l : adj)
        sum += l.size();
    return sum;
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortPairs(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {

    uint32_t sum = 0;

    for (auto sourceVec : adj) {

        std::sort(sourceVec.begin(), sourceVec.end(), sortPairs);

        uint32_t prevTarget = 0;
        uint32_t prevLabel = 0;
        bool first = true;

        for (const auto &labelTgtPair : sourceVec) {
            if (first || !(prevTarget == labelTgtPair.second && prevLabel == labelTgtPair.first)) {
                first = false;
                sum++;
                prevTarget = labelTgtPair.second;
                prevLabel = labelTgtPair.first;
            }
        }
    }

    return sum;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
}

void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    //adj[from].emplace_back(std::make_pair(edgeLabel, to));
    adj[from].emplace_back(edgeLabel, to);
}

bool SimpleGraph::edgeExists(uint32_t source, uint32_t target, uint32_t label) const {
    if (source >= adj.size()) return false;
    for (const auto &edge : adj[source]) {
        if (edge.first == label && edge.second == target) {
            return true;
        }
    }
    return false;
}

bool SimpleGraph::isEmpty() const {
    for (const auto &edges : adj) {
        if (!edges.empty()) {
            return false; // Found an edge, graph is not empty
        }
    }
    return true; // No edges found, graph is empty
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {
    std::ifstream graphFile(fileName);
    if (!graphFile.is_open()) {
        throw std::runtime_error("Could not open file: " + fileName);
    }

    std::string headerLine;
    if (!std::getline(graphFile, headerLine)) {
        throw std::runtime_error("Empty file or read error in header");
    }

    // The header is "noNodes,noEdges,noLabels"
    uint32_t noNodes = 0, noEdges = 0, noLabels = 0;
    int matched = std::sscanf(headerLine.c_str(), "%u,%u,%u",
                              &noNodes, &noEdges, &noLabels);
    if (matched != 3) {
        throw std::runtime_error("Invalid header format in: " + headerLine);
    }

    setNoVertices(noNodes);
    setNoLabels(noLabels);

    // Approximate out-degree for each node to reduce re-allocations
    if (noNodes > 0 && noEdges > 0) {
        uint32_t approxOutDegree = (noEdges / noNodes) + 1;
        for (auto &edgeList : adj) {
            edgeList.reserve(approxOutDegree);
        }
    }

    // Read edges
    uint32_t subject, predicate, object;
    char dot;
    while (true) {
        graphFile >> subject >> predicate >> object >> dot;
        if (!graphFile.good()) {
            // Either EOF or read error => break
            break;
        }
        if (!edgeExists(subject, object, predicate)) {
            addEdge(subject, object, predicate);
        }
    }

    graphFile.close();
}