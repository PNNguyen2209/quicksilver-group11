#include "SimpleGraph.h"
#include "SimpleEstimator.h"

#include <numeric>
SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){
    graph = g;
}

void SimpleEstimator::prepare() {

    labelFrequencies.resize(graph->getNoLabels(), 0);
    nodeDegrees.resize(graph->getNoVertices());

    reverseIndex.clear();
    reverseIndex.resize(graph->getNoVertices());

    for (uint32_t source = 0; source < graph->getNoVertices(); ++source) {
        // Out-degree is simply how many edges in graph->adj[source]
        nodeDegrees[source].outDegree = static_cast<uint32_t>(graph->adj[source].size());

        // Count label frequencies, build reverseIndex
        for (const auto &labelTgt : graph->adj[source]) {
            uint32_t label  = labelTgt.first;
            uint32_t target = labelTgt.second;
            labelFrequencies[label]++;

            // Populate reverseIndex for target
            reverseIndex[target].push_back({label, source});
        }
    }

    for (uint32_t v = 0; v < graph->getNoVertices(); ++v) {
        nodeDegrees[v].inDegree = static_cast<uint32_t>(reverseIndex[v].size());
    }

    uint32_t maxOutDegree = 0;
    uint32_t maxInDegree = 0;

    uint32_t totalOutDegrees = 0;
    uint32_t totalInDegrees = 0;

    std::vector<uint32_t> outDegrees(graph->getNoVertices(), 0);
    std::vector<uint32_t> inDegrees(graph->getNoVertices(), 0);

    for (uint32_t source = 0; source < graph->getNoVertices(); ++source) {
        uint32_t outDegree = graph->adj[source].size();
        outDegrees[source] = outDegree;

        totalOutDegrees += outDegree;

        if (outDegree > maxOutDegree) maxOutDegree = outDegree;

        for (const auto& labelTarget : graph->adj[source]) {
            inDegrees[labelTarget.second]++;
        }
    }

    for (uint32_t degree : inDegrees) {
        totalInDegrees += degree;
        if (degree > maxInDegree) maxInDegree = degree;
    }

    outDegreeHistogram.resize(maxOutDegree + 1, 0);
    inDegreeHistogram.resize(maxInDegree + 1, 0);

    for (uint32_t degree : outDegrees) {
        outDegreeHistogram[degree]++;
    }
    for (uint32_t degree : inDegrees) {
        inDegreeHistogram[degree]++;
    }

    avgOutDeg = static_cast<float>(totalOutDegrees) / static_cast<float>(graph->getNoVertices());
    avgInDeg = static_cast<float>(totalInDegrees) / static_cast<float>(graph->getNoVertices());
}

cardStat SimpleEstimator::estimate(Triple &q) {

    // cardStat result {0, 0, 0};

    return estimateConcat(q.src, q.path, q.trg);
}

cardStat SimpleEstimator::estimateLabelDir(Identifier &src, LabelDir &labelDir, Identifier &trg) {
    cardStat result {0, 0, 0};

    if (src == NO_IDENTIFIER && trg == NO_IDENTIFIER) {
        result.noPaths = labelFrequencies[labelDir.label];

        uint32_t sumLabelFreq = std::accumulate(labelFrequencies.begin(), labelFrequencies.end(), 0u);
        float labelProb = static_cast<float>(labelFrequencies[labelDir.label]) / static_cast<float>(sumLabelFreq);

        for (uint32_t i = 1; i < outDegreeHistogram.size(); ++i) {
            result.noOut += static_cast<uint32_t>(static_cast<float>(outDegreeHistogram[i] * i) * labelProb);
        }

        for (uint32_t i = 1; i < inDegreeHistogram.size(); ++i) {
            result.noIn += static_cast<uint32_t>(static_cast<float>(inDegreeHistogram[i] * i) * labelProb);
        }

        if (labelDir.reverse) {
            std::swap(result.noIn, result.noOut);
        }

    } else if ((trg == NO_IDENTIFIER && !labelDir.reverse) || (src == NO_IDENTIFIER && labelDir.reverse)) {
        auto sourceNode = labelDir.reverse ? trg : src;
        auto basePaths = nodeDegrees[sourceNode].outDegree;
        result.noPaths = basePaths;

        for (const auto& neighbor : graph->adj[src]) {
            if (neighbor.first == labelDir.label) {
                result.noPaths++;
            }
        }

        result.noPaths *= avgOutDeg;
        result.noIn = result.noPaths;
        result.noOut = result.noPaths > 0 ? 1 : 0;

        if (labelDir.reverse) {
            std::swap(result.noIn, result.noOut);
        }

    } else if ((src == NO_IDENTIFIER && !labelDir.reverse) || (trg == NO_IDENTIFIER && labelDir.reverse)) {
        auto node = labelDir.reverse ? src : trg;
        auto basePaths = nodeDegrees[node].inDegree;
        result.noPaths = basePaths;

        for (const auto &incoming : reverseIndex[node]) {
            if (incoming.first == labelDir.label) {
                result.noPaths++;
            }
        }

        result.noPaths *= avgInDeg;
        result.noOut = result.noPaths;
        result.noIn = result.noPaths > 0 ? 1 : 0;

        if (labelDir.reverse) {
            std::swap(result.noIn, result.noOut);
        }

    } else {
        if (labelDir.reverse) {
            std::swap(result.noIn, result.noOut);
        }

        for (const auto& neighbor : graph->adj[src]) {
            if (neighbor.first == labelDir.label && neighbor.second == trg) {
                result.noOut = 1;
            }
        }

        result.noIn = result.noOut;
        result.noPaths = result.noOut;
    }

    return result;
}

cardStat SimpleEstimator::estimateUnion(Identifier &src, std::vector<LabelDir> &labelDirs, Identifier &trg) {
    cardStat result {0, 0, 0};

    for (LabelDir &labelDir : labelDirs) {
        cardStat labelDirEst = estimateLabelDir(src, labelDir, trg);

        if (src == NO_IDENTIFIER)
            result.noOut += labelDirEst.noOut;

        result.noPaths += labelDirEst.noPaths;

        if (trg == NO_IDENTIFIER)
            result.noIn += labelDirEst.noIn;
    }

    if (src != NO_IDENTIFIER)
        result.noOut = 1;

    if (trg != NO_IDENTIFIER)
        result.noIn = 1;

    return result;
}

cardStat SimpleEstimator::estimateKleene(Identifier &src, PathEntry &pathEntry, Identifier &trg) {
    cardStat result = estimateUnion(src, pathEntry.labels, trg);

    if (pathEntry.kleene) {
        float discountFactor = 0.5f;    // Try 0.6 or 0.7
        result.noPaths = static_cast<uint32_t>(static_cast<float>(result.noPaths) / (1.0f - discountFactor));

        result.noIn = static_cast<uint32_t>(static_cast<float>(result.noIn) * avgOutDeg);
    }

    return result;
}

cardStat SimpleEstimator::estimateConcat(Identifier &src, std::vector<PathEntry> &pathEntries, Identifier &trg) {
    // Define a heuristic for sorting PathEntry
    auto heuristic = [&](const PathEntry &a, const PathEntry &b) {
        cardStat estA = estimateKleene(src, const_cast<PathEntry &>(a), trg);
        cardStat estB = estimateKleene(src, const_cast<PathEntry &>(b), trg);

        // Weighted combination of metrics
        float scoreA = estA.noPaths * 0.5f + estA.noIn * 0.25f + estA.noOut * 0.25f;
        float scoreB = estB.noPaths * 0.5f + estB.noIn * 0.25f + estB.noOut * 0.25f;

        return scoreA < scoreB; // Lower score is better
    };

    // Sort pathEntries using the heuristic
    std::sort(pathEntries.begin(), pathEntries.end(), heuristic);

    // Perform concatenation in the sorted order
    cardStat result = estimateKleene(src, pathEntries[0], trg);
    for (uint32_t i = 1; i < pathEntries.size(); ++i) {
        cardStat nextEstimate = estimateKleene(src, pathEntries[i], trg);

        if (src == NO_IDENTIFIER && trg == NO_IDENTIFIER) {
            uint32_t minNoIn = std::min(result.noIn, nextEstimate.noIn);
            uint32_t minNoOut = std::max(result.noOut, nextEstimate.noOut);
            float est1 = static_cast<float>(result.noPaths * nextEstimate.noPaths) / static_cast<float>(minNoIn);
            float est2 = static_cast<float>(result.noPaths * nextEstimate.noPaths) / static_cast<float>(minNoOut);
            result.noPaths = static_cast<uint32_t>(std::min(est1, est2));

            result.noIn = nextEstimate.noIn;
        } else {
            result.noPaths = result.noPaths + nextEstimate.noPaths;
            if (src != NO_IDENTIFIER)
                result.noIn = result.noPaths;
            if (trg != NO_IDENTIFIER)
                result.noOut = result.noPaths;
        }
    }

    return result;
}
