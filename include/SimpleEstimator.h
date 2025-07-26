#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

struct NodeDegree {
    uint32_t inDegree = 0;  // Number of incoming edges
    uint32_t outDegree = 0; // Number of outgoing edges
};


class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    // My fields
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> reverseIndex;

    struct NodeDegree {
        uint32_t inDegree = 0;
        uint32_t outDegree = 0;
    };

    std::vector<uint32_t> labelPaths;

    std::vector<NodeDegree> nodeDegrees;
    std::vector<uint32_t> labelFrequencies;
    std::vector<uint32_t> outDegreeHistogram;
    std::vector<uint32_t> inDegreeHistogram;
    std::vector<std::vector<uint32_t>> coOccurrenceMatrix;

    float avgOutDeg;
    float avgInDeg;


public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    uint32_t estimate(PathEntry &entry);
    cardStat estimate(Triple &q) override ;

    cardStat estimateLabelDir(Identifier &src, LabelDir &ld, Identifier &trg);
    cardStat estimateUnion(Identifier &src, std::vector<LabelDir> &labelDirs, Identifier &trg);
    cardStat estimateKleene(Identifier &src, PathEntry &pathEntry, Identifier &trg);
    cardStat estimateConcat(Identifier &src, std::vector<PathEntry> &pathEntries, Identifier &trg);

};


#endif //QS_SIMPLEESTIMATOR_H
