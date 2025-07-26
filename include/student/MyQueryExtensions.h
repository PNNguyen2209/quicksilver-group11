#ifndef MYQUERYEXTENSIONS_H
#define MYQUERYEXTENSIONS_H

#include "Query.h" // Include the baseline Query.h
#include <sstream>
#include <vector>
#include <string>

class MyQueryExtensions {
public:
    // Extended functionality for LabelDir
    static std::string labelDirToString(const LabelDir& labelDir);

    // Extended functionality for PathEntry
    static std::string pathEntryToString(const PathEntry& pathEntry, bool noKleene = false);

    // Extended functionality for Triple
    static std::string tripleToString(const Triple& triple, bool noEndpoints = false);
};

#endif // MYQUERYEXTENSIONS_H
