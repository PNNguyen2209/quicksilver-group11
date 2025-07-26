#include "MyQueryExtensions.h"

std::string MyQueryExtensions::labelDirToString(const LabelDir& labelDir) {
    std::ostringstream out;
    labelDir.appendTo(out);
    return out.str();
}

std::string MyQueryExtensions::pathEntryToString(const PathEntry& pathEntry, bool noKleene) {
    std::ostringstream out;
    bool parens = (!noKleene && pathEntry.kleene) || pathEntry.labels.size() > 1;
    if (parens) {
        out << '(';
    }

    for (std::size_t i = 0, max = pathEntry.labels.size(); i < max; i++) {
        if (i != 0) {
            out << '|';
        }
        pathEntry.labels[i].appendTo(out);
    }

    if (parens) {
        out << ')';
    }
    if (pathEntry.kleene && !noKleene) {
        out << '+';
    }
    return out.str();
}

std::string MyQueryExtensions::tripleToString(const Triple& triple, bool noEndpoints) {
    std::ostringstream out;
    if (!noEndpoints) {
        if (triple.src == NO_IDENTIFIER) {
            out << '*';
        } else {
            out << triple.src;
        }
        out << ',';
    }

    for (std::size_t i = 0, max = triple.path.size(); i < max; i++) {
        if (i != 0) {
            out << '/';
        }
        out << pathEntryToString(triple.path[i]);
    }

    if (!noEndpoints) {
        out << ',';
        if (triple.trg == NO_IDENTIFIER) {
            out << '*';
        } else {
            out << triple.trg;
        }
    }
    return out.str();
}
