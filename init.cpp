#include "init.h"
#include "mutation-ratios.h"

#include <boost/range/numeric.hpp>


RunConfig prepareConfig(int argc, char *argv[]) {
    return {
        readRepertoire("IGHJRepertoire.fasta"),
                readRepertoire("IGHDRepertoire.fasta"),
                boost::accumulate(
                    readRepertoire("IGHVRepertoire.fasta"),
                    seq_map_t(),
                    [](seq_map_t &map, FastaSequence f) -> seq_map_t & {
            map.insert({f.name(), f});
            return map;
        }
        ),
        readMutationProbsFile("Mutation_spectrum.txt")
    };
}


std::vector<FastaSequence> readRepertoire(const char *repo_path) {
    std::vector<FastaSequence> ret;
    FastaReader fr;
    std::ifstream is(repo_path);
    assert(is.is_open());

    fr.toupper(true);

    while (!is.eof()) {
        auto seq = fr.next(is);
        ret.push_back(*seq.get());
    }

    return ret;
}
