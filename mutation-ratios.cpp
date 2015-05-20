#include "common.h"
#include "mutation-ratios.h"

namespace nucleotides {
    char gap = '-';
    char any = '*';
}

namespace nt = nucleotides;

// compare to MutationSpectrum.getTNProbability()
double fetch_mutation_ratio(RunConfig config, std::string sequence,
                            std::size_t index, char to_nt) {
    enum tri_nt_pos {
        left = 0,
        from,
        to,
        right
    };
    auto tri_nt = get_tri_nucleotide(sequence, index);
    auto mutated_tri_nt = tri_nt.insert(2, 1, to_nt);

    // when we support gaps, this should cause probability = 0
    assert(to_nt != nt::gap);

    assert(sequence.find(nt::gap) == std::string::npos); // don't support gaps yet

    std::vector<std::string> keys = {mutated_tri_nt};
    if (to_nt == nt::any) {
        std::transform(TRACK.begin(), TRACK.end(), std::back_inserter(keys),
                       [=](char t) -> std::string {
//            return tri_nt.replace();
        });
    }

    return 1;
}

mutation_probs_t readMutationProbsFile(const char *fileName) {
    mutation_probs_t ret;

    std::ifstream is(fileName, std::ios::in);
    assert (is.is_open());

    while (is.good()) {
        std::string str; double prob;
        is >> str; is >> prob;

        std::transform(str.begin(), str.end(), str.begin(), toupper);

        // hopefully, the don't change the whitespace
        ret[{str[0], str[2], str[5], str[7]}] = prob;
}

return ret;
}
