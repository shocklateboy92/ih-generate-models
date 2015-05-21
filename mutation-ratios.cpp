#include "common.h"
#include "mutation-ratios.h"
#include <functional>

namespace nucleotides {
    char gap = '-';
    char any = '*';
}

namespace nt = nucleotides;
enum tri_nt_pos {
    left = 0,
    from,
    to,
    right
};

// compare to MutationSpectrum.getTNProbability()
double fetch_mutation_ratio(RunConfig config, std::string sequence,
                            std::size_t index, char to_nt) {
    auto tri_nt = get_tri_nucleotide(sequence, index);
    auto mutated_tri_nt = tri_nt.insert(tri_nt_pos::to, 1, to_nt);

    // when we support gaps, this should cause probability = 0
    assert(to_nt != nt::gap);

    assert(sequence.find(nt::gap) == std::string::npos); // don't support gaps yet

    double sum = 0;
    std::size_t matches = 0;
    for (char c : TRACK) {
        seq_t key(mutated_tri_nt);
        std::replace(key.begin(), key.end(), nt::any, c);

        if (key[tri_nt_pos::from] != key[tri_nt_pos::to]) {
            sum += config.mutation_probs.at(key);
            matches++;
        }
    }

    return sum / static_cast<double>(matches);
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
