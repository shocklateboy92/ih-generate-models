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

//namespace args = std::placeholders;

//double lookup_tri_nt(RunConfig config, seq_t tri_nt, char to_nt) {
//    std::cerr << tri_nt;
//    std::replace(tri_nt.begin(), tri_nt.end(), nt::any, to_nt);
//    std::cerr << " looking up " << tri_nt << std::endl;
//    return config.mutation_probs.at(tri_nt);
//}

// compare to MutationSpectrum.getTNProbability()
double fetch_mutation_ratio(RunConfig config, std::string sequence,
                            std::size_t index, char to_nt) {
    auto tri_nt = get_tri_nucleotide(sequence, index);
    auto mutated_tri_nt = tri_nt.insert(tri_nt_pos::to, 1, to_nt);

    // when we support gaps, this should cause probability = 0
    assert(to_nt != nt::gap);

    assert(sequence.find(nt::gap) == std::string::npos); // don't support gaps yet

    // if the to_nt is nt::any, we need to look up all permutations
    //    std::vector<double> probs(TRACK.size());
    //    if (tri_nt.find(nt::any) != std::string::npos) {
    //        probs.reserve(TRACK.size());
    //        std::copy_if(TRACK.begin(), TRACK.end(), probs.begin(),
    //                     std::bind(&lookup_tri_nt,
    //                               config,
    //                               mutated_tri_nt,
    //                               args::_1));
    //    }

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
