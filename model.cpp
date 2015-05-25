//
// Created by Lasath on 3/24/2015.
//

#include "common.h"
#include "mutation-ratios.h"
#include "mutability.h"

#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/assign.hpp>

using namespace boost::adaptors;
using namespace boost::assign;

BlastResult parseBlastOutput(const pugi::xpath_node &node) {
    // Only considering the most likely V-gene for now
    auto hit = node.node().child("Iteration_hits").child("Hit");
    auto hsps = hit.child("Hit_hsps").child("Hsp");
    return {
        // The current version of BLAST seems to prepend "lcl|" to Id
        std::string(hit.child_value("Hit_id")).substr(4),
        hsps.child_value("Hsp_hseq"),
        std::stoul(hsps.child_value("Hsp_hit-from")),
        hsps.child_value("Hsp_qseq")
    };
}

decltype(get_penta_nucleotide) get_penta_nucleotide = _get_n_nucleotide<5>;
decltype(get_tri_nucleotide) get_tri_nucleotide = _get_n_nucleotide<3>;

template <std::size_t N>
std::string _get_n_nucleotide(std::string seq_string, int nucl_pos) {
    static const std::size_t padding = 2;
    std::stringstream ss;

    ss << "**" << seq_string << "**";
    return ss.str().substr(nucl_pos + padding - N/2,  N);
}

double MIN_MUTATION_PROB = 0.02l;
double EXP_DECAY_INDEX_CONST = -0.0023999999999999998L;

StateInfo createState(const RunConfig &config, const SequenceInfo &input,
                      double exp_decay_prob, seq_t fstr, std::size_t i) {

        double mutability_score = fetch_mutability_score(fstr, i);

        double mutation_prob =
                exp_decay_prob * mutability_score * input.a_score;

        auto probs = transform(TRACK, [&](char c) -> double {

            double mutation_ratio =
                    fetch_mutation_ratio(config, fstr, i, c);

            return c == fstr[i]
                    ? 1 - ((mutation_prob - MIN_MUTATION_PROB) * 3)
                    : mutation_prob * mutation_ratio + MIN_MUTATION_PROB;
        });

        return {
            "V-" + std::to_string(i),
                    emission_probs_t(probs.begin(), probs.end())
        };

}

HiddenMarkovModel buildModel(const RunConfig &config, const SequenceInfo &input) {
    HiddenMarkovModel ret = {};

    // Use full V_gene from the repertoire, rather
    // than just the aligned segment BLAST spits out
    std::string full_v_seq = config.v_repo
            .at(input.blast_result.v_name)
            .c_str();
    // But we only care about the V sequence from the start of the aligned region
    std::string fstr = full_v_seq.substr(
                input.blast_result.v_match_start -1,
                full_v_seq.length());

    // Now make a state for each NT in the V-Gene
    double i = 0;
    for (nt_t nt : fstr) {
        ret.states.push_back(
                    createState(
                        config,
                        input,
                        std::exp(EXP_DECAY_INDEX_CONST * i),
                        fstr,
                        i)
                    );
        i++;
    }

    // Now, do the same for all possible D-genes
    for (auto seq : config.d_repo) {
        double i = 0;
        for (nt_t nt : seq_t(seq.c_str())) {
            ret.states.push_back(
                        createState(
                            config,
                            input,
                            std::exp(EXP_DECAY_INDEX_CONST * i),
                            seq.c_str(),
                            i)
                        );
            i++;
        }
    }

    return ret;
}
