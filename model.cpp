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
    auto t = transform(index(fstr, 0), [&](auto i) -> StateInfo {

        // for V gene
        // where does this constant come from? Ask Bruno
        double exp_decay_prob = std::exp(-0.0023999999999999998L * i.value());

        double mutability_score = fetch_mutability_score(fstr, i.index());

        double mutation_prob =
                exp_decay_prob * mutability_score * input.a_score;

        auto probs = transform(TRACK, [&](char c) -> double {

            double mutation_ratio =
                    fetch_mutation_ratio(config, fstr, i.index(), c);

            return c == i.value()
                    ? 1 - ((mutation_prob - MIN_MUTATION_PROB) * 3)
                    : mutation_prob * mutation_ratio + MIN_MUTATION_PROB;
        });

        return {
            "V-" + std::to_string(i.index()),
                    emission_probs_t(probs.begin(), probs.end())
        };
    });

    ret.states.assign(t.begin(), t.end());

    return ret;
}
