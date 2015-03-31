//
// Created by Lasath on 3/24/2015.
//

#include "common.h"

#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/assign.hpp>
#include <string>

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
        std::stoul(hsps.child_value("Hsp_hit-from"))
    };
}

double calculateAScore(const BlastResult &br) {
    return 7;
}

StateInfo create_v_state(const RunConfig &config, std::size_t index, char v) {
    auto probs = transform(TRACK, [v](char c) -> double {
        return c == v ? 1 : 0;
    });
    return {
        "V-" + std::to_string(index),
        emission_probs_t(probs.begin(), probs.end())
    };
}

HiddenMarkovModel buildModel(const RunConfig &config, const SequenceInfo &input) {
    HiddenMarkovModel ret = {};
    // Use full V_gene from the repertoire, rather
    // than just the aligned segment BLAST spits out
    std::string full_v_seq = config.v_repo.at(input.blast_result.v_name).c_str();
    // But we only care about the V sequence from the start of the aligned region
    std::string fstr = full_v_seq.substr(input.blast_result.v_match_start -1, full_v_seq.length());
    // Now make a state for each NT in the V-Gene
    auto t = transform(index(fstr, 0), [&config](auto i) {
        return create_v_state(config, i.index(), i.value());
    });

    ret.states.assign(t.begin(), t.end());

    return ret;
}
