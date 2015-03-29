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
    std::string full_v_seq = config.v_repo.at(input.blast_result.v_name).c_str();
//    std::string full_v_seq = input.blast_result.v_string;
    auto t = transform(index(full_v_seq, 0), [&config](auto i) {
        return create_v_state(config, i.index(), i.value());
    });

    ret.states.assign(t.begin(), t.end());

    return ret;
}
