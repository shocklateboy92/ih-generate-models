//
// Created by Lasath on 3/24/2015.
//

#include "common.h"

#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/assign.hpp>
#include <string>
#include <regex>

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

std::string get_penta_nucleotide(std::string seq_string, int nucl_pos) {
    std::stringstream ss;
    ss << "uu" << seq_string << "uu";
    return ss.str().substr(nucl_pos,  5);
}

//new coverage values are based on the germline frequency of the hotspots
//for 59982 4mers in all germline IGHV, 2284 RGYW, 2116 WRCY, 4349 WAN and 51233 Non-HS
double RGYW_Nucleotide_coverage = 0.0381D;
double WRCY_Nucleotide_coverage = 0.0353D;
double WAN_Nucleotide_coverage = 0.0725D;
double NO_HOTSPOT_Nucleotide_coverage = 0.854D;

//mutation probabilities calculated from Barington 2007 JI paper
//4592 muts total, 598 RGYW, 633 WAN, 558 WRCY, 2803 in non-hotspots
double RGYW_mutation_prob = 0.130D;
double WRCY_mutation_prob = 0.122D;
double WAN_mutation_prob = 0.138D;
double NO_HOTSPOT_mutation_prob = 0.610D;

double NO_HOTSPOT_MUTABILITY_SCORE = (NO_HOTSPOT_mutation_prob/NO_HOTSPOT_Nucleotide_coverage);
double WAN_MUTABILITY_SCORE = (WAN_mutation_prob/WAN_Nucleotide_coverage);
double NAN_MUTABILITY_SCORE = ((WAN_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGYW_MUTABILITY_SCORE = (RGYW_mutation_prob / RGYW_Nucleotide_coverage);
double NGYW_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGYN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGNN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE*0.25) + (NO_HOTSPOT_MUTABILITY_SCORE*0.75));
double WRCY_MUTABILITY_SCORE = (WRCY_mutation_prob/WRCY_Nucleotide_coverage);
double NRCY_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double NNCY_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE * 0.25)+(NO_HOTSPOT_MUTABILITY_SCORE*0.75);
double WRCN_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);

std::vector<std::pair<std::regex, double>> hotspots = {
    {"?[at]a??", WAN_MUTABILITY_SCORE},
//    {"?na??", NAN_MUTABILITY_SCORE},

//    {"??a??", NO_HOTSPOT_MUTABILITY_SCORE},

//    {"?[ag]?[ct][at]", RGYW_MUTABILITY_SCORE},
//    {"?n?yw", NGYW_MUTABILITY_SCORE},
//    {"?r?n?", RGNN_MUTABILITY_SCORE},
//    {"?r?yn", RGYN_MUTABILITY_SCORE}
};

StateInfo create_v_state(const RunConfig &config, std::size_t index, char v) {
    std::regex rgx("?[at]a??");

    // for V gene
    // where does this constant come from? Ask Bruno
    double exp_decay_prob = std::exp(-0.0023999999999999998L * index);

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
        std::string penta_nucleotide = get_penta_nucleotide(i.index(), fstr);

        return create_v_state(config, i.index(), i.value());
    });

    ret.states.assign(t.begin(), t.end());

    return ret;
}
