//
// Created by Lasath on 3/24/2015.
//

#ifndef _IH_GENERATE_MODELS_COMMON_H_
#define _IH_GENERATE_MODELS_COMMON_H_

#include <unordered_map>
#include <vector>
#include <boost/functional/hash.hpp>

#include <pugixml/pugixml.hpp>
#include <fastareader/fastareader.h>
#include <catch/catch.hpp>

#include "stdafx.h"

const std::array<char, 4> TRACK = {'A','G','C','T'};

using nt_t = char;
using seq_t = std::basic_string<nt_t>;
using seq_map_t = std::unordered_map<std::string, FastaSequence>;
using emission_probs_t = std::vector<double>;
using transitions_t = std::vector<std::unordered_map<std::size_t, double>>;
using range_t = std::pair<std::size_t, std::size_t>;
using probs_list_t = std::vector<double>;

using mutation_probs_t = std::unordered_map<
    std::string,
    double >;

struct StateInfo {
    std::string name;
    emission_probs_t emission_probs;
};

struct HiddenMarkovModel {
    std::vector<StateInfo> states;
    transitions_t transitions;
};

struct RunConfig {
    const std::vector<FastaSequence> j_repo;
    const std::vector<FastaSequence> d_repo;
    const seq_map_t v_repo;
    const mutation_probs_t mutation_probs;
};

struct BlastResult {
    const std::string v_name;
    const std::string v_match_string;
    const std::size_t v_match_start;
    const std::string v_input_string;
};

struct SequenceInfo {
    const std::string seq_name;
    const std::string seq_string;

    const BlastResult blast_result;
    const double a_score;
};

struct MutationProbabilites {

};

extern std::function<seq_t(seq_t, int)> get_penta_nucleotide;
extern std::function<seq_t(seq_t, int)> get_tri_nucleotide;

HiddenMarkovModel buildModel(const RunConfig &config, const SequenceInfo &input);
BlastResult parseBlastOutput(const pugi::xpath_node &node);

#endif //_IH_GENERATE_MODELS_COMMON_H_
