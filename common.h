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

constexpr std::size_t TRACK_LENGTH = 4;
using emission_probs_t = std::vector<double>;
const std::array<char, TRACK_LENGTH> TRACK = {'A','G','C','T'};
using seq_map_t = std::unordered_map<std::string, FastaSequence>;

struct StateInfo {
    std::string name;
    emission_probs_t emission_probs;
};

struct HiddenMarkovModel {
    std::vector<StateInfo> states;
    std::vector<std::unordered_map<std::size_t, double>> transitions;
};

struct RunConfig {
    const std::vector<FastaSequence> j_repo;
    const std::vector<FastaSequence> d_repo;
    const seq_map_t v_repo;
};

struct BlastResult {
    const std::string v_name;
    const std::string v_string;
    const std::size_t v_match_start;
};

struct SequenceInfo {
    const std::string seq_name;
    const std::string seq_string;

    const BlastResult blast_result;
    const double a_score;
};

struct MutationProbabilites {
    
};

HiddenMarkovModel buildModel(const RunConfig &config, const SequenceInfo &input);

#endif //_IH_GENERATE_MODELS_COMMON_H_
