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

struct StateInfo {
};

struct HiddenMarkovModel {
    std::vector<StateInfo> states;
    std::vector<std::unordered_map<std::size_t, double>> transitions;
};

struct RunConfig {
    std::vector<FastaSequence> j_repo;
    std::vector<FastaSequence> d_repo;
};

struct BlastResult {
    std::string v_name;
    std::string v_string;
    std::size_t v_match_start_index;
    std::string input_sequence;
};

struct SequenceInfo {
    std::string seq_name;
    std::string seq_string;

    std::unique_ptr<BlastResult> blast_result;
    double a_score;
};

HiddenMarkovModel buildModel(RunConfig config, BlastResult input);

#endif //_IH_GENERATE_MODELS_COMMON_H_
