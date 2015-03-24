//
// Created by Lasath on 3/24/2015.
//

#ifndef _IH_GENERATE_MODELS_COMMON_H_
#define _IH_GENERATE_MODELS_COMMON_H_

#include <unordered_map>
#include <vector>

struct StateInfo {

};

struct HiddenMarkovModel {
    std::vector<StateInfo> states;
//    std::unordered_map<std::pair<std::size_t, std::size_t>, double> transitions;
};

HiddenMarkovModel buildModel(void);

#endif //_IH_GENERATE_MODELS_COMMON_H_
