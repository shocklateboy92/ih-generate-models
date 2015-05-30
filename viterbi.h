#ifndef VITERBI_H
#define VITERBI_H

#include "common.h"

std::vector<StateInfo> do_viterbi(const HiddenMarkovModel &model);

#endif // VITERBI_H
