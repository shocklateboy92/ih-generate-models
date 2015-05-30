#include "viterbi.h"

#include <boost/multi_array.hpp>
#include "utils.hpp"

typedef std::vector<double> vec;
typedef boost::multi_array<double, 2> mat;
typedef std::vector<std::unordered_map<std::size_t, double>> s_mat;
typedef std::pair<std::size_t, double> s_val;

static const double SCALE = 1;

std::vector<StateInfo> do_viterbi(const HiddenMarkovModel &model) {
    std::size_t num_observations = model.input.seq_string.size();
    const auto &transitions = model.transitions;
    const auto &emissions = model.states;

    vec observations;
    for (char c : model.input.seq_string) {
        observations.push_back(std::find(TRACK.begin(), TRACK.end(), c) - TRACK.begin());
    }

    s_mat V(num_observations);

    // set state 0 as the initial state
    V[0][0] = 1.0l;

    // do viterbi iterative viterbi algorithm
    for (std::size_t t = 1; t < num_observations; t++) {
        for (s_val v : V[t-1]) {
            for (s_val transition : transitions[v.first]) {
                double v_next = v.second *
                        transition.second *
                        SCALE *
                        emissions[transition.first].emission_probs[observations[t-1]];
                if (V[t][transition.first] <= v_next) {
                    V[t][transition.first] = v_next;
                }
            }
        }
    }

    std::vector<StateInfo> most_likely_path;
    most_likely_path.reserve(num_observations);
    for (std::unordered_map<std::size_t, double> v : V) {
        auto max = *std::max_element(v.begin(), v.end(),
                                    [](s_val a, s_val b){
                                        return a.second < b.second;
                                    });
        most_likely_path.push_back(model.states[max.first]);
    }

    return most_likely_path;
}
