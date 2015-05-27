//
// Created by Lasath on 3/24/2015.
//

#include "common.h"
#include "mutation-ratios.h"
#include "mutability.h"

#include "utils.hpp"

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

template <std::size_t N>
std::string _get_n_nucleotide(std::string seq_string, int nucl_pos) {
    static const std::size_t padding = 2;
    std::stringstream ss;

    ss << "**" << seq_string << "**";
    return ss.str().substr(nucl_pos + padding - N/2,  N);
}

decltype(get_penta_nucleotide) get_penta_nucleotide = _get_n_nucleotide<5>;
decltype(get_tri_nucleotide) get_tri_nucleotide = _get_n_nucleotide<3>;

double MIN_MUTATION_PROB = 0.02l;
double EXP_DECAY_INDEX_CONST = -0.0023999999999999998L;

StateInfo createState(const RunConfig &config, const SequenceInfo &input,
                      double exp_decay_prob, const seq_t &fstr, std::size_t i) {

        double mutability_score = fetch_mutability_score(fstr, i);

        double mutation_prob =
                exp_decay_prob * mutability_score * input.a_score;

        return {
            "V-" + std::to_string(i),
             transform(TRACK, [&](char c) -> double {

                double mutation_ratio =
                        fetch_mutation_ratio(config, fstr, i, c);

                return c == fstr[i]
                        ? 1 - ((mutation_prob - MIN_MUTATION_PROB) * 3)
                        : mutation_prob * mutation_ratio + MIN_MUTATION_PROB;
            })
        };
}

void createStates(const RunConfig &config, const SequenceInfo &input,
                  std::function<double(double)> exp_decay_fn,
                  const seq_t &fstr, HiddenMarkovModel &ret)
{
    for (std::size_t i = 0; i < fstr.size(); i++) {
        ret.states.push_back(
                    createState(
                        config,
                        input,
                        exp_decay_fn(i),
                        fstr,
                        i)
                    );
    }
}

struct state_pos {
    enum {
        magic = 0,
        v_end_no_exo,
        v_end_exo,
        n1_start,
        n1_end,
        n2_start,
        n2_end,
        d_end_no_exo,
        d_end_exo
    };


};

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

    auto exp_decay_fn = [](double i) {return std::exp(EXP_DECAY_INDEX_CONST * i);};

    // Now make a state for each NT in the V-Gene
    createStates(config, input, exp_decay_fn, fstr, ret);

    range_t v_states = {0, ret.states.size() - 1};

    // Now, do the same for all possible D-genes
    for (auto seq : config.d_repo) {
        createStates(config, input,
                     [](double i) {
            return std::exp(EXP_DECAY_INDEX_CONST * i);
        },
        seq.c_str(), ret);
    }

    // Now, repeat for all the J genes as well
    for (auto seq : config.j_repo) {
        createStates(config, input,
                     [](double i) {
            return std::exp(EXP_DECAY_INDEX_CONST * i);
        },
        seq.c_str(), ret);
    }

    // Initialize the transitions
    ret.transitions = transitions_t(ret.states.size());

    // Using hardcoded values for now, will read from file later
    probs_list_t v_end_exo_probs = { 1.7849193927612356E-6, 4.337659583280163E-5, 6.481808081758074E-4, 0.005955808351520267, 0.03365031276420145, 0.11690722185950304, 0.24974535503543732, 0.32806306897882054, 0.2649848906871161 };

    for (std::size_t i = 0; i < v_states.second; i++) {
        assert(ret.transitions[i].empty());
        auto probs_offset = v_states.second - v_end_exo_probs.size();
        if (i > probs_offset) {
            double bail_prob = v_end_exo_probs.at(i - probs_offset - 1);
            ret.transitions[i] = {
                {
                    state_pos::v_end_exo,
                    bail_prob
                },
                {i + 1, 1 - bail_prob}
            };
        } else {
            ret.transitions[i] = {
                {i + 1, 1}
            };
        }
    }

    // We start at the beginning of V
    ret.transitions[state_pos::magic] = {
        {v_states.first, 1}
    };

    return ret;
}
