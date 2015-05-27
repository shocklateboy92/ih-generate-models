//
// Created by Lasath on 3/24/2015.
//

#include "common.h"
#include "mutation-ratios.h"
#include "mutability.h"

template <class C, class F>
auto my_transform(const C &c, F f) {
    std::vector<decltype(f(c[0]))> v(c.size());
    std::transform(c.begin(), c.end(), v.begin(), f);
    return v;
}

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
             my_transform(TRACK, [&](char c) -> double {

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

    return ret;
}
