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
        std::stoul(hsps.child_value("Hsp_hit-from")) -1,
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

struct state_pos {
    enum {
        magic = 0,
        v_end_no_exo,
        v_end_exo,
        n1_start,
        n1_end,
        n2_start,
        n2_end,
        end_state
    };
};

// Using hardcoded values for now, will read from file later
extern std::vector<std::vector<double>> d_start_exo_probs;
extern std::vector<std::vector<double>> d_end_exo_probs;
probs_list_t v_end_exo_probs = { 1.7849193927612356E-6, 4.337659583280163E-5, 6.481808081758074E-4, 0.005955808351520267, 0.03365031276420145, 0.11690722185950304, 0.24974535503543732, 0.32806306897882054, 0.2649848906871161 };

// Fill figure out how to do this later
static const double early_exit_prob = 0;

void createStates(const RunConfig &config, const SequenceInfo &input,
                  std::function<double(double)> exp_decay_fn,
                  probs_list_t start_exo_probs, probs_list_t end_exo_probs,
                  const seq_t &fstr, HiddenMarkovModel &ret)
{
    // TODO: abstract this into class
    std::size_t pos_start_no_exo = ret.states.size();
    std::size_t pos_start_exo = pos_start_no_exo + 1;
    std::size_t gene_start = pos_start_exo + 1;

    std::size_t pos_end_no_exo = gene_start + fstr.size();
    std::size_t pos_end_exo = pos_end_no_exo + 1;
    std::size_t gene_end = pos_end_no_exo - 1;

    ret.states.push_back({"start_no_exo"});
    ret.transitions.push_back({
                                  {state_pos::magic, early_exit_prob},
                                  {gene_start, 1 - early_exit_prob}
                              });
    ret.states.push_back({"start_exo"});
    ret.transitions.push_back({
                                  {state_pos::magic, early_exit_prob},
                                  {
                                      pos_start_no_exo,
                                      start_exo_probs.empty()
                                      ? 0
                                      : start_exo_probs.front()
                                  },
                              });
    for (std::size_t i = 0; i < start_exo_probs.size(); i++) {
        ret.transitions[pos_start_exo].insert({
                                                  gene_start + i,
                                                  start_exo_probs[i]
                                              });
    }
    // TODO: look at why the last has a 1 prob of going in
    // at line 1457 of VnDnJCnoC.java

    assert (ret.states.size() == gene_start);
    for (std::size_t i = 0; i < fstr.size(); i++) {
        ret.states.push_back(
                    createState(
                        config,
                        input,
                        exp_decay_fn(i),
                        fstr,
                        i)
                    );

        auto end_probs_start = fstr.size() - end_exo_probs.size();
        double end_exo_prob = i > end_probs_start
                ? end_exo_prob
                : 0;

        double forward_prob = 1 - early_exit_prob - end_exo_prob;

        ret.transitions.push_back({
                                      {gene_start + i + 1, forward_prob},
                                      {state_pos::magic, early_exit_prob},
                                      {pos_end_exo, end_exo_prob}
                                  });

        assert(ret.states.size() == ret.transitions.size());
    }

    ret.states.push_back({"end_no_exo"});
    ret.states.push_back({"end_exo"});
    assert (ret.states.size() == pos_end_exo + 1);

    ret.transitions.push_back({});
    ret.transitions.push_back({});
    assert(ret.states.size() == ret.transitions.size());
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
                input.blast_result.v_match_start,
                full_v_seq.length());

    auto exp_decay_fn = [](double i) {return std::exp(EXP_DECAY_INDEX_CONST * i);};

    // Initialize the transitions with the dot states built in
    ret.states = std::vector<StateInfo>(state_pos::end_state);
    ret.transitions = transitions_t(state_pos::end_state);

    std::size_t
            v_start_no_exo = state_pos::end_state,
            v_start_exo  = v_start_no_exo + 1,
            v_end_no_exo = v_start_exo + fstr.length() + 1,
            v_end_exo = v_end_no_exo + 1;

    // Now make a state for each NT in the V-Gene
    createStates(config, input, exp_decay_fn, {}, v_end_exo_probs, fstr, ret);

    // We start at the beginning of V
    ret.transitions[state_pos::magic] = {
        {state_pos::end_state, 1}
    };
    ret.transitions[v_end_no_exo] = {
        {state_pos::magic, early_exit_prob},
        {v_end_exo, 1 - early_exit_prob},
    };
    ret.transitions[v_end_exo] = {
        {state_pos::magic, early_exit_prob},
        {state_pos::n1_start, 1 - early_exit_prob}
    };
    ret.transitions[state_pos::n1_start] = {
        {state_pos::n1_end, 1}
    };

    assert (ret.states[v_start_no_exo].name == "start_no_exo");
    assert (ret.states[v_start_exo].name == "start_exo");
    assert (ret.states[v_end_no_exo].name == "end_no_exo");

    // Now, do the same for all possible D-genes
    for (std::size_t i = 0; i < config.d_repo.size(); i++) {
        std::size_t pos_start_no_exo = ret.states.size();
        std::size_t pos_start_exo = pos_start_no_exo + 1;
        std::size_t gene_start = pos_start_exo + 1;

        std::size_t pos_end_no_exo = gene_start + config.d_repo.at(i).size();
        std::size_t pos_end_exo = pos_end_no_exo + 1;
        std::size_t gene_end = pos_end_no_exo - 1;

        createStates(
            config,
            input,
            [](double i) {return std::exp(EXP_DECAY_INDEX_CONST * i);},
            d_start_exo_probs.at(i),
            d_end_exo_probs.at(i),
            config.d_repo.at(i).c_str(),
            ret
        );

        assert (ret.states[pos_start_exo].name == "start_exo");
        assert (ret.states[pos_end_no_exo].name == "end_no_exo");

        double num_genes = config.d_repo.at(i).size();
        ret.transitions[state_pos::n1_end].insert({
            {pos_start_exo,  (1 - early_exit_prob) / num_genes}
        });

        ret.transitions[pos_end_no_exo] = {
            {state_pos::magic, early_exit_prob},
            {state_pos::n2_start, 1 - early_exit_prob}
        };
        ret.transitions[pos_end_exo] = ret.transitions[pos_end_no_exo];
    }

    ret.transitions[state_pos::n2_start] = {
        {state_pos::n2_end, 1}
    };

    // FIXME: extract function
    // Now, repeat for all the J genes as well
    for (std::size_t i = 0; i < config.j_repo.size(); i++) {
        std::size_t pos_start_no_exo = ret.states.size();
        std::size_t pos_start_exo = pos_start_no_exo + 1;
        std::size_t gene_start = pos_start_exo + 1;

        std::size_t pos_end_no_exo = gene_start + config.j_repo.at(i).size();
        std::size_t pos_end_exo = pos_end_no_exo + 1;
        std::size_t gene_end = pos_end_no_exo - 1;

        createStates(
            config,
            input,
            [](double i) {return std::exp(EXP_DECAY_INDEX_CONST * i);},
            d_start_exo_probs.at(i),
            d_end_exo_probs.at(i),
            config.j_repo.at(i).c_str(),
            ret
        );

        assert (ret.states[pos_start_exo].name == "start_exo");
        assert (ret.states[pos_end_no_exo].name == "end_no_exo");

        double num_genes = config.d_repo.at(i).size();
        ret.transitions[state_pos::n2_end].insert({
            {pos_start_exo,  (1 - early_exit_prob) / num_genes}
        });

        ret.transitions[pos_end_no_exo] = {
            {state_pos::magic, 1},
        };
        ret.transitions[pos_end_exo] = ret.transitions[pos_end_no_exo];
    }


    return ret;
}

#ifndef SKIP_TESTS

bool check_choke(const HiddenMarkovModel &model, std::size_t cur_state,
                 std::size_t choke_state) {
    if (cur_state == choke_state || cur_state == state_pos::magic) {
        return true;
    }

    bool has_path_out = false;
    for (auto t : model.transitions[cur_state]) {
        if (t.second > 0) {
            has_path_out |= check_choke(model, t.first, choke_state);
        }
    }

    assert(has_path_out);
    REQUIRE(has_path_out);

    return has_path_out;
}

TEST_CASE("testing continuity of model", "[transitions]") {
    RunConfig config = prepareConfig(0, nullptr);
    SequenceInfo input = {
        "test",
        "<future>",
        BlastResult {
            "IGHV3-9*01(L1)",
            "GAGTCTGGGGGAGGCTTGGTACAGCCTGGCAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGGC",
            15,
            "GAGTCGGGGGGAGGCTGGGTACAGCCTGGCAGGTCCCTGAGACTCTCCTGTTCAGCCTCTGGACTCACCTTTGATGATTATGCCATGCACTGGGTCCGGC",
        },
        0.06
    };

    auto model = buildModel(config, input);

    // everyone gets to end of V
    check_choke(model, state_pos::end_state, state_pos::n1_start);

    // everyone gets to end of D
    check_choke(model, state_pos::end_state, state_pos::n2_start);

    check_choke(model, state_pos::end_state, state_pos::magic);
}

#endif
