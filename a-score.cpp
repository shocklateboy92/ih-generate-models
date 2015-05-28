#include "a-score.h"
#include "common.h"

constexpr int COMMON_AREA_START_POS = 129;
constexpr int COMMON_AREA_END_POS = 243;
constexpr int COMMON_AREA_NUCL_LENGTH = 115;

constexpr double A_SLOPE = 0.0064999999999999997l;
constexpr double B_ADDITION = 0.0015l;

std::size_t difference(std::string a, std::string b) {
    std::size_t diff = 0;
    for (auto ai = a.begin(), bi = b.begin();
         ai != a.end() && bi != b.end();
         ai++, bi++) {
        *ai == *bi && diff++;
    }
    return diff;
}

double calculateAScore(const BlastResult &br) {
    // currently only implementing the 'long' Vgene prob, assuming
    // that the commone region is included in the match
    // TODO: Implement 'short' Vgene prob
    assert(br.v_match_start < 128);

    // make sure v gene actually contains the common region
    assert(br.v_match_start + br.v_input_string.length() > COMMON_AREA_END_POS);
    assert(br.v_match_start <= COMMON_AREA_START_POS);

//    assert(br.v_input_string.length() == br.v_match_string.length());

    // TODO: deal with these properly - i.e. ignore them in the
    // 		 mutation count.
    assert(br.v_input_string.find('n') == std::string::npos);
    assert(br.v_input_string.find('x') == std::string::npos);

    double num_mutations = difference(br.v_input_string,
                                      br.v_match_string);
    double num_n_x_in_cr = 0; //TODO: Fix this

    double num_mutations_in_cr = difference(
                br.v_input_string.substr(COMMON_AREA_START_POS,
                                         COMMON_AREA_NUCL_LENGTH),
                br.v_match_string.substr(COMMON_AREA_START_POS,
                                         COMMON_AREA_NUCL_LENGTH)
                );

    double cr_to_v_gene_ratio = COMMON_AREA_NUCL_LENGTH / num_mutations;
    double avg_mutations_in_cr = num_mutations * cr_to_v_gene_ratio;

    double n_x_addition =
            (avg_mutations_in_cr / COMMON_AREA_NUCL_LENGTH)
            * A_SLOPE * num_n_x_in_cr;

    double prob_mutation = num_mutations_in_cr * A_SLOPE * B_ADDITION;

    return prob_mutation + n_x_addition;
}
