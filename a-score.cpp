#include "a-score.h"
#include "common.h"

constexpr int COMMON_AREA_START_POS = 129;
constexpr int COMMON_AREA_END_POS = 243;
constexpr int COMMON_AREA_NUCL_LENGTH = 115;

double calculateAScore(const BlastResult &br) {
    // currently only implementing the 'long' Vgene prob, assuming
    // that the commone region is included in the match
    assert(br.v_match_start < 128);

    // make sure v gene actually contains the common region
    assert(br.v_match_start + br.v_aligned_string.length() < COMMON_AREA_END_POS);
    assert(br.v_match_start <= COMMON_AREA_START_POS);

    assert(br.v_aligned_string.length() <= br.v_match_string.length());

    return 7;
}
