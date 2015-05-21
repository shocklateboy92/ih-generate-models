#define CATCH_CONFIG_MAIN

#include "common.h"
#include "mutation-ratios.h"

TEST_CASE("mutation ratio", "[mr]") {
    RunConfig config {
        .mutation_probs = readMutationProbsFile("Mutation_spectrum.txt")
    };

    REQUIRE(fetch_mutation_ratio(config, "AAT", 1, 'G') == 0.533465028);
}

TEST_CASE("testing get_penta_nucleotide", "[get_nt_5]") {
    REQUIRE(get_penta_nucleotide("AGTCGAC", 2) == "AGTCG");
    REQUIRE(get_penta_nucleotide("AGTCGAC", 0) == "**AGT");
    REQUIRE(get_penta_nucleotide("AGTCGAC", 5) == "CGAC*");
}

TEST_CASE("testing get_tri_nucleotide", "[get_nt_3]") {
    REQUIRE(get_tri_nucleotide("AGTGAC", 1) == "AGT");
    REQUIRE(get_tri_nucleotide("AGTGAC", 0) == "*AG");
    REQUIRE(get_tri_nucleotide("AGTGAC", 5) == "AC*");
}

TEST_CASE("testing mutation_spectrum.txt parsing") {
    mutation_probs_t probs = readMutationProbsFile("Mutation_spectrum.txt");

    REQUIRE(probs.at("AAGT") == 0.533465028);
    REQUIRE(probs.at("ACTT") == 0.41283902);
    REQUIRE(probs.at("AGCT") == 0.404121647);
}
