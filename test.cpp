#define CATCH_CONFIG_MAIN

#include "common.h"
#include "mutation-ratios.h"

TEST_CASE("mutation ratio", "[mr]") {
    RunConfig config {
        .mutation_probs = readMutationProbsFile("Mutation_spectrum.txt")
    };

    REQUIRE(fetch_mutation_ratio(config, "AAT", 1, 'G') == 0.533465028);
}
