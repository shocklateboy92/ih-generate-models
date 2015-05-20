#ifndef MUTATIONRATIOS_H
#define MUTATIONRATIOS_H

#include "common.h"

mutation_probs_t readMutationProbsFile(const char * fileName);

double fetch_mutation_ratio(RunConfig config, std::string sequence,
                            std::size_t index, char to_nt);

#endif // MUTATIONRATIOS_H

