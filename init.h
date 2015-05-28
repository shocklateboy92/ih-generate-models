#ifndef INIT_H
#define INIT_H

#include "common.h"

RunConfig prepareConfig(int argc, char *argv[]);

std::vector<FastaSequence> readRepertoire(const char *repo_path);

#endif // INIT_H
