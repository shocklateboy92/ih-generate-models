#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "lib/pugixml/pugixml.hpp"
#include "lib/fastareader/fastareader.h"
#include <assert.h>

using namespace std;

struct BlastResult {
    std::string v_name;

    double a_score;
};

struct RunConfig {
    std::vector<FastaSequence> j_repo;
    std::vector<FastaSequence> d_repo;
};

BlastResult parseBlastOutput(const pugi::xpath_node &node);

double calculateAScore(const BlastResult &br);

RunConfig prepareConfig(int argc, char *pString[]);

std::vector<FastaSequence> readRepertoire(const char *);

#include "common.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Incorrect Args!" << std::endl;
        return -2;
    }

    pugi::xml_document doc;
    std::cerr << "loading " << argv[1] << std::endl;
    auto result = doc.load_file(argv[1]);
    if (!result) {
        std::cerr << result.description() << endl;
        return -1;
    }

    auto nodes = doc.select_nodes("//BlastOutput/BlastOutput_iterations/Iteration");

    // Parse the Iteration nodes
    std::vector<BlastResult> results;
    results.reserve(nodes.size());
    std::transform(nodes.begin(), nodes.end(),
            results.begin(), parseBlastOutput);

    // Calculate A Score from Blast Results
    std::for_each(results.begin(), results.end(),
            [](BlastResult &br) {
                br.a_score = calculateAScore(br);
            });

    RunConfig globalConfig = prepareConfig(argc, argv);
    assert(globalConfig.j_repo.size() == 15);


    return 0;
}

RunConfig prepareConfig(int argc, char *argv[]) {
    return {
            readRepertoire("IGHJRepertoire.fasta"),
            readRepertoire("IGHDRepertoire.fasta")
    };
}

std::vector<FastaSequence> readRepertoire(const char *repo_path) {
    std::vector<FastaSequence> ret;
    FastaReader fr;
    std::ifstream is(repo_path);
    assert(is.is_open());

    while (!is.eof()) {
        auto seq = fr.next(is);
        ret.push_back(*seq.get());
    }

    return ret;
}

BlastResult parseBlastOutput(const pugi::xpath_node &node) {
    return {
            // Only considering the most likely V-gene for now
            node.node()
                    .child("Iteration_hits")
                    .child("Hit")
                    .child("Hit_id")
                    .value()
    };
}

double calculateAScore(const BlastResult &br) {
    return 7;
}
