#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include <assert.h>

#include "common.h"

BlastResult parseBlastOutput(const pugi::xpath_node &node);

double calculateAScore(const BlastResult &br);

RunConfig prepareConfig(int argc, char *pString[]);

std::vector<FastaSequence> readRepertoire(const char *);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Incorrect Args!" << std::endl;
        return -2;
    }

    pugi::xml_document doc;
    std::cerr << "loading " << argv[1] << std::endl;
    auto result = doc.load_file(argv[1]);
    if (!result) {
        std::cerr << result.description() << std::endl;
        return -1;
    }

    auto nodes = doc.select_nodes("//BlastOutput/BlastOutput_iterations/Iteration");

    //TODO: Replace this with Boost Range Adapters or pimped out std::bind
    // See: http://stackoverflow.com/questions/4302006/how-to-bind-a-constructor-in-c

    // Parse the Iteration nodes
    std::vector<BlastResult> results;
    results.reserve(nodes.size());
    std::transform(nodes.begin(), nodes.end(),
            results.begin(), parseBlastOutput);

    std::vector<SequenceInfo> results2;
    std::transform(results.begin(), results.end(),
    results2.begin(), [](BlastResult result1){
                return {
                        "name",
                        "data",
                        result1,
                        // Calculate A Score from Blast Results
                        calculateAScore(result1)
                };
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
    // Only considering the most likely V-gene for now
    auto hit = node.node().child("Iteration_hits").child("Hit");
    auto hsps = hit.child("Hit_hsps").child("Hsp");
    return {
            hit.child("Hit_id").value(),
            hsps.child_value("Hsp_qseq"),
            std::atol(hsps.child_value("Hsp_query-from"))
    };
}

double calculateAScore(const BlastResult &br) {
    return 7;
}
