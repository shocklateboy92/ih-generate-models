#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include <assert.h>
#include "common.h"

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/numeric.hpp>
//#include <boost/phoenix/phoenix.hpp>

BlastResult parseBlastOutput(const pugi::xpath_node &node);

double calculateAScore(const BlastResult &br);

RunConfig prepareConfig(int argc, char *pString[]);

std::vector<FastaSequence> readRepertoire(const char *);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Incorrect Args!" << std::endl;
        return -2;
    }

    RunConfig globalConfig = prepareConfig(argc, argv);
    assert(globalConfig.j_repo.size() == 15);

    pugi::xml_document doc;
    auto xml_parse_result = doc.load_file(argv[1]);
    if (!xml_parse_result) {
        std::cerr << "Failed to load '" << argv[1] << "': " <<
                xml_parse_result.description() << std::endl;
        return -1;
    }

    auto nodes = doc.select_nodes("//BlastOutput/BlastOutput_iterations/Iteration");

    using boost::adaptors::transform;

    #define BIND_CONFIG(x) std::bind(x, globalConfig, std::placeholders::_1)

    // Parse the Iteration nodes
    auto r1 = transform(nodes, parseBlastOutput);
    auto r2 = transform(r1, [](BlastResult result1) -> SequenceInfo {
        return {
                "name",
                "data",
                result1,
                // Calculate A Score from Blast Results
                calculateAScore(result1)
        };
    });
    auto r3 = transform(r2, BIND_CONFIG(buildModel));

    std::ifstream v_probs("v_probs.txt"), v_names("v_names.txt");
    assert(v_probs.is_open() && v_names.is_open());

    for (auto s : r2) {
        std::string ev_name;
        v_names >> ev_name;
        assert(s.blast_result.v_name == ev_name);
    }

    std::size_t matches = 0;
    for (const HiddenMarkovModel &model : r3) {
        for (auto s : model.states) {
            for (auto p : s.emission_probs) {
                double ep;
                v_probs >> ep;
                assert(p == std::round(ep));
                matches++;
            }
        }
    }

    return 0;
}

RunConfig prepareConfig(int argc, char *argv[]) {
    return {
            readRepertoire("IGHJRepertoire.fasta"),
            readRepertoire("IGHDRepertoire.fasta"),
            boost::accumulate(
                    readRepertoire("IGHVRepertoire.fasta"),
                    seq_map_t(),
                    [](seq_map_t &map, FastaSequence f) -> seq_map_t & {
                        map.insert({f.name(), f});
                        return map;
                    }
            )
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
