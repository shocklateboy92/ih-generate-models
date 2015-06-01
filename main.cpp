#include "common.h"
#include "a-score.h"
#include "mutation-ratios.h"
#include "init.h"
#include "viterbi.h"

#include <vector>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Incorrect Args!" << std::endl;
        return -2;
    }

    const RunConfig &globalConfig = prepareConfig(argc, argv);
    assert(globalConfig.j_repo.size() == 15);

    pugi::xml_document doc;
    auto xml_parse_result = doc.load_file(argv[1]);
    if (!xml_parse_result) {
        std::cerr << "Failed to load '" << argv[1] << "': " <<
                xml_parse_result.description() << std::endl;
        return -1;
    }

    auto nodes = doc.select_nodes("//BlastOutput/BlastOutput_iterations/Iteration");

//    using boost::adaptors::transform;

    #define BIND_CONFIG(x) std::bind(x, globalConfig, std::placeholders::_1)

    auto i1 = readRepertoire(argv[2]);
    assert (nodes.size() == i1.size());

//    std::vector<std::pair<const pugi::xpath_node&, FastaSequence>> pairs;
//    for (std::size_t i = 0; i < nodes.size(); i++) {
//        pairs.push_back({nodes[i], i1[i]});
//    }

//    // Parse the Iteration nodes
//    auto r1 = transform(pairs, [](auto b) {return std::make_pair(parseBlastOutput(b.first), b.second);});
//    auto r2 = transform(r1, [](auto b) -> SequenceInfo {
//        return {
//                b.second.name(),
//                b.second.c_str(),
//                b.first,
//                // Calculate A Score from Blast Results
//                calculateAScore(b.first)
//        };
//    });
//    auto r3 = transform(r2, BIND_CONFIG(buildModel));
//    auto r4 = transform(r3, do_viterbi);

//    for (auto v : r4) {
//        for (auto s : v)
//            std::cout << s.name << std::endl;
//    }
    tbb::concurrent_vector<std::vector<StateInfo>> vec(nodes.size());

    tbb::parallel_for(0ul, nodes.size(), [&](std::size_t i) {
        auto br = parseBlastOutput(nodes[i]);

        auto path = do_viterbi(
                    buildModel(
                        globalConfig, {
                            i1[i].name(),
                            i1[i].c_str(),
                            br,
                            calculateAScore(br),
                        })
                    );

        vec[i] = path;
    });

    for (auto path : vec) {
        for (auto v : path) {
            std::cout << v.name << " ";
        }
        std::endl(std::cout);
    }

    return 0;
}


