#include "common.h"
#include "a-score.h"
#include "mutation-ratios.h"
#include "init.h"
#include "viterbi.h"

#include <vector>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_scheduler_init.h>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Incorrect Args!" << std::endl;
        return -2;
    }

    auto num_threads = tbb::task_scheduler_init::automatic;
    if (argc == 4) {
        num_threads = atoi(argv[3]);
    }
    tbb::task_scheduler_init init(num_threads);

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

    auto i1 = readRepertoire(argv[2]);
    assert (nodes.size() == i1.size());

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


