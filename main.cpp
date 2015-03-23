#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib/pugixml/pugixml.hpp"
#include <assert.h>

using namespace std;

struct BlastResult {
    double a_score;
};

BlastResult parseBlastOutput(const pugi::xpath_node &node);
double calculateAScore(const BlastResult &br);

int main(int argc, char *argv[])
{
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
                  [](BlastResult &br){
        br.a_score = calculateAScore(br);
    });


    return 0;
}

BlastResult parseBlastOutput(const pugi::xpath_node &node) {
    return {};
}

double calculateAScore(const BlastResult &br) {
    return 7;
}
