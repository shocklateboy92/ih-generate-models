#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "lib/pugixml/pugixml.hpp"
#include <assert.h>

using namespace std;

struct BlastResult {

};

BlastResult parseBlastOutput(pugi::xml_node node);

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
    std::vector<BlastResult> results;
    results.reserve(nodes.size());

    for (auto n : nodes) {
        std::cout << n.node().name() << std::endl;
    }

    return 0;
}

BlastResult parseBlastOutput(pugi::xml_node node) {
    return {};
}
