#include "mutability.h"
#include <regex>

//new coverage values are based on the germline frequency of the hotspots
//for 59982 4mers in all germline IGHV, 2284 RGYW, 2116 WRCY, 4349 WAN and 51233 Non-HS
double RGYW_Nucleotide_coverage = 0.0381l;
double WRCY_Nucleotide_coverage = 0.0353l;
double WAN_Nucleotide_coverage = 0.0725l;
double NO_HOTSPOT_Nucleotide_coverage = 0.854l;

//mutation probabilities calculated from Barington 2007 JI paper
//4592 muts total, 598 RGYW, 633 WAN, 558 WRCY, 2803 in non-hotspots
double RGYW_mutation_prob = 0.130l;
double WRCY_mutation_prob = 0.122l;
double WAN_mutation_prob = 0.138l;
double NO_HOTSPOT_mutation_prob = 0.610l;

double NO_HOTSPOT_MUTABILITY_SCORE = (NO_HOTSPOT_mutation_prob/NO_HOTSPOT_Nucleotide_coverage);
double WAN_MUTABILITY_SCORE = (WAN_mutation_prob/WAN_Nucleotide_coverage);
double NAN_MUTABILITY_SCORE = ((WAN_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGYW_MUTABILITY_SCORE = (RGYW_mutation_prob / RGYW_Nucleotide_coverage);
double NGYW_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGYN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double RGNN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE*0.25) + (NO_HOTSPOT_MUTABILITY_SCORE*0.75));
double WRCY_MUTABILITY_SCORE = (WRCY_mutation_prob/WRCY_Nucleotide_coverage);
double NRCY_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
double NNCY_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE * 0.25)+(NO_HOTSPOT_MUTABILITY_SCORE*0.75);
double WRCN_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);

std::vector<std::pair<std::regex, double>> hotspots = {
    {std::regex(".[AT]A.."), WAN_MUTABILITY_SCORE},
    {std::regex(".NA.."), NAN_MUTABILITY_SCORE},

    {std::regex("..A.."), NO_HOTSPOT_MUTABILITY_SCORE},

    {std::regex(".[AG].[CT][AT]"), RGYW_MUTABILITY_SCORE},
    {std::regex(".N.YW"), NGYW_MUTABILITY_SCORE},
    {std::regex(".R.N."), RGNN_MUTABILITY_SCORE},
    {std::regex(".R.YN"), RGYN_MUTABILITY_SCORE},
    {std::regex("[AT][AG].[CT]."), WRCY_MUTABILITY_SCORE},
    {std::regex("[AT][AG].N."), WRCN_MUTABILITY_SCORE},
    {std::regex(".N.[CT]."), NNCY_MUTABILITY_SCORE},
    {std::regex("N[AG].[CT]."), NRCY_MUTABILITY_SCORE}
};
// TODO: Finish regexs for hotspots

double fetch_mutability_score(const seq_t &fstr, std::size_t index) {

    // fetch the current penta-nucleotide
    std::string penta_nucleotide = get_penta_nucleotide(fstr, index);

    // check for matching hotspots
    for (auto p : hotspots) {
        if (std::regex_match(penta_nucleotide, p.first)) {
            return p.second;
        }
    }

    return NO_HOTSPOT_MUTABILITY_SCORE;
}
