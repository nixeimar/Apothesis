#include "../src/config/configuration.h"
#include "../src/extLibs/nlohmann/json.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>

using Json = nlohmann::json;
namespace fs = std::filesystem;
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
int main(int argc, char** argv) {
    fs::path temp = fs::temp_directory_path() / ("apothesis-input-test-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        if (argc != 2) throw std::runtime_error("expected example JSON path");
        fs::create_directory(temp);
        Json original;
        std::ifstream(argv[1]) >> original;
        const auto file = temp / "input.json";
        auto load = [&](const Json& input) {
            std::ofstream(file) << input.dump(2);
            return Config::readJson(file.string());
        };
        auto invalid = [&](Json input, const std::string& expected) {
            try { load(input); }
            catch (const std::runtime_error& e) {
                check(std::string(e.what()).find(expected) != std::string::npos,
                      "wrong diagnostic: " + std::string(e.what())); return;
            }
            throw std::runtime_error("accepted invalid input: " + expected);
        };
        auto c = load(original);
        check(c.startTime == 5 && c.endTime == 505, "duration must be relative to start time");
        check(c.processes[1].equation == "O2 + 2* -> 2O*", "stoichiometry was not preserved");
        check(c.processes[1].modelArguments.back() == "all", "neighbor option lost");
        auto custom = original;
        custom["parameters"] = {{"nullable", nullptr}, {"integer", 9007199254740993LL}, {"label", "custom"}};
        custom["species"]["CO"]["parameters"]["binding"] = {{"enabled", true}, {"values", {1, 2.5, "x"}}};
        custom["lattice"]["parameters"]["custom"] = Json::array({false, nullptr});
        c = load(custom);
        check(c.document.at("parameters").at("integer").as<std::int64_t>() == 9007199254740993LL, "integer precision lost");
        check(c.document.at("species").at("CO").at("parameters").at("binding").at("enabled").as<bool>(), "nested boolean lost");
        check(c.document.at("processes").as<Config::Value::Array>()[1].at("parameters").as<Config::Value::Object>().empty(), "missing default parameter bag");
        auto duplicateEquation = original;
        duplicateEquation["processes"].push_back(duplicateEquation["processes"][0]);
        duplicateEquation["processes"][3]["id"] = "second_co_adsorption";
        check(load(duplicateEquation).processes.size() == 4, "distinct process ids must preserve duplicate equations");
        auto simple = original;
        simple["processes"][0]["rate"] = {{"model", "simple"}, {"parameters", {
            {"sticking_coefficient", 0.1}, {"molar_fraction", 0.2}, {"site_density", 1e19}}}};
        check(std::stod(load(simple).processes[0].modelArguments[4]) == 0.028, "species molecular weight not resolved");
        auto arrhenius = original;
        arrhenius["processes"][2]["rate"] = {{"model", "arrhenius"}, {"parameters", {{"prefactor", 1e13}, {"activation_energy", 50000}}}};
        check(load(arrhenius).processes[2].modelArguments.size() == 3, "Arrhenius mapping failed");
        auto bad = original; bad["schema_version"] = 2; invalid(bad, "schema_version");
        bad = original; bad["system"]["parameters"]["temperature"] = "hot"; invalid(bad, "temperature");
        bad = original; bad["simulation"]["duration"] = -1; invalid(bad, "duration");
        bad = original; bad["simulation"]["random_seed"] = 1.5; invalid(bad, "random_seed");
        bad = original; bad["lattice"]["parameters"]["dimensions"] = {10}; invalid(bad, "dimensions");
        bad = original; bad["lattice"]["type"] = "Unknown"; invalid(bad, "lattice.type");
        bad = original; bad["processes"][0]["reactants"][0]["species"] = "missing"; invalid(bad, "undefined species");
        bad = original; bad["processes"][1]["id"] = "co_adsorption"; invalid(bad, "duplicate process id");
        bad = original; bad["processes"][0]["rate"]["parameters"].erase("value"); invalid(bad, "required field");
        bad = original; bad["processes"][0]["rate"]["model"] = "unknown"; invalid(bad, "unsupported rate model");
        bad = original; bad["processes"][0]["parameters"] = 5; invalid(bad, "expected an object");
        bad = original; bad["processes"][0]["type"] = "reaction"; invalid(bad, "reaction reactants");
        bad = original; bad["system"]["temperature"] = 300; invalid(bad, "unknown field");
        bad = original; bad["species"]["CO*"] = Json::object(); invalid(bad, "name must not");
        bad = simple; bad["species"]["CO"]["parameters"].erase("molecular_weight"); invalid(bad, "molecular_weight");
        // Rectangular matrices and paths relative to the configuration file.
        auto matrix = original;
        auto& lp = matrix["lattice"]["parameters"];
        lp["dimensions"] = {4, 3}; lp.erase("initial_height");
        lp["initial_height_file"] = "heights.dat";
        lp["initial_species_file"] = "species.dat";
        std::ofstream(temp / "heights.dat") << "1 2 3 4\n5 6 7 8\n9 10 11 12\n";
        std::ofstream(temp / "species.dat") << "X X X CO*\nX X X X\nX X X X\n";
        check(load(matrix).heightsFile == (temp / "heights.dat").string(), "relative matrix path failed");
        std::ofstream(temp / "heights.dat") << "1 2\n";
        invalid(matrix, "each row");
        fs::remove(temp / "heights.dat"); invalid(matrix, "cannot open lattice matrix");
        std::ofstream(file) << "{\"schema_version\":1,\"schema_version\":1}";
        try { Config::readJson(file.string()); throw std::logic_error("duplicate keys accepted"); }
        catch (const std::runtime_error& e) { check(std::string(e.what()).find("duplicate key") != std::string::npos, "duplicate key diagnostic"); }
        std::ofstream(file) << "{ invalid json";
        try { Config::readJson(file.string()); throw std::logic_error("malformed JSON accepted"); }
        catch (const std::runtime_error& e) { check(std::string(e.what()).find("parse_error") != std::string::npos, "syntax diagnostic"); }
        fs::remove_all(temp);
        std::cout << "All JSON input tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        fs::remove_all(temp);
        std::cerr << e.what() << '\n'; return 1;
    }
}
