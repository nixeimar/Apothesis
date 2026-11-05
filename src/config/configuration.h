#ifndef APOTHESIS_CONFIGURATION_H
#define APOTHESIS_CONFIGURATION_H

#include <cstdint>
#include <map>
#include <optional>
#include "../lattice/atomistic_surface.h"
#include "../processes/atomistic_kmc.h"
#include "../processes/zgb_growth.h"
#include "../processes/zgb_3d.h"
#include <string>
#include <variant>
#include <utility>
#include <vector>

namespace Config {

// Independent of the file format: models can read custom, typed parameters here.
class Value {
public:
    using Array = std::vector<Value>;
    using Object = std::map<std::string, Value>;
    using Storage = std::variant<std::nullptr_t, bool, std::int64_t, double,
                                 std::string, Array, Object>;
    Value() : data(nullptr) {}
    template<class T> explicit Value(T value) : data(std::move(value)) {}
    template<class T> const T& as() const { return std::get<T>(data); }
    double number() const;
    const Value& at(const std::string& key) const { return as<Object>().at(key); }
    bool contains(const std::string& key) const;
private:
    Storage data;
};

struct Process {
    std::string id, type, equation;
    // Adapter for existing kinetic models; named parameters remain in document.
    std::vector<std::string> modelArguments;
};

struct Configuration {
    Value document;
    std::optional<ZGB::Options> zgb;
    std::optional<ZGB3D::Options> zgb3d;
    std::optional<Atomistic::Options> atomisticSurface;
    std::optional<Atomistic::RunOptions> atomisticSimulation;
    double temperature = 0, pressure = 0, startTime = 0, endTime = 0;
    double logInterval = 0, latticeInterval = 0;
    int seed = 0, x = 0, y = 0, height = 0;
    std::string latticeType, initialSpecies, heightsFile, speciesFile;
    std::vector<std::string> growth, coverage;
    std::vector<Process> processes;
};

Configuration readJson(const std::string& filename);
} // namespace Config
#endif
