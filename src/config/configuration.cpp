#include "configuration.h"
#include "../extLibs/nlohmann/json.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace Config {
namespace {
using Json = nlohmann::json;
[[noreturn]] void fail(const std::string& path, const std::string& message) {
    throw std::runtime_error(path + ": " + message);
}
void object(const Json& j, const std::string& path) {
    if (!j.is_object()) fail(path, "expected an object");
}
void keys(const Json& j, const std::string& path, std::initializer_list<std::string> allowed) {
    object(j, path);
    for (auto it = j.begin(); it != j.end(); ++it)
        if (std::find(allowed.begin(), allowed.end(), it.key()) == allowed.end())
            fail(path + "." + it.key(), "unknown field; custom properties belong in parameters");
    if (j.contains("description") && !j.at("description").is_string())
        fail(path + ".description", "expected a string");
    if (j.contains("parameters")) object(j.at("parameters"), path + ".parameters");
}
const Json& required(const Json& j, const std::string& key, const std::string& path) {
    object(j, path);
    if (!j.contains(key)) fail(path + "." + key, "required field is missing");
    return j.at(key);
}
std::string text(const Json& j, const std::string& path) {
    if (!j.is_string() || j.get<std::string>().empty()) fail(path, "expected a nonempty string");
    return j.get<std::string>();
}
double number(const Json& j, const std::string& path, bool positive = false) {
    if (!j.is_number()) fail(path, "expected a number");
    double n = j.get<double>();
    if (!std::isfinite(n) || (positive ? n <= 0 : n < 0))
        fail(path, positive ? "must be finite and greater than zero" : "must be finite and nonnegative");
    return n;
}
int integer(const Json& j, const std::string& path, bool positive = false) {
    double n = number(j, path, positive);
    if (!j.is_number_integer() || n > std::numeric_limits<int>::max()) fail(path, "expected an integer in range");
    return static_cast<int>(n);
}
Value convert(const Json& j) {
    if (j.is_object()) {
        Value::Object result;
        for (auto it = j.begin(); it != j.end(); ++it) result.emplace(it.key(), convert(it.value()));
        return Value(std::move(result));
    }
    if (j.is_array()) {
        Value::Array result;
        for (const auto& v : j) result.push_back(convert(v));
        return Value(std::move(result));
    }
    if (j.is_string()) return Value(j.get<std::string>());
    if (j.is_boolean()) return Value(j.get<bool>());
    if (j.is_number_unsigned() && j.get<std::uint64_t>() > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        fail("parameters", "integer exceeds signed 64-bit range");
    if (j.is_number_integer()) return Value(j.get<std::int64_t>());
    if (j.is_number()) return Value(j.get<double>());
    return Value();
}
std::string numericArgument(double n) {
    std::ostringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << n;
    return out.str();
}
struct Term { std::string label, phase; int coefficient; bool vacancy; };
std::vector<Term> terms(const Json& array, const Json& species, const std::string& path) {
    if (!array.is_array() || array.empty()) fail(path, "expected a nonempty array");
    std::vector<Term> result;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < array.size(); ++i) {
        const auto& j = array[i];
        auto p = path + "[" + std::to_string(i) + "]";
        keys(j, p, {"species", "phase", "coefficient", "vacancy"});
        Term term{"", "surface", 1, false};
        if (j.contains("coefficient")) term.coefficient = integer(j.at("coefficient"), p + ".coefficient", true);
        if (j.contains("vacancy")) {
            if (!j.at("vacancy").is_boolean() || !j.at("vacancy").get<bool>()) fail(p + ".vacancy", "must be true");
            if (j.contains("species") || j.contains("phase")) fail(p, "vacancy cannot also define species or phase");
            term.label = "*"; term.vacancy = true;
        } else {
            auto name = text(required(j, "species", p), p + ".species");
            if (!species.contains(name)) fail(p + ".species", "undefined species '" + name + "'");
            term.phase = text(required(j, "phase", p), p + ".phase");
            if (term.phase != "gas" && term.phase != "surface") fail(p + ".phase", "expected gas or surface");
            term.label = name + (term.phase == "surface" ? "*" : "");
        }
        if (!seen.insert(term.label).second) fail(p, "duplicate term; use coefficient instead");
        result.push_back(term);
    }
    return result;
}
std::string equationSide(const std::vector<Term>& terms) {
    std::string out;
    for (const auto& term : terms) {
        if (!out.empty()) out += " + ";
        if (term.coefficient != 1) out += std::to_string(term.coefficient);
        out += term.label;
    }
    return out;
}
void validateTopology(const std::string& type, const std::vector<Term>& r,
                      const std::vector<Term>& p, const std::string& path) {
    auto surface = [](const Term& t) { return !t.vacancy && t.phase == "surface"; };
    if (type == "adsorption") {
        if (r.size() != 2 || p.size() != 1 || !surface(p[0])) fail(path, "adsorption requires gas + vacancy -> surface species");
        const Term* gas = nullptr; const Term* vacancy = nullptr;
        for (const auto& t : r) { if (t.vacancy) vacancy = &t; else if (t.phase == "gas") gas = &t; }
        if (!gas || !vacancy || gas->coefficient != 1 || vacancy->coefficient != p[0].coefficient || p[0].coefficient > 2)
            fail(path, "supported adsorption requires one gas molecule and one or two matching surface sites");
    } else if (type == "desorption") {
        if (r.size() != 1 || !surface(r[0]) || r[0].coefficient != 1 || p.size() != 2)
            fail(path, "desorption requires surface species -> vacancy + gas");
        int vacancies = 0, gases = 0;
        for (const auto& t : p) { if (t.coefficient != 1) fail(path, "desorption coefficients must be one");
            if (t.vacancy) ++vacancies; else if (t.phase == "gas") ++gases; }
        if (vacancies != 1 || gases != 1) fail(path, "desorption requires vacancy + gas products");
    } else if (type == "diffusion") {
        if (r.size() != 1 || p.size() != 1 || !surface(r[0]) || !surface(p[0]) ||
            r[0].label != p[0].label || r[0].coefficient != 1 || p[0].coefficient != 1)
            fail(path, "diffusion requires the same single surface species on both sides");
    } else if (type == "reaction") {
        if (r.size() != 2 || p.size() > 2) fail(path, "current reactions support two reactant terms and at most two product terms");
        for (const auto& t : r) if (!surface(t) || t.coefficient != 1) fail(path, "reaction reactants must be single surface species");
        for (const auto& t : p) if (t.vacancy || t.coefficient != 1) fail(path, "reaction products must be species with coefficient one");
    } else fail(path + ".type", "unsupported process type '" + type + "'");
}
std::string elementFor(const Json& species,const std::string& name) {
    const std::string path="species."+name;
    if (!species.contains(name)) fail(path,"undefined species");
    const auto& parameters=required(species.at(name),"parameters",path);
    const auto element=text(required(parameters,"element",path+".parameters"),path+".parameters.element");
    const std::string elements=" H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og ";
    if (elements.find(" "+element+" ")==std::string::npos) fail(path+".parameters.element","expected a chemical element symbol");
    return element;
}
Atomistic::RunOptions atomisticRun(const Json& root) {
    Atomistic::RunOptions o;
    const auto& system=required(root,"system","input");
    keys(system,"system",{"parameters"});
    const auto& sp=required(system,"parameters","system");
    o.temperature=number(required(sp,"temperature","system.parameters"),"system.parameters.temperature",true);
    if (sp.contains("pressure")) number(sp.at("pressure"),"system.parameters.pressure");
    const auto& sim=required(root,"simulation","input");
    keys(sim,"simulation",{"start_time","duration","random_seed","max_events","parameters"});
    if (sim.contains("start_time")) o.start=number(sim.at("start_time"),"simulation.start_time");
    o.end=o.start+number(required(sim,"duration","simulation"),"simulation.duration",true);
    if (sim.contains("random_seed")) o.seed=integer(sim.at("random_seed"),"simulation.random_seed");
    if (sim.contains("max_events")) o.maxEvents=integer(sim.at("max_events"),"simulation.max_events");
    const auto& output=required(root,"output","input");
    keys(output,"output",{"log_interval","lattice_interval","parameters"});
    o.logInterval=number(required(output,"log_interval","output"),"output.log_interval",true);
    o.snapshotInterval=number(required(output,"lattice_interval","output"),"output.lattice_interval",true);
    const auto& species=required(root,"species","input");
    const auto& processes=required(root,"processes","input");
    if (!processes.is_array() || processes.empty()) fail("processes","expected a nonempty array");
    for (std::size_t i=0;i<processes.size();++i) {
        const auto path="processes["+std::to_string(i)+"]";
        const auto& j=processes[i];
        keys(j,path,{"id","type","description","reactants","products","rate","parameters","all_neighbors"});
        Atomistic::Process p;
        p.id=text(required(j,"id",path),path+".id");
        auto type=text(required(j,"type",path),path+".type");
        if (type!="adsorption" && type!="diffusion" && type!="desorption") fail(path+".type","atomistic backend supports adsorption, diffusion and desorption only");
        auto r=terms(required(j,"reactants",path),species,path+".reactants");
        auto products=terms(required(j,"products",path),species,path+".products");
        validateTopology(type,r,products,path);
        for (const auto& t : r) if (t.coefficient!=1) fail(path,"atomistic processes require coefficient one");
        for (const auto& t : products) if (t.coefficient!=1) fail(path,"atomistic processes require coefficient one");
        const auto label=type=="adsorption" ? products[0].label : r[0].label;
        p.species=label.substr(0,label.size()-1);
        p.element=elementFor(species,p.species);
        p.kind=type=="adsorption" ? Atomistic::Kind::Adsorption : (type=="diffusion" ? Atomistic::Kind::Diffusion : Atomistic::Kind::Desorption);
        p.minCoordination=type=="desorption" ? 0 : 1;
        for (const auto& t : r) if (!t.vacancy && t.phase=="gas" && elementFor(species,t.label)!=p.element) fail(path,"adsorption must conserve the chemical element");
        for (const auto& t : products) if (!t.vacancy && t.phase=="gas" && elementFor(species,t.label)!=p.element) fail(path,"desorption must conserve the chemical element");
        if (j.contains("all_neighbors") && (!j.at("all_neighbors").is_boolean() || j.at("all_neighbors").get<bool>()))
            fail(path+".all_neighbors","atomistic events enumerate individual sites/pairs; all_neighbors must be false or omitted");
        if (j.contains("parameters")) {
            const auto& cp=j.at("parameters");
            if (cp.contains("min_coordination")) p.minCoordination=integer(cp.at("min_coordination"),path+".parameters.min_coordination");
            if (cp.contains("max_coordination")) p.maxCoordination=integer(cp.at("max_coordination"),path+".parameters.max_coordination");
        }
        const auto& rate=required(j,"rate",path);
        keys(rate,path+".rate",{"model","parameters"});
        p.model=text(required(rate,"model",path+".rate"),path+".rate.model");
        const auto& rp=required(rate,"parameters",path+".rate");
        object(rp,path+".rate.parameters");
        auto get=[&](const std::string& key){return number(required(rp,key,path+".rate.parameters"),path+".rate.parameters."+key);};
        if (p.model=="constant") p.value=get("value");
        else if (p.model=="arrhenius") {
            p.prefactor=get("prefactor");p.activationEnergy=get("activation_energy");
            if (rp.contains("bond_energy")) p.bondEnergy=get("bond_energy");
            if (rp.contains("migration_energy")) fail(path+".rate.parameters.migration_energy","use activation_energy for the atomistic hop barrier");
        } else fail(path+".rate.model","atomistic rates support constant or arrhenius");
        o.processes.push_back(std::move(p));
    }
    Atomistic::validate(o);
    return o;
}
Atomistic::Options surfaceOptions(const Json& root) {
    const auto& lattice = required(root, "lattice", "input");
    keys(lattice, "lattice", {"representation", "type", "parameters"});
    const auto& p = required(lattice, "parameters", "lattice");
    object(p, "lattice.parameters");
    const std::string path = "lattice.parameters";
    Atomistic::Options o;
    o.crystal = text(required(lattice, "type", "lattice"), "lattice.type");
    o.latticeConstant = number(required(p, "lattice_constant", path), path + ".lattice_constant", true);
    if (p.contains("length_unit") && text(p.at("length_unit"), path + ".length_unit") != "angstrom")
        fail(path + ".length_unit", "only angstrom is supported");
    if (p.contains("orientation")) {
        auto orientation = text(p.at("orientation"), path + ".orientation");
        if (orientation != (o.crystal == "HCP" ? "0001" : "100"))
            fail(path + ".orientation", "supported orientations are cubic 100 and HCP 0001");
    }
    const auto& repeats = required(p, "repetitions", path);
    if (!repeats.is_array() || repeats.size() != 3) fail(path + ".repetitions", "expected [nx, ny, nz] unit cells");
    for (int axis=0;axis<3;++axis) o.repetitions[axis] = integer(repeats[axis], path + ".repetitions", true);
    if (p.contains("periodic")) {
        const auto& periodic = p.at("periodic");
        if (!periodic.is_array() || periodic.size() != 3) fail(path + ".periodic", "expected three booleans");
        for (int axis=0;axis<3;++axis) {
            if (!periodic[axis].is_boolean()) fail(path + ".periodic", "expected three booleans");
            o.periodic[axis] = periodic[axis].get<bool>();
        }
    }
    if (p.contains("origin")) {
        const auto& origin = p.at("origin");
        if (!origin.is_array() || origin.size() != 3) fail(path + ".origin", "expected [x, y, z]");
        for (int axis=0;axis<3;++axis) {
            if (!origin[axis].is_number() || !std::isfinite(origin[axis].get<double>())) fail(path + ".origin", "coordinates must be finite numbers");
            o.origin[axis] = origin[axis].get<double>();
        }
    }
    if (p.contains("vacuum")) o.vacuum = number(p.at("vacuum"), path + ".vacuum");
    if (p.contains("c_over_a")) o.cOverA = number(p.at("c_over_a"), path + ".c_over_a", true);
    if (p.contains("neighbor_cutoff")) o.neighborCutoff = number(p.at("neighbor_cutoff"), path + ".neighbor_cutoff", true);
    if (p.contains("occupied_unit_cells")) o.occupiedUnitCells = integer(p.at("occupied_unit_cells"), path + ".occupied_unit_cells");
    if (p.contains("fixed_unit_cells")) o.fixedUnitCells = integer(p.at("fixed_unit_cells"), path + ".fixed_unit_cells");
    if (p.contains("vacancies")) {
        if (!p.at("vacancies").is_array()) fail(path + ".vacancies", "expected an array of site IDs");
        for (const auto& id : p.at("vacancies")) o.vacancies.push_back(static_cast<std::size_t>(integer(id, path + ".vacancies")));
    }
    for (const auto& field : {"dimensions", "initial_height", "initial_height_file", "initial_species_file", "growth_species", "steps"})
        if (p.contains(field)) fail(path + "." + field, "height-field settings do not apply to atomistic surfaces");
    o.species = text(required(p, "initial_species", path), path + ".initial_species");
    const auto& species = required(root, "species", "input");
    object(species, "species");
    for (auto it=species.begin();it!=species.end();++it) keys(it.value(), "species." + it.key(), {"description", "parameters"});
    if (!species.contains(o.species)) fail(path + ".initial_species", "undefined species");
    o.element = elementFor(species,o.species);
    Atomistic::validate(o);
    return o;
}
} // namespace

double Value::number() const {
    if (std::holds_alternative<std::int64_t>(data)) return static_cast<double>(as<std::int64_t>());
    return as<double>();
}
bool Value::contains(const std::string& key) const {
    return std::holds_alternative<Object>(data) && as<Object>().count(key) != 0;
}

Configuration readJson(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) fail(filename, "cannot open input file");
    Json root;
    try {
        // Reject duplicate object keys instead of silently replacing their values.
        std::vector<std::set<std::string>> keyStack;
        root = Json::parse(file, [&keyStack](int, Json::parse_event_t event, Json& parsed) {
            if (event == Json::parse_event_t::object_start) keyStack.emplace_back();
            else if (event == Json::parse_event_t::object_end) keyStack.pop_back();
            else if (event == Json::parse_event_t::key && !keyStack.back().insert(parsed.get<std::string>()).second)
                fail("input", "duplicate key '" + parsed.get<std::string>() + "'");
            return true;
        });
    } catch (const Json::exception& e) { fail(filename, e.what()); }
    keys(root, "input", {"schema_version", "description", "parameters", "system", "simulation", "lattice", "species", "processes", "output"});
    if (integer(required(root, "schema_version", "input"), "schema_version") != 1) fail("schema_version", "only version 1 is supported");
    Configuration c;
    const auto& latticeDefinition = required(root, "lattice", "input");
    object(latticeDefinition, "lattice");
    if (latticeDefinition.contains("representation")) {
        const auto representation = text(latticeDefinition.at("representation"), "lattice.representation");
        if (representation == "atomistic") {
            c.atomisticSurface = surfaceOptions(root);
            if (root.contains("system") || root.contains("simulation") || root.contains("processes") || root.contains("output")) {
                if (c.atomisticSurface->periodic[2]) fail("lattice.parameters.periodic","atomistic gas access requires nonperiodic z");
                c.atomisticSimulation = atomisticRun(root);
            }
            c.document = convert(root);
            return c;
        }
        if (representation != "height_field") fail("lattice.representation", "expected height_field or atomistic");
    }
    const auto& system = required(root, "system", "input");
    keys(system, "system", {"parameters"});
    const auto& sys = required(system, "parameters", "system");
    c.temperature = number(required(sys, "temperature", "system.parameters"), "system.parameters.temperature", true);
    c.pressure = number(required(sys, "pressure", "system.parameters"), "system.parameters.pressure");
    const auto& sim = required(root, "simulation", "input");
    keys(sim, "simulation", {"start_time", "duration", "random_seed", "parameters"});
    if (sim.contains("start_time")) c.startTime = number(sim.at("start_time"), "simulation.start_time");
    c.endTime = c.startTime + number(required(sim, "duration", "simulation"), "simulation.duration", true);
    if (!std::isfinite(c.endTime) || c.endTime <= c.startTime) fail("simulation.duration", "end time is out of range");
    if (sim.contains("random_seed")) c.seed = integer(sim.at("random_seed"), "simulation.random_seed");
    const auto& species = required(root, "species", "input");
    object(species, "species");
    if (species.empty()) fail("species", "at least one species must be declared");
    for (auto it = species.begin(); it != species.end(); ++it) {
        const auto& name = it.key();
        if (name.empty() || name.find_first_of(" \t\r\n*+:->") != std::string::npos || (std::isdigit(static_cast<unsigned char>(name[0])) || name[0] == '.'))
            fail("species." + name, "name must not start with a digit or dot or contain whitespace or equation delimiters");
        keys(it.value(), "species." + name, {"description", "parameters"});
    }
    const auto& lattice = required(root, "lattice", "input");
    keys(lattice, "lattice", {"representation", "type", "parameters"});
    c.latticeType = text(required(lattice, "type", "lattice"), "lattice.type");
    if (c.latticeType != "SimpleCubic" && c.latticeType != "FCC" && c.latticeType != "HCP" && c.latticeType != "Diamond") fail("lattice.type", "unsupported lattice type");
    const auto& lp = required(lattice, "parameters", "lattice");
    object(lp, "lattice.parameters");
    const auto& dims = required(lp, "dimensions", "lattice.parameters");
    if (!dims.is_array() || dims.size() != 2) fail("lattice.parameters.dimensions", "expected [x, y]");
    c.x = integer(dims[0], "lattice.parameters.dimensions[0]", true);
    c.y = integer(dims[1], "lattice.parameters.dimensions[1]", true);
    if (c.x < 3 || c.y < 3 || static_cast<long long>(c.x) * c.y > std::numeric_limits<int>::max())
        fail("lattice.parameters.dimensions", "dimensions must be at least 3 and their product must fit an int");
    if ((c.latticeType == "HCP" && c.x % 2 != 0) ||
        (c.latticeType == "FCC" && (c.x % 2 != 0 || c.y % 2 != 0)))
        fail("lattice.parameters.dimensions", "HCP requires even x; FCC requires even x and y");
    c.initialSpecies = text(required(lp, "initial_species", "lattice.parameters"), "lattice.parameters.initial_species");
    if (!species.contains(c.initialSpecies)) fail("lattice.parameters.initial_species", "undefined species");
    auto path = [&](const Json& value, const std::string& field) {
        return (std::filesystem::path(filename).parent_path() / text(value, field)).lexically_normal().string();
    };
    if (lp.contains("initial_height_file")) {
        if (lp.contains("initial_height")) fail("lattice.parameters", "choose initial_height or initial_height_file");
        c.heightsFile = path(lp.at("initial_height_file"), "lattice.parameters.initial_height_file");
    } else c.height = integer(required(lp, "initial_height", "lattice.parameters"), "lattice.parameters.initial_height");
    if (lp.contains("initial_species_file")) c.speciesFile = path(lp.at("initial_species_file"), "lattice.parameters.initial_species_file");
    if ((!c.heightsFile.empty() || !c.speciesFile.empty()) && c.latticeType != "SimpleCubic") fail("lattice.parameters", "matrix files are supported only by SimpleCubic");
    if (lp.contains("steps")) fail("lattice.parameters.steps", "steps are not supported by JSON version 1");
    auto labels = [&](const Json& values, const std::string& field, bool growth) {
        if (!values.is_array()) fail(field, "expected an array");
        std::vector<std::string> result;
        for (const auto& v : values) {
            auto label = text(v, field);
            auto name = label;
            if (!growth && name.back() == '*') name.pop_back();
            if (!species.contains(name)) fail(field, "undefined species '" + name + "'");
            result.push_back(label + (growth ? "*" : ""));
        }
        return result;
    };
    if (lp.contains("growth_species")) c.growth = labels(lp.at("growth_species"), "lattice.parameters.growth_species", true);
    const auto& output = required(root, "output", "input");
    keys(output, "output", {"log_interval", "lattice_interval", "coverage", "parameters"});
    c.logInterval = number(required(output, "log_interval", "output"), "output.log_interval", true);
    c.latticeInterval = number(required(output, "lattice_interval", "output"), "output.lattice_interval", true);
    if (output.contains("coverage")) c.coverage = labels(output.at("coverage"), "output.coverage", false);
    const auto& processes = required(root, "processes", "input");
    if (!processes.is_array() || processes.empty()) fail("processes", "expected a nonempty array");
    std::set<std::string> ids;
    for (std::size_t i = 0; i < processes.size(); ++i) {
        const auto& proc = processes[i];
        auto p = "processes[" + std::to_string(i) + "]";
        keys(proc, p, {"id", "type", "description", "reactants", "products", "rate", "parameters", "all_neighbors"});
        Process model;
        model.id = text(required(proc, "id", p), p + ".id");
        if (!ids.insert(model.id).second) fail(p + ".id", "duplicate process id");
        model.type = text(required(proc, "type", p), p + ".type");
        auto r = terms(required(proc, "reactants", p), species, p + ".reactants");
        auto products = terms(required(proc, "products", p), species, p + ".products");
        validateTopology(model.type, r, products, p);
        model.equation = equationSide(r) + " -> " + equationSide(products);
        const auto& rate = required(proc, "rate", p);
        keys(rate, p + ".rate", {"model", "parameters"});
        auto rateModel = text(required(rate, "model", p + ".rate"), p + ".rate.model");
        const auto& rp = required(rate, "parameters", p + ".rate");
        object(rp, p + ".rate.parameters");
        model.modelArguments.push_back(rateModel);
        auto argument = [&](const std::string& key, bool positive = false) {
            auto field = p + ".rate.parameters." + key;
            double n = number(required(rp, key, p + ".rate.parameters"), field, positive);
            model.modelArguments.push_back(numericArgument(n)); return n;
        };
        if (rateModel == "constant") argument("value");
        else if (rateModel == "simple" && model.type == "adsorption") {
            if (argument("sticking_coefficient") > 1) fail(p + ".rate.parameters.sticking_coefficient", "must be at most 1");
            if (argument("molar_fraction") > 1) fail(p + ".rate.parameters.molar_fraction", "must be at most 1");
            argument("site_density", true);
            // Resolve molecular weight from the gas species unless explicitly overridden.
            if (rp.contains("molecular_weight")) argument("molecular_weight", true);
            else {
                std::string gas;
                for (const auto& t : r) if (!t.vacancy && t.phase == "gas") gas = t.label;
                const auto& sp = required(species.at(gas), "parameters", "species." + gas);
                model.modelArguments.push_back(numericArgument(number(required(sp, "molecular_weight", "species." + gas + ".parameters"), "species." + gas + ".parameters.molecular_weight", true)));
            }
        } else if (rateModel == "arrhenius" && model.type != "adsorption") {
            argument("prefactor"); argument("activation_energy");
            if (model.type == "diffusion") argument("migration_energy");
        } else fail(p + ".rate.model", "unsupported rate model for " + model.type);
        bool all = false;
        if (proc.contains("all_neighbors")) {
            if (!proc.at("all_neighbors").is_boolean()) fail(p + ".all_neighbors", "expected a boolean");
            all = proc.at("all_neighbors").get<bool>();
        }
        if (all && model.type == "reaction") fail(p + ".all_neighbors", "not supported for reactions");
        if (all) model.modelArguments.push_back("all");
        c.processes.push_back(std::move(model));
    }
    auto validateMatrix = [&](const std::string& filename, bool heights) {
        if (filename.empty()) return;
        std::ifstream matrix(filename);
        if (!matrix) fail(filename, "cannot open lattice matrix");
        std::string line;
        int rows = 0;
        while (std::getline(matrix, line)) {
            std::istringstream row(line);
            std::string token;
            int columns = 0;
            while (row >> token) {
                if (heights) {
                    std::size_t end = 0;
                    long long h = 0;
                    try { h = std::stoll(token, &end); }
                    catch (const std::exception&) { fail(filename, "invalid height '" + token + "'"); }
                    if (end != token.size() || h < 0 || h > std::numeric_limits<int>::max()) fail(filename, "height out of range");
                } else {
                    auto name = token;
                    if (!name.empty() && name.back() == '*') name.pop_back();
                    if (!species.contains(name)) fail(filename, "undefined species '" + name + "'");
                }
                ++columns;
            }
            if (!columns) continue;
            if (columns != c.x) fail(filename, "each row must contain exactly " + std::to_string(c.x) + " values");
            ++rows;
        }
        if (rows != c.y) fail(filename, "expected " + std::to_string(c.y) + " rows");
    };
    validateMatrix(c.heightsFile, true);
    validateMatrix(c.speciesFile, false);
    // Make optional parameter bags consistently available to model code.
    auto addParameters = [](Json& entity) {
        if (!entity.contains("parameters")) entity["parameters"] = Json::object();
    };
    addParameters(root);
    addParameters(root["simulation"]);
    addParameters(root["output"]);
    for (auto& entity : root["species"]) addParameters(entity);
    for (auto& entity : root["processes"]) addParameters(entity);
    c.document = convert(root);
    return c;
}
} // namespace Config
