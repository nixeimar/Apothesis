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
std::string elementFor(const Json& species,const std::string& name,bool allowPseudo=false) {
    const std::string path="species."+name;
    if (!species.contains(name)) fail(path,"undefined species");
    const auto& parameters=required(species.at(name),"parameters",path);
    const auto element=text(required(parameters,"element",path+".parameters"),path+".parameters.element");
    const std::string elements=" H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og ";
    if (!(allowPseudo && element=="X") && elements.find(" "+element+" ")==std::string::npos) fail(path+".parameters.element","expected a chemical element symbol");
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
Atomistic::Options surfaceOptions(const Json& root,bool allowPseudo=false) {
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
    o.element = elementFor(species,o.species,allowPseudo);
    Atomistic::validate(o);
    return o;
}

ZGB::Options zgbOptions(const Json& root,bool volume) {
    const auto& model=required(root,"model","input");
    keys(model,"model",{"type","parameters","description"});
    const auto& p=required(model,"parameters","model");
    object(p,"model.parameters");
    ZGB::Options o;
    if (p.contains("variant")) {
        auto variant=text(p.at("variant"),"model.parameters.variant");
        if (variant!="growth" && variant!="catalytic") fail("model.parameters.variant","expected growth or catalytic");
        o.growth=variant=="growth";
    }
    o.yAB=number(required(p,"y_ab","model.parameters"),"model.parameters.y_ab");
    if (p.contains("reaction_rate")) o.reactionRate=number(p.at("reaction_rate"),"model.parameters.reaction_rate",true);
    if (p.contains("pair_convention")) o.pairConvention=text(p.at("pair_convention"),"model.parameters.pair_convention");
    if (p.contains("reaction_selection")) o.reactionSelection=text(p.at("reaction_selection"),"model.parameters.reaction_selection");
    const auto& lattice=required(root,"lattice","input");
    keys(lattice,"lattice",{"representation","type","parameters"});
    if (text(required(lattice,"type","lattice"),"lattice.type")!="SimpleCubic") fail("lattice.type","ZGB requires SimpleCubic");
    const auto& lp=required(lattice,"parameters","lattice");object(lp,"lattice.parameters");
    if (!volume) {
        if (lattice.contains("representation") && text(lattice.at("representation"),"lattice.representation")!="height_field")
            fail("lattice.representation","use cheimarios_zgb_2023_3d for atomistic ZGB");
        const auto& dims=required(lp,"dimensions","lattice.parameters");
        if (!dims.is_array() || dims.size()!=2) fail("lattice.parameters.dimensions","expected [nx, ny]");
        o.x=integer(dims[0],"lattice.parameters.dimensions",true);o.y=integer(dims[1],"lattice.parameters.dimensions",true);
        o.initialHeight=integer(required(lp,"initial_height","lattice.parameters"),"lattice.parameters.initial_height");
        for (const auto& field:{"initial_height_file","initial_species_file","steps","repetitions","occupied_unit_cells"})
            if (lp.contains(field)) fail(std::string("lattice.parameters.")+field,"not supported by the predefined height-field ZGB model");
    } else {
        if (text(required(lattice,"representation","lattice"),"lattice.representation")!="atomistic") fail("lattice.representation","3D ZGB requires atomistic");
        auto geometry=surfaceOptions(root,true);
        o.x=geometry.repetitions[0];o.y=geometry.repetitions[1];
        if (geometry.occupiedUnitCells<1 || geometry.occupiedUnitCells>=geometry.repetitions[2])
            fail("lattice.parameters.occupied_unit_cells","3D ZGB requires occupied layers and empty growth headroom");
        o.initialHeight=geometry.occupiedUnitCells;
    }
    const auto& sim=required(root,"simulation","input");
    keys(sim,"simulation",{"start_time","duration","random_seed","max_events","parameters"});
    if (sim.contains("start_time")) o.start=number(sim.at("start_time"),"simulation.start_time");
    o.end=o.start+number(required(sim,"duration","simulation"),"simulation.duration",true);
    if (sim.contains("random_seed")) o.seed=integer(sim.at("random_seed"),"simulation.random_seed");
    if (sim.contains("max_events")) o.maxEvents=integer(sim.at("max_events"),"simulation.max_events");
    const auto& output=required(root,"output","input");
    if (volume) keys(output,"output",{"log_interval","lattice_interval","write_xyz","parameters"});
    else keys(output,"output",{"log_interval","parameters"});
    o.logInterval=number(required(output,"log_interval","output"),"output.log_interval",true);
    if (root.contains("processes") || root.contains("system")) fail("model","predefined ZGB kinetics do not accept system or processes sections");
    ZGB::validate(o);return o;
}
ZGB3D::Options zgb3dOptions(const Json& root) {
    ZGB3D::Options o;o.kinetics=zgbOptions(root,true);o.geometry=surfaceOptions(root,true);
    const auto& output=root.at("output");
    if (output.contains("lattice_interval")) o.snapshotInterval=number(output.at("lattice_interval"),"output.lattice_interval");
    if (output.contains("write_xyz")) {
        if (!output.at("write_xyz").is_boolean()) fail("output.write_xyz","expected a boolean");
        o.writeXYZ=output.at("write_xyz").get<bool>();
    }
    ZGB3D::validate(o);return o;
}


Atomistic::Surface customSurface(const Json& root,const std::map<std::string,std::string>& elements,double& siteVolume) {
    const auto& lattice=root.at("lattice");keys(lattice,"lattice",{"representation","type","parameters"});
    const auto& p=required(lattice,"parameters","lattice");object(p,"lattice.parameters");
    const auto path=std::string("lattice.parameters");Atomistic::Surface s;s.elements=elements;
    for(const auto& field:{"repetitions","vacancies","fixed_unit_cells","occupied_unit_cells","origin","lattice_constant","orientation","c_over_a","vacuum"}) if(p.contains(field)) fail(path+"."+field,"Custom geometry uses explicit positions, occupancy and cell");
    s.neighborCutoff=number(required(p,"neighbor_cutoff",path),path+".neighbor_cutoff",true);
    siteVolume=number(required(p,"site_volume",path),path+".site_volume",true);
    auto name=text(required(p,"initial_species",path),path+".initial_species");
    if(!elements.count(name)) fail(path+".initial_species","undefined species");s.element=elements.at(name);
    const auto& cell=required(p,"cell",path);
    if(!cell.is_array() || cell.size()!=3) fail(path+".cell","Custom geometry requires [Lx, Ly, Lz] orthorhombic box lengths");
    s.cell={{{{0,0,0}},{{0,0,0}},{{0,0,0}}}};
    for(int axis=0;axis<3;++axis) s.cell[axis][axis]=number(cell[axis],path+".cell",true);
    s.periodic={{false,false,false}};
    if(p.contains("periodic")) {
        const auto& periodic=p.at("periodic");if(!periodic.is_array() || periodic.size()!=3) fail(path+".periodic","expected three booleans");
        for(int axis=0;axis<3;++axis) {if(!periodic[axis].is_boolean()) fail(path+".periodic","expected three booleans");s.periodic[axis]=periodic[axis].get<bool>();}
    }
    if(p.contains("length_unit") && text(p.at("length_unit"),path+".length_unit")!="angstrom") fail(path+".length_unit","only angstrom is supported");
    const auto& sites=required(p,"sites",path);
    if(!sites.is_array() || sites.empty() || sites.size()>1000000) fail(path+".sites","expected 1 to 1000000 explicit candidate sites");
    std::set<Atomistic::Vec3> unique;
    for(const auto& j:sites) {
        keys(j,path+".sites",{"position","species","occupied","fixed"});
        const auto& xyz=required(j,"position",path+".sites");if(!xyz.is_array() || xyz.size()!=3) fail(path+".sites.position","expected [x,y,z]");
        Atomistic::Vec3 position;
        for(int axis=0;axis<3;++axis) {
            position[axis]=number(xyz[axis],path+".sites.position");
            if(position[axis]>s.cell[axis][axis] || (s.periodic[axis] && position[axis]>=s.cell[axis][axis])) fail(path+".sites.position","coordinates must lie inside the custom box (periodic upper faces are excluded)");
        }
        if(!unique.insert(position).second) fail(path+".sites.position","duplicate coordinates");
        auto species=j.contains("species") ? text(j.at("species"),path+".sites.species") : name;
        if(!elements.count(species)) fail(path+".sites.species","undefined species");
        const auto& occupied=required(j,"occupied",path+".sites");if(!occupied.is_boolean()) fail(path+".sites.occupied","expected a boolean");
        bool fixed=false;if(j.contains("fixed")) {if(!j.at("fixed").is_boolean()) fail(path+".sites.fixed","expected a boolean");fixed=j.at("fixed").get<bool>();}
        if(fixed && !occupied.get<bool>()) fail(path+".sites.fixed","fixed sites must be occupied");
        s.sites.push_back({s.sites.size(),position,species,occupied.get<bool>(),fixed,{}});
    }
    auto neighbors=Atomistic::neighborGraph(s,0,s.neighborCutoff);
    for(auto& site:s.sites) {
        site.neighbors=std::move(neighbors[site.id]);
        for(auto j:site.neighbors) {const auto d=Atomistic::displacement(s,site.position,s.sites[j].position);if(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]<1e-20) fail(path+".sites.position","coincident periodic sites");}
    }
    return s;
}
SurfaceKMC::Options explicitSurfaceOptions(const Json& root) {
    const auto& model=root.at("model");keys(model,"model",{"type","parameters","description"});
    const auto& lattice=required(root,"lattice","input");
    if (text(required(lattice,"representation","lattice"),"lattice.representation")!="atomistic") fail("lattice.representation","surface_kmc_3d requires atomistic");
    SurfaceKMC::Options o;
    o.surfaceSites="gas_exposed";
    if(model.contains("parameters")) {
        const auto& p=model.at("parameters");
        if(p.contains("surface_sites")) o.surfaceSites=text(p.at("surface_sites"),"model.parameters.surface_sites");
        if(p.contains("max_embeddings")) o.maxEmbeddings=integer(p.at("max_embeddings"),"model.parameters.max_embeddings",true);
        if(p.contains("max_search_steps")) o.maxSearchSteps=integer(p.at("max_search_steps"),"model.parameters.max_search_steps",true);
    }
    const auto& species=required(root,"species","input");object(species,"species");
    for (auto it=species.begin();it!=species.end();++it) {
        if (it.key().empty() || it.key().find_first_of(" \t\n\r\"")!=std::string::npos) fail("species","names must be nonempty tokens without quotes or whitespace");
        keys(it.value(),"species."+it.key(),{"description","parameters"});
        o.elements[it.key()]=elementFor(species,it.key(),true);
    }
    if(text(required(lattice,"type","lattice"),"lattice.type")=="Custom") {
        o.customGeometry=customSurface(root,o.elements,o.siteVolume);o.geometry.crystal="Custom";
    } else o.geometry=surfaceOptions(root,true);
    double temperature=0;
    if (root.contains("system")) {
        const auto& system=root.at("system");keys(system,"system",{"parameters"});
        const auto& p=required(system,"parameters","system");
        if (p.contains("temperature")) temperature=number(p.at("temperature"),"system.parameters.temperature",true);
    }
    const auto& sim=required(root,"simulation","input");keys(sim,"simulation",{"start_time","duration","random_seed","max_events","parameters"});
    if (sim.contains("start_time")) o.start=number(sim.at("start_time"),"simulation.start_time");
    o.end=o.start+number(required(sim,"duration","simulation"),"simulation.duration",true);
    if (sim.contains("random_seed")) o.seed=integer(sim.at("random_seed"),"simulation.random_seed");
    if (sim.contains("max_events")) o.maxEvents=integer(sim.at("max_events"),"simulation.max_events");
    const auto& output=required(root,"output","input");keys(output,"output",{"log_interval","lattice_interval","write_xyz","parameters"});
    o.logInterval=number(required(output,"log_interval","output"),"output.log_interval",true);
    if (output.contains("lattice_interval")) o.snapshotInterval=number(output.at("lattice_interval"),"output.lattice_interval");
    if (output.contains("write_xyz")) {if (!output.at("write_xyz").is_boolean()) fail("output.write_xyz","expected a boolean");o.writeXYZ=output.at("write_xyz").get<bool>();}
    const auto& processes=required(root,"processes","input");
    if (!processes.is_array() || processes.empty()) fail("processes","expected a nonempty array");
    for (std::size_t index=0;index<processes.size();++index) {
        const auto path="processes["+std::to_string(index)+"]";const auto& j=processes[index];
        keys(j,path,{"id","type","description","reactants","products","rate","parameters"});
        SurfaceKMC::Rule r;r.id=text(required(j,"id",path),path+".id");r.type=text(required(j,"type",path),path+".type");
        if (r.id.find_first_of(" \t\r\n\"")!=std::string::npos) fail(path+".id","expected a token without whitespace or quotes");
        std::map<int,std::string> before;std::map<int,bool> beforeSolid;std::map<int,SurfaceKMC::Product> after;
        int gasReactants=0,gasProducts=0;
        for (bool products:{false,true}) {
            const std::string side=products ? "products" : "reactants";const auto& list=required(j,side,path);
            if (!list.is_array() || list.empty()) fail(path+"."+side,"expected a nonempty array");
            for (const auto& term:list) {
                const auto tp=path+"."+side;keys(term,tp,{"species","phase","site","vacancy","coefficient"});
                if (term.contains("coefficient") && integer(term.at("coefficient"),tp+".coefficient",true)!=1) fail(tp+".coefficient","write one term per site slot; coefficients must be one");
                bool vacant=term.contains("vacancy");std::string name,phase;
                if (vacant) {
                    if (!term.at("vacancy").is_boolean() || !term.at("vacancy").get<bool>() || term.contains("species") || term.contains("phase")) fail(tp,"vacancy must be true without species or phase");
                } else {
                    name=text(required(term,"species",tp),tp+".species");
                    if (!o.elements.count(name)) fail(tp+".species","undefined species '"+name+"'");
                    phase=text(required(term,"phase",tp),tp+".phase");
                    if (phase!="gas" && phase!="surface" && phase!="solid") fail(tp+".phase","expected gas, surface or solid");
                    if (phase=="gas") {
                        if (term.contains("site")) fail(tp+".site","gas reservoirs do not occupy a site slot");
                        if (products) ++gasProducts;else ++gasReactants;
                        continue;
                    }
                    if (!products && phase=="solid" && o.surfaceSites=="column_tops") fail(tp+".phase","solid reactants require gas_exposed mode");
                }
                const int slot=integer(required(term,"site",tp),tp+".site");
                if (products) {if (!after.emplace(slot,SurfaceKMC::Product{name,phase=="solid"}).second) fail(tp+".site","duplicate site slot");}
                else {if (!before.emplace(slot,name).second) fail(tp+".site","duplicate site slot");beforeSolid[slot]=phase=="solid";}
            }
        }
        if(before.empty() || before.size()!=after.size()) fail(path,"reactants and products must define the same site slots");
        for(std::size_t slot=0;slot<before.size();++slot) {
            if(!before.count(static_cast<int>(slot)) || !after.count(static_cast<int>(slot))) fail(path,"site slots must be contiguous from zero and present on both sides");
            r.before.push_back(before.at(static_cast<int>(slot)));r.beforeSolid.push_back(beforeSolid.at(static_cast<int>(slot)));r.after.push_back(after.at(static_cast<int>(slot)));
        }
        const auto& p=required(j,"parameters",path);object(p,path+".parameters");
        r.counting=text(required(p,"counting",path+".parameters"),path+".parameters.counting");
        auto parseLink=[&](const Json& link,std::size_t first,std::size_t second) {
            SurfaceKMC::Link l;l.first=first;l.second=second;
            if(link.contains("neighbors")) l.neighbors=text(link.at("neighbors"),path+".parameters.neighbors");
            if(link.contains("same_height")) {if(!link.at("same_height").is_boolean()) fail(path+".parameters.same_height","expected a boolean");l.sameHeight=link.at("same_height").get<bool>();}
            if(l.neighbors=="distance") {
                if(link.contains("min_distance")) l.minDistance=number(link.at("min_distance"),path+".parameters.min_distance");
                l.maxDistance=number(required(link,"max_distance",path+".parameters"),path+".parameters.max_distance",true);
            } else if(link.contains("min_distance") || link.contains("max_distance")) fail(path+".parameters","distance bounds require neighbors: distance");
            r.links.push_back(l);
        };
        if(p.contains("links")) {
            if(p.contains("neighbors") || p.contains("same_height") || p.contains("min_distance") || p.contains("max_distance")) fail(path+".parameters","put per-link geometry inside links, without competing shorthand conditions");
            const auto& links=p.at("links");if(!links.is_array()) fail(path+".parameters.links","expected an array");
            for(const auto& link:links) {
                keys(link,path+".parameters.links",{"sites","neighbors","same_height","min_distance","max_distance"});
                const auto& ends=required(link,"sites",path+".parameters.links");
                if(!ends.is_array() || ends.size()!=2) fail(path+".parameters.links.sites","expected two slot indices");
                parseLink(link,integer(ends[0],path+".parameters.links.sites"),integer(ends[1],path+".parameters.links.sites"));
            }
        } else if(r.before.size()==2) parseLink(p,0,1);
        else if(r.before.size()>2) fail(path+".parameters.links","multi-site processes require a connected links array");
        else if(p.contains("neighbors") || p.contains("same_height")) fail(path+".parameters","single-site events do not use neighbor conditions");
        if (r.type=="adsorption") {
            if (!gasReactants || gasProducts || std::any_of(r.before.begin(),r.before.end(),[](const std::string& a){return !a.empty();}) || std::any_of(r.after.begin(),r.after.end(),[](const SurfaceKMC::Product& a){return a.species.empty() || a.solid;})) fail(path,"adsorption requires gas + vacant slots -> occupied surface slots");
        } else if (r.type=="desorption") {
            if (gasReactants || !gasProducts || std::any_of(r.before.begin(),r.before.end(),[](const std::string& a){return a.empty();}) || std::any_of(r.after.begin(),r.after.end(),[](const SurfaceKMC::Product& a){return !a.species.empty();})) fail(path,"desorption requires occupied surface slots -> vacant slots + gas");
        } else if (r.type=="diffusion" && (gasReactants || gasProducts)) fail(path,"diffusion cannot exchange gas species");
        const auto& rate=required(j,"rate",path);keys(rate,path+".rate",{"model","parameters"});
        const auto modelName=text(required(rate,"model",path+".rate"),path+".rate.model");
        const auto& rp=required(rate,"parameters",path+".rate");object(rp,path+".rate.parameters");
        if (modelName=="constant") r.rate=number(required(rp,"value",path+".rate.parameters"),path+".rate.parameters.value");
        else if (modelName=="arrhenius") {
            if (!temperature) fail("system.parameters.temperature","required for Arrhenius rates");
            const double prefactor=number(required(rp,"prefactor",path+".rate.parameters"),path+".rate.parameters.prefactor");
            const double barrier=number(required(rp,"activation_energy",path+".rate.parameters"),path+".rate.parameters.activation_energy");
            if (rp.contains("bond_energy") || rp.contains("migration_energy")) fail(path+".rate.parameters","environment-dependent barriers are not implemented by surface_kmc_3d");
            r.rate=prefactor*std::exp(-barrier/(8.31446261815324*temperature));
        } else fail(path+".rate.model","supported models are constant and arrhenius");
        o.processes.push_back(std::move(r));
    }
    SurfaceKMC::validate(o);return o;
}

DiamondCVD::Options diamondOptions(const Json& root) {
    DiamondCVD::Options o;const auto& model=root.at("model");keys(model,"model",{"type","description","parameters"});
    const auto& mp=required(model,"parameters","model");object(mp,"model.parameters");
    const auto& lattice=required(root,"lattice","input");keys(lattice,"lattice",{"representation","type","parameters"});
    if(text(required(lattice,"type","lattice"),"lattice.type")!="Custom" || text(required(lattice,"representation","lattice"),"lattice.representation")!="atomistic") fail("lattice","diamond_cvd_1997 requires explicit oriented Custom carbon candidate coordinates");
    const auto& species=required(root,"species","input");object(species,"species");
    std::map<std::string,std::string> elements;
    for(auto it=species.begin();it!=species.end();++it) {keys(it.value(),"species."+it.key(),{"description","parameters"});elements[it.key()]=elementFor(species,it.key());}
    if(std::none_of(elements.begin(),elements.end(),[](const auto& entry){return entry.second=="H";})) fail("species","diamond CVD requires a declared hydrogen element species");
    o.lattice=customSurface(root,elements,o.siteVolume);
    for(const auto& site:o.lattice.sites) if(elements.at(site.species)!="C") fail("lattice.parameters.sites.species","initial diamond candidate sites must be carbon");
    if(mp.contains("site_volume") && number(mp.at("site_volume"),"model.parameters.site_volume",true)!=o.siteVolume) fail("model.parameters.site_volume","must match the lattice site_volume");
    if(mp.contains("initial_passivation")) {if(!mp.at("initial_passivation").is_boolean()) fail("model.parameters.initial_passivation","expected a boolean");o.passivate=mp.at("initial_passivation").get<bool>();}
    if(mp.contains("ch_bond_length")) o.chLength=number(mp.at("ch_bond_length"),"model.parameters.ch_bond_length",true);
    const auto& system=required(root,"system","input");keys(system,"system",{"parameters"});const auto& sp=required(system,"parameters","system");
    o.temperature=number(required(sp,"temperature","system.parameters"),"system.parameters.temperature",true);
    o.gasTemperature=number(required(sp,"gas_temperature","system.parameters"),"system.parameters.gas_temperature",true);
    const auto& pressures=required(sp,"partial_pressures_torr","system.parameters");object(pressures,"system.parameters.partial_pressures_torr");
    for(auto it=pressures.begin();it!=pressures.end();++it) o.pressuresTorr[it.key()]=number(it.value(),"system.parameters.partial_pressures_torr."+it.key());
    for(const auto& gas:{"H","H2","CH3","C2H2"}) if(!o.pressuresTorr.count(gas)) fail("system.parameters.partial_pressures_torr",std::string("required reservoir ")+gas);
    const auto& sim=required(root,"simulation","input");keys(sim,"simulation",{"start_time","duration","random_seed","max_events","parameters"});
    if(sim.contains("start_time")) o.start=number(sim.at("start_time"),"simulation.start_time");
    o.end=o.start+number(required(sim,"duration","simulation"),"simulation.duration",true);
    if(sim.contains("random_seed")) o.seed=integer(sim.at("random_seed"),"simulation.random_seed");
    if(sim.contains("max_events")) o.maxEvents=integer(sim.at("max_events"),"simulation.max_events");
    const auto& output=required(root,"output","input");keys(output,"output",{"log_interval","lattice_interval","write_xyz","parameters"});
    o.logInterval=number(required(output,"log_interval","output"),"output.log_interval",true);
    if(output.contains("lattice_interval")) o.snapshotInterval=number(output.at("lattice_interval"),"output.lattice_interval");
    if(output.contains("write_xyz")) {if(!output.at("write_xyz").is_boolean()) fail("output.write_xyz","expected a boolean");o.writeXYZ=output.at("write_xyz").get<bool>();}
    const auto& list=required(root,"processes","input");if(!list.is_array() || list.empty()) fail("processes","explicit diamond process channels are required");
    auto signedNumber=[&](const Json& j,const std::string& path) {if(!j.is_number() || !std::isfinite(j.get<double>())) fail(path,"expected a finite number");return j.get<double>();};
    for(std::size_t index=0;index<list.size();++index) {
        const auto path="processes["+std::to_string(index)+"]";const auto& j=list[index];keys(j,path,{"id","type","event","description","parameters","rate"});
        if(text(required(j,"type",path),path+".type")!="reaction") fail(path+".type","covalent channels use reaction type");
        DiamondCVD::Process p;p.id=text(required(j,"id",path),path+".id");p.event=text(required(j,"event",path),path+".event");
        const auto& cp=required(j,"parameters",path);object(cp,path+".parameters");
        if(cp.contains("domain")) {if(!cp.at("domain").is_string()) fail(path+".parameters.domain","expected a string");p.domain=cp.at("domain").get<std::string>();}
        if(cp.contains("hydrogens")) {if(!cp.at("hydrogens").is_number_integer() || cp.at("hydrogens").get<double>()< -1 || cp.at("hydrogens").get<double>()>4) fail(path+".parameters.hydrogens","expected an integer -1..4");p.hydrogens=cp.at("hydrogens").get<int>();}
        const auto& rate=required(j,"rate",path);keys(rate,path+".rate",{"model","parameters"});const auto rateModel=text(required(rate,"model",path+".rate"),path+".rate.model");const auto& rp=required(rate,"parameters",path+".rate");object(rp,path+".rate.parameters");
        if(rateModel=="constant") p.rate=number(required(rp,"value",path+".rate.parameters"),path+".rate.parameters.value");
        else if(rateModel=="thermochemical") {
            auto get=[&](const std::string& name){return signedNumber(required(rp,name,path+".rate.parameters"),path+".rate.parameters."+name);};
            const auto& reverse=required(rp,"reverse",path+".rate.parameters");if(!reverse.is_boolean()) fail(path+".rate.parameters.reverse","expected a boolean");
            p.rate=DiamondCVD::rateCoefficient(get("A"),get("temperature_power"),get("activation_energy_kcal_mol"),get("reaction_enthalpy_kcal_mol"),get("reaction_entropy_cal_mol_K"),get("standard_concentration_mol_cm3"),o.temperature,reverse.get<bool>());
            if(rp.contains("gas_reactant")) {
                if(!rp.at("gas_reactant").is_string()) fail(path+".rate.parameters.gas_reactant","expected a reservoir name or empty string");
                auto gas=rp.at("gas_reactant").get<std::string>();if(!gas.empty()) {if(!o.pressuresTorr.count(gas)) fail(path+".rate.parameters.gas_reactant","undefined gas reservoir");p.rate*=o.pressuresTorr.at(gas)/760/(82.057366080960*o.gasTemperature);}
            }
        } else fail(path+".rate.model","diamond rates support constant or thermochemical");
        o.processes.push_back(std::move(p));
    }
    DiamondCVD::validate(o);return o;
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
    keys(root, "input", {"schema_version", "description", "parameters", "model", "system", "simulation", "lattice", "species", "processes", "output"});
    if (integer(required(root, "schema_version", "input"), "schema_version") != 1) fail("schema_version", "only version 1 is supported");
    Configuration c;
    if (root.contains("model")) {
        const auto& model=root.at("model");
        auto type=text(required(model,"type","model"),"model.type");
        if (type=="diamond_cvd_1997") c.diamondCVD=diamondOptions(root);
        else if (type=="surface_kmc_3d") {c.surfaceKMC=explicitSurfaceOptions(root);if(!c.surfaceKMC->customGeometry) c.atomisticSurface=c.surfaceKMC->geometry;}
        else if (type=="cheimarios_zgb_2023") c.zgb=zgbOptions(root,false);
        else if (type=="cheimarios_zgb_2023_3d") {
            c.zgb3d=zgb3dOptions(root);c.atomisticSurface=c.zgb3d->geometry;
        } else fail("model.type","unsupported predefined model");
        c.document=convert(root);return c;
    }
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
