#include "surface_kmc.h"
#include "../extLibs/nlohmann/json.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>

namespace SurfaceKMC {
void validate(const Options& o) {
    Atomistic::validate(o.geometry);
    const auto& g=o.geometry;
    if (g.crystal!="SimpleCubic" || !g.periodic[0] || !g.periodic[1] || g.periodic[2] ||
        g.occupiedUnitCells<1 || g.occupiedUnitCells>=g.repetitions[2])
        throw std::runtime_error("surface_kmc_3d requires a SimpleCubic slab, lateral PBC, open z and occupied layers with growth headroom");
    if (g.neighborCutoff && (g.neighborCutoff<g.latticeConstant || g.neighborCutoff>=std::sqrt(2.0)*g.latticeConstant))
        throw std::runtime_error("surface_kmc_3d requires the first SC neighbor shell");
    if (!std::isfinite(o.start) || o.start<0 || !std::isfinite(o.end) || o.end<=o.start ||
        !std::isfinite(o.logInterval) || o.logInterval<=0 || !std::isfinite(o.snapshotInterval) || o.snapshotInterval<0)
        throw std::runtime_error("simulation/output: invalid times or intervals");
    if (o.processes.empty()) throw std::runtime_error("processes: at least one explicit process is required");
    std::set<std::string> ids;
    for (const auto& p:o.processes) {
        if (p.id.empty() || !ids.insert(p.id).second) throw std::runtime_error("processes: IDs must be nonempty and unique");
        if (p.type!="adsorption" && p.type!="desorption" && p.type!="diffusion" && p.type!="reaction") throw std::runtime_error("processes.type: unsupported process type");
        if (p.before.empty() || p.before.size()>2 || p.after.size()!=p.before.size()) throw std::runtime_error("processes: use one or two explicit binding-site slots");
        if (p.counting!="sites" && p.counting!="undirected_pairs" && p.counting!="reactive_sites") throw std::runtime_error("processes.parameters.counting: unsupported counting convention");
        if ((p.before.size()==1)!=(p.counting=="sites")) throw std::runtime_error("processes.parameters.counting: sites requires one slot; pair counting requires two");
        if (!std::isfinite(p.rate) || p.rate<0) throw std::runtime_error("processes.rate: must be finite and nonnegative");
        bool changed=false;int solids=0;
        for (std::size_t i=0;i<p.before.size();++i) {
            if (!p.before[i].empty() && !o.elements.count(p.before[i])) throw std::runtime_error("processes.reactants: undefined species");
            if (!p.after[i].species.empty() && !o.elements.count(p.after[i].species)) throw std::runtime_error("processes.products: undefined species");
            if (p.after[i].solid && p.after[i].species.empty()) throw std::runtime_error("processes.products: solid product needs a species");
            solids+=p.after[i].solid;changed|=p.before[i]!=p.after[i].species || p.after[i].solid;
        }
        if (!changed || solids>1) throw std::runtime_error("processes: an event must change state and may incorporate at most one solid unit");
        if (p.type=="diffusion" && (p.before.size()!=2 || solids || p.before[0].empty() || !p.before[1].empty() || !p.after[0].species.empty() || p.after[1].species!=p.before[0]))
            throw std::runtime_error("processes.diffusion: requires surface species at slot 0 and vacancy at slot 1, exchanged in products");
    }
}
void Simulation::Catalog::set(std::size_t id,bool present) {
    auto& p=positions.at(id);
    if (present && p<0) {p=static_cast<long long>(ids.size());ids.push_back(id);}
    else if (!present && p>=0) {auto last=ids.back();ids[static_cast<std::size_t>(p)]=last;positions[last]=p;ids.pop_back();p=-1;}
}
Simulation::Simulation(Options options):m_options(std::move(options)) {
    validate(m_options);m_solid=Atomistic::build(m_options.geometry);m_solid.elements=m_options.elements;
    const auto n=m_solid.sites.size();m_ads.resize(n);m_binding.assign(n,false);m_deposited.assign(n,false);
    m_neighbors.resize(n);m_incident.resize(n);m_edges.resize(2*n);m_up.assign(n,n);m_down.assign(n,n);
    const double a=m_options.geometry.latticeConstant,tol=a*1e-7;
    for (const auto& site:m_solid.sites) {
        std::array<std::size_t,4> neighbors{{n,n,n,n}};
        for (auto j:site.neighbors) {
            Atomistic::Vec3 d;
            for (int k=0;k<3;++k) {d[k]=m_solid.sites[j].position[k]-site.position[k];if (m_solid.periodic[k]) d[k]-=std::round(d[k]/m_solid.cell[k][k])*m_solid.cell[k][k];}
            if (std::abs(d[2])<tol) {
                if (std::abs(d[0]-a)<tol && std::abs(d[1])<tol) neighbors[0]=j;
                else if (std::abs(d[0]+a)<tol && std::abs(d[1])<tol) neighbors[1]=j;
                else if (std::abs(d[1]-a)<tol && std::abs(d[0])<tol) neighbors[2]=j;
                else if (std::abs(d[1]+a)<tol && std::abs(d[0])<tol) neighbors[3]=j;
                else throw std::runtime_error("invalid lateral first neighbor");
            } else if (std::abs(d[2]-a)<tol) m_up[site.id]=j;
            else if (std::abs(d[2]+a)<tol) m_down[site.id]=j;
            else throw std::runtime_error("invalid vertical first neighbor");
        }
        if (std::find(neighbors.begin(),neighbors.end(),n)!=neighbors.end()) throw std::runtime_error("surface_kmc_3d requires four lateral neighbors");
        m_neighbors[site.id]=neighbors;m_edges[2*site.id]={{site.id,neighbors[0]}};m_edges[2*site.id+1]={{site.id,neighbors[2]}};
        m_incident[site.id]={{2*site.id,2*neighbors[1],2*site.id+1,2*neighbors[3]+1}};
    }
    for (const auto& bottom:m_solid.sites) if (m_down[bottom.id]==n) {
        bool empty=false,solid=false;
        for (auto id=bottom.id;id!=n;id=m_up[id]) {
            const auto& site=m_solid.sites[id];
            if (site.occupied) {if (empty) throw std::runtime_error("initial surface has a buried vacancy or unsupported overhang");solid=true;++m_solidCount;}
            else if (!empty) {if (!solid || site.fixed) throw std::runtime_error("initial surface needs a supported, unfixed binding site");m_binding[id]=true;empty=true;}
        }
        if (!empty) throw std::runtime_error("initial surface needs empty growth headroom");
    }
    for (const auto& rule:m_options.processes) m_catalogs.emplace_back(rule.counting=="undirected_pairs" ? 2*n : n);
    m_counts.assign(m_catalogs.size(),0);m_time=m_options.start;m_initialMean=meanHeight();
    if (!m_options.seed) {std::random_device rd;m_options.seed=(static_cast<std::uint64_t>(rd())<<32)^rd();}
    m_rng.seed(m_options.seed);
    for (std::size_t id=0;id<n;++id) refresh(id);
}
bool Simulation::matches(const Rule& p,std::size_t a,std::size_t b) const {
    return m_binding[a] && m_ads[a]==p.before[0] && (p.before.size()==1 || (m_binding[b] && m_ads[b]==p.before[1]));
}
bool Simulation::eligible(const Rule& p,std::size_t id) const {
    if (p.counting=="sites") return matches(p,id,id);
    if (p.counting=="undirected_pairs") {const auto& e=m_edges[id];return matches(p,e[0],e[1]) || matches(p,e[1],e[0]);}
    for (auto j:m_neighbors[id]) if (matches(p,id,j) || matches(p,j,id)) return true;
    return false;
}
void Simulation::refresh(std::size_t id) {
    for (std::size_t i=0;i<m_catalogs.size();++i) {
        const auto& p=m_options.processes[i];auto& c=m_catalogs[i];
        if (p.counting=="sites") c.set(id,eligible(p,id));
        else if (p.counting=="undirected_pairs") {for (auto e:m_incident[id]) c.set(e,eligible(p,e));}
        else {c.set(id,eligible(p,id));for (auto j:m_neighbors[id]) c.set(j,eligible(p,j));}
    }
}
void Simulation::apply(std::size_t process,std::size_t first,std::size_t second) {
    const auto& p=m_options.processes.at(process);m_ads.at(first);
    if (p.before.size()==2) {
        m_ads.at(second);
        if (std::find(m_neighbors[first].begin(),m_neighbors[first].end(),second)==m_neighbors[first].end()) throw std::runtime_error("process requires lateral same-height first neighbors");
    }
    if (!matches(p,first,second)) throw std::runtime_error("process reactant pattern does not match binding sites");
    // Validate the whole event before mutating any site.
    for (std::size_t slot=0;slot<p.after.size();++slot) if (p.after[slot].solid) {
        const auto id=slot ? second : first;const auto& site=m_solid.sites[id];
        if (site.fixed || site.occupied || m_down[id]==m_solid.sites.size() || !m_solid.sites[m_down[id]].occupied) throw std::runtime_error("solid incorporation requires first-neighbor support");
    }
    for (std::size_t slot=0;slot<p.after.size();++slot) {
        const auto id=slot ? second : first;const auto& product=p.after[slot];
        m_ads[id]=product.solid ? "" : product.species;
        if (product.solid) {
            auto& site=m_solid.sites[id];site.occupied=true;site.species=product.species;m_deposited[id]=true;++m_solidCount;m_binding[id]=false;
            if (m_up[id]==m_solid.sites.size()) m_capacity=true;
            else {m_binding[m_up[id]]=true;refresh(m_up[id]);}
        }
    }
    refresh(first);if (p.before.size()==2) refresh(second);
    ++m_counts[process];++m_steps;
}
std::size_t Simulation::eligibleCount(std::size_t process) const {return m_catalogs.at(process).ids.size();}
bool Simulation::catalogConsistent() const {
    for (std::size_t i=0;i<m_catalogs.size();++i) {
        const auto& c=m_catalogs[i];
        for (std::size_t id=0;id<c.positions.size();++id) if ((c.positions[id]>=0)!=eligible(m_options.processes[i],id)) return false;
        for (std::size_t j=0;j<c.ids.size();++j) if (c.positions[c.ids[j]]!=static_cast<long long>(j)) return false;
        if (std::count_if(c.positions.begin(),c.positions.end(),[](long long pos){return pos>=0;})!=static_cast<long long>(c.ids.size())) return false;
    }
    return true;
}
double Simulation::meanHeight() const {return m_solidCount/(static_cast<double>(m_options.geometry.repetitions[0])*m_options.geometry.repetitions[1]);}
double Simulation::timeAverageHeight() const {return m_time>m_options.start ? static_cast<double>(m_integral/(m_time-m_options.start)) : m_initialMean;}
double Simulation::growthRate() const {return m_time>m_options.start ? (meanHeight()-m_initialMean)/(m_time-m_options.start) : 0;}
void Simulation::run(const std::function<void(const Simulation&,bool)>& observe) {
    if (observe) observe(*this,false);
    while (m_time<m_options.end) {
        if (m_capacity) {m_stop="capacity_limit";break;}
        std::vector<double> propensities;double total=0;
        for (std::size_t i=0;i<m_catalogs.size();++i) {double r=m_options.processes[i].rate*m_catalogs[i].ids.size();propensities.push_back(r);total+=r;}
        if (!std::isfinite(total)) throw std::runtime_error("total process rate overflow");
        if (!total) {m_stop="absorbing";break;}
        if (m_options.maxEvents && m_steps>=m_options.maxEvents) {m_stop="event_limit";break;}
        double u=0;do {u=std::generate_canonical<double,53>(m_rng);} while (!u);
        const double dt=-std::log(u)/total;
        if (dt>m_options.end-m_time) {m_integral+=meanHeight()*(m_options.end-m_time);m_time=m_options.end;break;}
        if (!(dt>0) || m_time+dt<=m_time) throw std::runtime_error("time step below floating-point resolution");
        double pick=std::generate_canonical<double,53>(m_rng)*total;
        std::size_t process=0;while (process+1<propensities.size() && pick>=propensities[process]) {pick-=propensities[process];++process;}
        // Guard rounding at a channel boundary against selecting a zero-rate tail.
        if (!propensities[process]) {for (process=propensities.size();process-->0;) if (propensities[process]) break;}
        const auto& p=m_options.processes[process];const auto& c=m_catalogs[process];
        auto first=c.ids[std::uniform_int_distribution<std::size_t>(0,c.ids.size()-1)(m_rng)],second=first;
        if (p.counting=="undirected_pairs") {const auto e=m_edges[first];first=e[0];second=e[1];}
        else if (p.counting=="reactive_sites") {
            std::array<std::size_t,4> partners{};std::size_t count=0;
            for (auto j:m_neighbors[first]) if (matches(p,first,j) || matches(p,j,first)) partners[count++]=j;
            second=partners[std::uniform_int_distribution<std::size_t>(0,count-1)(m_rng)];
        }
        if (p.before.size()==2) {
            const bool forward=matches(p,first,second),reverse=matches(p,second,first);
            const bool symmetric=p.before[0]==p.before[1] && p.after[0].species==p.after[1].species && p.after[0].solid==p.after[1].solid;
            if (!forward || (reverse && !symmetric && std::uniform_int_distribution<int>(0,1)(m_rng))) std::swap(first,second);
        }
        m_integral+=meanHeight()*dt;m_time+=dt;apply(process,first,second);
        if (observe) observe(*this,false);
    }
    if (m_stop.empty()) m_stop=m_capacity ? "capacity_limit" : "end_time";
    if (observe) observe(*this,true);
}
void Simulation::writeXYZ(std::ostream& out) const {
    std::size_t particles=m_solidCount;for (const auto& a:m_ads) particles+=!a.empty();
    out<<std::setprecision(17)<<particles<<"\nLattice=\"";
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) out<<(i||j ? " " : "")<<m_solid.cell[i][j];
    out<<"\" Properties=species:S:1:pos:R:3:site_id:I:1:fixed:L:1:chemical_species:S:1:role:S:1 pbc=\"T T F\" Time="<<m_time<<'\n';
    for (const auto& site:m_solid.sites) {
        if (!site.occupied && m_ads[site.id].empty()) continue;
        const auto& species=site.occupied ? site.species : m_ads[site.id];
        out<<m_options.elements.at(species)<<' '<<site.position[0]<<' '<<site.position[1]<<' '<<site.position[2]<<' '<<site.id<<' '<<(site.fixed ? "T" : "F")<<' '<<species<<' '<<(site.occupied ? (m_deposited[site.id] ? "film" : "substrate") : "adsorbate")<<'\n';
    }
    if (!out) throw std::runtime_error("failed to write surface XYZ output");
}
void runToFiles(const Options& options) {
    Simulation simulation(options);std::ofstream log("surface.log"),trajectory;
    if (!log) throw std::runtime_error("cannot open surface.log");
    if (options.writeXYZ) {trajectory.open("surface_trajectory.xyz");std::ofstream initial("surface_initial.xyz");if (!initial || !trajectory) throw std::runtime_error("cannot open surface XYZ outputs");simulation.writeXYZ(initial);}
    log<<std::setprecision(17)<<"time\tmean_height\ttime_average_height\tmean_growth_rate";
    for (const auto& rule:options.processes) log<<'\t'<<rule.id;
    log<<'\n';double lastLog=-std::numeric_limits<double>::infinity(),lastFrame=lastLog;
    simulation.run([&](const Simulation& s,bool final) {
        if ((final && s.time()!=lastLog) || s.time()-lastLog>=options.logInterval) {
            log<<s.time()<<'\t'<<s.meanHeight()<<'\t'<<s.timeAverageHeight()<<'\t'<<s.growthRate();for (auto count:s.counts()) log<<'\t'<<count;log<<'\n';lastLog=s.time();
        }
        if (options.writeXYZ && ((final && s.time()!=lastFrame) || !std::isfinite(lastFrame) || (options.snapshotInterval>0 && s.time()-lastFrame>=options.snapshotInterval))) {s.writeXYZ(trajectory);lastFrame=s.time();}
    });
    nlohmann::json counts=nlohmann::json::object(),ads=nlohmann::json::object(),film=nlohmann::json::object(),rates=nlohmann::json::object();
    for (std::size_t i=0;i<options.processes.size();++i) {counts[options.processes[i].id]=simulation.counts()[i];rates[options.processes[i].id]=options.processes[i].rate;}
    for (const auto& name:simulation.adsorbates()) if (!name.empty()) ads[name]=ads.value(name,std::uint64_t(0))+1;
    for (const auto& site:simulation.solid().sites) if (site.occupied) film[site.species]=film.value(site.species,std::uint64_t(0))+1;
    nlohmann::json summary={{"model","surface_kmc_3d"},{"representation","explicit_3d_occupancy"},{"repetitions",options.geometry.repetitions},{"start_time",options.start},{"requested_end_time",options.end},{"time",simulation.time()},{"seed",simulation.options().seed},{"stop_reason",simulation.stopReason()},{"capacity_limited",simulation.capacityLimited()},{"poisoned",simulation.stopReason()=="absorbing"},{"events",simulation.steps()},{"solid_units",simulation.solid().atomCount()},{"mean_height",simulation.meanHeight()},{"time_average_height",simulation.timeAverageHeight()},{"mean_growth_rate",simulation.growthRate()},{"process_counts",counts},{"process_rates",rates},{"adsorbate_counts",ads},{"solid_species_counts",film}};
    std::ofstream result("surface_summary.json");result<<summary.dump(2)<<'\n';result.close();log.close();
    if (options.writeXYZ) {std::ofstream final("surface_final.xyz");simulation.writeXYZ(final);final.close();trajectory.close();if (!final || !trajectory) throw std::runtime_error("failed to write surface geometry");}
    if (!result || !log) throw std::runtime_error("failed to write surface summary/log");
}
}
