#include "atomistic_kmc.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>

namespace Atomistic {
namespace {
// The fixed substrate is the anchor. Test a proposed removal/hop without
// mutating the surface; all occupied components must reach a fixed atom.
bool connectedToSubstrate(const Surface& surface,std::size_t source,std::size_t destination) {
    std::vector<bool> seen(surface.sites.size(),false);
    std::vector<std::size_t> queue;
    auto occupied=[&](std::size_t id) {
        return id!=source && (id==destination || surface.sites[id].occupied);
    };
    for (const auto& site:surface.sites) if (site.fixed && occupied(site.id)) {
        seen[site.id]=true;queue.push_back(site.id);
    }
    if (queue.empty()) return true; // No fixed substrate was requested.
    for (std::size_t i=0;i<queue.size();++i)
        for (auto j:surface.sites[queue[i]].neighbors) if (!seen[j] && occupied(j)) {
            seen[j]=true;queue.push_back(j);
        }
    for (const auto& site:surface.sites) if (occupied(site.id) && !seen[site.id]) return false;
    return true;
}
// Articulation points in the occupied graph plus a virtual root connected to
// every fixed atom. Iterative DFS avoids recursion on large candidate graphs.
std::vector<bool> substrateBridges(const Surface& surface) {
    const auto n=surface.sites.size();
    std::vector<std::size_t> roots;
    for (const auto& s:surface.sites) if (s.fixed && s.occupied) roots.push_back(s.id);
    std::vector<bool> bridge(n+1,false);
    if (roots.empty()) return bridge;
    std::vector<long long> discovery(n+1,-1),low(n+1),parent(n+1,-1);
    std::vector<std::size_t> next(n+1,0),stack{n};
    long long clock=0;discovery[n]=low[n]=clock++;
    while (!stack.empty()) {
        auto u=stack.back();
        const auto degree=u==n ? roots.size() : surface.sites[u].neighbors.size()+(surface.sites[u].fixed ? 1 : 0);
        if (next[u]<degree) {
            const auto index=next[u]++;
            const auto v=u==n ? roots[index] : (index<surface.sites[u].neighbors.size() ? surface.sites[u].neighbors[index] : n);
            if ((v!=n && !surface.sites[v].occupied) || static_cast<long long>(v)==parent[u]) continue;
            if (discovery[v]<0) {
                parent[v]=static_cast<long long>(u);discovery[v]=low[v]=clock++;stack.push_back(v);
            } else low[u]=std::min(low[u],discovery[v]);
        } else {
            stack.pop_back();
            if (parent[u]>=0) {
                auto p=static_cast<std::size_t>(parent[u]);low[p]=std::min(low[p],low[u]);
                if (p!=n && low[u]>=discovery[p]) bridge[p]=true;
            }
        }
    }
    return bridge;
}
} // namespace
void validate(const RunOptions& o) {
    if (!std::isfinite(o.start) || !std::isfinite(o.end) || o.start<0 || o.end<=o.start ||
        !std::isfinite(o.temperature) || o.temperature<=0 ||
        !std::isfinite(o.logInterval) || o.logInterval<=0 || !std::isfinite(o.snapshotInterval) || o.snapshotInterval<=0)
        throw std::runtime_error("atomistic simulation: positive temperature, duration and output intervals are required");
    if (o.processes.empty()) throw std::runtime_error("processes: at least one atomistic process is required");
    std::set<std::string> ids;
    for (const auto& p : o.processes) {
        if (p.id.empty() || !ids.insert(p.id).second || p.species.empty() || p.element.empty())
            throw std::runtime_error("processes: unique IDs and an atomic species with element are required");
        if (p.model!="constant" && p.model!="arrhenius") throw std::runtime_error("processes: unsupported atomistic rate model");
        for (double v : {p.value,p.prefactor,p.activationEnergy,p.bondEnergy})
            if (!std::isfinite(v) || v<0) throw std::runtime_error("processes: rate parameters must be finite and nonnegative");
        if (p.kind!=Kind::Desorption && p.minCoordination<1)
            throw std::runtime_error("processes: adsorption and diffusion require min_coordination >= 1 occupied first-nearest neighbor");
        if (p.minCoordination<0 || p.maxCoordination<p.minCoordination) throw std::runtime_error("processes: invalid coordination bounds");
    }
}
Simulation::Simulation(Surface surface, RunOptions options)
    : m_surface(std::move(surface)), m_options(std::move(options)), m_time(m_options.start) {
    validate(m_options);
    if (m_surface.periodic[2]) throw std::runtime_error("atomistic simulation: gas access requires nonperiodic z");
    if (m_surface.sites.empty()) throw std::runtime_error("atomistic simulation: no candidate sites");
    for (const auto& site : m_surface.sites)
        if (site.occupied && m_surface.occupiedCoordination(site.id)==0)
            throw std::runtime_error("atomistic initial state: occupied site " + std::to_string(site.id) +
                                     " has no occupied first-nearest neighbor bond");
    if (!connectedToSubstrate(m_surface,m_surface.sites.size(),m_surface.sites.size()))
        throw std::runtime_error("atomistic initial state: occupied atoms must be connected to the fixed substrate");
    if (m_options.seed==0) { std::random_device device; m_options.seed=(static_cast<std::uint64_t>(device())<<32)^device(); }
    m_rng.seed(m_options.seed);
    m_counts.resize(m_options.processes.size(),0);
    for (const auto& p : m_options.processes) {
        auto existing=m_surface.elements.find(p.species);
        if (existing!=m_surface.elements.end() && existing->second!=p.element) throw std::runtime_error("processes: conflicting species elements");
        m_surface.elements[p.species]=p.element;
    }
}
std::vector<bool> Simulation::accessible() const {
    std::vector<bool> gas(m_surface.sites.size(),false);
    double top=-std::numeric_limits<double>::infinity();
    for (const auto& s : m_surface.sites) top=std::max(top,s.position[2]);
    const double tolerance=1e-9*std::max(1.0,m_surface.neighborCutoff);
    std::queue<std::size_t> pending;
    for (const auto& s : m_surface.sites)
        if (!s.occupied && s.position[2]>=top-tolerance) {gas[s.id]=true;pending.push(s.id);}
    while (!pending.empty()) {
        auto id=pending.front();pending.pop();
        for (auto j : m_surface.sites[id].neighbors)
            if (!m_surface.sites[j].occupied && !gas[j]) {gas[j]=true;pending.push(j);}
    }
    return gas;
}
double Simulation::rate(const Process& p,std::size_t coordination) const {
    if (p.model=="constant") return p.value;
    // Molar energies in J/mol; no implicit multiplication by coordination.
    const double barrier=p.activationEnergy+p.bondEnergy*static_cast<double>(coordination);
    return p.prefactor*std::exp(-barrier/(8.31446261815324*m_options.temperature));
}
std::vector<Event> Simulation::events() const {
    auto gas=accessible();
    double top=-std::numeric_limits<double>::infinity();
    for (const auto& s : m_surface.sites) top=std::max(top,s.position[2]);
    auto exposed=[&](const Site& s) {
        if (s.position[2]>=top-1e-9*std::max(1.0,m_surface.neighborCutoff)) return true;
        for (auto j : s.neighbors) if (gas[j]) return true;
        return false;
    };
    std::vector<std::size_t> coordination(m_surface.sites.size());
    for (const auto& s : m_surface.sites) coordination[s.id]=m_surface.occupiedCoordination(s.id);
    const auto bridges=substrateBridges(m_surface);
    // Removing a source may strand a different atom. Only its occupied
    // neighbors can lose a bond; a simultaneous hop can replace that bond.
    auto preservesBonds=[&](std::size_t source,std::size_t destination) {
        for (auto j : m_surface.sites[source].neighbors) {
            const auto& neighbor=m_surface.sites[j];
            if (!neighbor.occupied || coordination[j]>1) continue;
            const bool replacement=destination!=source &&
                std::find(neighbor.neighbors.begin(),neighbor.neighbors.end(),destination)!=neighbor.neighbors.end();
            if (!replacement) return false;
        }
        if (!bridges[source]) return true;
        if (destination==source) return false;
        return connectedToSubstrate(m_surface,source,destination);
    };
    std::vector<Event> result;
    auto add=[&](std::size_t p,std::size_t source,std::size_t target,std::size_t coordination) {
        double k=rate(m_options.processes[p],coordination);
        if (!std::isfinite(k) || k<0) throw std::runtime_error("atomistic event rate is not finite");
        if (k>0) result.push_back({p,source,target,k});
    };
    for (std::size_t p=0;p<m_options.processes.size();++p) {
        const auto& process=m_options.processes[p];
        auto allowed=[&](std::size_t coordination) {return coordination>=static_cast<std::size_t>(process.minCoordination) && coordination<=static_cast<std::size_t>(process.maxCoordination);};
        for (const auto& s : m_surface.sites) {
            if (s.fixed) continue;
            if (process.kind==Kind::Adsorption) {
                if (!s.occupied && gas[s.id]) {
                    auto n=coordination[s.id];
                    if (allowed(n)) add(p,s.id,s.id,n);
                }
            } else if (s.occupied && s.species==process.species && exposed(s)) {
                const auto n=coordination[s.id];
                if (process.kind==Kind::Desorption) {
                    if (allowed(n) && preservesBonds(s.id,s.id)) add(p,s.id,s.id,n);
                } else {
                    for (auto j : s.neighbors) {
                        const auto& target=m_surface.sites[j];
                        if (target.fixed || target.occupied || !gas[j]) continue;
                        // The departing source does not support the destination.
                        const auto targetCoordination=coordination[j]-1;
                        if (allowed(targetCoordination) && preservesBonds(s.id,j)) add(p,s.id,j,n);
                    }
                }
            }
        }
    }
    return result;
}
void Simulation::execute(const Event& event) {
    const auto& p=m_options.processes[event.process];
    auto& source=m_surface.sites[event.source];
    auto& target=m_surface.sites[event.target];
    if (p.kind==Kind::Adsorption) { target.occupied=true;target.species=p.species; }
    else if (p.kind==Kind::Desorption) {source.occupied=false;source.species.clear();}
    else {target.occupied=true;target.species=source.species;source.occupied=false;source.species.clear();}
    ++m_counts[event.process];++m_steps;
}
void Simulation::apply(const Event& event) {
    const auto candidates=events();
    auto found=std::find_if(candidates.begin(),candidates.end(),[&](const Event& e){return e.process==event.process && e.source==event.source && e.target==event.target;});
    if (found==candidates.end()) throw std::runtime_error("atomistic event is stale or ineligible");
    execute(*found);
}
void Simulation::run(const std::function<void(const Simulation&,bool)>& observe) {
    if (observe) observe(*this,false);
    while (m_time<m_options.end) {
        if (m_options.maxEvents && m_steps>=m_options.maxEvents) {m_stopReason="event_limit";break;}
        const auto candidates=events();
        double total=0;
        for (const auto& e : candidates) total+=e.rate;
        if (!std::isfinite(total)) throw std::runtime_error("atomistic total rate overflow");
        if (total==0) {m_time=m_options.end;m_stopReason="no_events";break;}
        double u=0;
        do {u=std::generate_canonical<double,53>(m_rng);} while (u==0);
        // Waiting time uses the pre-event total escape rate.
        const double dt=-std::log(u)/total;
        if (dt>m_options.end-m_time) {m_time=m_options.end;break;}
        if (!(dt>0) || m_time+dt<=m_time) throw std::runtime_error("atomistic time step is below floating-point resolution");
        const double threshold=std::generate_canonical<double,53>(m_rng)*total;
        double cumulative=0;
        const Event* selected=&candidates.back();
        for (const auto& e : candidates) {cumulative+=e.rate;if (threshold<cumulative) {selected=&e;break;}}
        m_time+=dt;
        execute(*selected);
        if (observe) observe(*this,false);
    }
    if (m_stopReason.empty()) m_stopReason="end_time";
    if (observe) observe(*this,true);
}
void runToFiles(Surface surface,RunOptions options) {
    Simulation simulation(std::move(surface),std::move(options));
    std::ofstream log("atomistic.log"),trajectory("atomistic_trajectory.xyz");
    if (!log || !trajectory) throw std::runtime_error("cannot open atomistic log or trajectory");
    log<<std::setprecision(17)<<"# seed "<<simulation.options().seed<<"\ntime\tatoms\tevents";
    for (const auto& p : simulation.options().processes) log<<'\t'<<p.id;
    log<<'\n';
    simulation.surface().writeXYZ("atomistic_initial.xyz",simulation.time());
    double lastLog=-std::numeric_limits<double>::infinity(),lastFrame=lastLog;
    simulation.run([&](const Simulation& s,bool final) {
        if ((final && s.time()!=lastLog) || s.time()-lastLog>=s.options().logInterval) {
            log<<s.time()<<'\t'<<s.surface().atomCount()<<'\t'<<s.steps();
            for (auto count:s.counts()) log<<'\t'<<count;
            log<<'\n';lastLog=s.time();
        }
        if ((final && s.time()!=lastFrame) || s.time()-lastFrame>=s.options().snapshotInterval) {
            s.surface().writeXYZFrame(trajectory,s.time());lastFrame=s.time();
        }
    });
    log<<"# stop "<<simulation.stopReason()<<'\n';
    simulation.surface().writeXYZ("atomistic_final.xyz",simulation.time());
    log.close();trajectory.close();
    if (!log || !trajectory) throw std::runtime_error("failed to write atomistic output");
}
} // namespace Atomistic
