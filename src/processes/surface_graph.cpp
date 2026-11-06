#include "surface_kmc.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <stdexcept>

namespace SurfaceKMC {
namespace {
Atomistic::Vec3 cross(const Atomistic::Vec3& a,const Atomistic::Vec3& b) {
    return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
double dot(const Atomistic::Vec3& a,const Atomistic::Vec3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
}
void Simulation::initializeGraph() {
    const auto n=m_solid.sites.size();
    if (!n) throw std::runtime_error("geometry must contain candidate sites");
    const auto normal=cross(m_solid.cell[0],m_solid.cell[1]);m_area=std::sqrt(dot(normal,normal));
    const double volume=dot(normal,m_solid.cell[2]);
    if (!std::isfinite(volume) || volume<=0 || !std::isfinite(m_area) || m_area<=0 || std::abs(normal[0])+std::abs(normal[1])>1e-8*m_area)
        throw std::runtime_error("surface cell must have a positive volume and its first two vectors in the xy plane");
    m_siteVolume=m_options.siteVolume;
    if (!m_options.customGeometry && m_siteVolume==0) m_siteVolume=(volume-m_options.geometry.vacuum*m_area)/n;
    if (!std::isfinite(m_siteVolume) || m_siteVolume<=0) throw std::runtime_error("custom geometry requires positive site_volume for thickness normalization");
    m_top=-std::numeric_limits<double>::infinity();std::vector<bool> occupied(n,false);
    for (const auto& site:m_solid.sites) {
        if (site.id>=n || &site!=&m_solid.sites[site.id]) throw std::runtime_error("site IDs must be contiguous and stable");
        if (!m_options.elements.count(site.species)) throw std::runtime_error("geometry contains an undefined species");
        if (site.fixed && !site.occupied) throw std::runtime_error("fixed substrate sites must be occupied");
        occupied[site.id]=site.occupied;m_solidCount+=site.occupied;m_top=std::max(m_top,site.position[2]);
    }
    if (!m_solidCount || !graphSupported(occupied)) throw std::runtime_error("initial solids require first-neighbor bonds and connectivity to the fixed substrate when present");
    m_relations.resize(m_options.processes.size());
    std::map<std::tuple<std::string,double,double,bool>,std::shared_ptr<const Graph>> cached;
    for (std::size_t i=0;i<m_options.processes.size();++i) for (const auto& link:m_options.processes[i].links) {
        const bool lateral=link.sameHeight || link.neighbors=="lateral_first";
        const auto key=std::make_tuple(link.neighbors,link.minDistance,link.maxDistance,lateral);
        auto found=cached.find(key);
        if (found==cached.end()) {
            auto graph=std::make_shared<Graph>(n);
            if (link.neighbors=="distance") *graph=Atomistic::neighborGraph(m_solid,link.minDistance,link.maxDistance,lateral);
            else for (const auto& site:m_solid.sites) for (auto j:site.neighbors) {
                if (!lateral || std::abs(site.position[2]-m_solid.sites[j].position[2])<=m_solid.neighborCutoff*1e-7) (*graph)[site.id].push_back(j);
            }
            found=cached.emplace(key,graph).first;
        }
        m_relations[i].push_back(found->second);
    }
    m_counts.assign(m_options.processes.size(),0);m_time=m_options.start;m_initialMean=meanHeight();
    if (!m_options.seed) {std::random_device rd;m_options.seed=(static_cast<std::uint64_t>(rd())<<32)^rd();}
    m_rng.seed(m_options.seed);rebuildGraph();
}
void Simulation::graphAccess(const std::vector<bool>& occupied,std::vector<bool>& binding,std::vector<bool>& exposed) const {
    const auto n=occupied.size();std::vector<bool> gas(n,false);std::vector<std::size_t> pending;
    for (const auto& site:m_solid.sites) if (!occupied[site.id] && site.position[2]>=m_top-m_solid.neighborCutoff*1e-7) {gas[site.id]=true;pending.push_back(site.id);}
    for (std::size_t k=0;k<pending.size();++k) for (auto j:m_solid.sites[pending[k]].neighbors) if (!occupied[j] && !gas[j]) {gas[j]=true;pending.push_back(j);}
    binding.assign(n,false);exposed.assign(n,false);
    for (const auto& site:m_solid.sites) {
        bool support=false,touchesGas=site.position[2]>=m_top-m_solid.neighborCutoff*1e-7;
        for (auto j:site.neighbors) {support=support || occupied[j];touchesGas=touchesGas || gas[j];}
        binding[site.id]=!site.fixed && !occupied[site.id] && gas[site.id] && support;
        exposed[site.id]=!site.fixed && occupied[site.id] && touchesGas;
    }
}
bool Simulation::graphSupported(const std::vector<bool>& occupied) const {
    std::vector<bool> anchored(occupied.size(),false);std::vector<std::size_t> pending;
    for (const auto& site:m_solid.sites) if (occupied[site.id]) {
        bool bonded=false;for (auto j:site.neighbors) bonded=bonded || occupied[j];
        if (!bonded) return false;
        if (site.fixed) {anchored[site.id]=true;pending.push_back(site.id);}
    }
    if (pending.empty()) return true;
    for (std::size_t k=0;k<pending.size();++k) for (auto j:m_solid.sites[pending[k]].neighbors) if (occupied[j] && !anchored[j]) {anchored[j]=true;pending.push_back(j);}
    for (const auto& site:m_solid.sites) if (occupied[site.id] && !anchored[site.id]) return false;
    return true;
}
bool Simulation::graphSlot(const Rule& p,std::size_t slot,std::size_t id) const {
    if (!p.beforeSolid.empty() && p.beforeSolid[slot]) return m_exposed[id] && m_solid.sites[id].species==p.before[slot];
    return m_binding[id] && m_ads[id]==p.before[slot];
}
bool Simulation::graphTransition(const Rule& p,const std::vector<std::size_t>& sites) const {
    bool occupancyChanged=false;
    for (std::size_t slot=0;slot<sites.size();++slot) {
        const auto& site=m_solid.sites[sites[slot]];
        if (site.fixed) return false;
        occupancyChanged|=site.occupied!=p.after[slot].solid;
    }
    if (!occupancyChanged) return true;
    // Validate the complete simultaneous update, including substrate anchoring
    // and adsorbates that could lose support or become sealed inside a cavity.
    std::vector<bool> occupied(m_solid.sites.size());
    for (const auto& site:m_solid.sites) occupied[site.id]=site.occupied;
    for (std::size_t slot=0;slot<sites.size();++slot) occupied[sites[slot]]=p.after[slot].solid;
    if (!graphSupported(occupied)) return false;
    std::vector<bool> binding,exposed;graphAccess(occupied,binding,exposed);
    for (std::size_t id=0;id<m_ads.size();++id) {
        bool adsorbed=!m_ads[id].empty();
        for (std::size_t slot=0;slot<sites.size();++slot) if (sites[slot]==id) adsorbed=!p.after[slot].solid && !p.after[slot].species.empty();
        if (adsorbed && !binding[id]) return false;
    }
    return true;
}
std::vector<std::vector<Simulation::EventGroup>> Simulation::enumerateGraph() const {
    const auto n=m_solid.sites.size();std::vector<std::vector<EventGroup>> catalogs(m_options.processes.size());
    std::size_t steps=0,embeddings=0;
    for (std::size_t process=0;process<m_options.processes.size();++process) {
        const auto& p=m_options.processes[process];const auto slots=p.before.size();
        std::vector<std::size_t> assignment(slots,n);std::vector<bool> used(n,false);
        std::map<std::vector<std::size_t>,EventGroup> groups;
        std::function<void(std::size_t)> search=[&](std::size_t depth) {
            if (depth==slots) {
                if (!graphTransition(p,assignment)) return;
                if (++embeddings>m_options.maxEmbeddings) throw std::runtime_error("matching exceeds model.parameters.max_embeddings; increase the explicit limit or reduce the pattern/system");
                auto add=[&](std::vector<std::size_t> key) {auto& a=groups[std::move(key)].assignments;if(std::find(a.begin(),a.end(),assignment)==a.end()) a.push_back(assignment);};
                if (p.counting=="reactive_sites") {add({assignment[0]});add({assignment[1]});}
                else if (p.counting=="undirected_pairs" || p.counting=="unique_site_sets") {auto key=assignment;std::sort(key.begin(),key.end());add(std::move(key));}
                else add(assignment);
                return;
            }
            std::size_t slot=slots;const std::vector<std::size_t>* candidates=nullptr;
            if (!depth) slot=0;
            else for (std::size_t i=0;i<p.links.size();++i) {
                const auto& link=p.links[i];std::size_t pending=slots,assigned=n;
                if (assignment[link.first]==n && assignment[link.second]!=n) {pending=link.first;assigned=assignment[link.second];}
                if (assignment[link.second]==n && assignment[link.first]!=n) {pending=link.second;assigned=assignment[link.first];}
                if (pending!=slots && (!candidates || (*m_relations[process][i])[assigned].size()<candidates->size())) {slot=pending;candidates=&(*m_relations[process][i])[assigned];}
            }
            if (slot==slots) throw std::logic_error("disconnected matching pattern");
            auto trySite=[&](std::size_t id) {
                if (++steps>m_options.maxSearchSteps) throw std::runtime_error("matching exceeds model.parameters.max_search_steps; increase the explicit limit or reduce the pattern/system");
                if (used[id] || !graphSlot(p,slot,id)) return;
                for (std::size_t i=0;i<p.links.size();++i) {
                    const auto& link=p.links[i];std::size_t other=slots;
                    if(link.first==slot) other=link.second;else if(link.second==slot) other=link.first;
                    if(other!=slots && assignment[other]!=n) {const auto& neighbors=(*m_relations[process][i])[id];if (!std::binary_search(neighbors.begin(),neighbors.end(),assignment[other])) return;}
                }
                assignment[slot]=id;used[id]=true;search(depth+1);used[id]=false;assignment[slot]=n;
            };
            if (candidates) {for(auto id:*candidates) trySite(id);} else for(std::size_t id=0;id<n;++id) trySite(id);
        };
        search(0);for(auto& entry:groups) catalogs[process].push_back(std::move(entry.second));
    }
    return catalogs;
}
void Simulation::rebuildGraph() {
    std::vector<bool> occupied(m_solid.sites.size());for(const auto& site:m_solid.sites) occupied[site.id]=site.occupied;
    graphAccess(occupied,m_binding,m_exposed);m_graphCatalogs=enumerateGraph();
}
void Simulation::apply(std::size_t process,const std::vector<std::size_t>& sites) {
    if (!graphMode()) {if(sites.size()==1) apply(process,sites[0]);else if(sites.size()==2) apply(process,sites[0],sites[1]);else throw std::runtime_error("multi-site events require gas_exposed mode");return;}
    const auto& p=m_options.processes.at(process);
    if(sites.size()!=p.before.size()) throw std::runtime_error("event site vector must match the rule slot count");
    std::set<std::size_t> unique;
    for(std::size_t slot=0;slot<sites.size();++slot) {
        if(sites[slot]>=m_solid.sites.size() || !unique.insert(sites[slot]).second || !graphSlot(p,slot,sites[slot])) throw std::runtime_error("event does not match distinct eligible reactant sites");
    }
    for(std::size_t i=0;i<p.links.size();++i) {const auto& link=p.links[i];const auto& neighbors=(*m_relations[process][i])[sites[link.first]];if(!std::binary_search(neighbors.begin(),neighbors.end(),sites[link.second])) throw std::runtime_error("event violates a geometric link");}
    if(!graphTransition(p,sites)) throw std::runtime_error("event would break solid support, detach the film, or bury an adsorbate");
    for(std::size_t slot=0;slot<sites.size();++slot) {
        auto& site=m_solid.sites[sites[slot]];const auto& product=p.after[slot];
        if(!site.occupied && product.solid) {++m_solidCount;m_deposited[site.id]=true;if(site.position[2]>=m_top-m_solid.neighborCutoff*1e-7) m_capacity=true;}
        else if(site.occupied && !product.solid) {--m_solidCount;m_deposited[site.id]=false;}
        site.occupied=product.solid;if(product.solid) site.species=product.species;
        m_ads[site.id]=product.solid ? "" : product.species;
    }
    ++m_counts[process];++m_steps;rebuildGraph();
}
void Simulation::runGraph(const std::function<void(const Simulation&,bool)>& observe) {
    if(observe) observe(*this,false);
    while(m_time<m_options.end) {
        if(m_capacity) {m_stop="capacity_limit";break;}
        std::vector<double> rates;double total=0;
        for(std::size_t i=0;i<m_graphCatalogs.size();++i) {rates.push_back(m_options.processes[i].rate*m_graphCatalogs[i].size());total+=rates.back();}
        if(!std::isfinite(total)) throw std::runtime_error("total process rate overflow");
        if(!total) {m_stop="absorbing";break;}
        if(m_options.maxEvents && m_steps>=m_options.maxEvents) {m_stop="event_limit";break;}
        double u=0;do {u=std::generate_canonical<double,53>(m_rng);} while(!u);
        const double dt=-std::log(u)/total;
        if(dt>m_options.end-m_time) {m_integral+=meanHeight()*(m_options.end-m_time);m_time=m_options.end;break;}
        if(!(dt>0) || m_time+dt<=m_time) throw std::runtime_error("time step below floating-point resolution");
        double pick=std::generate_canonical<double,53>(m_rng)*total;std::size_t process=0;
        while(process+1<rates.size() && pick>=rates[process]) {pick-=rates[process];++process;}
        if(!rates[process]) for(process=rates.size();process-->0;) if(rates[process]) break;
        const auto& catalog=m_graphCatalogs[process];const auto& group=catalog[std::uniform_int_distribution<std::size_t>(0,catalog.size()-1)(m_rng)];
        auto sites=group.assignments[std::uniform_int_distribution<std::size_t>(0,group.assignments.size()-1)(m_rng)];
        m_integral+=meanHeight()*dt;m_time+=dt;apply(process,sites);
        if(observe) observe(*this,false);
    }
    if(m_stop.empty()) m_stop=m_capacity ? "capacity_limit" : "end_time";
    if(observe) observe(*this,true);
}
}
