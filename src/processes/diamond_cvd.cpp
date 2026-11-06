#include "diamond_cvd.h"
#include "../extLibs/nlohmann/json.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>
#include <numeric>

namespace DiamondCVD {
namespace {
bool has(const std::vector<std::size_t>& v,std::size_t x) {return std::find(v.begin(),v.end(),x)!=v.end();}
}
double rateCoefficient(double A,double n,double E,double H,double S,double correction,double T,bool reverse) {
    if(!std::isfinite(A) || A<0 || !std::isfinite(n) || !std::isfinite(E) || E<0 || !std::isfinite(H) || !std::isfinite(S) || !std::isfinite(correction) || correction<=0 || !std::isfinite(T) || T<=0) throw std::runtime_error("invalid thermochemical rate parameters");
    if(!A) return 0;
    const double logRate=std::log(A)+n*std::log(T)-E/(.00198720425864083*T)+(reverse ? std::log(correction)-S/1.98720425864083+H/(.00198720425864083*T) : 0);
    const double rate=std::exp(logRate);if(!std::isfinite(rate)) throw std::runtime_error("thermochemical rate overflow");return rate;
}
void validate(const Options& o) {
    if(o.lattice.sites.empty() || o.lattice.periodic[2]) throw std::runtime_error("diamond CVD requires a fixed carbon lattice with open z");
    if(!std::isfinite(o.temperature) || o.temperature<=0 || !std::isfinite(o.gasTemperature) || o.gasTemperature<=0 || !std::isfinite(o.siteVolume) || o.siteVolume<=0 || !std::isfinite(o.start) || o.start<0 || !std::isfinite(o.end) || o.end<=o.start || !std::isfinite(o.logInterval) || o.logInterval<=0 || !std::isfinite(o.snapshotInterval) || o.snapshotInterval<0 || !std::isfinite(o.chLength) || o.chLength<=0) throw std::runtime_error("invalid diamond simulation temperatures, times or geometry scale");
    std::set<std::string> ids;
    const std::set<std::string> kinds={"hydrogen_addition","hydrogen_removal","adsorb_carbon","desorb_hydrocarbon","append_methyl","dimer_opening","dimer_closing","carbon_incorporation"};
    if(o.processes.empty()) throw std::runtime_error("explicit diamond processes are required");
    for(const auto& p:o.processes) {
        if(p.id.empty() || !ids.insert(p.id).second || !kinds.count(p.event) || !std::isfinite(p.rate) || p.rate<0) throw std::runtime_error("invalid or duplicate diamond process");
        if((p.event=="hydrogen_addition" || p.event=="hydrogen_removal") && p.domain!="diamond" && p.domain!="C1" && p.domain!="C2") throw std::runtime_error("hydrogen processes require diamond, C1 or C2 domain");
        if(p.event=="adsorb_carbon" && (p.domain!="C1" && p.domain!="C2")) throw std::runtime_error("adsorption domain must be C1 or C2");
        if(p.event=="adsorb_carbon" && ((p.domain=="C1" && (p.hydrogens<0 || p.hydrogens>3)) || (p.domain=="C2" && p.hydrogens!=2))) throw std::runtime_error("C1 adsorption supports 0..3 hydrogens; acetylene adsorption requires C2H2");
        if(p.event=="desorb_hydrocarbon" && (p.domain!="C1" && p.domain!="C2")) throw std::runtime_error("desorption domain must be C1 or C2");
    }
}
Simulation::Simulation(Options options):m_options(std::move(options)) {
    validate(m_options);const auto& s=m_options.lattice;const auto n=s.sites.size();
    m_atoms.resize(5*n);m_bonds.resize(n);m_frozenVacancy.assign(n,false);m_counts.assign(m_options.processes.size(),0);
    m_top=-std::numeric_limits<double>::infinity();m_initialTop=m_top;
    for(const auto& site:s.sites) {
        if(site.id>=n || site.neighbors.size()>4 || (site.fixed && !site.occupied)) throw std::runtime_error("diamond graph must have stable IDs and at most four first neighbors");
        m_top=std::max(m_top,site.position[2]);
        if(site.occupied) {m_atoms[site.id].kind=Kind::Diamond;++m_diamondCount;m_initialTop=std::max(m_initialTop,site.position[2]);}
    }
    if(!count(Kind::Diamond)) throw std::runtime_error("diamond substrate is empty");
    for(const auto& site:s.sites) if(carbon(site.id)) for(auto j:site.neighbors) if(carbon(j)) m_bonds[site.id].push_back(j);
    auto a=s.cell[0],b=s.cell[1];Atomistic::Vec3 normal{{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};m_area=std::sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.0));
    if(!std::isfinite(m_area) || m_area<=0 || std::abs(normal[0])+std::abs(normal[1])>1e-8*m_area) throw std::runtime_error("diamond surface requires in-plane first two cell vectors");
    gasAccess(false);
    if(m_options.passivate) for(auto i:m_surface) if(m_atoms[i].kind==Kind::Diamond) for(auto j:s.sites[i].neighbors) if(m_atoms[j].kind==Kind::Empty && m_gas[j]) {m_atoms[hydrogenSlot(i,j)]={Kind::Hydrogen,i,false};}
    m_time=m_options.start;m_preciseTime=m_time;m_initialHeight=meanHeight();if(!m_options.seed) {std::random_device rd;m_options.seed=(static_cast<std::uint64_t>(rd())<<32)^rd();}m_rng.seed(m_options.seed);
    rebuild();if(!consistent()) throw std::runtime_error("invalid initial diamond bond network");
}
bool Simulation::carbon(std::size_t id) const {return id<m_options.lattice.sites.size() && (m_atoms[id].kind==Kind::Diamond || m_atoms[id].kind==Kind::Hydrocarbon);}
std::size_t Simulation::hydrogenSlot(std::size_t parent,std::size_t target) const {
    const auto& neighbors=m_options.lattice.sites[parent].neighbors;auto found=std::find(neighbors.begin(),neighbors.end(),target);
    if(found==neighbors.end()) throw std::logic_error("H decoration requires a carbon bond direction");
    return m_options.lattice.sites.size()+4*parent+static_cast<std::size_t>(found-neighbors.begin());
}
std::vector<std::size_t> Simulation::hydrogens(std::size_t id) const {
    std::vector<std::size_t> result;for(auto j:m_options.lattice.sites[id].neighbors) {auto h=hydrogenSlot(id,j);if(m_atoms[h].kind==Kind::Hydrogen) result.push_back(h);}return result;
}
bool Simulation::free(std::size_t id,std::size_t target) const {
    if(has(m_bonds[id],target)) return false;
    if(m_atoms[hydrogenSlot(id,target)].kind==Kind::Hydrogen) return false;
    for(const auto& d:m_dimers) if(d.shared==target && (d.first==id || d.second==id)) return false;
    return true;
}
bool Simulation::exposed(std::size_t id) const {
    if(id>=m_options.lattice.sites.size()) return !m_atoms[id].trapped && exposed(m_atoms[id].parent);
    if(m_atoms[id].trapped || m_frozenVacancy[id]) return false;
    for(auto j:m_options.lattice.sites[id].neighbors) if(m_gas[j]) return true;
    return m_options.lattice.sites[id].position[2]>=m_top-1e-7;
}
void Simulation::gasAccess(bool freeze) {
    const auto& sites=m_options.lattice.sites;std::vector<std::size_t> pending;m_gas.assign(sites.size(),false);
    // H decorations do not form an artificial impermeable carbon-plane barrier.
    // Carbon occupancy defines the void network used to identify sealed defects.
    for(const auto& site:sites) if(!carbon(site.id) && site.position[2]>=m_top-1e-7) {m_gas[site.id]=true;pending.push_back(site.id);}
    for(std::size_t k=0;k<pending.size();++k) for(auto j:sites[pending[k]].neighbors) if(!carbon(j) && !m_gas[j]) {m_gas[j]=true;pending.push_back(j);}
    m_surface.clear();for(const auto& site:sites) if(carbon(site.id) && !site.fixed) {for(auto j:site.neighbors) if(m_gas[j]) {m_surface.push_back(site.id);break;}}
    if(freeze) {
        std::vector<bool> visited(sites.size(),false);
        for(const auto& start:sites) if(!carbon(start.id) && !m_gas[start.id] && !visited[start.id]) {
            std::vector<std::size_t> cavity{start.id};visited[start.id]=true;bool stable=true;
            for(std::size_t k=0;k<cavity.size();++k) for(auto j:sites[cavity[k]].neighbors) {
                if(m_atoms[j].kind==Kind::Hydrocarbon) stable=false;
                if(!carbon(j) && !visited[j]) {visited[j]=true;cavity.push_back(j);}
            }
            if(stable) for(auto id:cavity) if(sites[id].position[2]>m_initialTop+1e-7) m_frozenVacancy[id]=true;
        }
        for(std::size_t h=sites.size();h<m_atoms.size();++h) if(m_atoms[h].kind==Kind::Hydrogen && m_atoms[m_atoms[h].parent].kind==Kind::Diamond && !has(m_surface,m_atoms[h].parent) ) {
            const auto parent=m_atoms[h].parent;const auto target=sites[parent].neighbors[(h-sites.size())%4];
            const auto d=Atomistic::displacement(m_options.lattice,sites[parent].position,sites[target].position);
            const double length=std::sqrt(std::inner_product(d.begin(),d.end(),d.begin(),0.0));
            if(sites[parent].position[2]+m_options.chLength*d[2]/length<=m_initialTop+1e-7) continue;
            bool stable=true;for(auto j:sites[m_atoms[h].parent].neighbors) if(m_atoms[j].kind==Kind::Hydrocarbon) stable=false;
            if(stable) m_atoms[h].trapped=true;
        }
    }
}
bool Simulation::sixRing(std::size_t a,std::size_t b) const {
    // Diamond's first-neighbor graph has no shorter carbon rings; a length-five
    // existing bond path plus the new edge completes the six-carbon ring.
    std::vector<int> distance(m_atoms.size(),-1);std::vector<std::size_t> pending{a};distance[a]=0;
    for(std::size_t k=0;k<pending.size();++k) {auto i=pending[k];if(distance[i]>=5) continue;for(auto j:m_bonds[i]) if(distance[j]<0) {distance[j]=distance[i]+1;pending.push_back(j);}}
    return distance[b]==5;
}
void Simulation::rebuild() {
    m_events.clear();const auto& sites=m_options.lattice.sites;const auto n=sites.size();
    struct Group {std::vector<std::size_t> carbons,H;int anchors=0;};
    std::vector<Group> groups;std::vector<std::size_t> groupOf(n,n);
    for(auto start:m_surface) if(m_atoms[start].kind==Kind::Hydrocarbon && groupOf[start]==n) {
        Group g;g.carbons={start};groupOf[start]=groups.size();
        for(std::size_t k=0;k<g.carbons.size();++k) {auto i=g.carbons[k];for(auto h:hydrogens(i)) g.H.push_back(h);for(auto j:m_bonds[i]) {if(m_atoms[j].kind==Kind::Diamond) ++g.anchors;else if(m_atoms[j].kind==Kind::Hydrocarbon && groupOf[j]==n) {groupOf[j]=groups.size();g.carbons.push_back(j);}}}
        groups.push_back(std::move(g));
    }
    auto empty=[&](std::size_t id) {if(m_atoms[id].kind!=Kind::Empty || m_frozenVacancy[id] || !m_gas[id]) return false;for(const auto& d:m_dimers) if(d.shared==id) return false;return true;};
    auto domain=[&](std::size_t id) {if(m_atoms[id].kind==Kind::Diamond) return std::string("diamond");if(groupOf[id]==n) return std::string();return "C"+std::to_string(groups[groupOf[id]].carbons.size());};
    auto emit=[&](std::size_t process,Event event,double factor=1.0) {event.process=process;event.rate=m_options.processes[process].rate*factor;if(event.rate>0) m_events.push_back(std::move(event));};
    auto hydrogenPlacements=[&](const std::vector<std::size_t>& carbons,const std::vector<std::size_t>& reserved,int H,const std::vector<std::size_t>& erase) {
        std::vector<std::pair<std::size_t,std::size_t>> available;
        for(auto c:carbons) for(auto target:sites[c].neighbors) if(!has(reserved,target)) {auto h=hydrogenSlot(c,target);if(m_atoms[h].kind==Kind::Empty || has(erase,h)) available.emplace_back(h,c);}
        std::vector<std::vector<std::pair<std::size_t,std::size_t>>> result;std::vector<std::pair<std::size_t,std::size_t>> choice;
        std::function<void(std::size_t)> pick=[&](std::size_t index) {if(static_cast<int>(choice.size())==H) {result.push_back(choice);return;}for(std::size_t k=index;k<available.size();++k) {bool duplicate=false;for(auto h:choice) duplicate|=h.first==available[k].first;if(duplicate) continue;choice.push_back(available[k]);pick(k+1);choice.pop_back();}};pick(0);return result;
    };
    for(std::size_t p=0;p<m_options.processes.size();++p) {
        const auto& rule=m_options.processes[p];if(!rule.rate) continue;
        if(rule.event=="hydrogen_addition" || rule.event=="hydrogen_removal") {
            for(auto c:m_surface) if(domain(c)==rule.domain) {
                if(rule.event=="hydrogen_removal") {for(auto h:hydrogens(c)) if(exposed(h)) {Event e;e.erase={h};emit(p,e);}}
                else for(auto target:sites[c].neighbors) if(free(c,target)) {Event e;e.hydrogens={{hydrogenSlot(c,target),c}};emit(p,e);}
            }
        } else if(rule.event=="adsorb_carbon" || rule.event=="append_methyl") {
            for(auto c:m_surface) if((rule.event=="adsorb_carbon" && m_atoms[c].kind==Kind::Diamond) || (rule.event=="append_methyl" && domain(c)=="C1")) for(auto j:sites[c].neighbors) if(empty(j) && free(c,j)) {
                std::vector<Event> choices;
                if(rule.domain=="C1" || rule.event=="append_methyl") {
                    const int H=rule.event=="append_methyl" ? 3 : rule.hydrogens;
                    // A methyl has three C-H bonds; boundary-incomplete sites
                    // cannot silently accept a truncated molecule.
                    if(sites[j].neighbors.size()!=4) {m_capacity=true;continue;}
                    for(auto placement:hydrogenPlacements({j},{c,j},H,{})) {Event e;e.carbons={j};e.hydrogens=std::move(placement);e.bonds={{c,j}};choices.push_back(std::move(e));}
                } else for(auto k:sites[j].neighbors) if(k!=c && empty(k) && sites[j].neighbors.size()==4 && sites[k].neighbors.size()==4) {
                    // Sigma-lattice C2H2 embedding: one H on each carbon.
                    for(auto hj:sites[j].neighbors) if(hj!=c && hj!=k) for(auto hk:sites[k].neighbors) if(hk!=j) {Event e;e.carbons={j,k};e.hydrogens={{hydrogenSlot(j,hj),j},{hydrogenSlot(k,hk),k}};e.bonds={{c,j},{j,k}};choices.push_back(std::move(e));}
                }
                for(auto& e:choices) emit(p,std::move(e),1.0/choices.size());
            }
        } else if(rule.event=="desorb_hydrocarbon") {
            for(const auto& g:groups) if("C"+std::to_string(g.carbons.size())==rule.domain && static_cast<int>(g.H.size())==rule.hydrogens && g.anchors==1) {
                bool accessible=true;for(auto h:g.H) accessible&=!m_atoms[h].trapped;
                if(accessible) {Event e;e.erase=g.carbons;e.erase.insert(e.erase.end(),g.H.begin(),g.H.end());emit(p,e);}
            }
        } else if(rule.event=="carbon_incorporation") {
            // Ring completion joins radical hydrocarbons to film or another
            // anchored hydrocarbon; relabel the entire ring-connected molecule.
            for(auto a:m_surface) if(carbon(a)) for(auto b:sites[a].neighbors) if(b>a && carbon(b) && !has(m_bonds[a],b) && free(a,b) && free(b,a) && (m_atoms[a].kind==Kind::Hydrocarbon || m_atoms[b].kind==Kind::Hydrocarbon) && sixRing(a,b)) {
                Event e;e.bonds={{a,b}};std::set<std::size_t> molecules;
                if(groupOf[a]!=n) molecules.insert(groupOf[a]);if(groupOf[b]!=n) molecules.insert(groupOf[b]);int anchors=0;
                for(auto g:molecules) {anchors+=groups[g].anchors;e.diamond.insert(e.diamond.end(),groups[g].carbons.begin(),groups[g].carbons.end());}
                if(m_atoms[a].kind==Kind::Diamond || m_atoms[b].kind==Kind::Diamond) ++anchors;
                if(anchors>=2) emit(p,e);
            }
            // Insertion/rotation into an opened dimer site (Table II, row 12).
            for(std::size_t x=0;x<m_surface.size();++x) {auto a=m_surface[x];if(m_atoms[a].kind!=Kind::Diamond) continue;
                for(auto shared:sites[a].neighbors) if(empty(shared) && free(a,shared)) for(auto b:sites[shared].neighbors) if(b>a && m_atoms[b].kind==Kind::Diamond && free(b,shared)) {
                    for(auto attached:{a,b}) for(auto old:m_bonds[attached]) if(m_atoms[old].kind==Kind::Hydrocarbon && groupOf[old]!=n) {
                        const auto& g=groups[groupOf[old]];if(g.carbons.size()>2 || g.anchors!=1) continue;
                        const auto H=g.H.size();const auto occupied=m_bonds[old].size()+hydrogens(old).size();if(occupied>=4) continue;
                        std::vector<std::size_t> erased=g.carbons;erased.insert(erased.end(),g.H.begin(),g.H.end());
                        bool trapped=false;for(auto h:g.H) trapped|=m_atoms[h].trapped;if(trapped) continue;
                        std::vector<Event> choices;
                        if(g.carbons.size()==1 && H<=2) for(auto h:hydrogenPlacements({shared},{a,b,shared},static_cast<int>(H),erased)) {Event e;e.erase=erased;e.carbons={shared};e.diamond={shared};e.bonds={{a,shared},{b,shared}};e.hydrogens=std::move(h);choices.push_back(std::move(e));}
                        else if(g.carbons.size()==2 && H<=3) for(auto rest:sites[shared].neighbors) if(rest!=a && rest!=b && (empty(rest) || has(erased,rest))) for(auto h:hydrogenPlacements({rest},{a,b,shared,rest},static_cast<int>(H),erased)) {Event e;e.erase=erased;e.carbons={shared,rest};e.diamond={shared};e.bonds={{a,shared},{b,shared},{shared,rest}};e.hydrogens=std::move(h);choices.push_back(std::move(e));}
                        for(auto& e:choices) emit(p,std::move(e),1.0/choices.size());
                    }
                }
            }
        } else if(rule.event=="dimer_closing") {
            for(auto a:m_surface) if(m_atoms[a].kind==Kind::Diamond && m_bonds[a].size()<=2) for(auto shared:sites[a].neighbors) if(empty(shared) && free(a,shared)) for(auto b:sites[shared].neighbors) if(b>a && m_atoms[b].kind==Kind::Diamond && m_bonds[b].size()<=2 && free(b,shared)) {
                std::size_t da=m_bonds[a].size(),db=m_bonds[b].size();for(auto d:m_dimers) {da+=(d.first==a || d.second==a);db+=(d.first==b || d.second==b);}if(da>2 || db>2) continue;
                Event e;e.close={{a,b,shared}};emit(p,e);
            }
        } else if(rule.event=="dimer_opening") {
            for(const auto& d:m_dimers) {
                bool radical=false;for(auto a:{d.first,d.second}) for(auto c:sites[a].neighbors) if(m_atoms[c].kind==Kind::Hydrocarbon && m_bonds[c].size()+hydrogens(c).size()<4) radical=true;
                if(radical) {Event e;e.open={d};emit(p,e);}
            }
        }
    }
}
void Simulation::apply(const Event& e) {
    bool carbonChange=false;
    for(auto i:e.erase) {if(m_atoms[i].kind==Kind::Diamond) --m_diamondCount;if(carbon(i)) {carbonChange=true;for(auto j:m_bonds[i]) {auto& back=m_bonds[j];back.erase(std::remove(back.begin(),back.end(),i),back.end());}m_bonds[i].clear();}m_atoms[i]=Atom{};}
    for(auto i:e.carbons) {m_atoms[i]={Kind::Hydrocarbon,0,false};carbonChange=true;}
    for(auto h:e.hydrogens) m_atoms[h.first]={Kind::Hydrogen,h.second,false};
    for(auto bond:e.bonds) {m_bonds[bond.first].push_back(bond.second);m_bonds[bond.second].push_back(bond.first);}
    for(auto i:e.diamond) {if(m_atoms[i].kind!=Kind::Diamond) ++m_diamondCount;m_atoms[i].kind=Kind::Diamond;}
    for(auto d:e.open) m_dimers.erase(std::remove_if(m_dimers.begin(),m_dimers.end(),[&](const Dimer& x){return x.first==d.first && x.second==d.second && x.shared==d.shared;}),m_dimers.end());
    m_dimers.insert(m_dimers.end(),e.close.begin(),e.close.end());
    ++m_steps;++m_counts[e.process];
    if(carbonChange || !e.diamond.empty()) gasAccess(!e.diamond.empty());
    rebuild();
}
void Simulation::execute(std::size_t index) {auto event=m_events.at(index);apply(event);}
std::size_t Simulation::count(Kind kind) const {return static_cast<std::size_t>(std::count_if(m_atoms.begin(),m_atoms.end(),[&](const Atom& a){return a.kind==kind;}));}
std::size_t Simulation::trappedHydrogens() const {return static_cast<std::size_t>(std::count_if(m_atoms.begin(),m_atoms.end(),[](const Atom& a){return a.kind==Kind::Hydrogen && a.trapped;}));}
std::size_t Simulation::trappedVacancies() const {return static_cast<std::size_t>(std::count(m_frozenVacancy.begin(),m_frozenVacancy.end(),true));}
double Simulation::meanHeight() const {return m_diamondCount*m_options.siteVolume/m_area;}
double Simulation::growthRate() const {return m_time>m_options.start ? (meanHeight()-m_initialHeight)/(m_time-m_options.start) : 0;}
double Simulation::timeAverageHeight() const {return m_time>m_options.start ? static_cast<double>(m_integral/(m_preciseTime-m_options.start)) : m_initialHeight;}
bool Simulation::consistent() const {
    const auto n=m_options.lattice.sites.size();std::vector<bool> anchored(n,false);std::vector<std::size_t> pending;
    if(count(Kind::Diamond)!=m_diamondCount) return false;
    for(std::size_t i=0;i<n;++i) {
        if(carbon(i)) {
            if(m_bonds[i].empty()) return false;
            std::set<std::size_t> unique;std::size_t valence=m_bonds[i].size()+hydrogens(i).size();for(const auto& d:m_dimers) valence+=(d.first==i || d.second==i);
            if(valence>4) return false;
            for(auto j:m_bonds[i]) if(!carbon(j) || !has(m_options.lattice.sites[i].neighbors,j) || !has(m_bonds[j],i) || !unique.insert(j).second || m_atoms[hydrogenSlot(i,j)].kind==Kind::Hydrogen) return false;
            if(m_options.lattice.sites[i].fixed) {anchored[i]=true;pending.push_back(i);}
        } else if(!m_bonds[i].empty() || !hydrogens(i).empty()) return false;
        if(m_options.lattice.sites[i].fixed && m_atoms[i].kind!=Kind::Diamond) return false;
    }
    for(std::size_t i=n;i<m_atoms.size();++i) if(m_atoms[i].kind!=Kind::Empty) {
        if(m_atoms[i].kind!=Kind::Hydrogen || m_atoms[i].parent>=n || !carbon(m_atoms[i].parent) || (i-n)/4!=m_atoms[i].parent || (i-n)%4>=m_options.lattice.sites[m_atoms[i].parent].neighbors.size()) return false;
    }
    std::set<std::pair<std::size_t,std::size_t>> dimers;
    for(const auto& d:m_dimers) {
        if(d.first>=n || d.second>=n || d.shared>=n || m_atoms[d.first].kind!=Kind::Diamond || m_atoms[d.second].kind!=Kind::Diamond || !has(m_options.lattice.sites[d.first].neighbors,d.shared) || !has(m_options.lattice.sites[d.second].neighbors,d.shared) || carbon(d.shared) || !dimers.emplace(d.first,d.second).second) return false;
    }
    if(!pending.empty()) {for(std::size_t k=0;k<pending.size();++k) for(auto j:m_bonds[pending[k]]) if(!anchored[j]) {anchored[j]=true;pending.push_back(j);}for(std::size_t i=0;i<n;++i) if(carbon(i) && !anchored[i]) return false;}
    return true;
}
void Simulation::run(const std::function<void(const Simulation&,bool)>& observe) {
    if(observe) observe(*this,false);
    while(m_preciseTime<m_options.end) {
        if(m_capacity) {m_stop="capacity_limit";break;}
        if(m_options.maxEvents && m_steps>=m_options.maxEvents) {m_stop="event_limit";break;}
        double total=0;for(const auto& e:m_events) total+=e.rate;if(!std::isfinite(total)) throw std::runtime_error("diamond total rate overflow");if(!total) {m_stop="absorbing";break;}
        double u=0;do {u=std::generate_canonical<double,53>(m_rng);}while(!u);const double dt=-std::log(u)/total;
        if(dt>m_options.end-m_preciseTime) {m_integral+=meanHeight()*(m_options.end-m_preciseTime);m_preciseTime=m_options.end;m_time=m_options.end;break;}
        if(!(dt>0) || !std::isfinite(dt)) throw std::runtime_error("invalid diamond time step");
        const double pick=std::generate_canonical<double,53>(m_rng)*total;double sum=0;std::size_t selected=m_events.size()-1;for(std::size_t i=0;i<m_events.size();++i) {sum+=m_events[i].rate;if(pick<sum) {selected=i;break;}}
        m_integral+=meanHeight()*dt;
        const long double increment=static_cast<long double>(dt)-m_timeCorrection;
        const long double advanced=m_preciseTime+increment;
        m_timeCorrection=(advanced-m_preciseTime)-increment;m_preciseTime=advanced;m_time=static_cast<double>(m_preciseTime);
        execute(selected);if(observe) observe(*this,false);
    }
    if(m_stop.empty()) m_stop=m_capacity ? "capacity_limit" : "end_time";if(observe) observe(*this,true);
}
void Simulation::writeXYZ(std::ostream& out) const {
    const auto& lattice=m_options.lattice;out<<std::setprecision(17)<<count(Kind::Diamond)+count(Kind::Hydrocarbon)+count(Kind::Hydrogen)<<"\nLattice=\"";
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) out<<(i||j ? " " : "")<<lattice.cell[i][j];
    out<<"\" Properties=species:S:1:pos:R:3:site_id:I:1:fixed:L:1:chemical_state:S:1:parent:I:1:trapped:L:1 pbc=\"";for(int i=0;i<3;++i) out<<(i ? " " : "")<<(lattice.periodic[i] ? "T" : "F");out<<"\" Time="<<m_time<<'\n';
    for(std::size_t i=0;i<m_atoms.size();++i) if(m_atoms[i].kind!=Kind::Empty) {
        const auto& a=m_atoms[i];const bool H=a.kind==Kind::Hydrogen;
        auto r=H ? lattice.sites[lattice.sites[a.parent].neighbors[(i-lattice.sites.size())%4]].position : lattice.sites[i].position;
        if(H) {auto d=Atomistic::displacement(lattice,lattice.sites[a.parent].position,r);double length=std::sqrt(std::inner_product(d.begin(),d.end(),d.begin(),0.0));for(int k=0;k<3;++k) {r[k]=lattice.sites[a.parent].position[k]+m_options.chLength*d[k]/length;if(lattice.periodic[k]) r[k]-=std::floor(r[k]/lattice.cell[k][k])*lattice.cell[k][k];}}
        out<<(H ? "H" : "C")<<' '<<r[0]<<' '<<r[1]<<' '<<r[2]<<' '<<i<<' '<<(!H && lattice.sites[i].fixed ? "T" : "F")<<' '<<(H ? "hydrogen" : a.kind==Kind::Diamond ? "diamond" : "hydrocarbon")<<' '<<(H ? static_cast<long long>(a.parent) : -1)<<' '<<(a.trapped ? "T" : "F")<<'\n';
    }
    if(!out) throw std::runtime_error("failed to write diamond XYZ");
}
void runToFiles(const Options& options) {
    Simulation simulation(options);const auto initialDiamond=simulation.count(Kind::Diamond);std::ofstream log("diamond.log"),trajectory;if(!log) throw std::runtime_error("cannot open diamond.log");
    if(options.writeXYZ) {trajectory.open("diamond_trajectory.xyz");std::ofstream initial("diamond_initial.xyz");if(!initial || !trajectory) throw std::runtime_error("cannot open diamond XYZ outputs");simulation.writeXYZ(initial);}
    log<<std::setprecision(17)<<"time\tdiamond_C\thydrocarbon_C\tH\ttrapped_H\ttrapped_vacancies\tgrowth_rate_um_hr\tevents\n";double last=-std::numeric_limits<double>::infinity(),frame=last;
    simulation.run([&](const Simulation& s,bool final) {if((final && s.time()!=last) || s.time()-last>=options.logInterval) {log<<s.time()<<'\t'<<s.count(Kind::Diamond)<<'\t'<<s.count(Kind::Hydrocarbon)<<'\t'<<s.count(Kind::Hydrogen)<<'\t'<<s.trappedHydrogens()<<'\t'<<s.trappedVacancies()<<'\t'<<s.growthRate()*.36<<'\t'<<s.steps()<<'\n';last=s.time();}if(options.writeXYZ && ((final && s.time()!=frame) || !std::isfinite(frame) || (options.snapshotInterval>0 && s.time()-frame>=options.snapshotInterval))) {s.writeXYZ(trajectory);frame=s.time();}});
    nlohmann::json counts=nlohmann::json::object(),rates=nlohmann::json::object();for(std::size_t i=0;i<options.processes.size();++i) {counts[options.processes[i].id]=simulation.counts()[i];rates[options.processes[i].id]=options.processes[i].rate;}
    const auto deposited=simulation.count(Kind::Diamond)-initialDiamond;
    const double defectDenominator=static_cast<double>(deposited+simulation.trappedHydrogens()+simulation.trappedVacancies());
    nlohmann::json summary={{"model","diamond_cvd_1997"},{"initial_diamond_carbon_atoms",initialDiamond},{"deposited_diamond_carbon_atoms",deposited},{"trapped_hydrogen_fraction",defectDenominator ? simulation.trappedHydrogens()/defectDenominator : 0},{"trapped_vacancy_fraction",defectDenominator ? simulation.trappedVacancies()/defectDenominator : 0},{"temperature_K",options.temperature},{"gas_temperature_K",options.gasTemperature},{"partial_pressures_Torr",options.pressuresTorr},{"time",simulation.time()},{"events",simulation.steps()},{"stop_reason",simulation.stopReason()},{"diamond_carbon_atoms",simulation.count(Kind::Diamond)},{"hydrocarbon_carbon_atoms",simulation.count(Kind::Hydrocarbon)},{"hydrogen_atoms",simulation.count(Kind::Hydrogen)},{"trapped_hydrogen_atoms",simulation.trappedHydrogens()},{"trapped_vacancies",simulation.trappedVacancies()},{"dimers",simulation.dimers().size()},{"mean_height_angstrom",simulation.meanHeight()},{"time_average_height_angstrom",simulation.timeAverageHeight()},{"mean_growth_rate_um_hr",simulation.growthRate()*.36},{"process_counts",counts},{"process_rates_s_inverse",rates},{"seed",simulation.options().seed},{"numerical_reproduction_validated",false}};
    std::ofstream result("diamond_summary.json");result<<summary.dump(2)<<'\n';result.close();log.close();
    if(options.writeXYZ) {std::ofstream final("diamond_final.xyz"),bonds("diamond_bonds.tsv");simulation.writeXYZ(final);bonds<<"site1\tsite2\tkind\n";for(std::size_t i=0;i<simulation.bonds().size();++i) for(auto j:simulation.bonds()[i]) if(j>i) bonds<<i<<'\t'<<j<<"\tcarbon_sigma\n";for(const auto& d:simulation.dimers()) bonds<<d.first<<'\t'<<d.second<<"\tdimer\n";for(std::size_t i=0;i<simulation.atoms().size();++i) if(simulation.atoms()[i].kind==Kind::Hydrogen) bonds<<simulation.atoms()[i].parent<<'\t'<<i<<"\tcarbon_hydrogen\n";final.close();bonds.close();trajectory.close();if(!final || !bonds || !trajectory) throw std::runtime_error("failed to write diamond geometry");}
    if(!result || !log) throw std::runtime_error("failed to write diamond results");
}
}
