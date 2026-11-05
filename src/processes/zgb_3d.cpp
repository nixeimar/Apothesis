#include "zgb_3d.h"
#include "../extLibs/nlohmann/json.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ZGB3D {
void validate(const Options& o) {
    Atomistic::validate(o.geometry);ZGB::validate(o.kinetics);
    const auto& g=o.geometry;
    if (g.crystal!="SimpleCubic" || !g.periodic[0] || !g.periodic[1] || g.periodic[2])
        throw std::runtime_error("ZGB 3D requires SimpleCubic, lateral periodic boundaries and open z");
    if (g.repetitions[0]!=o.kinetics.x || g.repetitions[1]!=o.kinetics.y)
        throw std::runtime_error("ZGB 3D lateral dimensions do not match geometry");
    if (g.neighborCutoff && (g.neighborCutoff<g.latticeConstant || g.neighborCutoff>=std::sqrt(2.0)*g.latticeConstant))
        throw std::runtime_error("ZGB 3D neighbor_cutoff must include only the first SC neighbor shell");
    if (!std::isfinite(o.snapshotInterval) || o.snapshotInterval<0)
        throw std::runtime_error("output.lattice_interval must be finite and nonnegative");
}
void Simulation::IndexedSet::set(std::size_t id,bool present) {
    auto& position=positions.at(id);
    if (present && position<0) {position=static_cast<long long>(ids.size());ids.push_back(id);}
    else if (!present && position>=0) {
        auto replacement=ids.back();ids[static_cast<std::size_t>(position)]=replacement;
        positions[replacement]=position;ids.pop_back();position=-1;
    }
}
Simulation::Simulation(Options options):m_options(std::move(options)) {
    validate(m_options);
    m_solid=Atomistic::build(m_options.geometry);
    const auto n=m_solid.sites.size();
    m_adsorbates.assign(n,ZGB::Adsorbate::Vacant);m_binding.assign(n,false);
    m_lateral.resize(n);m_incident.resize(n);m_edges.resize(2*n);
    m_up.assign(n,n);m_down.assign(n,n);
    m_vacancies=IndexedSet(n);m_pairs=IndexedSet(2*n);m_reactions=IndexedSet(2*n);m_reactiveSites=IndexedSet(n);
    const auto a=m_options.geometry.latticeConstant,tolerance=a*1e-7;
    // Derive axial neighbors from explicit first-neighbor geometry, not a height map.
    for (const auto& site:m_solid.sites) {
        std::array<std::size_t,4> lateral{{n,n,n,n}};
        for (auto j:site.neighbors) {
            Atomistic::Vec3 d;
            for (int k=0;k<3;++k) {
                d[k]=m_solid.sites[j].position[k]-site.position[k];
                if (m_solid.periodic[k]) d[k]-=std::round(d[k]/m_solid.cell[k][k])*m_solid.cell[k][k];
            }
            if (std::abs(d[2])<tolerance) {
                if (std::abs(d[0]-a)<tolerance && std::abs(d[1])<tolerance) lateral[0]=j;
                else if (std::abs(d[0]+a)<tolerance && std::abs(d[1])<tolerance) lateral[1]=j;
                else if (std::abs(d[1]-a)<tolerance && std::abs(d[0])<tolerance) lateral[2]=j;
                else if (std::abs(d[1]+a)<tolerance && std::abs(d[0])<tolerance) lateral[3]=j;
                else throw std::runtime_error("ZGB 3D found a nonaxial lateral neighbor");
            } else if (std::abs(d[2]-a)<tolerance) m_up[site.id]=j;
            else if (std::abs(d[2]+a)<tolerance) m_down[site.id]=j;
            else throw std::runtime_error("ZGB 3D found an invalid vertical neighbor");
        }
        if (std::find(lateral.begin(),lateral.end(),n)!=lateral.end()) throw std::runtime_error("ZGB 3D requires four lateral first neighbors");
        m_lateral[site.id]=lateral;
        m_edges[2*site.id]={{site.id,lateral[0]}};m_edges[2*site.id+1]={{site.id,lateral[2]}};
        m_incident[site.id]={{2*site.id,2*lateral[1],2*site.id+1,2*lateral[3]+1}};
    }
    // The paper's surface is a supported slab: reject buried cavities and overhangs.
    // Walking vertical graph links initializes one real binding node per surface chain.
    for (const auto& bottom:m_solid.sites) if (m_down[bottom.id]==n) {
        bool foundEmpty=false,foundSolid=false;
        for (auto id=bottom.id;id!=n;id=m_up[id]) {
            auto& site=m_solid.sites[id];
            if (site.occupied) {
                if (foundEmpty) throw std::runtime_error("ZGB 3D initial solid contains a buried vacancy or unsupported overhang");
                foundSolid=true;++m_solidCount;
            } else if (!foundEmpty) {
                if (!foundSolid || site.fixed) throw std::runtime_error("ZGB 3D needs a nonempty substrate and an unfixed binding site");
                m_binding[id]=true;foundEmpty=true;
            }
        }
        if (!foundEmpty) throw std::runtime_error("ZGB 3D needs empty crystal layers above every initial surface site");
    }
    if (!m_options.kinetics.growth) {
        double level=0;bool first=true;
        for (std::size_t id=0;id<n;++id) if (m_binding[id]) {
            if (first) {level=m_solid.sites[id].position[2];first=false;}
            else if (std::abs(m_solid.sites[id].position[2]-level)>tolerance)
                throw std::runtime_error("ZGB 3D catalytic reference requires a flat initial binding surface");
        }
    }
    const auto columns=static_cast<double>(m_options.kinetics.x)*m_options.kinetics.y;
    m_initialMean=m_solidCount/columns;m_time=m_options.kinetics.start;
    if (!m_options.kinetics.seed) {std::random_device rd;m_options.kinetics.seed=(static_cast<std::uint64_t>(rd())<<32)^rd();}
    m_rng.seed(m_options.kinetics.seed);
    for (std::size_t id=0;id<n;++id) refresh(id);
}
bool Simulation::pairEligible(std::size_t edge) const {
    const auto& e=m_edges[edge];
    return m_binding[e[0]] && m_binding[e[1]] && m_adsorbates[e[0]]==ZGB::Adsorbate::Vacant && m_adsorbates[e[1]]==ZGB::Adsorbate::Vacant;
}
bool Simulation::reactionEligible(std::size_t edge) const {
    const auto& e=m_edges[edge];
    if (!m_binding[e[0]] || !m_binding[e[1]]) return false;
    return (m_adsorbates[e[0]]==ZGB::Adsorbate::AB && m_adsorbates[e[1]]==ZGB::Adsorbate::B) ||
           (m_adsorbates[e[1]]==ZGB::Adsorbate::AB && m_adsorbates[e[0]]==ZGB::Adsorbate::B);
}
void Simulation::refresh(std::size_t id) {
    if (id==m_solid.sites.size()) return;
    m_vacancies.set(id,m_binding[id] && m_adsorbates[id]==ZGB::Adsorbate::Vacant);
    for (auto e:m_incident[id]) {m_pairs.set(e,pairEligible(e));m_reactions.set(e,reactionEligible(e));}
    auto active=[&](std::size_t i) {
        bool eligible=false;for (auto e:m_incident[i]) if (reactionEligible(e)) {eligible=true;break;}
        m_reactiveSites.set(i,eligible);
    };
    active(id);for (auto j:m_lateral[id]) active(j);
}
double Simulation::pairRate() const {
    return (1-m_options.kinetics.yAB)*(m_options.kinetics.pairConvention=="zgb_attempts" ? 0.5 : 2.0);
}
void Simulation::apply(ZGB::Event event,std::size_t first,std::size_t second) {
    if (!m_binding.at(first)) throw std::runtime_error("ZGB 3D event requires an exposed binding node");
    if (event==ZGB::Event::AdsorbAB) {
        if (m_adsorbates[first]!=ZGB::Adsorbate::Vacant) throw std::runtime_error("AB requires a vacant binding node");
        m_adsorbates[first]=ZGB::Adsorbate::AB;++m_adsAB;
    } else {
        if (!m_binding.at(second) || std::find(m_lateral[first].begin(),m_lateral[first].end(),second)==m_lateral[first].end())
            throw std::runtime_error("ZGB 3D pair requires lateral, same-height first-neighbor binding nodes");
        if (event==ZGB::Event::AdsorbB2) {
            if (m_adsorbates[first]!=ZGB::Adsorbate::Vacant || m_adsorbates[second]!=ZGB::Adsorbate::Vacant)
                throw std::runtime_error("B2 requires two vacant binding nodes");
            m_adsorbates[first]=m_adsorbates[second]=ZGB::Adsorbate::B;++m_adsB2;
        } else {
            auto precursor=first;
            if (m_adsorbates[first]==ZGB::Adsorbate::B && m_adsorbates[second]==ZGB::Adsorbate::AB) precursor=second;
            else if (!(m_adsorbates[first]==ZGB::Adsorbate::AB && m_adsorbates[second]==ZGB::Adsorbate::B))
                throw std::runtime_error("reaction requires AB* and B*");
            if (m_options.kinetics.growth) {
                auto& unit=m_solid.sites[precursor];
                if (unit.fixed || unit.occupied || m_down[precursor]==m_solid.sites.size() || !m_solid.sites[m_down[precursor]].occupied)
                    throw std::runtime_error("solid incorporation requires an empty, supported, unfixed 3D site");
                unit.occupied=true;unit.species="AB2";++m_solidCount;
                m_binding[precursor]=false;
                const auto above=m_up[precursor];
                if (above==m_solid.sites.size()) m_capacity=true;
                else {
                    if (m_solid.sites[above].occupied || m_solid.sites[above].fixed) throw std::runtime_error("invalid 3D growth headroom");
                    m_binding[above]=true;refresh(above);
                }
            }
            m_adsorbates[first]=m_adsorbates[second]=ZGB::Adsorbate::Vacant;++m_reacted;
        }
        refresh(second);
    }
    refresh(first);++m_steps;
}
void Simulation::run(const std::function<void(const Simulation&,bool)>& observe) {
    const auto& k=m_options.kinetics;
    if (observe) observe(*this,false);
    while (m_time<k.end) {
        if (m_capacity) {m_stop="capacity_limit";break;}
        const double a=k.yAB*m_vacancies.ids.size(),b=pairRate()*m_pairs.ids.size();
        const double r=k.reactionRate*(k.reactionSelection=="sites" ? m_reactiveSites.ids.size() : m_reactions.ids.size());
        const double total=a+b+r;
        if (!std::isfinite(total)) throw std::runtime_error("ZGB 3D total rate overflow");
        if (!total) {m_stop="absorbing";break;}
        if (k.maxEvents && m_steps>=k.maxEvents) {m_stop="event_limit";break;}
        double u=0;do {u=std::generate_canonical<double,53>(m_rng);} while (!u);
        const double dt=-std::log(u)/total;
        const double mean=m_solidCount/(static_cast<double>(k.x)*k.y);
        if (dt>k.end-m_time) {m_integral+=mean*(k.end-m_time);m_time=k.end;break;}
        if (!(dt>0) || m_time+dt<=m_time) throw std::runtime_error("ZGB 3D time step below floating-point resolution");
        const double pick=std::generate_canonical<double,53>(m_rng)*total;
        auto choose=[&](const IndexedSet& set) {return set.ids[std::uniform_int_distribution<std::size_t>(0,set.ids.size()-1)(m_rng)];};
        m_integral+=mean*dt;m_time+=dt;
        if (pick<a || (b==0 && r==0)) apply(ZGB::Event::AdsorbAB,choose(m_vacancies));
        else if (pick<a+b || r==0) {auto e=m_edges[choose(m_pairs)];apply(ZGB::Event::AdsorbB2,e[0],e[1]);}
        else if (k.reactionSelection=="pairs") {auto e=m_edges[choose(m_reactions)];apply(ZGB::Event::React,e[0],e[1]);}
        else {
            auto center=choose(m_reactiveSites);std::array<std::size_t,4> partners{};std::size_t count=0;
            for (auto edge:m_incident[center]) if (reactionEligible(edge)) {
                const auto& e=m_edges[edge];partners[count++]=e[0]==center ? e[1] : e[0];
            }
            apply(ZGB::Event::React,center,partners[std::uniform_int_distribution<std::size_t>(0,count-1)(m_rng)]);
        }
        if (observe) observe(*this,false);
    }
    if (m_stop.empty()) m_stop=m_capacity ? "capacity_limit" : "end_time";
    if (observe) observe(*this,true);
}
ZGB::Statistics Simulation::statistics() const {
    ZGB::Statistics s;const auto& k=m_options.kinetics;
    const double columns=static_cast<double>(k.x)*k.y;
    s.time=m_time;s.events=m_steps;s.adsorbedAB=m_adsAB;s.adsorbedB2=m_adsB2;s.reactions=m_reacted;
    s.coverageAB=(m_adsAB-m_reacted)/columns;s.coverageB=(2*m_adsB2-m_reacted)/columns;s.coverageVacant=m_vacancies.ids.size()/columns;
    s.meanHeight=m_solidCount/columns;
    s.timeAverageHeight=m_time>k.start ? static_cast<double>(m_integral/(m_time-k.start)) : m_initialMean;
    if (m_time>k.start) s.growthRate=(s.meanHeight-m_initialMean)/(m_time-k.start);
    s.vacantSites=m_vacancies.ids.size();s.vacantPairs=m_pairs.ids.size();s.reactionPairs=m_reactions.ids.size();s.reactiveSites=m_reactiveSites.ids.size();
    s.poisoned=!m_capacity && k.yAB*s.vacantSites+pairRate()*s.vacantPairs+k.reactionRate*(k.reactionSelection=="sites" ? s.reactiveSites : s.reactionPairs)==0;
    s.stopReason=m_stop;
    const auto n=m_solid.sites.size();
    for (const auto& bottom:m_solid.sites) if (m_down[bottom.id]==n) {
        double occupied=0;for (auto id=bottom.id;id!=n;id=m_up[id]) occupied+=m_solid.sites[id].occupied ? 1 : 0;
        s.rms+=(occupied-s.meanHeight)*(occupied-s.meanHeight)/columns;
    }
    s.rms=std::sqrt(s.rms);return s;
}
bool Simulation::catalogConsistent() const {
    const auto n=m_solid.sites.size();std::size_t solids=0;
    for (const auto& site:m_solid.sites) {
        solids+=site.occupied ? 1 : 0;
        if (site.occupied && m_binding[site.id]) return false;
        if (!m_binding[site.id] && m_adsorbates[site.id]!=ZGB::Adsorbate::Vacant) return false;
        if (m_binding[site.id] && (m_down[site.id]==n || !m_solid.sites[m_down[site.id]].occupied)) return false;
        if ((m_vacancies.positions[site.id]>=0)!=(m_binding[site.id] && m_adsorbates[site.id]==ZGB::Adsorbate::Vacant)) return false;
        bool active=false;for (auto e:m_incident[site.id]) if (reactionEligible(e)) active=true;
        if ((m_reactiveSites.positions[site.id]>=0)!=active) return false;
    }
    if (solids!=m_solidCount) return false;
    for (std::size_t e=0;e<m_edges.size();++e)
        if ((m_pairs.positions[e]>=0)!=pairEligible(e) || (m_reactions.positions[e]>=0)!=reactionEligible(e)) return false;
    for (const auto* set:{&m_vacancies,&m_pairs,&m_reactions,&m_reactiveSites}) {
        for (std::size_t i=0;i<set->ids.size();++i) if (set->positions[set->ids[i]]!=static_cast<long long>(i)) return false;
        if (std::count_if(set->positions.begin(),set->positions.end(),[](long long p){return p>=0;})!=static_cast<long long>(set->ids.size())) return false;
    }
    return true;
}
void Simulation::writeXYZ(std::ostream& out) const {
    out<<std::setprecision(17)<<m_solidCount+(m_adsAB-m_reacted)+(2*m_adsB2-m_reacted)<<"\nLattice=\"";
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) out<<(i||j ? " " : "")<<m_solid.cell[i][j];
    out<<"\" Properties=species:S:1:pos:R:3:site_id:I:1:fixed:L:1:chemical_state:S:1:role:S:1 pbc=\"T T F\" Time="<<m_time<<'\n';
    for (const auto& site:m_solid.sites) {
        const auto ads=m_adsorbates[site.id];
        if (!site.occupied && ads==ZGB::Adsorbate::Vacant) continue;
        const auto state=site.occupied ? (site.species=="AB2" ? "AB2_solid" : "substrate") : (ads==ZGB::Adsorbate::AB ? "AB_adsorbate" : "B_adsorbate");
        const auto element=site.occupied && site.species!="AB2" ? m_solid.element : "X";
        out<<element<<' '<<site.position[0]<<' '<<site.position[1]<<' '<<site.position[2]<<' '<<site.id<<' '<<(site.fixed ? "T" : "F")<<' '<<state<<' '<<(site.occupied ? "solid" : "adsorbate")<<'\n';
    }
    if (!out) throw std::runtime_error("failed to write ZGB 3D XYZ frame");
}
void runToFiles(const Options& options) {
    Simulation simulation(options);
    std::ofstream log("zgb3d.log"),trajectory;
    if (!log) throw std::runtime_error("cannot open zgb3d.log");
    if (options.writeXYZ) {
        trajectory.open("zgb3d_trajectory.xyz");std::ofstream initial("zgb3d_initial.xyz");
        if (!initial || !trajectory) throw std::runtime_error("cannot open ZGB 3D XYZ outputs");
        simulation.writeXYZ(initial);
    }
    log<<std::setprecision(17)<<"time\ttheta_AB\ttheta_B\ttheta_vacant\tmean_height\ttime_average_height\tmean_growth_rate\tsolid_units\tAB_adsorptions\tB2_adsorptions\treactions\tevents\n";
    double lastLog=-std::numeric_limits<double>::infinity(),lastFrame=lastLog;
    simulation.run([&](const Simulation& model,bool final) {
        const auto time=model.time();
        if ((final && time!=lastLog) || time-lastLog>=options.kinetics.logInterval) {
            const auto s=model.statistics();
            log<<s.time<<'\t'<<s.coverageAB<<'\t'<<s.coverageB<<'\t'<<s.coverageVacant<<'\t'<<s.meanHeight<<'\t'<<s.timeAverageHeight<<'\t'<<s.growthRate<<'\t'<<model.solid().atomCount()<<'\t'<<s.adsorbedAB<<'\t'<<s.adsorbedB2<<'\t'<<s.reactions<<'\t'<<s.events<<'\n';lastLog=time;
        }
        if (options.writeXYZ && ((final && time!=lastFrame) || !std::isfinite(lastFrame) || (options.snapshotInterval>0 && time-lastFrame>=options.snapshotInterval))) {
            model.writeXYZ(trajectory);lastFrame=time;
        }
    });
    const auto s=simulation.statistics();const auto& k=simulation.options().kinetics;
    nlohmann::json summary={{"model","cheimarios_zgb_2023_3d"},{"representation","explicit_3d_occupancy"},{"growth",k.growth},{"y_ab",k.yAB},{"reaction_rate",k.reactionRate},{"pair_convention",k.pairConvention},{"reaction_selection",k.reactionSelection},{"repetitions",options.geometry.repetitions},{"seed",k.seed},{"time",s.time},{"start_time",k.start},{"requested_end_time",k.end},{"stop_reason",s.stopReason},{"poisoned",s.poisoned},{"capacity_limited",simulation.capacityLimited()},{"solid_units",simulation.solid().atomCount()},{"coverage_ab",s.coverageAB},{"coverage_b",s.coverageB},{"coverage_vacant",s.coverageVacant},{"mean_height",s.meanHeight},{"time_average_height",s.timeAverageHeight},{"mean_growth_rate",s.growthRate},{"rms",s.rms},{"ab_adsorptions",s.adsorbedAB},{"b2_adsorptions",s.adsorbedB2},{"reactions",s.reactions},{"events",s.events}};
    std::ofstream result("zgb3d_summary.json");result<<summary.dump(2)<<'\n';result.close();log.close();
    if (options.writeXYZ) {
        std::ofstream final("zgb3d_final.xyz");if (!final) throw std::runtime_error("cannot open zgb3d_final.xyz");simulation.writeXYZ(final);final.close();trajectory.close();
        if (!final || !trajectory) throw std::runtime_error("failed to write ZGB 3D geometry output");
    }
    if (!result || !log) throw std::runtime_error("failed to write ZGB 3D summary or log");
}
} // namespace ZGB3D
