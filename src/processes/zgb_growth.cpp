#include "zgb_growth.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include "../extLibs/nlohmann/json.hpp"

namespace ZGB {
void validate(const Options& o) {
    if (o.x<3 || o.y<3 || static_cast<long long>(o.x)*o.y>1000000)
        throw std::runtime_error("lattice.parameters.dimensions: square-grid dimensions must be at least 3, with at most 1000000 columns");
    if (o.initialHeight<0) throw std::runtime_error("lattice.parameters.initial_height: must be nonnegative");
    if (!std::isfinite(o.yAB) || o.yAB<0 || o.yAB>1) throw std::runtime_error("model.parameters.y_ab: must be between 0 and 1");
    if (!std::isfinite(o.reactionRate) || o.reactionRate<0) throw std::runtime_error("model.parameters.reaction_rate: must be finite and nonnegative");
    if (!std::isfinite(o.start) || o.start<0 || !std::isfinite(o.end) || o.end<=o.start || !std::isfinite(o.logInterval) || o.logInterval<=0)
        throw std::runtime_error("simulation: duration and output.log_interval must be positive and finite");
    if (o.reactionSelection!="sites" && o.reactionSelection!="pairs")
        throw std::runtime_error("model.parameters.reaction_selection: expected sites or pairs");
    if (o.pairConvention!="zgb_attempts" && o.pairConvention!="ordered_pairs")
        throw std::runtime_error("model.parameters.pair_convention: expected zgb_attempts or ordered_pairs");
}
void Simulation::DenseSet::set(std::size_t id,bool present) {
    auto& p=positions.at(id);
    if (present && p<0) {p=static_cast<long long>(ids.size());ids.push_back(id);}
    else if (!present && p>=0) {
        auto replacement=ids.back();ids[static_cast<std::size_t>(p)]=replacement;
        positions[replacement]=p;ids.pop_back();p=-1;
    }
}
Simulation::Simulation(Options options,std::vector<Column> initial)
    : m_options(options),m_vacancies(0),m_pairs(0),m_reactions(0),m_reactiveSites(0),m_time(options.start),m_initialMean(0) {
    validate(m_options);
    const auto n=static_cast<std::size_t>(m_options.x)*m_options.y;
    m_columns=initial.empty() ? std::vector<Column>(n,{m_options.initialHeight,Adsorbate::Vacant}) : std::move(initial);
    if (m_columns.size()!=n) throw std::runtime_error("ZGB initial surface size does not match dimensions");
    long double initialSum=0;
    for (const auto& c:m_columns) {
        if (c.height<0 || (c.adsorbate!=Adsorbate::Vacant && c.adsorbate!=Adsorbate::AB && c.adsorbate!=Adsorbate::B))
            throw std::runtime_error("invalid ZGB initial column");
        initialSum+=c.height;
    }
    m_initialMean=static_cast<double>(initialSum/n);
    m_vacancies=DenseSet(n);m_pairs=DenseSet(2*n);m_reactions=DenseSet(2*n);m_reactiveSites=DenseSet(n);
    m_neighbors.resize(n);m_incident.resize(n);m_edges.resize(2*n);
    for (int y=0;y<m_options.y;++y) for (int x=0;x<m_options.x;++x) {
        const auto i=static_cast<std::size_t>(y*m_options.x+x);
        const auto east=static_cast<std::size_t>(y*m_options.x+(x+1)%m_options.x);
        const auto west=static_cast<std::size_t>(y*m_options.x+(x+m_options.x-1)%m_options.x);
        const auto south=static_cast<std::size_t>(((y+1)%m_options.y)*m_options.x+x);
        const auto north=static_cast<std::size_t>(((y+m_options.y-1)%m_options.y)*m_options.x+x);
        m_neighbors[i]={{east,west,south,north}};
        m_edges[2*i]={{i,east}};m_edges[2*i+1]={{i,south}};
        m_incident[i]={{2*i,2*west,2*i+1,2*north+1}};
    }
    for (std::size_t i=0;i<n;++i) refresh(i);
    if (!m_options.seed) {std::random_device rd;m_options.seed=(static_cast<std::uint64_t>(rd())<<32)^rd();}
    m_rng.seed(m_options.seed);
}
bool Simulation::adjacent(std::size_t a,std::size_t b) const {
    const auto& neighbors=m_neighbors.at(a);
    return std::find(neighbors.begin(),neighbors.end(),b)!=neighbors.end();
}
bool Simulation::pairEligible(std::size_t edge) const {
    const auto& e=m_edges[edge];const auto& a=m_columns[e[0]];const auto& b=m_columns[e[1]];
    return a.adsorbate==Adsorbate::Vacant && b.adsorbate==Adsorbate::Vacant && (!m_options.growth || a.height==b.height);
}
bool Simulation::reactionEligible(std::size_t edge) const {
    const auto& e=m_edges[edge];const auto& a=m_columns[e[0]];const auto& b=m_columns[e[1]];
    return ((a.adsorbate==Adsorbate::AB && b.adsorbate==Adsorbate::B) ||
            (b.adsorbate==Adsorbate::AB && a.adsorbate==Adsorbate::B)) && (!m_options.growth || a.height==b.height);
}
void Simulation::refresh(std::size_t site) {
    m_vacancies.set(site,m_columns[site].adsorbate==Adsorbate::Vacant);
    for (auto e:m_incident[site]) {m_pairs.set(e,pairEligible(e));m_reactions.set(e,reactionEligible(e));}
    auto updateReactiveSite=[&](std::size_t id) {
        bool active=false;
        for (auto e:m_incident[id]) if (reactionEligible(e)) {active=true;break;}
        m_reactiveSites.set(id,active);
    };
    updateReactiveSite(site);
    for (auto j:m_neighbors[site]) updateReactiveSite(j);
}
double Simulation::b2PairRate() const {
    return (1-m_options.yAB)*(m_options.pairConvention=="zgb_attempts" ? 0.5 : 2.0);
}
void Simulation::apply(Event event,std::size_t first,std::size_t second) {
    auto& a=m_columns.at(first);
    if (event==Event::AdsorbAB) {
        if (a.adsorbate!=Adsorbate::Vacant) throw std::runtime_error("AB adsorption requires a vacant site");
        a.adsorbate=Adsorbate::AB;++m_adsAB;
    } else {
        auto& b=m_columns.at(second);
        if (!adjacent(first,second)) throw std::runtime_error("ZGB pair must be nearest neighbors");
        if (m_options.growth && a.height!=b.height) throw std::runtime_error("growth-model pairs must have the same height");
        if (event==Event::AdsorbB2) {
            if (a.adsorbate!=Adsorbate::Vacant || b.adsorbate!=Adsorbate::Vacant) throw std::runtime_error("B2 adsorption requires two vacant sites");
            a.adsorbate=b.adsorbate=Adsorbate::B;++m_adsB2;
        } else {
            if (!((a.adsorbate==Adsorbate::AB && b.adsorbate==Adsorbate::B) || (b.adsorbate==Adsorbate::AB && a.adsorbate==Adsorbate::B)))
                throw std::runtime_error("reaction requires AB* and B*");
            auto& precursor=a.adsorbate==Adsorbate::AB ? a : b;
            if (m_options.growth) {
                if (precursor.height==std::numeric_limits<int>::max()) throw std::runtime_error("ZGB height overflow");
                ++precursor.height;
            }
            a.adsorbate=b.adsorbate=Adsorbate::Vacant;++m_reacted;
        }
        refresh(second);
    }
    refresh(first);++m_steps;
}
void Simulation::run(const std::function<void(const Simulation&,bool)>& observe) {
    if (observe) observe(*this,false);
    while (m_time<m_options.end) {
        const double a=m_options.yAB*m_vacancies.ids.size(),b=b2PairRate()*m_pairs.ids.size(),r=m_options.reactionRate*(m_options.reactionSelection=="sites" ? m_reactiveSites.ids.size() : m_reactions.ids.size());
        const double total=a+b+r;
        if (!std::isfinite(total)) throw std::runtime_error("ZGB total rate overflow");
        if (total==0) {m_stop="absorbing";break;}
        if (m_options.maxEvents && m_steps>=m_options.maxEvents) {m_stop="event_limit";break;}
        double u=0;do {u=std::generate_canonical<double,53>(m_rng);} while (!u);
        const double dt=-std::log(u)/total;
        const double currentMean=m_initialMean+(m_options.growth ? m_reacted/static_cast<double>(m_columns.size()) : 0);
        if (dt>m_options.end-m_time) {
            m_heightIntegral+=currentMean*(m_options.end-m_time);m_time=m_options.end;break;
        }
        if (!(dt>0) || m_time+dt<=m_time) throw std::runtime_error("ZGB time step below floating-point resolution");
        const double pick=std::generate_canonical<double,53>(m_rng)*total;
        const auto choose=[&](const DenseSet& set) {return set.ids[std::uniform_int_distribution<std::size_t>(0,set.ids.size()-1)(m_rng)];};
        m_heightIntegral+=currentMean*dt;
        m_time+=dt;
        if (pick<a || (b==0 && r==0)) apply(Event::AdsorbAB,choose(m_vacancies));
        else if (pick<a+b || r==0) {auto e=m_edges[choose(m_pairs)];apply(Event::AdsorbB2,e[0],e[1]);}
        else if (m_options.reactionSelection=="pairs") {auto e=m_edges[choose(m_reactions)];apply(Event::React,e[0],e[1]);}
        else {
            const auto center=choose(m_reactiveSites);
            std::array<std::size_t,4> partners{};std::size_t count=0;
            for (auto edge:m_incident[center]) if (reactionEligible(edge)) {
                const auto& e=m_edges[edge];partners[count++]=e[0]==center ? e[1] : e[0];
            }
            const auto partner=partners[std::uniform_int_distribution<std::size_t>(0,count-1)(m_rng)];
            apply(Event::React,center,partner);
        }
        if (observe) observe(*this,false);
    }
    if (m_stop.empty()) m_stop="end_time";
    if (observe) observe(*this,true);
}
Statistics Simulation::statistics() const {
    Statistics s;s.time=m_time;s.adsorbedAB=m_adsAB;s.adsorbedB2=m_adsB2;s.reactions=m_reacted;s.events=m_steps;
    const double n=static_cast<double>(m_columns.size());
    long double sumHeight=0;
    std::size_t countAB=0,countB=0,countVacant=0;
    for (const auto& c:m_columns) {
        sumHeight+=c.height;
        if (c.adsorbate==Adsorbate::AB) ++countAB;
        else if (c.adsorbate==Adsorbate::B) ++countB;
        else ++countVacant;
    }
    s.meanHeight=static_cast<double>(sumHeight/n);
    s.coverageAB=countAB/n;s.coverageB=countB/n;s.coverageVacant=countVacant/n;
    for (const auto& c:m_columns) s.rms+=(c.height-s.meanHeight)*(c.height-s.meanHeight)/n;
    s.rms=std::sqrt(s.rms);
    s.timeAverageHeight=m_time>m_options.start ? static_cast<double>(m_heightIntegral/(m_time-m_options.start)) : m_initialMean;
    if (m_time>m_options.start) s.growthRate=(s.meanHeight-m_initialMean)/(m_time-m_options.start);
    s.vacantSites=m_vacancies.ids.size();s.vacantPairs=m_pairs.ids.size();s.reactionPairs=m_reactions.ids.size();s.reactiveSites=m_reactiveSites.ids.size();
    s.poisoned=m_options.yAB*s.vacantSites+b2PairRate()*s.vacantPairs+m_options.reactionRate*(m_options.reactionSelection=="sites" ? s.reactiveSites : s.reactionPairs)==0;
    s.stopReason=m_stop;
    return s;
}
bool Simulation::catalogConsistent() const {
    for (std::size_t i=0;i<m_columns.size();++i)
        if ((m_vacancies.positions[i]>=0)!=(m_columns[i].adsorbate==Adsorbate::Vacant)) return false;
    for (std::size_t e=0;e<m_edges.size();++e)
        if ((m_pairs.positions[e]>=0)!=pairEligible(e) || (m_reactions.positions[e]>=0)!=reactionEligible(e)) return false;
    for (std::size_t i=0;i<m_columns.size();++i) {
        bool active=false;for (auto e:m_incident[i]) if (reactionEligible(e)) active=true;
        if ((m_reactiveSites.positions[i]>=0)!=active) return false;
    }
    for (const auto* set:{&m_vacancies,&m_pairs,&m_reactions,&m_reactiveSites}) {
        for (std::size_t i=0;i<set->ids.size();++i) if (set->positions[set->ids[i]]!=static_cast<long long>(i)) return false;
        if (std::count_if(set->positions.begin(),set->positions.end(),[](long long p){return p>=0;})!=static_cast<long long>(set->ids.size())) return false;
    }
    return true;
}
void runToFiles(Options options) {
    Simulation simulation(options);
    std::ofstream log("zgb.log");if (!log) throw std::runtime_error("cannot open zgb.log");
    log<<std::setprecision(17)<<"time\ttheta_AB\ttheta_B\ttheta_vacant\tmean_height\trms\tmean_growth_rate\ttime_average_height\tAB_adsorptions\tB2_adsorptions\treactions\tevents\n";
    double last=-std::numeric_limits<double>::infinity();
    simulation.run([&](const Simulation& model,bool final) {
        // Avoid scanning the whole lattice after every event.
        const auto time=model.time();
        if ((final && time!=last) || time-last>=options.logInterval) {
            auto s=model.statistics();
            log<<s.time<<'\t'<<s.coverageAB<<'\t'<<s.coverageB<<'\t'<<s.coverageVacant<<'\t'<<s.meanHeight<<'\t'<<s.rms<<'\t'<<s.growthRate<<'\t'<<s.timeAverageHeight<<'\t'<<s.adsorbedAB<<'\t'<<s.adsorbedB2<<'\t'<<s.reactions<<'\t'<<s.events<<'\n';last=time;
        }
    });
    const auto s=simulation.statistics();
    nlohmann::json summary={{"model","cheimarios_zgb_2023"},{"growth",options.growth},{"y_ab",options.yAB},{"reaction_rate",options.reactionRate},{"pair_convention",options.pairConvention},{"reaction_selection",options.reactionSelection},{"dimensions",{options.x,options.y}},{"initial_height",options.initialHeight},{"seed",simulation.options().seed},{"start_time",options.start},{"requested_end_time",options.end},{"time",s.time},{"stop_reason",s.stopReason},{"poisoned",s.poisoned},{"coverage_ab",s.coverageAB},{"coverage_b",s.coverageB},{"coverage_vacant",s.coverageVacant},{"mean_height",s.meanHeight},{"rms",s.rms},{"mean_growth_rate",s.growthRate},{"time_average_height",s.timeAverageHeight},{"ab_adsorptions",s.adsorbedAB},{"b2_adsorptions",s.adsorbedB2},{"reactions",s.reactions},{"events",s.events}};
    std::ofstream result("zgb_summary.json"),surface("zgb_final.csv");
    if (!result || !surface) throw std::runtime_error("cannot open ZGB summary or surface output");
    result<<summary.dump(2)<<'\n';surface<<"x,y,height,adsorbate\n";
    const auto& columns=simulation.surface();
    for (int y=0;y<options.y;++y) for (int x=0;x<options.x;++x) {
        const auto& c=columns[static_cast<std::size_t>(y*options.x+x)];
        surface<<x<<','<<y<<','<<c.height<<','<<(c.adsorbate==Adsorbate::AB ? "AB" : c.adsorbate==Adsorbate::B ? "B" : "vacant")<<'\n';
    }
    log.close();result.close();surface.close();
    if (!log || !result || !surface) throw std::runtime_error("failed to write ZGB output");
}
} // namespace ZGB
