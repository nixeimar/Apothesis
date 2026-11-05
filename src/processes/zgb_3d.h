#ifndef APOTHESIS_ZGB_3D_H
#define APOTHESIS_ZGB_3D_H
#include "zgb_growth.h"
#include "../lattice/atomistic_surface.h"
#include <iosfwd>

namespace ZGB3D {
struct Options {
    Atomistic::Options geometry;
    ZGB::Options kinetics;
    double snapshotInterval=10;
    bool writeXYZ=true;
};
void validate(const Options& options);
// Explicit volume occupancy is authoritative. No ZGB::Simulation or column
// height array is used to advance this backend.
class Simulation {
public:
    explicit Simulation(Options options);
    void apply(ZGB::Event event,std::size_t first,std::size_t second=0);
    void run(const std::function<void(const Simulation&,bool final)>& observe={});
    const Atomistic::Surface& solid() const {return m_solid;}
    const std::vector<ZGB::Adsorbate>& adsorbates() const {return m_adsorbates;}
    const Options& options() const {return m_options;}
    bool isBinding(std::size_t id) const {return m_binding.at(id);}
    double time() const {return m_time;}
    bool capacityLimited() const {return m_capacity;}
    ZGB::Statistics statistics() const;
    bool catalogConsistent() const;
    void writeXYZ(std::ostream& out) const;
private:
    struct IndexedSet {
        std::vector<std::size_t> ids;
        std::vector<long long> positions;
        explicit IndexedSet(std::size_t size=0):positions(size,-1){}
        void set(std::size_t id,bool present);
    };
    Options m_options;
    Atomistic::Surface m_solid;
    std::vector<ZGB::Adsorbate> m_adsorbates;
    std::vector<bool> m_binding;
    std::vector<std::array<std::size_t,4>> m_lateral,m_incident;
    std::vector<std::array<std::size_t,2>> m_edges;
    std::vector<std::size_t> m_up,m_down;
    IndexedSet m_vacancies,m_pairs,m_reactions,m_reactiveSites;
    std::mt19937_64 m_rng;
    double m_time=0,m_initialMean=0;
    long double m_integral=0;
    std::uint64_t m_adsAB=0,m_adsB2=0,m_reacted=0,m_steps=0;
    std::size_t m_solidCount=0;
    bool m_capacity=false;
    std::string m_stop;
    bool pairEligible(std::size_t edge) const;
    bool reactionEligible(std::size_t edge) const;
    void refresh(std::size_t id);
    double pairRate() const;
};
void runToFiles(const Options& options);
} // namespace ZGB3D
#endif
