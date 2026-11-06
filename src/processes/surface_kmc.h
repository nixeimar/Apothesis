#ifndef APOTHESIS_SURFACE_KMC_H
#define APOTHESIS_SURFACE_KMC_H
#include "../lattice/atomistic_surface.h"
#include <functional>
#include <random>
#include <map>
#include <memory>
#include <optional>
#include <cmath>

namespace SurfaceKMC {
// Empty species denotes a vacancy. Solid products update mapped crystal nodes;
// exposure follows either the general geometry or the selected column model.
struct Product { std::string species; bool solid=false; };
struct Link {
    std::size_t first=0,second=1;
    std::string neighbors="first";
    double minDistance=0,maxDistance=0;
    bool sameHeight=false;
};
struct Rule {
    std::string id,type,counting="sites";
    std::vector<std::string> before;
    std::vector<Product> after;
    std::vector<bool> beforeSolid;
    std::vector<Link> links;
    double rate=0;
};
struct Options {
    Atomistic::Options geometry;
    std::optional<Atomistic::Surface> customGeometry;
    std::string surfaceSites="column_tops";
    double siteVolume=0;
    std::size_t maxEmbeddings=100000,maxSearchSteps=1000000;
    std::map<std::string,std::string> elements;
    std::vector<Rule> processes;
    double start=0,end=100,logInterval=1,snapshotInterval=10;
    std::uint64_t seed=1234,maxEvents=0;
    bool writeXYZ=true;
};
void validate(const Options& options);
class Simulation {
public:
    explicit Simulation(Options options);
    void run(const std::function<void(const Simulation&,bool)>& observe={});
    void apply(std::size_t process,std::size_t first,std::size_t second=0);
    void apply(std::size_t process,const std::vector<std::size_t>& sites);
    bool graphMode() const {return m_options.surfaceSites=="gas_exposed";}
    double siteVolume() const {return graphMode() ? m_siteVolume : std::pow(m_options.geometry.latticeConstant,3);}
    std::string heightUnit() const {return graphMode() ? "angstrom" : "layers";}
    const Atomistic::Surface& solid() const {return m_solid;}
    const std::vector<std::string>& adsorbates() const {return m_ads;}
    const Options& options() const {return m_options;}
    const std::vector<std::uint64_t>& counts() const {return m_counts;}
    std::size_t eligibleCount(std::size_t process) const;
    bool catalogConsistent() const;
    double time() const {return m_time;}
    double meanHeight() const;
    double timeAverageHeight() const;
    double growthRate() const;
    bool capacityLimited() const {return m_capacity;}
    std::uint64_t steps() const {return m_steps;}
    const std::string& stopReason() const {return m_stop;}
    void writeXYZ(std::ostream& out) const;
private:
    struct Catalog {
        std::vector<std::size_t> ids;
        std::vector<long long> positions;
        explicit Catalog(std::size_t n):positions(n,-1){}
        void set(std::size_t id,bool present);
    };
    Options m_options;
    Atomistic::Surface m_solid;
    std::vector<std::string> m_ads;
    std::vector<bool> m_binding,m_deposited;
    std::vector<std::array<std::size_t,4>> m_neighbors,m_incident;
    std::vector<std::array<std::size_t,2>> m_edges;
    std::vector<std::size_t> m_up,m_down;
    std::vector<Catalog> m_catalogs;
    std::vector<std::uint64_t> m_counts;
    std::mt19937_64 m_rng;
    std::size_t m_solidCount=0;
    std::uint64_t m_steps=0;
    double m_time=0,m_initialMean=0;
    long double m_integral=0;
    bool m_capacity=false;
    std::string m_stop;
    using Graph=std::vector<std::vector<std::size_t>>;
    struct EventGroup {std::vector<std::vector<std::size_t>> assignments;};
    std::vector<std::vector<std::shared_ptr<const Graph>>> m_relations;
    std::vector<std::vector<EventGroup>> m_graphCatalogs;
    std::vector<bool> m_exposed;
    double m_area=0,m_siteVolume=0,m_top=0;
    void initializeGraph();
    void graphAccess(const std::vector<bool>& occupied,std::vector<bool>& binding,std::vector<bool>& exposed) const;
    bool graphSupported(const std::vector<bool>& occupied) const;
    bool graphSlot(const Rule& rule,std::size_t slot,std::size_t id) const;
    bool graphTransition(const Rule& rule,const std::vector<std::size_t>& sites) const;
    std::vector<std::vector<EventGroup>> enumerateGraph() const;
    void rebuildGraph();
    void runGraph(const std::function<void(const Simulation&,bool)>& observe);
    bool matches(const Rule& rule,std::size_t first,std::size_t second) const;
    bool eligible(const Rule& rule,std::size_t candidate) const;
    void refresh(std::size_t site);
};
void runToFiles(const Options& options);
}
#endif
