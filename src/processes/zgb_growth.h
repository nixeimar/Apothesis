#ifndef APOTHESIS_ZGB_GROWTH_H
#define APOTHESIS_ZGB_GROWTH_H
#include <array>
#include <cstdint>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace ZGB {
enum class Adsorbate { Vacant, AB, B };
enum class Event { AdsorbAB, AdsorbB2, React };
struct Column { int height=10; Adsorbate adsorbate=Adsorbate::Vacant; };
struct Options {
    int x=100,y=100,initialHeight=10;
    double yAB=0.4,reactionRate=1e5,start=0,end=100,logInterval=1;
    bool growth=true;
    // zgb_attempts: each vacant directed neighbor attempt has rate (1-yAB)/4.
    // ordered_pairs is provided to audit the normalization ambiguity in Table 1.
    std::string pairConvention="zgb_attempts";
    std::string reactionSelection="sites";
    std::uint64_t seed=1234,maxEvents=0;
};
struct Statistics {
    double time=0,coverageAB=0,coverageB=0,coverageVacant=0,meanHeight=0,rms=0,growthRate=0,timeAverageHeight=0;
    std::uint64_t adsorbedAB=0,adsorbedB2=0,reactions=0,events=0;
    std::size_t vacantSites=0,vacantPairs=0,reactionPairs=0,reactiveSites=0;
    bool poisoned=false;
    std::string stopReason;
};
void validate(const Options& options);
class Simulation {
public:
    explicit Simulation(Options options,std::vector<Column> initial={});
    void run(const std::function<void(const Simulation&,bool final)>& observe={});
    void apply(Event event,std::size_t first,std::size_t second=0);
    const std::vector<Column>& surface() const {return m_columns;}
    const Options& options() const {return m_options;}
    double time() const {return m_time;}
    Statistics statistics() const;
    // Independent full recomputation for validation of incremental catalogs.
    bool catalogConsistent() const;
private:
    struct DenseSet {
        std::vector<std::size_t> ids;
        std::vector<long long> positions;
        explicit DenseSet(std::size_t capacity):positions(capacity,-1){}
        void set(std::size_t id,bool present);
    };
    Options m_options;
    std::vector<Column> m_columns;
    std::vector<std::array<std::size_t,4>> m_neighbors,m_incident;
    std::vector<std::array<std::size_t,2>> m_edges;
    DenseSet m_vacancies,m_pairs,m_reactions,m_reactiveSites;
    std::mt19937_64 m_rng;
    double m_time,m_initialMean;
    long double m_heightIntegral=0;
    std::uint64_t m_adsAB=0,m_adsB2=0,m_reacted=0,m_steps=0;
    std::string m_stop;
    bool adjacent(std::size_t a,std::size_t b) const;
    bool pairEligible(std::size_t edge) const;
    bool reactionEligible(std::size_t edge) const;
    void refresh(std::size_t site);
    double b2PairRate() const;
};
void runToFiles(Options options);
} // namespace ZGB
#endif
