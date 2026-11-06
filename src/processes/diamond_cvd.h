#ifndef APOTHESIS_DIAMOND_CVD_H
#define APOTHESIS_DIAMOND_CVD_H
#include "../lattice/atomistic_surface.h"
#include <random>
#include <map>
#include <functional>

namespace DiamondCVD {
struct Process {
    std::string id,event,domain;
    double rate=0;
    int hydrogens=-1;
};
struct Options {
    Atomistic::Surface lattice;
    std::vector<Process> processes;
    double temperature=1200,gasTemperature=1200,start=0,end=1,siteVolume=0,logInterval=.1,snapshotInterval=.5,chLength=1.09;
    std::uint64_t seed=1234,maxEvents=200000;
    bool passivate=true,writeXYZ=true;
    std::map<std::string,double> pressuresTorr;
};
double rateCoefficient(double A,double n,double E,double H,double S,double correction,double temperature,bool reverse);
void validate(const Options& options);
enum class Kind {Empty,Diamond,Hydrocarbon,Hydrogen};
struct Atom {Kind kind=Kind::Empty;std::size_t parent=0;bool trapped=false;};
struct Dimer {std::size_t first,second,shared;};
struct Event {
    std::size_t process=0;double rate=0;
    std::vector<std::size_t> erase,carbons,diamond;
    std::vector<std::pair<std::size_t,std::size_t>> hydrogens,bonds;
    std::vector<Dimer> close,open;
};
class Simulation {
public:
    explicit Simulation(Options options);
    const Options& options() const {return m_options;}
    const std::vector<Atom>& atoms() const {return m_atoms;}
    const std::vector<std::vector<std::size_t>>& bonds() const {return m_bonds;}
    const std::vector<Dimer>& dimers() const {return m_dimers;}
    const std::vector<Event>& events() const {return m_events;}
    const std::vector<std::uint64_t>& counts() const {return m_counts;}
    double time() const {return m_time;}
    std::uint64_t steps() const {return m_steps;}
    const std::string& stopReason() const {return m_stop;}
    std::size_t count(Kind kind) const;
    std::size_t trappedHydrogens() const;
    std::size_t trappedVacancies() const;
    double meanHeight() const;
    double growthRate() const;
    double timeAverageHeight() const;
    bool consistent() const;
    void execute(std::size_t event);
    void run(const std::function<void(const Simulation&,bool)>& observe={});
    void writeXYZ(std::ostream& out) const;
private:
    Options m_options;
    std::vector<Atom> m_atoms;
    std::vector<std::vector<std::size_t>> m_bonds;
    std::vector<Dimer> m_dimers;
    std::vector<Event> m_events;
    std::vector<bool> m_gas,m_frozenVacancy;
    std::vector<std::size_t> m_surface;
    std::vector<std::uint64_t> m_counts;
    std::mt19937_64 m_rng;
    std::uint64_t m_steps=0;
    std::size_t m_diamondCount=0;
    double m_time=0,m_initialHeight=0,m_initialTop=0,m_top=0,m_area=0;
    long double m_integral=0,m_preciseTime=0,m_timeCorrection=0;
    bool m_capacity=false;
    std::string m_stop;
    bool carbon(std::size_t id) const;
    bool free(std::size_t carbon,std::size_t target) const;
    bool exposed(std::size_t id) const;
    std::size_t hydrogenSlot(std::size_t parent,std::size_t target) const;
    std::vector<std::size_t> hydrogens(std::size_t carbon) const;
    void gasAccess(bool freeze);
    void rebuild();
    void apply(const Event& event);
    bool sixRing(std::size_t first,std::size_t second) const;
};
void runToFiles(const Options& options);
}
#endif
