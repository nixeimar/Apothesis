#ifndef APOTHESIS_ATOMISTIC_KMC_H
#define APOTHESIS_ATOMISTIC_KMC_H
#include "../lattice/atomistic_surface.h"
#include <cstdint>
#include <functional>
#include <map>
#include <random>

namespace Atomistic {
enum class Kind { Adsorption, Diffusion, Desorption };
struct Process {
    std::string id, species, element, model = "constant";
    Kind kind = Kind::Adsorption;
    double value = 0, prefactor = 0, activationEnergy = 0, bondEnergy = 0;
    int minCoordination = 1, maxCoordination = 1000000;
};
struct RunOptions {
    double start = 0, end = 0, temperature = 0, logInterval = 1, snapshotInterval = 1;
    std::uint64_t seed = 0, maxEvents = 0;
    std::vector<Process> processes;
};
struct Event {
    std::size_t process, source, target;
    double rate;
};
class Simulation {
public:
    Simulation(Surface surface, RunOptions options);
    std::vector<Event> events() const;
    // Revalidates the event; stale or illegal transitions are rejected.
    void apply(const Event& event);
    void run(const std::function<void(const Simulation&, bool final)>& observe = {});
    const Surface& surface() const { return m_surface; }
    const RunOptions& options() const { return m_options; }
    const std::vector<std::uint64_t>& counts() const { return m_counts; }
    double time() const { return m_time; }
    std::uint64_t steps() const { return m_steps; }
    const std::string& stopReason() const { return m_stopReason; }
private:
    Surface m_surface;
    RunOptions m_options;
    std::vector<std::uint64_t> m_counts;
    std::mt19937_64 m_rng;
    double m_time;
    std::uint64_t m_steps = 0;
    std::string m_stopReason;
    std::vector<bool> accessible() const;
    void execute(const Event& event);
    double rate(const Process& process, std::size_t coordination) const;
};
void validate(const RunOptions& options);
void runToFiles(Surface surface, RunOptions options);
} // namespace Atomistic
#endif
