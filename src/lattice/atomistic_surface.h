#ifndef APOTHESIS_ATOMISTIC_SURFACE_H
#define APOTHESIS_ATOMISTIC_SURFACE_H
#include <array>
#include <cstddef>
#include <string>
#include <map>
#include <iosfwd>
#include <vector>

namespace Atomistic {
using Vec3 = std::array<double, 3>;
struct Options {
    std::string crystal, species, element;
    double latticeConstant = 0, cOverA = 1.632993161855452;
    std::array<int, 3> repetitions{{1, 1, 1}};
    std::array<bool, 3> periodic{{true, true, false}};
    Vec3 origin{{0, 0, 0}};
    double vacuum = 0, neighborCutoff = 0;
    int fixedUnitCells = 0;
    int occupiedUnitCells = -1; // -1 fills every generated unit cell.
    std::vector<std::size_t> vacancies;
};
struct Site {
    std::size_t id;
    Vec3 position;
    std::string species;
    bool occupied = true, fixed = false;
    // Geometrical candidate sites, including vacant sites; IDs remain stable.
    std::vector<std::size_t> neighbors;
};
struct Surface {
    std::array<Vec3, 3> cell;
    std::array<bool, 3> periodic;
    std::vector<Site> sites;
    std::string element;
    std::map<std::string,std::string> elements;
    double neighborCutoff;
    std::size_t atomCount() const;
    std::size_t occupiedCoordination(std::size_t id) const;
    void writeXYZ(const std::string& filename, double time = 0) const;
    void writeXYZFrame(std::ostream& out, double time) const;
};
std::size_t validate(const Options& options);
Surface build(const Options& options);
// Geometric distance graphs for process-specific shells. Distances use the
// periodic cell, including the skew HCP basis, independently of occupancy.
Vec3 displacement(const Surface& surface,const Vec3& from,const Vec3& to);
std::vector<std::vector<std::size_t>> neighborGraph(const Surface& surface,
    double minDistance,double maxDistance,bool sameHeight=false);
} // namespace Atomistic
#endif
