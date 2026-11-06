#include "atomistic_surface.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace Atomistic {
namespace {
using Basis = std::vector<Vec3>;
Basis basis(const std::string& crystal) {
    if (crystal == "SimpleCubic") return {{{0,0,0}}};
    if (crystal == "BCC") return {{{0,0,0}}, {{0.5,0.5,0.5}}};
    if (crystal == "FCC") return {{{0,0,0}}, {{0,0.5,0.5}}, {{0.5,0,0.5}}, {{0.5,0.5,0}}};
    if (crystal == "HCP") return {{{0,0,0}}, {{2.0/3,1.0/3,0.5}}};
    if (crystal == "Diamond") {
        Basis b = basis("FCC");
        for (std::size_t i = 0; i < 4; ++i)
            b.push_back({{b[i][0]+0.25,b[i][1]+0.25,b[i][2]+0.25}});
        return b;
    }
    throw std::runtime_error("lattice.type: unsupported atomistic crystal '" + crystal + "'");
}
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec3 cartesian(const std::array<Vec3,3>& cell, const Vec3& f) {
    Vec3 r{{0,0,0}};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) r[j] += f[i]*cell[i][j];
    return r;
}
std::array<Vec3,3> cellFor(const Options& o) {
    const double a = o.latticeConstant;
    std::array<Vec3,3> cell{{{{a*o.repetitions[0],0,0}},{{0,a*o.repetitions[1],0}},{{0,0,a*o.repetitions[2]+o.vacuum}}}};
    if (o.crystal == "HCP") {
        cell[1] = {{-0.5*a*o.repetitions[1], std::sqrt(3.0)/2*a*o.repetitions[1],0}};
        cell[2][2] = a*o.cOverA*o.repetitions[2]+o.vacuum;
    }
    return cell;
}
double cutoffFor(const Options& o) {
    if (o.neighborCutoff > 0) return o.neighborCutoff;
    double distance = o.latticeConstant;
    if (o.crystal == "FCC") distance /= std::sqrt(2.0);
    if (o.crystal == "BCC") distance *= std::sqrt(3.0)/2;
    if (o.crystal == "Diamond") distance *= std::sqrt(3.0)/4;
    if (o.crystal == "HCP") distance *= std::max(1.0,std::sqrt(1.0/3+o.cOverA*o.cOverA/4));
    return distance*1.01;
}
// Fractional cell binning also works for the nonorthogonal HCP cell.
void connect(Surface& surface, const std::vector<Vec3>& fractions) {
    const double volume = dot(surface.cell[0],cross(surface.cell[1],surface.cell[2]));
    std::array<int,3> bins;
    for (int axis=0;axis<3;++axis) {
        auto normal = cross(surface.cell[(axis+1)%3],surface.cell[(axis+2)%3]);
        double planeHeight = volume/std::sqrt(dot(normal,normal));
        bins[axis] = std::max(1,static_cast<int>(std::min(1000000.0,std::floor(planeHeight/surface.neighborCutoff))));
    }
    std::map<std::array<int,3>,std::vector<std::size_t>> buckets;
    auto binOf = [&](const Vec3& f) {
        std::array<int,3> key;
        for (int axis=0;axis<3;++axis)
            key[axis] = std::min(bins[axis]-1, std::max(0,static_cast<int>(std::floor(f[axis]*bins[axis]))));
        return key;
    };
    for (std::size_t i=0;i<fractions.size();++i) buckets[binOf(fractions[i])].push_back(i);
    const double cutoff2 = surface.neighborCutoff*surface.neighborCutoff;
    for (std::size_t i=0;i<surface.sites.size();++i) {
        const auto key = binOf(fractions[i]);
        std::set<std::array<int,3>> candidates;
        for (int x=-1;x<=1;++x) for (int y=-1;y<=1;++y) for (int z=-1;z<=1;++z) {
            std::array<int,3> k{{key[0]+x,key[1]+y,key[2]+z}};
            bool valid = true;
            for (int axis=0;axis<3;++axis) {
                if (surface.periodic[axis]) k[axis] = (k[axis]+bins[axis])%bins[axis];
                else if (k[axis]<0 || k[axis]>=bins[axis]) valid=false;
            }
            if (valid) candidates.insert(k);
        }
        for (const auto& k : candidates) {
            auto bucket = buckets.find(k);
            if (bucket == buckets.end()) continue;
            for (auto j : bucket->second) {
                if (j<=i) continue;
                Vec3 delta;
                for (int axis=0;axis<3;++axis) {
                    delta[axis] = fractions[j][axis]-fractions[i][axis];
                    if (surface.periodic[axis]) delta[axis] -= std::round(delta[axis]);
                }
                double best = std::numeric_limits<double>::infinity();
                // Search adjacent images for skew cells rather than assuming
                // fractional rounding always yields the shortest image.
                for (int x=-1;x<=1;++x) for (int y=-1;y<=1;++y) for (int z=-1;z<=1;++z) {
                    if ((!surface.periodic[0] && x) || (!surface.periodic[1] && y) || (!surface.periodic[2] && z)) continue;
                    auto r = cartesian(surface.cell,{{delta[0]+x,delta[1]+y,delta[2]+z}});
                    best = std::min(best,dot(r,r));
                }
                if (best<=cutoff2*(1+1e-12)) {
                    surface.sites[i].neighbors.push_back(j);
                    surface.sites[j].neighbors.push_back(i);
                }
            }
        }
    }
    for (auto& site : surface.sites) std::sort(site.neighbors.begin(),site.neighbors.end());
}
} // namespace
std::size_t validate(const Options& o) {
    const auto b = basis(o.crystal);
    if (!std::isfinite(o.latticeConstant) || o.latticeConstant<=0 || !std::isfinite(o.cOverA) || o.cOverA<=0)
        throw std::runtime_error("lattice.parameters: lattice_constant and c_over_a must be positive and finite");
    if (!std::isfinite(o.vacuum) || o.vacuum<0 || !std::isfinite(o.neighborCutoff) || o.neighborCutoff<0)
        throw std::runtime_error("lattice.parameters: vacuum and neighbor_cutoff must be nonnegative and finite");
    if (o.periodic[2] && o.vacuum>0) throw std::runtime_error("lattice.parameters: vacuum requires a nonperiodic z direction");
    std::size_t n = b.size();
    for (int i=0;i<3;++i) {
        if (o.repetitions[i]<=0 || static_cast<std::size_t>(o.repetitions[i])>1000000/n)
            throw std::runtime_error("lattice.parameters.repetitions: positive dimensions required; limit is 1000000 sites");
        n *= static_cast<std::size_t>(o.repetitions[i]);
        if (!std::isfinite(o.origin[i])) throw std::runtime_error("lattice.parameters.origin: coordinates must be finite");
    }
    if (o.species.empty() || o.element.empty()) throw std::runtime_error("lattice.parameters.initial_species: atomic species and element are required");
    if (o.fixedUnitCells<0 || o.fixedUnitCells>o.repetitions[2])
        throw std::runtime_error("lattice.parameters.fixed_unit_cells: must be between zero and repetitions[2]");
    std::set<std::size_t> unique;
    for (auto id : o.vacancies) if (id>=n || !unique.insert(id).second)
        throw std::runtime_error("lattice.parameters.vacancies: site IDs must be unique and in range");
    if (o.occupiedUnitCells < -1 || o.occupiedUnitCells>o.repetitions[2] ||
        (o.occupiedUnitCells>=0 && o.fixedUnitCells>o.occupiedUnitCells))
        throw std::runtime_error("lattice.parameters.occupied_unit_cells: must cover fixed cells and fit within repetitions[2]");
    auto cell = cellFor(o);
    auto cutoff = cutoffFor(o);
    const double volume = dot(cell[0],cross(cell[1],cell[2]));
    if (!std::isfinite(volume) || volume<=0 || !std::isfinite(cutoff))
        throw std::runtime_error("lattice.parameters: cell dimensions or cutoff are out of range");
    for (int i=0;i<3;++i) {
        auto normal = cross(cell[(i+1)%3],cell[(i+2)%3]);
        double height = volume/std::sqrt(dot(normal,normal));
        if (o.periodic[i] && 2*cutoff>=height)
            throw std::runtime_error("lattice.parameters.repetitions: periodic cell heights must exceed twice the neighbor cutoff");
    }
    return n;
}
Surface build(const Options& o) {
    const auto n = validate(o);
    Surface s;
    s.cell=cellFor(o); s.periodic=o.periodic; s.element=o.element; s.neighborCutoff=cutoffFor(o);
    s.elements[o.species]=o.element;
    s.sites.reserve(n);
    std::vector<Vec3> fractions;
    fractions.reserve(n);
    const auto b = basis(o.crystal);
    const double zSize = o.latticeConstant*(o.crystal=="HCP" ? o.cOverA : 1.0)*o.repetitions[2];
    for (int z=0;z<o.repetitions[2];++z) for (int y=0;y<o.repetitions[1];++y)
        for (int x=0;x<o.repetitions[0];++x) for (const auto& f : b) {
            Vec3 fractional{{(x+f[0])/o.repetitions[0],(y+f[1])/o.repetitions[1],
                ((z+f[2])*zSize/o.repetitions[2]+o.vacuum/2)/(zSize+o.vacuum)}};
            auto r=cartesian(s.cell,fractional);
            for (int axis=0;axis<3;++axis) r[axis]+=o.origin[axis];
            s.sites.push_back({s.sites.size(),r,o.species,o.occupiedUnitCells<0 || z<o.occupiedUnitCells,z<o.fixedUnitCells,{}});
            fractions.push_back(fractional);
        }
    for (auto id : o.vacancies) s.sites[id].occupied=false;
    if (o.crystal=="SimpleCubic" && s.neighborCutoff>=o.latticeConstant*(1-1e-12) &&
        s.neighborCutoff<std::sqrt(2.0)*o.latticeConstant*(1-1e-12)) {
        // Exact first-neighbor graph avoids a spatial search over large SC slabs.
        const auto nx=o.repetitions[0],ny=o.repetitions[1],nz=o.repetitions[2];
        for (int z=0;z<nz;++z) for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
            const auto id=static_cast<std::size_t>((z*ny+y)*nx+x);
            for (int axis=0;axis<3;++axis) for (int step:{-1,1}) {
                std::array<int,3> c{{x,y,z}};c[axis]+=step;
                const int size=o.repetitions[axis];
                if (o.periodic[axis]) c[axis]=(c[axis]+size)%size;
                else if (c[axis]<0 || c[axis]>=size) continue;
                s.sites[id].neighbors.push_back(static_cast<std::size_t>((c[2]*ny+c[1])*nx+c[0]));
            }
            std::sort(s.sites[id].neighbors.begin(),s.sites[id].neighbors.end());
        }
    } else connect(s,fractions);
    return s;
}

Vec3 displacement(const Surface& s,const Vec3& from,const Vec3& to) {
    const double volume=dot(s.cell[0],cross(s.cell[1],s.cell[2]));
    if (!std::isfinite(volume) || volume<=0) throw std::runtime_error("invalid geometric cell");
    Vec3 delta{{to[0]-from[0],to[1]-from[1],to[2]-from[2]}},f;
    for (int k=0;k<3;++k) {f[k]=dot(delta,cross(s.cell[(k+1)%3],s.cell[(k+2)%3]))/volume;if(s.periodic[k]) f[k]-=std::round(f[k]);}
    Vec3 best=cartesian(s.cell,f);double distance=dot(best,best);
    for (int x=-1;x<=1;++x) for (int y=-1;y<=1;++y) for (int z=-1;z<=1;++z) {
        if ((!s.periodic[0] && x) || (!s.periodic[1] && y) || (!s.periodic[2] && z)) continue;
        auto d=cartesian(s.cell,{{f[0]+x,f[1]+y,f[2]+z}});auto r=dot(d,d);
        if (r<distance) {distance=r;best=d;}
    }
    return best;
}
std::vector<std::vector<std::size_t>> neighborGraph(const Surface& source,double minimum,double maximum,bool sameHeight) {
    if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum<0 || maximum<=minimum) throw std::runtime_error("neighbor distance interval requires 0 <= min_distance < max_distance");
    const double volume=dot(source.cell[0],cross(source.cell[1],source.cell[2]));
    if (!std::isfinite(volume) || volume<=0) throw std::runtime_error("invalid geometric cell");
    for (int k=0;k<3;++k) {
        const auto normal=cross(source.cell[(k+1)%3],source.cell[(k+2)%3]);
        const double height=volume/std::sqrt(dot(normal,normal));
        if (source.periodic[k] && 2*maximum>=height) throw std::runtime_error("periodic cell heights must exceed twice the process neighbor distance");
    }
    Surface s=source;s.neighborCutoff=maximum;
    std::vector<Vec3> fractions;fractions.reserve(s.sites.size());
    for (auto& site:s.sites) {
        site.neighbors.clear();Vec3 f;
        for (int k=0;k<3;++k) {f[k]=dot(site.position,cross(s.cell[(k+1)%3],s.cell[(k+2)%3]))/volume;if(s.periodic[k]) f[k]-=std::floor(f[k]);}
        fractions.push_back(f);
    }
    connect(s,fractions);
    std::vector<std::vector<std::size_t>> result(s.sites.size());
    for (const auto& site:s.sites) for (auto j:site.neighbors) {
        auto d=displacement(s,site.position,s.sites[j].position);const double r2=dot(d,d);
        if (r2>=minimum*minimum*(1-1e-12) && (!sameHeight || std::abs(d[2])<=maximum*1e-7)) result[site.id].push_back(j);
    }
    return result;
}
std::size_t Surface::atomCount() const {
    return static_cast<std::size_t>(std::count_if(sites.begin(),sites.end(),[](const Site& s){return s.occupied;}));
}
std::size_t Surface::occupiedCoordination(std::size_t id) const {
    std::size_t count=0;
    for (auto j : sites.at(id).neighbors) if (sites[j].occupied) ++count;
    return count;
}
void Surface::writeXYZ(const std::string& filename,double time) const {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error(filename+": cannot open XYZ output");
    writeXYZFrame(out,time);
    out.close();
    if (!out) throw std::runtime_error(filename+": failed to write XYZ output");
}
void Surface::writeXYZFrame(std::ostream& out,double time) const {
    out << std::setprecision(17) << atomCount() << "\nLattice=\"";
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) out << (i||j ? " " : "") << cell[i][j];
    out << "\" Properties=species:S:1:pos:R:3:site_id:I:1:fixed:L:1 pbc=\"";
    for (int i=0;i<3;++i) out << (i ? " " : "") << (periodic[i] ? "T" : "F");
    out << "\" Time=" << time << "\n";
    for (const auto& s : sites) if (s.occupied)
        out << (elements.count(s.species) ? elements.at(s.species) : element) << ' ' << s.position[0] << ' ' << s.position[1] << ' ' << s.position[2]
            << ' ' << s.id << ' ' << (s.fixed ? "T" : "F") << '\n';
    if (!out) throw std::runtime_error("failed to write XYZ frame");
}
} // namespace Atomistic
