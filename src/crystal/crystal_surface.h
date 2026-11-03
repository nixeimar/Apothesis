#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace Crystal {
using SiteId = std::size_t;
struct Vec3 { double x, y, z; };
enum class Face { FCC111, FCC100, HCP0001 };
enum class SiteType { Top, Bridge, FccHollow, HcpHollow, Fourfold };
const char* name(SiteType type);
struct Definition {
    Face face = Face::FCC111;
    double a = 1.0; // conventional cubic a for FCC; basal a for HCP
    double c = 1.632993161855452; // used only for HCP
    int nx = 4, ny = 4, layers = 3;
};
struct Atom { Vec3 position; int layer; };
struct AtomImage { std::size_t atom; int sx, sy; };
struct Site {
    Vec3 position; // geometric centre in top-layer plane; not relaxed adsorbate z
    SiteType type;
    std::vector<AtomImage> support;
};
struct Edge {
    SiteId from, to;
    int sx, sy; // image of destination relative to source
    double distance;
};
struct Change {
    // Conservative invalidation set: all sites. Safe for arbitrary rate dependencies.
    std::vector<SiteId> affected;
};
class CrystalSurface {
public:
    explicit CrystalSurface(const Definition&);
    const std::vector<Atom>& atoms() const { return atoms_; }
    const std::vector<Site>& sites() const { return sites_; }
    const std::vector<Edge>& transitions() const { return edges_; }
    Vec3 cellA() const { return cellA_; }
    Vec3 cellB() const { return cellB_; }
    double spacing() const { return spacing_; }

    // No site is occupiable until explicitly enabled by a physical model.
    void enable(SiteType type);
    bool enabled(SiteId id) const;
    // User-selected hard-core distance, not an inferred chemical footprint.
    // Must be configured before adsorbates are present.
    void setExclusionDistance(double distance);
    const std::vector<SiteId>& exclusions(SiteId id) const;
    // Adds directed transitions and their reverse, retaining periodic images.
    // Choose pairs/cutoff from the physical transition catalogue.
    void connect(SiteType from, SiteType to, double cutoff);
    bool canAdsorb(SiteId id) const;
    bool canHop(std::size_t edge) const;
    Change adsorb(SiteId id, const std::string& species);
    Change desorb(SiteId id);
    Change hop(std::size_t edge);
    const std::string& occupant(SiteId id) const;
    void writeCSV(const std::string& prefix) const;
private:
    Definition d_;
    double spacing_;
    Vec3 cellA_, cellB_;
    std::vector<Atom> atoms_;
    std::vector<Site> sites_;
    std::vector<Edge> edges_;
    std::vector<bool> enabled_;
    std::vector<std::string> occupancy_;
    std::vector<std::vector<SiteId>> exclusions_;
    bool available(SiteId id, SiteId ignore) const;
    Change changed() const;
};
} // namespace Crystal
