#include "crystal_surface.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace Crystal {
namespace {
double norm(Vec3 v) { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
Vec3 delta(Vec3 a, Vec3 b, Vec3 u, Vec3 v, int x, int y) {
    return {b.x+x*u.x+y*v.x-a.x,b.y+x*u.y+y*v.y-a.y,b.z-a.z};
}
bool positive(double x) { return std::isfinite(x) && x>0; }
}
const char* name(SiteType t) {
    switch(t) {
    case SiteType::Top:return "top";
    case SiteType::Bridge:return "bridge";
    case SiteType::FccHollow:return "fcc_hollow";
    case SiteType::HcpHollow:return "hcp_hollow";
    case SiteType::Fourfold:return "fourfold";
    }
    throw std::invalid_argument("unknown site type");
}
CrystalSurface::CrystalSurface(const Definition& d):d_(d) {
    if(d.face!=Face::FCC111 && d.face!=Face::FCC100 && d.face!=Face::HCP0001)
        throw std::invalid_argument("unsupported crystal face");
    if(!positive(d.a) || (d.face==Face::HCP0001 && !positive(d.c)) ||
       d.nx<1 || d.ny<1 || d.layers<3)
        throw std::invalid_argument("positive lattice parameters/cell sizes and >=3 layers required");
    const std::size_t n=static_cast<std::size_t>(d.nx)*d.ny;
    if(n>1000000 || static_cast<std::size_t>(d.layers)>1000000/n)
        throw std::invalid_argument("slab exceeds one million atoms");
    const bool square=d.face==Face::FCC100;
    spacing_=d.face==Face::HCP0001?d.a:d.a/std::sqrt(2.0);
    Vec3 u{spacing_,0,0}, v{square?0:spacing_/2,
                                  square?spacing_:spacing_*std::sqrt(3.0)/2,0};
    cellA_={u.x*d.nx,u.y*d.nx,0};
    cellB_={v.x*d.ny,v.y*d.ny,0};
    const double dz=d.face==Face::HCP0001?d.c/2:
                    (square?d.a/2:d.a/std::sqrt(3.0));
    for(int l=0;l<d.layers;++l) {
        const double shift=square?(l%2)*0.5:
                           (d.face==Face::FCC111?(l%3)/3.0:(l%2)/3.0);
        for(int j=0;j<d.ny;++j) for(int i=0;i<d.nx;++i)
            atoms_.push_back({{(i+shift)*u.x+(j+shift)*v.x,
                               (i+shift)*u.y+(j+shift)*v.y,-l*dz},l});
    }
    struct Motif { double x,y; SiteType type; };
    std::vector<Motif> motifs={{0,0,SiteType::Top},{0.5,0,SiteType::Bridge},
                              {0,0.5,SiteType::Bridge}};
    if(square) motifs.push_back({0.5,0.5,SiteType::Fourfold});
    else {
        motifs.push_back({0.5,0.5,SiteType::Bridge});
        // Second substrate layer is at (1/3,1/3): hcp hollow.
        motifs.push_back({1.0/3,1.0/3,SiteType::HcpHollow});
        motifs.push_back({2.0/3,2.0/3,SiteType::FccHollow});
    }
    for(int j=0;j<d.ny;++j) for(int i=0;i<d.nx;++i) for(auto m:motifs) {
        Site s{{(i+m.x)*u.x+(j+m.y)*v.x,(i+m.x)*u.y+(j+m.y)*v.y,0},m.type,{}};
        // Enumerate local top-layer atoms; preserve image shifts even in a 1x1 cell.
        double best=std::numeric_limits<double>::infinity();
        for(int y=j-1;y<=j+2;++y) for(int x=i-1;x<=i+2;++x) {
            Vec3 p{x*u.x+y*v.x,x*u.y+y*v.y,0};
            best=std::min(best,norm({p.x-s.position.x,p.y-s.position.y,0}));
        }
        for(int y=j-1;y<=j+2;++y) for(int x=i-1;x<=i+2;++x) {
            Vec3 p{x*u.x+y*v.x,x*u.y+y*v.y,0};
            if(std::abs(norm({p.x-s.position.x,p.y-s.position.y,0})-best)<1e-8*spacing_) {
                int wx=(x%d.nx+d.nx)%d.nx, wy=(y%d.ny+d.ny)%d.ny;
                s.support.push_back({static_cast<std::size_t>(wy*d.nx+wx),
                                     (x-wx)/d.nx,(y-wy)/d.ny});
            }
        }
        sites_.push_back(s);
    }
    enabled_.assign(sites_.size(),false);
    occupancy_.resize(sites_.size());
    exclusions_.resize(sites_.size());
}
void CrystalSurface::enable(SiteType t) {
    for(SiteId i=0;i<sites_.size();++i) if(sites_[i].type==t) enabled_[i]=true;
}
bool CrystalSurface::enabled(SiteId i) const { return enabled_.at(i); }
const std::string& CrystalSurface::occupant(SiteId i) const { return occupancy_.at(i); }
const std::vector<SiteId>& CrystalSurface::exclusions(SiteId i) const { return exclusions_.at(i); }
void CrystalSurface::setExclusionDistance(double r) {
    if(!std::isfinite(r) || r<0 || r>=0.5*std::min(norm(cellA_),norm(cellB_)))
        throw std::invalid_argument("exclusion distance must be >=0 and < half shortest cell vector");
    for(const auto& s:occupancy_) if(!s.empty()) throw std::logic_error("configure exclusions before occupancy");
    for(auto& e:exclusions_) e.clear();
    // O(N^2), intended as a reference implementation for modest surface cells.
    for(SiteId i=0;i<sites_.size();++i) for(SiteId j=i+1;j<sites_.size();++j) {
        double best=std::numeric_limits<double>::infinity();
        for(int x=-1;x<=1;++x) for(int y=-1;y<=1;++y)
            best=std::min(best,norm(delta(sites_[i].position,sites_[j].position,cellA_,cellB_,x,y)));
        if(best<r) { exclusions_[i].push_back(j); exclusions_[j].push_back(i); }
    }
}
void CrystalSurface::connect(SiteType a, SiteType b, double cutoff) {
    if(!positive(cutoff) || cutoff>spacing_*(1+1e-10))
        throw std::invalid_argument("transition cutoff must be positive and <= nearest substrate spacing");
    for(SiteId i=0;i<sites_.size();++i) if(sites_[i].type==a)
    for(SiteId j=0;j<sites_.size();++j) if(sites_[j].type==b)
    for(int x=-1;x<=1;++x) for(int y=-1;y<=1;++y) {
        double dist=norm(delta(sites_[i].position,sites_[j].position,cellA_,cellB_,x,y));
        if(dist<1e-9*spacing_ || dist>cutoff*(1+1e-10)) continue;
        auto add=[&](SiteId from, SiteId to,int sx,int sy) {
            for(const auto& e:edges_) if(e.from==from && e.to==to && e.sx==sx && e.sy==sy) return;
            edges_.push_back({from,to,sx,sy,dist});
        };
        add(i,j,x,y); add(j,i,-x,-y);
    }
}
bool CrystalSurface::available(SiteId i, SiteId ignore) const {
    if(!enabled_.at(i) || (i!=ignore && !occupancy_.at(i).empty())) return false;
    for(auto j:exclusions_.at(i)) if(j!=ignore && !occupancy_[j].empty()) return false;
    return true;
}
bool CrystalSurface::canAdsorb(SiteId i) const { return available(i,sites_.size()); }
bool CrystalSurface::canHop(std::size_t k) const {
    const auto& e=edges_.at(k);
    // A periodic self-image edge carries geometry but is not a state-changing event.
    return e.from!=e.to && !occupancy_[e.from].empty() && available(e.to,e.from);
}
Change CrystalSurface::changed() const {
    Change c;
    for(SiteId i=0;i<sites_.size();++i) c.affected.push_back(i);
    return c;
}
Change CrystalSurface::adsorb(SiteId i,const std::string& s) {
    if(s.empty() || !canAdsorb(i)) throw std::logic_error("adsorption is not allowed");
    occupancy_[i]=s; return changed();
}
Change CrystalSurface::desorb(SiteId i) {
    if(occupancy_.at(i).empty()) throw std::logic_error("site is empty");
    occupancy_[i].clear(); return changed();
}
Change CrystalSurface::hop(std::size_t k) {
    if(!canHop(k)) throw std::logic_error("hop is not allowed");
    const auto e=edges_[k];
    occupancy_[e.to]=occupancy_[e.from]; occupancy_[e.from].clear(); return changed();
}
void CrystalSurface::writeCSV(const std::string& p) const {
    std::ofstream a(p+"_atoms.csv"),s(p+"_sites.csv"),e(p+"_edges.csv");
    if(!a || !s || !e) throw std::runtime_error("cannot open CSV output");
    a<<std::setprecision(17)<<"id,x,y,z,layer\n";
    s<<std::setprecision(17)<<"id,type,x,y,z,enabled\n";
    e<<std::setprecision(17)<<"from,to,sx,sy,distance\n";
    for(std::size_t i=0;i<atoms_.size();++i) {auto q=atoms_[i];a<<i<<','<<q.position.x<<','<<q.position.y<<','<<q.position.z<<','<<q.layer<<'\n';}
    for(std::size_t i=0;i<sites_.size();++i) {auto q=sites_[i];s<<i<<','<<name(q.type)<<','<<q.position.x<<','<<q.position.y<<','<<q.position.z<<','<<enabled_[i]<<'\n';}
    for(auto q:edges_) e<<q.from<<','<<q.to<<','<<q.sx<<','<<q.sy<<','<<q.distance<<'\n';
    a.flush();s.flush();e.flush();
    if(!a || !s || !e) throw std::runtime_error("CSV write failed");
}
}
