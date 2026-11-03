#include "crystal_surface.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace Crystal;
void check(bool b) {if(!b) throw std::runtime_error("test failed");}
template<class F> void rejects(F f) {bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw);}
int main() {
 try {
    for(auto face:{Face::FCC111,Face::FCC100,Face::HCP0001})
    for(auto dims:{std::pair<int,int>{1,1},{2,3},{3,2},{4,4}}) {
        Definition d;d.face=face;d.nx=dims.first;d.ny=dims.second;
        CrystalSurface s(d);
        std::size_t n=d.nx*d.ny;
        check(s.atoms().size()==3*n);
        check(s.sites().size()==n*(face==Face::FCC100?4:6));
        for(SiteId i=0;i<s.sites().size();++i) {
            auto q=s.sites()[i];
            std::size_t expected=q.type==SiteType::Top?1:q.type==SiteType::Bridge?2:q.type==SiteType::Fourfold?4:3;
            check(q.support.size()==expected);
            check(!s.canAdsorb(i));
            // Each support image has the correct geometric distance.
            double expectedR=q.type==SiteType::Top?0:q.type==SiteType::Bridge?s.spacing()/2:
                q.type==SiteType::Fourfold?s.spacing()/std::sqrt(2.0):s.spacing()/std::sqrt(3.0);
            for(auto image:q.support) {
                auto p=s.atoms().at(image.atom).position;
                p.x+=image.sx*s.cellA().x+image.sy*s.cellB().x;
                p.y+=image.sx*s.cellA().y+image.sy*s.cellB().y;
                check(std::abs(std::hypot(p.x-q.position.x,p.y-q.position.y)-expectedR)<1e-9);
            }
            // Hcp hollow must have a second-layer atom directly below.
            if(q.type==SiteType::HcpHollow) {
                bool below=false;
                for(std::size_t a=n;a<2*n;++a) {
                    auto p=s.atoms()[a].position;
                    if(std::hypot(p.x-q.position.x,p.y-q.position.y)<1e-9) below=true;
                }
                check(below);
            }
        }
        if(face==Face::FCC111) {
            for(auto q:s.sites()) if(q.type==SiteType::FccHollow) {
                bool below=false;
                for(std::size_t a=2*n;a<3*n;++a) {
                    auto p=s.atoms()[a].position;
                    if(std::hypot(p.x-q.position.x,p.y-q.position.y)<1e-9) below=true;
                }
                check(below);
            }
        }
        s.connect(SiteType::Top,SiteType::Top,s.spacing());
        check(s.transitions().size()==n*(face==Face::FCC100?4:6));
        auto count=s.transitions().size();
        s.connect(SiteType::Top,SiteType::Top,s.spacing());
        check(count==s.transitions().size());
        for(auto e:s.transitions()) {
            bool reverse=false;
            for(auto r:s.transitions()) if(r.from==e.to&&r.to==e.from&&r.sx==-e.sx&&r.sy==-e.sy) reverse=true;
            check(reverse);
        }
    }
    Definition d;
    CrystalSurface s(d);
    s.enable(SiteType::Top);s.enable(SiteType::Bridge);
    s.setExclusionDistance(0.6*s.spacing());
    s.connect(SiteType::Top,SiteType::Bridge,0.51*s.spacing());
    auto change=s.adsorb(0,"A");
    check(change.affected.size()==s.sites().size());
    for(auto i:s.exclusions(0)) check(!s.canAdsorb(i));
    rejects([&]{s.adsorb(0,"B");});
    rejects([&]{s.setExclusionDistance(0.1);});
    // Top at origin excludes bridge across the periodic x boundary.
    SiteId wrappedBridge=3*6+1; // last x cell's x-directed bridge, y=0
    check(!s.canAdsorb(wrappedBridge));
    rejects([&]{s.adsorb(wrappedBridge,"B");});
    rejects([&]{s.canAdsorb(s.sites().size());});
    bool hopped=false;
    for(std::size_t i=0;i<s.transitions().size();++i) if(s.transitions()[i].from==0 && s.canHop(i)) {
        auto to=s.transitions()[i].to;s.hop(i);check(s.occupant(0).empty());check(s.occupant(to)=="A");
        s.desorb(to);hopped=true;break;
    }
    check(hopped);
    rejects([&]{s.desorb(0);});
    rejects([&]{s.adsorb(0,"");});
    rejects([&]{s.connect(SiteType::Top,SiteType::Top,-1);});
    rejects([&]{Definition bad;bad.a=0;CrystalSurface invalid(bad);});
    rejects([&]{Definition bad;bad.layers=2;CrystalSurface invalid(bad);});
    std::cout<<"All crystal-surface geometry, image, occupancy and transition tests passed.\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
