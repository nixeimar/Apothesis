#include "crystal_surface.h"
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace Crystal;
int main(int argc,char** argv) {
    try {
        Definition d;
        std::string face=argc>1?argv[1]:"fcc111";
        if(face=="fcc100") d.face=Face::FCC100;
        else if(face=="hcp0001") d.face=Face::HCP0001;
        else if(face!="fcc111") throw std::invalid_argument("use fcc111, fcc100 or hcp0001");
        CrystalSurface s(d);
        s.enable(SiteType::Top);
        s.enable(d.face==Face::FCC100?SiteType::Fourfold:SiteType::FccHollow);
        // Toy model: top and selected hollow sites metastable; bridges disabled.
        s.setExclusionDistance(0.6*s.spacing());
        s.connect(SiteType::Top,d.face==Face::FCC100?SiteType::Fourfold:SiteType::FccHollow,
                  0.8*s.spacing());
        s.writeCSV(argc>2?argv[2]:face);
        // Direct-method continuous-time KMC demonstration, dimensionless toy rates.
        // Rebuild all eligible events after each change (no stale event cache).
        struct Event {int kind; std::size_t id; double rate;};
        std::mt19937 rng(12345);
        double time=0;
        for(int step=0;step<1000;++step) {
            std::vector<Event> events;
            double total=0;
            auto add=[&](int k,std::size_t i,double r){events.push_back({k,i,r});total+=r;};
            for(SiteId i=0;i<s.sites().size();++i) {
                if(s.canAdsorb(i)) add(0,i,1.0);
                if(!s.occupant(i).empty()) add(1,i,0.2);
            }
            for(std::size_t i=0;i<s.transitions().size();++i) if(s.canHop(i)) add(2,i,2.0);
            if(total==0) break;
            auto uniform=[&](){return (static_cast<double>(rng())+0.5)/(static_cast<double>(rng.max())+1.0);};
            time-=std::log(uniform())/total;
            double target=uniform()*total;
            Event chosen=events.back();
            for(auto ev:events) {target-=ev.rate;if(target<0){chosen=ev;break;}}
            if(chosen.kind==0) s.adsorb(chosen.id,"X");
            else if(chosen.kind==1) s.desorb(chosen.id);
            else s.hop(chosen.id);
        }
        int occupied=0;
        for(SiteId i=0;i<s.sites().size();++i) {
            if(!s.occupant(i).empty()) {
                ++occupied;
                for(auto j:s.exclusions(i)) if(!s.occupant(j).empty())
                    throw std::logic_error("exclusion invariant violated");
            }
        }
        std::cout<<face<<": "<<s.atoms().size()<<" atoms, "<<s.sites().size()
                 <<" sites, "<<s.transitions().size()<<" directed edges; toy KMC t="
                 <<time<<", occupied="<<occupied<<"\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
