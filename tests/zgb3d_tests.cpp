#include "../src/config/configuration.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

void check(bool value,const std::string& message) {if (!value) throw std::runtime_error(message);}
ZGB3D::Options options() {
    ZGB3D::Options o;
    o.geometry.crystal="SimpleCubic";o.geometry.species="substrate";o.geometry.element="X";
    o.geometry.latticeConstant=1;o.geometry.repetitions={{4,4,12}};o.geometry.occupiedUnitCells=2;o.geometry.fixedUnitCells=1;
    o.kinetics.x=4;o.kinetics.y=4;o.kinetics.initialHeight=2;o.kinetics.end=100;o.kinetics.seed=1234;
    return o;
}
int main(int argc,char** argv) {
    try {
        auto o=options();ZGB3D::Simulation model(o);
        check(model.solid().sites.size()==192 && model.solid().atomCount()==32,"3D volume or initial atom inventory missing");
        check(model.statistics().vacantSites==16 && model.statistics().vacantPairs==32,"initial top binding catalog");
        for (std::size_t i=0;i<192;++i) check(model.isBinding(i)==(i>=32 && i<48),"binding nodes must lie immediately above solid");
        auto initialPositions=model.solid().sites;
        model.apply(ZGB::Event::AdsorbAB,32);
        check(model.solid().atomCount()==32 && !model.solid().sites[32].occupied,"adsorption must not create a solid particle");
        check(model.adsorbates()[32]==ZGB::Adsorbate::AB,"AB chemical state missing");
        model.apply(ZGB::Event::AdsorbB2,33,34);
        check(model.adsorbates()[33]==ZGB::Adsorbate::B && model.adsorbates()[34]==ZGB::Adsorbate::B,"B2 did not dissociate on two 3D nodes");
        model.apply(ZGB::Event::React,33,32);
        check(model.solid().atomCount()==33 && model.solid().sites[32].occupied && model.solid().sites[32].species=="AB2","reaction must incorporate a real 3D solid unit at the precursor site");
        check(!model.isBinding(32) && model.isBinding(48) && model.isBinding(33),"binding site did not move above deposited particle");
        check(model.adsorbates()[32]==ZGB::Adsorbate::Vacant && model.adsorbates()[33]==ZGB::Adsorbate::Vacant,"reaction did not consume both adsorbates");
        try {model.apply(ZGB::Event::AdsorbB2,48,33);throw std::logic_error("different-height 3D pair accepted");} catch(const std::runtime_error&) {}
        try {model.apply(ZGB::Event::AdsorbB2,35,38);throw std::logic_error("diagonal 3D pair accepted");} catch(const std::runtime_error&) {}
        check(model.catalogConsistent(),"3D catalogs stale after incorporation");
        std::ostringstream xyz;model.writeXYZ(xyz);
        check(xyz.str().find("AB2_solid solid")!=std::string::npos && xyz.str().find("B_adsorbate adsorbate")!=std::string::npos,"XYZ must distinguish solid and adsorbate particles");
        check(xyz.str().find("chemical_state:S:1:role:S:1")!=std::string::npos,"XYZ lacks chemical-role metadata");
        // Periodic lateral pairs are recognized in 3D, with no column-engine calls.
        ZGB3D::Simulation periodic(o);periodic.apply(ZGB::Event::AdsorbAB,32);periodic.apply(ZGB::Event::AdsorbB2,35,39);
        periodic.apply(ZGB::Event::React,32,35);check(periodic.solid().sites[32].occupied,"periodic reaction failed");
        // Finite z capacity must not be called poisoning.
        auto small=o;small.geometry.repetitions[2]=3;
        ZGB3D::Simulation capped(small);capped.apply(ZGB::Event::AdsorbAB,32);capped.apply(ZGB::Event::AdsorbB2,33,34);capped.apply(ZGB::Event::React,32,33);capped.run();
        check(capped.capacityLimited() && capped.statistics().stopReason=="capacity_limit" && !capped.statistics().poisoned,"finite growth capacity was mistaken for poisoning");
        ZGB3D::Simulation kinetic(o);std::size_t observed=0;double lastTime=0,lastHeight=2,integral=0;
        kinetic.run([&](const ZGB3D::Simulation& state,bool) {
            check(state.catalogConsistent(),"incremental 3D catalogs inconsistent");
            const auto s=state.statistics();const auto& volume=state.solid();
            check(volume.atomCount()==32+s.reactions,"3D solid inventory balance");
            check(std::abs(s.meanHeight-(2+s.reactions/16.0))<1e-12,"mean height must derive from 3D solid inventory");
            check(std::abs(s.coverageAB-(static_cast<double>(s.adsorbedAB)-s.reactions)/16)<1e-12,"3D AB balance");
            check(std::abs(s.coverageB-(2.0*s.adsorbedB2-s.reactions)/16)<1e-12,"3D B balance");
            std::size_t products=0,ab=0,b=0;
            for (const auto& node:volume.sites) {
                if (node.occupied && node.species=="AB2") ++products;
                if (state.adsorbates()[node.id]==ZGB::Adsorbate::AB) ++ab;
                if (state.adsorbates()[node.id]==ZGB::Adsorbate::B) ++b;
            }
            check(ab+products==s.adsorbedAB && ab+b+2*products==s.adsorbedAB+2*s.adsorbedB2,"symbolic A/B atom inventory not conserved");

            integral+=(s.time-lastTime)*lastHeight;
            if (s.time>0) check(std::abs(s.timeAverageHeight-integral/s.time)<1e-10,"3D residence-time height integral");
            lastTime=s.time;lastHeight=s.meanHeight;
            for (const auto& site:volume.sites) {
                check(site.position==initialPositions[site.id].position,"fixed-site coordinates changed");
                if (site.fixed) check(site.occupied && site.species=="substrate","fixed substrate changed");
                if (site.occupied && !site.fixed) check(volume.occupiedCoordination(site.id)>=1,"unbonded solid particle");
                if (state.adsorbates()[site.id]!=ZGB::Adsorbate::Vacant) {
                    check(!site.occupied && state.isBinding(site.id),"adsorbate overlaps solid or is not exposed");
                    check(volume.occupiedCoordination(site.id)>=1,"adsorbate lacks a solid bond");
                }
            }
            // Project only for a test oracle: compare eligible catalogs to the
            // independent paper implementation, never use it to advance 3D.
            if (!state.capacityLimited()) {
                std::vector<ZGB::Column> projection(16,{0,ZGB::Adsorbate::Vacant});
                for (const auto& site:volume.sites) {
                    const auto xy=site.id%16;
                    if (site.occupied) ++projection[xy].height;
                    if (state.isBinding(site.id)) projection[xy].adsorbate=state.adsorbates()[site.id];
                }
                auto referenceOptions=o.kinetics;ZGB::Simulation oracle(referenceOptions,projection);auto expected=oracle.statistics();
                check(s.vacantSites==expected.vacantSites && s.vacantPairs==expected.vacantPairs && s.reactionPairs==expected.reactionPairs && s.reactiveSites==expected.reactiveSites,"3D eligibility differs from same-height four-neighbor rules");
            }
            ++observed;
        });
        check(observed>2 && kinetic.statistics().reactions>0,"3D evolution was not exercised");
        ZGB3D::Simulation repeat(o);repeat.run();check(repeat.statistics().events==kinetic.statistics().events && repeat.time()==kinetic.time(),"3D seed reproducibility");
        for (std::size_t id=0;id<192;++id) check(repeat.solid().sites[id].occupied==kinetic.solid().sites[id].occupied && repeat.adsorbates()[id]==kinetic.adsorbates()[id],"seeded 3D final state mismatch");
        auto reference=o;reference.kinetics.growth=false;ZGB3D::Simulation catalytic(reference);catalytic.run();
        check(catalytic.solid().atomCount()==32,"3D catalytic reference deposited solid");
        auto cavity=o;cavity.geometry.vacancies={0};
        try {ZGB3D::Simulation invalid(cavity);throw std::logic_error("buried substrate vacancy accepted");} catch(const std::runtime_error&) {}
        auto tiny=o;tiny.kinetics.end=1e-12;ZGB3D::Simulation horizon(tiny);horizon.run();check(horizon.statistics().events==0 && horizon.time()==tiny.kinetics.end,"3D horizon overshoot");
        check(argc==2,"expected 3D input example");auto config=Config::readJson(argv[1]);
        check(config.zgb3d.has_value() && !config.zgb.has_value(),"3D input selected the height backend");
        std::cout<<"All explicit 3D ZGB tests passed.\n";return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
