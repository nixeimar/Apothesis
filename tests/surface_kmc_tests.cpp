#include "../src/config/configuration.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
using namespace SurfaceKMC;
void check(bool ok,const char* msg) {if (!ok) throw std::runtime_error(msg);}
int main(int argc,char** argv) {
 try {
    check(argc==2,"expected explicit example");auto config=Config::readJson(argv[1]);check(config.surfaceKMC.has_value() && !config.zgb3d,"generic dispatch");
    auto o=*config.surfaceKMC;o.geometry.repetitions={{4,4,12}};o.geometry.occupiedUnitCells=2;
    Simulation s(o);check(s.eligibleCount(0)==16 && s.eligibleCount(1)==32 && s.eligibleCount(2)==0,"initial catalogs");
    s.apply(1,32,33);s.apply(0,34);check(s.eligibleCount(2)==2,"reactive center counting");
    s.apply(2,34,33);check(s.solid().sites[34].occupied && s.solid().sites[34].species=="AB2","explicit solid product");
    check(s.adsorbates()[33].empty() && s.eligibleCount(0)==15,"vacancies move to next real crystal site");
    check(s.catalogConsistent(),"post growth catalog");
    try {s.apply(1,33,34);throw std::logic_error("unequal height event accepted");}catch(const std::runtime_error&) {}
    // Independent paper-model oracle verifies same-height eligibility after
    // each generic event, using the actual 3D occupancy only in this test.
    Simulation paper(o);paper.run([&](const Simulation& state,bool) {
        ZGB::Options reference;reference.x=4;reference.y=4;
        std::vector<ZGB::Column> columns(16,{0,ZGB::Adsorbate::Vacant});
        for (const auto& site:state.solid().sites) {
            auto& c=columns[site.id%16];c.height+=site.occupied;
            const auto& ads=state.adsorbates()[site.id];
            if(ads=="AB") c.adsorbate=ZGB::Adsorbate::AB;
            else if(ads=="B") c.adsorbate=ZGB::Adsorbate::B;
        }
        ZGB::Simulation oracle(reference,columns);auto expected=oracle.statistics();
        check(state.eligibleCount(0)==expected.vacantSites && state.eligibleCount(1)==expected.vacantPairs && state.eligibleCount(2)==expected.reactiveSites,"explicit patterns differ from paper eligibility");
    });
    auto cap=o;cap.geometry.repetitions[2]=3;Simulation capacity(cap);
    capacity.apply(0,32);capacity.apply(1,33,37);capacity.apply(2,32,33);capacity.run();
    check(capacity.capacityLimited() && capacity.stopReason()=="capacity_limit","finite headroom not reported");
    // Rename every chemical species and define extra, independently counted channels.
    for (auto& rule:o.processes) {
        for (auto& before:rule.before) if (!before.empty()) before="renamed_"+before;
        for (auto& after:rule.after) if (!after.species.empty()) after.species="renamed_"+after.species;
    }
    o.elements={{"substrate","X"},{"renamed_AB","X"},{"renamed_B","X"},{"renamed_AB2","X"}};
    Rule des;des.id="desorb_precursor";des.type="desorption";des.before={"renamed_AB"};des.after={{"",false}};des.rate=2;o.processes.push_back(des);
    Rule hop;hop.id="hop_precursor";hop.type="diffusion";hop.counting="undirected_pairs";hop.before={"renamed_AB",""};hop.after={{"",false},{"renamed_AB",false}};hop.rate=3;o.processes.push_back(hop);
    Rule duplicate=o.processes[0];duplicate.id="second_adsorption_channel";duplicate.rate=0;o.processes.push_back(duplicate);
    Simulation renamed(o);renamed.apply(0,32);renamed.apply(4,32,33);check(renamed.adsorbates()[32].empty() && renamed.adsorbates()[33]=="renamed_AB","generic diffusion");
    renamed.apply(3,33);check(renamed.adsorbates()[33].empty() && renamed.solid().atomCount()==32,"generic desorption");
    check(renamed.catalogConsistent() && renamed.counts()[0]==1 && renamed.counts()[3]==1 && renamed.counts()[4]==1 && renamed.counts()[5]==0,"independent process IDs");
    Simulation evolution(o);double last=0,height=2,integral=0;std::size_t observations=0;
    evolution.run([&](const Simulation& state,bool) {
        check(state.catalogConsistent(),"incremental catalog mismatch");
        integral+=(state.time()-last)*height;
        if(state.time()>0) check(std::abs(state.timeAverageHeight()-integral/state.time())<1e-9,"residence time integral");
        last=state.time();height=state.meanHeight();++observations;
        check(state.solid().atomCount()==32+state.counts()[2],"solid incorporation balance");
        std::size_t ab=0,b=0;
        for (const auto& a:state.adsorbates()) {ab+=a=="renamed_AB";b+=a=="renamed_B";}
        check(ab+state.counts()[2]+state.counts()[3]==state.counts()[0]+state.counts()[5],"precursor inventory");
        check(b+state.counts()[2]==2*state.counts()[1],"dissociated inventory");
        for (const auto& site:state.solid().sites) {
            if (site.occupied || !state.adsorbates()[site.id].empty()) check(state.solid().occupiedCoordination(site.id)>0,"unsupported particle");
            if (site.fixed) check(site.occupied && site.species=="substrate","fixed substrate mutation");
        }
    });
    check(observations>10 && evolution.counts()[2]>0,"chemistry not exercised");
    Simulation repeat(o);repeat.run();check(repeat.time()==evolution.time() && repeat.counts()==evolution.counts(),"seed reproducibility");
    auto clock=o;clock.processes.resize(2);clock.maxEvents=1;Simulation timed(clock);std::mt19937_64 rng(clock.seed);
    const double expected=-std::log(std::generate_canonical<double,53>(rng))/(16*.4+32*.3);timed.run();check(std::abs(timed.time()-expected)<1e-14,"explicit process propensities / pre-event clock");
    auto horizon=o;horizon.end=1e-12;Simulation shortRun(horizon);shortRun.run();check(shortRun.steps()==0 && shortRun.time()==horizon.end,"horizon overshoot");
    std::ostringstream xyz;evolution.writeXYZ(xyz);check(xyz.str().find("renamed_AB2")!=std::string::npos,"generic species output");
    std::cout<<"All explicit surface process tests passed.\n";return 0;
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
