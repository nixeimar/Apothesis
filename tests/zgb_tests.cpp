#include "../src/config/configuration.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void check(bool ok,const std::string& message) {if (!ok) throw std::runtime_error(message);}
int main(int argc,char** argv) {
    try {
        ZGB::Options o;o.x=4;o.y=4;o.end=100;o.seed=1234;
        ZGB::Simulation initial(o);
        auto s=initial.statistics();check(s.vacantSites==16 && s.vacantPairs==32 && s.reactionPairs==0,"initial catalogs must count each undirected edge once");
        check(initial.catalogConsistent(),"initial catalog invalid");
        auto terrace=std::vector<ZGB::Column>(16);terrace[5].height=11;
        ZGB::Simulation heightPairs(o,terrace);
        check(heightPairs.statistics().vacantPairs==28,"raising one empty column must remove exactly its four lateral pair interactions");
        // A center and its four unlike, equal-height neighbors are five active
        // sites, four undirected bonds, and eight directed pair orientations.
        auto cross=std::vector<ZGB::Column>(16,{10,ZGB::Adsorbate::B});cross[5].adsorbate=ZGB::Adsorbate::AB;
        auto crossOptions=o;crossOptions.maxEvents=1;
        ZGB::Simulation centerSelection(crossOptions,cross);
        check(centerSelection.statistics().reactionPairs==4 && centerSelection.statistics().reactiveSites==5,"four-neighbor center interaction counting");
        std::mt19937_64 crossRng(o.seed);
        const double centerTime=-std::log(std::generate_canonical<double,53>(crossRng))/(5*o.reactionRate);
        centerSelection.run();check(std::abs(centerSelection.time()-centerTime)<1e-14,"site-based reaction propensity");
        try {initial.apply(ZGB::Event::AdsorbB2,0,5);throw std::logic_error("diagonal pair accepted");} catch(const std::runtime_error&) {}
        initial.apply(ZGB::Event::AdsorbB2,0,1);
        check(initial.surface()[0].adsorbate==ZGB::Adsorbate::B && initial.surface()[1].adsorbate==ZGB::Adsorbate::B,"dissociation must occupy two sites");
        check(initial.statistics().vacantSites==14,"B2 vacancy balance");
        initial.apply(ZGB::Event::AdsorbAB,2);
        initial.apply(ZGB::Event::React,1,2); // precursor is second; only it grows
        check(initial.surface()[1].height==10 && initial.surface()[2].height==11,"only the precursor column must grow");
        check(initial.surface()[1].adsorbate==ZGB::Adsorbate::Vacant && initial.surface()[2].adsorbate==ZGB::Adsorbate::Vacant,"reaction must free both adsorption sites");
        check(initial.catalogConsistent(),"catalog stale after growth");
        try {initial.apply(ZGB::Event::AdsorbB2,1,2);throw std::logic_error("unequal-height dissociation accepted");} catch(const std::runtime_error&) {}
        initial.apply(ZGB::Event::AdsorbAB,1);initial.apply(ZGB::Event::AdsorbAB,2);
        try {initial.apply(ZGB::Event::React,1,2);throw std::logic_error("invalid reaction accepted");} catch(const std::runtime_error&) {}
        auto unequal=std::vector<ZGB::Column>(16);
        unequal[0]={10,ZGB::Adsorbate::AB};unequal[1]={11,ZGB::Adsorbate::B};
        ZGB::Simulation blocked(o,unequal);
        check(blocked.statistics().reactionPairs==0,"unequal-height reaction eligible");
        try {blocked.apply(ZGB::Event::React,0,1);throw std::logic_error("unequal-height reaction accepted");} catch(const std::runtime_error&) {}
        o.growth=false;ZGB::Simulation catalytic(o,unequal);
        catalytic.apply(ZGB::Event::React,0,1);
        check(catalytic.surface()[0].height==10 && catalytic.surface()[1].height==11,"catalytic reference must not grow");
        check(catalytic.catalogConsistent(),"reference catalog invalid");
        o.growth=true;
        // Periodic edges must be included and refreshed.
        ZGB::Simulation periodic(o);periodic.apply(ZGB::Event::AdsorbAB,0);periodic.apply(ZGB::Event::AdsorbB2,3,7);
        check(periodic.statistics().reactionPairs==1,"periodic nearest-neighbor reaction missing");periodic.apply(ZGB::Event::React,0,3);
        check(periodic.catalogConsistent(),"periodic edge not refreshed");
        // Inventory and deposited material balance after every event.
        ZGB::Simulation run(o);std::size_t observations=0;double lastTime=0,lastHeight=10,integral=0;
        run.run([&](const ZGB::Simulation& state,bool) {
            check(state.catalogConsistent(),"incremental catalog differs from full recomputation");
            auto stat=state.statistics();
            integral+=(stat.time-lastTime)*lastHeight;
            if (stat.time>0) check(std::abs(stat.timeAverageHeight-integral/stat.time)<1e-10,"exact time-weighted height average");
            lastTime=stat.time;lastHeight=stat.meanHeight;
            check(std::abs(stat.coverageAB-(static_cast<double>(stat.adsorbedAB)-stat.reactions)/16)<1e-12,"AB balance");
            check(std::abs(stat.coverageB-(2.0*stat.adsorbedB2-stat.reactions)/16)<1e-12,"B balance");
            check(std::abs(stat.meanHeight-(10+stat.reactions/16.0))<1e-12,"solid incorporation balance");
            check(std::abs(stat.coverageAB+stat.coverageB+stat.coverageVacant-1)<1e-12,"coverage balance");++observations;
        });check(observations>2,"no events validated");
        ZGB::Simulation repeat(o);repeat.run();check(repeat.statistics().events==run.statistics().events && repeat.time()==run.time(),"seed reproducibility");
        // Initial escape rate on a flat empty lattice is exactly N for standard ZGB attempts.
        o.maxEvents=1;ZGB::Simulation timed(o);std::mt19937_64 rng(o.seed);
        double expected=-std::log(std::generate_canonical<double,53>(rng))/16;timed.run();
        check(std::abs(timed.time()-expected)<1e-14,"B2 pair normalization or pre-event clock is wrong");
        auto ordered=o;ordered.pairConvention="ordered_pairs";
        ZGB::Simulation orderedClock(ordered);std::mt19937_64 orderedRng(ordered.seed);
        const double orderedTime=-std::log(std::generate_canonical<double,53>(orderedRng))/(16*(ordered.yAB+4*(1-ordered.yAB)));
        orderedClock.run();check(std::abs(orderedClock.time()-orderedTime)<1e-14,"ordered-pair audit normalization");
        o.end=1e-12;ZGB::Simulation horizon(o);horizon.run();check(horizon.statistics().events==0 && horizon.time()==o.end,"event executed beyond horizon");
        o.end=100;o.maxEvents=0;o.yAB=1;ZGB::Simulation abOnly(o);abOnly.run();s=abOnly.statistics();
        check(s.poisoned && s.coverageAB==1 && s.reactions==0 && s.meanHeight==10,"AB endpoint");
        o.yAB=0;ZGB::Simulation bOnly(o);bOnly.run();s=bOnly.statistics();check(s.poisoned && s.coverageAB==0 && s.reactions==0,"B2 endpoint");
        check(argc==2,"expected input fixture");auto config=Config::readJson(argv[1]);check(config.zgb.has_value(),"named model not parsed");
        check(config.zgb->x==100 && config.zgb->initialHeight==10 && config.zgb->reactionRate==1e5,"paper defaults");
        o.yAB=1.1;try {ZGB::validate(o);throw std::logic_error("invalid y accepted");} catch(const std::runtime_error&) {}
        std::cout<<"All ZGB growth tests passed.\n";return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
