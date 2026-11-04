#include "../src/config/configuration.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

void check(bool ok,const std::string& message) {if (!ok) throw std::runtime_error(message);}
Atomistic::RunOptions runOptions(Atomistic::Kind kind) {
    Atomistic::RunOptions o;o.start=2;o.end=3;o.temperature=600;o.seed=1234;
    Atomistic::Process p;p.id="test";p.kind=kind;p.species="Cu";p.element="Cu";p.value=1;
    p.minCoordination=kind==Atomistic::Kind::Desorption ? 0 : 1;
    o.processes={p};return o;
}
int main(int argc,char** argv) {
    try {
        Atomistic::Options geometry;
        geometry.crystal="SimpleCubic";geometry.species="Cu";geometry.element="Cu";
        geometry.latticeConstant=2;geometry.repetitions={{4,4,4}};geometry.occupiedUnitCells=2;geometry.fixedUnitCells=1;
        auto surface=Atomistic::build(geometry);
        auto original=surface;
        check(surface.atomCount()==32,"initial atom count");
        Atomistic::Simulation adsorption(surface,runOptions(Atomistic::Kind::Adsorption));
        auto candidates=adsorption.events();check(candidates.size()==16,"only supported sites in first growth plane may adsorb");
        auto event=candidates[0];adsorption.apply(event);
        check(adsorption.surface().atomCount()==33,"adsorption must add one atom");
        check(adsorption.surface().sites[event.target].species=="Cu","adsorption species");
        try {adsorption.apply(event);throw std::logic_error("stale event accepted");} catch(const std::runtime_error&) {}
        Atomistic::Simulation diffusion(adsorption.surface(),runOptions(Atomistic::Kind::Diffusion));
        auto moves=diffusion.events();check(!moves.empty(),"supported diffusion missing");
        auto move=std::find_if(moves.begin(),moves.end(),[&](const Atomistic::Event& e){return e.source==event.target;});
        check(move!=moves.end(),"deposited atom cannot diffuse");
        auto hop=*move;diffusion.apply(hop);
        check(diffusion.surface().atomCount()==33,"diffusion must conserve atoms");
        check(!diffusion.surface().sites[hop.source].occupied && diffusion.surface().sites[hop.target].occupied,"diffusion occupancy");
        for (const auto& s : diffusion.surface().sites)
            check(s.position==original.sites[s.id].position,"events must not move crystal coordinates");
        Atomistic::Simulation desorption(diffusion.surface(),runOptions(Atomistic::Kind::Desorption));
        auto removals=desorption.events();
        auto removal=std::find_if(removals.begin(),removals.end(),[&](const Atomistic::Event& e){return e.source==hop.target;});
        check(removal!=removals.end(),"desorption missing");desorption.apply(*removal);
        check(desorption.surface().atomCount()==32,"desorption must remove one atom");
        for (const auto& s : desorption.surface().sites) if (s.fixed) check(s.occupied,"fixed atom changed");
        auto frozenGeometry=geometry;
        frozenGeometry.occupiedUnitCells=4;frozenGeometry.fixedUnitCells=4;frozenGeometry.vacancies={48};
        for (auto kind : {Atomistic::Kind::Adsorption,Atomistic::Kind::Diffusion,Atomistic::Kind::Desorption}) {
            Atomistic::Simulation frozen(Atomistic::build(frozenGeometry),runOptions(kind));
            check(frozen.events().empty(),"fixed site accepted an event");
        }
        // A mobile supporting atom must not strand a singly bonded neighbor.
        auto supported=surface;
        supported.sites[32].occupied=true; // bonded only to site 16 below it
        supported.sites[17].occupied=false; // a supported hop destination for 16
        for (auto kind : {Atomistic::Kind::Desorption,Atomistic::Kind::Diffusion}) {
            Atomistic::Simulation supportCheck(supported,runOptions(kind));
            for (const auto& e : supportCheck.events())
                check(e.source!=16,"event would remove another atom's last bond");
        }
        // A simultaneous hop can replace a dependent atom's last bond.
        // Two equilateral triangles share the source-destination edge.
        Atomistic::Surface triangles;
        triangles.periodic={{false,false,false}};triangles.neighborCutoff=1.01;triangles.element="Cu";
        triangles.sites={
            {0,{{0,0,0}},"Cu",true,false,{1,2,3}},
            {1,{{0.5,std::sqrt(3.0)/2,0}},"Cu",true,false,{0,2}},
            {2,{{1,0,0}},"Cu",false,false,{0,1,3}},
            {3,{{0.5,-std::sqrt(3.0)/2,0}},"Cu",true,true,{0,2}}
        };
        Atomistic::Simulation rescue(triangles,runOptions(Atomistic::Kind::Diffusion));
        auto rescueEvents=rescue.events();
        auto rescued=std::find_if(rescueEvents.begin(),rescueEvents.end(),[](const Atomistic::Event& e){return e.source==0 && e.target==2;});
        check(rescued!=rescueEvents.end(),"valid bond-preserving hop was rejected");
        rescue.apply(*rescued);
        for (const auto& s:rescue.surface().sites) if (s.occupied)
            check(rescue.surface().occupiedCoordination(s.id)>=1,"hop left an isolated atom");
        // A detached dimer has one bond per atom, but is no longer film.
        Atomistic::Surface branch;
        branch.periodic={{false,false,false}};branch.neighborCutoff=1.01;branch.element="Cu";
        branch.sites={
            {0,{{0,0,0}},"Cu",true,true,{1,4,5,6}},
            {1,{{1,0,0}},"Cu",true,false,{0,2,5,6}},
            {2,{{2,0,0}},"Cu",true,false,{1,3,6}},
            {3,{{3,0,0}},"Cu",true,false,{2}},
            {4,{{0,1,0}},"Cu",true,true,{0}},
            {5,{{1,1,0}},"Cu",false,false,{0,1}},
            {6,{{2,1,0}},"Cu",false,false,{0,1,2}}
        };
        Atomistic::Simulation anchorRemoval(branch,runOptions(Atomistic::Kind::Desorption));
        for (const auto& e:anchorRemoval.events()) check(e.source!=1,"desorption detached a bonded cluster");
        Atomistic::Simulation anchorHop(branch,runOptions(Atomistic::Kind::Diffusion));
        bool foundBridgeHop=false;
        for (const auto& e:anchorHop.events()) {
            check(!(e.source==1 && e.target==5),"diffusion detached a bonded cluster");
            if (e.source==1 && e.target==6) foundBridgeHop=true;
        }
        check(foundBridgeHop,"cluster-preserving replacement bridge hop rejected");
        auto detached=branch;detached.sites[1].occupied=false;
        try {Atomistic::Simulation bad(detached,runOptions(Atomistic::Kind::Adsorption));throw std::logic_error("initial floating bonded cluster accepted");}
        catch(const std::runtime_error& e) {check(std::string(e.what()).find("fixed substrate")!=std::string::npos,"initial connectivity diagnostic");}

        auto isolated=surface;isolated.sites[48].occupied=true;
        try {Atomistic::Simulation bad(isolated,runOptions(Atomistic::Kind::Adsorption));throw std::logic_error("isolated initial atom accepted");}
        catch(const std::runtime_error& e) {check(std::string(e.what()).find("no occupied first-nearest")!=std::string::npos,"initial bond diagnostic");}
        for (auto kind : {Atomistic::Kind::Adsorption,Atomistic::Kind::Diffusion}) {
            auto unsupported=runOptions(kind);unsupported.processes[0].minCoordination=0;
            try {Atomistic::validate(unsupported);throw std::logic_error("zero bond minimum accepted");}
            catch(const std::runtime_error& e) {check(std::string(e.what()).find("min_coordination")!=std::string::npos,"bond minimum diagnostic");}
        }
        // Closed bulk vacancies cannot adsorb, and buried atoms cannot desorb.
        geometry.occupiedUnitCells=4;geometry.vacancies={21};
        Atomistic::Simulation sealed(Atomistic::build(geometry),runOptions(Atomistic::Kind::Adsorption));
        check(sealed.events().empty(),"adsorption into sealed bulk vacancy");
        Atomistic::Simulation buried(Atomistic::build(geometry),runOptions(Atomistic::Kind::Desorption));
        for (auto e : buried.events()) check(e.source>=48,"buried atom desorption");
        // Rate units and optional bond-counting barrier.
        auto options=runOptions(Atomistic::Kind::Adsorption);
        options.processes[0].model="arrhenius";options.processes[0].prefactor=100;
        options.processes[0].activationEnergy=10000;options.processes[0].bondEnergy=2000;
        Atomistic::Simulation activated(surface,options);
        check(std::abs(activated.events()[0].rate-100*std::exp(-12000/(8.31446261815324*600)))<1e-12,"Arrhenius barrier");
        // Sample the clock before changing the number of eligible sites.
        auto timedOptions=runOptions(Atomistic::Kind::Adsorption);timedOptions.maxEvents=1;
        std::mt19937_64 rng(timedOptions.seed);
        const double expectedTime=timedOptions.start-std::log(std::generate_canonical<double,53>(rng))/16;
        Atomistic::Simulation timed(surface,timedOptions);timed.run();
        check(timed.steps()==1 && std::abs(timed.time()-expectedTime)<1e-14,"waiting time must use the pre-event total rate");
        // Species identity is preserved and written with its own chemical element.
        auto nickelOptions=runOptions(Atomistic::Kind::Adsorption);
        nickelOptions.processes[0].species="Ni";nickelOptions.processes[0].element="Ni";
        Atomistic::Simulation nickel(surface,nickelOptions);auto nickelEvent=nickel.events()[0];nickel.apply(nickelEvent);
        check(nickel.surface().sites[nickelEvent.target].species=="Ni" && nickel.surface().elements.at("Ni")=="Ni","deposited species identity");
        Atomistic::Simulation copperOnly(nickel.surface(),runOptions(Atomistic::Kind::Desorption));
        for (const auto& e:copperOnly.events()) check(e.source!=nickelEvent.target,"desorption selected the wrong species");
        // A sampled event beyond the horizon must not be executed.
        options=runOptions(Atomistic::Kind::Adsorption);options.end=options.start+1e-9;
        Atomistic::Simulation horizon(surface,options);horizon.run();
        check(horizon.steps()==0 && horizon.time()==options.end && horizon.surface().atomCount()==32,"end-time overshoot");
        options=runOptions(Atomistic::Kind::Adsorption);options.processes[0].value=0;
        Atomistic::Simulation zero(surface,options);zero.run();
        check(zero.steps()==0 && zero.time()==options.end && zero.stopReason()=="no_events","zero-rate termination");
        check(argc==2,"expected example input");const auto config=Config::readJson(argv[1]);
        check(config.atomisticSimulation.has_value(),"atomistic simulation not parsed");
        auto a=Atomistic::Simulation(Atomistic::build(*config.atomisticSurface),*config.atomisticSimulation);
        auto b=Atomistic::Simulation(Atomistic::build(*config.atomisticSurface),*config.atomisticSimulation);
        std::size_t checkedStates=0;
        a.run([&](const Atomistic::Simulation& state,bool) {
            for (const auto& site:state.surface().sites) if (site.occupied)
                check(state.surface().occupiedCoordination(site.id)>=1,"trajectory contains an isolated atom");
            std::vector<bool> anchored(state.surface().sites.size(),false);
            std::vector<std::size_t> pending;
            for (const auto& site:state.surface().sites) if (site.fixed && site.occupied) {
                anchored[site.id]=true;pending.push_back(site.id);
            }
            for (std::size_t i=0;i<pending.size();++i)
                for (auto j:state.surface().sites[pending[i]].neighbors)
                    if (state.surface().sites[j].occupied && !anchored[j]) {anchored[j]=true;pending.push_back(j);}
            for (const auto& site:state.surface().sites) if (site.occupied)
                check(anchored[site.id],"trajectory contains a floating bonded cluster");
            ++checkedStates;
        });
        check(checkedStates>=a.steps()+1,"trajectory bond checks not executed");
        b.run();check(a.counts()==b.counts() && a.time()==b.time(),"seeded reproducibility");
        for (std::size_t i=0;i<a.surface().sites.size();++i) check(a.surface().sites[i].occupied==b.surface().sites[i].occupied,"seeded final state");
        check(a.steps()>0,"no KMC events");
        check(a.surface().atomCount()==128+a.counts()[0]-a.counts()[2],"KMC atom balance");
        check(a.time()<=config.atomisticSimulation->end,"KMC horizon exceeded");
        std::cout<<"All atomistic KMC tests passed; event counts:";
        for (auto n:a.counts()) std::cout<<' '<<n;
        std::cout<<'\n';return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
