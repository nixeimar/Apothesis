#include "../src/config/configuration.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <set>
using namespace SurfaceKMC;
void check(bool yes,const char* message) {if(!yes) throw std::runtime_error(message);}
Options base() {
    Options o;o.surfaceSites="gas_exposed";o.geometry.crystal="SimpleCubic";o.geometry.species="S";o.geometry.element="X";
    o.geometry.latticeConstant=1;o.geometry.repetitions={{4,4,5}};o.geometry.occupiedUnitCells=2;o.geometry.fixedUnitCells=1;
    o.elements={{"S","X"},{"A","X"},{"B","X"},{"C","X"},{"F","X"}};o.end=2;o.maxEvents=8;
    Rule a;a.id="ads";a.type="adsorption";a.before={""};a.after={{"A",false}};a.rate=2;o.processes={a};return o;
}
Link first(std::size_t a,std::size_t b) {Link link;link.first=a;link.second=b;link.sameHeight=true;return link;}
int main(int argc,char** argv) {
 try {
    for(const auto& crystal:{"SimpleCubic","BCC","FCC","HCP","Diamond"}) {
        auto o=base();o.geometry.crystal=crystal;Simulation s(o);
        check(s.eligibleCount(0)>0,"no supported binding sites on generated crystal");
        s.run([&](const Simulation& state,bool) {check(state.catalogConsistent(),"graph catalog inconsistent");});
        check(s.steps()>0 && s.solid().atomCount()==Atomistic::build(o.geometry).atomCount(),"adsorption changed solid occupancy");
    }
    auto o=base();Rule triple;triple.id="three_site";triple.type="reaction";triple.before={"","",""};triple.after={{"A",false},{"B",false},{"C",false}};triple.counting="embeddings";triple.links={first(0,1),first(1,2)};triple.rate=.3;
    o.processes.push_back(triple);Simulation embeddings(o);check(embeddings.eligibleCount(1)==192,"wrong number of ordered three-site embeddings");
    auto squareOptions=base();Rule square=triple;square.id="square";square.before={"","","",""};square.after={{"A",false},{"A",false},{"A",false},{"A",false}};square.counting="unique_site_sets";square.links={first(0,1),first(1,2),first(2,3),first(3,0)};squareOptions.processes.push_back(square);Simulation squares(squareOptions);check(squares.eligibleCount(1)==24,"four-site cycle constraints/symmetry count");
    auto triangleOptions=base();triangleOptions.geometry.crystal="HCP";Rule triangle=triple;triangle.counting="unique_site_sets";triangle.links.push_back(first(0,2));triangleOptions.processes.push_back(triangle);Simulation triangles(triangleOptions);check(triangles.eligibleCount(1)==32,"HCP basal triangle count");
    auto sets=o;sets.processes[1].counting="unique_site_sets";Simulation unique(sets);check(unique.eligibleCount(1)==96,"symmetric site sets double counted");
    try{unique.apply(1,std::vector<std::size_t>{32,33,32});throw std::logic_error("repeated site accepted");}catch(const std::runtime_error&){}
    check(unique.steps()==0 && unique.adsorbates()[32].empty(),"invalid multi-site event partially applied");
    unique.apply(1,std::vector<std::size_t>{32,33,34});check(unique.adsorbates()[32]=="A" && unique.adsorbates()[33]=="B" && unique.adsorbates()[34]=="C","multi-site products not assigned by slot");
    Rule grow=triple;grow.id="multi_growth";grow.before={"A","B","C"};grow.after={{"F",true},{"F",true},{"F",true}};
    sets.processes.push_back(grow);Simulation grown(sets);grown.apply(1,std::vector<std::size_t>{32,33,34});grown.apply(2,std::vector<std::size_t>{32,33,34});
    check(grown.solid().atomCount()==35 && grown.catalogConsistent(),"simultaneous solid products or exposure update failed");
    // A process-specific second-neighbor distance shell, not the support graph.
    auto diagonal=base();Rule pair;pair.id="diagonal_pair";pair.type="reaction";pair.before={"",""};pair.after={{"B",false},{"B",false}};pair.counting="undirected_pairs";
    Link shell;shell.neighbors="distance";shell.minDistance=1.4;shell.maxDistance=1.43;shell.sameHeight=true;pair.links={shell};diagonal.processes.push_back(pair);
    Simulation secondShell(diagonal);check(secondShell.eligibleCount(1)==32,"second shell count");secondShell.apply(1,32,37);
    Simulation axial(diagonal);try{axial.apply(1,32,33);throw std::logic_error("first neighbor accepted by second shell");}catch(const std::runtime_error&){}
    Simulation wrap(diagonal);wrap.apply(1,32,47);check(wrap.catalogConsistent(),"periodic diagonal match failed");
    // Same-height is optional: a geometrical link can connect a step edge.
    auto step=diagonal;step.geometry.vacancies={16};step.processes[1].links[0].sameHeight=false;Simulation rough(step);rough.apply(1,16,33);
    step.processes[1].links[0].sameHeight=true;Simulation flatOnly(step);try{flatOnly.apply(1,16,33);throw std::logic_error("height restriction ignored");}catch(const std::runtime_error&){}
    // Cavity/overhang topology is accepted; sealed vacancies are not gas binding sites.
    auto cavity=base();cavity.geometry.occupiedUnitCells=3;cavity.geometry.vacancies={16};Simulation sealed(cavity);
    try{sealed.apply(0,16);throw std::logic_error("adsorption inside a sealed cavity");}catch(const std::runtime_error&){}
    // Solid desorption and solid diffusion use the same matcher and safety checks.
    auto removal=base();Rule des;des.id="solid_desorption";des.type="desorption";des.before={"S"};des.beforeSolid={true};des.after={{"",false}};des.rate=1;removal.processes.push_back(des);
    Simulation eroded(removal);eroded.apply(1,16);check(eroded.solid().atomCount()==31 && eroded.catalogConsistent(),"solid reactant desorption");
    auto hopping=base();Rule hop;hop.id="solid_hop";hop.type="diffusion";hop.before={"S",""};hop.beforeSolid={true,false};hop.after={{"",false},{"S",true}};hop.counting="embeddings";Link vertical;hop.links={vertical};hopping.processes.push_back(hop);
    Simulation unsupported(hopping);try{unsupported.apply(1,16,32);throw std::logic_error("unsupported solid hop accepted");}catch(const std::runtime_error&){}
    hopping.processes[1].links[0]=shell;hopping.processes[1].links[0].sameHeight=false;
    Simulation moved(hopping);moved.apply(1,16,33);check(moved.solid().atomCount()==32 && !moved.solid().sites[16].occupied && moved.solid().sites[33].occupied,"simultaneous supported solid hop");
    Simulation adsorbed(removal);adsorbed.apply(0,32);try{adsorbed.apply(1,16);throw std::logic_error("adsorbate lost its support");}catch(const std::runtime_error&){}
    check(adsorbed.steps()==1 && adsorbed.solid().sites[16].occupied,"rejected removal partially applied");
    auto clock=base();clock.maxEvents=1;Simulation timed(clock);std::mt19937_64 rng(clock.seed);
    double expected=-std::log(std::generate_canonical<double,53>(rng))/32;timed.run();check(std::abs(timed.time()-expected)<1e-14,"graph KMC clock/counting");
    auto evolving=sets;evolving.maxEvents=20;evolving.end=10;Simulation run(evolving);double last=0,height=run.meanHeight(),integral=0;
    run.run([&](const Simulation& state,bool) {check(state.catalogConsistent(),"multi-site evolving catalog");integral+=(state.time()-last)*height;if(state.time()>0) check(std::abs(state.timeAverageHeight()-integral/state.time())<1e-10,"graph residence-time integral");last=state.time();height=state.meanHeight();});
    Simulation repeat(evolving);repeat.run();check(run.counts()==repeat.counts() && run.time()==repeat.time(),"graph seed reproducibility");
    auto limit=sets;limit.maxEmbeddings=10;try{Simulation big(limit);throw std::logic_error("embedding limit ignored");}catch(const std::runtime_error&){}
    auto bad=sets;bad.processes[1].links.pop_back();try{Simulation disconnected(bad);throw std::logic_error("disconnected pattern accepted");}catch(const std::runtime_error&){}
    auto small=diagonal;small.geometry.repetitions[0]=2;try{Simulation ambiguous(small);throw std::logic_error("periodic distance ambiguity accepted");}catch(const std::runtime_error&){}
    check(argc==2,"expected custom coordinate example");auto custom=Config::readJson(argv[1]);check(custom.surfaceKMC && custom.surfaceKMC->customGeometry,"custom geometry parser");
    Simulation bridge(*custom.surfaceKMC);try{bridge.apply(1,2);throw std::logic_error("bonded cluster detached from substrate");}catch(const std::runtime_error&){}
    check(bridge.solid().atomCount()==5,"failed bridge removal changed state");
    Simulation customRun(*custom.surfaceKMC);customRun.run();check(customRun.steps()>0,"custom coordinate process did not run");
    std::cout<<"All generalized geometry and matching tests passed.\n";return 0;
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
