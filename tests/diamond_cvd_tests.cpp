#include "../src/config/configuration.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <algorithm>
using namespace DiamondCVD;
void check(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
std::size_t find(const Simulation& s,const std::string& id) {
    for(std::size_t i=0;i<s.events().size();++i) if(s.options().processes[s.events()[i].process].id==id) return i;
    return std::numeric_limits<std::size_t>::max();
}
void execute(Simulation& s,const std::string& id) {auto i=find(s,id);check(i!=std::numeric_limits<std::size_t>::max(),("missing event "+id).c_str());s.execute(i);check(s.consistent(),"bond/valence consistency after forced event");}
int main(int argc,char** argv) {
 try {
    check(argc==2,"expected full diamond input");auto c=Config::readJson(argv[1]);check(c.diamondCVD.has_value(),"diamond dispatch");auto o=*c.diamondCVD;
    check(o.processes.size()==22,"all twelve Table II rows must produce 22 directional channels");
    Simulation initial(o);check(initial.count(Kind::Diamond)==40 && initial.count(Kind::Hydrogen)==8,"(111) substrate/passivation counts");check(initial.events().size()==16,"initial surface allows only H abstraction and desorption");
    const double expected=1.3e14*std::exp(-7.3/(.00198720425864083*1200))*(.04/760/(82.057366080960*1200));
    check(std::abs(o.processes[0].rate/expected-1)<1e-12,"Table II forward units/gas concentration");
    const double reverse=1.016e-5*1e13*std::exp(32.8/1.98720425864083)*std::exp(-96.9/(.00198720425864083*1200));
    check(std::abs(o.processes[3].rate/reverse-1)<1e-12,"reverse thermochemistry sign/standard concentration");
    execute(initial,"r1_forward");
    for(const auto& id:{"r3_reverse","r4_forward","r5_forward"}) {
        double sum=0;for(const auto& e:initial.events()) if(o.processes[e.process].id==id) sum+=e.rate;
        auto p=std::find_if(o.processes.begin(),o.processes.end(),[&](const Process& p){return p.id==id;});check(p!=o.processes.end() && std::abs(sum/p->rate-1)<1e-12,"adsorption geometry multiplicity incorrectly multiplies a site rate");
    }
    auto methyl=initial;execute(methyl,"r4_forward");check(methyl.count(Kind::Hydrocarbon)==1 && methyl.count(Kind::Hydrogen)==10,"actual methyl atom inventory");
    execute(methyl,"r6_forward");check(methyl.count(Kind::Hydrogen)==9,"hydrocarbon hydrogen abstraction");execute(methyl,"r7_forward");check(methyl.count(Kind::Hydrogen)==10,"hydrocarbon repassivation");execute(methyl,"r4_reverse");check(methyl.count(Kind::Hydrocarbon)==0 && methyl.count(Kind::Hydrogen)==7,"methyl desorption inventory");
    auto acetylene=initial;execute(acetylene,"r5_forward");check(acetylene.count(Kind::Hydrocarbon)==2 && acetylene.count(Kind::Hydrogen)==9,"C2H2 atomic representation");execute(acetylene,"r8_forward");execute(acetylene,"r9_forward");execute(acetylene,"r5_reverse");check(acetylene.count(Kind::Hydrocarbon)==0,"acetylene desorption");
    // Every enumerated event in representative molecular states must be a
    // complete valid transition, not merely a species-name match.
    auto molecule=initial;execute(molecule,"r4_forward");for(std::size_t i=0;i<molecule.events().size();++i) {auto copy=molecule;copy.execute(i);check(copy.consistent(),"invalid enumerated molecular transition");}
    // Figure 2(c): an acetylene plus a methyl complete a three-carbon bridge
    // on (111). Force only the preparatory abstraction steps, then execute the
    // actual six-ring incorporation event from the catalog.
    auto bareOptions=o;bareOptions.passivate=false;Simulation bare(bareOptions);bool bridge=false;
    for(std::size_t first=0;first<bare.events().size() && !bridge;++first) if(o.processes[bare.events()[first].process].id=="r5_forward") {
        auto c2=bare;c2.execute(first);
        while(find(c2,"r8_forward")!=std::numeric_limits<std::size_t>::max()) c2.execute(find(c2,"r8_forward"));
        for(std::size_t second=0;second<c2.events().size() && !bridge;++second) if(o.processes[c2.events()[second].process].id=="r4_forward") {
            auto combined=c2;combined.execute(second);
            while(find(combined,"r6_forward")!=std::numeric_limits<std::size_t>::max()) combined.execute(find(combined,"r6_forward"));
            for(std::size_t k=0;k<combined.events().size();++k) if(o.processes[combined.events()[k].process].id=="r12_forward" && combined.events()[k].carbons.empty() && combined.events()[k].diamond.size()==3) {
                combined.execute(k);check(combined.consistent() && combined.count(Kind::Diamond)==43 && combined.count(Kind::Hydrocarbon)==0,"(111) three-carbon ring incorporation");bridge=true;break;
            }
        }
    }
    check(bridge,"CH3/C2H2 three-carbon bridge was not exercised");
    auto dimer=o;Atomistic::Options geometry;geometry.crystal="Diamond";geometry.species="C";geometry.element="C";geometry.latticeConstant=3.567;geometry.repetitions={{4,4,4}};geometry.occupiedUnitCells=2;geometry.fixedUnitCells=1;dimer.lattice=Atomistic::build(geometry);dimer.passivate=false;
    Simulation reconstructed(dimer);execute(reconstructed,"r11_reverse");check(reconstructed.dimers().size()==1,"(100) dimer closing");check(find(reconstructed,"r11_forward")==std::numeric_limits<std::size_t>::max(),"dimer opening without adjacent radical hydrocarbon");
    bool insertion=false;
    for(std::size_t i=0;i<reconstructed.events().size() && !insertion;++i) if(dimer.processes[reconstructed.events()[i].process].id=="r4_forward") {
        auto trial=reconstructed;trial.execute(i);
        auto abs=find(trial,"r6_forward");if(abs==std::numeric_limits<std::size_t>::max()) continue;trial.execute(abs);
        auto opening=find(trial,"r11_forward");if(opening==std::numeric_limits<std::size_t>::max()) continue;trial.execute(opening);check(trial.dimers().empty(),"dimer opening event");
        for(std::size_t k=0;k<trial.events().size();++k) if(dimer.processes[trial.events()[k].process].id=="r12_forward" && !trial.events()[k].carbons.empty()) {
            auto grown=trial;auto C=grown.count(Kind::Diamond)+grown.count(Kind::Hydrocarbon),H=grown.count(Kind::Hydrogen);grown.execute(k);check(grown.consistent(),"dimer insertion valence/connectivity");check(grown.count(Kind::Diamond)+grown.count(Kind::Hydrocarbon)==C && grown.count(Kind::Hydrogen)==H,"dimer insertion must conserve C/H atoms");check(grown.count(Kind::Diamond)>trial.count(Kind::Diamond),"dimer insertion did not incorporate carbon");insertion=true;break;
        }
    }
    check(insertion,"dimer insertion/rotation path was not exercised");
    // Direct clock and atomic material balances over a stochastic trajectory.
    auto runOptions=o;runOptions.maxEvents=2000;runOptions.end=1;Simulation run(runOptions);double last=0,height=run.meanHeight(),integral=0;
    const auto initialC=run.count(Kind::Diamond),initialH=run.count(Kind::Hydrogen);
    run.run([&](const Simulation& s,bool) {
        check(s.consistent(),"trajectory contains broken valence or unsupported carbon");
        integral+=(s.time()-last)*height;if(s.time()>0) check(std::abs(s.timeAverageHeight()-integral/s.time())<1e-9,"diamond residence time integral");last=s.time();height=s.meanHeight();
        auto count=[&](const std::string& id){for(std::size_t i=0;i<s.options().processes.size();++i) if(s.options().processes[i].id==id) return static_cast<long long>(s.counts()[i]);return 0LL;};
        const auto C=static_cast<long long>(initialC)+count("r3_reverse")+count("r4_forward")+2*count("r5_forward")+count("r10_forward")-count("r3_forward")-count("r4_reverse")-2*count("r5_reverse");
        check(C==static_cast<long long>(s.count(Kind::Diamond)+s.count(Kind::Hydrocarbon)),"stochastic carbon balance");
        const auto H=static_cast<long long>(initialH)-count("r1_forward")+count("r1_reverse")+count("r2_forward")-count("r2_reverse")-2*count("r3_forward")+2*count("r3_reverse")+3*count("r4_forward")-3*count("r4_reverse")+2*count("r5_forward")-2*count("r5_reverse")-count("r6_forward")+count("r6_reverse")+count("r7_forward")-count("r7_reverse")-count("r8_forward")+count("r8_reverse")+count("r9_forward")-count("r9_reverse")+3*count("r10_forward");
        check(H==static_cast<long long>(s.count(Kind::Hydrogen)),"stochastic hydrogen balance");
    });check(run.steps()>10,"no diamond kinetics exercised");
    Simulation repeated(runOptions);repeated.run();check(run.time()==repeated.time() && run.counts()==repeated.counts(),"diamond seed reproducibility");
    auto clock=o;clock.maxEvents=1;Simulation timed(clock);double total=0;for(const auto& e:timed.events()) total+=e.rate;std::mt19937_64 rng(o.seed);double dt=-std::log(std::generate_canonical<double,53>(rng))/total;timed.run();check(std::abs(timed.time()-dt)<1e-14,"diamond pre-event KMC clock");
    auto stiff=o;stiff.start=1;stiff.end=2;stiff.maxEvents=1000;
    for(auto& p:stiff.processes) p.rate=(p.id=="r1_forward" || p.id=="r2_forward") ? 1e20 : 0;
    Simulation compensated(stiff);compensated.run();check(compensated.steps()==1000 && compensated.stopReason()=="event_limit" && compensated.consistent(),"stiff compensated clock lost or rejected valid fast events");
    std::cout<<"All diamond chemistry, dimer and conservation tests passed.\n";return 0;
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
