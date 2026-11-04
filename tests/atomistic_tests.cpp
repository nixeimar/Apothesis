#include "../src/config/configuration.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <chrono>

void check(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }
bool close(double a,double b) { return std::abs(a-b)<1e-9; }
int main(int argc, char** argv) {
    namespace fs = std::filesystem;
    fs::path file = fs::temp_directory_path() / ("apothesis-atoms-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".xyz");
    try {
        check(argc==2,"expected atomistic example path");
        auto config=Config::readJson(argv[1]);
        check(config.atomisticSurface.has_value(),"atomistic input not selected");
        auto example=Atomistic::build(*config.atomisticSurface);
        check(example.sites.size()==576 && example.atomCount()==576,"FCC atom count");
        check(close(example.cell[2][2],4*3.615+15),"vacuum not included in cell");
        check(close(example.sites[0].position[2],7.5),"vacuum must be split above and below slab");
        check(std::count_if(example.sites.begin(),example.sites.end(),[](const Atomistic::Site& s){return s.fixed;})==144,"fixed substrate count");
        example.writeXYZ(file.string());
        std::ifstream xyz(file);
        std::string line;
        std::getline(xyz,line); check(line=="576","XYZ atom count");
        std::getline(xyz,line); check(line.find("pbc=\"T T F\"")!=std::string::npos,"XYZ periodic metadata");
        int count=0; while(std::getline(xyz,line)) ++count;
        check(count==576,"XYZ rows");
        fs::remove(file);
        Atomistic::Options options;
        options.latticeConstant=2; options.species="Cu"; options.element="Cu";
        options.repetitions={{4,4,4}}; options.periodic={{true,true,true}};
        for (const auto& item : {std::make_pair("SimpleCubic",6),{"BCC",8},{"FCC",12},{"Diamond",4},{"HCP",12}}) {
            options.crystal=item.first;
            auto surface=Atomistic::build(options);
            for (const auto& site : surface.sites) {
                check(site.neighbors.size()==static_cast<std::size_t>(item.second),options.crystal+": wrong bulk coordination at "+std::to_string(site.id));
                check(surface.occupiedCoordination(site.id)==site.neighbors.size(),"occupied coordination mismatch");
                for (auto j : site.neighbors) {
                    check(j!=site.id,"self neighbor");
                    check(std::binary_search(surface.sites[j].neighbors.begin(),surface.sites[j].neighbors.end(),site.id),"asymmetric neighbor graph");
                }
            }
        }
        options.crystal="SimpleCubic";
        options.periodic={{true,true,false}};
        auto slab=Atomistic::build(options);
        check(slab.sites[0].neighbors.size()==5,"open bottom boundary");
        check(slab.sites[16].neighbors.size()==6,"interior coordination");
        auto expected=slab.sites[0].position;
        options.origin={{-5,4,-3}};
        auto moved=Atomistic::build(options);
        for (int a=0;a<3;++a) check(close(moved.sites[0].position[a],expected[a]+options.origin[a]),"origin translation");
        for (std::size_t i=0;i<slab.sites.size();++i) check(slab.sites[i].neighbors==moved.sites[i].neighbors,"translation changed topology");
        options.vacancies={16};
        auto defective=Atomistic::build(options);
        check(defective.atomCount()==63 && defective.sites.size()==64,"vacancies must preserve candidate sites and IDs");
        check(defective.occupiedCoordination(0)==4,"vacancy not reflected in coordination");
        auto invalid=[&](const Atomistic::Options& bad,const std::string& field) {
            try { Atomistic::validate(bad); }
            catch(const std::runtime_error& e) {check(std::string(e.what()).find(field)!=std::string::npos,"unexpected diagnostic: "+std::string(e.what()));return;}
            throw std::runtime_error("invalid input accepted: "+field);
        };
        auto bad=options; bad.vacancies={1000};invalid(bad,"vacancies");
        bad=options;bad.vacancies={1,1};invalid(bad,"vacancies");
        bad=options;bad.vacancies.clear();bad.repetitions={{1,1,1}};invalid(bad,"periodic cell");
        bad=options;bad.repetitions={{1000000,1000000,1000000}};invalid(bad,"limit");
        bad=options;bad.periodic[2]=true;bad.vacuum=10;invalid(bad,"vacuum");
        bad=options;bad.latticeConstant=0;invalid(bad,"lattice_constant");
        std::cout<<"All atomistic surface tests passed.\n";
        return 0;
    } catch(const std::exception& e) { fs::remove(file);std::cerr<<e.what()<<'\n';return 1; }
}
