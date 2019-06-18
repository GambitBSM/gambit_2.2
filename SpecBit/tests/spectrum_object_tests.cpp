//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Standalone executable to test basic Spectrum
///  and SpectrumContents classes 
///
///  Not built with GAMBIT cmake system. 
///  Compile with 'make -f MakeSpecBitTests' from
///  GAMBIT root directory. May require customisation
///  for your system.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (benjamin.farmer@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************

#include "gambit/Utils/standalone_error_handlers.hpp"
#include "gambit/Utils/standalone_utils.hpp"
#include "gambit/Utils/static_members.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Models/spectrum_contents.hpp"
#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"
#include "gambit/Elements/spectrum.hpp"
#include "gambit/Elements/slhaea_helpers.hpp"

using namespace Gambit;

int main(int argc, char* argv[])
{
    std::cout<<"Creating MSSM Spectrum Contents object"<<std::endl; 
    SpectrumContents::MSSM mssm_contents; 

    std::cout<<"Creating template SLHAea object corresponding to required spectrum contents: mssm_template.slha"<<std::endl;
    // Create template SLHAea object complying with MSSM contents definition
    // (Won't include stuff the Spectrum objects don't care about, like MODSEL blocks and whatnot)
    SLHAstruct mssm_slha = mssm_contents.create_template_SLHAea();

    std::cout<<"Writing SLHAea object to file for inspection: mssm_template.slha"<<std::endl;
    std::ofstream ofs("mssm_template.slha");
    ofs << mssm_slha;
    ofs.close();

    std::cout<<"Creating Spectrum object from template MSSM SLHAea object and corresponding Contents object"<<std::endl;
    Spectrum mssm_spec(mssm_slha,mssm_contents,100); // Last parameter is scale at which running parameters are defined. Could try to infer from certain blocks, but I think better to explictly specify it.

    std::cout<<"Writing SLHA-compliant outputs from MSSM Spectrum object"<<std::endl;
    mssm_spec.writeSLHAfile("mssm_compliant_template.slha1", 1);
    mssm_spec.writeSLHAfile("mssm_compliant_template.slha2", 2);

    // Do the same thing with a couple of test SLHA files from SoftSUSY
    std::cout<<"Reading SoftSUSY-generated SLHA1 spectrum and wrapping it"<<std::endl;
    SLHAstruct ss1;
    std::ifstream ifs1("SpecBit/tests/softsusy_output.slha1");
    ifs1 >> ss1;
    ifs1.close();
    double Q = SLHAea_get_scale(ss1,"GAUGE"); // Determine scale of spectrum
    mssm_spec = Spectrum(ss1,mssm_contents,Q);
    std::cout<<"Writing SLHA-compliant outputs from MSSM Spectrum object"<<std::endl;
    mssm_spec.writeSLHAfile("mssm_compliant_softsusy1.slha1", 1);
    mssm_spec.writeSLHAfile("mssm_compliant_softsusy1.slha2", 2);

    std::cout<<"Reading SoftSUSY-generated SLHA2 spectrum and wrapping it"<<std::endl;
    SLHAstruct ss2;
    std::ifstream ifs2("SpecBit/tests/softsusy_output.slha2");
    ifs2 >> ss2;
    ifs2.close();
    Q = SLHAea_get_scale(ss2,"GAUGE"); // Determine scale of spectrum
    mssm_spec = Spectrum(ss2,mssm_contents,Q);
    std::cout<<"Writing SLHA-compliant outputs from MSSM Spectrum object"<<std::endl;
    mssm_spec.writeSLHAfile("mssm_compliant_softsusy2.slha1", 1);
    mssm_spec.writeSLHAfile("mssm_compliant_softsusy2.slha2", 2);

    std::cout<<"Explicity call a few of the getters for good measure..."<<std::endl;
    for(int i=1;i<=3;i++)
    {
       std::cout<<"mssm_spec.get(Par::Pole_Mass,\"e-\","<<i<<"): "<<mssm_spec.get(Par::Pole_Mass,"e-",i)<<std::endl;
       std::cout<<"mssm_spec.get(Par::Pole_Mass,\"e+\","<<i<<"): "<<mssm_spec.get(Par::Pole_Mass,"e+",i)<<std::endl;
    }
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"e-\"): "<<mssm_spec.get(Par::Pole_Mass,"e-")<<std::endl;
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"mu-\"): "<<mssm_spec.get(Par::Pole_Mass,"mu-")<<std::endl;
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"tau-\"): "<<mssm_spec.get(Par::Pole_Mass,"tau-")<<std::endl;
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"e+\"): "<<mssm_spec.get(Par::Pole_Mass,"e+")<<std::endl;
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"mu+\"): "<<mssm_spec.get(Par::Pole_Mass,"mu+")<<std::endl;
    std::cout<<"mssm_spec.get(Par::Pole_Mass,\"tau+\"): "<<mssm_spec.get(Par::Pole_Mass,"tau+")<<std::endl;
    std::cout<<"mssm_spec.get(Par::dimensionless,\"sinW2\"): "<<mssm_spec.get(Par::dimensionless,"sinW2")<<std::endl;

    std::cout<<"End of tests!"<<std::endl;
}
