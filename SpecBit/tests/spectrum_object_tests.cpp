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
#include "gambit/Models/SpectrumContents/spectrum_contents.hpp"
#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"
#include "gambit/Elements/spectrum.hpp"

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
    Spectrum mssm_spec(mssm_slha,mssm_contents);
    std::cout<<"End of tests!"<<std::endl;
}
