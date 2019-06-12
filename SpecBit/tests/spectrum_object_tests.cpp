//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Standalone executable to test basic Spectrum
///  and SpectrumContents classes 
///
///  Not built with GAMBIT cmake system. 
///  Compile like this:
///
///  > 
//    -IModels/include
//    -IElements/include
//    -IUtils/include
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

#include "gambit/Models/SpectrumContents/spectrum_contents.hpp"
#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"
#include "gambit/Elements/spectrum.hpp"

using namespace Gambit;

int main(int argc, char* argv[])
{
    // Create MSSM spectrum contents object
    SpectrumContents::MSSM mssm(); 

    // Create template SLHA file complying with MSSM contents definition
    // (Won't include stuff the Spectrum objects don't care about, like MODSEL blocks and whatnot)
    mssm.create_template_SLHA_file("mssm_template.slha");
}
