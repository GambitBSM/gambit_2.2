"""
Master module for all DarkBit related routines.
"""

import numpy as np

from setup import *
from files import *

def sort_annihilations(dm, three_fields, four_fields):
    """
    Sorts BSM vertices, returns DM+DM -> X + Y
    """

    if not isinstance(dm, Particle):
        GumError("\n\nDM not passed over as an instance of class Particle.")

    dm_dm_to_x_y = []

    # Add all 4-pt vertices that are DM + DM -> X + Y
    for i in range(0, len(four_fields)):

        if dm.is_sc and four_fields[i].count(dm.PDG_code) == 2:
            products = [f for f in four_fields[i] if f not in {dm.PDG_code}]
            dm_dm_to_x_y.append(products)

        elif not dm.is_sc and four_fields[i].count(dm.PDG_code) == 1 and \
                four_fields[i].count(dm.Conjugate.PDG_code) == 1:
            products = [f for f in four_fields[i] if
                        f not in {dm.PDG_code, dm.Conjugate.PDG_code}]
            dm_dm_to_x_y.append(products)

    v_with_dm = []
    v_with_dmbar = []
    v_with_2dm = []
    v_with_dmdmbar = []
    s_channel_propagators = []

    # Separate all vertices
    for i in range(0, len(three_fields)):

        # Self conjugate DM.
        if dm.is_sc():

            if three_fields[i].count(dm.PDG_code) == 1:
                v_with_dm.append(three_fields[i])

            elif three_fields[i].count(dm.PDG_code) == 2:
                v_with_2dm.append(three_fields[i])
                # Add s-channel propagator: DM+DM -> prop
                product = [f for f in three_fields[i] if f != dm.PDG_code]
                s_channel_propagators.append(product[0])

        # Not self-conjugate DM.
        elif not dm.is_sc():

            if three_fields[i].count(dm.PDG_code) == 1:

                if three_fields[i].count(dm.Conjugate.PDG_code) == 0:
                    v_with_dm.append(three_fields[i])

                elif three_fields[i].count(dm.Conjugate.PDG_code) == 1:
                    v_with_dmdmbar.append(three_fields[i])
                    # Add s-channel propagator: DM+DMbar -> prop
                    product = [f for f in three_fields[i] if
                               abs(f) != dm.PDG_code]
                    s_channel_propagators.append(product[0])

            elif three_fields[i].count(dm.Conjugate.PDG_code) == 1:
                v_with_dmbar.append(three_fields[i])

    # s-channel: sticking together
    # DM + DMbar -> propagator & propagator -> stuff + otherstuff

    # See what each propagator can decay into.
    for i in range(0, len(s_channel_propagators)):

        propagator_pdg = s_channel_propagators[i]

        for j in range(0, len(three_fields)):

            curr_vertex = three_fields[j][:]

            if curr_vertex.count(propagator_pdg) >= 1:
                # Copy the vertex & remove one instance of the propagator.
                curr_vertex.remove(
                    curr_vertex[curr_vertex.index(propagator_pdg)])
                dm_dm_to_x_y.append(curr_vertex)

    # t/u-channel: sticking together
    # DM + stuff -> propagator & propagator -> DMbar + otherstuff

    t_channel_propagators = []
    t_channel_potential = []

    # self-conj:
    if dm.is_sc():

        dm_verts = v_with_dm + v_with_2dm
        verts = dm_verts[:]

        # Remove one instance of DM from each vertex -- the incoming DM.
        for i in range(0, len(verts)):
            verts[i].remove(verts[i][verts[i].index(dm.PDG_code)])
            t_channel_potential.append(verts[i])

        vertsbar = verts

    elif not dm.is_sc():

        verts = v_with_dm[:]
        vertsbar = v_with_dmbar[:]

        # Remove one instance of DM from each vertex -- the incoming DM.
        for i in range(0, len(verts)):
            verts[i].remove(verts[i][verts[i].index(dm.PDG_code)])
            t_channel_potential.append(verts[i])

        # Remove DMbar from each vertex too -- also incoming.
        for i in range(0, len(vertsbar)):
            vertsbar[i].remove(
                vertsbar[i][vertsbar[i].index(dm.Conjugate.PDG_code)])

    # Remove any duplicate t-channel props.
    t_channel_propagators = list(set(t_channel_propagators))

    for i in range(0, len(t_channel_potential)):

        for j in range(0, len(vertsbar)):

            curr_vertex = vertsbar[j][:]

            if dm.spinX2 != 1:
                for k in range(2):
                    # If the mediator shows up in the DMbar ->
                    # add the two products and the mediator to respective lists.
                    if t_channel_potential[i][k] in curr_vertex:
                        ind = curr_vertex.index(t_channel_potential[i][k])
                        dm_dm_to_x_y.append([t_channel_potential[i][k - 1],
                                             curr_vertex[ind - 1]])
                        t_channel_propagators.append(t_channel_potential[i][k])

            # If DM is fermionic -> also allowed fermionic t-channel prop
            # TODO -- this requires continuous fermion chain
            else:
                pass

    propagators = s_channel_propagators + t_channel_propagators
    propagators[:] = list(set(propagators))

    # Remove all duplicates within the list of annihilation products.
    ann_products = map(list,
                       sorted(set(map(tuple, dm_dm_to_x_y)), reverse=True))
    # Remove any DM-DM self-interactions
    if [dm.PDG_code, dm.Conjugate.PDG_code] in ann_products:
        ann_products.remove([dm.PDG_code, dm.Conjugate.PDG_code])
    return np.array(ann_products), propagators

    
def xsecs(dm, ann_products, gambit_pdg_dict, gambit_model_name,
          calchep_pdg_dict):
    """
    Writes all entries for <sigma v> within the Process Catalogue, 
    utilising CalcHEP.
    """
    
    # DM mass parameter
    gb_id = pdg_to_particle(dm.PDG_code, gambit_pdg_dict)
    dm_mass = "m" + gb_id.replace("\"", "")
    
    towrite_class = (
            "double sv(std::vector<str> channel, double v_rel)\n"
            "{\n"
            "/// Returns sigma*v for a given channel.\n"
            "double GeV2tocm3s1 = gev2cm2*s2cm;\n\n"
    )

    out1g = np.array([pdg_to_particle(x, gambit_pdg_dict) for x in ann_products[:,0]])
    out2g = np.array([pdg_to_particle(x, gambit_pdg_dict) for x in ann_products[:,1]])
    dm_chep = pdg_to_particle(dm.PDG_code, calchep_pdg_dict)
    dm_chepc = pdg_to_particle(dm.Conjugate.PDG_code, calchep_pdg_dict)

    # Add each channel individually for annihilation cross sections
    for i in np.arange(len(ann_products)):
        towrite_class += (
                "if (channel == '{0}, {1}')"
                " return BEreq::CH_Sigma_V({2},"
                " std::vector<str> {{'{3}', '{4}'}},"
                " std::vector<str> {{'{0}', '{1}'}},"
                " QCD_coupling, v_rel, tbl)*GeV2tocm3s1; \n"
        ).format(out1g[i], out2g[i], gambit_model_name, dm_chep, dm_chepc)

    towrite_class += (
            "else return 0;\n"
            "}\n\n"
    )
    
    channels = ', '.join("\'{}, {}'".format(*t) for t in zip(out1g, out2g))
    p1 = ', '.join("\'{}'".format(*t) for t in zip(out1g))
    p2 = ', '.join("\'{}'".format(*t) for t in zip(out2g))

    towrite_pc = (
            "// Instantiate new {0} object.\n"
            "auto pc = boost::make_shared<{0}>(&catalog, &tbl);\n"
            "\n"
            "// Populate annihilation channel list and add "
            "thresholds to threshold list.\n"
            "process_ann.resonances_thresholds.threshold_energy"
            ".push_back(2*{1});\n"
            "auto channels = \n"
            "  daFunk::vec<string>({2});\n"
            "auto p1 = \n"
            "  daFunk::vec<string>({3});\n"
            "auto p2 = \n"
            "  daFunk::vec<string>({4});\n"
            "\n"
            "for (unsigned int i = 0; i < channels.size(); ++i)\n"
            "{{\n"
            "double mtot_final = \n"
            "catalog.getParticleProperty(p1.[i]).mass + \n"
            "catalog.getParticleProperty(p2.[i]).mass;  \n"
            "if ({1}*2 > mtot_final*0.5)\n"
            "{{\n"
            "daFunk::Funk kinematicFunction = daFunk::funcM("
            "pc, &{0}::sv, channel[i], daFunk::var('v'));\n"
            "TH_Channel new_channel("
            "daFunk::vec<string>(p1[i], p2[i]), kinematicFunction);\n"
            "process_ann.channelList.push_back(new_channel);\n"
            "}}\n"
            "if ({1}*2 > mtot_final)\n"
            "{{\n"
            "process_ann.resonances_thresholds.threshold_energy.\n"
            "push_back(mtot_final);\n"
            "}}\n"
            "}}\n"
            "\n"
    ).format(gambit_model_name, dm_mass, channels, p1, p2)       
    
    return towrite_class, towrite_pc      

def proc_cat(dm, sv, ann_products, propagators, gambit_pdg_dict, 
             gambit_model_name, calchep_pdg_dict, model_specific_particles, 
             exclude_decays):
    """
    Writes all entries for the Process Catalogue for DarkBit.
    """

    gb_id = pdg_to_particle(dm.PDG_code, gambit_pdg_dict)
    gb_conj = pdg_to_particle(dm.Conjugate.PDG_code, gambit_pdg_dict)

    towrite = (
            "class {0}\n"
            "{{\n"
            "public:\n"
            "/// Initialize {0} object (branching ratios etc)\n"
            "{0}(TH_ProcessCatalog* const catalog, DecayTable* const tbl);\n"
            "~{0}();\n\n"
    ).format(gambit_model_name)
    
    if sv:
        sv_class, sv_src = xsecs(dm, ann_products, gambit_pdg_dict, 
                                 gambit_model_name, calchep_pdg_dict)
        towrite += sv_class

    towrite += (
            "\n"
            "}};\n\n"
            "void TH_ProcessCatalog_{0}(DarkBit::TH_ProcessCatalog &result)\n"
            "{{\n"
            "using namespace Pipes::TH_ProcessCatalog_{0};\n"
            "using std::vector;\n"
            "using std::string;\n\n"
            "// Initialize empty catalog, main annihilation process\n"
            "TH_ProcessCatalog catalog;\n"
            "TH_Process process_ann(\"{1}\", \"{2}\");"
            "\n"
    ).format(gambit_model_name, gb_id, gb_conj)

    # Add flag for (non-)self-conjugate DM to rescale spectra properly
    if not dm.is_sc():
        towrite += (
                "\n"
                "// Explicitly state that Dirac DM is not self-conjugate to add"
                " extra \n// factors of 1/2 where necessary\n"
                "process_ann.isSelfConj = false;\n\n"
        )
    else:
        towrite += (
                "\n"
                "// Explicitly state that DM is self-conjugate\n"
                "process_ann.isSelfConj = true;\n\n"
        )

    towrite += add_SM_macros()

    # Add the new BSM particles to the Process Catalog
    dm_mass = "m" + gb_id.replace("\"", "")

    towrite += (
            "// {0}-specific masses\n"
            "double {1} = spec.get(Par::Pole_Mass, '{2}'));\n"
            "addParticle('{2}', {1}, {3});\n"
    ).format(gambit_model_name, dm_mass, gb_id, dm.spinX2)

    for i in np.arange(len(model_specific_particles)):
        if model_specific_particles[i].PDG_code != dm.PDG_code:
          towrite += (
                  "addParticle('{0}', spec.get(Par::Pole_Mass, '{1}'), {2});\n"
          ).format(pdg_to_particle(model_specific_particles[i].PDG_code, gambit_pdg_dict),
                   pdg_to_particle(model_specific_particles[i].PDG_code, gambit_pdg_dict),
                   str(model_specific_particles[i].spinX2))
                   
    towrite += (
            "\n"
            "// Get rid of convenience macros\n"
            "#undef getSMmass\n"
            "#undef addParticle\n"
            "\n"
            "// Import decay table from DecayBit\n"
            "const DecayTable* tbl = &(*Dep::decay_rates);\n"
            "\n"
            "// Set of imported decays\n"
            "std::set<string> importedDecays;\n"
            "\n"
            "// Minimum branching ratio to include\n"
            "double minBranching = runOptions->getValueOrDef<double>(0.0,"
            " \"ProcessCatalog_MinBranching\");\n"
            "\n"
            "// Import relevant decays\n"
            "using DarkBit_utils::ImportDecays;\n"
            "\n"
            "      *** TODO: excludeDecays?? ***\n"
            "      *** And ImportDecays!! For all propagators? ***\n"
            "      *** and all final states?? Which to exclude...?\n"
            "\n"
    )

    if sv:
        towrite += sv_src

    for i in np.arange(len(propagators)):
        if abs(propagators[i]) != abs(dm.PDG_code):
            towrite += (
                    "if (spec.get(Par::Pole_Mass, '{0}') >= 2*{1}) "
                    "process_ann.resonances_thresholds.resonances.\n    "
                    "push_back(TH_Resonance((spec.get(Par::Pole_Mass, '{0}'), "
                    "tbl.at('{0}').width_in_GeV)));\n"
            ).format(pdg_to_particle(propagators[i], gambit_pdg_dict),
                 dm_mass)

    towrite += (
            "\n"
            "catalog.processList.push_back(process_ann);\n\n"
            "// Validate\n"
            "catalog.validate();\n\n"
            "result = catalog;\n"
            "} // function TH_ProcessCatalog\n"
    )

    return towrite


def write_darkbit_src(dm, pc, sv, dd, ann_products, propagators, 
                      gambit_pdg_dict, gambit_model_name, calchep_pdg_dict,
                      model_specific_particles, exclude_decays = []):
    """
    Collects all source for DarkBit: process catalogue, direct detection...
    """

    gb_id = pdg_to_particle(dm.PDG_code, gambit_pdg_dict)
    gb_conj = pdg_to_particle(dm.Conjugate.PDG_code, gambit_pdg_dict)

    if not isinstance(dm, Particle):
        print("DM not passed over as an instance of class Particle.")
        exit()

    intro_message = (
            "///  Implementation of {0}\n"
            "///  DarkBit routines."
    ).format(gambit_model_name)

    towrite = blame_gum(intro_message)

    towrite += (
            "#include \"gambit/Elements/gambit_module_headers.hpp\"\n"
            "#include \"gambit/DarkBit/DarkBit_rollcall.hpp\"\n"
            "#include \"gambit/Utils/ascii_table_reader.hpp\"\n"
            "#include \"boost/make_shared.hpp\"\n"
            "#include \"gambit/DarkBit/DarkBit_utils.hpp\"\n"
            "\n"
            "namespace Gambit\n"
            "{\n"
            "namespace DarkBit\n"
            "{\n"
    )

    if pc:
        towrite += proc_cat(dm, sv, ann_products, propagators,
                            gambit_pdg_dict, gambit_model_name,
                            calchep_pdg_dict, model_specific_particles,
                            exclude_decays)

    if dd:
        towrite += write_direct_detection(gambit_model_name)

    towrite += write_dm_id(gambit_model_name, gb_id)

    towrite += (
            "} //namespace DarkBit\n\n"
            "} //namespace Gambit\n"
            "\n"
    )

    return indent(towrite)

def write_dm_id(model_name, dm_id):
    """
    Returns entry for DarkMatter_ID in DarkBit.
    """

    towrite = (
            "\n"
            "void DarkMatter_ID_{0}(std::string& result)"
            "{{ result = \"{1}\"; }}"
            "\n\n"
    ).format(model_name, dm_id)

    return towrite;

def add_SM_macros():
    """
    Adds Standard Model macros to the Process Catalogue.
    """

    towrite = (
            "\n"
            "// Import particle masses \n"
            "\n"
            "// Convenience macros\n"
            "#define getSMmass(Name, spinX2) "
            "catalog.particleProperties.insert(std::pair<string, "
            "TH_ParticleProperty> (Name, TH_ParticleProperty"
            "(SM.get(Par::Pole_Mass,Name), spinX2)));\n"
            "#define addParticle(Name, Mass, spinX2) "
            "catalog.particleProperties.insert(std::pair<string, "
            "TH_ParticleProperty> (Name, "
            "TH_ParticleProperty(Mass, spinX2)));\n"
            "\n"
            "// Import Spectrum objects\n"
            "const Spectrum& spec = *Dep::SingletDM_spectrum;\n"
            "const SubSpectrum& SM = spec.get_LE();\n"
            "\n"
            "// Get SM pole masses\n"
            "getSMmass(\"e-_1\",     1)\n"
            "getSMmass(\"e+_1\",     1)\n"
            "getSMmass(\"e-_2\",     1)\n"
            "getSMmass(\"e+_2\",     1)\n"
            "getSMmass(\"e-_3\",     1)\n"
            "getSMmass(\"e+_3\",     1)\n"
            "getSMmass(\"Z0\",       2)\n"
            "getSMmass(\"W+\",       2)\n"
            "getSMmass(\"W-\",       2)\n"
            "getSMmass(\"g\",        2)\n"
            "getSMmass(\"gamma\",    2)\n"
            "getSMmass(\"u_3\",      1)\n"
            "getSMmass(\"ubar_3\",   1)\n"
            "getSMmass(\"d_3\",      1)\n"
            "getSMmass(\"dbar_3\",   1)\n"
            "\n"
            "// Pole masses not available for the light quarks.\n"
            "addParticle(\"u_1\"   , SMI.mU,  1) // mu(2 GeV)^MS-bar\n"
            "addParticle(\"ubar_1\", SMI.mU,  1) // mu(2 GeV)^MS-bar\n"
            "addParticle(\"d_1\"   , SMI.mD,  1) // md(2 GeV)^MS-bar\n"
            "addParticle(\"dbar_1\", SMI.mD,  1) // md(2 GeV)^MS-bar\n"
            "addParticle(\"u_2\"   , SMI.mCmC,1) // mc(mc)^MS-bar\n"
            "addParticle(\"ubar_2\", SMI.mCmC,1) // mc(mc)^MS-bar\n"
            "addParticle(\"d_2\"   , SMI.mS,  1) // ms(2 GeV)^MS-bar\n"
            "addParticle(\"dbar_2\", SMI.mS,  1) // ms(2 GeV)^MS-bar\n"
            "\n"
            "// Masses for neutrino flavour eigenstates. Set to zero.\n"
            "// (presently not required)\n"
            "addParticle(\"nu_e\",     0.0, 1)\n"
            "addParticle(\"nubar_e\",  0.0, 1)\n"
            "addParticle(\"nu_mu\",    0.0, 1)\n"
            "addParticle(\"nubar_mu\", 0.0, 1)\n"
            "addParticle(\"nu_tau\",   0.0, 1)\n"
            "addParticle(\"nubar_tau\",0.0, 1)\n"
            "\n"
            "// Meson masses\n"
            "addParticle(\"pi0\",   meson_masses.pi0,       0)\n"
            "addParticle(\"pi+\",   meson_masses.pi_plus,   0)\n"
            "addParticle(\"pi-\",   meson_masses.pi_minus,  0)\n"
            "addParticle(\"eta\",   meson_masses.eta,       0)\n"
            "addParticle(\"rho0\",  meson_masses.rho0,      1)\n"
            "addParticle(\"rho+\",  meson_masses.rho_plus,  1)\n"
            "addParticle(\"rho-\",  meson_masses.rho_minus, 1)\n"
            "addParticle(\"omega\", meson_masses.omega,     1)\n"
            "\n"
    )

    return towrite

def write_direct_detection(model_name):
    """
    Writes direct detection bits in DarkBit... TODO
    """

    towrite = (
            "\n"
            "/// Direct detection couplings.\n"
            "void DD_couplings_{0}(DM_nucleon_couplings & result)\n"
            "{{\n"
            "using namespace Pipes::DD_couplings_{0};\n"
            "const Spectrum& spec = *Dep::{0}_spectrum;\n"
            "\n"
            "  *** TODO *** \n"
            "result.gps = ();\n"
            "result.gns = ();\n"
            "result.gpa = ();\n"
            "result.gna = ();\n"
            "\n"
            "logger() << LogTags::debug << '{0} DD couplings:' << std::endl;\n"
            "logger() << ' gps = ' << result.gps << std::endl;\n"
            "logger() << ' gns = ' << result.gns << std::endl;\n"
            "logger() << ' gpa = ' << result.gpa << std::endl;\n"
            "logger() << ' gna = ' << result.gna << EOM;\n"
            "\n"
            "}}\n"
            "\n"
    ).format(model_name)

    return towrite

def write_darkbit_rollcall(model_name, pc, dd):
    """
    Writes the rollcall header entries for new DarkBit entry.
    """

    if pc:
        pro_cat = dumb_indent(4, (
                "#define FUNCTION TH_ProcessCatalog_{0}\n"
                "  START_FUNCTION(DarkBit::TH_ProcessCatalog)\n"
                "  DEPENDENCY(decay_rates, DecayTable)\n"
                "  DEPENDENCY({0}_spectrum, Spectrum)\n"
                "  BACKEND_REQ(CH_Sigma_V, (), double, (str&, std::vector<str>&, "
                "std::vector<str>&, double&, double&, const DecayTable&))\n"
                "  ALLOW_MODELS({0})\n"
                "#undef FUNCTION\n"
        ).format(model_name))
    else:
        pro_cat = None

    if dd:
        dir_det = dumb_indent(4, (
                "#define FUNCTION DD_couplings_{0}\n"
                "START_FUNCTION(DM_nucleon_couplings)\n"
                "DEPENDENCY({0}_spectrum, Spectrum)\n"
                "*** TODO *** \n"
                "#undef FUNCTION\n"
        ).format(model_name))
    else:
        dir_det = None

    dm_id = dumb_indent(4, (
            "#define FUNCTION DarkMatter_ID_{0}\n"
            "START_FUNCTION(std::string)\n"
            "ALLOW_MODELS({0})\n"
            "#undef FUNCTION\n"
    ).format(model_name))

    return pro_cat, dir_det, dm_id


