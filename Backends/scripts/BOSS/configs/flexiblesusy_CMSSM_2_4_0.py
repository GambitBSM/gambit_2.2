###################################
#                                 #
#  Configuration module for BOSS  #
#                                 #
###################################


# ~~~~~ CASTXML options ~~~~~

# See CastXML documentation for details on these options:
#
#   https://github.com/CastXML/CastXML/blob/master/doc/manual/castxml.1.rst
#

#
# *** Special note for OS X ***
#
# BOSS will most likely fail if 'g++' points to the Clang compiler.
# Install GNU g++ and point the castxml_cc variable below the GNU
# g++ executable.
#

castxml_cc_id  = 'gnu'         # Reference compiler: 'gnu', 'gnu-c', 'msvc', 'msvc-c'
castxml_cc     = 'g++'         # Name a specific compiler: 'g++', 'cl', ...
castxml_cc_opt = '-std=c++11'  # Additional option string passed to the compiler in castxml_cc (e.g. '-m32')


# ~~~~~ GAMBIT-specific options ~~~~~

gambit_backend_name    = 'FlexibleSUSY_CMSSM'
gambit_backend_version = '2.4.0'
gambit_base_namespace  = ''


# ~~~~~ Information about the external code ~~~~~

# Use either absolute paths or paths relative to the main BOSS directory.
input_files   = [
'../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/src/lowe.h',
'../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/src/spectrum_generator_settings.hpp',
'../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/src/spectrum_generator_problems.hpp',
'../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/src/command_line_options.hpp',
#'../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/models/CMSSM/CMSSM_slha_io.hpp'
]
include_paths = ['../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/src', '../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/models/CMSSM']
base_paths    = ['../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/']

header_files_to = '../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/models/CMSSM'
src_files_to    = '../../../Backends/installed/flexiblesusy/'+gambit_backend_version+'/CMSSM/models/CMSSM'



load_classes = [
    'softsusy::QedQcd',
    'flexiblesusy::CMSSM_slha_io', 
    'flexiblesusy::CMSSM_scales',
    'flexiblesusy::Spectrum_generator_settings',
    'flexiblesusy::Spectrum_generator_problems',
    'flexiblesusy::Command_line_options::get_parameter_value',
    'flexiblesusy::Command_line_options::starts_with'
]

load_functions = [
]

ditch = [
    'softsusy::QedQcd::displayMass',
    'softsusy::QedQcd::displayAlphas',
    'softsusy::QedQcd::display_input_parameter_names',
    'flexiblesusy::Command_line_options::get_parameter_value'
]


auto_detect_stdlib_paths = False


load_parent_classes    = False
wrap_inherited_members = False


header_extension = '.h'
source_extension = '.cpp'

indent = 4


# ~~~~~ Information about other known types ~~~~~

# Dictionary key: type name
# Dictionary value: header file with containing type declaration.
#
# Example:
#   known_classes = {"SomeNamespace::KnownClassOne" : "path_to_header/KnownClassOne.hpp",
#                    "AnotherNamespace::KnownClassTwo" : "path_to_header/KnownClassTwo.hpp" }

known_classes = {}


# ~~~~~ Pragma directives for the inclusion of BOSSed classes in GAMBIT ~~~~~

# The listed pragma directives will be added before/after including the
# the BOSS-generated headers in GAMBIT.
#
# Example:
#   pragmas_begin = [
#       '#pragma GCC diagnostic push',
#       '#pragma GCC diagnostic ignored "-Wdeprecated-declarations"',
#   ]
#
#   pragmas_end = [
#       '#pragma GCC diagnostic pop'
#   ]

pragmas_begin = [
       '#pragma GCC diagnostic push',
       '#pragma GCC diagnostic ignored "-Wdeprecated-declarations"',
]
pragmas_end = [
       '#pragma GCC diagnostic pop'
]


