#!/bin/python
#
# GAMBIT: Global and Modular BSM Inference Tool
#*********************************************
# \file
#
#  Model harvesting script.
#  Generates model_rollcall.hpp
#  
#  This script identifies all the headers that 
#  define GAMBIT models, and includes them
#  in model_rollcall unless asked not to.
#
#*********************************************
#
#  Authors (add name and date if you modify):
#
#  \author Pat Scott 
#          (patscott@physics.mcgill.ca)
#    \date 2014 Nov
#
#*********************************************
import os
execfile("./Utils/scripts/harvesting_tools.py")

def main(argv):

    model_headers=set([])
    model_type_headers = set([])
    exclude_models=set([])

    # Handle command line options
    verbose = False
    try:
        opts, args = getopt.getopt(argv,"vx:",["verbose","exclude-models="])
    except getopt.GetoptError:
        print 'Usage: model_harvestor.py [flags]'
        print ' flags:'
        print '        -v                   : More verbose output'  
        print '        -x model1,model2,... : Exclude model1, model2, etc.' 
        sys.exit(2)
    for opt, arg in opts:
      if opt in ('-v','--verbose'):
        verbose = True
        print 'model_harvester.py: verbose=True'
      elif opt in ('-x','--exclude-models'):
        exclude_models.update(neatsplit(",",arg))

    # Get list of models to include in models_rollcall.hpp
    model_headers.update(retrieve_generic_headers(verbose,"./Models/include/models","model",exclude_models))   
    # Get lists of model type header files
    #model_type_headers.update(retrieve_type_headers(verbose,"./Models/include/model_types",exclude_models))

    print "\nModel headers identified:"
    for h in model_headers:
        print ' ',h
    print "Model type headers identified:"
    for h in model_type_headers:
        print ' ',h
    print

    # Generate a c++ header containing all the frontend headers we have just harvested.
    towrite = "\
//   GAMBIT: Global and Modular BSM Inference Tool\n\
//   *********************************************\n\
///  \\file                                       \n\
///                                               \n\
///  Compile-time registration of GAMBIT models.  \n\
///                                               \n\
///  This file was automatically generated by     \n\
///  model_harvester.py. Do not modify.           \n\
///                                               \n\
///  Do not add to this if you want to add a new  \n\
///  model -- just add your model header to       \n\
///  Models/include/models and rest assured that  \n\
///  model_harvester.py will make sure it ends    \n\
///  up here.                                     \n\
///                                               \n\
///  *********************************************\n\
///                                               \n\
///  Authors (add name and date if you modify):   \n\
///                                               \n\
///  \\author The GAMBIT Collaboration            \n\
///  \date "+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+"\n\
///                                               \n\
///  *********************************************\n\
                                                  \n\
#ifndef __model_rollcall_hpp__                    \n\
#define __model_rollcall_hpp__                    \n\
                                                  \n\
// Include the model macro definitions            \n\
#include \"model_macros.hpp\"                     \n\
                                                  \n\
// Automatically-generated list of models.        \n"

    for h in model_headers:
        towrite+='#include \"models/{0}\"\n'.format(h)
    towrite+="\n#endif // defined __model_rollcall_hpp__\n"
    
    with open("./Models/include/model_rollcall.hpp","w") as f:
        f.write(towrite)

    if verbose: print "Generated model_rollcall.hpp.\n" 


# Handle command line arguments (verbosity)
if __name__ == "__main__":
   main(sys.argv[1:])
