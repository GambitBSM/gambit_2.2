from __future__ import print_function
import numpy as np

#
# This will run when this backend is first loaded, i.e. during GAMBIT startup
#
prefix = "salami_gambit:"
print(prefix, "Starting up...")

# This needs to be run from within GAMBIT to avoid the hardcoded absolute path
import salami
KPREDS = {}
SLHA = None


#
# Import SLHA content as string
#
def import_slha_string(slha_string):

    global SLHA
    print("Reading SLHA input...")
    import pyslha
    SLHA = pyslha.readSLHA(slha_string)
    print("...done")


# #
# # Set parameters
# #
# def set_parameters(params):

#     for k,v in params.items():
#         xsec.set_parameter(k, v)


#
# Return the cross-section for a given process, identified by a PID pair
#
def get_xsection(energy, proc, xsec_lo):

    pids_to_prosids = { [ 1000022,  1000022] : "11",
                        [ 1000022,  1000023] : "12",
                        [ 1000022,  1000025] : "13",
                        [ 1000022,  1000035] : "14",
                        [ 1000023,  1000023] : "22",
                        [ 1000023,  1000025] : "23",
                        [ 1000023,  1000035] : "24",
                        [ 1000025,  1000025] : "33",
                        [ 1000025,  1000035] : "34",
                        [ 1000035,  1000035] : "44",
                        [ 1000022,  1000024] : "15",
                        [ 1000022,  1000034] : "16",
                        [ 1000023,  1000024] : "25",
                        [ 1000023,  1000034] : "26",
                        [ 1000025,  1000024] : "35",
                        [ 1000025,  1000034] : "36",
                        [ 1000035,  1000024] : "45",
                        [ 1000035,  1000034] : "46",
                        [ 1000024,  1000024] : "55",
                        [ 1000024,  1000034] : "56",
                        [ 1000034,  1000034] : "66",
                        [-1000024, -1000024] : "77",
                        [-1000024, -1000034] : "78",
                        [-1000034, -1000034] : "88",
                        [ 1000024, -1000024] : "57",
                        [ 1000024, -1000034] : "58",
                        [ 1000034, -1000024] : "67",
                        [ 1000034, -1000034] : "68" }

    # Load process
    known_process = True
    print("Load process:", proc)
    try:
        energy = str(int(energy))
        pidkey = list(sorted(proc, key=abs))
        proskey = pids_to_prosids[pidkey]
        key = "{}_{}".format(energy, proskey)
        if proskey not in KPREDS:
            KPREDS[key] = salami.KPred(energy, proskey)
        kpred = KPREDS[key]
    except KeyError:
        print("Unknown process", proc, "-- setting cross-section to 0")
        known_process = False

    # Evaluate cross-sections and errors (factor of 1000 for Prospino pb -> Gambit fb conversion)
    result_dict = {}
    result_dict["central"] = 1000*kpred.predict_xsec(slha, xsec_lo) if known_process else 0.0
    # result_dict["regdown_rel"] = 0.0
    # result_dict["regup_rel"] = 0.0
    # result_dict["scaledown_rel"] = 0.0
    # result_dict["scaleup_rel"] = 0.0
    # result_dict["pdfdown_rel"] = 0.0
    # result_dict["pdfup_rel"] = 0.0
    # result_dict["alphasdown_rel"] = 0.0
    # result_dict["alphasup_rel"] = 0.0
    result_dict["tot_err_down"] = 0.0
    result_dict["tot_err_up"] = 0.0

    return result_dict



# # ================ TESTING ================

# params_in = {'par1': 1.0, 'par2': 2.0}
# set_parameters(params_in)

# flags_in = {'flag1': True, 'flag2': False}
# set_flags(flags_in)

# xsec_fb((1000021,1000021), {'LO_xsec_fb': 78.9}, {})

# xsec_err_fb((1000021,1000021), {}, {'alphas_error': True})
