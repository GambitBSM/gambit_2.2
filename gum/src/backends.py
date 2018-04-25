"""
Master file containing all routines for modifying Backend interfaces.
Think CalcHEP, MadGraph, SPheno.
"""

def add_calchep_switch(model_name):
    """
    Adds an 'if ModelInUse()' switch to the CalcHEP frontend to make GAMBIT 
    point to the correct CalcHEP files.
    """
