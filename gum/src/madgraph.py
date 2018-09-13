"""
Master routines for all MadGraph related routines.
"""

from distutils.dir_util import copy_tree

def call_madgraph(message):
    """
    Calls MadGraph.
    """

def copy_madgraph_files(mgdir, model_name):
    """
    Copies all MadGraph files into the contrib/MadGraph directory.
    """

    model_name.strip('/')

    target = './contrib/MadGraph/models/' + model_name

    copy_tree(mgdir, target)

    print("MadGraph files moved to correct directory.")
