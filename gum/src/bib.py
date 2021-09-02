#  GUM: GAMBIT Universal Model Machine
#  ***********************************
#  \file
#
#  Contains all routines for generating
#  a list of bibtex tags for all used
#  software. 
#
#  *************************************
#
#  \author Christopher Chang
#  \date 2021
#
#  **************************************

def write_bib_file(bibtags,bibcontents,output_dir,name):
  """
  Write a <model>.bib file and a <model>_bibtags.txt file
  containing necessary citations when using GUM.
  """

  open(output_dir + "/" + name + "_bibtags.txt", 'w').write(bibtags)
  open(output_dir + "/" + name + ".bib", 'w').write(bibcontents)
  
  print("Generated a bib file in the Outputs folder")

  return


def generate_bib_tags(output_opts,gum_math):
  """
  For each used piece of software relevant for gum,
  generate the list of bib tags and bib entries.
  """
  
  bibtags = "Please use these bib tags when referencing the software used in GUM: \n\n"
  contents = "% Bib file for the software used with GUM.\n\n"
  
  # GUM. TODO: Needs updating once published
  bibtags += "GUM: GUM \n \n"
  
  contents += (
  "@misc{GUM,\n"
  "  author = \"Bloor, Sanjay and Gonzalo, Tomas E. and Scott, Pat and others\",\n"
  "  title = \"{The GAMBIT Universal Model Machine: from Lagrangians to Likelihoods}\",\n"
  "  year = \"2021\",\n"
  "  note = \"{in preparation}\"\n"
  "}\n"
  )
  
  # FeynRules
  if gum_math == 'feynrules':

    bibtags += "FeynRules: Alloul:2013bka,Christensen:2008py,Christensen:2010wz,Christensen:2009jx,2012CoPhC.183.1201D \n\n"
    
    # FeynRules 2.0 Manual
    contents += (
    "@article{Alloul:2013bka,\n"
    "author         = \"Alloul, Adam and Christensen, Neil D. and Degrande,Céline and Duhr, Claude and Fuks, Benjamin\",\n"
    "title          = \"{FeynRules  2.0 - A complete toolbox for tree-level\n"
    "phenomenology}\",\n"
    "journal        = \"Comput. Phys. Commun.\",\n"
    "volume         = \"185\",\n"
    "year           = \"2014\",\n"
    "pages          = \"2250-2300\",\n"
    "doi            = \"10.1016/j.cpc.2014.04.012\",\n"
    "eprint         = \"1310.1921\",\n"
    "archivePrefix  = \"arXiv\",\n"
    "primaryClass   = \"hep-ph\",\n"
    "reportNumber   = \"CERN-PH-TH-2013-239, MCNET-13-14, IPPP-13-71,DCPT-13-142, PITT-PACC-1308\",\n"
    "SLACcitation   = \"%%CITATION = ARXIV:1310.1921;%%\"\n"
    "}\n" 
    )
    
    # FeynRules made easy
    contents += (
    "@article{Christensen:2008py,\n"
    "  author         = \"Christensen, Neil D. and Duhr, Claude\",\n"
    "  title          = \"{FeynRules - Feynman rules made easy}\",\n"
    "  journal        = \"Comput. Phys. Commun.\",\n"
    "  volume         = \"180\",\n"
    "  year           = \"2009\",\n"
    "  pages          = \"1614-1641\",\n"
    "  doi            = \"10.1016/j.cpc.2009.02.018\",\n"
    "  eprint         = \"0806.4194\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"MSUHEP-080616, CP3-08-20\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:0806.4194;%%\"\n"
    "}\n"
    )
    
    # Whizard Interface
    contents += (
    "@article{Christensen:2010wz,\n"
    "  author         = \"Christensen, Neil D. and Duhr, Claude and Fuks, Benjamin\n"
    "                    and Reuter, Jurgen and Speckner, Christian\",\n"
    "  title          = \"{Introducing an interface between WHIZARD and FeynRules}\",\n"
    "  journal        = \"Eur. Phys. J.\",\n"
    "  volume         = \"C72\",\n"
    "  year           = \"2012\",\n"
    "  pages          = \"1990\",\n"
    "  doi            = \"10.1140/epjc/s10052-012-1990-5\",\n"
    "  eprint         = \"1010.3251\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"FR-PHENO-2010-030, IPHC-PHENO-10-03, IPPP-10-80,\n"
    "                    DCPT-10-160, MADPH-10-1562, EDINBURGH-2010-25,\n"
    "                    DESY-12-020, PITT-PACC-12002\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1010.3251;%%\"\n"
    "}\n"
    )
    # CalcHEP, FeynArts, Sherpa interfaces
    contents += (
    "@article{Christensen:2009jx,\n"
    "  author         = \"Christensen, Neil D. and de Aquino, Priscila and\n"
    "                    Degrande, Celine and Duhr, Claude and Fuks, Benjamin and\n"
    "                    Herquet, Michel and Maltoni, Fabio and Schumann, Steffen\",\n"
    "  title          = \"{A Comprehensive approach to new physics simulations}\",\n"
    "  journal        = \"Eur. Phys. J.\",\n"
    "  volume         = \"C71\",\n"
    "  year           = \"2011\",\n"
    "  pages          = \"1541\",\n"
    "  doi            = \"10.1140/epjc/s10052-011-1541-5\",\n"
    "  eprint         = \"0906.2474\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"CP3-09-24, HD-THEP-09-11, IPHC-PHENO-09-01,\n"
    "                    MSUHEP-090612, NIKHEF-2009-009\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:0906.2474;%%\"\n"
    "}\n"
    )
    # Universal FeynRules Output
    contents += (
    "@ARTICLE{2012CoPhC.183.1201D,\n"
    "  author = {{Degrande}, C{\\'e}line and {Duhr}, Claude and {Fuks}, Benjamin and {Grellscheid}, David and {Mattelaer}, Olivier and {Reiter}, Thomas},\n"
    "  title = \"{UFO - The Universal FEYNRULES Output}\",\n"
    "  journal = {Computer Physics Communications},\n"
    "  keywords = {High Energy Physics - Phenomenology},\n"
    "  year = 2012,\n"
    "  month = jun,\n"
    "  volume = {183},\n"
    "  number = {6},\n"
    "  pages = {1201-1214},\n"
    "  doi = {10.1016/j.cpc.2012.01.022},\n"
    "  archivePrefix = {arXiv},\n"
    "  eprint = {1108.2040},\n"
    "  primaryClass = {hep-ph},\n"
    "  adsurl = {https://ui.adsabs.harvard.edu/abs/2012CoPhC.183.1201D},\n"
    "  adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n"
    "}\n"
    )

  
  # SARAH
  if gum_math == 'sarah':

    bibtags += "SARAH: Staub:2013tta,2012pb,Staub:2008uz,Staub:2010jh,Staub:2015kfa\n\n"
    
    # SARAH 4
    contents += (
    "@article{Staub:2013tta,\n"
    "  author         = \"Staub, Florian\",\n"
    "  title          = \"{SARAH 4 : A tool for (not only SUSY) model builders}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"185\",\n"
    "  year           = \"2014\",\n"
    "  pages          = \"1773-1790\",\n"
    "  doi            = \"10.1016/j.cpc.2014.02.018\",\n"
    "  eprint         = \"1309.7223\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"BONN-TH-2013-17\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1309.7223;%%\"\n"
    "}\n"
    )
    
    # SARAH 3.2
    contents += (
    "@article{Staub:2012pb,\n"
    "  author         = \"Staub, Florian\",\n"
    "  title          = \"{SARAH 3.2: Dirac Gauginos, UFO output, and more}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"184\",\n"
    "  year           = \"2013\",\n"
    "  pages          = \"1792-1809\",\n"
    "  doi            = \"10.1016/j.cpc.2013.02.019\",\n"
    "  eprint         = \"1207.0906\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"BONN-TH-2012-17\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1207.0906;%%\"\n"
    "}\n"
    )
    
    # SARAH
    contents += (
    "@article{Staub:2008uz,\n"
    "  author         = \"Staub, F.\",\n"
    "  title          = \"{SARAH}\",\n"
    "  year           = \"2008\",\n"
    "  eprint         = \"0806.0538\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:0806.0538;%%\"\n"
    "}\n"
    )
    
    # Automatic Calculation of Supersymmetric Renormalization Group Equations and Self Energies
    contents += (
    "@article{Staub:2010jh,\n"
    "  author         = \"Staub, Florian\",\n"
    "  title          = \"{Automatic Calculation of Supersymmetric\n"
    "                     Renormalization Group Equations and Self Energies}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"182\",\n"
    "  year           = \"2011\",\n"
    "  pages          = \"808-833\",\n"
    "  doi            = \"10.1016/j.cpc.2010.11.030\",\n"
    "  eprint         = \"1002.0840\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1002.0840;%%\"\n"
    "}\n"
    )
    
    # Exploring new models in detail with SARAH
    contents += (
    "@article{Staub:2015kfa,\n"
    "  author         = \"Staub, Florian\",\n"
    "  title          = \"{Exploring new models in all detail with SARAH}\",\n"
    "  journal        = \"Adv. High Energy Phys.\",\n"
    "  volume         = \"2015\",\n"
    "  year           = \"2015\",\n"
    "  pages          = \"840780\",\n"
    "  doi            = \"10.1155/2015/840780\",\n"
    "  eprint         = \"1503.04200\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"CERN-PH-TH-2015-051\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1503.04200;%%\"\n"
    "}\n"
    )

  # CalcHEP
  if output_opts.ch:

    bibtags += "CalcHEP: Pukhov:2004ca,Belyaev:2012qa \n\n"

    # Calchep 2.3    
    contents += (
    "@article{Pukhov:2004ca,\n"
    "  author         = \"Pukhov, A.\",\n"
    "  title          = \"{CalcHEP 2.3: MSSM, structure functions, event\n"
    "                    generation, batchs, and generation of matrix elements for\n"
    "                    other packages}\",\n"
    "  year           = \"2004\",\n"
    "  eprint         = \"hep-ph/0412191\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0412191;%%\"\n"
    "}\n"
    )
    
    # Calchep 3.4
    contents += (
    "@article{Belyaev:2012qa,\n"
    "  author         = \"Belyaev, Alexander and Christensen, Neil D. and Pukhov,\n"
    "                    Alexander\",\n"
    "  title          = \"{CalcHEP 3.4 for collider physics within and beyond the\n"
    "                    Standard Model}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"184\",\n"
    "  year           = \"2013\",\n"
    "  pages          = \"1729-1769\",\n"
    "  doi            = \"10.1016/j.cpc.2013.01.014\",\n"
    "  eprint         = \"1207.6082\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"PITT-PACC-1209\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1207.6082;%%\"\n"
    "}\n"
    )
  
  # Vevacious
  if output_opts.vev:

    bibtags += "Vevacious: Camargo-Molina:2013qva \n\n"
    
    contents += (
    "@article{Camargo-Molina:2013qva,\n"
    "  author         = \"Camargo-Molina, J. E. and O'Leary, B. and Porod, W. and\n"
    "                    Staub, F.\",\n"
    "  title          = \"{\\textsf{Vevacious}: A Tool For Finding The Global\n"
    "                    Minima Of One-Loop Effective Potentials With Many\n"
    "                    Scalars}\",\n"
    "  journal        = \"Eur. Phys. J.\",\n"
    "  volume         = \"C73\",\n"
    "  year           = \"2013\",\n"
    "  number         = \"10\",\n"
    "  pages          = \"2588\",\n"
    "  doi            = \"10.1140/epjc/s10052-013-2588-2\",\n"
    "  eprint         = \"1307.1477\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1307.1477;%%\"\n"
    "}\n"
    )
    
  # MadGraph
  if (output_opts.ufo and output_opts.pythia):
  
    bibtags += "MadGraph: Stelzer:1994ta,Maltoni:2002qb,Alwall:2007st,Alwall:2011uj \n\n"
  
    # Generation of Helicity Amplitudes
    contents += (
    "@article{Stelzer:1994ta,\n"
    "  author         = \"Stelzer, T. and Long, W. F.\",\n"
    "  title          = \"{Automatic generation of tree level helicity amplitudes}\",\n"
    "  journal        = \"Comput. Phys. Commun.\",\n"
    "  volume         = \"81\",\n"
    "  year           = \"1994\",\n"
    "  pages          = \"357-371\",\n"
    "  doi            = \"10.1016/0010-4655(94)90084-1\",\n"
    "  eprint         = \"hep-ph/9401258\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"MAD-PH-813\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/9401258;%%\"\n"
    "}\n"
    )
    
    # Event generation with MadEvent
    contents += (
    "@article{Maltoni:2002qb,\n"
    "  author         = \"Maltoni, Fabio and Stelzer, Tim\",\n"
    "  title          = \"{MadEvent: Automatic event generation with MadGraph}\",\n"
    "  journal        = \"\jhep\",\n"
    "  volume         = \"02\",\n"
    "  year           = \"2003\",\n"
    "  pages          = \"027\",\n"
    "  doi            = \"10.1088/1126-6708/2003/02/027\",\n"
    "  eprint         = \"hep-ph/0208156\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0208156;%%\"\n"
    "}\n"
    )
    
    # Web Generation
    contents += (
    "@article{Alwall:2007st,\n"
    "  author         = \"Alwall, Johan and Demin, Pavel and de Visscher, Simon and\n"
    "                    Frederix, Rikkert and Herquet, Michel and Maltoni, Fabio\n"
    "                    and Plehn, Tilman and Rainwater, David L. and Stelzer,\n"
    "                    Tim\",\n"
    "  title          = \"{MadGraph/MadEvent v4: The New Web Generation}\",\n"
    "  journal        = \"\jhep\",\n"
    "  volume         = \"09\",\n"
    "  year           = \"2007\",\n"
    "  pages          = \"028\",\n"
    "  doi            = \"10.1088/1126-6708/2007/09/028\",\n"
    "  eprint         = \"0706.2334\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"SLAC-PUB-12603, CP3-07-17\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:0706.2334;%%\"\n"
    "}\n"
    )
    
    
    # MadGrah 5: Going Beyond
    contents += (
    "@article{Alwall:2011uj,\n"
    "  author         = \"Alwall, Johan and Herquet, Michel and Maltoni, Fabio and\n"
    "  Mattelaer, Olivier and Stelzer, Tim\",\n"
    "  title          = \"{MadGraph 5 : Going Beyond}\",\n"
    "  journal        = \"\jhep\",\n"
    "  volume         = \"06\",\n"
    "  year           = \"2011\",\n"
    "  pages          = \"128\",\n"
    "  doi            = \"10.1007/JHEP06(2011)128\",\n"
    "  eprint         = \"1106.0522\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"FERMILAB-PUB-11-448-T\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1106.0522;%%\"\n"
    "}\n"
    )
    
  # Pythia
  if output_opts.pythia:
  
    bibtags += "Pythia: Sjostrand:2014zea \n\n"
    
    # Pythia 8.2 Manual
    contents += (
    "@article{Sjostrand:2014zea,\n"
    "  author         = \"Sjostrand, Torbjorn and Ask, Stefan and Christiansen,\n"
    "                    Jesper R. and Corke, Richard and Desai, Nishita and Ilten,\n"
    "                    Philip and Mrenna, Stephen and Prestel, Stefan and\n"
    "                    Rasmussen, Christine O. and Skands, Peter Z.\",\n"
    "  title          = \"{An Introduction to PYTHIA 8.2}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"191\",\n"
    "  year           = \"2015\",\n"
    "  pages          = \"159-177\",\n"
    "  doi            = \"10.1016/j.cpc.2015.01.024\",\n"
    "  eprint         = \"1410.3012\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LU-TP-14-36, MCNET-14-22, CERN-PH-TH-2014-190,\n"
    "                    FERMILAB-PUB-14-316-CD, DESY-14-178, SLAC-PUB-16122\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1410.3012;%%\"\n"
    "}\n"
    )
    
  
  # Micromegas
  if output_opts.mo:

    bibtags += "Micromegas: micromegas,Belanger:2013oya,Belanger:2010gh,Belanger:2008sj,Belanger:2006is,Belanger:2004yn,Belanger:2001fz \n\n" 
    
    # Micromegas 4.1
    contents += (
    "@ARTICLE{micromegas,\n"
    "  author = {{B{\\'e}langer}, G. and {Boudjema}, F. and {Pukhov}, A. and {Semenov}, A.},\n"
    "  title = \"{micrOMEGAs4.1: Two dark matter candidates}\",\n"
    "  journal = {\cpc},\n"
    "  archivePrefix = \"arXiv\",\n"
    "  eprint = {1407.6129},\n"
    "  primaryClass = \"hep-ph\",\n"
    "  keywords = {Dark matter, Relic density, Indirect detection, MSSM, Beyond standard model},\n"
    "  year = 2015,\n"
    "  month = jul,\n"
    "  volume = 192,\n"
    "  pages = {322-329},\n"
    "  doi = {10.1016/j.cpc.2015.03.003},\n"
    "  adsurl = {http://adsabs.harvard.edu/abs/2015CoPhC.192..322B},\n"
    "  adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n"
    "}\n"
    )
    
    
    # Micromegas 3
    contents += (
    "@article{Belanger:2013oya,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Pukhov, A. and Semenov,\n"
    "                    A.\",\n"
    "  title          = \"{micrOMEGAs 3: A program for calculating dark matter\n"
    "                    observables}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"185\",\n"
    "  year           = \"2014\",\n"
    "  pages          = \"960-985\",\n"
    "  doi            = \"10.1016/j.cpc.2013.10.016\",\n"
    "  eprint         = \"1305.0237\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LAPTH-023-13\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1305.0237;%%\"\n"
    "}\n"
    )
    
    # Micromegas 2.4
    contents += (
    "@article{Belanger:2010gh,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Brun, P. and Pukhov, A.\n"
    "                    and Rosier-Lees, S. and Salati, P. and Semenov, A.\",\n"
    "  title          = \"{Indirect search for dark matter with micrOMEGAs2.4}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"182\",\n"
    "  year           = \"2011\",\n"
    "  pages          = \"842-856\",\n"
    "  doi            = \"10.1016/j.cpc.2010.11.033\",\n"
    "  eprint         = \"1004.1092\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"IRFU-10-24, LAPTH-012-10.\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1004.1092;%%\"\n"
    "}\n"
    )
    
    # Micromegas 2.2
    contents += (
    "@article{Belanger:2008sj,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Pukhov, A. and Semenov,\n"
    "                    A.\",\n"
    "  title          = \"{Dark matter direct detection rate in a generic model\n"
    "                    with micrOMEGAs 2.2}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"180\",\n"
    "  year           = \"2009\",\n"
    "  pages          = \"747-767\",\n"
    "  doi            = \"10.1016/j.cpc.2008.11.019\",\n"
    "  eprint         = \"0803.2360\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LAPTH-1237-08\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:0803.2360;%%\"\n"
    "}\n"
    )
    
    # Micromegas 2.0
    contents += (
    "@article{Belanger:2006is,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Pukhov, A. and Semenov,\n"
    "                    A.\",\n"
    "  title          = \"{MicrOMEGAs 2.0: A Program to calculate the relic density\n"
    "                    of dark matter in a generic model}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"176\",\n"
    "  year           = \"2007\",\n"
    "  pages          = \"367-382\",\n"
    "  doi            = \"10.1016/j.cpc.2006.11.008\",\n"
    "  eprint         = \"hep-ph/0607059\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LAPTH-1152-06\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0607059;%%\"\n"
    "}\n"
    )
    
    # Micromegas 1.3
    contents += (
    "@article{Belanger:2004yn,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Pukhov, A. and Semenov,\n"
    "                    A.\",\n"
    "  title          = \"{micrOMEGAs: Version 1.3}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"174\",\n"
    "  year           = \"2006\",\n"
    "  pages          = \"577-604\",\n"
    "  doi            = \"10.1016/j.cpc.2005.12.005\",\n"
    "  eprint         = \"hep-ph/0405253\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LAPTH-1044\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0405253;%%\"\n"
    "}\n"
    )
    
    # Micromegas 1.0
    contents += (
    "@article{Belanger:2001fz,\n"
    "  author         = \"B{\\'e}langer, G. and Boudjema, F. and Pukhov, A. and Semenov,\n"
    "                    A.\",\n"
    "  title          = \"{MicrOMEGAs: A Program for calculating the relic density\n"
    "                    in the MSSM}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"149\",\n"
    "  year           = \"2002\",\n"
    "  pages          = \"103-120\",\n"
    "  doi            = \"10.1016/S0010-4655(02)00596-9\",\n"
    "  eprint         = \"hep-ph/0112278\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"LAPTH-881-01\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0112278;%%\"\n"
    "}\n"
    )
    
  # Spheno
  if output_opts.spheno:

    bibtags += "Spheno: Porod:2003um, Porod:2011nf \n\n"
    
    # Spheno
    contents += (
    "@article{Porod:2003um,\n"
    "  author         = \"Porod, Werner\",\n"
    "  title          = \"{SPheno, a program for calculating supersymmetric\n"
    "                    spectra, SUSY particle decays and SUSY particle production\n"
    "                    at $e^+e^-$ colliders}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"153\",\n"
    "  year           = \"2003\",\n"
    "  pages          = \"275-315\",\n"
    "  doi            = \"10.1016/S0010-4655(03)00222-4\",\n"
    "  eprint         = \"hep-ph/0301101\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  reportNumber   = \"ZU-TH-01-03\",\n"
    "  SLACcitation   = \"%%CITATION = HEP-PH/0301101;%%\"\n"
    "}\n"
    )
    
    # Spheno 3.1
    contents += (
    "@article{Porod:2011nf,\n"
    "  author         = \"Porod, W. and Staub, F.\",\n"
    "  title          = \"{SPheno 3.1: Extensions including flavour, CP-phases and\n"
    "                    models beyond the MSSM}\",\n"
    "  journal        = \"\cpc\",\n"
    "  volume         = \"183\",\n"
    "  year           = \"2012\",\n"
    "  pages          = \"2458-2469\",\n"
    "  doi            = \"10.1016/j.cpc.2012.05.021\",\n"
    "  eprint         = \"1104.1573\",\n"
    "  archivePrefix  = \"arXiv\",\n"
    "  primaryClass   = \"hep-ph\",\n"
    "  SLACcitation   = \"%%CITATION = ARXIV:1104.1573;%%\"\n"
    "}\n"
    )
    
  return bibtags, contents
