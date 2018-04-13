SetDirectory["~/.Mathematica/Applications/SARAH-4.9.1"];
<< SARAH`;
Start["SingletDM"];
ModelOutput[EWSB];
MakeVertexList[EWSB];
MakeCHep[SLHAinput -> False, UseRunningCoupling -> True, CalculateMasses -> True];

Quit[];
