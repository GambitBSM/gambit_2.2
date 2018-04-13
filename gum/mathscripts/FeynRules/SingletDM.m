$FeynRulesPath = SetDirectory["~/.Mathematica/Applications/feynrules-current"];
<< FeynRules`;
LoadModel["Models/SingletDM/SingletDM.fr"];
LoadRestriction["Models/DiagonalCKM.rst"];
CheckHermiticity[LTotal];
FeynmanGauge = False;
WriteCHOutput[LTotal, CHAutoWidths -> False];
Quit[];
