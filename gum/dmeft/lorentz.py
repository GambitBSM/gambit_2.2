# This file was automatically created by FeynRules 2.3.24
# Mathematica version: 11.0.1 for Linux x86 (64-bit) (September 21, 2016)
# Date: Wed 13 Feb 2019 17:51:47


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot
try:
   import form_factors as ForFac 
except ImportError:
   pass


UUS1 = Lorentz(name = 'UUS1',
               spins = [ -1, -1, 1 ],
               structure = '1')

UUV1 = Lorentz(name = 'UUV1',
               spins = [ -1, -1, 3 ],
               structure = 'P(3,2) + P(3,3)')

SSS1 = Lorentz(name = 'SSS1',
               spins = [ 1, 1, 1 ],
               structure = '1')

FFS1 = Lorentz(name = 'FFS1',
               spins = [ 2, 2, 1 ],
               structure = 'ProjM(2,1)')

FFS2 = Lorentz(name = 'FFS2',
               spins = [ 2, 2, 1 ],
               structure = 'ProjM(2,1) - ProjP(2,1)')

FFS3 = Lorentz(name = 'FFS3',
               spins = [ 2, 2, 1 ],
               structure = 'ProjP(2,1)')

FFS4 = Lorentz(name = 'FFS4',
               spins = [ 2, 2, 1 ],
               structure = 'ProjM(2,1) + ProjP(2,1)')

FFV1 = Lorentz(name = 'FFV1',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,1)')

FFV2 = Lorentz(name = 'FFV2',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFV3 = Lorentz(name = 'FFV3',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) - 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV4 = Lorentz(name = 'FFV4',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV5 = Lorentz(name = 'FFV5',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 4*Gamma(3,2,-1)*ProjP(-1,1)')

VSS1 = Lorentz(name = 'VSS1',
               spins = [ 3, 1, 1 ],
               structure = 'P(1,2) - P(1,3)')

VVS1 = Lorentz(name = 'VVS1',
               spins = [ 3, 3, 1 ],
               structure = 'Metric(1,2)')

VVV1 = Lorentz(name = 'VVV1',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

SSSS1 = Lorentz(name = 'SSSS1',
                spins = [ 1, 1, 1, 1 ],
                structure = '1')

FFFF1 = Lorentz(name = 'FFFF1',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(2,1)*Gamma5(4,3)')

FFFF2 = Lorentz(name = 'FFFF2',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,1)*Gamma(-1,4,3)')

FFFF3 = Lorentz(name = 'FFFF3',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(-2,1)*Gamma(-1,2,-2)*Gamma(-1,4,3)')

FFFF4 = Lorentz(name = 'FFFF4',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(-2,3)*Gamma(-1,2,1)*Gamma(-1,4,-2)')

FFFF5 = Lorentz(name = 'FFFF5',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(-3,3)*Gamma5(-2,1)*Gamma(-1,2,-2)*Gamma(-1,4,-3)')

FFFF6 = Lorentz(name = 'FFFF6',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(4,3)*Identity(2,1)')

FFFF7 = Lorentz(name = 'FFFF7',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(2,1)*Identity(4,3)')

FFFF8 = Lorentz(name = 'FFFF8',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Identity(2,1)*Identity(4,3)')

FFFF9 = Lorentz(name = 'FFFF9',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Sigma(-1,-2,2,1)*Sigma(-1,-2,4,3)')

FFFF10 = Lorentz(name = 'FFFF10',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma5(-3,1)*Sigma(-1,-2,2,-3)*Sigma(-1,-2,4,3)')

FFVV1 = Lorentz(name = 'FFVV1',
                spins = [ 2, 2, 3, 3 ],
                structure = '-4*Epsilon(3,4,-1,-2)*P(-2,4)*P(-1,3)*Gamma5(2,1) + 4*Epsilon(3,4,-1,-2)*P(-2,3)*P(-1,4)*Gamma5(2,1)')

FFVV2 = Lorentz(name = 'FFVV2',
                spins = [ 2, 2, 3, 3 ],
                structure = '-4*Epsilon(3,4,-1,-2)*P(-2,4)*P(-1,3)*Identity(2,1) + 4*Epsilon(3,4,-1,-2)*P(-2,3)*P(-1,4)*Identity(2,1)')

FFVV3 = Lorentz(name = 'FFVV3',
                spins = [ 2, 2, 3, 3 ],
                structure = 'P(3,4)*P(4,3)*Gamma5(2,1) - P(-1,3)*P(-1,4)*Gamma5(2,1)*Metric(3,4)')

FFVV4 = Lorentz(name = 'FFVV4',
                spins = [ 2, 2, 3, 3 ],
                structure = 'P(3,4)*P(4,3)*Identity(2,1) - P(-1,3)*P(-1,4)*Identity(2,1)*Metric(3,4)')

VVSS1 = Lorentz(name = 'VVSS1',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Metric(1,2)')

VVVV1 = Lorentz(name = 'VVVV1',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVV2 = Lorentz(name = 'VVVV2',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVV3 = Lorentz(name = 'VVVV3',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVV4 = Lorentz(name = 'VVVV4',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVV5 = Lorentz(name = 'VVVV5',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

FFVVV1 = Lorentz(name = 'FFVVV1',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = '-4*Epsilon(3,4,5,-1)*P(-1,3)*Gamma5(2,1) - 4*Epsilon(3,4,5,-1)*P(-1,4)*Gamma5(2,1) - 4*Epsilon(3,4,5,-1)*P(-1,5)*Gamma5(2,1)')

FFVVV2 = Lorentz(name = 'FFVVV2',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = '-4*Epsilon(3,4,5,-1)*P(-1,3)*Identity(2,1) - 4*Epsilon(3,4,5,-1)*P(-1,4)*Identity(2,1) - 4*Epsilon(3,4,5,-1)*P(-1,5)*Identity(2,1)')

FFVVV3 = Lorentz(name = 'FFVVV3',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'P(5,3)*Gamma5(2,1)*Metric(3,4) - P(5,4)*Gamma5(2,1)*Metric(3,4) - P(4,3)*Gamma5(2,1)*Metric(3,5) + P(4,5)*Gamma5(2,1)*Metric(3,5) + P(3,4)*Gamma5(2,1)*Metric(4,5) - P(3,5)*Gamma5(2,1)*Metric(4,5)')

FFVVV4 = Lorentz(name = 'FFVVV4',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'P(5,3)*Identity(2,1)*Metric(3,4) - P(5,4)*Identity(2,1)*Metric(3,4) - P(4,3)*Identity(2,1)*Metric(3,5) + P(4,5)*Identity(2,1)*Metric(3,5) + P(3,4)*Identity(2,1)*Metric(4,5) - P(3,5)*Identity(2,1)*Metric(4,5)')

FFVVVV1 = Lorentz(name = 'FFVVVV1',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Gamma5(2,1)*Metric(3,6)*Metric(4,5) - Gamma5(2,1)*Metric(3,5)*Metric(4,6)')

FFVVVV2 = Lorentz(name = 'FFVVVV2',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Identity(2,1)*Metric(3,6)*Metric(4,5) - Identity(2,1)*Metric(3,5)*Metric(4,6)')

FFVVVV3 = Lorentz(name = 'FFVVVV3',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Gamma5(2,1)*Metric(3,6)*Metric(4,5) - Gamma5(2,1)*Metric(3,4)*Metric(5,6)')

FFVVVV4 = Lorentz(name = 'FFVVVV4',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Gamma5(2,1)*Metric(3,5)*Metric(4,6) - Gamma5(2,1)*Metric(3,4)*Metric(5,6)')

FFVVVV5 = Lorentz(name = 'FFVVVV5',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Identity(2,1)*Metric(3,6)*Metric(4,5) - Identity(2,1)*Metric(3,4)*Metric(5,6)')

FFVVVV6 = Lorentz(name = 'FFVVVV6',
                  spins = [ 2, 2, 3, 3, 3, 3 ],
                  structure = 'Identity(2,1)*Metric(3,5)*Metric(4,6) - Identity(2,1)*Metric(3,4)*Metric(5,6)')

