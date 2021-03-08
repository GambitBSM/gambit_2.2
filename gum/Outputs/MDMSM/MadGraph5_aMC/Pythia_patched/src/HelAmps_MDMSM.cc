//==========================================================================
// This file has been automatically generated for Pythia 8 by
// MadGraph5_aMC@NLO v. 2.8.3.2, 2021-02-02
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_MDMSM.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace Pythia8_MDMSM 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = sqrt(0.5); 
  sqs = sqrt(0.5/3); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], sqrt(pt2 + p[3] * p[3])); 
  pt = min(pp, sqrt(pt2)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(std::labs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
    if (pp == 0.000)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[std::abs(ip)]; 
      fo[3] = ip * nsf * sqm[std::abs(ip)]; 
      fo[4] = im * nsf * sqm[std::abs(im)]; 
      fo[5] = ip * sqm[std::abs(im)]; 
    }
    else
    {
      pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/sqrt(2.0 * pp * pp3); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.00)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])); 
    if (pp == 0.0)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = sqrt(0.5); 
  hel = double(nhel); 
  nsvahl = nsv * std::abs(hel); 
  pt2 = (p[1] * p[1]) + (p[2] * p[2]); 
  pp = min(p[0], sqrt(pt2 + (p[3] * p[3]))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - std::abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2])); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void VVVVS1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (OM3 * P3[0] * (-cI * (TMP18 * TMP30) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[2] * TMP29) + cI * (V1[2] * TMP30)));
  V3[3] = denom * S5[2] * (OM3 * P3[1] * (-cI * (TMP18 * TMP30) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[3] * TMP29) + cI * (V1[3] * TMP30)));
  V3[4] = denom * S5[2] * (OM3 * P3[2] * (-cI * (TMP18 * TMP30) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[4] * TMP29) + cI * (V1[4] * TMP30)));
  V3[5] = denom * S5[2] * (OM3 * P3[3] * (-cI * (TMP18 * TMP30) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[5] * TMP29) + cI * (V1[5] * TMP30)));
}


void VVVS1_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> S4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP20; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0] + S4[0]; 
  V1[1] = +V2[1] + V3[1] + S4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP20 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S4[2] * (OM1 * P1[0] * (TMP27 * (-cI * (TMP20) + cI *
      (TMP19)) + (-cI * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 *
      (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI * (TMP23) + cI * (TMP24))
      + V3[2] * (-cI * (TMP26) + cI * (TMP22)))));
  V1[3] = denom * S4[2] * (OM1 * P1[1] * (TMP27 * (-cI * (TMP20) + cI *
      (TMP19)) + (-cI * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 *
      (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI * (TMP23) + cI * (TMP24))
      + V3[3] * (-cI * (TMP26) + cI * (TMP22)))));
  V1[4] = denom * S4[2] * (OM1 * P1[2] * (TMP27 * (-cI * (TMP20) + cI *
      (TMP19)) + (-cI * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 *
      (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI * (TMP23) + cI * (TMP24))
      + V3[4] * (-cI * (TMP26) + cI * (TMP22)))));
  V1[5] = denom * S4[2] * (OM1 * P1[3] * (TMP27 * (-cI * (TMP20) + cI *
      (TMP19)) + (-cI * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 *
      (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI * (TMP23) + cI * (TMP24))
      + V3[5] * (-cI * (TMP26) + cI * (TMP22)))));
}


void FFV3_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * (-2. * cI) * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + (+1./2. * (M1 * (+2. * (F2[4] * (-1./2.) * (V3[2] + V3[5]))
      - F2[5] * (V3[3] + cI * (V3[4])))) + F2[3] * (P1[0] * (V3[3] + cI *
      (V3[4])) + (P1[1] * (-1.) * (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI *
      (V3[2] + V3[5])) + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * (-2. * cI) * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) +
      (P1[1] * (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P1[3] * (+cI * (V3[4]) - V3[3])))) + (+1./2. * (M1 * (F2[5] * (V3[5] -
      V3[2]) + 2. * (F2[4] * 1./2. * (+cI * (V3[4]) - V3[3])))) + F2[3] *
      (P1[0] * (-1.) * (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) +
      (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * cI * (F2[4] * (P1[0] * (-1.) * (V3[2] + V3[5]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[2]
      + V3[5])))) + (F2[5] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (-1.) * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] *
      (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3]
      * (V3[2] - V3[5])))) + M1 * (F2[2] * 2. * (+cI * (V3[4]) - V3[3]) + 2. *
      (F2[3] * (V3[2] + V3[5])))));
}


void FFV4C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * (F2[2] * (P1[0] * (-1.) * (V3[2] + V3[5]) + (P1[1]
      * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] *
      (V3[2] + V3[5])))) + (+1./2. * (M1 * (F2[5] * (V3[3] - cI * (V3[4])) + 2.
      * (F2[4] * 1./2. * (V3[5] - V3[2])))) + F2[3] * (P1[0] * (+cI * (V3[4]) -
      V3[3]) + (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI *
      (V3[5])) + P1[3] * (V3[3] - cI * (V3[4])))))));
  F1[3] = denom * 2. * cI * (F2[2] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M1 * (+2. * (F2[4] * 1./2. *
      (V3[3] + cI * (V3[4]))) - F2[5] * (V3[2] + V3[5]))) + F2[3] * (P1[0] *
      (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI
      * (V3[3])) + P1[3] * (V3[5] - V3[2]))))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + (F2[5] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + M1 * (F2[2] * 2. * (V3[2] + V3[5]) + 2. * (F2[3]
      * (V3[3] - cI * (V3[4]))))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P1[3] * (V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] + V3[5]) +
      (P1[1] * (+cI * (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI *
      (V3[3])) - P1[3] * (V3[2] + V3[5])))) + M1 * (F2[2] * 2. * (V3[3] + cI *
      (V3[4])) + 2. * (F2[3] * (V3[2] - V3[5])))));
}


void FFV5C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP10; 
  std::complex<double> TMP16; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP10 = (-1.) * (F1[2] * (F2[4] * (P3[0] - P3[3]) + F2[5] * (+cI * (P3[2]) -
      P3[1])) + F1[3] * (F2[4] * (-1.) * (P3[1] + cI * (P3[2])) + F2[5] *
      (P3[0] + P3[3])));
  TMP16 = (-1.) * (F1[4] * (F2[2] * (P3[0] + P3[3]) + F2[3] * (P3[1] - cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] + cI * (P3[2])) + F2[3] * (P3[0] -
      P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * 4. * cI * (OM3 * 1./4. * P3[0] * (TMP10 + 4. * (TMP16)) +
      (+1./4. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * (-4. * cI) * (OM3 * - 1./4. * P3[1] * (TMP10 + 4. * (TMP16))
      + (-1./4. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * 4. * cI * (OM3 * 1./4. * P3[2] * (TMP10 + 4. * (TMP16)) +
      (-1./4. * cI * (F2[5] * F1[2]) + 1./4. * cI * (F2[4] * F1[3]) - cI *
      (F2[2] * F1[5]) + cI * (F2[3] * F1[4])));
  V3[5] = denom * 4. * cI * (OM3 * 1./4. * P3[3] * (TMP10 + 4. * (TMP16)) +
      (+1./4. * (F2[4] * F1[2]) - 1./4. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void FFS3_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) - F1[2] * M2));
  F2[3] = denom * cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] * (-1.)
      * (P2[0] + P2[3]) + F1[3] * M2));
  F2[4] = denom * - cI * S3[2] * (F1[2] * (-1.) * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) + F1[4] * M2));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) - F1[5] * M2));
}


void VVS2P0_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP22; 
  std::complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (P2[0] * TMP22) + cI * (TMP19 * V2[2])); 
  V1[3] = denom * S3[2] * (-cI * (P2[1] * TMP22) + cI * (TMP19 * V2[3])); 
  V1[4] = denom * S3[2] * (-cI * (P2[2] * TMP22) + cI * (TMP19 * V2[4])); 
  V1[5] = denom * S3[2] * (-cI * (P2[3] * TMP22) + cI * (TMP19 * V2[5])); 
}


void FFS4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP2; 
}


void FFV5_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP13; 
  std::complex<double> TMP8; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP13 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP8 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-4. * cI) * (OM3 * - 1./4. * P3[0] * (TMP8 + 4. * (TMP13)) +
      (+1./4. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * (-4. * cI) * (OM3 * - 1./4. * P3[1] * (TMP8 + 4. * (TMP13)) +
      (-1./4. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * 4. * cI * (OM3 * 1./4. * P3[2] * (TMP8 + 4. * (TMP13)) +
      (+1./4. * cI * (F2[5] * F1[2]) - 1./4. * cI * (F2[4] * F1[3]) - cI *
      (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
  V3[5] = denom * 4. * cI * (OM3 * 1./4. * P3[3] * (TMP8 + 4. * (TMP13)) +
      (+1./4. * (F2[4] * F1[2]) - 1./4. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void FFS5C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (+cI * (TMP1 + TMP2)); 
}


void VVVV5_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  vertex = COUP * 1./2. * (-2. * cI * (TMP27 * TMP29) + cI * (TMP25 * TMP30 +
      TMP21 * TMP35));
}


void VVVVS1P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (-cI * (V1[2] * TMP27) + cI * (V2[2] * TMP25)); 
  V4[3] = denom * S5[2] * (-cI * (V1[3] * TMP27) + cI * (V2[3] * TMP25)); 
  V4[4] = denom * S5[2] * (-cI * (V1[4] * TMP27) + cI * (V2[4] * TMP25)); 
  V4[5] = denom * S5[2] * (-cI * (V1[5] * TMP27) + cI * (V2[5] * TMP25)); 
}


void VVVV3P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP29; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI * (V2[2] * TMP29) + cI * (TMP21 * V4[2])); 
  V3[3] = denom * (-cI * (V2[3] * TMP29) + cI * (TMP21 * V4[3])); 
  V3[4] = denom * (-cI * (V2[4] * TMP29) + cI * (TMP21 * V4[4])); 
  V3[5] = denom * (-cI * (V2[5] * TMP29) + cI * (TMP21 * V4[5])); 
}


void VVVVS2P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (-cI * (V3[2] * TMP29) + cI * (V1[2] * TMP35)); 
  V2[3] = denom * S5[2] * (-cI * (V3[3] * TMP29) + cI * (V1[3] * TMP35)); 
  V2[4] = denom * S5[2] * (-cI * (V3[4] * TMP29) + cI * (V1[4] * TMP35)); 
  V2[5] = denom * S5[2] * (-cI * (V3[5] * TMP29) + cI * (V1[5] * TMP35)); 
}


void FFS1C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * S3[2] * (F1[4] * (P2[0] + P2[3]) + (F1[5] * (P2[1] +
      cI * (P2[2])) - F1[2] * M2));
  F2[3] = denom * cI * S3[2] * (F1[4] * (+cI * (P2[2]) - P2[1]) + (F1[5] *
      (P2[3] - P2[0]) + F1[3] * M2));
  F2[4] = denom * cI * S3[2] * (F1[2] * (P2[3] - P2[0]) + (F1[3] * (P2[1] + cI
      * (P2[2])) + F1[4] * M2));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (+cI * (P2[2]) - P2[1]) + (F1[3] *
      (P2[0] + P2[3]) - F1[5] * M2));
}


void VVVV3_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP24; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (-cI * (TMP17 * TMP35) + cI * (TMP24 * TMP29))
      + (-cI * (V3[2] * TMP29) + cI * (V1[2] * TMP35)));
  V2[3] = denom * (OM2 * P2[1] * (-cI * (TMP17 * TMP35) + cI * (TMP24 * TMP29))
      + (-cI * (V3[3] * TMP29) + cI * (V1[3] * TMP35)));
  V2[4] = denom * (OM2 * P2[2] * (-cI * (TMP17 * TMP35) + cI * (TMP24 * TMP29))
      + (-cI * (V3[4] * TMP29) + cI * (V1[4] * TMP35)));
  V2[5] = denom * (OM2 * P2[3] * (-cI * (TMP17 * TMP35) + cI * (TMP24 * TMP29))
      + (-cI * (V3[5] * TMP29) + cI * (V1[5] * TMP35)));
}


void FFV4C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * (-cI) * (F1[2] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] -
      cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[3] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * (-1.) *
      (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * 2. * (V3[2] + V3[5]) + 2. *
      (F1[5] * (V3[3] + cI * (V3[4]))))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + (F1[3] * (P2[0] * (-1.) * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[4] * 2. * (V3[3] - cI * (V3[4])) + 2. *
      (F1[5] * (V3[2] - V3[5])))));
  F2[4] = denom * (-2. * cI) * (F1[4] * (P2[0] * (-1.) * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + (+1./2. * (M2 * (+2. * (F1[2] * 1./2. * (V3[2] -
      V3[5])) - F1[3] * (V3[3] + cI * (V3[4])))) + F1[5] * (P2[0] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[1] * (V3[2] - V3[5]) + (P2[2] * (-cI *
      (V3[5]) + cI * (V3[2])) + P2[3] * (V3[3] + cI * (V3[4])))))));
  F2[5] = denom * (-2. * cI) * (F1[4] * (P2[0] * (+cI * (V3[4]) - V3[3]) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (+cI * (V3[4]) - V3[3])))) + (+1./2. * (M2 * (F1[3] * (V3[2] +
      V3[5]) + 2. * (F1[2] * 1./2. * (+cI * (V3[4]) - V3[3])))) + F1[5] *
      (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] *
      (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] - V3[2]))))));
}


void FFS2C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * F1[2] * M2 * S3[2]; 
  F2[3] = denom * cI * F1[3] * M2 * S3[2]; 
  F2[4] = denom * cI * S3[2] * (F1[2] * (P2[3] - P2[0]) + F1[3] * (P2[1] + cI *
      (P2[2])));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (+cI * (P2[2]) - P2[1]) + F1[3] *
      (P2[0] + P2[3]));
}

void FFS2_4C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFS2C1_1(F1, S3, COUP1, M2, W2, F2); 
  FFS4C1_1(F1, S3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFS4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * - cI * TMP2 * S3[2]; 
}


void SSSS1_2(std::complex<double> S1[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> S2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  S2[0] = +S1[0] + S3[0] + S4[0]; 
  S2[1] = +S1[1] + S3[1] + S4[1]; 
  P2[0] = -S2[0].real(); 
  P2[1] = -S2[1].real(); 
  P2[2] = -S2[1].imag(); 
  P2[3] = -S2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  S2[2] = denom * cI * S4[2] * S3[2] * S1[2]; 
}


void VVV1_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP20; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP20 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (TMP27 * (-cI * (TMP20) + cI * (TMP19)) + (-cI
      * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 * (-cI * (P2[0]) + cI
      * (P3[0])) + (V2[2] * (-cI * (TMP23) + cI * (TMP24)) + V3[2] * (-cI *
      (TMP26) + cI * (TMP22)))));
  V1[3] = denom * (OM1 * P1[1] * (TMP27 * (-cI * (TMP20) + cI * (TMP19)) + (-cI
      * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 * (-cI * (P2[1]) + cI
      * (P3[1])) + (V2[3] * (-cI * (TMP23) + cI * (TMP24)) + V3[3] * (-cI *
      (TMP26) + cI * (TMP22)))));
  V1[4] = denom * (OM1 * P1[2] * (TMP27 * (-cI * (TMP20) + cI * (TMP19)) + (-cI
      * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 * (-cI * (P2[2]) + cI
      * (P3[2])) + (V2[4] * (-cI * (TMP23) + cI * (TMP24)) + V3[4] * (-cI *
      (TMP26) + cI * (TMP22)))));
  V1[5] = denom * (OM1 * P1[3] * (TMP27 * (-cI * (TMP20) + cI * (TMP19)) + (-cI
      * (TMP22 * TMP24) + cI * (TMP23 * TMP26))) + (TMP27 * (-cI * (P2[3]) + cI
      * (P3[3])) + (V2[5] * (-cI * (TMP23) + cI * (TMP24)) + V3[5] * (-cI *
      (TMP26) + cI * (TMP22)))));
}


void VVS1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * - cI * TMP21 * S3[2]; 
}


void FFS2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * - cI * TMP1 * S3[2]; 
}

void FFS2_4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFS2_0(F1, F2, S3, COUP1, vertex); 
  FFS4_0(F1, F2, S3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFV1C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * (-cI) * (F1[2] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] -
      cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[3] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * (-1.) *
      (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * (V3[2] + V3[5]) + F1[5] *
      (V3[3] + cI * (V3[4])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) + P2[3] *
      (V3[3] - cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (-1.) * (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[4] * (+cI * (V3[4]) - V3[3]) + F1[5] *
      (V3[5] - V3[2]))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (+cI *
      (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI * (V3[3])) - P2[3] *
      (V3[2] + V3[5])))) + (F1[5] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[2] * (V3[5] - V3[2]) + F1[3] *
      (V3[3] + cI * (V3[4])))));
  F2[5] = denom * (-cI) * (F1[4] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + (F1[5] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[2] * (+cI * (V3[4]) - V3[3]) + F1[3] * (V3[2] +
      V3[5]))));
}


void VVSS1P0_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + S3[0] + S4[0]; 
  V1[1] = +V2[1] + S3[1] + S4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * - cI * V2[2] * S4[2] * S3[2]; 
  V1[3] = denom * - cI * V2[3] * S4[2] * S3[2]; 
  V1[4] = denom * - cI * V2[4] * S4[2] * S3[2]; 
  V1[5] = denom * - cI * V2[5] * S4[2] * S3[2]; 
}


void VVVVS2P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (-cI * (TMP27 * V4[2]) + cI * (V2[2] * TMP35)); 
  V1[3] = denom * S5[2] * (-cI * (TMP27 * V4[3]) + cI * (V2[3] * TMP35)); 
  V1[4] = denom * S5[2] * (-cI * (TMP27 * V4[4]) + cI * (V2[4] * TMP35)); 
  V1[5] = denom * S5[2] * (-cI * (TMP27 * V4[5]) + cI * (V2[5] * TMP35)); 
}


void VSS1P0_1(std::complex<double> S2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> denom; 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  V1[0] = +S2[0] + S3[0]; 
  V1[1] = +S2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S2[2] * S3[2] * (-cI * (P2[0]) + cI * (P3[0])); 
  V1[3] = denom * S2[2] * S3[2] * (-cI * (P2[1]) + cI * (P3[1])); 
  V1[4] = denom * S2[2] * S3[2] * (-cI * (P2[2]) + cI * (P3[2])); 
  V1[5] = denom * S2[2] * S3[2] * (-cI * (P2[3]) + cI * (P3[3])); 
}


void VVSS1_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> S4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> denom; 
  S4[0] = +V1[0] + V2[0] + S3[0]; 
  S4[1] = +V1[1] + V2[1] + S3[1]; 
  P4[0] = -S4[0].real(); 
  P4[1] = -S4[1].real(); 
  P4[2] = -S4[1].imag(); 
  P4[3] = -S4[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  S4[2] = denom * cI * TMP21 * S3[2]; 
}


void FFV3C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP10; 
  std::complex<double> TMP16; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP10 = (-1.) * (F1[2] * (F2[4] * (P3[0] - P3[3]) + F2[5] * (+cI * (P3[2]) -
      P3[1])) + F1[3] * (F2[4] * (-1.) * (P3[1] + cI * (P3[2])) + F2[5] *
      (P3[0] + P3[3])));
  TMP16 = (-1.) * (F1[4] * (F2[2] * (P3[0] + P3[3]) + F2[3] * (P3[1] - cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] + cI * (P3[2])) + F2[3] * (P3[0] -
      P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-2. * cI) * (OM3 * 1./2. * P3[0] * (+2. * (TMP16) - TMP10) +
      (-1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * 2. * cI * (OM3 * 1./2. * P3[1] * (TMP10 - 2. * (TMP16)) +
      (+1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * (-2. * cI) * (OM3 * 1./2. * P3[2] * (+2. * (TMP16) - TMP10) +
      (+1./2. * cI * (F2[5] * F1[2]) - 1./2. * cI * (F2[4] * F1[3]) - cI *
      (F2[2] * F1[5]) + cI * (F2[3] * F1[4])));
  V3[5] = denom * (-2. * cI) * (OM3 * 1./2. * P3[3] * (+2. * (TMP16) - TMP10) +
      (-1./2. * (F2[4] * F1[2]) + 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void VVVVS2_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP27; 
  std::complex<double> TMP31; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (OM1 * P1[0] * (-cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[2]) + cI * (V2[2] * TMP35)));
  V1[3] = denom * S5[2] * (OM1 * P1[1] * (-cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[3]) + cI * (V2[3] * TMP35)));
  V1[4] = denom * S5[2] * (OM1 * P1[2] * (-cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[4]) + cI * (V2[4] * TMP35)));
  V1[5] = denom * S5[2] * (OM1 * P1[3] * (-cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[5]) + cI * (V2[5] * TMP35)));
}


void FFS5_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) - F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1));
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) + F2[4] * M1));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) - F2[5] * M1));
}


void FFS4C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * - cI * TMP2 * S3[2]; 
}


void VVVV4P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (-cI * (V2[2] * TMP25) + cI * (V3[2] * TMP21)); 
  V4[3] = denom * (-cI * (V2[3] * TMP25) + cI * (V3[3] * TMP21)); 
  V4[4] = denom * (-cI * (V2[4] * TMP25) + cI * (V3[4] * TMP21)); 
  V4[5] = denom * (-cI * (V2[5] * TMP25) + cI * (V3[5] * TMP21)); 
}


void VVS2_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP19; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP19 * TMP21) + cI * (TMP17 * TMP22)); 
}


void FFV1_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * (-cI) * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * (-1.) *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * (-1.) * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) - P2[3] * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI *
      (V3[4])) + F1[3] * (V3[2] - V3[5]))));
}


void SSS1_3(std::complex<double> S1[], std::complex<double> S2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> denom; 
  S3[0] = +S1[0] + S2[0]; 
  S3[1] = +S1[1] + S2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * S2[2] * S1[2]; 
}


void VVVVS3_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (OM1 * P1[0] * (-cI * (TMP22 * TMP35) + cI * (TMP23 *
      TMP30)) + (-cI * (V3[2] * TMP30) + cI * (V2[2] * TMP35)));
  V1[3] = denom * S5[2] * (OM1 * P1[1] * (-cI * (TMP22 * TMP35) + cI * (TMP23 *
      TMP30)) + (-cI * (V3[3] * TMP30) + cI * (V2[3] * TMP35)));
  V1[4] = denom * S5[2] * (OM1 * P1[2] * (-cI * (TMP22 * TMP35) + cI * (TMP23 *
      TMP30)) + (-cI * (V3[4] * TMP30) + cI * (V2[4] * TMP35)));
  V1[5] = denom * S5[2] * (OM1 * P1[3] * (-cI * (TMP22 * TMP35) + cI * (TMP23 *
      TMP30)) + (-cI * (V3[5] * TMP30) + cI * (V2[5] * TMP35)));
}


void VVVV2_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  vertex = COUP * (-1.) * (-2. * cI * (TMP21 * TMP35) + cI * (TMP27 * TMP29 +
      TMP25 * TMP30));
}


void VVSS1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> denom; 
  S3[0] = +V1[0] + V2[0] + S4[0]; 
  S3[1] = +V1[1] + V2[1] + S4[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP21 * S4[2]; 
}


void VVVV5P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * 1./2. * (-2. * cI * (TMP27 * V4[2]) + cI * (V3[2] * TMP30 +
      V2[2] * TMP35));
  V1[3] = denom * 1./2. * (-2. * cI * (TMP27 * V4[3]) + cI * (V3[3] * TMP30 +
      V2[3] * TMP35));
  V1[4] = denom * 1./2. * (-2. * cI * (TMP27 * V4[4]) + cI * (V3[4] * TMP30 +
      V2[4] * TMP35));
  V1[5] = denom * 1./2. * (-2. * cI * (TMP27 * V4[5]) + cI * (V3[5] * TMP30 +
      V2[5] * TMP35));
}


void VSS1_0(std::complex<double> V1[], std::complex<double> S2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  vertex = COUP * S2[2] * S3[2] * (-cI * (TMP17) + cI * (TMP18)); 
}


void FFS3C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] + P2[3]) + (F1[5] * (P2[1] + cI
      * (P2[2])) + F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (+cI * (P2[2]) - P2[1]) + (F1[5] *
      (P2[3] - P2[0]) - F1[3] * M2));
  F2[4] = denom * cI * S3[2] * (F1[2] * (P2[3] - P2[0]) + (F1[3] * (P2[1] + cI
      * (P2[2])) - F1[4] * M2));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (+cI * (P2[2]) - P2[1]) + (F1[3] *
      (P2[0] + P2[3]) + F1[5] * M2));
}


void FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP7; 
  TMP7 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP7; 
}

void FFV2_5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2_0(F1, F2, V3, COUP1, vertex); 
  FFV5_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV2_3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2_0(F1, F2, V3, COUP1, vertex); 
  FFV3_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV2_4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2_0(F1, F2, V3, COUP1, vertex); 
  FFV4_0(F1, F2, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVVV1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * (-cI * (TMP27 * TMP29) + cI * (TMP25 * TMP30)); 
}


void FFS1_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) + F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] *
      (-1.) * (P2[0] + P2[3]) - F1[3] * M2));
  F2[4] = denom * - cI * S3[2] * (F1[2] * (-1.) * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) - F1[4] * M2));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) + F1[5] * M2));
}


void FFV2C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] - V3[5]) + (P2[1] * (+cI *
      (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI * (V3[3])) + P2[3] *
      (V3[2] - V3[5])))) + F1[3] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + F1[3] * (P2[0] * (-1.) * (V3[2] + V3[5]) + (P2[1]
      * (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] *
      (V3[2] + V3[5])))));
  F2[4] = denom * cI * M2 * (F1[2] * (V3[5] - V3[2]) + F1[3] * (V3[3] + cI *
      (V3[4])));
  F2[5] = denom * - cI * M2 * (F1[2] * (+cI * (V3[4]) - V3[3]) + F1[3] * (V3[2]
      + V3[5]));
}

void FFV2_5C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_1(F1, V3, COUP1, M2, W2, F2); 
  FFV5C1_1(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_3C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_1(F1, V3, COUP1, M2, W2, F2); 
  FFV3C1_1(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_4C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_1(F1, V3, COUP1, M2, W2, F2); 
  FFV4C1_1(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFS5_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + (F1[5] * (+cI *
      (P2[2]) - P2[1]) + F1[2] * M2));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) + (F1[5] *
      (-1.) * (P2[0] + P2[3]) - F1[3] * M2));
  F2[4] = denom * - cI * S3[2] * (F1[2] * (-1.) * (P2[0] + P2[3]) + (F1[3] *
      (+cI * (P2[2]) - P2[1]) - F1[4] * M2));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + (F1[3] *
      (P2[0] - P2[3]) + F1[5] * M2));
}


void VVVV4P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V3[2] * TMP30) + cI * (V2[2] * TMP35)); 
  V1[3] = denom * (-cI * (V3[3] * TMP30) + cI * (V2[3] * TMP35)); 
  V1[4] = denom * (-cI * (V3[4] * TMP30) + cI * (V2[4] * TMP35)); 
  V1[5] = denom * (-cI * (V3[5] * TMP30) + cI * (V2[5] * TMP35)); 
}


void VVVV4_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP25; 
  std::complex<double> TMP32; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (-cI * (TMP17 * TMP35) + cI * (TMP25 * TMP32))
      + (-cI * (TMP25 * V4[2]) + cI * (V1[2] * TMP35)));
  V2[3] = denom * (OM2 * P2[1] * (-cI * (TMP17 * TMP35) + cI * (TMP25 * TMP32))
      + (-cI * (TMP25 * V4[3]) + cI * (V1[3] * TMP35)));
  V2[4] = denom * (OM2 * P2[2] * (-cI * (TMP17 * TMP35) + cI * (TMP25 * TMP32))
      + (-cI * (TMP25 * V4[4]) + cI * (V1[4] * TMP35)));
  V2[5] = denom * (OM2 * P2[3] * (-cI * (TMP17 * TMP35) + cI * (TMP25 * TMP32))
      + (-cI * (TMP25 * V4[5]) + cI * (V1[5] * TMP35)));
}


void FFV4_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * 2. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3]
      - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] -
      V3[2])))) + (+1./2. * (M1 * (F2[5] * (V3[3] + cI * (V3[4])) + 2. * (F2[4]
      * 1./2. * (V3[2] + V3[5])))) + F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) +
      (P1[1] * (-1.) * (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] +
      V3[5])) + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * 2. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (+1./2. * (M1 * (F2[5] * (V3[2] - V3[5]) + 2. *
      (F2[4] * 1./2. * (V3[3] - cI * (V3[4]))))) + F2[3] * (P1[0] * (-1.) *
      (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI
      * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 2. * (V3[5] - V3[2]) + 2. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * 2. * (+cI * (V3[4]) - V3[3]) + 2. * (F2[3] *
      (V3[2] + V3[5])))));
}


void VVVS1P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> S4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V2[0] = +V1[0] + V3[0] + S4[0]; 
  V2[1] = +V1[1] + V3[1] + S4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S4[2] * (TMP25 * (-cI * (P3[0]) + cI * (P1[0])) + (V1[2] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[2] * (-cI * (TMP17) + cI * (TMP18))));
  V2[3] = denom * S4[2] * (TMP25 * (-cI * (P3[1]) + cI * (P1[1])) + (V1[3] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[3] * (-cI * (TMP17) + cI * (TMP18))));
  V2[4] = denom * S4[2] * (TMP25 * (-cI * (P3[2]) + cI * (P1[2])) + (V1[4] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[4] * (-cI * (TMP17) + cI * (TMP18))));
  V2[5] = denom * S4[2] * (TMP25 * (-cI * (P3[3]) + cI * (P1[3])) + (V1[5] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[5] * (-cI * (TMP17) + cI * (TMP18))));
}


void FFV2P0_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (-1.) * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + F1[3] * (P2[0] * (V3[2] - V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) + P2[3] * (V3[2] - V3[5])))));
  F2[4] = denom * - cI * M2 * (F1[2] * (-1.) * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}


void FFS4P0_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + F1[5] * (+cI *
      (P2[2]) - P2[1]));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) - F1[5] *
      (P2[0] + P2[3]));
  F2[4] = denom * cI * F1[4] * M2 * S3[2]; 
  F2[5] = denom * cI * F1[5] * M2 * S3[2]; 
}


void VVV1P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP26; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V3[0] = +V1[0] + V2[0]; 
  V3[1] = +V1[1] + V2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (TMP21 * (-cI * (P1[0]) + cI * (P2[0])) + (V1[2] * (-cI *
      (TMP26) + cI * (TMP22)) + V2[2] * (-cI * (TMP17) + cI * (TMP18))));
  V3[3] = denom * (TMP21 * (-cI * (P1[1]) + cI * (P2[1])) + (V1[3] * (-cI *
      (TMP26) + cI * (TMP22)) + V2[3] * (-cI * (TMP17) + cI * (TMP18))));
  V3[4] = denom * (TMP21 * (-cI * (P1[2]) + cI * (P2[2])) + (V1[4] * (-cI *
      (TMP26) + cI * (TMP22)) + V2[4] * (-cI * (TMP17) + cI * (TMP18))));
  V3[5] = denom * (TMP21 * (-cI * (P1[3]) + cI * (P2[3])) + (V1[5] * (-cI *
      (TMP26) + cI * (TMP22)) + V2[5] * (-cI * (TMP17) + cI * (TMP18))));
}


void VVVVS1_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP34; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (OM4 * P4[0] * (-cI * (TMP25 * TMP34) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[2] * TMP27) + cI * (V2[2] * TMP25)));
  V4[3] = denom * S5[2] * (OM4 * P4[1] * (-cI * (TMP25 * TMP34) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[3] * TMP27) + cI * (V2[3] * TMP25)));
  V4[4] = denom * S5[2] * (OM4 * P4[2] * (-cI * (TMP25 * TMP34) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[4] * TMP27) + cI * (V2[4] * TMP25)));
  V4[5] = denom * S5[2] * (OM4 * P4[3] * (-cI * (TMP25 * TMP34) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[5] * TMP27) + cI * (V2[5] * TMP25)));
}


void VVVS1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  vertex = COUP * S4[2] * (TMP21 * (-cI * (TMP23) + cI * (TMP24)) + (TMP25 *
      (-cI * (TMP26) + cI * (TMP22)) + TMP27 * (-cI * (TMP17) + cI * (TMP18))));
}


void FFV3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP11; 
  std::complex<double> TMP12; 
  TMP11 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP12 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  vertex = COUP * (-cI * (TMP11) + 2. * cI * (TMP12)); 
}


void VVVV2_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (-2. * cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29 + TMP18 * TMP30)) + (-cI * (V2[2] * TMP29 + V1[2] * TMP30) + 2. *
      cI * (TMP21 * V4[2])));
  V3[3] = denom * (OM3 * P3[1] * (-2. * cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29 + TMP18 * TMP30)) + (-cI * (V2[3] * TMP29 + V1[3] * TMP30) + 2. *
      cI * (TMP21 * V4[3])));
  V3[4] = denom * (OM3 * P3[2] * (-2. * cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29 + TMP18 * TMP30)) + (-cI * (V2[4] * TMP29 + V1[4] * TMP30) + 2. *
      cI * (TMP21 * V4[4])));
  V3[5] = denom * (OM3 * P3[3] * (-2. * cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29 + TMP18 * TMP30)) + (-cI * (V2[5] * TMP29 + V1[5] * TMP30) + 2. *
      cI * (TMP21 * V4[5])));
}


void VVVV3_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (-cI * (TMP21 * TMP36) + cI * (TMP26 * TMP29))
      + (-cI * (V2[2] * TMP29) + cI * (TMP21 * V4[2])));
  V3[3] = denom * (OM3 * P3[1] * (-cI * (TMP21 * TMP36) + cI * (TMP26 * TMP29))
      + (-cI * (V2[3] * TMP29) + cI * (TMP21 * V4[3])));
  V3[4] = denom * (OM3 * P3[2] * (-cI * (TMP21 * TMP36) + cI * (TMP26 * TMP29))
      + (-cI * (V2[4] * TMP29) + cI * (TMP21 * V4[4])));
  V3[5] = denom * (OM3 * P3[3] * (-cI * (TMP21 * TMP36) + cI * (TMP26 * TMP29))
      + (-cI * (V2[5] * TMP29) + cI * (TMP21 * V4[5])));
}


void FFS3_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] + cI
      * (P1[2])) + F2[2] * M1));
  F1[3] = denom * - cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) - F2[3] * M1));
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) - F2[4] * M1));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) + F2[5] * M1));
}


void FFS3C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] - P1[3]) + (F2[5] * (+cI *
      (P1[2]) - P1[1]) - F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (P1[1] + cI * (P1[2])) + (F2[5] * (-1.)
      * (P1[0] + P1[3]) + F2[3] * M1));
  F1[4] = denom * - cI * S3[2] * (F2[2] * (-1.) * (P1[0] + P1[3]) + (F2[3] *
      (+cI * (P1[2]) - P1[1]) + F2[4] * M1));
  F1[5] = denom * cI * S3[2] * (F2[2] * (P1[1] + cI * (P1[2])) + (F2[3] *
      (P1[0] - P1[3]) - F2[5] * M1));
}


void FFV5_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 4. * (V3[2] - V3[5]) + 4. * (F1[5]
      * (+cI * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (-1.) * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) + P2[3] * (V3[2] - V3[5])))) + M2 * (F1[4] * (-4.) * (V3[3] + cI
      * (V3[4])) + 4. * (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * (-4. * cI) * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + (+1./4. * (M2 * (F1[3] * (+cI * (V3[4]) - V3[3]) + 4. *
      (F1[2] * (-1./4.) * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI *
      (V3[4])) + (P2[1] * (-1.) * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] +
      V3[5])) + P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * (-4. * cI) * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) -
      P2[3] * (V3[3] + cI * (V3[4]))))) + (+1./4. * (M2 * (F1[3] * (V3[5] -
      V3[2]) + 4. * (F1[2] * (-1./4.) * (V3[3] + cI * (V3[4]))))) + F1[5] *
      (P2[0] * (-1.) * (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) +
      (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void FFS1_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + (F2[5] * (P1[1] +
      cI * (P1[2])) - F2[2] * M1));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + (F2[5] *
      (P1[3] - P1[0]) + F2[3] * M1));
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + (F2[3] * (P1[1] + cI
      * (P1[2])) + F2[4] * M1));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + (F2[3] *
      (P1[0] + P1[3]) - F2[5] * M1));
}


void FFS5C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * S3[2] * (F2[4] * (P1[0] - P1[3]) + (F2[5] * (+cI *
      (P1[2]) - P1[1]) + F2[2] * M1));
  F1[3] = denom * - cI * S3[2] * (F2[4] * (P1[1] + cI * (P1[2])) + (F2[5] *
      (-1.) * (P1[0] + P1[3]) - F2[3] * M1));
  F1[4] = denom * - cI * S3[2] * (F2[2] * (-1.) * (P1[0] + P1[3]) + (F2[3] *
      (+cI * (P1[2]) - P1[1]) - F2[4] * M1));
  F1[5] = denom * cI * S3[2] * (F2[2] * (P1[1] + cI * (P1[2])) + (F2[3] *
      (P1[0] - P1[3]) + F2[5] * M1));
}


void VVVV5_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * 1./2. * (OM3 * - P3[0] * (-2. * cI * (TMP26 * TMP29) + cI *
      (TMP18 * TMP30 + TMP21 * TMP36)) + (-2. * cI * (V2[2] * TMP29) + cI *
      (V1[2] * TMP30 + TMP21 * V4[2])));
  V3[3] = denom * 1./2. * (OM3 * - P3[1] * (-2. * cI * (TMP26 * TMP29) + cI *
      (TMP18 * TMP30 + TMP21 * TMP36)) + (-2. * cI * (V2[3] * TMP29) + cI *
      (V1[3] * TMP30 + TMP21 * V4[3])));
  V3[4] = denom * 1./2. * (OM3 * - P3[2] * (-2. * cI * (TMP26 * TMP29) + cI *
      (TMP18 * TMP30 + TMP21 * TMP36)) + (-2. * cI * (V2[4] * TMP29) + cI *
      (V1[4] * TMP30 + TMP21 * V4[4])));
  V3[5] = denom * 1./2. * (OM3 * - P3[3] * (-2. * cI * (TMP26 * TMP29) + cI *
      (TMP18 * TMP30 + TMP21 * TMP36)) + (-2. * cI * (V2[5] * TMP29) + cI *
      (V1[5] * TMP30 + TMP21 * V4[5])));
}


void VVVV4P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (-cI * (TMP25 * V4[2]) + cI * (V1[2] * TMP35)); 
  V2[3] = denom * (-cI * (TMP25 * V4[3]) + cI * (V1[3] * TMP35)); 
  V2[4] = denom * (-cI * (TMP25 * V4[4]) + cI * (V1[4] * TMP35)); 
  V2[5] = denom * (-cI * (TMP25 * V4[5]) + cI * (V1[5] * TMP35)); 
}


void VVVV3P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (-cI * (V3[2] * TMP29) + cI * (V1[2] * TMP35)); 
  V2[3] = denom * (-cI * (V3[3] * TMP29) + cI * (V1[3] * TMP35)); 
  V2[4] = denom * (-cI * (V3[4] * TMP29) + cI * (V1[4] * TMP35)); 
  V2[5] = denom * (-cI * (V3[5] * TMP29) + cI * (V1[5] * TMP35)); 
}


void FFS1C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * S3[2] * (F2[4] * (P1[0] - P1[3]) + (F2[5] * (+cI *
      (P1[2]) - P1[1]) + F2[2] * M1));
  F1[3] = denom * - cI * S3[2] * (F2[4] * (P1[1] + cI * (P1[2])) + (F2[5] *
      (-1.) * (P1[0] + P1[3]) - F2[3] * M1));
  F1[4] = denom * - cI * S3[2] * (F2[2] * (-1.) * (P1[0] + P1[3]) + (F2[3] *
      (+cI * (P1[2]) - P1[1]) - F2[4] * M1));
  F1[5] = denom * cI * S3[2] * (F2[2] * (P1[1] + cI * (P1[2])) + (F2[3] *
      (P1[0] - P1[3]) + F2[5] * M1));
}


void VVVVS1_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP23; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP31; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (OM1 * P1[0] * (-cI * (TMP23 * TMP30) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[2]) + cI * (V3[2] * TMP30)));
  V1[3] = denom * S5[2] * (OM1 * P1[1] * (-cI * (TMP23 * TMP30) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[3]) + cI * (V3[3] * TMP30)));
  V1[4] = denom * S5[2] * (OM1 * P1[2] * (-cI * (TMP23 * TMP30) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[4]) + cI * (V3[4] * TMP30)));
  V1[5] = denom * S5[2] * (OM1 * P1[3] * (-cI * (TMP23 * TMP30) + cI * (TMP27 *
      TMP31)) + (-cI * (TMP27 * V4[5]) + cI * (V3[5] * TMP30)));
}


void VVVVS3P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (-cI * (V1[2] * TMP30) + cI * (TMP21 * V4[2])); 
  V3[3] = denom * S5[2] * (-cI * (V1[3] * TMP30) + cI * (TMP21 * V4[3])); 
  V3[4] = denom * S5[2] * (-cI * (V1[4] * TMP30) + cI * (TMP21 * V4[4])); 
  V3[5] = denom * S5[2] * (-cI * (V1[5] * TMP30) + cI * (TMP21 * V4[5])); 
}


void FFV3_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP13; 
  std::complex<double> TMP8; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP13 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP8 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (OM3 * 1./2. * P3[0] * (TMP8 - 2. * (TMP13)) +
      (-1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * 2. * cI * (OM3 * 1./2. * P3[1] * (TMP8 - 2. * (TMP13)) +
      (+1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * (-2. * cI) * (OM3 * 1./2. * P3[2] * (+2. * (TMP13) - TMP8) +
      (-1./2. * cI * (F2[5] * F1[2]) + 1./2. * cI * (F2[4] * F1[3]) - cI *
      (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
  V3[5] = denom * (-2. * cI) * (OM3 * 1./2. * P3[3] * (+2. * (TMP13) - TMP8) +
      (-1./2. * (F2[4] * F1[2]) + 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void FFV4C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP14; 
  std::complex<double> TMP15; 
  TMP15 = (-1.) * (F1[4] * (F2[2] * (V3[2] + V3[5]) + F2[3] * (V3[3] - cI *
      (V3[4]))) + F1[5] * (F2[2] * (V3[3] + cI * (V3[4])) + F2[3] * (V3[2] -
      V3[5])));
  TMP14 = (-1.) * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + F1[3] * (F2[4] * (-1.) * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])));
  vertex = COUP * (-1.) * (+cI * (TMP14) + 2. * cI * (TMP15)); 
}


void VVVV3_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (OM4 * P4[0] * (-cI * (TMP21 * TMP37) + cI * (TMP27 * TMP33))
      + (-cI * (V1[2] * TMP27) + cI * (V3[2] * TMP21)));
  V4[3] = denom * (OM4 * P4[1] * (-cI * (TMP21 * TMP37) + cI * (TMP27 * TMP33))
      + (-cI * (V1[3] * TMP27) + cI * (V3[3] * TMP21)));
  V4[4] = denom * (OM4 * P4[2] * (-cI * (TMP21 * TMP37) + cI * (TMP27 * TMP33))
      + (-cI * (V1[4] * TMP27) + cI * (V3[4] * TMP21)));
  V4[5] = denom * (OM4 * P4[3] * (-cI * (TMP21 * TMP37) + cI * (TMP27 * TMP33))
      + (-cI * (V1[5] * TMP27) + cI * (V3[5] * TMP21)));
}


void FFS2C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * - cI * TMP1 * S3[2]; 
}

void FFS2_4C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFS2C1_0(F2, F1, S3, COUP1, vertex); 
  FFS4C1_0(F2, F1, S3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void VVVVS2_5(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, double M5, double W5, std::complex<double> S5[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P5[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P5[0] * P5[0]) - (P5[1] * P5[1]) - (P5[2] * P5[2]) - (P5[3] *
      P5[3]) - M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP21 * TMP35) + cI * (TMP27 * TMP29)); 
}


void FFS4_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI
      * (P1[2])));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + F2[5] *
      (P1[3] - P1[0]));
  F1[4] = denom * cI * F2[4] * M1 * S3[2]; 
  F1[5] = denom * cI * F2[5] * M1 * S3[2]; 
}


void SSSS1_1(std::complex<double> S2[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> S1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  S1[0] = +S2[0] + S3[0] + S4[0]; 
  S1[1] = +S2[1] + S3[1] + S4[1]; 
  P1[0] = -S1[0].real(); 
  P1[1] = -S1[1].real(); 
  P1[2] = -S1[1].imag(); 
  P1[3] = -S1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  S1[2] = denom * cI * S4[2] * S3[2] * S2[2]; 
}


void VVV1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  vertex = COUP * (TMP21 * (-cI * (TMP23) + cI * (TMP24)) + (TMP25 * (-cI *
      (TMP26) + cI * (TMP22)) + TMP27 * (-cI * (TMP17) + cI * (TMP18))));
}


void VVS1_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (V2[2]) + cI * (P1[0] * OM1 * TMP22)); 
  V1[3] = denom * S3[2] * (-cI * (V2[3]) + cI * (P1[1] * OM1 * TMP22)); 
  V1[4] = denom * S3[2] * (-cI * (V2[4]) + cI * (P1[2] * OM1 * TMP22)); 
  V1[5] = denom * S3[2] * (-cI * (V2[5]) + cI * (P1[3] * OM1 * TMP22)); 
}


void FFS2_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * F2[2] * M1 * S3[2]; 
  F1[3] = denom * cI * F2[3] * M1 * S3[2]; 
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + F2[3] * (P1[1] + cI *
      (P1[2])));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] *
      (P1[0] + P1[3]));
}

void FFS2_4_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFS2_1(F2, S3, COUP1, M1, W1, F1); 
  FFS4_1(F2, S3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void FFS5C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * S3[2] * (F1[4] * (P2[0] + P2[3]) + (F1[5] * (P2[1] +
      cI * (P2[2])) - F1[2] * M2));
  F2[3] = denom * cI * S3[2] * (F1[4] * (+cI * (P2[2]) - P2[1]) + (F1[5] *
      (P2[3] - P2[0]) + F1[3] * M2));
  F2[4] = denom * cI * S3[2] * (F1[2] * (P2[3] - P2[0]) + (F1[3] * (P2[1] + cI
      * (P2[2])) + F1[4] * M2));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (+cI * (P2[2]) - P2[1]) + (F1[3] *
      (P2[0] + P2[3]) - F1[5] * M2));
}


void FFV1C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP5; 
  TMP5 = (-1.) * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + (F1[3] * (F2[4] * (-1.) * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])) + (F1[4] * (F2[2] * (V3[2] + V3[5]) + F2[3] * (V3[3] -
      cI * (V3[4]))) + F1[5] * (F2[2] * (V3[3] + cI * (V3[4])) + F2[3] * (V3[2]
      - V3[5])))));
  vertex = COUP * - cI * TMP5; 
}


void VVS1P0_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * - cI * V2[2] * S3[2]; 
  V1[3] = denom * - cI * V2[3] * S3[2]; 
  V1[4] = denom * - cI * V2[4] * S3[2]; 
  V1[5] = denom * - cI * V2[5] * S3[2]; 
}


void VVVV3P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP27 * V4[2]) + cI * (V2[2] * TMP35)); 
  V1[3] = denom * (-cI * (TMP27 * V4[3]) + cI * (V2[3] * TMP35)); 
  V1[4] = denom * (-cI * (TMP27 * V4[4]) + cI * (V2[4] * TMP35)); 
  V1[5] = denom * (-cI * (TMP27 * V4[5]) + cI * (V2[5] * TMP35)); 
}


void VVVS1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> S4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0] + S4[0]; 
  V1[1] = +V2[1] + V3[1] + S4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S4[2] * (TMP27 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[2] * (-cI * (TMP26) + cI * (TMP22))));
  V1[3] = denom * S4[2] * (TMP27 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[3] * (-cI * (TMP26) + cI * (TMP22))));
  V1[4] = denom * S4[2] * (TMP27 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[4] * (-cI * (TMP26) + cI * (TMP22))));
  V1[5] = denom * S4[2] * (TMP27 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] *
      (-cI * (TMP23) + cI * (TMP24)) + V3[5] * (-cI * (TMP26) + cI * (TMP22))));
}


void VVVV2_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP34; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (OM4 * P4[0] * (-2. * cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33 + TMP25 * TMP34)) + (-cI * (V1[2] * TMP27 + V2[2] * TMP25) + 2. *
      cI * (V3[2] * TMP21)));
  V4[3] = denom * (OM4 * P4[1] * (-2. * cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33 + TMP25 * TMP34)) + (-cI * (V1[3] * TMP27 + V2[3] * TMP25) + 2. *
      cI * (V3[3] * TMP21)));
  V4[4] = denom * (OM4 * P4[2] * (-2. * cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33 + TMP25 * TMP34)) + (-cI * (V1[4] * TMP27 + V2[4] * TMP25) + 2. *
      cI * (V3[4] * TMP21)));
  V4[5] = denom * (OM4 * P4[3] * (-2. * cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33 + TMP25 * TMP34)) + (-cI * (V1[5] * TMP27 + V2[5] * TMP25) + 2. *
      cI * (V3[5] * TMP21)));
}


void FFS2C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP1; 
}

void FFS2_4C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> S3[])
{
  std::complex<double> Stmp[3]; 
  int i; 
  FFS2C1_3(F2, F1, COUP1, M3, W3, S3); 
  FFS4C1_3(F2, F1, COUP2, M3, W3, Stmp); 
  i = 2; 
  while (i < 3)
  {
    S3[i] = S3[i] + Stmp[i]; 
    i++; 
  }
}

void FFV5C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP14; 
  std::complex<double> TMP15; 
  TMP15 = (-1.) * (F1[4] * (F2[2] * (V3[2] + V3[5]) + F2[3] * (V3[3] - cI *
      (V3[4]))) + F1[5] * (F2[2] * (V3[3] + cI * (V3[4])) + F2[3] * (V3[2] -
      V3[5])));
  TMP14 = (-1.) * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + F1[3] * (F2[4] * (-1.) * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])));
  vertex = COUP * (-1.) * (+cI * (TMP14) + 4. * cI * (TMP15)); 
}


void VVVV1_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP34; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (OM4 * P4[0] * (-cI * (TMP25 * TMP34) + cI * (TMP27 * TMP33))
      + (-cI * (V1[2] * TMP27) + cI * (V2[2] * TMP25)));
  V4[3] = denom * (OM4 * P4[1] * (-cI * (TMP25 * TMP34) + cI * (TMP27 * TMP33))
      + (-cI * (V1[3] * TMP27) + cI * (V2[3] * TMP25)));
  V4[4] = denom * (OM4 * P4[2] * (-cI * (TMP25 * TMP34) + cI * (TMP27 * TMP33))
      + (-cI * (V1[4] * TMP27) + cI * (V2[4] * TMP25)));
  V4[5] = denom * (OM4 * P4[3] * (-cI * (TMP25 * TMP34) + cI * (TMP27 * TMP33))
      + (-cI * (V1[5] * TMP27) + cI * (V2[5] * TMP25)));
}


void VVVVS2_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    S5[], std::complex<double> COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  vertex = COUP * S5[2] * (-cI * (TMP27 * TMP29) + cI * (TMP21 * TMP35)); 
}


void VVV1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP20; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP26; 
  std::complex<double> TMP28; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0]; 
  V3[1] = +V1[1] + V2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP28 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP20 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (TMP21 * (-cI * (TMP28) + cI * (TMP20)) + (-cI
      * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 * (-cI * (P1[0]) + cI
      * (P2[0])) + (V1[2] * (-cI * (TMP26) + cI * (TMP22)) + V2[2] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V3[3] = denom * (OM3 * P3[1] * (TMP21 * (-cI * (TMP28) + cI * (TMP20)) + (-cI
      * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 * (-cI * (P1[1]) + cI
      * (P2[1])) + (V1[3] * (-cI * (TMP26) + cI * (TMP22)) + V2[3] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V3[4] = denom * (OM3 * P3[2] * (TMP21 * (-cI * (TMP28) + cI * (TMP20)) + (-cI
      * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 * (-cI * (P1[2]) + cI
      * (P2[2])) + (V1[4] * (-cI * (TMP26) + cI * (TMP22)) + V2[4] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V3[5] = denom * (OM3 * P3[3] * (TMP21 * (-cI * (TMP28) + cI * (TMP20)) + (-cI
      * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 * (-cI * (P1[3]) + cI
      * (P2[3])) + (V1[5] * (-cI * (TMP26) + cI * (TMP22)) + V2[5] * (-cI *
      (TMP17) + cI * (TMP18)))));
}


void VVS1_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + S3[0]; 
  V2[1] = +V1[1] + S3[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S3[2] * (-cI * (V1[2]) + cI * (P2[0] * TMP17 * OM2)); 
  V2[3] = denom * S3[2] * (-cI * (V1[3]) + cI * (P2[1] * TMP17 * OM2)); 
  V2[4] = denom * S3[2] * (-cI * (V1[4]) + cI * (P2[2] * TMP17 * OM2)); 
  V2[5] = denom * S3[2] * (-cI * (V1[5]) + cI * (P2[3] * TMP17 * OM2)); 
}


void FFS4C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP2; 
}


void FFS2_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * F1[2] * M2 * S3[2]; 
  F2[3] = denom * cI * F1[3] * M2 * S3[2]; 
  F2[4] = denom * - cI * S3[2] * (F1[2] * (-1.) * (P2[0] + P2[3]) + F1[3] *
      (+cI * (P2[2]) - P2[1]));
  F2[5] = denom * cI * S3[2] * (F1[2] * (P2[1] + cI * (P2[2])) + F1[3] * (P2[0]
      - P2[3]));
}

void FFS2_4_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFS2_2(F1, S3, COUP1, M2, W2, F2); 
  FFS4_2(F1, S3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void VVS2_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP19; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  V2[0] = +V1[0] + S3[0]; 
  V2[1] = +V1[1] + S3[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S3[2] * (-cI * (P1[0] * TMP17) + cI * (V1[2] * TMP19)); 
  V2[3] = denom * S3[2] * (-cI * (P1[1] * TMP17) + cI * (V1[3] * TMP19)); 
  V2[4] = denom * S3[2] * (-cI * (P1[2] * TMP17) + cI * (V1[4] * TMP19)); 
  V2[5] = denom * S3[2] * (-cI * (P1[3] * TMP17) + cI * (V1[5] * TMP19)); 
}


void FFV1C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP6; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP6 = (-1.) * (F1[2] * (F2[4] * (P3[0] - P3[3]) + F2[5] * (+cI * (P3[2]) -
      P3[1])) + (F1[3] * (F2[4] * (-1.) * (P3[1] + cI * (P3[2])) + F2[5] *
      (P3[0] + P3[3])) + (F1[4] * (F2[2] * (P3[0] + P3[3]) + F2[3] * (P3[1] -
      cI * (P3[2]))) + F1[5] * (F2[2] * (P3[1] + cI * (P3[2])) + F2[3] * (P3[0]
      - P3[3])))));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (-F2[4] * F1[2] - F2[5] * F1[3] - F2[2] * F1[4] -
      F2[3] * F1[5] - P3[0] * OM3 * TMP6);
  V3[3] = denom * (-cI) * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] -
      F2[4] * F1[3] - P3[1] * OM3 * TMP6);
  V3[4] = denom * (-cI) * (-cI * (F2[4] * F1[3] + F2[3] * F1[4]) + cI * (F2[5]
      * F1[2] + F2[2] * F1[5]) - P3[2] * OM3 * TMP6);
  V3[5] = denom * (-cI) * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] -
      F2[3] * F1[5] - P3[3] * OM3 * TMP6);
}


void FFV1_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * (-1.) * (V3[2] +
      V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * (-cI) * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void VVVV3P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (-cI * (V1[2] * TMP27) + cI * (V3[2] * TMP21)); 
  V4[3] = denom * (-cI * (V1[3] * TMP27) + cI * (V3[3] * TMP21)); 
  V4[4] = denom * (-cI * (V1[4] * TMP27) + cI * (V3[4] * TMP21)); 
  V4[5] = denom * (-cI * (V1[5] * TMP27) + cI * (V3[5] * TMP21)); 
}


void SSS1_0(std::complex<double> S1[], std::complex<double> S2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  vertex = COUP * - cI * S3[2] * S2[2] * S1[2]; 
}


void VVVS1_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> S4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  double P4[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  S4[0] = +V1[0] + V2[0] + V3[0]; 
  S4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -S4[0].real(); 
  P4[1] = -S4[1].real(); 
  P4[2] = -S4[1].imag(); 
  P4[3] = -S4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  S4[2] = denom * (TMP21 * (-cI * (TMP24) + cI * (TMP23)) + (TMP25 * (-cI *
      (TMP22) + cI * (TMP26)) + TMP27 * (-cI * (TMP18) + cI * (TMP17))));
}


void VVVVS3_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    S5[], std::complex<double> COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  vertex = COUP * S5[2] * (-cI * (TMP25 * TMP30) + cI * (TMP21 * TMP35)); 
}


void VVVV1P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (-cI * (V3[2] * TMP29) + cI * (TMP25 * V4[2])); 
  V2[3] = denom * (-cI * (V3[3] * TMP29) + cI * (TMP25 * V4[3])); 
  V2[4] = denom * (-cI * (V3[4] * TMP29) + cI * (TMP25 * V4[4])); 
  V2[5] = denom * (-cI * (V3[5] * TMP29) + cI * (TMP25 * V4[5])); 
}


void VVSS1_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + S3[0] + S4[0]; 
  V2[1] = +V1[1] + S3[1] + S4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S3[2] * S4[2] * (-cI * (V1[2]) + cI * (P2[0] * TMP17 * OM2)); 
  V2[3] = denom * S3[2] * S4[2] * (-cI * (V1[3]) + cI * (P2[1] * TMP17 * OM2)); 
  V2[4] = denom * S3[2] * S4[2] * (-cI * (V1[4]) + cI * (P2[2] * TMP17 * OM2)); 
  V2[5] = denom * S3[2] * S4[2] * (-cI * (V1[5]) + cI * (P2[3] * TMP17 * OM2)); 
}


void VSS1_1(std::complex<double> S2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP20; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  V1[0] = +S2[0] + S3[0]; 
  V1[1] = +S2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP20 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S2[2] * S3[2] * (OM1 * P1[0] * (-cI * (TMP20) + cI * (TMP19))
      + (-cI * (P2[0]) + cI * (P3[0])));
  V1[3] = denom * S2[2] * S3[2] * (OM1 * P1[1] * (-cI * (TMP20) + cI * (TMP19))
      + (-cI * (P2[1]) + cI * (P3[1])));
  V1[4] = denom * S2[2] * S3[2] * (OM1 * P1[2] * (-cI * (TMP20) + cI * (TMP19))
      + (-cI * (P2[2]) + cI * (P3[2])));
  V1[5] = denom * S2[2] * S3[2] * (OM1 * P1[3] * (-cI * (TMP20) + cI * (TMP19))
      + (-cI * (P2[3]) + cI * (P3[3])));
}


void VVVS1P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP26; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V3[0] = +V1[0] + V2[0] + S4[0]; 
  V3[1] = +V1[1] + V2[1] + S4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S4[2] * (TMP21 * (-cI * (P1[0]) + cI * (P2[0])) + (V1[2] *
      (-cI * (TMP26) + cI * (TMP22)) + V2[2] * (-cI * (TMP17) + cI * (TMP18))));
  V3[3] = denom * S4[2] * (TMP21 * (-cI * (P1[1]) + cI * (P2[1])) + (V1[3] *
      (-cI * (TMP26) + cI * (TMP22)) + V2[3] * (-cI * (TMP17) + cI * (TMP18))));
  V3[4] = denom * S4[2] * (TMP21 * (-cI * (P1[2]) + cI * (P2[2])) + (V1[4] *
      (-cI * (TMP26) + cI * (TMP22)) + V2[4] * (-cI * (TMP17) + cI * (TMP18))));
  V3[5] = denom * S4[2] * (TMP21 * (-cI * (P1[3]) + cI * (P2[3])) + (V1[5] *
      (-cI * (TMP26) + cI * (TMP22)) + V2[5] * (-cI * (TMP17) + cI * (TMP18))));
}


void FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP8; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP8 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP8); 
  V3[3] = denom * (-cI) * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP8);
  V3[4] = denom * (-cI) * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2]
      * OM3 * TMP8);
  V3[5] = denom * (-cI) * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP8); 
}

void FFV2_5_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV5_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_3_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV3_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2_3(F1, F2, COUP1, M3, W3, V3); 
  FFV4_3(F1, F2, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVVV1_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP23; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP31; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-cI * (TMP23 * TMP30) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[2]) + cI * (V3[2] * TMP30)));
  V1[3] = denom * (OM1 * P1[1] * (-cI * (TMP23 * TMP30) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[3]) + cI * (V3[3] * TMP30)));
  V1[4] = denom * (OM1 * P1[2] * (-cI * (TMP23 * TMP30) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[4]) + cI * (V3[4] * TMP30)));
  V1[5] = denom * (OM1 * P1[3] * (-cI * (TMP23 * TMP30) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[5]) + cI * (V3[5] * TMP30)));
}


void VVVVS2_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (OM3 * P3[0] * (-cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[2] * TMP29) + cI * (TMP21 * V4[2])));
  V3[3] = denom * S5[2] * (OM3 * P3[1] * (-cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[3] * TMP29) + cI * (TMP21 * V4[3])));
  V3[4] = denom * S5[2] * (OM3 * P3[2] * (-cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[4] * TMP29) + cI * (TMP21 * V4[4])));
  V3[5] = denom * S5[2] * (OM3 * P3[3] * (-cI * (TMP21 * TMP36) + cI * (TMP26 *
      TMP29)) + (-cI * (V2[5] * TMP29) + cI * (TMP21 * V4[5])));
}


void VVVVS1P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (-cI * (V2[2] * TMP29) + cI * (V1[2] * TMP30)); 
  V3[3] = denom * S5[2] * (-cI * (V2[3] * TMP29) + cI * (V1[3] * TMP30)); 
  V3[4] = denom * S5[2] * (-cI * (V2[4] * TMP29) + cI * (V1[4] * TMP30)); 
  V3[5] = denom * S5[2] * (-cI * (V2[5] * TMP29) + cI * (V1[5] * TMP30)); 
}


void FFV2C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP9; 
  TMP9 = (-1.) * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + F1[3] * (F2[4] * (-1.) * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])));
  vertex = COUP * - cI * TMP9; 
}

void FFV2_5C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2C1_0(F2, F1, V3, COUP1, vertex); 
  FFV5C1_0(F2, F1, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV2_3C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2C1_0(F2, F1, V3, COUP1, vertex); 
  FFV3C1_0(F2, F1, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}
void FFV2_4C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP1, std::complex<double>
    COUP2, std::complex<double> & vertex)
{
  std::complex<double> tmp; 
  FFV2C1_0(F2, F1, V3, COUP1, vertex); 
  FFV4C1_0(F2, F1, V3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

void FFS5_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (+cI * (TMP1 + TMP2)); 
}


void VVVV1P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (-cI * (V1[2] * TMP27) + cI * (V2[2] * TMP25)); 
  V4[3] = denom * (-cI * (V1[3] * TMP27) + cI * (V2[3] * TMP25)); 
  V4[4] = denom * (-cI * (V1[4] * TMP27) + cI * (V2[4] * TMP25)); 
  V4[5] = denom * (-cI * (V1[5] * TMP27) + cI * (V2[5] * TMP25)); 
}


void VVVV4_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-cI * (TMP22 * TMP35) + cI * (TMP23 * TMP30))
      + (-cI * (V3[2] * TMP30) + cI * (V2[2] * TMP35)));
  V1[3] = denom * (OM1 * P1[1] * (-cI * (TMP22 * TMP35) + cI * (TMP23 * TMP30))
      + (-cI * (V3[3] * TMP30) + cI * (V2[3] * TMP35)));
  V1[4] = denom * (OM1 * P1[2] * (-cI * (TMP22 * TMP35) + cI * (TMP23 * TMP30))
      + (-cI * (V3[4] * TMP30) + cI * (V2[4] * TMP35)));
  V1[5] = denom * (OM1 * P1[3] * (-cI * (TMP22 * TMP35) + cI * (TMP23 * TMP30))
      + (-cI * (V3[5] * TMP30) + cI * (V2[5] * TMP35)));
}


void FFV4_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP11; 
  std::complex<double> TMP12; 
  TMP11 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP12 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  vertex = COUP * (-1.) * (+cI * (TMP11) + 2. * cI * (TMP12)); 
}


void VVS2_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP22; 
  std::complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  V1[0] = +V2[0] + S3[0]; 
  V1[1] = +V2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * (-cI * (P2[0] * TMP22) + cI * (TMP19 * V2[2])); 
  V1[3] = denom * S3[2] * (-cI * (P2[1] * TMP22) + cI * (TMP19 * V2[3])); 
  V1[4] = denom * S3[2] * (-cI * (P2[2] * TMP22) + cI * (TMP19 * V2[4])); 
  V1[5] = denom * S3[2] * (-cI * (P2[3] * TMP22) + cI * (TMP19 * V2[5])); 
}


void FFV2P0_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])));
  F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
      - V3[2]));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (-1.) * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] *
      (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3]
      * (V3[2] - V3[5])))));
}


void VVV1P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V2[0] = +V1[0] + V3[0]; 
  V2[1] = +V1[1] + V3[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (TMP25 * (-cI * (P3[0]) + cI * (P1[0])) + (V1[2] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[2] * (-cI * (TMP17) + cI * (TMP18))));
  V2[3] = denom * (TMP25 * (-cI * (P3[1]) + cI * (P1[1])) + (V1[3] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[3] * (-cI * (TMP17) + cI * (TMP18))));
  V2[4] = denom * (TMP25 * (-cI * (P3[2]) + cI * (P1[2])) + (V1[4] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[4] * (-cI * (TMP17) + cI * (TMP18))));
  V2[5] = denom * (TMP25 * (-cI * (P3[3]) + cI * (P1[3])) + (V1[5] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[5] * (-cI * (TMP17) + cI * (TMP18))));
}


void VVVVS1_5(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, double M5, double W5, std::complex<double> S5[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P5[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P5[0] * P5[0]) - (P5[1] * P5[1]) - (P5[2] * P5[2]) - (P5[3] *
      P5[3]) - M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP25 * TMP30) + cI * (TMP27 * TMP29)); 
}


void VVVS1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP20; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  std::complex<double> TMP26; 
  std::complex<double> TMP28; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + S4[0]; 
  V3[1] = +V1[1] + V2[1] + S4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP28 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP20 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S4[2] * (OM3 * P3[0] * (TMP21 * (-cI * (TMP28) + cI *
      (TMP20)) + (-cI * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 *
      (-cI * (P1[0]) + cI * (P2[0])) + (V1[2] * (-cI * (TMP26) + cI * (TMP22))
      + V2[2] * (-cI * (TMP17) + cI * (TMP18)))));
  V3[3] = denom * S4[2] * (OM3 * P3[1] * (TMP21 * (-cI * (TMP28) + cI *
      (TMP20)) + (-cI * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 *
      (-cI * (P1[1]) + cI * (P2[1])) + (V1[3] * (-cI * (TMP26) + cI * (TMP22))
      + V2[3] * (-cI * (TMP17) + cI * (TMP18)))));
  V3[4] = denom * S4[2] * (OM3 * P3[2] * (TMP21 * (-cI * (TMP28) + cI *
      (TMP20)) + (-cI * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 *
      (-cI * (P1[2]) + cI * (P2[2])) + (V1[4] * (-cI * (TMP26) + cI * (TMP22))
      + V2[4] * (-cI * (TMP17) + cI * (TMP18)))));
  V3[5] = denom * S4[2] * (OM3 * P3[3] * (TMP21 * (-cI * (TMP28) + cI *
      (TMP20)) + (-cI * (TMP18 * TMP22) + cI * (TMP17 * TMP26))) + (TMP21 *
      (-cI * (P1[3]) + cI * (P2[3])) + (V1[5] * (-cI * (TMP26) + cI * (TMP22))
      + V2[5] * (-cI * (TMP17) + cI * (TMP18)))));
}


void VVVVS1P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (-cI * (V3[2] * TMP29) + cI * (TMP25 * V4[2])); 
  V2[3] = denom * S5[2] * (-cI * (V3[3] * TMP29) + cI * (TMP25 * V4[3])); 
  V2[4] = denom * S5[2] * (-cI * (V3[4] * TMP29) + cI * (TMP25 * V4[4])); 
  V2[5] = denom * S5[2] * (-cI * (V3[5] * TMP29) + cI * (TMP25 * V4[5])); 
}


void VVVV2_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> TMP32; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (-2. * cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29 + TMP25 * TMP32)) + (-cI * (V3[2] * TMP29 + TMP25 * V4[2]) + 2. *
      cI * (V1[2] * TMP35)));
  V2[3] = denom * (OM2 * P2[1] * (-2. * cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29 + TMP25 * TMP32)) + (-cI * (V3[3] * TMP29 + TMP25 * V4[3]) + 2. *
      cI * (V1[3] * TMP35)));
  V2[4] = denom * (OM2 * P2[2] * (-2. * cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29 + TMP25 * TMP32)) + (-cI * (V3[4] * TMP29 + TMP25 * V4[4]) + 2. *
      cI * (V1[4] * TMP35)));
  V2[5] = denom * (OM2 * P2[3] * (-2. * cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29 + TMP25 * TMP32)) + (-cI * (V3[5] * TMP29 + TMP25 * V4[5]) + 2. *
      cI * (V1[5] * TMP35)));
}


void FFV5C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * (-cI) * (F1[2] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] -
      cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[3] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * (-1.) *
      (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * 4. * (V3[2] + V3[5]) + 4. *
      (F1[5] * (V3[3] + cI * (V3[4]))))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + (F1[3] * (P2[0] * (-1.) * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[4] * 4. * (V3[3] - cI * (V3[4])) + 4. *
      (F1[5] * (V3[2] - V3[5])))));
  F2[4] = denom * (-4. * cI) * (F1[4] * (P2[0] * (-1.) * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + (+1./4. * (M2 * (+4. * (F1[2] * 1./4. * (V3[2] -
      V3[5])) - F1[3] * (V3[3] + cI * (V3[4])))) + F1[5] * (P2[0] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[1] * (V3[2] - V3[5]) + (P2[2] * (-cI *
      (V3[5]) + cI * (V3[2])) + P2[3] * (V3[3] + cI * (V3[4])))))));
  F2[5] = denom * (-4. * cI) * (F1[4] * (P2[0] * (+cI * (V3[4]) - V3[3]) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (+cI * (V3[4]) - V3[3])))) + (+1./4. * (M2 * (F1[3] * (V3[2] +
      V3[5]) + 4. * (F1[2] * 1./4. * (+cI * (V3[4]) - V3[3])))) + F1[5] *
      (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] *
      (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] - V3[2]))))));
}


void VVSS1_1(std::complex<double> V2[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + S3[0] + S4[0]; 
  V1[1] = +V2[1] + S3[1] + S4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S3[2] * S4[2] * (-cI * (V2[2]) + cI * (P1[0] * OM1 * TMP22)); 
  V1[3] = denom * S3[2] * S4[2] * (-cI * (V2[3]) + cI * (P1[1] * OM1 * TMP22)); 
  V1[4] = denom * S3[2] * S4[2] * (-cI * (V2[4]) + cI * (P1[2] * OM1 * TMP22)); 
  V1[5] = denom * S3[2] * S4[2] * (-cI * (V2[5]) + cI * (P1[3] * OM1 * TMP22)); 
}


void FFS3_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * S3[2] * (-cI * (TMP1) + cI * (TMP2)); 
}


void VSS1_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> S2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> denom; 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  S2[0] = +V1[0] + S3[0]; 
  S2[1] = +V1[1] + S3[1]; 
  P2[0] = -S2[0].real(); 
  P2[1] = -S2[1].real(); 
  P2[2] = -S2[1].imag(); 
  P2[3] = -S2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  S2[2] = denom * S3[2] * (-cI * (TMP18) + cI * (TMP17)); 
}


void FFS3C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP2) + cI * (TMP1)); 
}


void FFV5_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3]
      - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] -
      V3[2])))) + (+1./4. * (M1 * (F2[5] * (V3[3] + cI * (V3[4])) + 4. * (F2[4]
      * 1./4. * (V3[2] + V3[5])))) + F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) +
      (P1[1] * (-1.) * (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] +
      V3[5])) + P1[3] * (V3[3] + cI * (V3[4])))))));
  F1[3] = denom * 4. * cI * (F2[2] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (+1./4. * (M1 * (F2[5] * (V3[2] - V3[5]) + 4. *
      (F2[4] * 1./4. * (V3[3] - cI * (V3[4]))))) + F2[3] * (P1[0] * (-1.) *
      (V3[2] + V3[5]) + (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI
      * (V3[3])) + P1[3] * (V3[2] + V3[5]))))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * 4. * (V3[5] - V3[2]) + 4. *
      (F2[3] * (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * 4. * (+cI * (V3[4]) - V3[3]) + 4. * (F2[3] *
      (V3[2] + V3[5])))));
}


void VVVV1_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> TMP32; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (-cI * (TMP25 * TMP32) + cI * (TMP24 * TMP29))
      + (-cI * (V3[2] * TMP29) + cI * (TMP25 * V4[2])));
  V2[3] = denom * (OM2 * P2[1] * (-cI * (TMP25 * TMP32) + cI * (TMP24 * TMP29))
      + (-cI * (V3[3] * TMP29) + cI * (TMP25 * V4[3])));
  V2[4] = denom * (OM2 * P2[2] * (-cI * (TMP25 * TMP32) + cI * (TMP24 * TMP29))
      + (-cI * (V3[4] * TMP29) + cI * (TMP25 * V4[4])));
  V2[5] = denom * (OM2 * P2[3] * (-cI * (TMP25 * TMP32) + cI * (TMP24 * TMP29))
      + (-cI * (V3[5] * TMP29) + cI * (TMP25 * V4[5])));
}


void FFS1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP0; 
  TMP0 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * - cI * TMP0 * S3[2]; 
}


void VVVVS1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (-cI * (TMP27 * V4[2]) + cI * (V3[2] * TMP30)); 
  V1[3] = denom * S5[2] * (-cI * (TMP27 * V4[3]) + cI * (V3[3] * TMP30)); 
  V1[4] = denom * S5[2] * (-cI * (TMP27 * V4[4]) + cI * (V3[4] * TMP30)); 
  V1[5] = denom * S5[2] * (-cI * (TMP27 * V4[5]) + cI * (V3[5] * TMP30)); 
}


void FFV2C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP10; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP10 = (-1.) * (F1[2] * (F2[4] * (P3[0] - P3[3]) + F2[5] * (+cI * (P3[2]) -
      P3[1])) + F1[3] * (F2[4] * (-1.) * (P3[1] + cI * (P3[2])) + F2[5] *
      (P3[0] + P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (-F2[4] * F1[2] - F2[5] * F1[3] - P3[0] * OM3 *
      TMP10);
  V3[3] = denom * (-cI) * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 *
      TMP10);
  V3[4] = denom * (-cI) * (-cI * (F2[4] * F1[3]) + cI * (F2[5] * F1[2]) - P3[2]
      * OM3 * TMP10);
  V3[5] = denom * (-cI) * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 *
      TMP10);
}

void FFV2_5C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2C1_3(F2, F1, COUP1, M3, W3, V3); 
  FFV5C1_3(F2, F1, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_3C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2C1_3(F2, F1, COUP1, M3, W3, V3); 
  FFV3C1_3(F2, F1, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}
void FFV2_4C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> V3[])
{
  std::complex<double> Vtmp[6]; 
  int i; 
  FFV2C1_3(F2, F1, COUP1, M3, W3, V3); 
  FFV4C1_3(F2, F1, COUP2, M3, W3, Vtmp); 
  i = 2; 
  while (i < 6)
  {
    V3[i] = V3[i] + Vtmp[i]; 
    i++; 
  }
}

void VVVV5_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> TMP32; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * 1./2. * (OM2 * - P2[0] * (-2. * cI * (TMP24 * TMP29) + cI *
      (TMP25 * TMP32 + TMP17 * TMP35)) + (-2. * cI * (V3[2] * TMP29) + cI *
      (TMP25 * V4[2] + V1[2] * TMP35)));
  V2[3] = denom * 1./2. * (OM2 * - P2[1] * (-2. * cI * (TMP24 * TMP29) + cI *
      (TMP25 * TMP32 + TMP17 * TMP35)) + (-2. * cI * (V3[3] * TMP29) + cI *
      (TMP25 * V4[3] + V1[3] * TMP35)));
  V2[4] = denom * 1./2. * (OM2 * - P2[2] * (-2. * cI * (TMP24 * TMP29) + cI *
      (TMP25 * TMP32 + TMP17 * TMP35)) + (-2. * cI * (V3[4] * TMP29) + cI *
      (TMP25 * V4[4] + V1[4] * TMP35)));
  V2[5] = denom * 1./2. * (OM2 * - P2[3] * (-2. * cI * (TMP24 * TMP29) + cI *
      (TMP25 * TMP32 + TMP17 * TMP35)) + (-2. * cI * (V3[5] * TMP29) + cI *
      (TMP25 * V4[5] + V1[5] * TMP35)));
}


void VVVV4P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI * (V1[2] * TMP30) + cI * (TMP21 * V4[2])); 
  V3[3] = denom * (-cI * (V1[3] * TMP30) + cI * (TMP21 * V4[3])); 
  V3[4] = denom * (-cI * (V1[4] * TMP30) + cI * (TMP21 * V4[4])); 
  V3[5] = denom * (-cI * (V1[5] * TMP30) + cI * (TMP21 * V4[5])); 
}


void VVVV4_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP34; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * (OM4 * P4[0] * (-cI * (TMP21 * TMP37) + cI * (TMP25 * TMP34))
      + (-cI * (V2[2] * TMP25) + cI * (V3[2] * TMP21)));
  V4[3] = denom * (OM4 * P4[1] * (-cI * (TMP21 * TMP37) + cI * (TMP25 * TMP34))
      + (-cI * (V2[3] * TMP25) + cI * (V3[3] * TMP21)));
  V4[4] = denom * (OM4 * P4[2] * (-cI * (TMP21 * TMP37) + cI * (TMP25 * TMP34))
      + (-cI * (V2[4] * TMP25) + cI * (V3[4] * TMP21)));
  V4[5] = denom * (OM4 * P4[3] * (-cI * (TMP21 * TMP37) + cI * (TMP25 * TMP34))
      + (-cI * (V2[5] * TMP25) + cI * (V3[5] * TMP21)));
}


void FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP13; 
  std::complex<double> TMP8; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP13 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP8 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[0] * (TMP8 + 2. * (TMP13)) +
      (+1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[1] * (TMP8 + 2. * (TMP13)) +
      (-1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP8 + 2. * (TMP13)) +
      (+1./2. * cI * (F2[5] * F1[2]) - 1./2. * cI * (F2[4] * F1[3]) - cI *
      (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
  V3[5] = denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP8 + 2. * (TMP13)) +
      (+1./2. * (F2[4] * F1[2]) - 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void VVVVS2P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (-cI * (V1[2] * TMP27) + cI * (V3[2] * TMP21)); 
  V4[3] = denom * S5[2] * (-cI * (V1[3] * TMP27) + cI * (V3[3] * TMP21)); 
  V4[4] = denom * S5[2] * (-cI * (V1[4] * TMP27) + cI * (V3[4] * TMP21)); 
  V4[5] = denom * S5[2] * (-cI * (V1[5] * TMP27) + cI * (V3[5] * TMP21)); 
}


void VVSS1P0_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> S4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + S3[0] + S4[0]; 
  V2[1] = +V1[1] + S3[1] + S4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * - cI * V1[2] * S4[2] * S3[2]; 
  V2[3] = denom * - cI * V1[3] * S4[2] * S3[2]; 
  V2[4] = denom * - cI * V1[4] * S4[2] * S3[2]; 
  V2[5] = denom * - cI * V1[5] * S4[2] * S3[2]; 
}


void FFS1C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP0; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP0 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP0; 
}


void VVV1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP26; 
  std::complex<double> TMP27; 
  std::complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP27 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[2] * (-cI * (TMP26) + cI * (TMP22))));
  V1[3] = denom * (TMP27 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[3] * (-cI * (TMP26) + cI * (TMP22))));
  V1[4] = denom * (TMP27 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[4] * (-cI * (TMP26) + cI * (TMP22))));
  V1[5] = denom * (TMP27 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP23) + cI * (TMP24)) + V3[5] * (-cI * (TMP26) + cI * (TMP22))));
}


void VVVVS1_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> TMP32; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (OM2 * P2[0] * (-cI * (TMP25 * TMP32) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[2] * TMP29) + cI * (TMP25 * V4[2])));
  V2[3] = denom * S5[2] * (OM2 * P2[1] * (-cI * (TMP25 * TMP32) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[3] * TMP29) + cI * (TMP25 * V4[3])));
  V2[4] = denom * S5[2] * (OM2 * P2[2] * (-cI * (TMP25 * TMP32) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[4] * TMP29) + cI * (TMP25 * V4[4])));
  V2[5] = denom * S5[2] * (OM2 * P2[3] * (-cI * (TMP25 * TMP32) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[5] * TMP29) + cI * (TMP25 * V4[5])));
}


void VVVVS3P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (-cI * (TMP25 * V4[2]) + cI * (V1[2] * TMP35)); 
  V2[3] = denom * S5[2] * (-cI * (TMP25 * V4[3]) + cI * (V1[3] * TMP35)); 
  V2[4] = denom * S5[2] * (-cI * (TMP25 * V4[4]) + cI * (V1[4] * TMP35)); 
  V2[5] = denom * S5[2] * (-cI * (TMP25 * V4[5]) + cI * (V1[5] * TMP35)); 
}


void FFV3_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 2. * (V3[5] - V3[2]) + 2. * (F1[5]
      * (V3[3] - cI * (V3[4]))))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (-1.) * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) + P2[3] * (V3[2] - V3[5])))) + M2 * (F1[4] * 2. * (V3[3] + cI *
      (V3[4])) - 2. * (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3]
      + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (+1./2. * (M2 * (F1[3] * (V3[3] - cI * (V3[4])) + 2. * (F1[2]
      * 1./2. * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) +
      (P2[1] * (-1.) * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * 2. * cI * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M2 * (F1[3] * (V3[2] - V3[5]) + 2.
      * (F1[2] * 1./2. * (V3[3] + cI * (V3[4]))))) + F1[5] * (P2[0] * (-1.) *
      (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI
      * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


void FFV4C1_3(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP10; 
  std::complex<double> TMP16; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP10 = (-1.) * (F1[2] * (F2[4] * (P3[0] - P3[3]) + F2[5] * (+cI * (P3[2]) -
      P3[1])) + F1[3] * (F2[4] * (-1.) * (P3[1] + cI * (P3[2])) + F2[5] *
      (P3[0] + P3[3])));
  TMP16 = (-1.) * (F1[4] * (F2[2] * (P3[0] + P3[3]) + F2[3] * (P3[1] - cI *
      (P3[2]))) + F1[5] * (F2[2] * (P3[1] + cI * (P3[2])) + F2[3] * (P3[0] -
      P3[3])));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * 2. * cI * (OM3 * 1./2. * P3[0] * (TMP10 + 2. * (TMP16)) +
      (+1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] *
      F1[5]));
  V3[3] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[1] * (TMP10 + 2. * (TMP16))
      + (-1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] *
      F1[5]));
  V3[4] = denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP10 + 2. * (TMP16)) +
      (-1./2. * cI * (F2[5] * F1[2]) + 1./2. * cI * (F2[4] * F1[3]) - cI *
      (F2[2] * F1[5]) + cI * (F2[3] * F1[4])));
  V3[5] = denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP10 + 2. * (TMP16)) +
      (+1./2. * (F2[4] * F1[2]) - 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] +
      F2[3] * F1[5]));
}


void VVVVS2_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (OM4 * P4[0] * (-cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[2] * TMP27) + cI * (V3[2] * TMP21)));
  V4[3] = denom * S5[2] * (OM4 * P4[1] * (-cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[3] * TMP27) + cI * (V3[3] * TMP21)));
  V4[4] = denom * S5[2] * (OM4 * P4[2] * (-cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[4] * TMP27) + cI * (V3[4] * TMP21)));
  V4[5] = denom * S5[2] * (OM4 * P4[3] * (-cI * (TMP21 * TMP37) + cI * (TMP27 *
      TMP33)) + (-cI * (V1[5] * TMP27) + cI * (V3[5] * TMP21)));
}


void FFS3_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * (-cI * (TMP2) + cI * (TMP1)); 
}


void VVS2P0_2(std::complex<double> V1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP19; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  V2[0] = +V1[0] + S3[0]; 
  V2[1] = +V1[1] + S3[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S3[2] * (-cI * (P1[0] * TMP17) + cI * (V1[2] * TMP19)); 
  V2[3] = denom * S3[2] * (-cI * (P1[1] * TMP17) + cI * (V1[3] * TMP19)); 
  V2[4] = denom * S3[2] * (-cI * (P1[2] * TMP17) + cI * (V1[4] * TMP19)); 
  V2[5] = denom * S3[2] * (-cI * (P1[3] * TMP17) + cI * (V1[5] * TMP19)); 
}


void FFS4_2(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + F1[5] * (+cI *
      (P2[2]) - P2[1]));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) - F1[5] *
      (P2[0] + P2[3]));
  F2[4] = denom * cI * F1[4] * M2 * S3[2]; 
  F2[5] = denom * cI * F1[5] * M2 * S3[2]; 
}


void SSSS1_0(std::complex<double> S1[], std::complex<double> S2[],
    std::complex<double> S3[], std::complex<double> S4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  vertex = COUP * - cI * S4[2] * S3[2] * S2[2] * S1[2]; 
}


void FFS1_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP0; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP0 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP0; 
}


void VVVVS3_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP30; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (OM3 * P3[0] * (-cI * (TMP21 * TMP36) + cI * (TMP18 *
      TMP30)) + (-cI * (V1[2] * TMP30) + cI * (TMP21 * V4[2])));
  V3[3] = denom * S5[2] * (OM3 * P3[1] * (-cI * (TMP21 * TMP36) + cI * (TMP18 *
      TMP30)) + (-cI * (V1[3] * TMP30) + cI * (TMP21 * V4[3])));
  V3[4] = denom * S5[2] * (OM3 * P3[2] * (-cI * (TMP21 * TMP36) + cI * (TMP18 *
      TMP30)) + (-cI * (V1[4] * TMP30) + cI * (TMP21 * V4[4])));
  V3[5] = denom * S5[2] * (OM3 * P3[3] * (-cI * (TMP21 * TMP36) + cI * (TMP18 *
      TMP30)) + (-cI * (V1[5] * TMP30) + cI * (TMP21 * V4[5])));
}


void FFS5C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * - S3[2] * (+cI * (TMP1 + TMP2)); 
}


void VVVV5_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP31; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * 1./2. * (OM1 * - P1[0] * (-2. * cI * (TMP27 * TMP31) + cI *
      (TMP23 * TMP30 + TMP22 * TMP35)) + (-2. * cI * (TMP27 * V4[2]) + cI *
      (V3[2] * TMP30 + V2[2] * TMP35)));
  V1[3] = denom * 1./2. * (OM1 * - P1[1] * (-2. * cI * (TMP27 * TMP31) + cI *
      (TMP23 * TMP30 + TMP22 * TMP35)) + (-2. * cI * (TMP27 * V4[3]) + cI *
      (V3[3] * TMP30 + V2[3] * TMP35)));
  V1[4] = denom * 1./2. * (OM1 * - P1[2] * (-2. * cI * (TMP27 * TMP31) + cI *
      (TMP23 * TMP30 + TMP22 * TMP35)) + (-2. * cI * (TMP27 * V4[4]) + cI *
      (V3[4] * TMP30 + V2[4] * TMP35)));
  V1[5] = denom * 1./2. * (OM1 * - P1[3] * (-2. * cI * (TMP27 * TMP31) + cI *
      (TMP23 * TMP30 + TMP22 * TMP35)) + (-2. * cI * (TMP27 * V4[5]) + cI *
      (V3[5] * TMP30 + V2[5] * TMP35)));
}


void VVVVS2P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP29; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0] + S5[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1] + S5[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * S5[2] * (-cI * (V2[2] * TMP29) + cI * (TMP21 * V4[2])); 
  V3[3] = denom * S5[2] * (-cI * (V2[3] * TMP29) + cI * (TMP21 * V4[3])); 
  V3[4] = denom * S5[2] * (-cI * (V2[4] * TMP29) + cI * (TMP21 * V4[4])); 
  V3[5] = denom * S5[2] * (-cI * (V2[5] * TMP29) + cI * (TMP21 * V4[5])); 
}


void VVVV3_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  vertex = COUP * (-cI * (TMP27 * TMP29) + cI * (TMP21 * TMP35)); 
}


void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] +
      F2[3] * F1[5]);
  V3[3] = denom * (-cI) * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] -
      F2[4] * F1[3]);
  V3[4] = denom * (-cI) * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4]
      * F1[3] + F2[3] * F1[4]));
  V3[5] = denom * (-cI) * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] -
      F2[3] * F1[5]);
}


void FFS1C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP0; 
  TMP0 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * - cI * TMP0 * S3[2]; 
}


void VVVVS3P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M1, double W1, std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0] + S5[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1] + S5[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * S5[2] * (-cI * (V3[2] * TMP30) + cI * (V2[2] * TMP35)); 
  V1[3] = denom * S5[2] * (-cI * (V3[3] * TMP30) + cI * (V2[3] * TMP35)); 
  V1[4] = denom * S5[2] * (-cI * (V3[4] * TMP30) + cI * (V2[4] * TMP35)); 
  V1[5] = denom * S5[2] * (-cI * (V3[5] * TMP30) + cI * (V2[5] * TMP35)); 
}


void VVVV2P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-1.) * (-2. * cI * (V2[2] * TMP35) + cI * (TMP27 * V4[2] +
      V3[2] * TMP30));
  V1[3] = denom * (-1.) * (-2. * cI * (V2[3] * TMP35) + cI * (TMP27 * V4[3] +
      V3[3] * TMP30));
  V1[4] = denom * (-1.) * (-2. * cI * (V2[4] * TMP35) + cI * (TMP27 * V4[4] +
      V3[4] * TMP30));
  V1[5] = denom * (-1.) * (-2. * cI * (V2[5] * TMP35) + cI * (TMP27 * V4[5] +
      V3[5] * TMP30));
}


void FFS2C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * F2[2] * M1 * S3[2]; 
  F1[3] = denom * cI * F2[3] * M1 * S3[2]; 
  F1[4] = denom * - cI * S3[2] * (F2[2] * (-1.) * (P1[0] + P1[3]) + F2[3] *
      (+cI * (P1[2]) - P1[1]));
  F1[5] = denom * cI * S3[2] * (F2[2] * (P1[1] + cI * (P1[2])) + F2[3] * (P1[0]
      - P1[3]));
}

void FFS2_4C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFS2C1_2(F2, S3, COUP1, M1, W1, F1); 
  FFS4C1_2(F2, S3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV1P0_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP27 * V4[2]) + cI * (V3[2] * TMP30)); 
  V1[3] = denom * (-cI * (TMP27 * V4[3]) + cI * (V3[3] * TMP30)); 
  V1[4] = denom * (-cI * (TMP27 * V4[4]) + cI * (V3[4] * TMP30)); 
  V1[5] = denom * (-cI * (TMP27 * V4[5]) + cI * (V3[5] * TMP30)); 
}


void FFV3C1_1(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] - V3[5]) + (P2[1] * (+cI *
      (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI * (V3[3])) + P2[3] *
      (V3[2] - V3[5])))) + (F1[3] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * 2. * (V3[2] + V3[5]) + 2. *
      (F1[5] * (V3[3] + cI * (V3[4]))))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + (F1[3] * (P2[0] * (-1.) * (V3[2] + V3[5]) +
      (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[4] * 2. * (+cI * (V3[4]) - V3[3]) + 2. *
      (F1[5] * (V3[5] - V3[2])))));
  F2[4] = denom * 2. * cI * (F1[4] * (P2[0] * (-1.) * (V3[2] + V3[5]) + (P2[1]
      * (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] *
      (V3[2] + V3[5])))) + (+1./2. * (M2 * (F1[3] * (V3[3] + cI * (V3[4])) + 2.
      * (F1[2] * 1./2. * (V3[5] - V3[2])))) + F1[5] * (P2[0] * (-1.) * (V3[3] +
      cI * (V3[4])) + (P2[1] * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI *
      (V3[2])) + P2[3] * (V3[3] + cI * (V3[4])))))));
  F2[5] = denom * 2. * cI * (F1[4] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1]
      * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P2[3] *
      (+cI * (V3[4]) - V3[3])))) + (+1./2. * (M2 * (+2. * (F1[2] * 1./2. *
      (V3[3] - cI * (V3[4]))) - F1[3] * (V3[2] + V3[5]))) + F1[5] * (P2[0] *
      (V3[5] - V3[2]) + (P2[1] * (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI
      * (V3[3])) + P2[3] * (V3[5] - V3[2]))))));
}


void SSSS1_3(std::complex<double> S1[], std::complex<double> S2[],
    std::complex<double> S4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> denom; 
  S3[0] = +S1[0] + S2[0] + S4[0]; 
  S3[1] = +S1[1] + S2[1] + S4[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * S4[2] * S2[2] * S1[2]; 
}


void VVVVS3_5(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, double M5, double W5, std::complex<double> S5[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P5[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  S5[0] = +V1[0] + V2[0] + V3[0] + V4[0]; 
  S5[1] = +V1[1] + V2[1] + V3[1] + V4[1]; 
  P5[0] = -S5[0].real(); 
  P5[1] = -S5[1].real(); 
  P5[2] = -S5[1].imag(); 
  P5[3] = -S5[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P5[0] * P5[0]) - (P5[1] * P5[1]) - (P5[2] * P5[2]) - (P5[3] *
      P5[3]) - M5 * (M5 - cI * W5));
  S5[2] = denom * (-cI * (TMP21 * TMP35) + cI * (TMP25 * TMP30)); 
}


void VVS1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> denom; 
  S3[0] = +V1[0] + V2[0]; 
  S3[1] = +V1[1] + V2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP21; 
}


void VVVVS1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    S5[], std::complex<double> COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  vertex = COUP * S5[2] * (-cI * (TMP27 * TMP29) + cI * (TMP25 * TMP30)); 
}


void FFS4C1_2(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * S3[2] * (F2[4] * (P1[0] - P1[3]) + F2[5] * (+cI *
      (P1[2]) - P1[1]));
  F1[3] = denom * - cI * S3[2] * (F2[4] * (P1[1] + cI * (P1[2])) - F2[5] *
      (P1[0] + P1[3]));
  F1[4] = denom * cI * F2[4] * M1 * S3[2]; 
  F1[5] = denom * cI * F2[5] * M1 * S3[2]; 
}


void FFS2_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP1; 
  std::complex<double> denom; 
  S3[0] = +F1[0] + F2[0]; 
  S3[1] = +F1[1] + F2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP1; 
}

void FFS2_4_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
    W3, std::complex<double> S3[])
{
  std::complex<double> Stmp[3]; 
  int i; 
  FFS2_3(F1, F2, COUP1, M3, W3, S3); 
  FFS4_3(F1, F2, COUP2, M3, W3, Stmp); 
  i = 2; 
  while (i < 3)
  {
    S3[i] = S3[i] + Stmp[i]; 
    i++; 
  }
}

void FFV1C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * (-cI) * (F2[2] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3] *
      (V3[2] + V3[5])))) + (F2[3] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + M1 * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI *
      (V3[4]) - V3[3]))));
  F1[3] = denom * cI * (F2[2] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[4] * (V3[3] + cI * (V3[4])) - F2[5] * (V3[2] +
      V3[5]))));
  F1[4] = denom * cI * (F2[4] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] + cI
      * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[5] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * (-1.) * (V3[2] +
      V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] - cI *
      (V3[4]))))) + M1 * (F2[2] * (-1.) * (V3[2] + V3[5]) + F2[3] * (+cI *
      (V3[4]) - V3[3]))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P1[3] * (V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] + V3[5]) +
      (P1[1] * (+cI * (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI *
      (V3[3])) - P1[3] * (V3[2] + V3[5])))) + M1 * (F2[2] * (V3[3] + cI *
      (V3[4])) + F2[3] * (V3[2] - V3[5]))));
}


void VVV1_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP19; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP28; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V2[0] = +V1[0] + V3[0]; 
  V2[1] = +V1[1] + V3[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP28 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (TMP25 * (-cI * (TMP19) + cI * (TMP28)) + (-cI
      * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 * (-cI * (P3[0]) + cI
      * (P1[0])) + (V1[2] * (-cI * (TMP23) + cI * (TMP24)) + V3[2] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V2[3] = denom * (OM2 * P2[1] * (TMP25 * (-cI * (TMP19) + cI * (TMP28)) + (-cI
      * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 * (-cI * (P3[1]) + cI
      * (P1[1])) + (V1[3] * (-cI * (TMP23) + cI * (TMP24)) + V3[3] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V2[4] = denom * (OM2 * P2[2] * (TMP25 * (-cI * (TMP19) + cI * (TMP28)) + (-cI
      * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 * (-cI * (P3[2]) + cI
      * (P1[2])) + (V1[4] * (-cI * (TMP23) + cI * (TMP24)) + V3[4] * (-cI *
      (TMP17) + cI * (TMP18)))));
  V2[5] = denom * (OM2 * P2[3] * (TMP25 * (-cI * (TMP19) + cI * (TMP28)) + (-cI
      * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 * (-cI * (P3[3]) + cI
      * (P1[3])) + (V1[5] * (-cI * (TMP23) + cI * (TMP24)) + V3[5] * (-cI *
      (TMP17) + cI * (TMP18)))));
}


void FFV1_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP3; 
  TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP3; 
}


void SSS1_1(std::complex<double> S2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> S1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  S1[0] = +S2[0] + S3[0]; 
  S1[1] = +S2[1] + S3[1]; 
  P1[0] = -S1[0].real(); 
  P1[1] = -S1[1].real(); 
  P1[2] = -S1[1].imag(); 
  P1[3] = -S1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  S1[2] = denom * cI * S3[2] * S2[2]; 
}


void VVVV1P0_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI * (V2[2] * TMP29) + cI * (V1[2] * TMP30)); 
  V3[3] = denom * (-cI * (V2[3] * TMP29) + cI * (V1[3] * TMP30)); 
  V3[4] = denom * (-cI * (V2[4] * TMP29) + cI * (V1[4] * TMP30)); 
  V3[5] = denom * (-cI * (V2[5] * TMP29) + cI * (V1[5] * TMP30)); 
}


void VVVVS3P0_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> denom; 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (-cI * (V2[2] * TMP25) + cI * (V3[2] * TMP21)); 
  V4[3] = denom * S5[2] * (-cI * (V2[3] * TMP25) + cI * (V3[3] * TMP21)); 
  V4[4] = denom * S5[2] * (-cI * (V2[4] * TMP25) + cI * (V3[4] * TMP21)); 
  V4[5] = denom * S5[2] * (-cI * (V2[5] * TMP25) + cI * (V3[5] * TMP21)); 
}


void VVVV2P0_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP25; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * (-1.) * (-2. * cI * (V1[2] * TMP35) + cI * (V3[2] * TMP29 +
      TMP25 * V4[2]));
  V2[3] = denom * (-1.) * (-2. * cI * (V1[3] * TMP35) + cI * (V3[3] * TMP29 +
      TMP25 * V4[3]));
  V2[4] = denom * (-1.) * (-2. * cI * (V1[4] * TMP35) + cI * (V3[4] * TMP29 +
      TMP25 * V4[4]));
  V2[5] = denom * (-1.) * (-2. * cI * (V1[5] * TMP35) + cI * (V3[5] * TMP29 +
      TMP25 * V4[5]));
}


void FFV2_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (-1.) * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + F1[3] * (P2[0] * (V3[2] - V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) + P2[3] * (V3[2] - V3[5])))));
  F2[4] = denom * - cI * M2 * (F1[2] * (-1.) * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]));
  F2[5] = denom * cI * M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] * (V3[2] -
      V3[5]));
}

void FFV2_5_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_2(F1, V3, COUP1, M2, W2, F2); 
  FFV5_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_3_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_2(F1, V3, COUP1, M2, W2, F2); 
  FFV3_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_4_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M2, double
    W2, std::complex<double> F2[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_2(F1, V3, COUP1, M2, W2, F2); 
  FFV4_2(F1, V3, COUP2, M2, W2, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F2[i] = F2[i] + Ftmp[i]; 
    i++; 
  }
}

void FFV3C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * (-2. * cI) * (F2[2] * (P1[0] * (-1.) * (V3[2] + V3[5]) +
      (P1[1] * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3]
      * (V3[2] + V3[5])))) + (+1./2. * (M1 * (F2[5] * (+cI * (V3[4]) - V3[3]) +
      2. * (F2[4] * 1./2. * (V3[2] - V3[5])))) + F2[3] * (P1[0] * (+cI *
      (V3[4]) - V3[3]) + (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) +
      cI * (V3[5])) + P1[3] * (V3[3] - cI * (V3[4])))))));
  F1[3] = denom * (-2. * cI) * (F2[2] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4]))
      + (P1[1] * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./2. * (M1 * (F2[5] * (V3[2] + V3[5]) + 2.
      * (F2[4] * (-1./2.) * (V3[3] + cI * (V3[4]))))) + F2[3] * (P1[0] * (V3[5]
      - V3[2]) + (P1[1] * (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI *
      (V3[3])) + P1[3] * (V3[5] - V3[2]))))));
  F1[4] = denom * cI * (F2[4] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] + cI
      * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[5] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] * (-1.) * (V3[2] +
      V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] - cI *
      (V3[4]))))) + M1 * (F2[2] * 2. * (V3[2] + V3[5]) + 2. * (F2[3] * (V3[3] -
      cI * (V3[4]))))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (-1.) * (V3[2] + V3[5]) +
      (P1[1] * (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[2] * 2. * (V3[3] + cI * (V3[4])) + 2. *
      (F2[3] * (V3[2] - V3[5])))));
}


void VVVVS2_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP24; 
  std::complex<double> TMP29; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (OM2 * P2[0] * (-cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[2] * TMP29) + cI * (V1[2] * TMP35)));
  V2[3] = denom * S5[2] * (OM2 * P2[1] * (-cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[3] * TMP29) + cI * (V1[3] * TMP35)));
  V2[4] = denom * S5[2] * (OM2 * P2[2] * (-cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[4] * TMP29) + cI * (V1[4] * TMP35)));
  V2[5] = denom * S5[2] * (OM2 * P2[3] * (-cI * (TMP17 * TMP35) + cI * (TMP24 *
      TMP29)) + (-cI * (V3[5] * TMP29) + cI * (V1[5] * TMP35)));
}


void FFS2P0_1(std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * F2[2] * M1 * S3[2]; 
  F1[3] = denom * cI * F2[3] * M1 * S3[2]; 
  F1[4] = denom * cI * S3[2] * (F2[2] * (P1[3] - P1[0]) + F2[3] * (P1[1] + cI *
      (P1[2])));
  F1[5] = denom * - cI * S3[2] * (F2[2] * (+cI * (P1[2]) - P1[1]) + F2[3] *
      (P1[0] + P1[3]));
}


void FFS5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * - S3[2] * (+cI * (TMP1 + TMP2)); 
}


void FFS4C1_1(std::complex<double> F1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * S3[2] * (F1[4] * (P2[0] + P2[3]) + F1[5] * (P2[1] + cI
      * (P2[2])));
  F2[3] = denom * cI * S3[2] * (F1[4] * (+cI * (P2[2]) - P2[1]) + F1[5] *
      (P2[3] - P2[0]));
  F2[4] = denom * cI * F1[4] * M2 * S3[2]; 
  F2[5] = denom * cI * F1[5] * M2 * S3[2]; 
}


void VVVVS3_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> S5[], std::complex<double>
    COUP, double M4, double W4, std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP34; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0] + S5[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1] + S5[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * S5[2] * (OM4 * P4[0] * (-cI * (TMP21 * TMP37) + cI * (TMP25 *
      TMP34)) + (-cI * (V2[2] * TMP25) + cI * (V3[2] * TMP21)));
  V4[3] = denom * S5[2] * (OM4 * P4[1] * (-cI * (TMP21 * TMP37) + cI * (TMP25 *
      TMP34)) + (-cI * (V2[3] * TMP25) + cI * (V3[3] * TMP21)));
  V4[4] = denom * S5[2] * (OM4 * P4[2] * (-cI * (TMP21 * TMP37) + cI * (TMP25 *
      TMP34)) + (-cI * (V2[4] * TMP25) + cI * (V3[4] * TMP21)));
  V4[5] = denom * S5[2] * (OM4 * P4[3] * (-cI * (TMP21 * TMP37) + cI * (TMP25 *
      TMP34)) + (-cI * (V2[5] * TMP25) + cI * (V3[5] * TMP21)));
}


void VVVV4_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP30; 
  std::complex<double> TMP35; 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  vertex = COUP * (-cI * (TMP25 * TMP30) + cI * (TMP21 * TMP35)); 
}


void SSS1_2(std::complex<double> S1[], std::complex<double> S3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> S2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  S2[0] = +S1[0] + S3[0]; 
  S2[1] = +S1[1] + S3[1]; 
  P2[0] = -S2[0].real(); 
  P2[1] = -S2[1].real(); 
  P2[2] = -S2[1].imag(); 
  P2[3] = -S2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  S2[2] = denom * cI * S3[2] * S1[2]; 
}


void VVS2_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP19; 
  std::complex<double> TMP21; 
  std::complex<double> TMP22; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  vertex = COUP * S3[2] * (-cI * (TMP17 * TMP22) + cI * (TMP19 * TMP21)); 
}


void FFV1_3(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP4; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP4 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] +
      F2[3] * F1[5] - P3[0] * OM3 * TMP4);
  V3[3] = denom * (-cI) * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] -
      F2[4] * F1[3] - P3[1] * OM3 * TMP4);
  V3[4] = denom * (-cI) * (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) + cI * (F2[4]
      * F1[3] + F2[3] * F1[4]) - P3[2] * OM3 * TMP4);
  V3[5] = denom * (-cI) * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] -
      F2[3] * F1[5] - P3[3] * OM3 * TMP4);
}


void FFV3C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP14; 
  std::complex<double> TMP15; 
  TMP15 = (-1.) * (F1[4] * (F2[2] * (V3[2] + V3[5]) + F2[3] * (V3[3] - cI *
      (V3[4]))) + F1[5] * (F2[2] * (V3[3] + cI * (V3[4])) + F2[3] * (V3[2] -
      V3[5])));
  TMP14 = (-1.) * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + F1[3] * (F2[4] * (-1.) * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])));
  vertex = COUP * (-cI * (TMP14) + 2. * cI * (TMP15)); 
}


void FFV5C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * 4. * cI * (F2[2] * (P1[0] * (-1.) * (V3[2] + V3[5]) + (P1[1]
      * (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] *
      (V3[2] + V3[5])))) + (+1./4. * (M1 * (F2[5] * (V3[3] - cI * (V3[4])) + 4.
      * (F2[4] * 1./4. * (V3[5] - V3[2])))) + F2[3] * (P1[0] * (+cI * (V3[4]) -
      V3[3]) + (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI *
      (V3[5])) + P1[3] * (V3[3] - cI * (V3[4])))))));
  F1[3] = denom * 4. * cI * (F2[2] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + (+1./4. * (M1 * (+4. * (F2[4] * 1./4. *
      (V3[3] + cI * (V3[4]))) - F2[5] * (V3[2] + V3[5]))) + F2[3] * (P1[0] *
      (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI
      * (V3[3])) + P1[3] * (V3[5] - V3[2]))))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + (F2[5] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + M1 * (F2[2] * 4. * (V3[2] + V3[5]) + 4. * (F2[3]
      * (V3[3] - cI * (V3[4]))))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (-1.) * (V3[3] + cI * (V3[4])) +
      (P1[1] * (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[5]) + cI * (V3[2])) +
      P1[3] * (V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (V3[2] + V3[5]) +
      (P1[1] * (+cI * (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI *
      (V3[3])) - P1[3] * (V3[2] + V3[5])))) + M1 * (F2[2] * 4. * (V3[3] + cI *
      (V3[4])) + 4. * (F2[3] * (V3[2] - V3[5])))));
}


void VVVS1_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> S4[], std::complex<double> COUP, double M2, double W2,
    std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> TMP19; 
  std::complex<double> TMP23; 
  std::complex<double> TMP24; 
  std::complex<double> TMP25; 
  std::complex<double> TMP28; 
  std::complex<double> denom; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V2[0] = +V1[0] + V3[0] + S4[0]; 
  V2[1] = +V1[1] + V3[1] + S4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP28 = (P2[0] * P3[0] - P2[1] * P3[1] - P2[2] * P3[2] - P2[3] * P3[3]); 
  TMP19 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP24 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S4[2] * (OM2 * P2[0] * (TMP25 * (-cI * (TMP19) + cI *
      (TMP28)) + (-cI * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 *
      (-cI * (P3[0]) + cI * (P1[0])) + (V1[2] * (-cI * (TMP23) + cI * (TMP24))
      + V3[2] * (-cI * (TMP17) + cI * (TMP18)))));
  V2[3] = denom * S4[2] * (OM2 * P2[1] * (TMP25 * (-cI * (TMP19) + cI *
      (TMP28)) + (-cI * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 *
      (-cI * (P3[1]) + cI * (P1[1])) + (V1[3] * (-cI * (TMP23) + cI * (TMP24))
      + V3[3] * (-cI * (TMP17) + cI * (TMP18)))));
  V2[4] = denom * S4[2] * (OM2 * P2[2] * (TMP25 * (-cI * (TMP19) + cI *
      (TMP28)) + (-cI * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 *
      (-cI * (P3[2]) + cI * (P1[2])) + (V1[4] * (-cI * (TMP23) + cI * (TMP24))
      + V3[4] * (-cI * (TMP17) + cI * (TMP18)))));
  V2[5] = denom * S4[2] * (OM2 * P2[3] * (TMP25 * (-cI * (TMP19) + cI *
      (TMP28)) + (-cI * (TMP18 * TMP24) + cI * (TMP17 * TMP23))) + (TMP25 *
      (-cI * (P3[3]) + cI * (P1[3])) + (V1[5] * (-cI * (TMP23) + cI * (TMP24))
      + V3[5] * (-cI * (TMP17) + cI * (TMP18)))));
}


void VVVVS3_2(std::complex<double> V1[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> S5[], std::complex<double>
    COUP, double M2, double W2, std::complex<double> V2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM2; 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP25; 
  std::complex<double> TMP32; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./(M2 * M2); 
  V2[0] = +V1[0] + V3[0] + V4[0] + S5[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1] + S5[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP32 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  V2[2] = denom * S5[2] * (OM2 * P2[0] * (-cI * (TMP17 * TMP35) + cI * (TMP25 *
      TMP32)) + (-cI * (TMP25 * V4[2]) + cI * (V1[2] * TMP35)));
  V2[3] = denom * S5[2] * (OM2 * P2[1] * (-cI * (TMP17 * TMP35) + cI * (TMP25 *
      TMP32)) + (-cI * (TMP25 * V4[3]) + cI * (V1[3] * TMP35)));
  V2[4] = denom * S5[2] * (OM2 * P2[2] * (-cI * (TMP17 * TMP35) + cI * (TMP25 *
      TMP32)) + (-cI * (TMP25 * V4[4]) + cI * (V1[4] * TMP35)));
  V2[5] = denom * S5[2] * (OM2 * P2[3] * (-cI * (TMP17 * TMP35) + cI * (TMP25 *
      TMP32)) + (-cI * (TMP25 * V4[5]) + cI * (V1[5] * TMP35)));
}


void VVVV2_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP23; 
  std::complex<double> TMP27; 
  std::complex<double> TMP30; 
  std::complex<double> TMP31; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP23 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-2. * cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31 + TMP23 * TMP30)) + (-cI * (TMP27 * V4[2] + V3[2] * TMP30) + 2. *
      cI * (V2[2] * TMP35)));
  V1[3] = denom * (OM1 * P1[1] * (-2. * cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31 + TMP23 * TMP30)) + (-cI * (TMP27 * V4[3] + V3[3] * TMP30) + 2. *
      cI * (V2[3] * TMP35)));
  V1[4] = denom * (OM1 * P1[2] * (-2. * cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31 + TMP23 * TMP30)) + (-cI * (TMP27 * V4[4] + V3[4] * TMP30) + 2. *
      cI * (V2[4] * TMP35)));
  V1[5] = denom * (OM1 * P1[3] * (-2. * cI * (TMP22 * TMP35) + cI * (TMP27 *
      TMP31 + TMP23 * TMP30)) + (-cI * (TMP27 * V4[5] + V3[5] * TMP30) + 2. *
      cI * (V2[5] * TMP35)));
}


void VVVV3_1(std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, double M1, double W1,
    std::complex<double> V1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM1; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP27; 
  std::complex<double> TMP31; 
  std::complex<double> TMP35; 
  std::complex<double> denom; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./(M1 * M1); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP35 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP31 = (P1[0] * V4[2] - P1[1] * V4[3] - P1[2] * V4[4] - P1[3] * V4[5]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP22 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-cI * (TMP22 * TMP35) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[2]) + cI * (V2[2] * TMP35)));
  V1[3] = denom * (OM1 * P1[1] * (-cI * (TMP22 * TMP35) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[3]) + cI * (V2[3] * TMP35)));
  V1[4] = denom * (OM1 * P1[2] * (-cI * (TMP22 * TMP35) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[4]) + cI * (V2[4] * TMP35)));
  V1[5] = denom * (OM1 * P1[3] * (-cI * (TMP22 * TMP35) + cI * (TMP27 * TMP31))
      + (-cI * (TMP27 * V4[5]) + cI * (V2[5] * TMP35)));
}


void VVSS1_0(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> S3[], std::complex<double> S4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * - cI * TMP21 * S4[2] * S3[2]; 
}


void VSS1_3(std::complex<double> V1[], std::complex<double> S2[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP17; 
  std::complex<double> TMP18; 
  std::complex<double> denom; 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  S3[0] = +V1[0] + S2[0]; 
  S3[1] = +V1[1] + S2[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP17 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  S3[2] = denom * S2[2] * (-cI * (TMP18) + cI * (TMP17)); 
}


void FFS3C1_0(std::complex<double> F2[], std::complex<double> F1[],
    std::complex<double> S3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP2; 
  TMP2 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  vertex = COUP * S3[2] * (-cI * (TMP1) + cI * (TMP2)); 
}


void FFV5_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
    & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP11; 
  std::complex<double> TMP12; 
  TMP11 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP12 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  vertex = COUP * (-1.) * (+cI * (TMP11) + 4. * cI * (TMP12)); 
}


void FFV2_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])));
  F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
      - V3[2]));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))));
  F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
      (-1.) * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] *
      (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3]
      * (V3[2] - V3[5])))));
}

void FFV2_5_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_1(F2, V3, COUP1, M1, W1, F1); 
  FFV5_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_3_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_1(F2, V3, COUP1, M1, W1, F1); 
  FFV3_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_4_1(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2_1(F2, V3, COUP1, M1, W1, F1); 
  FFV4_1(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV1_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP26; 
  std::complex<double> TMP29; 
  std::complex<double> TMP30; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP29 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP26 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (-cI * (TMP18 * TMP30) + cI * (TMP26 * TMP29))
      + (-cI * (V2[2] * TMP29) + cI * (V1[2] * TMP30)));
  V3[3] = denom * (OM3 * P3[1] * (-cI * (TMP18 * TMP30) + cI * (TMP26 * TMP29))
      + (-cI * (V2[3] * TMP29) + cI * (V1[3] * TMP30)));
  V3[4] = denom * (OM3 * P3[2] * (-cI * (TMP18 * TMP30) + cI * (TMP26 * TMP29))
      + (-cI * (V2[4] * TMP29) + cI * (V1[4] * TMP30)));
  V3[5] = denom * (OM3 * P3[3] * (-cI * (TMP18 * TMP30) + cI * (TMP26 * TMP29))
      + (-cI * (V2[5] * TMP29) + cI * (V1[5] * TMP30)));
}


void SSSS1_4(std::complex<double> S1[], std::complex<double> S2[],
    std::complex<double> S3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> S4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P4[4]; 
  std::complex<double> denom; 
  S4[0] = +S1[0] + S2[0] + S3[0]; 
  S4[1] = +S1[1] + S2[1] + S3[1]; 
  P4[0] = -S4[0].real(); 
  P4[1] = -S4[1].real(); 
  P4[2] = -S4[1].imag(); 
  P4[3] = -S4[0].imag(); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  S4[2] = denom * cI * S3[2] * S2[2] * S1[2]; 
}


void FFV2C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * M1 * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4])
      - V3[3]));
  F1[3] = denom * cI * M1 * (F2[4] * (V3[3] + cI * (V3[4])) - F2[5] * (V3[2] +
      V3[5]));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] - V3[5]) + (P1[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3] *
      (V3[2] - V3[5])))) + F2[5] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + F2[5] * (P1[0] * (-1.) * (V3[2] + V3[5]) +
      (P1[1] * (V3[3] - cI * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3]
      * (V3[2] + V3[5])))));
}

void FFV2_5C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_2(F2, V3, COUP1, M1, W1, F1); 
  FFV5C1_2(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_3C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_2(F2, V3, COUP1, M1, W1, F1); 
  FFV3C1_2(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}
void FFV2_4C1_2(std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, double M1, double
    W1, std::complex<double> F1[])
{
  std::complex<double> Ftmp[6]; 
  int i; 
  FFV2C1_2(F2, V3, COUP1, M1, W1, F1); 
  FFV4C1_2(F2, V3, COUP2, M1, W1, Ftmp); 
  i = 2; 
  while (i < 6)
  {
    F1[i] = F1[i] + Ftmp[i]; 
    i++; 
  }
}

void VVVV5_4(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V3[], std::complex<double> COUP, double M4, double W4,
    std::complex<double> V4[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM4; 
  double P4[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP25; 
  std::complex<double> TMP27; 
  std::complex<double> TMP33; 
  std::complex<double> TMP34; 
  std::complex<double> TMP37; 
  std::complex<double> denom; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./(M4 * M4); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP34 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP27 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP37 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP33 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  denom = COUP/((P4[0] * P4[0]) - (P4[1] * P4[1]) - (P4[2] * P4[2]) - (P4[3] *
      P4[3]) - M4 * (M4 - cI * W4));
  V4[2] = denom * 1./2. * (OM4 * - P4[0] * (-2. * cI * (TMP27 * TMP33) + cI *
      (TMP25 * TMP34 + TMP21 * TMP37)) + (-2. * cI * (V1[2] * TMP27) + cI *
      (V2[2] * TMP25 + V3[2] * TMP21)));
  V4[3] = denom * 1./2. * (OM4 * - P4[1] * (-2. * cI * (TMP27 * TMP33) + cI *
      (TMP25 * TMP34 + TMP21 * TMP37)) + (-2. * cI * (V1[3] * TMP27) + cI *
      (V2[3] * TMP25 + V3[3] * TMP21)));
  V4[4] = denom * 1./2. * (OM4 * - P4[2] * (-2. * cI * (TMP27 * TMP33) + cI *
      (TMP25 * TMP34 + TMP21 * TMP37)) + (-2. * cI * (V1[4] * TMP27) + cI *
      (V2[4] * TMP25 + V3[4] * TMP21)));
  V4[5] = denom * 1./2. * (OM4 * - P4[3] * (-2. * cI * (TMP27 * TMP33) + cI *
      (TMP25 * TMP34 + TMP21 * TMP37)) + (-2. * cI * (V1[5] * TMP27) + cI *
      (V2[5] * TMP25 + V3[5] * TMP21)));
}


void VVVV4_3(std::complex<double> V1[], std::complex<double> V2[],
    std::complex<double> V4[], std::complex<double> COUP, double M3, double W3,
    std::complex<double> V3[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  std::complex<double> TMP18; 
  std::complex<double> TMP21; 
  std::complex<double> TMP30; 
  std::complex<double> TMP36; 
  std::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +V1[0] + V2[0] + V4[0]; 
  V3[1] = +V1[1] + V2[1] + V4[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP21 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP30 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP36 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]); 
  TMP18 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
      P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (OM3 * P3[0] * (-cI * (TMP21 * TMP36) + cI * (TMP18 * TMP30))
      + (-cI * (V1[2] * TMP30) + cI * (TMP21 * V4[2])));
  V3[3] = denom * (OM3 * P3[1] * (-cI * (TMP21 * TMP36) + cI * (TMP18 * TMP30))
      + (-cI * (V1[3] * TMP30) + cI * (TMP21 * V4[3])));
  V3[4] = denom * (OM3 * P3[2] * (-cI * (TMP21 * TMP36) + cI * (TMP18 * TMP30))
      + (-cI * (V1[4] * TMP30) + cI * (TMP21 * V4[4])));
  V3[5] = denom * (OM3 * P3[3] * (-cI * (TMP21 * TMP36) + cI * (TMP18 * TMP30))
      + (-cI * (V1[5] * TMP30) + cI * (TMP21 * V4[5])));
}


void FFV4_2(std::complex<double> F1[], std::complex<double> V3[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> F2[])
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] *
      P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * 2. * (V3[2] - V3[5]) + 2. * (F1[5]
      * (+cI * (V3[4]) - V3[3])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (-1.) * (V3[2] + V3[5]) + (P2[2] * (-1.) * (+cI * (V3[2] + V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] - V3[5]) +
      (P2[1] * (+cI * (V3[4]) - V3[3]) + (P2[2] * (-1.) * (V3[4] + cI *
      (V3[3])) + P2[3] * (V3[2] - V3[5])))) + M2 * (F1[4] * (-2.) * (V3[3] + cI
      * (V3[4])) + 2. * (F1[5] * (V3[2] + V3[5])))));
  F2[4] = denom * (-2. * cI) * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + (+1./2. * (M2 * (F1[3] * (+cI * (V3[4]) - V3[3]) + 2. *
      (F1[2] * (-1./2.) * (V3[2] + V3[5])))) + F1[5] * (P2[0] * (V3[3] - cI *
      (V3[4])) + (P2[1] * (-1.) * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] +
      V3[5])) + P2[3] * (V3[3] - cI * (V3[4])))))));
  F2[5] = denom * (-2. * cI) * (F1[4] * (P2[0] * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) -
      P2[3] * (V3[3] + cI * (V3[4]))))) + (+1./2. * (M2 * (F1[3] * (V3[5] -
      V3[2]) + 2. * (F1[2] * (-1./2.) * (V3[3] + cI * (V3[4]))))) + F1[5] *
      (P2[0] * (-1.) * (V3[2] + V3[5]) + (P2[1] * (V3[3] - cI * (V3[4])) +
      (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[2] + V3[5]))))));
}


}  // end namespace $(namespace)s_MDMSM

