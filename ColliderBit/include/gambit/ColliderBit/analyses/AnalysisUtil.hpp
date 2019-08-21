//
// Created by dsteiner on 31/07/18.
//

#include "gambit/ColliderBit/Utils.hpp"
#include "HEPUtils/Event.h"

namespace Gambit
{
namespace ColliderBit
{
class AnalysisUtil
{

public:
  static bool sortParticlesByPt(HEPUtils::Particle *particle1, HEPUtils::Particle *particle2)
  {
    return (particle1->pT() > particle2->pT());
  }

  static bool sortJetsByPt(HEPUtils::Jet *jet1, HEPUtils::Jet *jet2)
  {
    return (jet1->pT() > jet2->pT());
  }

  static std::vector<HEPUtils::Jet*> filterPtEta(std::vector<HEPUtils::Jet*> jets, double pT, double absEta)
  {
    std::vector<HEPUtils::Jet*> outJets;
    for (HEPUtils::Jet* jet : jets)
    {
      if (jet->pT() > pT && jet->abseta() < absEta)
      {
        outJets.push_back(jet);
      }
    }
    return outJets;
  }

  static std::vector<HEPUtils::Particle*> filterPtEta(std::vector<HEPUtils::Particle*> particles, double pT, double absEta)
  {
    std::vector<HEPUtils::Particle*> outParticles;
    for (HEPUtils::Particle* particle : particles)
    {
      if (particle->pT() > pT && particle->abseta() < absEta)
      {
        outParticles.push_back(particle);
      }
    }
    return outParticles;
  }

  static std::vector<HEPUtils::Jet*> jetLeptonOverlapRemoval(std::vector<HEPUtils::Jet*> jets, std::vector<HEPUtils::Particle*> leptons, double dR)
  {
    std::vector<HEPUtils::Jet*> outJets;
    for (HEPUtils::Jet* jet : jets)
    {
      bool overlap = false;
      for (HEPUtils::Particle* lepton : leptons)
      {
        double dRJetElectron = lepton->mom().deltaR_eta(jet->mom());
        if (fabs(dRJetElectron) < dR)
        {
          overlap = true;
        }
      }
      if (overlap) continue;
      outJets.push_back(jet);
    }
    return outJets;
  }


  static std::vector<HEPUtils::Particle*> leptonJetOverlapRemoval(std::vector<HEPUtils::Particle*> leptons, std::vector<HEPUtils::Jet*> jets, double dR)
  {
    std::vector<HEPUtils::Particle*> outLeptons;
    for (HEPUtils::Particle* lepton : leptons)
    {
      bool overlap = false;
      for (HEPUtils::Jet* jet : jets)
      {
        double dRLeptonJet = lepton->mom().deltaR_eta(jet->mom());
        if (fabs(dRLeptonJet) < dR)
        {
          overlap = true;
        }
      }
      if (overlap) continue;
      outLeptons.push_back(lepton);
    }
    return outLeptons;
  }

  static bool isSingleParticleTriggered(std::vector<HEPUtils::Particle*> particles, double pTrequirement)
  {
    for (HEPUtils::Particle* particle : particles)
    {
      if (particle->pT() > pTrequirement)
      {
        return true;
      }
    }
    return false;
  }

  static bool isMultipleParticleTriggered(std::vector<HEPUtils::Particle*> particles, std::vector<double> pTrequirements)
  {
    size_t numTriggers = 0;
    for (HEPUtils::Particle* particle : particles)
    {
      for (double pTrequirement : pTrequirements)
      {
        if (particle->pT() > pTrequirement)
        {
          numTriggers++;
          break;
        }
      }
    }
    return numTriggers >= pTrequirements.size();
  }

  static std::vector<HEPUtils::Particle*> getSortedLeptons(const std::vector<std::vector<HEPUtils::Particle*>> allLeptons)
  {
    std::vector<HEPUtils::Particle*> leptons;
    for (std::vector<HEPUtils::Particle*> setOfLeptons : allLeptons)
    {
      leptons.insert(leptons.end(), setOfLeptons.begin(), setOfLeptons.end());
    }
    std::sort(leptons.begin(), leptons.end(), sortParticlesByPt);
    return leptons;
  }

  static double dot2D(HEPUtils::P4 mom1, HEPUtils::P4 mom2)
  {
    return mom1.px() * mom2.px() + mom1.py() * mom2.py();
  }

  static std::vector<HEPUtils::Jet*> filterMaxEta(const std::vector<HEPUtils::Jet*>& jets, double maxAbsEta)
  {
    std::vector<HEPUtils::Jet*> outJets;
    for (HEPUtils::Jet* jet : jets)
    {
      if (jet->abseta() < maxAbsEta)
      {
        outJets.push_back(jet);
      }
    }
    return outJets;
  }


  static bool muonFilter7TeV(const std::vector<HEPUtils::Particle*>& muons)
  {
    double effProduct = 1.0;
    for (HEPUtils::Particle* muon : muons)
    {
      if (muon->abseta() < 1.05)
      {
        effProduct *= (1 - 0.7);
      }
      else
      {
        effProduct *= (1 - 0.95);
      }
    }
    double efficiency = 1 - effProduct;
    return random_bool(efficiency);
  }

  static bool oppositeSign(HEPUtils::Particle* a, HEPUtils::Particle* b)
  {
    return a->pid() * b->pid() < 0;
  }
};
}
}