///
///  \author Rose Kudzman-Blais
///  \date 2017 May
///
///  *********************************************


#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/analyses/Perf_Plot.hpp"

using namespace std;

namespace Gambit {
  namespace ColliderBit {

    class Analysis_CMS_13TeV_1LEPbb_36invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSRA, _numSRB; 
      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      size_t NCUTS;

      Perf_Plot* plots;
      ofstream cutflowFile;
      string analysisRunName;

    public:

      Analysis_CMS_13TeV_1LEPbb_36invfb() {

        _numSRA=0;
        _numSRB=0;

        NCUTS=10;
        set_luminosity(35.9);

        for (size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVector_str.push_back("");
        }

        analysisRunName = "CMS_13TeV_1LEPbb_36invfb_test";
        vector<const char*> variables = {"met","mCT","mbb"};
        plots = new Perf_Plot(analysisRunName, &variables);
	
      }


      void analyze(const HEPUtils::Event* event) {
	HEPUtilsAnalysis::analyze(event);

        // Missing energy
        double met = event->met();

        // Baseline objects
        vector<HEPUtils::Particle*> baselineLeptons;
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          baselineElectrons.push_back(electron);
	  baselineLeptons.push_back(electron);
        }
        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
          baselineMuons.push_back(muon);
	  baselineLeptons.push_back(muon);
        }
        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          baselineJets.push_back(jet);
        }

        // Signal objects
        vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Particle*> signalElectrons;
        vector<HEPUtils::Particle*> signalMuons;
        vector<HEPUtils::Jet*> signalJets;   
        vector<HEPUtils::Jet*> signalBJets;

	// Electrons
	for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
	  if (baselineElectrons.at(iEl)->pT() > 30. && fabs(baselineElectrons.at(iEl)->eta()) < 1.44) {
	      signalElectrons.push_back(baselineElectrons.at(iEl));
	      signalLeptons.push_back(baselineElectrons.at(iEl));
	  }
	} 
       
        //Muons
        for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
          if (baselineMuons.at(iMu)->pT() > 25. && fabs(baselineMuons.at(iMu)->eta()) < 2.1) {
	      signalMuons.push_back(baselineMuons.at(iMu));
	      signalLeptons.push_back(baselineMuons.at(iMu));
          }
        }

	//Jets 
        const vector<double>  a = {0,10.};
        const vector<double>  b = {0,10000.};
        const vector<double> c = {0.65};
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);

	for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
           if (baselineJets.at(iJet)->pT() > 30. && fabs(baselineJets.at(iJet)->eta()) < 2.40) {
            bool hasTag=has_tag(_eff2d, baselineJets.at(iJet)->eta(), baselineJets.at(iJet)->pT()); 
	    signalJets.push_back(baselineJets.at(iJet));                
	    if (baselineJets.at(iJet)->btag() && hasTag)signalBJets.push_back(baselineJets.at(iJet));
	  }
         }

        int nSignalLeptons = signalElectrons.size() + signalMuons.size();
        int nSignalJets = signalJets.size();
        int nSignalBJets = signalBJets.size();

	//Preselection
	bool preselection=false;	
 
        //2nd lepton veto
        bool lepton2_veto=true;
	if (nSignalLeptons>0) {
	  for (size_t iLe=0;iLe<baselineLeptons.size();iLe++) {
	    if (baselineLeptons.at(iLe)->pT()>5 && baselineLeptons.at(iLe)!=signalLeptons.at(0))lepton2_veto=false;
	  }
	} 

	//Tau veto
	bool tau_veto=true;
	for (HEPUtils::Particle* tau : event->taus()) {
	  if (tau->pT()>20)tau_veto=false;
	}	

        if (nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2)preselection=true;

        //Signal regions
        double mCT=0;
        double mbb=0;
	double mT=0;
	
	if (nSignalBJets>=2) {
	  mCT=sqrt(2*signalBJets.at(0)->pT()*signalBJets.at(1)->pT()*(1+cos(signalBJets.at(0)->phi()-signalBJets.at(1)->phi())));
          mbb = (signalBJets.at(0)->mom()+signalBJets.at(1)->mom()).m();
	}
	if (signalElectrons.size()>0) {
	  mT = sqrt(2*signalElectrons.at(0)->pT()*met*(1-cos(signalElectrons.at(0)->phi()-event->missingmom().phi())));
	}
        if (signalMuons.size()>0) {
          mT = sqrt(2*signalMuons.at(0)->pT()*met*(1-cos(signalMuons.at(0)->phi()-event->missingmom().phi()))); 
        }

        if (preselection && mbb>90 && mbb<150 && mCT>170. && met>125. && mT>150.) {
          //SRA
          if (met>125 && met<200)_numSRA++;
          //SRB
          if (met>200)_numSRB++;   
	}

        if (preselection) {
          vector<double> variables = {met, mCT, mbb};
          plots->fill(&variables);
        }             

        cutFlowVector_str[0] = "All events";
        cutFlowVector_str[1] = ">= 1 signal lepton; > 50 GeV";
        cutFlowVector_str[2] = "2nd lepton veto";
        cutFlowVector_str[3] = "Tau veto";
        cutFlowVector_str[4] = "2 jets";
        cutFlowVector_str[5] = "2 bjets";
        cutFlowVector_str[6] = "90 < mbb < 150 GeV";
        cutFlowVector_str[7] = "mCT > 170 GeV";
        cutFlowVector_str[8] = "MET > 125 GeV";
        cutFlowVector_str[9] = "mT > 150 GeV";

        for (size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||

	     (j==1 && nSignalLeptons>=1 && met>=50) ||

	     (j==2 && nSignalLeptons>=1 && met>=50 && lepton2_veto) ||

	     (j==3 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto) ||

	     (j==4 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2) ||

	     (j==5 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2) ||

	     (j==6 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150) ||

	     (j==7 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170.) ||            
 
             (j==8 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170. && met>125.) ||
             
             (j==9 && nSignalLeptons>=1 && met>=50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170. && met>125. && mT>150.) )

            cutFlowVector[j]++;
	}

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_CMS_13TeV_1LEPbb_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_1LEPbb_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (size_t j = 0; j < NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }
        _numSRA += specificOther->_numSRA;
        _numSRB += specificOther->_numSRB;
      }


      void collect_results() {

        string path = "ColliderBit/results/cutflow_";
        path.append(analysisRunName);
        path.append(".txt");
        cutflowFile.open(path.c_str()); 
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile << "CUT FLOW: CMS 1 lepton, 2 bjets paper "<<endl;
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------"<<endl;

        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str[j].c_str() << setw(20) << cutFlowVector[j] << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile.close();

        plots->createFile();

        SignalRegionData results_SRA;
        results_SRA.analysis_name = "Analysis_CMS_13TeV_1LEPbb_36invfb";
        results_SRA.sr_label = "SRA";
        results_SRA.n_observed = 11.;
        results_SRA.n_background = 7.5; 
        results_SRA.background_sys = 2.5;
        results_SRA.signal_sys = 0.; 
        results_SRA.n_signal = _numSRA;
        add_result(results_SRA);

        SignalRegionData results_SRB;
        results_SRB.analysis_name = "Analysis_CMS_13TeV_1LEPbb_36invfb";
        results_SRB.sr_label = "SRB";
        results_SRB.n_observed = 7.;
        results_SRB.n_background = 8.7; 
        results_SRB.background_sys = 2.2;
        results_SRB.signal_sys = 0.; 
        results_SRB.n_signal = _numSRB;
        add_result(results_SRB);

      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_1LEPbb_36invfb)


  }
}
