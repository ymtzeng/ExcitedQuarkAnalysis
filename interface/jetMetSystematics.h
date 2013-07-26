#ifndef __JETMETSYSTEMATICS__
#define __JETMETSYSTEMATICS__

#include <string>
#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

class EvtInfoBranches;
class JetInfoBranches;
class LepInfoBranches;

using std::string;

class jetMetSystematics {

   public:

      jetMetSystematics(const edm::ParameterSet&,
                        EvtInfoBranches &evt, JetInfoBranches &jet, LepInfoBranches &lep);

      void scale();

      enum Systematic{
         NOSYS,
         JES,
         JER,
         MET,
         UNC,//unclustered energy
         nSystematics
      };
      static const string sType[];

   private:
      double jesScale(const int index, const double sigma);
      double jerScale(const int index, const double sigma);
      void   setType (std::string s);

//       const int nEtaBins;
//       const int nPtBins;
      std::string _jecFile;
      std::string stype;
      Systematic _type;

      EvtInfoBranches* _evt;
      JetInfoBranches* _jets;
      LepInfoBranches* _leps;

//       double *etaBins, *ptBins;
//       double *jesPlus, *jesMinus;

//       double c_sw;// = 0.015;
//       double e_pu;// = 0.75;
//       //double e_pu123;// = 0.2;
//       double ja;//   = 0.8;
//       double avgpu;// = 2.2;
//       double c_bjes1;// = 0.02;
//       double c_bjes2;// = 0.03;

      double _scale;

      double _minDRljet;

      //JES scaling info
      JetCorrectionUncertainty *_jesSigma;
      std::string _jesUncType;

      //JER scaling info
      std::vector<double> _jerEta;//define max edge of eta bin
      std::vector<double> _jerNominal;//measured data/MC ratio
      std::vector<double> _jerSigmaSym, _jerSigmaNeg, _jerSigmaPos;//uncertainties on ratio
      double _jerMinGenJetPt;

      bool   _debug;

};

#endif
