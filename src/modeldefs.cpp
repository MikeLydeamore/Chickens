#include <math.h>
#include <map>
#include <vector>
#include <Rcpp.h>
#include <algorithm>
#include "MarkovChainSimulator/MarkovChain/MarkovChain.hpp"

typedef std::map<std::string, double> stringmap;

class WithinPatchParameters
{
public:
  stringmap mX0;
  std::vector<double> mN;
  std::vector<double> mDelta;
  double mY;
  double mX;
  stringmap mAlpha;
  stringmap mBeta;
  double mSigma;
  double mGamma;
  double mNEgg;
  double mQ; 
  double mW;
  double mK;
  double mBr;
  
  WithinPatchParameters(stringmap x0, std::vector<double> n, std::vector<double> delta, double y, double x, stringmap alpha, stringmap beta, double sigma, double gamma, double nEgg, double q, double w, double K, double br) :
    mX0(x0), mN(n), mDelta(delta), mY(y), mX(x), mAlpha(alpha), mBeta(beta), mSigma(sigma), mGamma(gamma), mNEgg(nEgg), mQ(q), mW(w), mK(K), mBr(br) {}
  
  WithinPatchParameters() {}
};

class ModelChickenFlu {
  
private:   
  static double eggHatchingRate(state_values states, parameter_map parameters) {
    //Calculate n_1*E and total sum of death rates.
    std::string patchname;
    if (parameters["Sc"] == 1)
      patchname = "Sc";
    if (parameters["Ns"] == 1)
      patchname = "Ns";
    if (parameters["Es"] == 1)
      patchname = "Es";
    if (parameters["Bs"] == 1)
      patchname = "Bs";
    
    double population_size = 0;
    for (state_values::iterator it = states.begin() ; it != states.end() ; it++) 
    {
      std::string s = it->first;
      std::string delimeter = ".";
      std::string current_patch = s.substr(0, s.find(delimeter));
      s.erase(0, s.find(delimeter) + delimeter.length());
      std::string state = s.substr(0, s.find(delimeter));
      size_t n = std::count(s.begin(), s.end(), '.');
      if (n > 0 && current_patch == patchname) 
      {
        population_size = population_size + it->second;
      }
    }
    
    double n1E = parameters["n1"]*states[patchname+".E"];
    double totaldeaths = 0;
    
    std::vector<std::string> disease_states = {"S","E","I"};
    for (std::string disease_state : disease_states)
    {
      totaldeaths += states[patchname+".Ch."+disease_state] * parameters["delta1"];
      totaldeaths += states[patchname+".eG."+disease_state] * parameters["delta2"];
      totaldeaths += states[patchname+".lG."+disease_state] * parameters["alphaLG"];
      totaldeaths += states[patchname+".lG."+disease_state] * (1-parameters["alphaLG"])*parameters["delta3"];
      totaldeaths += states[patchname+".He."+disease_state] * parameters["alphaHe"];
      totaldeaths += states[patchname+".He."+disease_state] * (1-parameters["alphaHe"])*parameters["delta3"];
      totaldeaths += states[patchname+".Rs."+disease_state] * parameters["alphaRs"];
      totaldeaths += states[patchname+".Rs."+disease_state] * (1-parameters["alphaRs"])*parameters["delta4"];
    }
    
    double alpha = 0;
    double rho = 0;
    if (population_size + n1E - totaldeaths <= parameters["K"])
    {
      alpha = 1;
      rho = parameters["K"] - population_size - n1E + totaldeaths;
    }
    else
    {
      rho = 0;
      alpha = (parameters["K"] - population_size + totaldeaths)/n1E;
      if (alpha > 1)
        alpha = 1;
      if (alpha < 0 || std::isinf(alpha))
        alpha = 0;
    }
    
    if (parameters["returnrho"])
    {
      if (parameters["Ch.S"] == 1)
        return (parameters["y"]*rho); 
      else
        return ((1-parameters["y"])*rho);
    }
    
    if (parameters["returnhatching"])
      return (alpha*n1E);
    
    if (parameters["returnsold"])
      return ((1-alpha)*n1E);
  }
  
  std::map<std::string, WithinPatchParameters> mPatchParams; 
  std::vector<std::string> mPatchNames;
  
public:
  
  ModelChickenFlu(std::vector<std::string> patchNames, std::map<std::string, WithinPatchParameters> patchParams) :
  mPatchParams(patchParams), mPatchNames(patchNames) {}
  
  void setupModel(MarkovChain &rChain) {
    
    const std::vector<std::string> disease_states = {"S", "E", "I"};
    const std::vector<std::string> demographic_states = {"Ch","eG","lG","He","Rs"};
    
    std::vector<std::string> sub_infectious = {"Ch.I", "eG.I", "lG.I", "He.I", "Rs.I"};
    std::vector<std::string> population_infectious;
    for (std::string patchName : mPatchNames)
    {
      for (std::string sub : sub_infectious)
        population_infectious.push_back(patchName+"."+sub);
    }
    
    
    std::vector<std::string> population_states;
    for (std::string patchName : mPatchNames)
    {
      //States!
      rChain.addState(patchName+".E", mPatchParams[patchName].mX0["E"]);
      
      
      std::vector<std::string> within_patch_population_states;
      double initial_size = 0;
      
      for (std::string demographic_state : demographic_states) 
      {
        for (std::string disease_state : disease_states)    
        {
          std::string state_name = patchName + "." + demographic_state + "." + disease_state;
          
          rChain.addState(state_name, mPatchParams[patchName].mX0[demographic_state+"."+disease_state]);
          
          population_states.push_back(state_name);
          within_patch_population_states.push_back(state_name);
          initial_size += mPatchParams[patchName].mX0[demographic_state+"."+disease_state];
        }
      }
      
      //Aging transitions
      parameter_map hatchingParameters;
      hatchingParameters[patchName] = 1;
      hatchingParameters["n1"]=mPatchParams[patchName].mN[1];
      hatchingParameters["delta1"]=mPatchParams[patchName].mDelta[1];
      hatchingParameters["delta2"]=mPatchParams[patchName].mDelta[2];
      hatchingParameters["delta3"]=mPatchParams[patchName].mDelta[3];
      hatchingParameters["delta4"]=mPatchParams[patchName].mDelta[4];
      hatchingParameters["alphaLG"]=mPatchParams[patchName].mAlpha["lG"];
      hatchingParameters["alphaHe"]=mPatchParams[patchName].mAlpha["He"];
      hatchingParameters["alphaRs"]=mPatchParams[patchName].mAlpha["Rs"];
      hatchingParameters["returnhatching"]=1;
      hatchingParameters["K"]=initial_size;
      hatchingParameters["y"]=mPatchParams[patchName].mY;
      hatchingParameters["Ch.S"]=1;
      
      rChain.addTransition(new TransitionCustom(patchName+".E", patchName+".Ch.S", hatchingParameters, *eggHatchingRate));

      
      hatchingParameters["returnsold"]=1; hatchingParameters["returnhatching"]=0;
      rChain.addTransition(new TransitionCustomToVoid(patchName+".E", hatchingParameters, *eggHatchingRate));
      
      hatchingParameters["returnrho"]=1; hatchingParameters["returnsold"]=0;
      TransitionCustomFromVoid* import_chicks = new TransitionCustomFromVoid(patchName+".Ch.S", hatchingParameters, *eggHatchingRate);
      import_chicks->addCounter(patchName+".importedChicks");
      rChain.addTransition(import_chicks);

      hatchingParameters["Ch.S"]=0;
      TransitionCustomFromVoid* import_hens = new TransitionCustomFromVoid(patchName+".He.S", hatchingParameters, *eggHatchingRate);
      import_hens->addCounter(patchName+".importedHens");
      rChain.addTransition(import_hens);
      
      
      //Need (1-y) into Hens
      for (std::string disease_state : disease_states) 
      {
        //Ageing
        rChain.addTransition(new TransitionIndividual(patchName+".Ch."+disease_state, patchName+".eG."+disease_state, mPatchParams[patchName].mN[1]));
        rChain.addTransition(new TransitionIndividual(patchName+".eG."+disease_state, patchName+".lG."+disease_state, 2*mPatchParams[patchName].mN[2]));
        rChain.addTransition(new TransitionIndividual(patchName+".lG."+disease_state, patchName+".He."+disease_state, 2*mPatchParams[patchName].mX*mPatchParams[patchName].mN[2]));
        rChain.addTransition(new TransitionIndividual(patchName+".lG."+disease_state, patchName+".Rs."+disease_state, 2*(1-mPatchParams[patchName].mX)*mPatchParams[patchName].mN[2]));
        
        //Death
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".Ch."+disease_state, mPatchParams[patchName].mDelta[1]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".eG."+disease_state, mPatchParams[patchName].mDelta[2]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".lG."+disease_state, mPatchParams[patchName].mAlpha["lG"]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".lG."+disease_state, (1-mPatchParams[patchName].mAlpha["lG"])*mPatchParams[patchName].mDelta[3]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".He."+disease_state, mPatchParams[patchName].mAlpha["He"]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".He."+disease_state, (1-mPatchParams[patchName].mAlpha["He"])*mPatchParams[patchName].mDelta[3]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".Rs."+disease_state, mPatchParams[patchName].mAlpha["Rs"]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".Rs."+disease_state, (1-mPatchParams[patchName].mAlpha["Rs"])*mPatchParams[patchName].mDelta[4]));
        
      }

      //New egg laying rate from discussion on 10/08/18
      std::vector<std::string> governing_states = {patchName+".He.S", patchName+".He.E"};
      rChain.addTransition(new TransitionIndividualFromVoid(patchName+".E", mPatchParams[patchName].mNEgg * mPatchParams[patchName].mBr / 365.0, governing_states));
      
      std::vector<std::string> infected_states = {"Ch.I", "eG.I", "lG.I", "He.I", "Rs.I"};
      for (int i = 0 ; i < infected_states.size() ; i++)
      {
        infected_states[i] = patchName+"."+infected_states[i];
      }
      
      for (std::string demographic_state : demographic_states)
      {
        //Within patch:
        TransitionMassActionByPopulation* infection_transition = new TransitionMassActionByPopulation(patchName +"." + demographic_state+".S", patchName + "." + demographic_state+".E", mPatchParams[patchName].mBeta[patchName], within_patch_population_states, infected_states);
        rChain.addTransition(infection_transition);
        TransitionIndividual* incidence_transition = new TransitionIndividual(patchName + "." + demographic_state+".E", patchName + "." + demographic_state+".I", mPatchParams[patchName].mSigma);
        incidence_transition->addCounter(patchName+".infection");
        rChain.addTransition(incidence_transition);
        rChain.addTransition(new TransitionIndividualToVoid(patchName + "." + demographic_state+".I", mPatchParams[patchName].mGamma));
      }
    }

    //Between patches:
    for (std::string patchName : mPatchNames)
    {
      for (std::string demographic_state : demographic_states)
      {
        for (std::string other_patch : mPatchNames)
        {
          if (other_patch != patchName)
          {
            std::vector<std::string> infected_states = {"Ch.I", "eG.I", "lG.I", "He.I", "Rs.I"};
            for (std::string& s : infected_states)
              s = other_patch + "." + s;
            
            const std::vector<std::string> disease_states = {"S", "E", "I"};
            const std::vector<std::string> demographic_states = {"Ch","eG","lG","He","Rs"};
            std::vector<std::string> denominator_states;

            for (std::string demo_state : demographic_states)
            {
              for (std::string dis_state : disease_states)
              {
                denominator_states.push_back(other_patch+"."+demo_state+"."+dis_state);
              }
            }

            TransitionMassActionByPopulation* infection_transition = new TransitionMassActionByPopulation(patchName+"."+demographic_state+".S", patchName+"."+demographic_state+".E", mPatchParams[patchName].mBeta[other_patch], denominator_states, infected_states);
            rChain.addTransition(infection_transition);
          }
        }
      }
    }
  }
};