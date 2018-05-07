#include <math.h>
#include <map>
#include <vector>
#include <Rcpp.h>
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
  
  WithinPatchParameters(stringmap x0, std::vector<double> n, std::vector<double> delta, double y, double x, stringmap alpha, stringmap beta, double sigma, double gamma, double nEgg, double q, double w, double K) :
    mX0(x0), mN(n), mDelta(delta), mY(y), mX(x), mAlpha(alpha), mBeta(beta), mSigma(sigma), mGamma(gamma), mNEgg(nEgg), mQ(q), mW(w), mK(K) {}
  
  WithinPatchParameters() {}
};

class ModelChickenFlu {
  
private:   
  static double eggLayingRate(state_values states, parameter_map parameters) {
    double population_size = 0;
    //Get patch name:
    std::string patchname;
    if (parameters["Sc"] == 1)
      patchname = "Sc";
    if (parameters["Ns"] == 1)
      patchname = "Ns";
    if (parameters["Es"] == 1)
      patchname = "Es";
    if (parameters["Bs"] == 1)
      patchname = "Bs";
    
    for (state_values::iterator it = states.begin() ; it != states.end() ; it++) 
    {
      std::string s = it->first;
      std::string delimeter = ".";
      std::string current_patch = s.substr(0, s.find(delimeter));
      s.erase(0, s.find(delimeter) + delimeter.length());
      std::string state = s.substr(0, s.find(delimeter));
      if (state != "E" && current_patch == patchname) 
      {
        population_size = population_size + it->second;
      }
    }
    if (population_size == 0)
    {
      return (0);
    }
    
    double f = ((double) states[patchname+".He.S"])/population_size;
    double mu = (f * parameters["n_egg"] * parameters["q"]); // div 365 for years -> days
    double b = log(mu + 1);
    
    if (patchname != "Sc")
      return (b);
    
    double mu_bar = b * population_size * (1-( ((double) population_size)/parameters["K"]));
    if (mu_bar < 0)
      mu_bar = 0;
    if (parameters["ES"]>0) {
      return ( parameters["w"] * mu_bar);
    }
    else {
      return (mu_bar);
    }
  }
  
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
      if (state != "E" && current_patch == patchname) 
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
      if (alpha < 0)
        alpha = 0;
    }
    
    if (parameters["returnrho"])
    {
      return (rho); 
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
      
      rChain.addTransition(new TransitionCustom(patchName+".E", patchName+".Ch.S", hatchingParameters, *eggHatchingRate));
      
      hatchingParameters["returnsold"]=1; hatchingParameters["returnhatching"]=0;
      rChain.addTransition(new TransitionCustomToVoid(patchName+".E", hatchingParameters, *eggHatchingRate));
      
      hatchingParameters["returnrho"]=1; hatchingParameters["returnsold"]=0;
      rChain.addTransition(new TransitionCustomFromVoid(patchName+".Ch.S", hatchingParameters, *eggHatchingRate));
      
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
      
      parameter_map eggParameters;
      eggParameters["w"] = mPatchParams[patchName].mW;
      eggParameters["K"] = mPatchParams[patchName].mK;
      eggParameters["n_egg"] = mPatchParams[patchName].mNEgg;
      eggParameters["q"] = mPatchParams[patchName].mQ;
      eggParameters[patchName] = 1;
      
      rChain.addTransition(new TransitionCustomFromVoid(patchName+".E", eggParameters, *eggLayingRate));
      
      std::vector<std::string> infected_states = {"Ch.I", "eG.I", "lG.I", "He.I", "Rs.I"};
      for (int i = 0 ; i < infected_states.size() ; i++)
      {
        infected_states[i] = patchName+"."+infected_states[i];
      }
      
      for (std::string demographic_state : demographic_states)
      {
        //Within patch:
        rChain.addTransition(new TransitionMassActionByPopulation(patchName +"." + demographic_state+".S", patchName + "." + demographic_state+".E", mPatchParams[patchName].mBeta[patchName], within_patch_population_states, infected_states));
        //Betwen patches:
        for (std::string other_patch : mPatchNames)
        {
          if (other_patch != patchName)
          {
            rChain.addTransition(new TransitionMassActionByPopulation(patchName+"."+demographic_state+".S", patchName+"."+demographic_state+".E", mPatchParams[patchName].mBeta[other_patch], population_states, population_infectious));
          }
        }
        rChain.addTransition(new TransitionIndividual(patchName + "." + demographic_state+".E", patchName + "." + demographic_state+".I", mPatchParams[patchName].mSigma));
        rChain.addTransition(new TransitionIndividualToVoid(patchName + "." + demographic_state+".I", mPatchParams[patchName].mGamma));
      }
    }
  }
};