#include <math.h>
#include <map>
#include <vector>
#include "MarkovChainSimulator/MarkovChain/MarkovChain.hpp"

typedef std::map<std::string, double> stringmap;

class WithinPatchParameters
{
public:
  stringmap mX0;
  std::vector<double> mN;
  std::vector<double> mDelta;
  double mRho;
  double mY;
  double mX;
  stringmap mAlpha;
  stringmap mBeta;
  double mSigma;
  double mGamma;
  
  WithinPatchParameters(stringmap x0, std::vector<double> n, std::vector<double> delta, double rho, double y, double x, stringmap alpha, stringmap beta, double sigma, double gamma) :
    mX0(x0), mN(n), mDelta(delta), mRho(rho), mY(y), mX(x), mAlpha(alpha), mBeta(beta), mSigma(sigma), mGamma(gamma) {}
  
  WithinPatchParameters() {}
};

class ModelChickenFlu {
  
private:   
  static double eggLayingRate(state_values states, parameter_map parameters) {
    double population_size = 0;
    for (typename state_values::iterator it = states.begin() ; it != states.end() ; it++) {
      if (it->first != "E" || it->first != "ES") {
        population_size = population_size + it->second;
      }
    }
    if (population_size == 0)
    {
      return (0);
    }
    double n_egg = 8.9;
    double q = 4;
    double f = ((double) states["He.S"])/population_size;
    double mu = (f * population_size * n_egg * q) / population_size; // div 365 for years -> days
    double b = log(mu + 1);
    
    double mu_bar = b * population_size * (1-( ((double) population_size)/parameters["K"]));
    if (parameters["ES"]>0) {
      return ( parameters["w"] * mu_bar);
    }
    else {
      return (mu_bar);
    }
    
    //if nonscavenging, return b*N;
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
      rChain.addState(patchName+".E", mPatchParams[patchName].mX0[patchName+".E"]);
      
      
      std::vector<std::string> within_patch_population_states;
      
      for (std::string demographic_state : demographic_states) 
      {
        for (std::string disease_state : disease_states) 
        {
          std::string state_name = patchName + "." + demographic_state + "." + disease_state;
          
          rChain.addState(state_name, mPatchParams[patchName].mX0[state_name]);
          
          population_states.push_back(state_name);
          within_patch_population_states.push_back(state_name);
        }
      }
      
      //Aging transitions
      rChain.addTransition(new TransitionIndividual(patchName+".E", patchName+".Ch.S", mPatchParams[patchName].mN[0]));
      rChain.addTransition(new TransitionIndividualToVoid("E", mPatchParams[patchName].mDelta[0]));
      
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
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".lG."+disease_state, (1-mPatchParams[patchName].mAlpha["lG"]*mPatchParams[patchName].mDelta[3])));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".He."+disease_state, mPatchParams[patchName].mAlpha["He"]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".He."+disease_state, (1-mPatchParams[patchName].mAlpha["He"]*mPatchParams[patchName].mDelta[3])));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".Rs."+disease_state, mPatchParams[patchName].mAlpha["Rs"]));
        rChain.addTransition(new TransitionIndividualToVoid(patchName+".Rs."+disease_state, (1-mPatchParams[patchName].mAlpha["Rs"])*mPatchParams[patchName].mDelta[4]));
        
        //Import rates
      }
      
      parameter_map eggParameters;
      eggParameters["w"] = 1.0/32.1;
      eggParameters["K"] = 1100;
      rChain.addTransition(new TransitionCustomFromVoid("E", eggParameters, *eggLayingRate));
      
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