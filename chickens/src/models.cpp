#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <string>
#include "MarkovChainSimulator/MarkovChain/MarkovChain.hpp"
#include "MarkovChainSimulator/MarkovChain/MarkovChain.cpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.hpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.cpp"
#include "modeldefs.cpp"
using namespace Rcpp;

class SerialiserR : public Serialiser {
  
private:
  std::map<std::string, std::vector<double>> mResults;
  std::vector<double> mSerialiseTimes;
  state_values mLastState;
  double mLastT;
  bool shouldInterpolate = true;
public:
  
  SerialiserGASReducedScabies(std::vector<double> serialiseTimes) : mSerialiseTimes(serialiseTimes) {}
  
  void setShouldInterpolate(bool status)
  {
    shouldInterpolate = status;
  }
  
  virtual void serialise(double t, state_values states)
  {
    double next_time = *mSerialiseTimes.begin();
    while (t > next_time && !mSerialiseTimes.empty())
    {
      state_values interpolated_states;
      if (shouldInterpolate)
      {
        for (auto &p : states)
        {
          double slope = ( (double) p.second - mLastState[p.first]) / (t - mLastT);
          interpolated_states[p.first] = slope * (next_time - mLastT) + mLastState[p.first];
        }
      } else
      {
        for (auto &p : states)
        {
          interpolated_states[p.first] = mLastState[p.first];
        }
      }
      mResults["t"].push_back(next_time);
      for (auto& state : interpolated_states)
      {
        mResults[state.first].push_back(state.second);
      }
      mSerialiseTimes.erase(mSerialiseTimes.begin());
      next_time = *mSerialiseTimes.begin();
    }
    mLastState = states;
    mLastT = t;
    
  }
  
  virtual void serialiseHeader(state_values states)
  {
  }
  
  virtual void serialiseFinally(double t, state_values states)
  {
    serialise(t, states);
  }
  
  List getResults()
  {
    List list(mResults.size());
    CharacterVector namevec;
    int i = 0;
    for (auto& result : mResults)
    {
      namevec.push_back(result.first);
      list[i] = result.second;
      i++;
    }
    list.attr("names") = namevec;
    return list;
  }
  
};

std::map<std::string, double> convertListToMap(List list)
{
  std::map<std::string, double> ret;
  CharacterVector names = list.names();
  for (int i = 0 ; i < names.size() ; i++)
  {
    ret[as<std::string>(names[i])] = as<double>(list[i]);
  }
  
  return (ret);
}

// [[Rcpp::export(.chickens_model)]]
List chickens_model(List parameters_patch, NumericMatrix betas, double max_time, double dt, int solver_type) {
  //parameters_patch contains the within-patch parameters
  //betas is the mixing matrix, which is named.
  
  std::vector<std::string> patchNames = as<std::vector<std::string>>(rownames(betas));
  //Can access each set of within-patch parameters using patchNames now.
  
  std::map<std::string, WithinPatchParameters> param_map;
  
  int i = 0;
  for (std::string patchName : patchNames)
  {
    List sublist = parameters_patch[patchName];
    stringmap x0 = convertListToMap(sublist["x0"]);
    std::vector<double> n = as<std::vector<double>>(sublist["n"]);
    std::vector<double> delta = as<std::vector<double>>(sublist["delta"]);
    double rho = as<double>(sublist["rho"]);
    double y = as<double>(sublist["y"]);
    double x = as<double>(sublist["x"]);
    stringmap alpha = convertListToMap(sublist["alpha"]);
    stringmap beta;
    double j = 0;
    for (std::string p : patchNames)
    {
      beta[p] = betas(i, j);
      j++;
    }
    double sigma = as<double>(sublist["sigma"]);
    double gamma = as<double>(sublist["gamma"]);
    
    param_map[patchName] = WithinPatchParameters(x0, n, delta, rho, y, x, alpha, beta, sigma, gamma);
    i++;
  }
  
  std::vector<double> serialiser_times(max_time/dt + 1);
  double n = {-1 * dt};
  std::generate(serialiser_times.begin(), serialiser_times.end(), [&n, dt] { return n+=dt;});
  SerialiserR serialiser(serialiser_times)
    
  MarkovChain chain;
  chain.setSerialiser(serialiser);
  chain.setMaxTime(max_time);
  
  ModelChickenFlu model = ModelChickenFlu(patchNames, param_map);
  model.setupModel(chain);
  
  chain.solve(solver_type);
  
  return (serialiser.getResults());
}
