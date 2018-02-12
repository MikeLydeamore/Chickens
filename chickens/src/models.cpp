#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <string>
#include "MarkovChainSimulator/MarkovChain/MarkovChain.hpp"
#include "MarkovChainSimulator/MarkovChain/MarkovChain.cpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.hpp"
#include "MarkovChainSimulator/MarkovChain/Serialiser.cpp"
#include "modeldefs.cpp"
using namespace Rcpp;

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
List chickens_model(List parameters_patch, NumericMatrix betas) {
  //parameters_patch contains the within-patch parameters
  //betas is the mixing matrix, which is named.
  
  std::vector<std::string> patchNames = as<std::vector<std::string>>(rownames(betas));
  //Can access each set of within-patch parameters using patchNames now.
  
  std::map<std::string, WithinPatchParameters> param_map;
  
  CharacterVector names = parameters_patch.names();
  for (int i = 0 ; i < names.size() ; i++)
  {
    std::cout << names[i] << std::endl;
  }
  int i = 0;
  for (std::string patchName : patchNames)
  {
    List sublist = parameters_patch[patchName];
    stringmap x0 = convertListToMap(sublist["x0"]);
    std::vector<double> n = as<std::vector<double>>(sublist["n"]);
    std::vector<double> delta = as<std::vector<double>>(sublist["delta"]);
    double rho = as<double> (sublist["rho"]);
    double y = as<double> (sublist["y"]);
    double x = as<double> (sublist["x"]);
    stringmap alpha = convertListToMap(sublist["alpha"]);
    stringmap beta;
    double j = 0;
    for (std::string p : patchNames)
    {
      beta[p] = betas(i, j);
      j++;
    }
    double sigma = as<double> (sublist["sigma"]);
    double gamma = as<double> (sublist["gamma"]);
    
    param_map[patchName] = WithinPatchParameters(x0, n, delta, rho, y, x, alpha, beta, sigma, gamma);
    i++;
  }
  
  return (List::create(Named("vec") = 0));
}
