#include <iostream>
#include <utility>
#include "StateValues.h"

class Transition {

protected:
  std::string mSource_state;
  std::string mDestination_state;
  std::vector<std::string> mGoverning_states;
  parameter_map mParameters;
  std::vector<std::string> mCounters;

  double (*mpGetActualRate)(state_values pStates, parameter_map parameters);

  int mTransition_type = 0;

  virtual void incrementCounters(state_values &rStates)
  {
    if (mCounters.size() > 0)
    {
      for (std::string counter : this->mCounters)
      {
        rStates[counter] += 1;
      }
    }
  }

public:
  
  Transition() {};
  
  Transition(std::string source_state, std::string destination_state, parameter_map parameters, double (*getActualRate)(state_values pStates, parameter_map parameters), std::vector<std::string> governing_states = {}) :
    mSource_state(source_state), mDestination_state(destination_state), mParameters(parameters), mpGetActualRate(getActualRate), mGoverning_states(governing_states) 
    {
      if (mGoverning_states.size() == 0) {
        mGoverning_states = {mDestination_state};
      }
    }

    Transition(std::string source_state, std::string destination_state, double parameter, std::vector<std::string> governing_states = {}) :
    mSource_state(source_state), mDestination_state(destination_state), mGoverning_states(governing_states)
    {
      if (mGoverning_states.size() == 0) {
        mGoverning_states = {mDestination_state};
      }

      mParameters["parameter"] = parameter;
    }

    
  virtual void setStates(std::string source_state, std::string destination_state) {
    mSource_state = source_state;
    mDestination_state = destination_state;
  }

  virtual void do_transition (double t, state_values &rStates) = 0;

  virtual double getRate(state_values states) = 0;

  virtual std::string getSourceState() const 
  {
    return (mSource_state);
  }

  virtual void setSourceState(std::string source_state) 
  {
    mSource_state = source_state;
  }

  virtual std::string getDestinationState() const 
  {
    return (mDestination_state);
  }

  virtual void setDestinationState(std::string destination_state) 
  {
    mDestination_state = destination_state;
  }

  double getSingleParameter() const {
    return (mParameters.begin()->second);
  }

  void setSingleParameter(double parameter) {
    mParameters["parameter"] = parameter;
  }

  std::vector<std::string> getGoverningStates() const {
    return (mGoverning_states); 
  }

  void setGoverningStates(std::vector<std::string> governing_states) {
    mGoverning_states = governing_states;
  }

  virtual void addCounter(std::string counter)
  {
    mCounters.push_back(counter);
  }

  virtual std::vector<std::string> getCounters()
  {
    return (mCounters);
  }

  virtual int getNumCounters()
  {
    return (mCounters.size());
  }
  

};


class TransitionIndividual : public Transition 
{
public:
  TransitionIndividual(std::string source_state, std::string destination_state, double parameter)
    : Transition(source_state, destination_state, parameter, {})
    {}

  virtual double getRate(state_values states)
  {
    //std::cout << this->mSource_state << std::endl;
    return (this->mParameters["parameter"] * states[this->mSource_state]);
  }

  virtual void do_transition(double t, state_values &rStates) {
    rStates[this->mSource_state] -= 1;
    rStates[this->mDestination_state] += 1;

    incrementCounters(rStates);
  }
};


class TransitionMassAction : public Transition
{
public:
  TransitionMassAction(std::string source_state, std::string destination_state, double parameter, std::vector<std::string> governing_states = {})
    : Transition(source_state, destination_state, parameter, governing_states)
    {}

  virtual double getRate(state_values states)
  {
    double mass = 0;
    for (std::vector<std::string>::iterator it = this->mGoverning_states.begin() ; it != this->mGoverning_states.end() ; it++)
    {
      mass += states[*it];
    }
    return (this->mParameters["parameter"] * states[this->mSource_state] * mass);
  }

  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mSource_state] -= 1;
    rStates[this->mDestination_state] += 1;

    incrementCounters(rStates);
  }
};


class TransitionIndividualToVoid : public TransitionIndividual
{
public:
  TransitionIndividualToVoid(std::string source_state, double parameter)
    : TransitionIndividual(source_state, "Void", parameter)
    {}
  
  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mSource_state] -= 1;

    incrementCounters(rStates);
  }
};

class TransitionIndividualFromVoid : public Transition
{
public:
  TransitionIndividualFromVoid(std::string destination_state, double parameter, std::vector<std::string> governing_states)
    : Transition("Void", destination_state, parameter, governing_states)
    {}
  
  virtual double getRate(state_values states)
  {
    double mass = 0;
    for (std::vector<std::string>::iterator it = this->mGoverning_states.begin() ; it != this->mGoverning_states.end() ; it++)
    {
      mass += states[*it];
    }
    return (this->mParameters["parameter"] * mass);
  }

  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mDestination_state] += 1;

    incrementCounters(rStates);
  }
};

class TransitionCustom : public Transition
{
public:
  TransitionCustom(std::string source_state, std::string destination_state, parameter_map parameters, double (*getActualRate)(state_values pStates, parameter_map parameters))
    : Transition(source_state, destination_state, parameters, getActualRate)
    {}

  virtual double getRate(state_values states)
  {
    return (this->mpGetActualRate(states, this->mParameters));
  }

  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mSource_state] -= 1;
    rStates[this->mDestination_state] += 1;

    incrementCounters(rStates);
  }
};

class TransitionCustomFromVoid : public TransitionCustom
{
public:
  TransitionCustomFromVoid(std::string destination_state, parameter_map parameters, double (*getActualRate)(state_values pStates, parameter_map parameters))
    : TransitionCustom("Void", destination_state, parameters, getActualRate)
    {}

  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mDestination_state] += 1;

    incrementCounters(rStates);
  }
};

class TransitionCustomToVoid : public TransitionCustom
{
public:
  TransitionCustomToVoid(std::string source_state, parameter_map parameters, double (*getActualRate)(state_values pStates, parameter_map parameters))
    : TransitionCustom(source_state, "Void", parameters, getActualRate)
    {}

  virtual void do_transition(double t, state_values &rStates)
  {
    rStates[this->mSource_state] -= 1;

    incrementCounters(rStates);
  }
};


class TransitionMassActionByPopulation : public TransitionMassAction
{
private:
  std::vector<std::string> mPopulationStates;
public:
  TransitionMassActionByPopulation(std::string source_state, std::string destination_state, double parameter, std::vector<std::string> population_states, std::vector<std::string> governing_states = {})
    : TransitionMassAction(source_state, destination_state, parameter, governing_states), mPopulationStates(population_states)
    {}

  virtual double getRate(state_values states)
  {
    double mass = 0;
    for (std::vector<std::string>::iterator it = this->mGoverning_states.begin() ; it != this->mGoverning_states.end() ; it++)
    {
      mass += states[*it];
    }

    double population_size = 0;
    for (std::vector<std::string>::iterator it = mPopulationStates.begin() ; it != mPopulationStates.end() ; it++)
    {
      population_size += states[*it];
    }
    if (population_size == 0)
    {
      return (0);
    }

    return ( (this->mParameters["parameter"] * states[this->mSource_state] * mass)/population_size );
  }
};


class TransitionConstant : public TransitionIndividual
{
public:
  TransitionConstant(std::string source_state, std::string destination_state, double parameter) : TransitionIndividual(source_state, destination_state, parameter) {}

  virtual double getRate(state_values states)
  {
    return (this->mParameters["parameter"]*((double) states[this->mSource_state] > 0));
  }
};