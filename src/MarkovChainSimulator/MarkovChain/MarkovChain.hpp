#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H

#include <iostream>
#include <utility>
#include <vector>
#include <assert.h>
#include <cmath>
#include <random>
#include <numeric>
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/operators.hpp>
#include <functional>
#include "StateValues.h"
#include "Transitions.cpp"
#include "Serialiser.hpp"

namespace pl = std::placeholders;

namespace Deterministic {
class State : boost::additive1<State,
                boost::additive2<State, double, 
                    boost::multiplicative2<State, double> > >
{
public:
    using Map = std::map<std::string, double>;
    State(Map const& map);
    State() = default;
    State(const State &p) = default;
    State &operator=(State const&a) = default;

    void addToKey(std::string key, double a);
    State &operator+=(const State &p);
    State &operator+=(double a);
    State &operator*=(double f);
    friend State abs(const State &p);
    friend State operator/(const State &p1, const State &p2);
    friend double vector_space_norm_inf_impl(State const& p);
    size_t size() const;
    void resize(State const& other);
    Map getMap() const;

private:
    Map mMap;
};
}

class MarkovChain
{
private:
    unsigned long mix(unsigned long a, unsigned long b, unsigned long c);
    double T_MAX = 500;
    std::string filename;
    Serialiser *mpSerialiser;
    bool debug = false;

    unsigned long seed = mix(clock(), time(NULL), getpid());

    void solveGillespie();
    using DeterministicStateType = Deterministic::State;
    void derivative(const DeterministicStateType p, DeterministicStateType &dpdt, const double t);
    void serialiserDeterministic(const DeterministicStateType &p, const double t);
    void solveDeterministic();
    void solveRK4();
    void solveRKD5();
    void solveForwardEuler();

protected:
    state_values states;
    std::vector<Transition* > transitions;

public:
    void setDebug();
    void setSerialiser(Serialiser *serialiser);
    const static int SOLVER_TYPE_GILLESPIE = -1;
    void addState(std::string state_name, double initial_value);
    void addTransition(Transition *transition);
    void setMaxTime(double newMaxTime);
    void solve(int solver_type);
    void setSeed(double seed);
};
#endif