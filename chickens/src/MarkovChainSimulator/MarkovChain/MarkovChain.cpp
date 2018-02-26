#include "MarkovChain.hpp"

namespace Deterministic
{
State::State(Map const &map) : mMap(map) {}

void State::addToKey(std::string key, double a)
{
    mMap[key] += a;
}

State &State::operator+=(const State &p)
{
    for (auto &p : p.mMap)
        mMap[p.first] += p.second;
    return *this;
}

State &State::operator+=(double a)
{
    for (auto &p : mMap)
        p.second += a;
    return *this;
}

State &State::operator*=(double f)
{
    for (auto &p : mMap)
        mMap[p.first] *= f;
    return *this;
}

State abs(const State &p)
{
    using std::abs;
    auto map = p.mMap;

    for (auto &e : map)
        e.second = abs(e.second);

    return map;
}

State operator/(const State &p1, const State &p2)
{
    auto map = p1.mMap;

    for (auto &e : map)
        e.second /= p2.mMap.at(e.first);

    return map;
}

double vector_space_norm_inf_impl(State const &p)
{
    double max = 0;
    using std::abs;
    for (auto &el : p.mMap)
        max = std::max(abs(el.second), max);
    return max;
}

size_t State::size() const { return mMap.size(); }

void State::resize(State const &other)
{
    for (auto &el : other.mMap)
        mMap[el.first] += 0; // inserts if non-existent
}

using Map = std::map<std::string, double>;
Map State::getMap() const
{
    return (mMap);
}
}

using DeterministicStateType = Deterministic::State;
namespace boost
{
namespace numeric
{
namespace odeint
{
template <>
struct vector_space_norm_inf<DeterministicStateType>
{
    typedef double result_type;
    double operator()(const DeterministicStateType &p) const { return vector_space_norm_inf_impl(p); }
};

template <>
struct is_resizeable<DeterministicStateType>
{
    typedef boost::true_type type;
    const static bool value = type::value;
};

template <>
struct same_size_impl<DeterministicStateType, DeterministicStateType>
{
    static bool same_size(const DeterministicStateType &v1, const DeterministicStateType &v2)
    {
        return v1.size() == v2.size();
    }
};

template <>
struct resize_impl<DeterministicStateType, DeterministicStateType>
{
    static void resize(DeterministicStateType &v1, const DeterministicStateType &v2)
    {
        v1.resize(v2);
    }
};
}
}
}

using namespace boost::numeric::odeint;

unsigned long MarkovChain::mix(unsigned long a, unsigned long b, unsigned long c)
{
    a = a - b;
    a = a - c;
    a = a ^ (c >> 13);
    b = b - c;
    b = b - a;
    b = b ^ (a << 8);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 13);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 12);
    b = b - c;
    b = b - a;
    b = b ^ (a << 16);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 5);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 3);
    b = b - c;
    b = b - a;
    b = b ^ (a << 10);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 15);
    return c;
}

void MarkovChain::solveGillespie()
{
    double t = 0;
    bool ended_infinite = false;

    int population_size = 0;
    for (typename state_values::iterator it = states.begin(); it != states.end(); it++)
    {
        assert(it->second >= 0);
        population_size += it->second;
    }

    typedef boost::uniform_real<> NumberDistribution;
    typedef boost::mt19937 RandomNumberGenerator;
    typedef boost::variate_generator<RandomNumberGenerator &,
                                     NumberDistribution>
        Generator;

    NumberDistribution distribution(0, 1);
    RandomNumberGenerator generator;
    Generator runif(generator, distribution);
    
    generator.seed(seed); // seed with the current time

    mpSerialiser->serialiseHeader(states);

    mpSerialiser->serialise(t, states);

    while (t < T_MAX && !ended_infinite)
    {

        int actual_size = 0;
        for (typename state_values::iterator it = states.begin(); it != states.end(); it++)
        {
            assert(it->second >= 0);
            actual_size += it->second;
        }
        //assert(actual_size == population_size);

        std::vector<double> rates(transitions.size());
        for (int i = 0; i < transitions.size(); i++)
        {
            rates[i] = transitions[i]->getRate(states);
            if (rates[i] < 0.0 || isnan(rates[i]))
            {
                std::cout << "Transition from " << transitions[i]->getSourceState() << " to " << transitions[i]->getDestinationState() << " has rate " << rates[i] << std::endl;
                exit(-1);
            }
            assert(rates[i] >= 0);
        }
        double rates_sum = std::accumulate(rates.begin(), rates.end(), (double)0.0);

        double event_time = -(1.0 / rates_sum) * log(runif());
        if (isinf(event_time))
        {
            ended_infinite = true;
            break;
        }
        t += event_time;

        std::vector<double> rates_normalised(rates.size());
        rates_normalised[0] = rates[0] / rates_sum;

        for (int i = 1; i < rates.size(); i++)
        {
            rates_normalised[i] = rates_normalised[i - 1] + (rates[i] / rates_sum);
        }

        int eventOccurred = 0;
        double u = runif();

        while (u > rates_normalised[eventOccurred])
        {
            eventOccurred++;
        }
        transitions[eventOccurred]->do_transition(t, states);
        mpSerialiser->serialise(t, states);
    }
    mpSerialiser->serialiseFinally(t, states);
}

void MarkovChain::derivative(const DeterministicStateType p, DeterministicStateType &dpdt, const double t)
{
    for (int i = 0; i < transitions.size(); i++)
    {
        dpdt.addToKey(transitions[i]->getSourceState(), -1 * transitions[i]->getRate(p.getMap()));
        dpdt.addToKey(transitions[i]->getDestinationState(), transitions[i]->getRate(p.getMap()));
    }
}

void MarkovChain::serialiserDeterministic(const DeterministicStateType &p, const double t)
{
    mpSerialiser->serialise(t, p.getMap());
}

void MarkovChain::solveDeterministic()
{
    DeterministicStateType x0(states);
    typedef runge_kutta_cash_karp54<DeterministicStateType, double, DeterministicStateType, double, vector_space_algebra>
        rkck54;

    //typedef controlled_runge_kutta< rkck54, double, DeterministicStateType, double, vector_space_algebra > ctrl_rkck54;
    mpSerialiser->serialiseHeader(x0.getMap());
    integrate_adaptive(make_controlled(1e-10, 1e-6, rkck54()), std::bind(&MarkovChain::derivative, *this, pl::_1, pl::_2, pl::_3), x0, 0.0,
                       T_MAX, 0.1,
                       std::bind(&MarkovChain::serialiserDeterministic, *this, pl::_1, pl::_2));

    mpSerialiser->serialiseFinally(T_MAX, x0.getMap());
}

void MarkovChain::solveRK4()
{
    DeterministicStateType y(states);
    mpSerialiser->serialiseHeader(y.getMap());
    double t = 0;
    double h = 1.0 / 120;
    if (T_MAX < 5)
        h = 1.0 / (5000 * T_MAX);
    while (t < T_MAX)
    {
        mpSerialiser->serialise(t, y.getMap());
        DeterministicStateType k_1;
        DeterministicStateType k_2;
        DeterministicStateType k_3;
        DeterministicStateType k_4;
        derivative(y, k_1, t);
        derivative(y + (h / 2) * k_1, k_2, t + (h / 2));
        derivative(y + (h / 2) * k_2, k_3, t + (h / 2));
        derivative(y + h * k_3, k_4, t + h);

        y += (h / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
        t += h;
    }

    mpSerialiser->serialiseFinally(t, y.getMap());
}

void MarkovChain::solveRKD5()
{
    DeterministicStateType y(states);
    mpSerialiser->serialiseHeader(y.getMap());
    double t = 0;
    double tol = 0.00001;

    double a21 = 1.0 / 5.0;
    double a31 = 3.0 / 40.0;
    double a32 = 9.0 / 40.0;
    double a41 = 44.0 / 45.0;
    double a42 = -56.0 / 15.0;
    double a43 = 32.0 / 9.0;
    double a51 = 19372.0 / 6561.0;
    double a52 = -25360.0 / 2187.0;
    double a53 = 64448.0 / 6561.0;
    double a54 = -212.0 / 729.0;
    double a61 = 9017.0 / 3168.0;
    double a62 = -355.0 / 33.0;
    double a63 = 46732.0 / 5247.0;
    double a64 = 49.0 / 176.0;
    double a65 = -5103.0 / 18656.0;
    double a71 = 35.0 / 384.0;
    double a72 = 0.0;
    double a73 = 500.0 / 1113.0;
    double a74 = 125.0 / 192.0;
    double a75 = -2187.0 / 6784.0;
    double a76 = 11.0 / 84.0;

    double c2 = 1.0 / 5.0;
    double c3 = 3.0 / 10.0;
    double c4 = 4.0 / 5.0;
    double c5 = 8.0 / 9.0;
    double c6 = 1;
    double c7 = 1;

    double b1 = 35.0 / 384.0;
    double b2 = 0;
    double b3 = 500.0 / 1113.0;
    double b4 = 125.0 / 192.0;
    double b5 = -2187.0 / 6784.0;
    double b6 = 11.0 / 84.0;
    double b7 = 0;

    double b1p = 5179.0 / 57600.0;
    double b2p = 0;
    double b3p = 7571.0 / 16695.0;
    double b4p = 393.0 / 640.0;
    double b5p = -92097.0 / 339200.0;
    double b6p = 187.0 / 2100.0;
    double b7p = 1.0 / 40.0;

    double h = 1;

    while (t < T_MAX)
    {
        DeterministicStateType k_1;
        DeterministicStateType k_2;
        DeterministicStateType k_3;
        DeterministicStateType k_4;
        DeterministicStateType k_5;
        DeterministicStateType k_6;
        DeterministicStateType k_7;

        derivative(y, k_1, t);
        derivative(y + h * (a21 * k_1), k_2, t + c2 * h);
        derivative(y + h * (a31 * k_1 + a32 * k_2), k_3, t + c3 * h);
        derivative(y + h * (a41 * k_1 + a42 * k_2 + a43 * k_3), k_4, t + c4 * h);
        derivative(y + h * (a51 * k_1 + a52 * k_2 + a53 * k_3 + a54 * k_4), k_5, t + c5 * h);
        derivative(y + h * (a61 * k_1 + a62 * k_2 + a63 * k_3 + a64 * k_4 + a65 * k_5), k_6, t + c6 * h);
        derivative(y + h * (a71 * k_1 + a72 * k_2 + a73 * k_3 + a74 * k_4 + a75 * k_5 + a76 * k_6), k_7, t + c7 * h);

        DeterministicStateType yn = y + h * (b1 * k_1 + b3 * k_3 + b4 * k_4 + b5 * k_5 + b6 * k_6);
        DeterministicStateType ynstar = y + h * (b1p * k_1 + b3p * k_3 + b4p * k_4 + b5p * k_5 + b6p * k_6 + b7p * k_7);
        DeterministicStateType diff = yn + (-1 * ynstar);
        double error = 0;
        for (auto &p : diff.getMap())
        {
            error += (p.second * p.second);
        }
        error = pow(error, 0.5);

        double delta = 0.84 * pow(tol / error, (1.0 / 5.0));

        if (error < tol)
        {
            mpSerialiser->serialise(t, y.getMap());
            t += h;
            y += h * (b1 * k_1 + b3 * k_3 + b4 * k_4 + b5 * k_5 + b6 * k_6);
        }

        if (delta <= 0.1)
        {
            h *= 0.1;
        }
        else if (delta >= 4.0)
        {
            h *= 4.0;
        }
        else
        {
            h *= delta;
        }
    }

    mpSerialiser->serialiseFinally(t, y.getMap());
}

void MarkovChain::solveForwardEuler()
{
    DeterministicStateType y0(states);

    double t = 0;
    double h = 1.0 / 120.0;
    mpSerialiser->serialiseHeader(y0.getMap());
    while (t < T_MAX)
    {
        mpSerialiser->serialise(t, y0.getMap());
        DeterministicStateType dpdt;
        derivative(y0, dpdt, t);
        //std::cout << y0.getMap()["S"] << std::endl;
        y0 += h * dpdt;

        t += h;
    }

    mpSerialiser->serialiseFinally(t, y0.getMap());
}

void MarkovChain::setDebug()
{
    debug = true;
}

void MarkovChain::setSerialiser(Serialiser *serialiser)
{
    mpSerialiser = serialiser;
}

void MarkovChain::addState(std::string state_name, double initial_value)
{
    states[state_name] = initial_value;
}

void MarkovChain::addTransition(Transition *transition)
{
    if (debug)
    {
        double two_decimals = round(100 / transition->getSingleParameter()) / 100;
        std::cout << ("Adding transition from " + transition->getSourceState() + " to " + transition->getDestinationState() + " at rate " + std::to_string(transition->getSingleParameter()) + " (1/") << std::setprecision(2) << std::fixed << two_decimals << ")" << std::endl;
    }

    transitions.push_back(transition);
}

void MarkovChain::setMaxTime(double newMaxTime)
{
    T_MAX = newMaxTime;
}

void MarkovChain::solve(int solver_type)
{
    if (solver_type == SOLVER_TYPE_GILLESPIE)
    {
        solveGillespie();
    }
    else
    {
        solveRKD5();
    }
}

void MarkovChain::setSeed(double newSeed)
{
    seed = newSeed;
}