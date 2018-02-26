#include "Serialiser.hpp"
#include "StateValues.h"

void Serialiser::serialise(double t, state_values states) 
{
    std::cout << t;
    for (typename state_values::iterator it = states.begin(); it != states.end(); it++)
    {
        std::cout << "," << it->second;
    }
    std::cout << "\n";
}


void Serialiser::serialiseHeader(state_values states) 
{
    std::cout << "t";
    for (typename state_values::iterator it = states.begin(); it != states.end(); it++)
    {
        std::cout << "," << it->first;
    }
    std::cout << "\n";
}

void Serialiser::serialiseFinally(double t, state_values states) {}

SerialiserFile::SerialiserFile(std::string filename) 
{
    mOutputfile.open(filename);
}


void SerialiserFile::serialise(double t, state_values states) 
{   
    mOutputfile << t;

    typename std::map<std::string, double>::iterator it;

    for(it=states.begin(); it != states.end(); it++) {
        mOutputfile << "," << it->second;
    }
    mOutputfile << "\n";
    mOutputfile.flush();
}


void SerialiserFile::serialiseHeader(state_values states) 
{
    mOutputfile << "t";
    for (typename state_values::iterator it = states.begin(); it != states.end(); it++)
    {
        mOutputfile << "," << it->first;
    }
    mOutputfile << "\n";
    mOutputfile.flush();
}


void SerialiserFile::serialiseFinally(double t, state_values states) {}



SerialiserFileFinalState::SerialiserFileFinalState(std::string filename) : SerialiserFile(filename) {}

    
void SerialiserFileFinalState::serialise(double t, state_values states) {}

void SerialiserFileFinalState::serialiseFinally(double t, state_values states) 
{
    SerialiserFile::serialise(t, states);
}


SerialiserPredefinedTimes::SerialiserPredefinedTimes(std::vector<double> serialiseTimes) : mSerialiseTimes(serialiseTimes) {}
    

void SerialiserPredefinedTimes::serialise(double t, state_values states) 
{
    double next_time = *mSerialiseTimes.begin();

    while (t > next_time) 
    {
        Serialiser::serialise(next_time, mLastState);
        mSerialiseTimes.erase(mSerialiseTimes.begin());
        next_time = *mSerialiseTimes.begin();
    }
    mLastState = states;
}


SerialiserPredefinedTimesFile::SerialiserPredefinedTimesFile(std::vector<double> serialiseTimes, std::string filename) : SerialiserFile(filename), mSerialiseTimes(serialiseTimes) {}

 
void SerialiserPredefinedTimesFile::serialise(double t, state_values states) {
    double next_time = *mSerialiseTimes.begin();
    while (t > next_time && !mSerialiseTimes.empty()) {
        //Interpolate between the points.
        state_values interpolated_states;
        for (auto &p : states)
        {
            double slope = ( (double) p.second - mLastState[p.first]) / (t - mLastT);
            interpolated_states[p.first] = slope * (next_time - mLastT) + mLastState[p.first];
        }
        SerialiserFile::serialise(next_time, interpolated_states);
        mSerialiseTimes.erase(mSerialiseTimes.begin());
        next_time = *mSerialiseTimes.begin();
    }
    mLastState = states;
    mLastT = t;
}


void SerialiserPredefinedTimesFile::serialiseFinally(double t, state_values states)
{
    serialise(t, states);
}