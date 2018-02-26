#ifndef SERIALISER_H
#define SERIALISER_H

#include "StateValues.h"
#include <iostream>
#include <fstream>
#include <vector>

class Serialiser
{
public:
    virtual void serialise(double t, state_values states);
    virtual void serialiseHeader(state_values states);
    virtual void serialiseFinally(double t, state_values states);
};

class SerialiserFile : public Serialiser
{
private:
    std::ofstream mOutputfile;

public:
    SerialiserFile(std::string filename);
    virtual void serialise(double t, state_values states);
    virtual void serialiseHeader(state_values states);
    virtual void serialiseFinally(double t, state_values states);
};


class SerialiserFileFinalState : public SerialiserFile
{
public:
    SerialiserFileFinalState(std::string filename);
    virtual void serialise(double t, state_values states);
    virtual void serialiseFinally(double t, state_values states);
};


class SerialiserPredefinedTimes : public Serialiser
{
private:
    std::vector<double> mSerialiseTimes;
    state_values mLastState;

public:
    SerialiserPredefinedTimes(std::vector<double> serialiseTimes);
    virtual void serialise(double t, state_values states);
};


class SerialiserPredefinedTimesFile : public SerialiserFile {
private:
    std::vector<double> mSerialiseTimes;
    state_values mLastState;
    double mLastT;

public:
    SerialiserPredefinedTimesFile(std::vector<double> serialiseTimes, std::string filename);
    virtual void serialise(double t, state_values states);
    virtual void serialiseFinally(double t, state_values states);
};

#endif