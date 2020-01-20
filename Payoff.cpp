#ifndef Payoff_cpp
#define Payoff_cpp


#include "Payoff.h"

PayOff::PayOff() {}

//payoff Call
PayOffCall::PayOffCall(const double& _K) :K(_K) {}

double PayOffCall::operator() (const double& S) const {
	return std::max(S - K, 0.0); // Standard European call pay-off
}
	
// Constructor with single strike parameter
PayOffPut::PayOffPut(const double& _K) {
	K = _K;
}

// Over-ridden operator() method, which turns PayOffPut into a function object
double PayOffPut::operator() (const double& S) const {
	return std::max(K - S, 0.0); // Standard European put pay-off
}


#endif






