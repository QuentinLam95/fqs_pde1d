#pragma once

#ifndef Payoff_hpp
#define Payoff_hpp

#include <algorithm>


class PayOff {
public:
	PayOff(); //Constructeur
	virtual ~PayOff() {}; //destructeur

	// functor pour le payoff
	virtual double operator() (const double& S) const = 0;
};

class PayOffCall : public PayOff {
private:
	double K; // Strike 

public:
	PayOffCall(const double& K_);
	virtual ~PayOffCall() {};

	virtual double operator() (const double& S) const; 
};

class PayOffPut : public PayOff {
private:
	double K; // Strike

public:
	PayOffPut(const double& K_);
	virtual ~PayOffPut() {};
	virtual double operator() (const double& S) const;
};

#endif



