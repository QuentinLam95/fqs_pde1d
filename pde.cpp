#ifndef pde_cpp
#define pde_cpp

#include "pde.hpp"
#include <math.h>

BS_PDE::BS_PDE(VanillaOption* _option, const std::string& _left_boundary_type , const std::string& _right_boundary_type)
	: option(_option), left_boundary_type(_left_boundary_type), right_boundary_type(_right_boundary_type) {}

std::string BS_PDE::get_right_boundary_type() const {
	return right_boundary_type ;
}

std::string BS_PDE::get_left_boundary_type() const {
	return left_boundary_type ;
}

double BS_PDE::diff_coeff(double t, double x, double v) const {
	double vol = option->sigma;
	return 0.5 * vol * vol ;  
}

double BS_PDE::conv_coeff(double t, double x, double v) const {
	return (option->r) - diff_coeff(double t, double x, double v) ;
}

double BS_PDE::zero_coeff(double t, double x, double v) const {
	return (option->r) ;
}

//maybe delete this
double BS_PDE::source_coeff(double t, double x, double v) const {
	return 0.0;
}

double BS_PDE::boundary_left(double t, double x, double v) const {
	return 0.0;
}

double BS_PDE::boundary_right(double t, double x, double v) const {
	if (left_boundary_type.compare("D") == 0) {
		return (x - (option->K) * exp(-(option->r) * ((option->T) - t))); //careful to use exp(x) when calling this function
	}
	else if ((left_boundary_type.compare("N") == 0)) {
		return 1.0; //valeur de la dérivée à s infini
	}
}

double BS_PDE::init_cond(double x)  {
	return option->pay_off->operator()(x);
}

std::vector<double> init_cond(std::vector<double> X)
{
	std::vector<double> res(X.size());
	std::transform(X.begin(), X.end(), res.begin(), [this](double x)->{ return init_cond(x);});
	return res;
}

double BS_PDE::standard_dev() {
	double vol = option->sigma;
	double maturity = option->T;
	return vol*sqrt(maturity)
}

#endif // !1
