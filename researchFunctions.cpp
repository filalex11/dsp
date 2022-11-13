#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> u1;
std::vector<double> u2;
const double f1 = 44100;
const double f2 = 1;
const double f3 = 2 * f1;
const double f4 = 2 * f2;
const double T = 1;
const double a = -T / 2;
const double b = T / 2;
const double dt = 0.001;
const double period = b - a;
const double pi = 3.1415926535898;

bool check_constant(const double number, const double reference) {
		return std::abs(number - reference) < 1e-6;
}

double find_scalar(const std::vector<double> &u1, const std::vector<double> &u2) {
	double result = 0;
	for (unsigned long int i = 0; i < u1.size(); ++i)
		result += u1[i] * u2[i];
	return (1 / period) * result;
}

double find_norm(const std::vector<double> &u) {
	return sqrt((1 / period) * find_scalar(u, u) * dt); 
}

bool is_orth(const std::vector<double> &u1, const std::vector<double> &u2) {
	return check_constant(find_scalar(u1, u2), 0);
}

bool is_basis(const std::vector<double> &u1, const std::vector<double> &u2) {
	return check_constant(find_scalar(u1, u2), 0) && check_constant(find_norm(u1), 1) && check_constant(find_norm(u2), 1);
}

void do_basis(std::vector<double> &u) {
	double norm = find_norm(u);
	for (unsigned long int i = 0; i < u.size(); ++i)
		u[i] = u[i] / norm;
}

int main() {
	for (double t = a; t <= b; t += dt) {
		u1.push_back(sin(2 * pi * f1 * t));
		u2.push_back(sin(2 * pi * f2 * t));
	}
	std::cout << "scalar " << find_scalar(u1, u2) << std::endl;
	std::cout << "norms " << find_norm(u1) << " " << find_norm(u2) << std::endl;
	std::cout << "is orth " << is_orth(u1, u2) << std::endl;
	std::cout << "is basis (before division by norm) " << is_basis(u1, u2) << std::endl;
	do_basis(u1); do_basis(u2);
	std::cout << "is basis (after division by norm) " << is_basis(u1, u2) << std::endl;
	u1.clear(); u2.clear();
	for (double t = a; t <= b; t += dt) {
		u1.push_back(sin(2 * pi * f3 * t));
		u2.push_back(sin(2 * pi * f4 * t));
	}
	std::cout << "is basis (after increasing frequencies) " << is_basis(u1, u2) << std::endl;
	u1.clear(); u2.clear();
	for (double t = - 2 * a; t <= 2 * b; t += dt) {
		u1.push_back(sin(2 * pi * f1 * t));
		u2.push_back(sin(2 * pi * f2 * t));
	}
	std::cout << "is basis (after increasing period) " << is_basis(u1, u2) << std::endl;
	u1.clear(); u2.clear();
	for (double t = - a / 2; t <= b / 2; t += dt) {
		u1.push_back(sin(2 * pi * f1 * t));
		u2.push_back(sin(2 * pi * f2 * t));
	}
	std::cout << "is basis (after decreasing period) " << is_basis(u1, u2) << std::endl;
	return 0;	
}
