#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <armadillo>

std::string scientificFormat(double d, const int width = 20, const int prec = 10);

std::string scientificFormat(const std::vector<double>& v, const int width = 20, const int prec = 10);

std::string scientificFormat(const arma::vec v, const int width = 20, const int prec = 10);

std::string scientificFormat(const std::vector<std::vector<double>>& v, const int width = 20, const int prec = 10);
