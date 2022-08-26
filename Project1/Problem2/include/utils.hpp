#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>

std::string scientific_format(const double d, const int width = 20, const int prec = 10);

std::string scientific_format(const std::vector<double>& v, const int width = 20, const int prec = 10);