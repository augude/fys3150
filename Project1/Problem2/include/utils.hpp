#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>

std::string scientificFormat(double d, const int width = 20, const int prec = 10);

std::string scientificFormat(const std::vector<double>& v, const int width = 20, const int prec = 10);