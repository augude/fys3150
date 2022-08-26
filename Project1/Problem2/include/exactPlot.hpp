#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <tuple>

std::tuple<std::vector<double>, std::vector<double>> calculateExact(int steps, int prec);

void printToFile(int steps, int prec);