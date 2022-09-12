#pragma once
#include <iostream>
#include <armadillo>
#include <map>
#include <math.h>

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);