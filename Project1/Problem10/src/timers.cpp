#include "../include/timers.hpp"

double time_general(int n){
      auto t1 = std::chrono::high_resolution_clock::now();

      int N = n;
      double h = 1.0/N; //step size
      double v0 = 0.0;  
      double vN = 0.0;

      std::vector<double> x;
      std::vector<double> g;
      for (double i = 0.0; i <= N; i++){
        //x takes steps values from 0 to 1. In total N + 1 points
        x.push_back(i/N);
      }
      for (int i = 1; i < N; i++){
          double value = h*h*100*exp(-10*x.at(i));
          g.push_back(value);
      }

      std::vector<double> a(N - 1, -1.0);
      std::vector<double> b(N - 1, 2);
      std::vector<double> c(N - 1, -1.0);
      std::vector <double> v = tridiagAlgo(a, b, c, g);

      auto t2 = std::chrono::high_resolution_clock::now();

      double duration = std::chrono::duration<double>(t2 - t1).count();

      return duration;
}

//double time_special(int n){
      //auto t1 = std::chrono::high_resolution_clock::now();

      //int N = n;
      //double h = 1.0/N; //step size
      //double v0 = 0.0;  
      //double vN = 0.0;

      //std::vector<double> x;
      //std::vector<double> g;
      //for (double i = 0.0; i <= N; i++){
        ////x takes steps values from 0 to 1. In total N + 1 points
        //x.push_back(i/N);
      //}
      //for (int i = 1; i < N; i++){
          //double value = h*h*100*exp(-10*x.at(i));
          //g.push_back(value);
      //}

      //std::vector <double> v = tridiagAlgo(a, b, c, g); //change to special

      //auto t2 = chrono::high_resolution_clock::now();

      //double duration = chrono::duration<double>(t2 - t1).count();

      //return duration;
//}
