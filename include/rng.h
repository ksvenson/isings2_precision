#ifndef QFERNG_H
#define QFERNG_H

#include <iostream>
#include <vector>
#include <random>
#include <ctime>

//Copied essentially from Evan's code, with some addition and deleting unnecessary functions
class QfeRng {

public:
  QfeRng(int seed = 12345678);
  double RandReal(double min = 0.0, double max = 1.0);
  double RandNormal(double mean = 0.0, double stddev = 1.0);
  int RandInt(int min, int max);
  bool RandBool();
  std::vector<int> getRandomSpinArray(int size);
  double SampleExponential(double rate);

  std::mt19937 gen;
};

QfeRng::QfeRng(int seed) {
  gen = std::mt19937(seed);
}

double QfeRng::RandReal(double min, double max) {
  std::uniform_real_distribution<double> dist(min, max);
  return dist(gen);
}

double QfeRng::RandNormal(double mean, double stddev) {
  std::normal_distribution<double> dist(mean, stddev);
  return dist(gen);
}

int QfeRng::RandInt(int min, int max) {
  std::uniform_int_distribution<int> dist(min, max);
  return dist(gen);
}

bool QfeRng::RandBool() {
  std::uniform_int_distribution<int> dist(0, 1);
  return (dist(gen) == 1);
}

std::vector<int> QfeRng::getRandomSpinArray(int size) {
    // Define a uniform distribution to generate 0 or 1
    std::uniform_int_distribution<int> dist(0, 1);

    // Create a vector to store the random integers (-1 or +1)
    std::vector<int> randomIntArray;

    // Generate random integers and store them in the vector
    for (int i = 0; i < size; ++i) {
        int value = dist(gen) == 0 ? -1 : 1;
        randomIntArray.push_back(value);
    }
    return randomIntArray;
}

double QfeRng::SampleExponential(double rate){
    std::exponential_distribution<double> dist(rate);
    return dist(gen);
}

#endif