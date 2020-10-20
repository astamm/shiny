#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_diameter_upper_bound(const NumericVector &distanceMatrix,
                                const unsigned int numberOfPoints,
                                const unsigned int kNearestNeighbor = 10)
{
  double maxDiameter = 0.0;
  NumericVector workVector(numberOfPoints - 1);
  for (unsigned int i = 0;i < numberOfPoints;++i)
  {
    unsigned int pos = 0;

    for (unsigned int j = 0;j < numberOfPoints;++j)
    {
      if (i == j)
        continue;

      if (i < j)
        workVector[pos] = distanceMatrix[numberOfPoints * i - (i + 1) * i / 2 + j - i - 1];
      else
        workVector[pos] = distanceMatrix[numberOfPoints * j - (j + 1) * j / 2 + i - j - 1];

      ++pos;
    }

    std::nth_element(workVector.begin(), workVector.begin() + kNearestNeighbor - 1, workVector.end());
    double currentDiameter = workVector[kNearestNeighbor - 1];

    if (currentDiameter > maxDiameter)
      maxDiameter = currentDiameter;
  }

  return maxDiameter;
}
