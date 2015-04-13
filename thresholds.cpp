#include "thresholds.h"

unsigned getThreshold(vector<vector<cv::Point>> &input)
{
  unsigned threshold = 0;

  for (vector<vector<cv::Point>>::iterator it = input.begin(); it != input.end(); ++it)
    threshold += it->size()*it->size();

  threshold /= input.size();
  
  return pow(threshold, 3.14159265359/5);
}

unsigned getThreshold2(vector<vector<cv::Point>> &input)
{
  const double gr = 1.61803398875;
  unsigned min = 0, max = 0;

  for (vector<vector<cv::Point>>::iterator it = input.begin(); it != input.end(); ++it)
  {
    if (min > it->size())
      min = it->size();
    if (max < it->size())
      max = it->size();
  }

  return (max - min) / gr;
}
