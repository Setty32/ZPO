#include "thresholds.h"

using namespace std;

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

/**
 *  If you want to have at least N elements in 2nd group, choose to ignore N-1
 *  thresholds. Default is 0.
 */
unsigned getThreshold3(vector<vector<cv::Point>> &input, unsigned min)
{
	if (input.size() == 0)
    return 0;

	vector<unsigned> sizes;
  for (vector<vector<cv::Point>>::iterator it = input.begin(); it != input.end(); ++it)
		sizes.push_back(it->size());

	sort(sizes.begin(), sizes.end());

	unsigned index_a = 0, index_min = 0,  diff_a = 0, diff_min = 0, tmp = 0;
	for (vector<unsigned>::iterator it = sizes.begin() + 1; it != sizes.end(); ++it)
	{
		tmp = *it - *(it - 1);
		if (diff_a < tmp)
		{
			diff_a = tmp;
			index_a = it - sizes.begin() - 1;
		}
	}

	if (index_a > sizes.size() - min - 1)
	{
		for (vector<unsigned>::iterator it = sizes.begin() + 1, end_it = sizes.end() - min; it != end_it; ++it)
		{
			tmp = *it - *(it - 1);
			if (diff_min < tmp)
			{
				diff_min = tmp;
				index_min = it - sizes.begin() - 1;
			}
		}
		
		tmp = (sizes[index_a] + sizes[index_min] + (diff_a + diff_min)/2)/2;
	}
	else
		tmp = sizes[index_a] + diff_a / 2;

	return tmp;
}
