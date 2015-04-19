#include "lines.h"

using namespace std;

void minMax(vector<cv::Point> &in, int &minX, int &maxX, int &minY, int &maxY)
{
	minX = INT_MAX;
	maxX = INT_MIN;
	minY = INT_MAX;
	maxY = INT_MIN;
	
	for(unsigned i = 0; i < in.size(); i++)
	{
		minX = min(in.at(i).x, minX);
		maxX = max(in.at(i).x, maxX);
		
		minY = min(in.at(i).y, minY);
		maxY = max(in.at(i).y, maxY);
	}
}

/**
 *  Compute average distance between line and contour
 */
double avgDistance(double k, double q, const vector<cv::Point> &contour, double x_limit, double delta, double init_value)
{
	unsigned cnt;
	double avg_dist, tmp;

	if (init_value < 0.0)
	{
		cnt = 0;
		avg_dist = 0.0;
	}
	else
	{
		cnt = 1;
		avg_dist = init_value;
	}

	for (unsigned j = 0; j < contour.size(); ++j)
	{
		if (contour[j].x < x_limit)
		{
			tmp = abs(k*contour[j].x + q - contour[j].y);

			if (tmp < delta)
			{
				avg_dist += tmp;
				++cnt;
			}
		}
	}

	if (cnt != 0)
		avg_dist /= cnt;

	return avg_dist;
}

/**
 *  Opperates under the assumption that input lines are relatively closely
 *  fitted already to one side of the contour or another. In reality only at 
 *  maximum two line groups usually exist. Lines are intersecting within 
 *  these groups and groups are in turn parallel to one another.
 */
vector<int> validLines(vector<cv::Point> &contour, vector<cv::Vec2f> &lines)
{
	int minX, minY, maxX, maxY;
	minMax(contour, minX, maxX, minY, maxY);

	bool done = false;
	vector<double> k, q;
	double x, y, fraction, integer, avg_dist;
	const double epsilon_k = 0.001, epsilon_q = 1.0, delta = 5.0; // limit of difference before equality
	vector<multimap<double, int>> line_groups;

	for (unsigned i = 0; i < lines.size(); ++i)
	{
		// check if theta isn't close to multiple of Pi
		fraction = modf(lines[i][1] / M_PI, &integer);
		if (fraction < epsilon_k)
		{
			// will never be used anyway 'cause we are interested only in horizontal lines
			k.push_back(numeric_limits<double>::max());
			q.push_back(numeric_limits<double>::min());
			continue;
		}

		// compute parameters of line's slope representation
		k.push_back(tan(lines[i][1] - M_PI/2));
		q.push_back((-k[i]*cos(lines[i][1]) + sin(lines[i][1]))*lines[i][0]);

		// compute average distance between line and closest part of contour
		avg_dist = avgDistance(k[i], q[i], contour);

		// put line into appropriate group
		for (unsigned m = 0; !done && m < line_groups.size(); ++m)
			for (multimap<double, int>::iterator it = line_groups[m].begin(); !done && it != line_groups[m].end(); ++it)
			{
				if (abs(k[i] - k[it->second]) < epsilon_k)
				{
					// it's the same line and is already in the structure
					if (abs(q[i] - q[it->second]) < epsilon_q)
						done = true;
					else
						continue;
				}
				else
				{
					x = (q[it->second] - q[i])/(k[i] - k[it->second]);
					y = k[i]*x + q[i];

					if (round(x) >= minX && round(x) <= maxX && round(y) >= minY && round(y) <= maxY)
					{
						// average difference on whole closest par of contour and its part before intersection
						avg_dist = avgDistance(k[i], q[i], contour, x, delta, avg_dist);
						line_groups[m].insert(make_pair(avg_dist, i));

						// do that for the line already in the structure as well
						avg_dist = avgDistance(k[it->second], q[it->second], contour, x, delta, it->first);
						line_groups[m].insert(make_pair(avg_dist, it->second));
						line_groups[m].erase(it);

						done = true;
					}
					// it's the same line otherwise and is already in the structure
				}
			}

		if (!done) // if you managed to get here, you are member of new group, congrats :)
		{
			multimap<double, int> tmp;
			tmp.insert(make_pair(avg_dist, i));
			line_groups.push_back(tmp);
		}
		else
			done = false;
	}

	multimap<double, int> tmp;
	// pick index of the best line from each group
	for (unsigned m = 0; m < line_groups.size(); ++m)
		tmp.insert(*(line_groups[m].begin()));

	vector<int> indexes;
	for (multimap<double, int>::iterator it = tmp.begin(); it != tmp.end(); ++it)
		indexes.push_back(it->second);

	return indexes;
}
