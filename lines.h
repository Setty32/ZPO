#ifndef LINES_H_
#define LINES_H_

#include <cmath>
#include <vector>
#include <map>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

void minMax(std::vector<cv::Point> &in, int &minX, int &maxX, int &minY, int &maxY);
double avgDistance(double k, double q, const std::vector<cv::Point> &contour, double x_limit = DBL_MAX, double delta = 5.0, double init_value = -1.0);
std::vector<int> validLines(std::vector<cv::Point> &contour, std::vector<cv::Vec2f> &lines);

#endif
