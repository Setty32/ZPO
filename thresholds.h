#ifndef THRESHOLDS_H_
#define THRESHOLDS_H_

#include <vector>
#include <cmath>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

unsigned getThreshold(std::vector<std::vector<cv::Point>> &input);
unsigned getThreshold2(std::vector<std::vector<cv::Point>> &input);
unsigned getThreshold3(std::vector<std::vector<cv::Point>> &input, unsigned min = 0);

#endif
