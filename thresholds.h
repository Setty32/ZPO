#include <vector>
#include <cmath>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

unsigned getThreshold(vector<vector<cv::Point>> &input);
unsigned getThreshold2(vector<vector<cv::Point>> &input);
