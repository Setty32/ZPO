#include<iostream>
#include<string>
#include<vector>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

int main(int argc, char* argv[]){
	string inFile = "";
	string outFile = "";

	if( argc > 1 ) inFile = std::string( argv[1] );
	if( argc > 2 ) outFile = std::string( argv[2] );

	cv::Mat src = cv::imread( inFile );
	cv::Mat dst;
	src.copyTo(dst);

	cv::threshold(src,dst,220,255,cv::THRESH_TOZERO);
	
	vector<vector<cv::Point> > contours;
	vector<cv::Vec4i> hierarchy;

	//cv::findContours(dst,
	cv::imshow("dst",dst);
	cv::waitKey(0);	

	return 0;
}
