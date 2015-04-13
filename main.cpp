#include<iostream>
#include<string>
#include<vector>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "thresholds.h"

using namespace std;

int main(int argc, char* argv[]){
	string inFile = "";
	string outFile = "";

	if( argc > 1 ) inFile = std::string( argv[1] );
	if( argc > 2 ) outFile = std::string( argv[2] );

	cv::Mat src = cv::imread( inFile );
	cv::Mat assist;
	cv::Mat src_gray;
	cv::cvtColor( src, src_gray, CV_BGR2GRAY );

	cv::Mat dst=cv::Mat::zeros(src_gray.size(),src_gray.type());;
	

	cv::threshold(src_gray,dst,220,255,cv::THRESH_TOZERO);
	
	
	vector<vector<cv::Point> > contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(dst, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE,cv::Point(0,0));

	cv::Mat draw = cv::Mat::zeros(dst.size(),CV_8UC3);

	cerr << contours.size() << endl;
	
	
	unsigned threshold = getThreshold(contours);
	for(unsigned i = 1; i < contours.size(); i++){
		//if(contours.at(i).size() > 100)
		if(contours.at(i).size() > threshold)
	        cv::drawContours( draw, contours, i, cv::Scalar(0,0,255),1,8,hierarchy);
	}
	
	cv::imshow("dst",draw);
	cv::waitKey(0);	

	return 0;
}
