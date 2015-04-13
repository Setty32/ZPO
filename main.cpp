#include<iostream>
#include<string>
#include<vector>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "thresholds.h"

using namespace std;

void minMax(vector<cv::Point> &in, int &minX, int &maxX, int &minY, int &maxY){
	minX = INT_MAX;
	maxX = INT_MIN;
	minY = INT_MAX;
	maxY = INT_MIN;
	
	
	for(unsigned i = 0; i < in.size(); i++){
		minX = min(in.at(i).x, minX);
		maxX = max(in.at(i).x, maxX);
		
		minY = min(in.at(i).y, minY);
		maxY = max(in.at(i).y, maxY);
	}
	
}

int main(int argc, char* argv[]){
	string inFile = "";
	string outFile = "";

	if( argc > 1 ) inFile = std::string( argv[1] );
	if( argc > 2 ) outFile = std::string( argv[2] );

	cv::Mat src = cv::imread( inFile );
	cv::Mat assist;
	cv::Mat src_gray;
	cv::cvtColor( src, src_gray, CV_BGR2GRAY );

	cv::Mat dst=cv::Mat::zeros(src_gray.size(),src_gray.type());
	cv::Mat finalResult = cv::Mat::zeros(src_gray.size(),src_gray.type());

	cv::threshold(src_gray,dst,200,255,cv::THRESH_TOZERO);		//JEDEN PARAMETR == PRAH
	
	vector<vector<cv::Point> > contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(dst, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE,cv::Point(0,0));

	cv::Mat draw = cv::Mat::zeros(dst.size(),CV_8UC3);

//	cerr << contours.size() << endl;
	
	
	vector<cv::Point3d> direct;
	int step = 50;
	int cushion = 3; //manipulation space
	unsigned threshold = getThreshold2(contours);

	for(unsigned i = 0; i < contours.size(); i++){
		if(contours.at(i).size() > threshold){
			draw = cv::Mat::zeros(dst.size(),CV_8UC3);
	        cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));
	        
	        cv::Mat hough;
			cv::cvtColor( draw, hough, CV_BGR2GRAY );
//	draw.convertTo(hough,CV_8UC1);

			vector<cv::Vec2f> lines;
			cv::HoughLines(hough,lines, 1, CV_PI/360, 250, 0, 0);	//JEDEN PARAMETR == PRAH
		
			  for( size_t i = 0; i < lines.size(); i++ )
			  {
				float rho = lines[i][0], theta = lines[i][1];
				cv::Point pt1, pt2;
				double a = cos(theta), b = sin(theta);
				double x0 = a*rho, y0 = b*rho;
				pt1.x = cvRound(x0 + 1000*(-b));
				pt1.y = cvRound(y0 + 1000*(a));
				pt2.x = cvRound(x0 - 1000*(-b));
				pt2.y = cvRound(y0 - 1000*(a));
				
				if (a < 0.5 && y0 > hough.rows * 0.1 && y0 < hough.rows * 0.9)
					cv::line( draw, pt1, pt2, cv::Scalar(255,255,0), 1, CV_AA);
			 
			  } //for lines
		
			int maxX, minX, minY, maxY;
			minMax(contours.at(i), minX, maxX, minY, maxY);
			
			minX -= cushion;
			maxX += cushion;
			minY -= cushion;
			maxY += cushion;
			
			vector<vector<int>> points_in_slices;
			vector<int> delta_of_slice;
			points_in_slices.resize(maxX - minX);
			delta_of_slice.resize(maxX - minX);
			
			for(int j = 0; j < contours.at(i).size(); j++){
				
				points_in_slices.at(contours.at(i).at(j).x - minX).push_back(j);
				
			}
			
			
			//process one field by vertical slices
			for(int sloupec = minX; sloupec < maxX; sloupec++){
				
				int closest = -1; //closest point from slice to line
				int minDistance = INT_MAX;
				for(int chn = 0; chn < points_in_slices.at(sloupec).size(); chn++){
					
					//z pt1 a pt2 --->
					int distance = minDistance, linePointDistance(points_in_slices.at(sloupec).at(chn), --------);
					
					if(abs(distance) < minDistance){
					 closest = chn;
					 delta_of_slice[sloupec] = distance;
					}
					
				}
				
								
				for(int radek = minY; radek < maxY; radek++){
					//correction
				}
			}
			
		
		//cv::imshow("dst",hough); 
		cv::imshow("draw",draw);
		cv::waitKey(0);	

		}
	}
	

	return 0;
}
