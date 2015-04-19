#include<iostream>
#include<string>
#include<vector>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "thresholds.h"
#include "lines.h"

#define filter_size 13

//sice vyfiltruje trochu sumu, ale pak pokud se neco prudceji zmeni
//tak chvili drzi a je tam patrny zlom..
void filter(std::vector<double> &input, std::vector<double> &output){
	
    double coeff[filter_size] = { 0.5, 0.5, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.5, 0.5};
	
    double weight_sum = 0;
	for(int i = 0; i < filter_size; i++){
		weight_sum += coeff[i];
	}
	
    int filterHalf = 6;
	
	for(unsigned i = 0; i < input.size(); i++){
		
        double result = 0.0;

        for(int j = - filterHalf; j <= filterHalf; j++){


            if(((int)i + j) < 0){
                result += coeff[j + filterHalf] * (double)input.at(0);
            } else if( (i + j) >= input.size() - 5){
                result += coeff[j + filterHalf] * (double)input.at(input.size() -5);
            } else {
                result += coeff[j + filterHalf] * (double)input.at(i + j);
            }
        }
        std::cerr<< input.at(i) <<" " <<  cvRound(result / weight_sum) << std::endl;
        output.at(i) = result / weight_sum;
	}
	
}


// http://geomalgorithms.com/a02-_lines.html
#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y)
#define d(u,v)      sqrt(dot(u-v,u-v))          // distance = norm of difference

double dist_Point_to_Line( cv::Point P, cv::Point L1, cv::Point L2)
{
     cv::Point v = L1 - L2;
     cv::Point w = P - L1;

     double c1 = dot(w,v);
     double c2 = dot(v,v);
     double b = c1 / c2;

     cv::Point Pb = L1 + b * v;
     if(b < 0)
		return d(P, Pb);
     else
		return -d(P,Pb);
}



using namespace std;
// PRESUNUL SOM DO LINES.CPP, POUZIL SOM TIEZ
/*void minMax(vector<cv::Point> &in, int &minX, int &maxX, int &minY, int &maxY){
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
	
}*/

int main(int argc, char* argv[]){
	string inFile = "";
	string outFile = "";

	if( argc > 1 ) inFile = std::string( argv[1] );
	if( argc > 2 ) outFile = std::string( argv[2] );

    cv::Mat src = cv::imread( inFile, CV_LOAD_IMAGE_COLOR);

    cerr<< "Input img sizes: " << src.cols << " " << src.rows << endl;
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
	int cushion = 3; //manipulation space
	unsigned threshold = getThreshold3(contours, 4); // we want at least 5 contours so we ignore 4 thresholds

	for(unsigned i = 0; i < contours.size(); i++){
		//if(contours.at(i).size() > 740 /*threshold*/){
		if(contours.at(i).size() > threshold){
			draw = cv::Mat::zeros(dst.size(),CV_8UC3);
	        cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));
	        
	        cv::Mat hough;
			cv::cvtColor( draw, hough, CV_BGR2GRAY );
//	draw.convertTo(hough,CV_8UC1);
			bool stop = true;
			
			vector<cv::Vec2f> lines;
			cv::HoughLines(hough,lines, 1, CV_PI/360, 250, 0, 0);	//JEDEN PARAMETR == PRAH
			draw = cv::Mat::zeros(dst.size(),CV_8UC3);

			std::vector<int> goodLines = validLines(contours[i], lines); // FUNKCIA NA VYBER INDEXOV TYCH NAJ LINII

			  //for( size_t j = 0; j < lines.size(); j++ ) // POVODNE
			  for (std::vector<int>::iterator it = goodLines.begin(); it != goodLines.end(); ++it)
			  {
				//float rho = lines[j][0], theta = lines[j][1]; // POVODNE
				float rho = lines[*it][0], theta = lines[*it][1];
				cv::Point pt1, pt2;
				double a = cos(theta), b = sin(theta);
				double x0 = a*rho, y0 = b*rho;
				pt1.x = cvRound(x0 + 1000*(-b));
				pt1.y = cvRound(y0 + 1000*(a));
				pt2.x = cvRound(x0 - 1000*(-b));
				pt2.y = cvRound(y0 - 1000*(a));
				
				if (a < 0.5 && y0 > hough.rows * 0.1 && y0 < hough.rows * 0.9){
					cv::line( draw, pt1, pt2, cv::Scalar(255,255,0), 1, CV_AA);
					stop = false;
			 }
			  } //for lines
			  
			if(lines.size() < 1 || stop){
				continue;
			}
			int maxX, minX, minY, maxY;
			minMax(contours.at(i), minX, maxX, minY, maxY);
			
			minX -= cushion;
			maxX += cushion;
			minY -= cushion;
			maxY += cushion;
			
			vector<vector<int> > points_in_slices;
            vector<double> delta_of_slice;
			points_in_slices.resize(maxX - minX);
			delta_of_slice.resize(maxX - minX);
			
			for(unsigned j = 0; j < contours.at(i).size(); j++){
				
				points_in_slices.at(contours.at(i).at(j).x - minX).push_back(j);
				
			}
			
			
			//process one field by vertical slices
			for(int slice = 5; slice < maxX - minX; slice++){
				
				int minDistance = INT_MAX;
				for(unsigned point = 0; point < points_in_slices.at(slice).size(); point++){
					
					//for(unsigned vec = 0; vec < lines.size(); vec++){ // POVODNE
			  	for (std::vector<int>::iterator it = goodLines.begin(); it != goodLines.end(); ++it){
						//double rho = lines[vec][0], theta = lines[vec][1]; // POVODNE
						double rho = lines[*it][0], theta = lines[*it][1];
						double a = cos(theta), b = sin(theta);
						double x0 = a*rho, y0 = b*rho;
						
						//contours.at(i).at(contourPoint).x - x0 / -b = k;
						//y = y0 + k *(a);
						
						unsigned contourPoint = points_in_slices.at(slice).at(point);
						//double m = dist_Point_to_Line(contours.at(i).at(contourPoint),cv::Point(cvRound(x0 - 1000*(-b)),cvRound(y0 - 1000*(a))), cv::Point(cvRound(x0 + 1500*(-b)),cvRound(y0 + 1500*(a))));
						
						double k = (contours.at(i).at(contourPoint).x - x0) / -b;
						double y = y0 + k * a;
						
						double m = y - contours.at(i).at(contourPoint).y;
						
						if(abs(m) < abs(minDistance)){
							delta_of_slice.at(slice) = m;
							minDistance = m;
						}
						
					}//for lines
		
					
				} //for points in slices
            }

            vector<double> newDeltaSlice(delta_of_slice.size());
            filter(delta_of_slice,newDeltaSlice);
            for(int slice = 5; slice < maxX - minX; slice++){
//				for(unsigned point = 0; point < points_in_slices.at(slice).size(); point++){
                    if(abs(newDeltaSlice.at(slice)) > 0){
						//vizualizace opravenych kontur
//						unsigned contourPoint = points_in_slices.at(slice).at(point);
//                        contours.at(i).at(contourPoint).y += newDeltaSlice.at(slice);
						
						int src_X = slice + minX;
						int height = maxY - minY;

                        float delta = newDeltaSlice.at(slice);
                        float fraction = newDeltaSlice.at(slice);
                        int wholePart = (int) fraction;
                        fraction = abs(fraction - wholePart);
                        float complement = 1 - fraction;

                        //cout << newDeltaSlice.at(slice) << " "<< wholePart <<" "<< fraction << " " << complement << endl;
						for(int offset = 0; offset < height; offset++){

                            //above line
                            if(delta > 0){
                                src_gray.at<uchar>(maxY + wholePart - offset, src_X) = src_gray.at<uchar>(maxY - offset, src_X);

                             //below line
                            } else if (delta < 0){
                                src_gray.at<uchar>(minY + wholePart + offset, src_X) = src_gray.at<uchar>(minY + offset, src_X);
                            }

                        }

                        //go through slice once again and fine adjust values..
                        for(int offset = 0; offset < height; offset++){

                            //above line
                            if(delta > 0){
                                float v1 = complement * src_gray.at<uchar>(maxY - offset, src_X);
                                float v2 = fraction * src_gray.at<uchar>(maxY - offset - 1, src_X);
                                float value = v1 + v2;
                                src_gray.at<uchar>(maxY - offset, src_X) = (uchar) value;

                             //below line
                            } else if (delta < 0){
                                float v1 = complement * src_gray.at<uchar>(minY + offset, src_X);
                                float v2 =  fraction * src_gray.at<uchar>(minY + offset + 1, src_X);
                                float value = v1 + v2;
                                src_gray.at<uchar>(minY + offset, src_X) = (uchar) value;
                            }

                        }
					}
//				}
//                cerr << delta_of_slice.at(slice) << endl;
				
			} // for all slices
			
                cerr << "------------------------" << endl;
            cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));
		
		//cv::imshow("dst",hough); 
                //cv::imshow("draw",draw);
                //cv::waitKey(0);

		}
		
	} //for all contours
	
	
	
//	draw = cv::Mat::zeros(dst.size(),CV_8UC3);
//	for(unsigned i = 0; i < contours.size(); i++){
//		if(contours.at(i).size() > 740 /*threshold*/)
//		cv::drawContours( src_gray, contours, i, cv::Scalar(0,0,255));
//	}
	
	//cv::imshow("draw",draw);
			cv::imshow("dst",src_gray);
		cv::waitKey(0);

	return 0;
}
