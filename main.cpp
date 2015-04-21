#include<iostream>
#include<string>
#include<vector>
#include <sys/stat.h>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "thresholds.h"
#include "lines.h"

//use this to build for presentation..
#define DEMO



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
        //std::cerr<< input.at(i) <<" " <<  cvRound(result / weight_sum) << std::endl;
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


int main(int argc, char* argv[]){
	string inFile = "";
	string outFile = "";

	if( argc > 1 ) inFile = std::string( argv[1] );
	if( argc > 2 ) outFile = std::string( argv[2] );

    if(argc == 1){
        cout << "need at least input file\n";
        return 1;
    }



    struct stat statbuf;

    if (stat("./out", &statbuf) != -1) {
       if (not S_ISDIR(statbuf.st_mode)) {

           const int dir_err = mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
           if (-1 == dir_err)
           {
               printf("Error creating directory!n");
               exit(1);
           }
       }
    }


    cv::Mat src = cv::imread( inFile, CV_LOAD_IMAGE_COLOR);

    //get rid of path at the beginning to keep only filename
     if(inFile.size() > 2){
         size_t pos = inFile.rfind("/");

         if(pos != string::npos)
            inFile.erase(0, pos + 1);
    }

    //cerr<< "Input img sizes: " << src.cols << " " << src.rows << endl;

    //convert input file to greyscale
    cv::Mat src_gray;
	cv::cvtColor( src, src_gray, CV_BGR2GRAY );

    if( src.cols < src.rows){
        cv::Mat tmp(src_gray);

        cv::transpose(tmp, src_gray);
        flip(src_gray, src_gray,1);
    }

    cv::Mat dst = cv::Mat::zeros(src_gray.size(),src_gray.type());
	
#ifdef DEMO
cv::imshow("original", src_gray);
cv::waitKey();

//helper metrices for visualization
cv::Mat tmp_draw = cv::Mat::zeros(dst.size(),CV_8UC3);
cv::Mat tmp_draw2 = cv::Mat::zeros(dst.size(),CV_8UC3);
cv::Mat tmp_draw3 = cv::Mat::zeros(dst.size(),CV_8UC1);
cv::Mat tmp_draw4 = cv::Mat::zeros(dst.size(),CV_8UC1);
#endif

    //isolate white areas; dst is binary image now
	cv::threshold(src_gray,dst,200,255,cv::THRESH_TOZERO);		//JEDEN PARAMETR == PRAH

#ifdef DEMO
    //after thresholding
cv::imshow("progress", dst);
cv::waitKey();
#endif

	vector<vector<cv::Point> > contours;
	vector<cv::Vec4i> hierarchy;
    vector<cv::Vec2f> lines;
    std::vector<int> goodLines;


    //find contours of these areas and only outer ones (no inner contours). Stored as points.
	cv::findContours(dst, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE,cv::Point(0,0));

	cv::Mat draw = cv::Mat::zeros(dst.size(),CV_8UC3);

//	cerr << contours.size() << endl;
	
	
    int cushion = 3;    //manipulation space to add to bounding box

    //want to keep at least 5 biggest contours by point count and cut off everything else..
    unsigned threshold = getThreshold3(contours, 4);

#ifdef DEMO

for(unsigned i = 0; i < contours.size(); i++){
    cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));
}

//show all contours
cv::imshow("progress", draw);
cv::waitKey();

draw = cv::Mat::zeros(dst.size(),CV_8UC3);

for(unsigned i = 0; i < contours.size(); i++){

    if(contours.at(i).size() < threshold)
        continue;

    cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));
}

//and all contours that pass threshold
cv::imshow("progress", draw);
cv::waitKey();


draw.copyTo(tmp_draw);
draw.copyTo(tmp_draw2);
#endif

	for(unsigned i = 0; i < contours.size(); i++){

        cv::Mat hough, tmp; //matrices for finding lines

        //so here we use only those 5+ biggest ones - can be more than 5 - in case of slight rotation of image, border contours are there too
        if(contours.at(i).size() < threshold)
            continue;


        tmp = cv::Mat::zeros(dst.size(),CV_8UC3);
        cv::drawContours( tmp, contours, i, cv::Scalar(0,0,255));


        cv::cvtColor( tmp, hough, CV_BGR2GRAY );

        bool stop = true;


        cv::HoughLines(hough,lines, 1, CV_PI/360, 250, 0, 0);	//JEDEN PARAMETR == PRAH
        tmp = cv::Mat::zeros(dst.size(),CV_8UC3);

#ifdef DEMO

for(unsigned pp = 0; pp < lines.size(); pp++){
    float rho = lines[pp][0], theta = lines[pp][1];
    cv::Point pt1, pt2;
    double a = cos(theta), b = sin(theta);
    double x0 = a*rho, y0 = b*rho;
    pt1.x = cvRound(x0 + 1000*(-b));
    pt1.y = cvRound(y0 + 1000*(a));
    pt2.x = cvRound(x0 - 1000*(-b));
    pt2.y = cvRound(y0 - 1000*(a));

    //all detected lines
    cv::line( tmp_draw, pt1, pt2, cv::Scalar(255,255,0), 1, CV_AA);
}

#endif

        //in case of multiple detected lines for any contour, select only the "best ones"
        goodLines = validLines(contours[i], lines);

        for (auto it = goodLines.begin(); it != goodLines.end(); ++it){

            float rho = lines[*it][0], theta = lines[*it][1];
            cv::Point pt1, pt2;
            double a = cos(theta), b = sin(theta);
            double x0 = a*rho, y0 = b*rho;
            pt1.x = cvRound(x0 + 1000*(-b));
            pt1.y = cvRound(y0 + 1000*(a));
            pt2.x = cvRound(x0 - 1000*(-b));
            pt2.y = cvRound(y0 - 1000*(a));

            //filters out fields with lines with more than 45 deg inclination and at top and bottmo 10% of space
            if (a < 0.5 && y0 > hough.rows * 0.1 && y0 < hough.rows * 0.9){
                #ifdef DEMO
                    //all good lines
                    cv::line( tmp_draw2, pt1, pt2, cv::Scalar(255,255,0), 1, CV_AA);

                #endif
                stop = false;       //basically is set only for border contours
            }
        } //for lines

        //if no lines or unfit go on..
        if(lines.size() < 1 || stop)
            continue;


        //find out bounding box of contour and inflate it a bit
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

        //reprocess contour points into slices - eg put them into vector at their X index
        for(unsigned j = 0; j < contours.at(i).size(); j++)
            points_in_slices.at(contours.at(i).at(j).x - minX).push_back(j);



        //choose to ignore first few slices - mostly vertical and unimportant
        for(int slice = 5; slice < maxX - minX; slice++){

            //for all points in each slice find minimum distance to any line and save it as delta of the whole slice
            int minDistance = INT_MAX;
            for(unsigned point = 0; point < points_in_slices.at(slice).size(); point++){

                for (auto it = goodLines.begin(); it != goodLines.end(); ++it){

                    double rho = lines[*it][0], theta = lines[*it][1];
                    double a = cos(theta), b = sin(theta);
                    double x0 = a*rho, y0 = b*rho;

                    //contours.at(i).at(contourPoint).x - x0 / -b = k;
                    //y = y0 + k *(a);

                    unsigned contourPoint = points_in_slices.at(slice).at(point);

                    double k = (contours.at(i).at(contourPoint).x - x0) / -b;
                    double y = y0 + k * a;

                    double m = y - contours.at(i).at(contourPoint).y;

                    if(abs(m) < abs(minDistance)){
                        delta_of_slice.at(slice) = m;
                        minDistance = m;
                    }

                }//for lines


            } //for points in slices
        } //for slice


        //now because there is inevitably some noise in contours and deltas due to thresholding, lack of sharpnes etc.,
        //run it through linear filter to get rid of some of this noise.
        vector<double> newDeltaSlice(delta_of_slice.size());
        filter(delta_of_slice,newDeltaSlice);


        //now to get to repairing image.. take each slice again
        for(int slice = 5; slice < maxX - minX; slice++){

        //v demu to vypada blbe bez casti s deltou 0
        #ifndef DEMO
            if(abs(newDeltaSlice.at(slice)) < 0.0001f) //floatove porovnavani
                continue;
        #endif


            int src_X = slice + minX;
            int height = maxY - minY;

            //breaking the delta to whole part and fraction + its complement
            float delta = newDeltaSlice.at(slice);
            float fraction = newDeltaSlice.at(slice);
            int wholePart = (int) fraction;
            fraction = abs(fraction - wholePart);
            float complement = 1 - fraction;

            //shift by whole part
            for(int offset = 0; offset < height; offset++){

                //above line
                if(delta > 0){
                    src_gray.at<uchar>(maxY + wholePart - offset, src_X) = src_gray.at<uchar>(maxY - offset, src_X);

                 //below line
                } else if (delta < 0){
                    src_gray.at<uchar>(minY + wholePart + offset, src_X) = src_gray.at<uchar>(minY + offset, src_X);
                }

            }
            #ifdef DEMO
                for(int offset = 0; offset < height; offset++){

                    //above line
                    if(delta > 0){
                        tmp_draw3.at<uchar>(maxY + wholePart - offset, src_X) = src_gray.at<uchar>(maxY - offset, src_X);
                        tmp_draw4.at<uchar>(maxY + wholePart - offset, src_X) = tmp_draw3.at<uchar>(maxY + wholePart - offset, src_X);

                     //below line
                    } else if (delta < 0){
                        tmp_draw3.at<uchar>(minY + wholePart + offset, src_X) = src_gray.at<uchar>(minY + offset, src_X);
                        tmp_draw4.at<uchar>(minY + wholePart + offset, src_X) = tmp_draw3.at<uchar>(minY + wholePart + offset, src_X);

                    }

                }

            #endif

            //go through slice once again and fine adjust values.. - linear interpolation
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

            #ifdef DEMO
                for(int offset = 0; offset < height; offset++){

                    //above line
                    if(delta > 0){

                        tmp_draw4.at<uchar>(maxY - offset, src_X) = src_gray.at<uchar>(maxY - offset, src_X);

                     //below line
                    } else if (delta < 0){

                        tmp_draw4.at<uchar>(minY + offset, src_X) = src_gray.at<uchar>(minY + offset, src_X);
                    }

                }

            #endif

        } // for all slices

        //cv::drawContours( draw, contours, i, cv::Scalar(0,0,255));




		
	} //for all contours
	
#ifdef DEMO

    char key = ' ';

    cv::imshow("progress",tmp_draw);
    cv::waitKey();

    cv::imshow("progress",tmp_draw2);
    cv::waitKey();

    while(key == ' '){
        cv::imshow("progress",tmp_draw3);
        cv::waitKey();

        cv::imshow("progress",tmp_draw4);
        key = cv::waitKey();
    }

    cv::imshow("result",src_gray);
    cv::waitKey();
#endif


	



    if(outFile == "")
        outFile = inFile;

    vector<int> params;
    params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    params.push_back(9);

    int code = cv::imwrite("./out/repaired_"+outFile, src_gray, params);


    return code;
}
