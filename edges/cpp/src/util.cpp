#include <opencv2/highgui/highgui.hpp>
#include "edge_detect.h"

using namespace std;
using namespace cv;

// function to calculate the numerical gradient of a matrix in X and Y directions
void gradient(Mat &I, Mat &gx, Mat &gy) {
    gx = Mat(I.rows, I.cols, CV_32FC1);
    gy = Mat(I.rows, I.cols, CV_32FC1);
    
    // compute rows of Gy 
    for(int i = 1; i < I.rows-1; i++) gy.row(i) = (I.row(i+1) - I.row(i-1)) * 0.5;
    gy.row(0) = I.row(1) - I.row(0);
    gy.row(gy.rows-1) = I.row(I.rows-1) - I.row(I.rows-2);

    // compute the columns of Gx
    for(int i = 1; i < I.cols-1; i++) gx.col(i) = (I.col(i+1) - I.col(i-1)) * 0.5;
    gx.col(0) = I.col(1) - I.col(0);
    gx.col(gx.cols-1) = I.col(I.cols-1) - I.col(I.cols-2);
}

// function to visualize a CV_32FC1 matrix
void vis_matrix(Mat &m, char *window_name) {
    double min_val, max_val;
    minMaxLoc(m, &min_val, &max_val);
    Mat m_show;
    convertScaleAbs(m-min_val, m_show, 255.0/(max_val - min_val));
    namedWindow(window_name, WINDOW_NORMAL);
    imshow(window_name, m_show);
}

// function to implement the elementwise signum function
Mat signum(Mat &src) {
    Mat dst = Mat::zeros(src.size(), CV_32FC1);
    add(dst, 1, dst, src>0);
    add(dst, -1, dst, src<0);
    //dst += ((src > 0) & 1);
    //dst -= ((src < 0) & 1);

    return dst;
} 
