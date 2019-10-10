/*******************************************************************************
* Structured Edge Detection Toolbox      Version 3.01
* Code written by Piotr Dollar and Larry Zitnick, 2014.
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/

// Modified by Samarth Brahmbhatt to remove MEX parts

#include <opencv2/highgui/highgui.hpp>

#include "contour/structuredforest.h"
#include "imgproc/image.h"
#include "edge_detect.h"

using namespace std;
using namespace cv;

// function to convert Eigen matrix to OpenCV Mat
template<typename _Tp, int _rows, int _cols, int _options, int _maxRows, int _maxCols>
void eigen2cv(Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols>& src, cv::Mat& dst) {
    Mat _src(src.rows(), src.cols(), DataType<_Tp>::type, (void*)src.data(), src.stride()*sizeof(_Tp));
    _src.copyTo(dst);
}

Mat sf_edges(const Mat& im, string st_path) {
    // use GOP library to get Structured Forest Edges
    Mat im_rgb;
    cvtColor(im, im_rgb, CV_BGR2RGB);
    Image8u im_8u(im_rgb.data, im_rgb.cols, im_rgb.rows, im_rgb.channels());

    // StructuredForestSettings sf_settings(stride, shrink, out_patch_size, feature_patch_size, patch_smooth, sim_smooth, sim_cells);
    StructuredForestSettings sf_settings(2, 2, 16, 32, 2, 8, 5);
    MultiScaleStructuredForest detector(1, -1, sf_settings);
    //detector.load("/home/samarth/Documents/MATLAB/edges/cpp/external/gop_1.3/data/sf.dat");
    detector.load(st_path.c_str());
    RMatrixXf im_gop_e = detector.detectAndFilter(im_8u);

    Mat im_e;
    eigen2cv(im_gop_e, im_e);

    return im_e;
}

Mat coarse_ori(Mat E) { // get coarse orientation, see Piotr Dollar's edgesDetect.m
    int r = 4;
    copyMakeBorder(E, E, r, r, r, r, BORDER_REFLECT);

    //f = [1:r r+1 r:-1:1]/(r+1)^2;
    Mat f = (r+1) * Mat::ones(1, 2*r+1, CV_32F);
    for(int i = 0; i < r; i++) {
        f.at<float>(0, i) = i+1;
        f.at<float>(0, f.cols-i-1) = i+1;
    }
    sepFilter2D(E, E, CV_32F, f, f.t());
    //vis_matrix(E, "E_conv");

    Mat Ox, Oy, Oxx, Oxy, Oyy;
    gradient(E, Ox, Oy);
    gradient(Ox, Oxx, Oxy);
    gradient(Oy, Oxy, Oyy);

    //vis_matrix(Ox, "Ox");
    //vis_matrix(Oy, "Oy");

    //O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
    /*
    Mat M, O;
    cartToPolar(Oxx + 1e-5, Oyy.mul(-signum(Oxy)), M, O);
    // convert from atan2 to atan
    add(O, -PI, O, Oxx<0);
    //O -= PI*((Oxx < 0) & 1);
    // mod pi
    add(O, PI, O, O<0);
    //O += PI*((O < 0) & 1);
    */
    Mat M, O;
    cartToPolar(Oxx + 1e-7, Oyy, M, O);
    //O -= PI * ((O > PI) & 1);
    add(O, -PI, O, O>PI);
    return O;
}

void edge_detect(const Mat &im, Mat &E, Mat &O, string st_path) {
    // get structured forest edges
    E = sf_edges(im, st_path);
    E.setTo(0, E < 0.5);
    //vis_matrix(E, "E");
    // get edge orientation
    O = coarse_ori(E);
    //vis_matrix(O, "O");
    // NMS on edges
    //E = edge_nms(E, O, 2, 0, 1, 4); 
    //vis_matrix(E, "E_nms");
}
