#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <algorithm>

#include "edge_detect.h"
#include "edge_boxes.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <exception>

// np_opencv_converter
#include "np_opencv_converter.hpp"

namespace py = boost::python;

using namespace std;
using namespace cv;

class edge_boxes_python {
    public:
        string st_path; // path to the trained structured forest for edge detection
        edge_boxes_python(string st_path_) : st_path(st_path_) {}
        
        Mat get_edge_boxes(const Mat &im) {
            Mat ime, grad_ori, ime_t, grad_ori_t;

            // setup and run EdgeBoxGenerator
            EdgeBoxGenerator edgeBoxGen; Boxes boxes;
            edgeBoxGen._alpha = 0.65; 
            edgeBoxGen._beta = 0.75;
            edgeBoxGen._eta = 1;
            edgeBoxGen._minScore = 0.01;
            edgeBoxGen._maxBoxes = 10000;
            edgeBoxGen._edgeMinMag = 0.1;
            edgeBoxGen._edgeMergeThr = 0.5;
            edgeBoxGen._clusterMinMag = 0.5;
            edgeBoxGen._maxAspectRatio = 3;
            edgeBoxGen._minBoxArea = 1000;
            edgeBoxGen._gamma = 2;
            edgeBoxGen._kappa = 1.5;

            double t = (double)getTickCount();
            edge_detect(im, ime, grad_ori, st_path);
            //vis_matrix(ime, "E");
            transpose(ime, ime_t);
            transpose(grad_ori, grad_ori_t);

            if(!(ime_t.isContinuous() && grad_ori_t.isContinuous())) {
                cout << "Matrices are not continuous, hence the Array struct will not work" << endl; 
            }

            arrayf E; E._x = ime_t.ptr<float>();
            arrayf O; O._x = grad_ori_t.ptr<float>();
            Size sz = ime.size();
            int h = sz.height; O._h=E._h=h;
            int w = sz.width; O._w=E._w=w;

            arrayf V;

            edgeBoxGen.generate( boxes, E, O, V );
            t = ((double)getTickCount() - t)/getTickFrequency();
            //cout << "Generated boxes, t = " << t*1000 << " ms" << endl;

            // create output bbs
            int n = (int) boxes.size();
            //cout << "Found " << n << " boxes" << endl;
            Mat bbs(n, 5, CV_32FC1);
            for(int i=0; i<n; i++) {
                bbs.at<float>(i, 0) = (float) boxes[i].c+1;
                bbs.at<float>(i, 1) = (float) boxes[i].r+1;
                bbs.at<float>(i, 2) = (float) boxes[i].w;
                bbs.at<float>(i, 3) = (float) boxes[i].h;
                bbs.at<float>(i, 4) = boxes[i].s;
            }

            return bbs;
        }
};

namespace fs { 
    namespace python {
        BOOST_PYTHON_MODULE(edge_boxes_python) {
            // Main types export
            fs::python::init_and_export_converters();
            py::scope scope = py::scope();
            //py::def("get_edge_boxes", &get_edge_boxes);
            py::class_<edge_boxes_python>("edge_boxes_python", py::init<std::string>())
                .def("get_edge_boxes", &edge_boxes_python::get_edge_boxes);
        }
    }
}
