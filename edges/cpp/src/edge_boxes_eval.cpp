#include <opencv2/highgui/highgui.hpp>
#include <boost/filesystem.hpp>

#include "edge_detect.h"
#include "edge_boxes.h"

using namespace std;
using namespace cv;
using namespace boost::filesystem;

int main(int argc, char **argv) {
    if(argc == 1) {
        cout << "Usage: ./edge_boxes test_image_path" << endl << "Produces results in test_image_path/../results_cpp/" << endl;
        return -1;
    }

    string path(argv[1]);
    cout << "Path is " << path << endl;

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

    int count = 0;
    for(directory_iterator i(path), end_iter; i != end_iter; i++, count++) {
        string filename = path + i->path().filename().string();
        Mat im = imread(filename), ime, grad_ori, ime_t, grad_ori_t;
        if(im.data == NULL) {
            cout << "Error reading image" << endl;
            return -1;
        }

        edge_detect(im, ime, grad_ori, string("/home/samarth/Documents/MATLAB/edges/cpp/external/gop_1.3/data/sf.dat"));
        //vis_matrix(ime, "E");
        transpose(ime, ime_t);
        transpose(grad_ori, grad_ori_t);

        if(!(ime_t.isContinuous() && grad_ori_t.isContinuous())) {
            cout << "Matrices are not continuous, hence the Array struct will not work" << endl; 
            return -1;
        }

        arrayf E; E._x = ime_t.ptr<float>();
        arrayf O; O._x = grad_ori_t.ptr<float>();
        Size sz = ime.size();
        int h = sz.height; O._h=E._h=h;
        int w = sz.width; O._w=E._w=w;

        arrayf V;
        edgeBoxGen.generate( boxes, E, O, V );

        // create output bbs
        int n = (int) boxes.size();
        //cout << "Found " << n << " boxes" << endl;
        float *bbs = new float[5 * n];
        for(int i=0; i<n; i++) {
            bbs[ i + 0*n ] = (float) boxes[i].c+1;
            bbs[ i + 1*n ] = (float) boxes[i].r+1;
            bbs[ i + 2*n ] = (float) boxes[i].w;
            bbs[ i + 3*n ] = (float) boxes[i].h;
            bbs[ i + 4*n ] = boxes[i].s;
        }

        // show the bbs
        int n_show = 25;
        Mat im_show = im.clone();
        for(int i = 0; i < n_show; i++) {
            Point p1(int(bbs[i+0*n]), int(bbs[i+1*n])), p2(int(bbs[i+0*n] + bbs[i+2*n]), int(bbs[i+1*n] + bbs[i+3*n]));
            Scalar color(rand()&255, rand()&255, rand()&255);
            rectangle(im_show, p1, p2, color, 2);
        }
        filename = path + string("../results_cpp/") +  i->path().filename().string();
        imwrite(filename, im_show);
        cout << "Finished image " << count << endl;

        delete []bbs;
    }

    return 0;
}
