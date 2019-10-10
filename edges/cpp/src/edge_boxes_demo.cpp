#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <algorithm>
#include <string>

#include "edge_detect.h"
#include "edge_boxes.h"

using namespace std;
using namespace cv;

void get_edge_boxes(Mat &im, vector<vector<float> > &bbs) {
    
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
    edge_detect(im, ime, grad_ori, string("/home/xcy/edges/cpp/external/gop_1.3/data/sf.dat"));
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
    cout << "Generated boxes, t = " << t*1000 << " ms" << endl;

    // create output bbs
    int n = (int) boxes.size();
    //cout << "Found " << n << " boxes" << endl;
    bbs.resize(n, vector<float>(5, 0));
    for(int i=0; i<n; i++) {
        bbs[i][0] = (float) boxes[i].c+1;
        bbs[i][1] = (float) boxes[i].r+1;
        bbs[i][2] = (float) boxes[i].w;
        bbs[i][3] = (float) boxes[i].h;
        bbs[i][4] = boxes[i].s;
    }
}

int main(int argc, char **argv) {
    if(argc != 3) {
        cout << "Usage: ./edge_boxes image_file num_boxes_to_show" << endl;
        return -1;
    }

    Mat im = imread(argv[1]);
    if(im.data == NULL) {
        cout << "Error reading image" << endl;
        return -1;
    }

    int n_show = stoi(argv[2], NULL);

    vector<vector<float> > bbs;
    get_edge_boxes(im, bbs);

    // show the bbs
    n_show = std::min(n_show, int(bbs.size()));
    Mat im_show = im.clone();
    
    int nn=4;
    
    for(int i = 0; i < n_show; i++) {
        Point p1(int(bbs[i][0]), int(bbs[i][1])), p2(int(bbs[i][0] + bbs[i][2]), int(bbs[i][1] + bbs[i][3]));
        //Scalar color(rand()&255, rand()&255, rand()&255);
        //rectangle(im_show, p1, p2, color, 2);
        
        Rect rect(int(bbs[i][0]), int(bbs[i][1]),int(bbs[i][2]), int(bbs[i][3]));
        Mat result_image(im_show, rect);
        
        string sss=argv[1];
        
        if (sss[sss.length()-5]=='.') nn=5;
        sss.insert(sss.length()-nn,".");
        sss.insert(sss.length()-nn,to_string(i+1));
        
        //imshow(sss, result_image);
        
        imwrite(sss,result_image);
        
        //output the w and h of the subpicture;
        
        string sout=".txt";
        sout.insert(0,sss);
        ofstream out(sout);
        out<<int(bbs[i][2])<<' '<<int(bbs[i][3]);
        
    }

    //imshow("Edge-Boxes", im_show);
    
    //while(waitKey(1) != 'q') {}
    return 0;
}
