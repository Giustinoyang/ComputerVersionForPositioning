#include <opencv2/opencv.hpp>
#define PI 3.14159265f

cv::Mat edge_nms(cv::Mat &E, cv::Mat &O, int r, int s, float m, int nThreads);
void edge_detect(const cv::Mat &im, cv::Mat &E, cv::Mat &O, std::string st_path);

// util functions
void vis_matrix(cv::Mat &m, char *window_name);
void gradient(cv::Mat &I, cv::Mat &gx, cv::Mat &gy);
cv::Mat signum(cv::Mat &src);
