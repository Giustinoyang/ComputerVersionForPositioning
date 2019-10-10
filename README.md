# ComputerVersionForPositioning
 2019 summer in Beijing


工具准备
使用edge_boxes初始化图像数据，获取特征子图。
编译edge_boxes需要numpy-opencv-converter，cvmatio工具。
编译cvmatio工具时，需要注释掉/scmatio/src/MatlabIO.cpp:523,526,527行。
疑似ubuntu64位系统兼容问题导致uint32_t,int64_t,uint64_t报错，修改edge_boxes源文件时注意避开相关数据类型。
编译上述工具需要boost1.71X，opencv2.x，gcc4.X，g++4.x，cudnn8.x等旧版本工具。
需要手动编译安装gflags，glog工具，ubuntu自动获取的版本疑似兼容问题。
以上大部分工具通过cmake编译安装时，需要加入c++11支持。
使用tensorflow获取特征子图的特征向量。
按官网教程安装tensorflow。
下载准备Inception-v3模型。

基本思路
在/edges/cpp/build/edge_data下放入数据图片，与rename.sh文件
运行rename.sh
在/edges/cpp/build下运行make_boxes.sh
源代码是edges/cpp/src/edge_boxes_demo.cpp 在edges/cpp/build下make编译

移动/edges/cpp/build/edge_data里生成的文件夹到/work/data
移除/work/bottleneck中所有文件
运行/work下1.py生成特征子图的特征向量
打开work.cpp修改NumOfTest的值为文件夹数量
执行
g++ work.cpp -o work -std=c++11 -fopenmp
编译文件
运行./work
命令行输出，每张图片在数据集中除去本身最相似的图像，总相似度，平均子相似度，计算时间，子相似度。

2.cpp为暴力解，3.cpp为宽搜解，4.cpp为预建图宽搜解，work.cpp为与建图宽搜减枝并且留有跳表结构的解
实践表明，100个特征子图的特征向量建图，使用跳表结构时间优化有限，留下结构，以应对更多的特征子图求解需求。


以相似度的和除以总子图数量，为判定图片之间相似度的标准，论文如是说，实践可行。
以每张互相最近临的子图的平均相似度，验证相似图片的准确性，我手动验证数据时，发现以子图平均相似度验证，可以基本确保准确率和召回率，标记部分图片集中没有的图片查询结果。
