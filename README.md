# ComputerVersionForPositioning
 2019 summer in Beijing


����׼��
ʹ��edge_boxes��ʼ��ͼ�����ݣ���ȡ������ͼ��
����edge_boxes��Ҫnumpy-opencv-converter��cvmatio���ߡ�
����cvmatio����ʱ����Ҫע�͵�/scmatio/src/MatlabIO.cpp:523,526,527�С�
����ubuntu64λϵͳ�������⵼��uint32_t,int64_t,uint64_t�����޸�edge_boxesԴ�ļ�ʱע��ܿ�����������͡�
��������������Ҫboost1.71X��opencv2.x��gcc4.X��g++4.x��cudnn8.x�Ⱦɰ汾���ߡ�
��Ҫ�ֶ����밲װgflags��glog���ߣ�ubuntu�Զ���ȡ�İ汾���Ƽ������⡣
���ϴ󲿷ֹ���ͨ��cmake���밲װʱ����Ҫ����c++11֧�֡�
ʹ��tensorflow��ȡ������ͼ������������
�������̳̰�װtensorflow��
����׼��Inception-v3ģ�͡�

����˼·
��/edges/cpp/build/edge_data�·�������ͼƬ����rename.sh�ļ�
����rename.sh
��/edges/cpp/build������make_boxes.sh
Դ������edges/cpp/src/edge_boxes_demo.cpp ��edges/cpp/build��make����

�ƶ�/edges/cpp/build/edge_data�����ɵ��ļ��е�/work/data
�Ƴ�/work/bottleneck�������ļ�
����/work��1.py����������ͼ����������
��work.cpp�޸�NumOfTest��ֵΪ�ļ�������
ִ��
g++ work.cpp -o work -std=c++11 -fopenmp
�����ļ�
����./work
�����������ÿ��ͼƬ�����ݼ��г�ȥ���������Ƶ�ͼ�������ƶȣ�ƽ�������ƶȣ�����ʱ�䣬�����ƶȡ�

2.cppΪ�����⣬3.cppΪ���ѽ⣬4.cppΪԤ��ͼ���ѽ⣬work.cppΪ�뽨ͼ���Ѽ�֦������������ṹ�Ľ�
ʵ��������100��������ͼ������������ͼ��ʹ������ṹʱ���Ż����ޣ����½ṹ����Ӧ�Ը����������ͼ�������


�����ƶȵĺͳ�������ͼ������Ϊ�ж�ͼƬ֮�����ƶȵı�׼����������˵��ʵ�����С�
��ÿ�Ż�������ٵ���ͼ��ƽ�����ƶȣ���֤����ͼƬ��׼ȷ�ԣ����ֶ���֤����ʱ����������ͼƽ�����ƶ���֤�����Ի���ȷ��׼ȷ�ʺ��ٻ��ʣ���ǲ���ͼƬ����û�е�ͼƬ��ѯ�����
