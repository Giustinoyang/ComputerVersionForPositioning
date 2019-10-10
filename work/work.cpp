#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdlib.h> 
#include <omp.h>
#include <sys/time.h>

using namespace std;

//设置测试图片数量。
#define NumOfTest 254   
//设置特征向量维度
#define NumOfVec 2048
//设置建图时，每个点的近邻最大数量
#define MostNN  30
//设定可相信的最小平均相似度
#define rangeOfError 0.83

//存储余弦距离
struct distcos{
    double dist; //存储距离
    int who;    //存储谁的距离
};

//建图的点
struct typeOfGraph{
    vector<int> son;    //存储点在图中的邻近点集
    int level;          //存储点在跳表中的层级
    vector<double> vec; //存储点的特征向量
};

//比大小
int minn(int x,int y){
    if (x>y) return y;
    else return x;
}
int maxx(int x,int y){
    if (x>y) return x;
    else return y;
}

//判定宽搜队列前size个有没有变化
bool find_same(vector<distcos> l, vector<distcos> r, int size){
    int i;
    for (i=0;i<size;i++){
        if(l[i].who!=r[i].who) return false;
    }
    
    return true;
}

//求余弦相似度
double cosine(vector<double>ln,vector<double>rn){

    int i,j;
    double coss=0.0;
    double sl=0.0;
    double sr=0.0;
    
    for (i=0;i<NumOfVec;i++){
        coss+=ln[i]*rn[i];
        sl+=ln[i]*ln[i];
        sr+=rn[i]*rn[i];
    }
    
    coss = coss/(sqrt(sl)*sqrt(sr));
    
    return coss;
}

//读取d图的第i个子图的特征向量
vector<double> init_read(int i, int d){

    
    string ss="/";
    ss.insert(0,to_string(d));
    ss.insert(0,"bottleneck/");
    string sm=".jpg.txt";
    sm.insert(0,".");
    sm.insert(0,to_string(d));
    string sf;
    int j,k;
    vector<double> n;
    double t;
    char c;
    
    sf=sm;
    sf.insert(sf.length()-8,to_string(i+1));
    sf.insert(0,ss);
    {
        ifstream in(sf);
        for (j=0;j<NumOfVec;j++){
            in>>t>>c;
            n.push_back(t);
        }
        in.close();
    }

    return n;
    
}

//读取d的子图的所有w和h
void init_read_wh(int d, vector<int> &w, vector<int> &h){

    int t1,t2;
    
    string ss="/";
    ss.insert(0,to_string(d));
    ss.insert(0,"data/");
    string sm=".jpg.txt";
    sm.insert(0,".");
    sm.insert(0,to_string(d));
    string sf;
    
    int i,j;
    
    for (i=0;i<100;i++){
    
        sf=sm;
        sf.insert(sf.length()-8,to_string(i+1));
        sf.insert(0,ss);
        ifstream in(sf);
        
        in>>t1>>t2;
        w.push_back(t1);
        h.push_back(t2);
        
        in.close();
    }
    
}

//获取该点该处于跳表的层数，如下公式
int level_get(int m){

    double m1=(double)(rand()%100000)/100000;
    
    double lm=log(m);
    
    int ans=int(-1*log(m1)*lm);
        
    if (ans>7) ans=7;
    
    return ans;
}

//向量排序sort的参考函数
bool sortfun(const distcos x,const distcos y){
    return x.dist<=y.dist;
}

//寻找第k近邻
vector<distcos> find_kNN(int k, vector<double>ln, vector<typeOfGraph> gr){

    bool b[110];//每个点只加入一次队列
    int i,j;
    for (i=0;i<110;i++)
        b[i]=true;
    
    vector<distcos> vc(0), vf(0);//队列向量，和队列向量的镜像
    distcos dd; //中介变量
    /*
    for (i=0;i<6;i++){
        dd.who=i;
        dd.dist =1-cosine(ln,gr[i].vec);
        vc.push_back(dd);
        b[i]=false;
    }*/
    
    
    dd.who=0;
    dd.dist =1-cosine(ln,gr[0].vec);
    vc.push_back(dd);
    
    int layer=0;    //跳表结构接口，开始宽搜层级设置，最大为7最小为0,
    
    while (layer>=0){
        j=0;        //队头指针
        vf=vc;
        while (j<vc.size()){
            int p = vc[j].who;
            if (p>=gr.size()) { //测试中，有时下面的sort函数的边界会不清晰，导致无效内存参与排序，如果队头元素是无效元素，则删除。
                vc.erase(vc.begin()+j);
                continue;
            }
            for (i=0;i<gr[p].son.size();i++){
                if (gr[gr[p].son[i]].level<layer) continue; //跳表结构判定，只在更高层级搜索。
                if (b[gr[p].son[i]]==false) continue;       //每个元素只需要进入队列一次。
                if (k==1) b[gr[p].son[i]]=false;            //求解最近临时，就算元素无法进入队列，但也同时证明了该点不是最优解，所以提前减枝。k>1时，无法确定不是解集之一。
                double ddd=1-cosine(ln,gr[gr[p].son[i]].vec);//余弦距离
                if(ddd<vc[j].dist) {            //如果比队头元素更近
                    dd.who = gr[p].son[i];
                    dd.dist = ddd;
                    b[dd.who]=false;        //该点以加入队列标记
                    vc.push_back(dd);
                }
            }
            
            j++;
            
            if (j==vf.size()) {     //宽搜每层的长度等于vc队列的镜像vf的长度
                sort(vc.begin(), vc.end(), sortfun);//对队列中点集进行排序
                int v_size=minn(vc.size(), 6);      //本算法中，实验后，队列留6值较佳，但是前几轮搜索可能不足6个需要处理。
                vc.resize(v_size);                  //置vc为v_size长
                if (find_same(vf, vc, v_size)){     //如果队列无变化，证明已经是本层的较优解点集。
                    vc.resize(v_size);
                    break;
                }
                vf.resize(0);//初始化vf
                vf=vc;      //vc的镜像vf
                //vf.reserve(v_size);
                j=0;
            }
        }
        layer--;
    }
    vc.resize(k);//处理vc为查询长度k以供返回
    return vc;
}

//在G中插入一个点
void insert_G(int i,int d, vector<typeOfGraph> &gr){

    int j,k;
    vector<distcos> vc;
    vector<double> vi;
    typeOfGraph tg;
    
    vi=init_read(i,d);
    tg.vec=vi;
        
    int level=level_get(6); //给这个点一个level
    tg.level=level;
    
    vc=find_kNN(6,vi,gr);   //查询图中离本点最近的6个点
        
    for (j=0;j<vc.size();j++){  //连接本点和最近的6个点
        if (gr[vc[j].who].son.size()>MostNN) continue;//每个点最多有有限个近邻，防止图过于复杂，搜索时间退化。
        tg.son.push_back(vc[j].who);
        gr[vc[j].who].son.push_back(i);
    }
    
    gr.push_back(tg);   //把本点加入点集

}

//建图
void init_Graph(int d, vector<typeOfGraph> &gr){

    int i,j,k,level;
    typeOfGraph tg;

    for (i=0;i<6;i++){  //初始化将前六个点，互相相连，也可以在不起用跳表结构时，提供图中随机的“高速通道”，加快搜索速度。
        level = level_get(6);
        if (i==0) level = 7;
        tg.level = level;
        tg.vec = init_read(i,d);
        for (j=0;j<6;j++){
            if (i==j) continue;
            tg.son.push_back(j);
        }
        gr.push_back(tg);
    }
    
    for (i=6;i<100;i++){//剩余点插入
        
        insert_G(i,d,gr);
    }
}

//求形状相似度
double shapesim(int wl,int hl,int wr,int hr){//如论文，公式如下
    double sh,ww,hh;
    
    ww = abs(wl-wr)/maxx(wl,wr);
    hh = abs(hl-hr)/maxx(hl,hr);
    
    sh = exp((ww+hh)/2);
    
    return sh;
}

//求总相似度过程
double find_sim(int l, vector<typeOfGraph> gl, int r, vector<typeOfGraph> gr, int& num){
    
    int i,j,k;
    distcos fl[110], fr[110];
    vector<distcos> v;
    vector<double> ln, rn;
    double sim;
    
    for (i=0;i<100;i++){ //寻找每个子图的最近临子图
    
        ln=init_read(i,l);
        rn=init_read(i,r);
        v=find_kNN(1,ln, gr);
        fl[i]=v[0];
        v=find_kNN(1,rn, gl);
        fr[i]=v[0];
    }
    
    vector<int> wl,wr,hl,hr;
    
    //读取子图的宽高信息
    init_read_wh(l,wl,hl);
    init_read_wh(r,wr,hr);
    
    double fsim=0.0;
    
    for (i=0;i<100;i++){
        int wh = fl[i].who;
        if (fr[wh].who == i){//参与计算相似度的子图，需要双向验证，总相似度如论文，公式如下。
            double shSim=shapesim(wl[i],hl[i],wr[wh],hr[wh]);
            fsim+=1-shSim*fl[i].dist;
            num++;
        }
    }
    
    fsim = fsim/100;
    
    return fsim; //返回总相似度

}

int main(){
    
    int n_test=0;
    
    int wrong_n=0; //记录可能匹配图像非同一特征地点的数量

    vector<typeOfGraph> gr[NumOfTest+10]; 
    
    struct timeval t_b, t_e;
    gettimeofday(&t_b, NULL);
    
    #pragma omp parallel for 
    for (n_test=1;n_test<=NumOfTest;n_test++){//提前对图像库建图
        
        init_Graph(n_test, gr[n_test]);
         
    }
    
    gettimeofday(&t_e, NULL);
    int tt=(t_e.tv_sec-t_b.tv_sec)*1000+(t_e.tv_usec-t_b.tv_usec)/1000;
    cout<<"build graph cost: "<<tt<<"ms"<<endl;
    
    
    for (n_test=1;n_test<=NumOfTest;n_test++){//对图像库中每张图片自测
    
        int i,max_n,max_num;
        double max_sim=0.0;
        struct timeval tt_b, tt_e;
        
        gettimeofday(&tt_b, NULL);
        
        vector<typeOfGraph> gl(0);
        init_Graph(n_test, gl);

            #pragma omp parallel for num_threads(32)
            for (i=1;i<=NumOfTest;i++){ 
            
                //gettimeofday(&t_b, NULL);
            
                if (i==n_test) continue;
                int num=0;
                double fs=find_sim(n_test,gl,i,gr[i],num);
            
                if (fs>max_sim) {   //寻找图像库中最相似的图片
                    max_sim=fs;
                    max_n=i;
                    max_num=num;
                }
                
            
                //gettimeofday(&t_e, NULL);
                //int tt=(t_e.tv_sec-t_b.tv_sec)*1000+(t_e.tv_usec-t_b.tv_usec)/1000;
                //cout<<"computing with "<<i<<" cost: "<<tt<<"ms"<<", and the sim is : "<<fs<<endl;
            
            }
        
        gettimeofday(&tt_e, NULL);
        int ttt=(tt_e.tv_sec-tt_b.tv_sec)*1000+(tt_e.tv_usec-tt_b.tv_usec)/1000;
        
        double max_ans=max_sim*100/max_num;
        cout<<n_test<<"-"<<max_n<<": "<<max_sim<<"/"<<max_ans<<"  ";
        
        cout<<"Cost: "<<ttt<<"ms"<<" and the sim num: "<<max_num;
        
        if (max_ans<rangeOfError) {
            wrong_n++;
            cout<<" MAYBE WRONG";
        }
        
        cout<<endl;
    
    }
    
    double rate_w=(double)wrong_n/NumOfTest;
    cout<<wrong_n<<":"<<rate_w<<endl;

}
