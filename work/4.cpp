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

struct distcos{
    double dist;
    int who;
};

struct typeOfGraph{
    vector<int> son;
    int level;
    vector<double> vec;
};

int maxx(int x,int y){
    if (x>y) return x;
    else return y;
}

bool find_same(vector<distcos> l, vector<distcos> r){
    int i;
    for (i=0;i<6;i++){
        if(l[i].who!=r[i].who) return false;
    }
    
    return true;
}

double cosine(vector<double>ln,vector<double>rn){

    int i,j;
    double coss=0.0;
    double sl=0.0;
    double sr=0.0;
    
    for (i=0;i<2048;i++){
        coss+=ln[i]*rn[i];
        sl+=ln[i]*ln[i];
        sr+=rn[i]*rn[i];
    }
    
    coss = coss/(sqrt(sl)*sqrt(sr));
    
    return coss;
}

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
        for (j=0;j<2048;j++){
            in>>t>>c;
            n.push_back(t);
        }
        in.close();
    }

    return n;
    
}

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

int level_get(int m){

    double m1=(double)(rand()%100000)/100000;
    
    double lm=log(m);
    
    int ans=int(-1*log(m1)*lm);
        
    if (ans>7) ans=7;
    
    return ans;
}

bool sortfun(const distcos x,const distcos y){
    return x.dist<=y.dist;
}

vector<distcos> find_kNN(int k, vector<double>ln, vector<typeOfGraph> gr){

    bool b[110];
    int i,j;
    for (i=0;i<110;i++)
        b[i]=true;
    
    vector<distcos> vc(0), vf(0);
    distcos dd;
    
    for (i=0;i<6;i++){
        dd.who=i;
        dd.dist =1-cosine(ln,gr[i].vec);
        vc.push_back(dd);
        b[i]=false;
    }
    
    j=0;
    vf=vc;
    
    while (j<vc.size()){
        int p = vc[j].who;
        if (p>=gr.size()) {
            vc.erase(vc.begin());
            j++;
            continue;
        }
        for (i=0;i<gr[p].son.size();i++){
            if (b[gr[p].son[i]]==false) continue;
            double ddd=1-cosine(ln,gr[gr[p].son[i]].vec);
            if(ddd<vc[j].dist) {
                dd.who = gr[p].son[i];
                dd.dist = ddd;
                b[dd.who]=false;
                vc.push_back(dd);
            }
        }
        j++;
        if (j==6) {
            sort(vc.begin(), vc.end(), sortfun);
            vc.reserve(6);
            if (find_same(vf, vc)){
                vc.resize(k);
                return vc;
            }
            vf.resize(0);
            vf=vc;
            j=0;
        }
    }
    cout<<"some problem"<<endl;
    return vc;
}

void insert_G(int i,int d, vector<typeOfGraph> &gr){

    int j,k;
    vector<distcos> vc;
    vector<double> vi;
    typeOfGraph tg;
    
    vi=init_read(i,d);
    tg.vec=vi;
        
    int level=level_get(6);
    tg.level=level;
    
    vc=find_kNN(6,vi,gr);
        
    for (j=0;j<vc.size();j++){
        if (gr[vc[j].who].son.size()>20) continue;
        tg.son.push_back(vc[j].who);
        gr[vc[j].who].son.push_back(i);
    }
    
    gr.push_back(tg);

}

void init_Graph(int d, vector<typeOfGraph> &gr){

    int i,j,k,level;
    typeOfGraph tg;
    
    for (i=0;i<6;i++){
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
    
    for (i=6;i<100;i++){
        
        insert_G(i,d,gr);
    }
}

double shapesim(int wl,int hl,int wr,int hr){
    double sh,ww,hh;
    
    ww = abs(wl-wr)/maxx(wl,wr);
    hh = abs(hl-hr)/maxx(hl,hr);
    
    sh = exp((ww+hh)/2);
    
    return sh;
}

double find_sim(int l, vector<typeOfGraph> gl,int r, vector<typeOfGraph> gr){
    
    int i,j,k;
    distcos fl[110], fr[110];
    vector<distcos> v;
    vector<double> ln, rn;
    double sim;
    

    //init_Graph(l, gl);
    //init_Graph(r, gr);
    
    for (i=0;i<100;i++){
    
        ln=init_read(i,l);
        rn=init_read(i,r);
        v=find_kNN(1,ln, gr);
        fl[i]=v[0];
        v=find_kNN(1,rn, gl);
        fr[i]=v[0];
    }
    
    vector<int> wl,wr,hl,hr;
    

    init_read_wh(l,wl,hl);
    init_read_wh(r,wr,hr);
    
    double fsim=0.0;
    
    for (i=0;i<100;i++){
        int wh = fl[i].who;
        if (fr[wh].who == i){
            double shSim=shapesim(wl[i],hl[i],wr[wh],hr[wh]);
            fsim+=1-shSim*fl[i].dist;
        }
    }
    
    fsim = fsim/100;
    
    return fsim;

}

int main(){
    
    int n_test=103;
    
    vector<typeOfGraph> gr[300];
    
    struct timeval t_b, t_e;
    gettimeofday(&t_b, NULL);
    
    #pragma omp parallel for 
    for (n_test=1;n_test<=289;n_test++){
        
        init_Graph(n_test, gr[n_test]);
         
    }
    
    gettimeofday(&t_e, NULL);
    int tt=(t_e.tv_sec-t_b.tv_sec)*1000+(t_e.tv_usec-t_b.tv_usec)/1000;
    cout<<"build graph cost: "<<tt<<"ms"<<endl;
    
    
    for (n_test=1;n_test<=289;n_test++){
    
    int i,max_n;
    double max_sim=0.0;
    struct timeval tt_b, tt_e;
    
    gettimeofday(&tt_b, NULL);
    
    vector<typeOfGraph> gl;
    init_Graph(n_test, gl);

        #pragma omp parallel for
        for (i=1;i<=289;i++){
        
            //struct timeval t_b, t_e;
            //gettimeofday(&t_b, NULL);
        
            if (i==n_test) continue;
            double fs=find_sim(n_test,gl,i,gr[i]);
        
            if (fs>max_sim) {
                max_sim=fs;
                max_n=i;
            }
        
            //gettimeofday(&t_e, NULL);
            //int tt=(t_e.tv_sec-t_b.tv_sec)*1000+(t_e.tv_usec-t_b.tv_usec)/1000;
            //cout<<"computing with "<<i<<" cost: "<<tt<<"ms"<<", and the sim is : "<<fs<<endl;
        
        }
    

    
    gettimeofday(&tt_e, NULL);
    int ttt=(tt_e.tv_sec-tt_b.tv_sec)*1000+(tt_e.tv_usec-tt_b.tv_usec)/1000;
    
    cout<<n_test<<"-"<<max_n<<": "<<max_sim<<"  ";
    
    cout<<"Cost: "<<ttt<<"ms"<<endl;
    
    }

}
