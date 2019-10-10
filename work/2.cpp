#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <omp.h>
#include <sys/time.h>

using namespace std;

struct distcos{
    double dist;
    int who;
};

int maxx(int x,int y){
    if (x>y) return x;
    else return y;
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
    //string ssl="bottleneck/l_s/";
    //string ssr="bottleneck/r_s/";
    string sm=".jpg.txt";
    sm.insert(0,".");
    sm.insert(0,to_string(d));
    //string sl="l.jpeg.txt";
    //string sr="r.jpeg.txt";
    string sf;
    //string sfl,sfr;
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
    
    /*
    if (d=='l'){
        sfl=sl;
        sfl.insert(0,to_string(i+1));
        sfl.insert(0,ssl);
        ifstream in(sfl);
        
        for (j=0;j<2048;j++){
            in>>t>>c;
            n.push_back(t);
        }
        
        in.close();
    }else {
        sfr=sr;
        sfr.insert(0,to_string(i+1));
        sfr.insert(0,ssr);
        ifstream in(sfr);
            
        for (k=0;k<2048;k++){
            in>>t>>c;
            n.push_back(t);
        }
        in.close();
    }
    */
    
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
    
    
    /*
    string ssl="data/l_s/";
    string ssr="data/r_s/";
    string sl="l.jpeg.txt";
    string sr="r.jpeg.txt";
    string sfl,sfr;
    
    {
    sfl=sl;
    sfl.insert(0,to_string(i+1));
    sfl.insert(0,ssl);
    ifstream iin(sfl);
        
    iin>>t1>>t2;
    wl.push_back(t1);
    hl.push_back(t2);
        
    iin.close();
    }
      
    {  
    sfr=sr;
    sfr.insert(0,to_string(i+1));
    sfr.insert(0,ssr);
    ifstream iin(sfr);
        
    iin>>t1>>t2;
    wr.push_back(t1);
    hr.push_back(t2);
    
    //cout<<t1<<' '<<t2<<endl;
    
    iin.close();
    }
    */
    
}

double shapesim(int wl,int hl,int wr,int hr){
    double sh,ww,hh;
    
    ww = abs(wl-wr)/maxx(wl,wr);
    hh = abs(hl-hr)/maxx(hl,hr);
    
    sh = exp((ww+hh)/2);
    
    return sh;
}

double find_sim(int l, int r){

    int i,j,k;
    distcos fl[110],fr[110];
    vector<double> ln,rn;
    double sim;

    for (i=0;i<100;i++)
        fr[i].dist=1;
    
    for (i=0;i<100;i++){
        
        ln=init_read(i,l);
        
        fl[i].dist=1;
        
        for (j=0;j<100;j++){
            
            rn=init_read(j,r);
            
            sim=1-cosine(ln,rn);
            
            if (sim<fl[i].dist) {
                fl[i].dist=sim;
                fl[i].who=j;
            }
            
            if (sim<fr[j].dist) {
                fr[j].dist=sim;
                fr[j].who=i;
            }
        }
    }
    
    vector<int> wl,wr,hl,hr;
    

    init_read_wh(l,wl,hl);
    init_read_wh(r,wr,hr);
    
    double fsim=0.0;
    
    for (i=0;i<100;i++){
        int wh = fl[i].who;
        //cout<<wh<<':'<<fl[i].dist<<endl;
        if (fr[wh].who == i){
            double shSim=shapesim(wl[i],hl[i],wr[wh],hr[wh]);
            //cout<<i<<','<<wh<<":"<<fl[i].dist<<endl;
            fsim+=1-shSim*fl[i].dist;
        }
    }
    
    fsim = fsim/100;
    
    return fsim;
    
}

int main(){
    
    int n_test=1;
    
    for (n_test=1;n_test<290;n_test++){
    
    int i,max_n;
    double max_sim=0.0;
    struct timeval tt_b, tt_e;
    
    gettimeofday(&tt_b, NULL);
    
        #pragma omp parallel for
        for (i=1;i<=289;i++){
        
            struct timeval t_b, t_e;
            //gettimeofday(&t_b, NULL);
        
            if (i==n_test) continue;
            double fs=find_sim(n_test,i);
        
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
    
    //double fs=find_sim(1,2);
    //cout<<fs<<endl;
    }
  
}
