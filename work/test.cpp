#include <iostream>
#include <cmath>
#include <stdlib.h> 
#include <vector>

using namespace std;

void work(vector<int> &f){

    int i;
    srand(time(0));
    
    int m=6;

    for (i=0;i<100;i++){
        
        double m1=(double)(rand()%100000)/100000;
    
        double lm=log(m);
    
        int ans=int(-1*log(m1)*lm);
        
        if (ans>7) ans=7;
        
        f[ans]++;
    }
}

int main(){
    
    vector<int> f;
    
    int i;
    
    for (i=0;i<8;i++)
        f.push_back(0);
    
    work(f);
    
    int ans=0;
    for (i=7;i>=0;i--){
        ans=ans+f[i];
        cout<<ans<<' ';
    }
    cout<<endl;
}
