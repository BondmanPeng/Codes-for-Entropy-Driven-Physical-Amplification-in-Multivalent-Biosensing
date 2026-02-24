#include "system.h"

using namespace std;

void sampling_theta(System& host_guest,vector<double>& data_theta,int& sampling_counts,int& repeat_counts);
double getAverage(vector<double>& sequence);

//simulation code for calculating theta as a function of rho_r, with fixed mu_g, kl, kr, nl, area, xi. The results are written in theta.txt. The other output files are for testing the code and will not be used in the paper.

int main(int argc, char* argv[]){
    vector<double> data_guest;
    ofstream theta_data,alpha_data;
    theta_data.open("theta.txt");
    theta_data.precision(20);
    double mu_g=stod(argv[1]);
    // double mu_g=log(pow(10,-5)/(M_PI/4.));
    double rho_r=stod(argv[2]);
    // double rho_r=12.73; 
    double kl=stod(argv[3]);
    double kr=stod(argv[4]);
    // double kl=3;
    // double kr=3;
    double nl=stod(argv[5]);
    // double nl=10;
    double area=stod(argv[6]);
    // double area=150;
    double xi=stod(argv[7]);
    // double xi=exp(-2);
    double xi_l=sqrt(xi); double xi_r=sqrt(xi);
    double fcnf=2;
    for(int q=0;q<25;q++){
        data_guest={};
        int sampling_counts=0,repeat_counts=0;
        // double mu_linker=-30-2.5*kl+0.5*q; 
        double mu_linker=-40+2*q;
        // double xi=exp(5); double xi_l=sqrt(xi); double xi_r=sqrt(xi);
        System host_guest;
        host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);
        host_guest.initialize_interactions();
        cout<<"generated_density: "<<host_guest.rec_substrate.Nr/area<<endl;
        int iteration=250000000;
        for(int i=0;i<iteration;++i){
            if(i!=0 && i%10000000==0){
                int sampling_counts=0,repeat_counts=0;
                // System host_guest;
                host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);   
                host_guest.initialize_interactions();     
                cout<<"generated_density: "<<host_guest.rec_substrate.Nr/area<<endl;
            }
            host_guest.MCmove();
            sampling_theta(host_guest,data_guest,sampling_counts,repeat_counts);
            // cout<<host_guest.Ng<<endl;
            
        }
        cout<<"one finished"<<endl;
        theta_data<<getAverage(data_guest)/host_guest.Nmax<<endl;
    }
}

void sampling_theta(System& host_guest, vector<double>& data_guest,int& sampling_counts,int& repeat_counts){
    sampling_counts++;
    repeat_counts++;
    if(sampling_counts>5000000){
        if(repeat_counts>5){
            data_guest.push_back((double)host_guest.Ng);
            repeat_counts=0;
            // cout<<"step_size: "<<host_guest.ds<<endl;
        }
        // if(sampling_counts%5000000==0){
        //     cout<<host_guest.Ng<<endl;
        // }
    }
}

double getAverage(vector<double>& sequence){ // calculate average value of a vector
    int length=sequence.size();
    double sum=0,average;
    for(int i =0; i<length; i++){
        sum=sum+(double) sequence[i];
    }
    average=sum/length;
    return average;

}