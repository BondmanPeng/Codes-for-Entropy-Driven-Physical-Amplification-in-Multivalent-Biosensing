#include "system.h"

using namespace std;

double sampling_average_theta(double mu_g,double rho_r,double mu_linker,double nl,double kr,double kl,double xi_r,double xi_l,double fcnf,double area);
void sampling_theta(System& host_guest,vector<double>& data_theta,int& sampling_counts,int& repeat_counts);
double getAverage(vector<double>& sequence);

int main(int argc, char* argv[]){
    ofstream theta_data;
    theta_data.open("adsorption.txt");
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
    double lower_limit=-50-3*kl;
    double upper_limit=0-2*kl;
    double middle=(upper_limit+lower_limit)/2.;
    double delta=fabs(upper_limit-lower_limit);
    double theta;
    while(delta>0.01){
        middle=(lower_limit+upper_limit)/2.;
        theta=sampling_average_theta(mu_g,rho_r,middle,nl,kr,kl,xi_r,xi_l,fcnf,area);
        if(theta>0.1){
            upper_limit=middle;
        }else{
            lower_limit=middle;
        }
        delta=fabs(upper_limit-lower_limit);
    }
    theta_data<<(lower_limit+upper_limit)/2.<<endl;
}

double sampling_average_theta(double mu_g,double rho_r,double mu_linker,double nl,double kr,double kl,double xi_r,double xi_l,double fcnf,double area){
    vector<double> data_guest={};
    int sampling_counts=0,repeat_counts=0;
        // double mu_linker=-25+0.5*q;
        // double xi=exp(5); double xi_l=sqrt(xi); double xi_r=sqrt(xi);
    System host_guest;
    host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);
    host_guest.initialize_interactions();
    cout<<"generated_density: "<<host_guest.rec_substrate.Nr/area<<endl;
    int iteration=100000000;
    for(int i=0;i<iteration;++i){
        // if(i!=0 && i%30000000==0){
        //     int sampling_counts=0,repeat_counts=0;
        //         // System host_guest;
        //     host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);   
        //     // host_guest.initialize_interactions();     
        //     cout<<"generated_density: "<<host_guest.rec_substrate.Nr/area<<endl;
        // }
        host_guest.MCmove();
        sampling_theta(host_guest,data_guest,sampling_counts,repeat_counts);
            // cout<<host_guest.Ng<<endl;
            
    }
    cout<<"one finished"<<endl;
    return getAverage(data_guest)/host_guest.Nmax;
}

void sampling_theta(System& host_guest, vector<double>& data_guest,int& sampling_counts,int& repeat_counts){
    sampling_counts++;
    repeat_counts++;
    if(sampling_counts>20000000){
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