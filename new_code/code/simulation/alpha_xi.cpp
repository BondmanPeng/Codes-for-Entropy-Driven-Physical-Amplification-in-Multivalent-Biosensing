#include "system.h"

using namespace std;

class parameters{
public:
    double mu_g,rho_r,kl,kr,nl,area,mu_linker,fcnf;
    double xi_l,xi_r;
};

void sampling_theta(System& host_guest,vector<double>& data_theta,int& sampling_counts,int& repeat_counts);
void sampling_MC(vector<double>& data_guest,parameters package);
double getAverage(vector<double>& sequence);



int main(int argc, char* argv[]){
    vector<double> data_guest;
    ofstream alpha_data;
    alpha_data.open("alpha.txt");
    alpha_data.precision(20);
    parameters package;
    package.mu_g=stod(argv[1]);
    // double mu_g=log(pow(10,-5)/(M_PI/4.));
    package.rho_r=stod(argv[2]);
    // double rho_r=12.73; 
    package.kl=stod(argv[3]);
    package.kr=stod(argv[4]);
    // double kl=3;
    // double kr=3;
    package.nl=stod(argv[5]);
    // double nl=10;
    package.area=stod(argv[6]);
    // double area=150;
    package.mu_linker=stod(argv[7]); 

    package.fcnf=2.;
    double interval=0.02;
    for(int q=0;q<30;q++){
        vector<double> thetas={};
        for(int p=0;p<2;p++){
            data_guest={};
            int sampling_counts=0,repeat_counts=0;
            double lnxi=(-15+q)+interval*(2.*p-1.);
            // double kf=-lnxik/2;
            package.xi_l=exp(lnxi/2.); package.xi_r=exp(lnxi/2.);
            // double xi=exp(5); double xi_l=sqrt(xi); double xi_r=sqrt(xi);
            
            int reshuffle=15;
            for(int i=0;i<reshuffle;++i){
                sampling_MC(data_guest,package);
                // cout<<host_guest.Ng<<endl;
            }
            cout<<"one finished"<<endl;
            double Nmax=package.area/(M_PI*pow(1./2.,2));
            thetas.push_back(getAverage(data_guest)/Nmax);
        }
        alpha_data<<-(log(thetas[1])-log(thetas[0]))/(2*interval)<<endl;
    }
}

void sampling_theta(System& host_guest, vector<double>& data_guest,int& sampling_counts,int& repeat_counts){
    sampling_counts++;
    repeat_counts++;
    if(sampling_counts>1000000){
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

void sampling_MC(vector<double>& data_guest,parameters package){
    int iteration=20000000;
    int sampling_counts=0,repeat_counts=0;
    System host_guest;
    host_guest.initialize_system(package.mu_g,package.rho_r,package.mu_linker,package.nl,package.kr,package.kl,package.xi_r,package.xi_l,package.fcnf,package.area);   
    host_guest.initialize_interactions(); 
    cout<<"generated_density: "<<host_guest.rec_substrate.Nr/package.area<<endl;  
    for(int i=0;i<iteration;++i){
        host_guest.MCmove();
        sampling_theta(host_guest,data_guest,sampling_counts,repeat_counts);
        // cout<<host_guest.Ng<<endl;
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