#include "system.h"

using namespace std;
void one_run(int iteration,double mu_g,double rho_r,double mu_linker,double nl,double kr,double kl,double xi_r,double xi_l, double fcnf,double area, vector<double>& data_guest,ofstream& configuration);
void sampling_theta(System& host_guest,vector<double>& data_theta,int& sampling_counts,int& repeat_counts);
double getAverage(vector<double>& sequence);
void record_configuration(System & any_system, ofstream& configuration);

int main(int argc, char* argv[]){
    vector<double> data_guest;

    ofstream alpha_data,configuration;
    alpha_data.open("alpha.txt");
    configuration.open("configuration.txt");
    alpha_data.precision(20);
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
    double interval=0.25;
    // double prev_theta, next_theta;
    for(int q=0;q<25;q++){
        vector<double> thetas={};
        for(int p=0;p<4;p++){
            data_guest={};
            int iteration=360000000;
            
            // double mu_linker=-30-2.5*kl+0.5*q; 
            double mu_linker;
            if (p<2){
                mu_linker=-40+2*q+interval*(p-2);
            }else{
                mu_linker=-40+2*q+interval*(p-1);
            }
            one_run(iteration,mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l, fcnf, area, data_guest,configuration);
            // double xi=exp(5); double xi_l=sqrt(xi); double xi_r=sqrt(xi);
            
            cout<<"one finished"<<endl;
            double Nmax=area/(M_PI*pow(1./2.,2));
            thetas.push_back(getAverage(data_guest)/Nmax);
        }
        alpha_data<<(-log(thetas[3])+8.*log(thetas[2])-8.*log(thetas[1])+log(thetas[0]))/(12*interval)<<endl;
        
    }
}
void one_run(int iteration,double mu_g,double rho_r,double mu_linker,double nl,double kr,double kl,double xi_r,double xi_l, double fcnf,double area, vector<double>& data_guest,ofstream& configuration){
    System host_guest;
    host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);
    host_guest.initialize_interactions();
    // cout<<"generated_density: "<<host_guest.rec_substrate.Nr/area<<endl;
    
    int sampling_counts=0,repeat_counts=0;
    for(int i=0;i<iteration;++i){
        // if(i!=0 && i%60000000==0){
        
        // System host_guest;
        // host_guest.initialize_system(mu_g,rho_r,mu_linker,nl,kr,kl,xi_r,xi_l,fcnf,area);   
        // host_guest.initialize_interactions();     
        
        // }
        host_guest.MCmove();
        if(i%10000000==0){
            record_configuration(host_guest,configuration);
            cout<<host_guest.Ng<<endl;
        }
        sampling_theta(host_guest,data_guest,sampling_counts,repeat_counts);
        // cout<<host_guest.Ng<<endl;
    }
}

void sampling_theta(System& host_guest, vector<double>& data_guest,int& sampling_counts,int& repeat_counts){
    sampling_counts++;
    repeat_counts++;
    if(sampling_counts>100000000){
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

void record_configuration(System & any_system, ofstream& configuration){
    configuration<<any_system.Ng<<endl;
    configuration<<" "<<endl;
    for(Guest* iter=any_system.guest_particles.headPtr; iter != nullptr; iter=iter->list_next_pointer){
        configuration<<"H"<<" "<<iter->boxs[0]*any_system.L<<" "<<iter->boxs[1]*any_system.L<<endl;
    }   
    configuration<<" "<<endl;

}