#include "system.h"

using namespace std;

double poisson_average(double average_nr,vector<double> container);// a function for calculating the average of a quantity over different valencies of the substrate plate, which follows a Poisson distribution with average nr. The input container is the quantity calculated for different valencies of the substrate plate, and the output is the average of this quantity over the Poisson distribution. This function will be used in the main function to calculate the average theta and alpha over different valencies of the substrate plate.

// codes for calculating theoretical results for theta and alpha as a function of mu_linker, with fixed kl, nl, kr, nr. The results are averaged over different valencies of the substrate plate, which follows a Poisson distribution with average nr. The results are written in theta.txt and alpha.txt, respectively. The other output files are for testing the code and will not be used in the paper.

int main(int argc, char* argv[]){
    ofstream theta_data,alpha_data,bridge_data,bridge_bonding_data,absorb_data,unbound_data,product_data;
    theta_data.open("theta.txt");
    alpha_data.open("alpha.txt");
    bridge_data.open("bridge_valency.txt");
    bridge_bonding_data.open("bonding_valency.txt");
    absorb_data.open("absorb_valency.txt");
    unbound_data.open("unbound_ligands.txt");
    product_data.open("product.txt");
    theta_data.precision(20);
    alpha_data.precision(20);
    double xi_l=exp(-stod(argv[1]));
            // double xi_l=exp(-0);
    double kl=stod(argv[2]);
            // double kl=3;
    double nl=stod(argv[3]);
            // double nl=10;
    double xi_r=exp(-stod(argv[4]));
            // double xi_r=exp(-0);
    double kr=stod(argv[5]);
    double nr=stod(argv[6]);
    for(int i=0;i<500;i++){
        vector<double> thetas, alphas;
        for(int j=1;j<35;j++){
            System host_guest;
            double xi_l=exp(-stod(argv[1]));
            // double xi_l=exp(-0);
            double kl=stod(argv[2]);
            // double kl=3;
            double nl=stod(argv[3]);
            // double nl=10;
            double xi_r=exp(-stod(argv[4]));
            // double xi_r=exp(-0);
            double kr=stod(argv[5]);
            // double kr=3;
            // double nr=stod(argv[6]);
            double nrp=j;
            // double nr=10;
            double mu_linker=-40+0.1*i;// the chemical potential of the linker, which is varied to calculate theta and alpha as a function of mu_linker. The range of mu_linker is from -40 to 10, with a step of 0.1. The value of mu_linker is chosen to cover the range of theta from 0 to 1, and the range of alpha from 0 to 1. The results will be used to plot theta and alpha as a function of mu_linker in the paper.
            
            double fcnf=2;
            double zg=pow(10,-5);
            double beta=1;
            host_guest.initialize_system(xi_l,kl,nl,xi_r,kr,nr,mu_linker,fcnf,zg,beta);
            host_guest.picard_solve_system();
            host_guest.picard_solve_ref_system();   
            host_guest.cal_Feff();
            host_guest.cal_q();
            host_guest.cal_dFeffdx();
            host_guest.cal_dqdx();
            host_guest.cal_dlnthetadx();
            host_guest.cal_alpha();
            // double block_l=host_guest.guest_part.bar_nl*host_guest.guest_part.xi_l;
            // double block_r=host_guest.substrate_plate.bar_nr*host_guest.substrate_plate.xi_r;
            // double block_l2=pow(host_guest.guest_part.bar_nl*host_guest.guest_part.xi_l+1,host_guest.guest_part.kl);
            // double block_r2=pow(host_guest.substrate_plate.bar_nr*host_guest.substrate_plate.xi_r+1,host_guest.substrate_plate.kr);
            // double block_l_ref=host_guest.substrate_plate.ref_bar_nr*host_guest.substrate_plate.xi_r;
            // double sumkql=exp(mu_linker-fcnf)*host_guest.guest_part.kl*block_l*pow(block_l+1.,host_guest.guest_part.kl-1.)*(pow(block_r+1.,host_guest.substrate_plate.kr)-1.);
            
            // double sumb=exp(mu_linker-fcnf)*(block_l2*block_r2-block_l2-block_r2+1);
            // double sumkml=exp(mu_linker)*host_guest.guest_part.kl*block_l_ref*pow(block_l_ref+1.,host_guest.guest_part.kl-1.);
            
            // double summl=exp(mu_linker)*(block_l_ref-1);
            thetas.push_back(host_guest.q*host_guest.zg/(1.+host_guest.q*host_guest.zg));
            alphas.push_back(host_guest.alpha); 
        }
        double avg_theta,avg_alpha;
        avg_theta=poisson_average(nr,thetas);
        avg_alpha=poisson_average(nr,alphas);
        theta_data<<avg_theta<<endl;
        alpha_data<<avg_alpha<<endl;
        // bridge_data<<sumb<<endl;
        // bridge_bonding_data<<sumkql/sumb<<endl;
        // absorb_data<<sumkml<<endl;
        // unbound_data<<host_guest.guest_part.bar_nl<<endl;
        // product_data<<host_guest.guest_part.bar_nl*exp(-f_l)<<endl;
    }
    
}

double poisson_average(double average_nr,vector<double> container){
    int length=container.size();
    double result=0;
    for(int i=1;i<=length;i++){
        result+=container[i-1]*pow(average_nr,i)*exp(-average_nr)/tgamma(i+1);
    }
    return result;
}