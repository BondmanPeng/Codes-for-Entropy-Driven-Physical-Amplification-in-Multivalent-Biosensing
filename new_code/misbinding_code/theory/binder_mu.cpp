#include "binder_system.h"

using namespace std;

double poisson_average(double average_nr,vector<double> container);

int main(int argc, char* argv[]){
    ofstream theta_data,alpha_data,unbound_data;
    ofstream bivalency,multivalency;
    theta_data.open("theta.txt");
    alpha_data.open("alpha.txt");
    unbound_data.open("unbound_ligands.txt");
    bivalency.open("bivalency.txt");
    multivalency.open("multivalency.txt");

    // bridge_bonding_data.open("bonding_valency.txt");
    // absorb_data.open("absorb_valency.txt");
    // unbound_data.open("unbound_ligands.txt");
    // product_data.open("product.txt");
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
        vector<double> bivalences, multivalences;
        vector<double> unbonds;
        for(int j=1;j<120;j++){
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

            double mu_m=stod(argv[7]);
            // double mu_lm=mu_m;
            // double mu_rm=mu_m;
            double xim=exp(-stod(argv[8]));
            double nrp=j;
            // double nr=10;
            double mu_linker=-40+0.1*i;
            
            double fcnf=2;
            double zg=pow(10,-5);
            double beta=1;
            
            host_guest.initialize_system(xi_l,xim,kl,nl,xi_r,xim,kr,nrp,mu_linker,mu_m,fcnf,zg,beta);
            host_guest.picard_solve_system();
            host_guest.picard_solve_ref_system();   
            host_guest.cal_Feff();
            host_guest.cal_q();
            host_guest.cal_dFeffdx();
            host_guest.cal_dqdx();
            host_guest.cal_dlnthetadx();
            host_guest.cal_alpha();
            
            double block_l=host_guest.guest_part.bar_nl*host_guest.guest_part.xi_l;
            double block_r=host_guest.substrate_plate.bar_nr*host_guest.substrate_plate.xi_r;
            double sumkml=exp(mu_linker)*host_guest.guest_part.kl*block_l*pow(block_l+1.,host_guest.guest_part.kl-1.);
            double sumkql=exp(mu_linker-fcnf)*host_guest.guest_part.kl*block_l*pow(block_l+1.,host_guest.guest_part.kl-1.)*(pow(block_r+1.,host_guest.substrate_plate.kr)-1.);
            double Eal=host_guest.guest_part.bar_nl*exp(mu_m)*host_guest.guest_part.xi_lm;
            double Eb= host_guest.guest_part.bar_nl*host_guest.substrate_plate.bar_nr*exp(mu_m-fcnf)*host_guest.guest_part.xi_lm*host_guest.substrate_plate.xi_rm;
            thetas.push_back(host_guest.q*host_guest.zg/(1.+host_guest.q*host_guest.zg));
            // bivalences.push_back(Eal+Eb);
            // multivalences.push_back(sumkml+sumkql);
            // unbonds.push_back(host_guest.guest_part.bar_nl);
            
            alphas.push_back(host_guest.alpha);
            
        }
        double avg_theta,avg_alpha;
        avg_theta=poisson_average(nr,thetas);
        avg_alpha=poisson_average(nr,alphas);
        theta_data<<avg_theta<<endl;
        alpha_data<<avg_alpha<<endl;
        double avg_bivalency,avg_multivalency,avg_unbonds;
        avg_bivalency=poisson_average(nr,bivalences);
        avg_multivalency=poisson_average(nr,multivalences);
        avg_unbonds=poisson_average(nr,unbonds);
        bivalency<<avg_bivalency<<endl;
        multivalency<<avg_multivalency<<endl;
        unbound_data<<avg_unbonds<<endl;


    }
    // cout<<"xi_l: "<<xi_l<<endl;
    // cout<<"kl: "<<kl<<endl;
    // cout<<"nl: "<<nl<<endl;
    // cout<<"xi_r: "<<xi_r<<endl;
    // cout<<"kr: "<<kr<<endl;
    // cout<<"nr: "<<nr<<endl;
    // cout<<"mu_m: "<<stod(argv[7])<<endl;
    // cout<<"xim: "<<exp(-stod(argv[8]))<<endl;
    
}

double poisson_average(double average_nr,vector<double> container){
    int length=container.size();
    double result=0;
    for(int i=1;i<=length;i++){
        result+=container[i-1]*pow(average_nr,i)*exp(-average_nr)/tgamma(i+1);
    }
    return result;
}
