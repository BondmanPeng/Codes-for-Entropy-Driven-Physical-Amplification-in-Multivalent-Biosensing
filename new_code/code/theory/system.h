#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

class Guest{// the guest part of the system, which has nl ligands, and each ligand can form kl bonds with the substrate plate. The variable xi_l is the Boltzmann weight of forming a bond between a ligand and the substrate plate, which is related to the binding free energy of the bond. The variable bar_nl is the average number of ligands that are not bound to the substrate plate, which is calculated from the self-consistent equations. The variable dndx is the derivative of bar_nl with respect to mu_linker, which is calculated from the self-consistent equations. The variable ref_bar_nl and drefndx are the corresponding variables for the reference system, which is used to calculate the effective free energy of binding.
public:
    double xi_l,kl,nl,bar_nl,ref_bar_nl;
    double dndx,drefndx;
    void initialize_guest(double given_xi_l,double given_kl,double given_nl){
        xi_l=given_xi_l; kl=given_kl; nl=given_nl; bar_nl=nl/2.; ref_bar_nl=nl/2.;
        dndx=0; drefndx=0;
    }
};

class Substrate{// the substrate part of the system, which has nr receptors, and each receptor can form kr bonds with the ligands on the guest part. The variable xi_r is the Boltzmann weight of forming a bond between a receptor and the ligands on the guest part, which is related to the binding free energy of the bond. The variable bar_nr is the average number of receptors that are not bound to the ligands on the guest part, which is calculated from the self-consistent equations. The variable dndx is the derivative of bar_nr with respect to mu_linker, which is calculated from the self-consistent equations. The variable ref_bar_nr and drefndx are the corresponding variables for the reference system, which is used to calculate the effective free energy of binding.
public:
    double xi_r,kr, nr,bar_nr, ref_bar_nr;
    double dndx,drefndx;
    double nrp;
    void initialize_substrate(double given_xi_r,double given_kr, double given_nr){
        xi_r=given_xi_r; kr=given_kr; nr=given_nr; bar_nr=nr/2.; ref_bar_nr=nr/2.;
        dndx=0; drefndx=0;
    }
};

class System{// the whole system, which consists of the guest part and the substrate part. The variable mu_linker is the chemical potential of the linker, which is varied to calculate theta and alpha as a function of mu_linker. The variable fcnf is the free energy cost of forming a bond between a ligand and a receptor, which is related to the non-specific interactions between the ligands and the receptors. The variable beta is the inverse temperature. The variable zg is the concentration of the guest in the bulk solution, which is related to the standard state of the system. The variable Feff is the effective free energy of binding, which is calculated from the self-consistent equations. The variable q is related to theta by q=theta/(zg*(1-theta)). The variables dMldx, dMrdx, dQdx are the derivatives of M_l, M_r, Q with respect to mu_linker, which are calculated from the self-consistent equations. The variables drefMldx and drefMrdx are the corresponding variables for the reference system, which are used to calculate the derivative of Feff with respect to mu_linker. The variable dFdx is the derivative of Feff with respect to mu_linker, which is used to find mu_linker for a given target theta. The variable dqdx is the derivative of q with respect to mu_linker, which is used to calculate dlnthetadx. The variable dlnthetadx is the derivative of ln(theta/(zg*(1-theta))) with respect to mu_linker, which is used to calculate alpha. The function initialize_system initializes all the variables in the system. The function solve_system solves the self-consistent equations for bar_nl and bar_nr using Newton's method. The function picard_solve_system solves the self-consistent equations for bar_nl and bar_nr using Picard iteration. The function solve_ref_system solves the self-consistent equations for bar_nl and bar_nr for the reference system using Newton's method. The function picard_solve_ref_system solves the self-consistent equations for bar_nl and bar_nr for the reference system using Picard iteration. The function cal_Feff calculates Feff from bar_nl and bar_nr. The function cal_q calculates q from Feff. The function cal_dFeffdx calculates dFdx from dndx, dMldx, dMrdx, dQdx, d
public:
    Guest guest_part;
    Substrate substrate_plate;
    double mu_linker, fcnf;
    double beta;
    double zg;
    double Feff,q;
    double dMldx,dMrdx,dQdx;
    double drefMldx,drefMrdx;
    double dFdx,dqdx,dlnthetadx,alpha;
    void initialize_system(double given_xi_l, double given_kl, double given_nl,double given_xi_r, double given_kr, double given_nr,double given_mu_linker, double given_fcnf,double given_zg,double given_beta){
        guest_part.initialize_guest(given_xi_l,given_kl,given_nl);
        substrate_plate.initialize_substrate(given_xi_r,given_kr,given_nr);
        mu_linker=given_mu_linker;
        fcnf=given_fcnf;
        zg=given_zg;
        beta=given_beta;
    }
    void solve_system(void);
    void picard_solve_system(void);
    void solve_ref_system(void);
    void picard_solve_ref_system(void);
    void cal_Feff(void);
    void cal_q(void);
    void cal_dFeffdx(void);
    void cal_dqdx(void);
    void cal_dlnthetadx(void);
    void cal_alpha(void);
    void find_mu(double target_theta);
    void picard_find_mu(double target_theta);
    void dichotomy_find_mu_ad(double target_theta, double lower_limit, double higher_limit);
    void dichotomy_find_mu_dp(double target_theta, double lower_limit, double upper_limit);
    double sampling_theta(double mu);
};

void System::solve_system(void){//solving self-consistent equations using Newton's method
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr,der_sumkml,der_sumkmr,der_sumkql,der_sumkqr;
    double block_l,block_r;
    double der_f_l,func_l,der_f_r,func_r;
    while(error>0.000000001){
        block_l=guest_part.bar_nl*guest_part.xi_l;
        block_r=substrate_plate.bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.);
        sumkql=exp(mu_linker-fcnf)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.)*(pow(block_r+1.,substrate_plate.kr)-1.);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.);
        sumkqr=exp(mu_linker-fcnf)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.)*(pow(block_l+1.,guest_part.kl)-1.);
        der_sumkml=sumkml/guest_part.bar_nl+sumkml*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkql=sumkql/guest_part.bar_nl+sumkql*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkmr=sumkmr/substrate_plate.bar_nr+sumkmr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);
        der_sumkqr=sumkqr/substrate_plate.bar_nr+sumkqr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);
        
        func_l=guest_part.bar_nl+sumkml+sumkql-guest_part.nl;
        der_f_l=1.+der_sumkml+der_sumkql;
        func_r=substrate_plate.bar_nr+sumkmr+sumkqr-substrate_plate.nr;
        der_f_r=1.+der_sumkmr+der_sumkqr;
        old_bar_nl=guest_part.bar_nl; old_bar_nr=substrate_plate.bar_nr;
        guest_part.bar_nl=guest_part.bar_nl-0.0001*func_l/der_f_l;
        substrate_plate.bar_nr=substrate_plate.bar_nr-0.0001*func_r/der_f_r;
        
        error=max(fabs(old_bar_nl-guest_part.bar_nl)/0.0001,fabs(old_bar_nr-substrate_plate.bar_nr)/0.0001);
    }
    double dfldnr,dfrdnl;
    dfldnr=exp(-fcnf)*sumkml*substrate_plate.kr*pow(block_r+1,substrate_plate.kr-1)*substrate_plate.xi_r;
    dfrdnl=exp(-fcnf)*sumkmr*guest_part.kl*pow(block_l+1,guest_part.kl-1)*guest_part.xi_l;
    double dfldmu,dfrdmu;
    dfldmu=sumkml+sumkql;
    dfrdmu=sumkmr+sumkqr;
    guest_part.dndx=((-dfldmu)*der_f_r-(-dfrdmu)*dfldnr)/(der_f_l*der_f_r-dfldnr*dfrdnl);
    substrate_plate.dndx=-((-dfldmu)*dfrdnl-(-dfrdmu)*der_f_l)/(der_f_l*der_f_r-dfldnr*dfrdnl);
    dMldx=exp(mu_linker)*(pow(block_l+1,guest_part.kl)-1)+sumkml*guest_part.dndx/guest_part.bar_nl; 
    dMrdx=exp(mu_linker)*(pow(block_r+1,substrate_plate.kr)-1)+sumkmr*substrate_plate.dndx/substrate_plate.bar_nr;
    dQdx=exp(mu_linker-fcnf)*(pow(block_l+1,guest_part.kl)*pow(block_r+1,substrate_plate.kr)-pow(block_l+1,guest_part.kl)-pow(block_r+1,substrate_plate.kr)+1)
                        +sumkql*guest_part.dndx/guest_part.bar_nl+sumkqr*substrate_plate.dndx/substrate_plate.bar_nr;
}

void System::picard_solve_system(void){//solving self-consistent equations using Picard iteration
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr,der_sumkml,der_sumkmr,der_sumkql,der_sumkqr;
    double block_l,block_r;
    double der_f_l,func_l,der_f_r,func_r;
    while(error>0.000000001){
        block_l=guest_part.bar_nl*guest_part.xi_l;
        block_r=substrate_plate.bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.);
        sumkql=exp(mu_linker-fcnf)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.)*(pow(block_r+1.,substrate_plate.kr)-1.);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.);
        sumkqr=exp(mu_linker-fcnf)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.)*(pow(block_l+1.,guest_part.kl)-1.);
        der_sumkml=sumkml/guest_part.bar_nl+sumkml*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkql=sumkql/guest_part.bar_nl+sumkql*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkmr=sumkmr/substrate_plate.bar_nr+sumkmr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);
        der_sumkqr=sumkqr/substrate_plate.bar_nr+sumkqr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);

        der_f_l=1.+der_sumkml+der_sumkql;
        der_f_r=1.+der_sumkmr+der_sumkqr;

        func_l=1+(sumkml+sumkql)/guest_part.bar_nl;
        func_r=1+(sumkmr+sumkqr)/substrate_plate.bar_nr;

        old_bar_nl=guest_part.bar_nl; old_bar_nr=substrate_plate.bar_nr;
        guest_part.bar_nl=guest_part.bar_nl*0.995+0.005*guest_part.nl/func_l;
        substrate_plate.bar_nr=substrate_plate.bar_nr*0.995+0.005*substrate_plate.nr/func_r;
        error=max(fabs(old_bar_nl-guest_part.bar_nl),fabs(old_bar_nr-substrate_plate.bar_nr))/0.005;

    }
    double dfldnr,dfrdnl;
    dfldnr=exp(-fcnf)*sumkml*substrate_plate.kr*pow(block_r+1,substrate_plate.kr-1)*substrate_plate.xi_r;
    dfrdnl=exp(-fcnf)*sumkmr*guest_part.kl*pow(block_l+1,guest_part.kl-1)*guest_part.xi_l;
    double dfldmu,dfrdmu;
    dfldmu=sumkml+sumkql;
    dfrdmu=sumkmr+sumkqr;
    guest_part.dndx=((-dfldmu)*der_f_r-(-dfrdmu)*dfldnr)/(der_f_l*der_f_r-dfldnr*dfrdnl);
    substrate_plate.dndx=-((-dfldmu)*dfrdnl-(-dfrdmu)*der_f_l)/(der_f_l*der_f_r-dfldnr*dfrdnl);
    dMldx=exp(mu_linker)*(pow(block_l+1,guest_part.kl)-1)+sumkml*guest_part.dndx/guest_part.bar_nl; 
    dMrdx=exp(mu_linker)*(pow(block_r+1,substrate_plate.kr)-1)+sumkmr*substrate_plate.dndx/substrate_plate.bar_nr;
    dQdx=exp(mu_linker-fcnf)*(pow(block_l+1,guest_part.kl)*pow(block_r+1,substrate_plate.kr)-pow(block_l+1,guest_part.kl)-pow(block_r+1,substrate_plate.kr)+1)
                        +sumkql*guest_part.dndx/guest_part.bar_nl+sumkqr*substrate_plate.dndx/substrate_plate.bar_nr;
}

void System::solve_ref_system(void){//solving self-consistent equations for the reference system using Newton's method
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,der_sumkml,der_sumkmr;
    double der_f_l,func_l,der_f_r,func_r;
    double block_l,block_r;
    while(error>0.0000000001){

        block_l=guest_part.ref_bar_nl*guest_part.xi_l;
        block_r=substrate_plate.ref_bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1,guest_part.kl-1);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1,substrate_plate.kr-1);
        der_sumkml=sumkml/guest_part.ref_bar_nl+sumkml*(guest_part.kl-1)*guest_part.xi_l/(block_l+1);
        der_sumkmr=sumkmr/substrate_plate.ref_bar_nr+sumkmr*(substrate_plate.kr-1)*substrate_plate.xi_r/(block_r+1);

        func_l=guest_part.ref_bar_nl+sumkml-guest_part.nl;
        der_f_l=1+der_sumkml;
        func_r=substrate_plate.ref_bar_nr+sumkmr-substrate_plate.nr;
        der_f_r=1+der_sumkmr;
        old_bar_nl=guest_part.ref_bar_nl; old_bar_nr=substrate_plate.ref_bar_nr;
        guest_part.ref_bar_nl=guest_part.ref_bar_nl-0.0001*func_l/der_f_l;
        substrate_plate.ref_bar_nr=substrate_plate.ref_bar_nr-0.0001*func_r/der_f_r;
        error=max(fabs(old_bar_nl-guest_part.ref_bar_nl)/0.0001,fabs(old_bar_nr-substrate_plate.ref_bar_nr)/0.0001);

    }
    double dfldmu,dfrdmu;
    dfldmu=sumkml;
    dfrdmu=sumkmr;
    guest_part.drefndx=(-dfldmu)/der_f_l;
    substrate_plate.drefndx=(-dfrdmu)/der_f_r;
    drefMldx=exp(mu_linker)*(pow(block_l+1,guest_part.kl)-1)+sumkml*guest_part.drefndx/guest_part.ref_bar_nl;
    drefMrdx=exp(mu_linker)*(pow(block_r+1,substrate_plate.kr)-1)+sumkmr*substrate_plate.drefndx/substrate_plate.ref_bar_nr;
}

void System::picard_solve_ref_system(void){//solving self-consistent equations for the reference system using Picard iteration
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,der_sumkml,der_sumkmr;
    double der_f_l,func_l,der_f_r,func_r;
    double block_l,block_r;
    while(error>0.0000000001){

        block_l=guest_part.ref_bar_nl*guest_part.xi_l;
        block_r=substrate_plate.ref_bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1,guest_part.kl-1);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1,substrate_plate.kr-1);
        der_sumkml=sumkml/guest_part.ref_bar_nl+sumkml*(guest_part.kl-1)*guest_part.xi_l/(block_l+1);
        der_sumkmr=sumkmr/substrate_plate.ref_bar_nr+sumkmr*(substrate_plate.kr-1)*substrate_plate.xi_r/(block_r+1);

        func_l=1+sumkml/guest_part.ref_bar_nl;
        func_r=1+sumkmr/substrate_plate.ref_bar_nr;
        der_f_l=1+der_sumkml;
        der_f_r=1+der_sumkmr;
        old_bar_nl=guest_part.ref_bar_nl; old_bar_nr=substrate_plate.ref_bar_nr;
        guest_part.ref_bar_nl=guest_part.ref_bar_nl*0.995+0.005*guest_part.nl/func_l;
        substrate_plate.ref_bar_nr=substrate_plate.ref_bar_nr*0.995+0.005*substrate_plate.nr/func_r;
        error=max(fabs(old_bar_nl-guest_part.ref_bar_nl),fabs(old_bar_nr-substrate_plate.ref_bar_nr))/0.005;
    }
    double dfldmu,dfrdmu;
    dfldmu=sumkml;
    dfrdmu=sumkmr;
    guest_part.drefndx=(-dfldmu)/der_f_l;
    substrate_plate.drefndx=(-dfrdmu)/der_f_r;
    drefMldx=exp(mu_linker)*(pow(block_l+1,guest_part.kl)-1)+sumkml*guest_part.drefndx/guest_part.ref_bar_nl;
    drefMrdx=exp(mu_linker)*(pow(block_r+1,substrate_plate.kr)-1)+sumkmr*substrate_plate.drefndx/substrate_plate.ref_bar_nr;
}

void System::cal_Feff(void){//calculating the effective free energy of binding from bar_nl and bar_nr
    double summl,summr;
    double sumb;
    double block_l,block_r;
    block_l=pow(guest_part.bar_nl*guest_part.xi_l+1,guest_part.kl);
    block_r=pow(substrate_plate.bar_nr*substrate_plate.xi_r+1,substrate_plate.kr);
    summl=exp(mu_linker)*(block_l-1);
    summr=exp(mu_linker)*(block_r-1);
    sumb=exp(mu_linker-fcnf)*(block_l*block_r-block_l-block_r+1);
    Feff= guest_part.nl*log(guest_part.bar_nl)-guest_part.bar_nl-summl
                +substrate_plate.nr*log(substrate_plate.bar_nr)-substrate_plate.bar_nr-summr
                -sumb;
    //set a reference
    double ref_summl,ref_summr;
    double ref_block_l,ref_block_r;
    ref_block_l=pow(guest_part.ref_bar_nl*guest_part.xi_l+1,guest_part.kl);
    ref_block_r=pow(substrate_plate.ref_bar_nr*substrate_plate.xi_r+1,substrate_plate.kr);
    ref_summl=exp(mu_linker)*(ref_block_l-1);
    ref_summr=exp(mu_linker)*(ref_block_r-1);

    Feff=Feff-guest_part.nl*log(guest_part.ref_bar_nl)+guest_part.ref_bar_nl+ref_summl
                -substrate_plate.nr*log(substrate_plate.ref_bar_nr)+substrate_plate.ref_bar_nr+ref_summr;

}

void System::cal_q(void){
    q=exp(-beta*Feff)-1.;
}

void System::cal_dFeffdx(void){//calculating the derivative of Feff with respect to mu_linker, which is used to find alpha
    dFdx=(guest_part.nl/guest_part.bar_nl-1)*guest_part.dndx-dMldx
            +(substrate_plate.nr/substrate_plate.bar_nr-1)*substrate_plate.dndx-dMrdx-dQdx
            -(guest_part.nl/guest_part.ref_bar_nl-1)*guest_part.drefndx+drefMldx
            -(substrate_plate.nr/substrate_plate.ref_bar_nr-1)*substrate_plate.drefndx+drefMrdx;
}

void System::cal_dqdx(void){
    dqdx=-exp(-Feff)*dFdx;
}

void System::cal_dlnthetadx(void){
    dlnthetadx=(1./q-zg/(1.+zg*q))*dqdx;
}

void System::cal_alpha(void){
    alpha=dlnthetadx;
}

void System::find_mu(double target_theta){
    double target_q=target_theta/(zg*(1-target_theta));
    double target_Feff=-log(1+target_q);
    double deltaF,step;
    double old_mu=mu_linker;
    double error=1;
    while(error>0.000001){
        old_mu=mu_linker;
        solve_system();
        solve_ref_system();
        cal_Feff();
        cal_dFeffdx();
        deltaF=Feff-target_Feff;
        step=deltaF/dFdx;
        // step=min(step,5.);
        // step=max(step,-5.);
        mu_linker=mu_linker-0.001*step;
        error=fabs(mu_linker-old_mu)/0.001;
        cout<<mu_linker<<endl;
    }
}

void System::picard_find_mu(double target_theta){//
    double target_q=target_theta/(zg*(1-target_theta));
    double target_Feff=-log(1+target_q);
    double deltaF,step;
    double old_mu=mu_linker;
    double error=1;
    while(error>0.000001){
        old_mu=mu_linker;
        picard_solve_system();
        picard_solve_ref_system();
        cal_Feff();
        cal_dFeffdx();
        deltaF=Feff-target_Feff;
        step=deltaF/dFdx;
        // step=min(step,5.);
        // step=max(step,-5.);
        mu_linker=mu_linker-0.001*step;
        error=fabs(mu_linker-old_mu)/0.001;
        cout<<mu_linker<<endl;
    }
}

void System::dichotomy_find_mu_ad(double target_theta, double lower_limit, double upper_limit){//find the adsorption threshold
    double middle=(lower_limit+upper_limit)/2;
    double delta=upper_limit-lower_limit;
    double theta;
    while (delta>0.000000000000001){
        middle=(lower_limit+upper_limit)/2.;
        theta=sampling_theta(middle);
        if(theta>target_theta){
            upper_limit=middle;
        }else{
            lower_limit=middle;
        }
        delta=fabs(upper_limit-lower_limit);
    }
    mu_linker=middle;
}

void System::dichotomy_find_mu_dp(double target_theta, double lower_limit, double upper_limit){
    double middle=(lower_limit+upper_limit)/2;
    double delta=upper_limit-lower_limit;
    double theta;
    while (delta>0.000000000000001){
        middle=(lower_limit+upper_limit)/2.;
        theta=sampling_theta(middle);
        if(theta<target_theta){
            upper_limit=middle;
        }else{
            lower_limit=middle;
        }
        delta=fabs(upper_limit-lower_limit);
    }
    mu_linker=middle;
}

double System::sampling_theta(double mu){
    vector<double> thetas={0};
    double average_nr=substrate_plate.nr;
    mu_linker=mu;
    for(int j=1;j<55;j++){
        substrate_plate.nr=j;
        picard_solve_system();
        picard_solve_ref_system();
        cal_Feff();
        cal_q();
        thetas.push_back(q*zg/(1.+q*zg));
        // cout<<"nr: "<< j<<" "<<"corresponding theta"<<thetas[j]<<endl;
        cout<<"nr: "<<j<< ", barnl: "<< guest_part.bar_nl<<", barnr: "<<substrate_plate.bar_nr<<endl;
        cout<<"slope: "<<log(guest_part.bar_nl+1)+log(substrate_plate.bar_nr+1)<<endl;
    }
    double result=0;
    for(int i=1;i<55;i++){
        result+=thetas[i]*pow(average_nr,i)*exp(-average_nr)/tgamma(i+1);
    }
    substrate_plate.nr=average_nr;
    return result;
}