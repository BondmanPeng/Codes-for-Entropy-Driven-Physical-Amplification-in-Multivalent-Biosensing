#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

class Guest{
public:
    double xi_l,kl,nl,bar_nl,ref_bar_nl,xi_lm;
    double dndx,drefndx;
    void initialize_guest(double given_xi_l,double given_xi_lm,double given_kl,double given_nl){
        xi_l=given_xi_l; xi_lm=given_xi_lm;
        kl=given_kl; nl=given_nl; bar_nl=nl/2.; ref_bar_nl=nl/2.;
        dndx=0; drefndx=0;
    }
};

class Substrate{
public:
    double xi_r,kr, nr,bar_nr, ref_bar_nr,xi_rm;
    double dndx,drefndx;
    double nrp;
    void initialize_substrate(double given_xi_r, double given_xi_rm,double given_kr, double given_nr){
        xi_r=given_xi_r; xi_rm=given_xi_rm; 
        kr=given_kr; nr=given_nr; bar_nr=nr/2.; ref_bar_nr=nr/2.;
        dndx=0; drefndx=0;
    }
};

class System{
public:
    Guest guest_part;
    Substrate substrate_plate;
    double mu_linker, fcnf;
    double mu_m;
    double beta;
    double zg;
    double Feff,q;
    double dMldx,dMrdx,dQdx,dEaldx,dEardx;
    double drefMldx,drefMrdx,drefEaldx,drefEardx;
    double dFdx,dqdx,dlnthetadx,alpha;
    void initialize_system(double given_xi_l, double given_xi_lm,double given_kl, double given_nl,double given_xi_r, double given_xi_rm,double given_kr, double given_nr,double given_mu_linker,double given_mu_m, double given_fcnf,double given_zg,double given_beta){
        guest_part.initialize_guest(given_xi_l,given_xi_lm,given_kl,given_nl);
        substrate_plate.initialize_substrate(given_xi_r,given_xi_rm,given_kr,given_nr);
        mu_linker=given_mu_linker;
        mu_m=given_mu_m;
        fcnf=given_fcnf;
        zg=given_zg;
        beta=given_beta;
    }
    void picard_solve_system(void);
    void picard_solve_ref_system(void);
    void cal_Feff(void);
    void cal_q(void);
    void cal_dFeffdx(void);
    void cal_dqdx(void);
    void cal_dlnthetadx(void);
    void cal_alpha(void);

};

void System::picard_solve_system(void){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr, Eal, Ear;
    double der_sumkml,der_sumkmr,der_sumkql,der_sumkqr;
    double block_l,block_r;
    double der_f_l,func_l,der_f_r,func_r;
    double der_Eal, der_Ear;
    while(error>0.000000001){
        block_l=guest_part.bar_nl*guest_part.xi_l;
        block_r=substrate_plate.bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.);
        sumkql=exp(mu_linker-fcnf)*guest_part.kl*block_l*pow(block_l+1.,guest_part.kl-1.)*(pow(block_r+1.,substrate_plate.kr)-1.);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.);
        sumkqr=exp(mu_linker-fcnf)*substrate_plate.kr*block_r*pow(block_r+1.,substrate_plate.kr-1.)*(pow(block_l+1.,guest_part.kl)-1.);
        Eal=guest_part.bar_nl*exp(mu_m)*guest_part.xi_lm;
        Ear=substrate_plate.bar_nr*exp(mu_m)*substrate_plate.xi_rm;
        // Eb=guest_part.bar_nl*substrate_plate.bar_nr*exp(mu_m-fcnf)*guest_part.xi_lm*substrate_plate.xi_rm;
        der_sumkml=sumkml/guest_part.bar_nl+sumkml*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkql=sumkql/guest_part.bar_nl+sumkql*(guest_part.kl-1.)*guest_part.xi_l/(block_l+1.);
        der_sumkmr=sumkmr/substrate_plate.bar_nr+sumkmr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);
        der_sumkqr=sumkqr/substrate_plate.bar_nr+sumkqr*(substrate_plate.kr-1.)*substrate_plate.xi_r/(block_r+1.);
        der_Eal=exp(mu_m)*guest_part.xi_lm; 
        der_Ear=exp(mu_m)* substrate_plate.xi_rm;
        // der_Ebl=exp(mu_linker-fcnf)*substrate_plate.bar_nr*guest_part.xi_lm*substrate_plate.xi_rm;
        // der_Ebr=exp(mu_linker-fcnf)*guest_part.bar_nl*guest_part.xi_lm*substrate_plate.xi_rm;

        der_f_l=1.+der_sumkml+der_sumkql+der_Eal;
        der_f_r=1.+der_sumkmr+der_sumkqr+der_Ear;

        func_l=1+(sumkml+sumkql+Eal)/guest_part.bar_nl;
        func_r=1+(sumkmr+sumkqr+Ear)/substrate_plate.bar_nr;

        old_bar_nl=guest_part.bar_nl; old_bar_nr=substrate_plate.bar_nr;
        guest_part.bar_nl=guest_part.bar_nl*0.99+0.01*guest_part.nl/func_l;
        substrate_plate.bar_nr=substrate_plate.bar_nr*0.99+0.01*substrate_plate.nr/func_r;
        error=max(fabs(old_bar_nl-guest_part.bar_nl),fabs(old_bar_nr-substrate_plate.bar_nr))/0.01;

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
    dEaldx=Eal*guest_part.dndx/guest_part.bar_nl;
    dEardx=Ear*substrate_plate.dndx/substrate_plate.bar_nr;
    // dEbdx=Eb*(guest_part.dndx/guest_part.bar_nl+substrate_plate.dndx/substrate_plate.bar_nr);
}

void System::picard_solve_ref_system(void){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,Eal,Ear;
    double der_sumkml,der_sumkmr;
    double der_Eal, der_Ear;
    double der_f_l,func_l;
    double der_f_r,func_r;
    double block_l,block_r;
    while(error>0.000000001){

        block_l=guest_part.ref_bar_nl*guest_part.xi_l;
        block_r=substrate_plate.ref_bar_nr*substrate_plate.xi_r;
        sumkml=exp(mu_linker)*guest_part.kl*block_l*pow(block_l+1,guest_part.kl-1);
        sumkmr=exp(mu_linker)*substrate_plate.kr*block_r*pow(block_r+1,substrate_plate.kr-1);
        Eal=guest_part.ref_bar_nl*exp(mu_m)*guest_part.xi_lm;
        Ear=substrate_plate.ref_bar_nr*exp(mu_m)*substrate_plate.xi_rm;
        der_sumkml=sumkml/guest_part.ref_bar_nl+sumkml*(guest_part.kl-1)*guest_part.xi_l/(block_l+1);
        der_sumkmr=sumkmr/substrate_plate.ref_bar_nr+sumkmr*(substrate_plate.kr-1)*substrate_plate.xi_r/(block_r+1);
        der_Eal=exp(mu_m)*guest_part.xi_lm; 
        der_Ear=exp(mu_m)*substrate_plate.xi_rm;

        func_l=1+(sumkml+Eal)/guest_part.ref_bar_nl;
        func_r=1+(sumkmr+Ear)/substrate_plate.ref_bar_nr;
        der_f_l=1+der_sumkml+der_Eal;
        der_f_r=1+der_sumkmr+der_Ear;
        old_bar_nl=guest_part.ref_bar_nl; old_bar_nr=substrate_plate.ref_bar_nr;
        guest_part.ref_bar_nl=guest_part.ref_bar_nl*0.99+0.01*guest_part.nl/func_l;
        substrate_plate.ref_bar_nr=substrate_plate.ref_bar_nr*0.99+0.01*substrate_plate.nr/func_r;
        error=max(fabs(old_bar_nl-guest_part.ref_bar_nl),fabs(old_bar_nr-substrate_plate.ref_bar_nr))/0.01;
    }
    double dfldmu,dfrdmu;
    dfldmu=sumkml;
    dfrdmu=sumkmr;
    guest_part.drefndx=(-dfldmu)/der_f_l;
    substrate_plate.drefndx=(-dfrdmu)/der_f_r;
    drefMldx=exp(mu_linker)*(pow(block_l+1,guest_part.kl)-1)+sumkml*guest_part.drefndx/guest_part.ref_bar_nl;
    drefMrdx=exp(mu_linker)*(pow(block_r+1,substrate_plate.kr)-1)+sumkmr*substrate_plate.drefndx/substrate_plate.ref_bar_nr;
    drefEaldx=Eal*guest_part.drefndx/guest_part.ref_bar_nl;
    drefEardx=Ear*substrate_plate.drefndx/substrate_plate.ref_bar_nr;

}

void System::cal_Feff(void){
    double summl,summr;
    double sumb;
    double block_l,block_r;
    double Eal, Ear;
    block_l=pow(guest_part.bar_nl*guest_part.xi_l+1,guest_part.kl);
    block_r=pow(substrate_plate.bar_nr*substrate_plate.xi_r+1,substrate_plate.kr);
    summl=exp(mu_linker)*(block_l-1);
    summr=exp(mu_linker)*(block_r-1);
    sumb=exp(mu_linker-fcnf)*(block_l*block_r-block_l-block_r+1);
    Eal=guest_part.bar_nl*exp(mu_m)*guest_part.xi_lm;
    Ear=substrate_plate.bar_nr*exp(mu_m)*substrate_plate.xi_rm;
    
    Feff= guest_part.nl*log(guest_part.bar_nl)-guest_part.bar_nl-summl
                +substrate_plate.nr*log(substrate_plate.bar_nr)-substrate_plate.bar_nr-summr
                -sumb-Eal-Ear;
    //set a reference
    double ref_summl,ref_summr;
    double ref_block_l,ref_block_r;
    double ref_Eal,ref_Ear;
    ref_block_l=pow(guest_part.ref_bar_nl*guest_part.xi_l+1,guest_part.kl);
    ref_block_r=pow(substrate_plate.ref_bar_nr*substrate_plate.xi_r+1,substrate_plate.kr);
    ref_summl=exp(mu_linker)*(ref_block_l-1);
    ref_summr=exp(mu_linker)*(ref_block_r-1);
    ref_Eal=guest_part.ref_bar_nl*exp(mu_m)*guest_part.xi_lm;
    ref_Ear=substrate_plate.ref_bar_nr*exp(mu_m)*substrate_plate.xi_rm;

    Feff=Feff-guest_part.nl*log(guest_part.ref_bar_nl)+guest_part.ref_bar_nl+ref_summl
                -substrate_plate.nr*log(substrate_plate.ref_bar_nr)+substrate_plate.ref_bar_nr+ref_summr+ref_Eal+ref_Ear;
    if(mu_linker-floor(mu_linker)<0.01){
        cout<<"linker mu: "<<mu_linker<<endl;
        cout<<"nr: "<<substrate_plate.nr<<" nl: "<<guest_part.nl<<endl;
        cout<<"Feff: "<<Feff<<endl;
    }
    
}

void System::cal_q(void){
    q=exp(-beta*Feff)-1.;
}

void System::cal_dFeffdx(void){
    dFdx=(guest_part.nl/guest_part.bar_nl-1)*guest_part.dndx-dMldx-dEaldx
            +(substrate_plate.nr/substrate_plate.bar_nr-1)*substrate_plate.dndx-dMrdx-dEardx-dQdx
            -(guest_part.nl/guest_part.ref_bar_nl-1)*guest_part.drefndx+drefMldx+drefEaldx
            -(substrate_plate.nr/substrate_plate.ref_bar_nr-1)*substrate_plate.drefndx+drefMrdx+drefEardx;
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