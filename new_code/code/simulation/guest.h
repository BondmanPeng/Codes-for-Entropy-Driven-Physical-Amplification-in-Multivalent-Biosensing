#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>

using namespace std;

mt19937 randNum(time(nullptr));

const double guest_size=1;

double randCV(void){
    return randNum()*(1.0/4294967296.0);
}

class Guest{
public:
    double r_eff=0.5;
    double r=0.5;
    double sigma=2*r;
    double boxs[2],backup_boxs[2];
    double nl,xi_l;
    double nr,backup_nr;
    double bar_nl,ref_bar_nl;
    double backup_bar_nl;

    double bar_nr,ref_bar_nr;
    double backup_bar_nr,backup_ref_bar_nr;

    double Feff,backup_Feff;
    double q, backup_q;
    // double dnl,drefnl;
    // double backup_dnl;
    // double dnr,drefnr;
    // double backup_dnr,backup_drefnr;

    Guest *list_pre_pointer, *list_next_pointer;
    Guest *cell_pre_pointer, *cell_next_pointer;
    void initialize_guest(double given_xi_l,double given_nl);
    void random_move(double& ds);
    void backup(void);
    void reject(void);
};

void Guest::initialize_guest(double given_xi_l,double given_nl){
        xi_l=given_xi_l; nl=given_nl; bar_nl=nl/2.; ref_bar_nl=nl/2.;
        nr=1; bar_nr=1; ref_bar_nr=1;
        Feff=0;q=0;
        // dnl=0; drefnl=0; 
        // dnr=0;drefnr=0;
        boxs[0]=randCV();
        boxs[1]=randCV();
        list_pre_pointer=nullptr; list_next_pointer=nullptr;
        cell_pre_pointer=nullptr; cell_next_pointer=nullptr;
}

void Guest::random_move(double& ds){
    for(int i=0; i<2;i++){
        double abs_pos=(boxs[i]+2*ds*(randCV()-0.5));
        boxs[i]=abs_pos-floor(abs_pos);
    }
}

void Guest::backup(void){
    backup_nr=nr;
    backup_bar_nl=bar_nl;
    // backup_dnl=dnl;
    backup_bar_nr=bar_nr;
    backup_ref_bar_nr=ref_bar_nr;
    backup_Feff=Feff;
    backup_q=q;
    // backup_dnr=dnr;
    // backup_drefnr=drefnr;

    backup_boxs[0]=boxs[0];
    backup_boxs[1]=boxs[1];
}

void Guest::reject(void){
    nr=backup_nr;
    bar_nl=backup_bar_nl;
    // dnl=backup_dnl;
    bar_nr=backup_bar_nr;
    ref_bar_nr=backup_ref_bar_nr;
    Feff=backup_Feff;
    q=backup_q;
    // dnr=backup_dnr;
    // drefnr=backup_drefnr;
    boxs[0]=backup_boxs[0];
    boxs[1]=backup_boxs[1];
}