#include "guest.h"

using namespace std;

double generate_a_Poisson_arrival_interval(double lambda){
    double U=randCV();
    return -(1./lambda)*log(U);
}

double generate_Poisson_arrival_counts(double lambda,double area){
    double summation=0;
    int counts=0;
    while(summation>=-lambda*area){
        counts++;
        summation+=log(randCV());
    }
    counts-=1;
    return counts;
}

class Receptor{
public:
    int index;
    double xi_r;
    double position[2];
    void initialize_a_receptor(double label_number,double binding_strength);
};

void Receptor::initialize_a_receptor(double label_number, double binding_strength){
    index=label_number;
    xi_r=binding_strength;
    position[0]=randCV();
    position[1]=randCV();
}

class Rec_cells{
public:
    vector<int> inside_receptors={};
    void take_in_one_receptor(int i);
};

void Rec_cells::take_in_one_receptor(int i){
    inside_receptors.push_back(i);
}

class Substrate{
public:
    double Nr;
    double rho_r;
    double area;
    double L;
    double xi_r;
    int cell_limit;
    double cell_length;
    vector<Receptor> receptors;
    vector<vector<Rec_cells>> receptor_cells;
    void initialize_receptors(double sigma, double given_rho_r, double given_area,double binding_strength);
};

void Substrate::initialize_receptors(double sigma, double given_rho_r, double given_area,double binding_strength){
    rho_r=given_rho_r;
    area=given_area;
    xi_r=binding_strength;
    L=sqrt(area);
    Nr=generate_Poisson_arrival_counts(rho_r,area);
    // Nr=round(12.73*area);
    cell_limit=floor(L/sigma);
    cell_length= 1./((double) cell_limit);
    vector<Receptor> initialized_receptors(Nr);
    vector<vector<Rec_cells>> initialized_receptor_cell(cell_limit,vector<Rec_cells>(cell_limit));
    receptors=initialized_receptors;
    receptor_cells=initialized_receptor_cell;
    for(int i=0;i<Nr;i++){
        receptors[i].initialize_a_receptor(i,binding_strength);
        int index_x,index_y;
        index_x=floor(receptors[i].position[0]/cell_length);
        index_y=floor(receptors[i].position[1]/cell_length);
        receptor_cells[index_x][index_y].take_in_one_receptor(i);
    }
}
