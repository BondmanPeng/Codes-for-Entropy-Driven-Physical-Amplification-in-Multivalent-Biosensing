#include "substrate.h"
#include "array_computation.h"

using namespace std;

const int max_Nr=35;

Guest* randomLink(Guest* headPointer,int N);

Guest* randomLink(Guest* headPointer,int N){
    int randomNumber;
    Guest* ptcPosition=headPointer;
    randomNumber=randNum()%N;
    for(int i =0; i<randomNumber;i++){
        ptcPosition=ptcPosition->list_next_pointer;
    }
    return ptcPosition;
}

class Cell_List {//cell list for storing particles in the same cell, which is used for checking the overlap between particles and for counting the number of neighbors of a particle. 
public:
    Guest *headPtr; // head pointer that points to the first particle in the cell
    int index; // size of the cell

    void insert(Guest *pNewPtr) { // a function to insert a particle pointed by pNewPtr in a chosen cell
    // cannot insert a particle in a corresponding cell automatically, but need to pick a cell first
    // the inserted particle is pointed by the head pointer
        if (headPtr != nullptr) {
            pNewPtr->cell_next_pointer = headPtr;
            headPtr->cell_pre_pointer=pNewPtr;
        } else {
            pNewPtr->cell_next_pointer = nullptr;
        }
        pNewPtr->cell_pre_pointer = nullptr;
        headPtr = pNewPtr;

    }

    void remove(Guest *pOldPtr) {// a function to remove a particle from a cell
    // need to pick a cell first to make it work
    // to cut the linked list and to link the two neighbor particle around the removed particle
        if (pOldPtr != headPtr){
           (pOldPtr->cell_pre_pointer)->cell_next_pointer = pOldPtr->cell_next_pointer;
        }else{
            headPtr = pOldPtr->cell_next_pointer;
        }

        if (pOldPtr->cell_next_pointer != nullptr) {
            (pOldPtr->cell_next_pointer)->cell_pre_pointer = pOldPtr->cell_pre_pointer;
        }
        pOldPtr->cell_next_pointer = nullptr;
        pOldPtr->cell_pre_pointer = nullptr;
    }
};

class Particle_List {
public:
    Guest *headPtr; // head pointer that points to the first particle in the cell
    int index; // size of the cell

    void insert(Guest *pNewPtr) { // a function to insert a particle pointed by pNewPtr in a chosen cell
    // cannot insert a particle in a corresponding cell automatically, but need to pick a cell first
    // the inserted particle is pointed by the head pointer
        if (headPtr != nullptr) {
            pNewPtr->list_next_pointer = headPtr;
            headPtr->list_pre_pointer=pNewPtr;
        } else {
            pNewPtr->list_next_pointer = nullptr;
        }
        pNewPtr->list_pre_pointer = nullptr;
        headPtr = pNewPtr;

    }

    void remove(Guest *pOldPtr) {// a function to remove a particle from a cell
    // need to pick a cell first to make it work
    // to cut the linked list and to link the two neighbor particle around the removed particle
        if (pOldPtr != headPtr){
           (pOldPtr->list_pre_pointer)->list_next_pointer = pOldPtr->list_next_pointer;
            
        }else{
            headPtr = pOldPtr->list_next_pointer;
        }

        if (pOldPtr->list_next_pointer != nullptr) {
            (pOldPtr->list_next_pointer)->list_pre_pointer = pOldPtr->list_pre_pointer;
        }
        pOldPtr->list_next_pointer = nullptr;
        pOldPtr->list_pre_pointer = nullptr;
    }
};

class System{// the system class that contains all the information of the system, including the parameters of the system, the particles in the system, and the functions to initialize the system, to check for the overlap between particles, to count the number of neighbors of a particle, and to solve the self-consistent equations for calculating Feff and theta.
public:    
    double Nmax,theta,alpha;
    int Ng;
    double mu_g,rho_r,mu_linker,nl;
    double kr,kl;
    double xi_r,xi_l;
    // double zg;
    double fcnf;

    double Nsteps=0,accSteps=0;

    double ds=0.01;
    
    double sigma=1.;
    double L;
    int cell_limit;
    double cell_length;
    
    double interactions[max_Nr];
    Particle_List guest_particles;
    vector<vector<Cell_List>> cells;

    Substrate rec_substrate;
    void initialize_system(double given_mu_g, double given_rho_r, double given_mu_linker, double given_nl, double given_kr,double given_kl,double given_xi_r, double given_xi_l, double given_fcnf,double given_area);
    void initialize_interactions(void);
    void check_for_Feff(Guest* the_particle);
    int check_overlap(Guest* the_particle);
    void count_nr(Guest* the_particle);
    void solve_system(Guest* the_particle);
    void solve_system_restore(double& bar_nl,double& bar_nr,double nr);
    void picard_solve_system(Guest* the_particle);
    void picard_solve_system_restore(double& bar_nl,double& bar_nr,double nr);
    void solve_ref_system(Guest* the_particle);
    void solve_ref_system_restore(double& ref_bar_nl,double& ref_bar_nr,double nr);
    void picard_solve_ref_system(Guest* the_particle);
    void picard_solve_ref_system_restore(double& ref_bar_nl,double& ref_bar_nr,double nr);
    void cal_Feff(Guest* the_particle);
    void update_q(Guest* the_particle);
    void extract_particle(Guest* the_particle);
    void put_particle(Guest* the_particle);
    void trans_move(void);
    void add_particles(void);
    void remove_particles(void);
    void MCmove(void);
    void adapt_step(void);
    // sampling
    void cal_alpha(void);

};

void System::initialize_system(double given_mu_g, double given_rho_r, double given_mu_linker, double given_nl, double given_kr,double given_kl,double given_xi_r, double given_xi_l, double given_fcnf,double given_area){
    alpha=0;
    L=sqrt(given_area);
    cell_limit=floor(L/sigma);
    cell_length= 1./((double) cell_limit);
    // initialize guest part: no guests inside
    guest_particles.headPtr=nullptr;
    vector<vector<Cell_List>> initialized_cells(cell_limit,vector<Cell_List>(cell_limit));
    cells=initialized_cells;
    for(int i =0;i<cell_limit;i++){// initialize the cell matrix
        for(int j=0;j<cell_limit;j++){
            cells[i][j].headPtr=nullptr;
        }
    }
    // initialize substrate part
    rec_substrate.initialize_receptors(sigma,given_rho_r,given_area,given_xi_r);
    // initialize interaction part
    Nmax=given_area/(M_PI*pow(sigma/2.,2));
    Ng=0; theta=0; 
    mu_g=given_mu_g; rho_r=given_rho_r; mu_linker=given_mu_linker; nl=given_nl;
    kr=given_kr; kl=given_kl; xi_r=given_xi_r; xi_l=given_xi_l;
    fcnf=given_fcnf;

}

void System::initialize_interactions(void){
    interactions[0]=0;
    for(int i=1;i<max_Nr;++i){
        double Feff=0;
        double nr=i;
        double bar_nl,bar_nr,ref_bar_nl,ref_bar_nr;
        bar_nl=nl/2.; ref_bar_nl=nl/2.;
        bar_nr=nr/2.; ref_bar_nr=nr/2.;
        picard_solve_system_restore(bar_nl,bar_nr,nr);
        picard_solve_ref_system_restore(ref_bar_nl,ref_bar_nr,nr);
        double summl,summr;
        double sumb;
        double block_l,block_r;
        block_l=pow(bar_nl*xi_l+1,kl);
        block_r=pow(bar_nr*xi_r+1,kr);
        summl=exp(mu_linker)*(block_l-1);
        summr=exp(mu_linker)*(block_r-1);
        sumb=exp(mu_linker-fcnf)*(block_l*block_r-block_l-block_r+1);
        Feff= nl*log(bar_nl)-bar_nl-summl
                    +nr*log(bar_nr)-bar_nr-summr
                    -sumb;
        //set a reference
        double ref_summl,ref_summr;
        double ref_block_l,ref_block_r;
        ref_block_l=pow(ref_bar_nl*xi_l+1,kl);
        ref_block_r=pow(ref_bar_nr*xi_r+1,kr);
        ref_summl=exp(mu_linker)*(ref_block_l-1);
        ref_summr=exp(mu_linker)*(ref_block_r-1);

        Feff=Feff-nl*log(ref_bar_nl)+ref_bar_nl+ref_summl
                    -nr*log(ref_bar_nr)+ref_bar_nr+ref_summr;
        interactions[i]=Feff;
        // cout<<Feff<<endl;
    }
}

void System::check_for_Feff(Guest* the_particle){
    if(the_particle->nr<max_Nr){
        int n_rec=the_particle->nr;
        the_particle->Feff=interactions[n_rec];
    }else{
        picard_solve_system(the_particle);
        picard_solve_ref_system(the_particle);
        cal_Feff(the_particle);
        cout<<"exception"<<endl;
        cout<<"nr="<<the_particle->nr<<endl;
    }
}

int System::check_overlap(Guest* the_particle){
    int neighborX,neighborY;
    int cell_x,cell_y;
    double distance;
    double dr[2],iter_box_adj[2];
    int periodic_adj[2];
    cell_x=floor(the_particle->boxs[0]/cell_length);
    cell_y=floor(the_particle->boxs[1]/cell_length);
    for(Guest* iter= the_particle->cell_next_pointer; iter !=nullptr; iter=iter->cell_next_pointer){
        Substract(dr,iter->boxs,the_particle->boxs);
        distance=Norm(dr)*L;
        if(distance<the_particle->sigma){
            return 1;
        }
    }
    for(Guest* iter= the_particle->cell_pre_pointer; iter !=nullptr; iter=iter->cell_pre_pointer){
        Substract(dr,iter->boxs,the_particle->boxs);
        distance=Norm(dr)*L;
        if(distance<the_particle->sigma){
            return 1;
        }
    }
    for(int q=-1; q<2;q++)// check its neighboring cells
    for(int n=-1; n<2;n++){
        if(q!=0 || n!=0){
            
            neighborX=cell_x+q;neighborY=cell_y+n;
            periodic_adj[0]=0;periodic_adj[1]=0;

            if(neighborX<0) periodic_adj[0]=-1;
            else if(neighborX>=rec_substrate.cell_limit) periodic_adj[0]=1;

            if(neighborY<0) periodic_adj[1]=-1;
            else if(neighborY>=rec_substrate.cell_limit) periodic_adj[1]=1;

            neighborX-=rec_substrate.cell_limit*periodic_adj[0]; neighborY-=rec_substrate.cell_limit*periodic_adj[1];
            
            Guest* ngbPtr= cells[neighborX][neighborY].headPtr;
            for(Guest* iter=ngbPtr;iter!=nullptr;iter=iter->cell_next_pointer){
                periodic_boundary(iter_box_adj,iter->boxs,periodic_adj);
                Substract(dr,iter_box_adj,the_particle->boxs);
                distance=Norm(dr)*L;
                if(distance<the_particle->sigma){
                    return 1;
                }
            }  
        }
    }
    return 0;
}

void System::count_nr(Guest* the_particle){
    int count=0;
    int neighborX,neighborY;
    int cell_x,cell_y;
    double distance;
    double dr[2],iter_box_adj[2];
    int periodic_adj[2];
    cell_x=floor(the_particle->boxs[0]/cell_length);
    cell_y=floor(the_particle->boxs[1]/cell_length);
    for(int q=-1; q<2;q++)// check its neighboring cells
    for(int n=-1; n<2;n++){
        neighborX=cell_x+q;neighborY=cell_y+n;
        periodic_adj[0]=0;periodic_adj[1]=0;

        if(neighborX<0) periodic_adj[0]=-1;
        else if(neighborX>=cell_limit) periodic_adj[0]=1;

        if(neighborY<0) periodic_adj[1]=-1;
        else if(neighborY>=cell_limit) periodic_adj[1]=1;

        neighborX-=cell_limit*periodic_adj[0]; neighborY-=rec_substrate.cell_limit*periodic_adj[1];

        int length=rec_substrate.receptor_cells[neighborX][neighborY].inside_receptors.size();
        for(int i=0;i<length;i++){
            int receptor_index=rec_substrate.receptor_cells[neighborX][neighborY].inside_receptors[i];
            periodic_boundary(iter_box_adj,rec_substrate.receptors[receptor_index].position,periodic_adj);
            Substract(dr,iter_box_adj,the_particle->boxs);
            distance=Norm(dr)*L;
            if(distance<the_particle->r_eff){
                count++;
            }
        }
    }
    the_particle->nr=count;
    the_particle->bar_nr=count/2.;
    the_particle->ref_bar_nr=count/2.;
}

void System::solve_system(Guest* the_particle){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr,der_sumkml,der_sumkmr,der_sumkql,der_sumkqr;
    double block_l,block_r;
    double der_f_l,func_l,der_f_r,func_r;
    // double grad_r,grad_l;
    while(error>0.00000001){
        block_l=the_particle->bar_nl*xi_l;
        block_r=the_particle->bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*block_l*pow(block_l+1.,kl-1.);
        sumkql=exp(mu_linker-fcnf)*kl*block_l*pow(block_l+1.,kl-1.)*(pow(block_r+1.,kr)-1.);
        sumkmr=exp(mu_linker)*kr*block_r*pow(block_r+1.,kr-1.);
        sumkqr=exp(mu_linker-fcnf)*kr*block_r*pow(block_r+1.,kr-1.)*(pow(block_l+1.,kl)-1.);
        der_sumkml=sumkml/the_particle->bar_nl+sumkml*(kl-1.)*the_particle->xi_l/(block_l+1.);
        der_sumkql=sumkql/the_particle->bar_nl+sumkql*(kl-1.)*the_particle->xi_l/(block_l+1.);
        der_sumkmr=sumkmr/the_particle->bar_nr+sumkmr*(kr-1.)*rec_substrate.xi_r/(block_r+1.);
        der_sumkqr=sumkqr/the_particle->bar_nr+sumkqr*(kr-1.)*rec_substrate.xi_r/(block_r+1.);
        
        func_l=the_particle->bar_nl+sumkml+sumkql-the_particle->nl;
        der_f_l=1.+der_sumkml+der_sumkql;
        func_r=the_particle->bar_nr+sumkmr+sumkqr-the_particle->nr;
        der_f_r=1.+der_sumkmr+der_sumkqr;
        old_bar_nl=the_particle->bar_nl; old_bar_nr=the_particle->bar_nr;
        the_particle->bar_nl=the_particle->bar_nl-0.0001*func_l/der_f_l;
        the_particle->bar_nr=the_particle->bar_nr-0.0001*func_r/der_f_r;
        error=max(fabs(old_bar_nl-the_particle->bar_nl),fabs(old_bar_nr-the_particle->bar_nr))/0.0001;
    }
}

void System::solve_system_restore(double& bar_nl,double& bar_nr,double nr){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr,der_sumkml,der_sumkmr,der_sumkql,der_sumkqr;
    double block_l,block_r;
    double der_f_l,func_l,der_f_r,func_r;
    // double grad_r,grad_l;
    while(error>0.00000001){
        block_l=bar_nl*xi_l;
        block_r=bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*block_l*pow(block_l+1.,kl-1.);
        sumkql=exp(mu_linker-fcnf)*kl*block_l*pow(block_l+1.,kl-1.)*(pow(block_r+1.,kr)-1.);
        sumkmr=exp(mu_linker)*kr*block_r*pow(block_r+1.,kr-1.);
        sumkqr=exp(mu_linker-fcnf)*kr*block_r*pow(block_r+1.,kr-1.)*(pow(block_l+1.,kl)-1.);
        der_sumkml=sumkml/bar_nl+sumkml*(kl-1.)*xi_l/(block_l+1.);
        der_sumkql=sumkql/bar_nl+sumkql*(kl-1.)*xi_l/(block_l+1.);
        der_sumkmr=sumkmr/bar_nr+sumkmr*(kr-1.)*xi_r/(block_r+1.);
        der_sumkqr=sumkqr/bar_nr+sumkqr*(kr-1.)*xi_r/(block_r+1.);
        
        func_l=bar_nl+sumkml+sumkql-nl;
        der_f_l=1.+der_sumkml+der_sumkql;
        func_r=bar_nr+sumkmr+sumkqr-nr;
        der_f_r=1.+der_sumkmr+der_sumkqr;
        old_bar_nl=bar_nl; old_bar_nr=bar_nr;
        bar_nl=bar_nl-0.0001*func_l/der_f_l;
        bar_nr=bar_nr-0.0001*func_r/der_f_r;
        error=max(fabs(old_bar_nl-bar_nl),fabs(old_bar_nr-bar_nr))/0.0001;
    }
}

void System::picard_solve_system(Guest* the_particle){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr;
    double block_l,block_r;
    double func_l,func_r;
    while(error>0.0000001){
        block_l=the_particle->bar_nl*xi_l;
        block_r=the_particle->bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*xi_l*pow(block_l+1.,kl-1.);
        sumkql=exp(mu_linker-fcnf)*kl*xi_l*pow(block_l+1.,kl-1.)*(pow(block_r+1.,kr)-1.);
        sumkmr=exp(mu_linker)*kr*xi_r*pow(block_r+1.,kr-1.);
        sumkqr=exp(mu_linker-fcnf)*kr*xi_r*pow(block_r+1.,kr-1.)*(pow(block_l+1.,kl)-1.);
        
        func_l=1+(sumkml+sumkql);
        func_r=1+(sumkmr+sumkqr);

        old_bar_nl=the_particle->bar_nl; old_bar_nr=the_particle->bar_nr;
        the_particle->bar_nl=the_particle->bar_nl*0.99+0.01*the_particle->nl/func_l;
        the_particle->bar_nr=the_particle->bar_nr*0.99+0.01*the_particle->nr/func_r;
        error=max(fabs(old_bar_nl-the_particle->bar_nl)/0.01,fabs(old_bar_nr-the_particle->bar_nr)/0.01);
    }
}

void System::picard_solve_system_restore(double& bar_nl, double& bar_nr,double nr){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,sumkql,sumkqr;
    double block_l,block_r;
    double func_l,func_r;
    while(error>0.0000001){
        block_l=bar_nl*xi_l;
        block_r=bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*xi_l*pow(block_l+1.,kl-1.);
        sumkql=exp(mu_linker-fcnf)*kl*xi_l*pow(block_l+1.,kl-1.)*(pow(block_r+1.,kr)-1.);
        sumkmr=exp(mu_linker)*kr*xi_r*pow(block_r+1.,kr-1.);
        sumkqr=exp(mu_linker-fcnf)*kr*xi_r*pow(block_r+1.,kr-1.)*(pow(block_l+1.,kl)-1.);
        
        func_l=1+(sumkml+sumkql);
        func_r=1+(sumkmr+sumkqr);

        old_bar_nl=bar_nl; old_bar_nr=bar_nr;
        bar_nl=bar_nl*0.99+0.01*nl/func_l;
        bar_nr=bar_nr*0.99+0.01*nr/func_r;
        error=max(fabs(old_bar_nl-bar_nl)/0.01,fabs(old_bar_nr-bar_nr)/0.01);
    }
}


void System::solve_ref_system(Guest* the_particle){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,der_sumkml,der_sumkmr;
    double der_f_l,func_l,der_f_r,func_r;
    double block_l,block_r;
    while(error>0.00000001){

        block_l=the_particle->ref_bar_nl*xi_l;
        block_r=the_particle->ref_bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*block_l*pow(block_l+1,kl-1);
        sumkmr=exp(mu_linker)*kr*block_r*pow(block_r+1,kr-1);
        der_sumkml=sumkml/the_particle->ref_bar_nl+sumkml*(kl-1)*xi_l/(block_l+1);
        der_sumkmr=sumkmr/the_particle->ref_bar_nr+sumkmr*(kr-1)*xi_r/(block_r+1);

        func_l=the_particle->ref_bar_nl+sumkml-the_particle->nl;
        der_f_l=1+der_sumkml;
        func_r=the_particle->ref_bar_nr+sumkmr-the_particle->nr;
        der_f_r=1+der_sumkmr;
        old_bar_nl=the_particle->ref_bar_nl; old_bar_nr=the_particle->ref_bar_nr;
        the_particle->ref_bar_nl=the_particle->ref_bar_nl-0.0001*func_l/der_f_l;
        the_particle->ref_bar_nr=the_particle->ref_bar_nr-0.0001*func_r/der_f_r;
        error=max(fabs(old_bar_nl-the_particle->ref_bar_nl),fabs(old_bar_nr-the_particle->ref_bar_nr))/0.0001;
    }
}
void System::solve_ref_system_restore(double& ref_bar_nl,double& ref_bar_nr,double nr){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr,der_sumkml,der_sumkmr;
    double der_f_l,func_l,der_f_r,func_r;
    double block_l,block_r;
    while(error>0.00000001){

        block_l=ref_bar_nl*xi_l;
        block_r=ref_bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*block_l*pow(block_l+1,kl-1);
        sumkmr=exp(mu_linker)*kr*block_r*pow(block_r+1,kr-1);
        der_sumkml=sumkml/ref_bar_nl+sumkml*(kl-1)*xi_l/(block_l+1);
        der_sumkmr=sumkmr/ref_bar_nr+sumkmr*(kr-1)*xi_r/(block_r+1);

        func_l=ref_bar_nl+sumkml-nl;
        der_f_l=1+der_sumkml;
        func_r=ref_bar_nr+sumkmr-nr;
        der_f_r=1+der_sumkmr;
        old_bar_nl=ref_bar_nl; old_bar_nr=ref_bar_nr;
        ref_bar_nl=ref_bar_nl-0.0001*func_l/der_f_l;
        ref_bar_nr=ref_bar_nr-0.0001*func_r/der_f_r;
        error=max(fabs(old_bar_nl-ref_bar_nl),fabs(old_bar_nr-ref_bar_nr));
    }
}

void System::picard_solve_ref_system(Guest* the_particle){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr;
    double func_l,func_r;
    double block_l,block_r;
    while(error>0.00000001){

        block_l=the_particle->ref_bar_nl*xi_l;
        block_r=the_particle->ref_bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*xi_l*pow(block_l+1,kl-1);
        sumkmr=exp(mu_linker)*kr*xi_r*pow(block_r+1,kr-1);

        func_l=1+sumkml;
        func_r=1+sumkmr;
        
        old_bar_nl=the_particle->ref_bar_nl; old_bar_nr=the_particle->ref_bar_nr;
        the_particle->ref_bar_nl=the_particle->ref_bar_nl*0.95+0.05*the_particle->nl/func_l;
        the_particle->ref_bar_nr=the_particle->ref_bar_nr*0.95+0.05*the_particle->nr/func_r;
        error=max(fabs(old_bar_nl-the_particle->ref_bar_nl)/0.05,fabs(old_bar_nr-the_particle->ref_bar_nr)/0.05);
    }
}

void System::picard_solve_ref_system_restore(double& ref_bar_nl, double& ref_bar_nr,double nr){
    double old_bar_nl, old_bar_nr;
    double error=1;
    double sumkml,sumkmr;
    double func_l,func_r;
    double block_l,block_r;
    while(error>0.00000001){

        block_l=ref_bar_nl*xi_l;
        block_r=ref_bar_nr*xi_r;
        sumkml=exp(mu_linker)*kl*xi_l*pow(block_l+1,kl-1);
        sumkmr=exp(mu_linker)*kr*xi_r*pow(block_r+1,kr-1);

        func_l=1+sumkml;
        func_r=1+sumkmr;
        
        old_bar_nl=ref_bar_nl; old_bar_nr=ref_bar_nr;
        ref_bar_nl=ref_bar_nl*0.99+0.01*nl/func_l;
        ref_bar_nr=ref_bar_nr*0.99+0.01*nr/func_r;
        error=max(fabs(old_bar_nl-ref_bar_nl)/0.01,fabs(old_bar_nr-ref_bar_nr)/0.01);
    }
}

void System::cal_Feff(Guest* the_particle){
    double summl,summr;
    double sumb;
    double block_l,block_r;
    block_l=pow(the_particle->bar_nl*xi_l+1,kl);
    block_r=pow(the_particle->bar_nr*xi_r+1,kr);
    summl=exp(mu_linker)*(block_l-1);
    summr=exp(mu_linker)*(block_r-1);
    sumb=exp(mu_linker-fcnf)*(block_l*block_r-block_l-block_r+1);
    the_particle->Feff= the_particle->nl*log(the_particle->bar_nl)-the_particle->bar_nl-summl
                +the_particle->nr*log(the_particle->bar_nr)-the_particle->bar_nr-summr
                -sumb;
    //set a reference
    double ref_summl,ref_summr;
    double ref_block_l,ref_block_r;
    ref_block_l=pow(the_particle->ref_bar_nl*xi_l+1,kl);
    ref_block_r=pow(the_particle->ref_bar_nr*xi_r+1,kr);
    ref_summl=exp(mu_linker)*(ref_block_l-1);
    ref_summr=exp(mu_linker)*(ref_block_r-1);

    the_particle->Feff=the_particle->Feff-the_particle->nl*log(the_particle->ref_bar_nl)+the_particle->ref_bar_nl+ref_summl
                -the_particle->nr*log(the_particle->ref_bar_nr)+the_particle->ref_bar_nr+ref_summr;

}

void System::update_q(Guest* the_particle){
    check_for_Feff(the_particle);
    the_particle->q=exp(-the_particle->Feff)-1.;
}

void System::extract_particle(Guest* the_particle){
    the_particle->backup();
    int cell_x,cell_y;
    cell_x=floor(the_particle->boxs[0]/cell_length);
    cell_y=floor(the_particle->boxs[1]/cell_length);
    cells[cell_x][cell_y].remove(the_particle);
}

void System::put_particle(Guest* the_particle){
    int cell_x,cell_y;
    cell_x=floor(the_particle->boxs[0]/cell_length);
    cell_y=floor(the_particle->boxs[1]/cell_length);
    cells[cell_x][cell_y].insert(the_particle);
}

void System::trans_move(void){
    Nsteps++;
    Guest* chosen_part=randomLink(guest_particles.headPtr,Ng);
    int collision;
    // chosen_part->backup();
    extract_particle(chosen_part);
    chosen_part->random_move(ds);
    put_particle(chosen_part);
    collision=check_overlap(chosen_part);
    if(collision ==1){
        int cell_x,cell_y;
        cell_x=floor(chosen_part->boxs[0]/cell_length);
        cell_y=floor(chosen_part->boxs[1]/cell_length);
        cells[cell_x][cell_y].remove(chosen_part);
        chosen_part->reject();
        put_particle(chosen_part);
    }else{
        count_nr(chosen_part);
        update_q(chosen_part);
        if(randCV()<chosen_part->q/chosen_part->backup_q){
            accSteps++;
        }else{
            int cell_x,cell_y;
            cell_x=floor(chosen_part->boxs[0]/cell_length);
            cell_y=floor(chosen_part->boxs[1]/cell_length);
            cells[cell_x][cell_y].remove(chosen_part);
            chosen_part->reject();
            put_particle(chosen_part);
        }
    }
}

void System::add_particles(void){
    Ng=Ng+1;

    Guest *newInserted;

    newInserted= new Guest;
    newInserted->initialize_guest(xi_l,nl);
    guest_particles.insert(newInserted);

    put_particle(newInserted);
    int collision=check_overlap(newInserted);
    if(collision == 1){
        guest_particles.remove(newInserted);
        extract_particle(newInserted);
        delete newInserted;
        Ng=Ng-1;
    }else{
        count_nr(newInserted);
        newInserted->bar_nr=newInserted->nr/2.;
        newInserted->ref_bar_nr=newInserted->bar_nr;
        update_q(newInserted);
        if(randCV()<rec_substrate.area*newInserted->q*exp(mu_g)/(Ng)){
            
        }else{
            guest_particles.remove(newInserted);
            extract_particle(newInserted);
            delete newInserted;
            Ng=Ng-1;
        }
    }
}

void System::remove_particles(void){

    Guest* here;

    if(Ng==0){
        //cout<<"No particle in box"<<endl;
    }else{
        here=randomLink(guest_particles.headPtr,Ng);
        guest_particles.remove(here);
        extract_particle(here);
    
        if(randCV()<Ng*exp(-1*mu_g)/(rec_substrate.area*here->q)){
            delete here;
            Ng=Ng-1;
        }else{
            guest_particles.insert(here);
            put_particle(here);
        }
    }
}

void System::MCmove(void){
    double indicator=randCV();
    Guest* here;
    if(indicator<0.2 && Ng>0){
        here=randomLink(guest_particles.headPtr,Ng);
        trans_move();
        adapt_step();
    }else{
        // indicator=randCV();
        if(indicator<0.6){
            remove_particles();
        }else{
            add_particles();
        }
    }
}

void System::adapt_step(void){
    if(Nsteps>200){
        double acceptance_rate;
        acceptance_rate=accSteps/Nsteps;
        if(acceptance_rate>0.6 && ds<0.09){
            ds=ds*1.01;
        }else if(acceptance_rate<0.3 && ds>0.00005){
            ds=ds/1.01;
        }
        Nsteps=0;accSteps=0;
    }
}

void System::cal_alpha(void){

}