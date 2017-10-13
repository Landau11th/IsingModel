#include "Ising.hpp"

using namespace Deng::Monte_Carlo::Ising;

Site_2D::Site_2D()
{
    _N_spin = 0;
    _N_neighbour = 0;
    neighbour = nullptr;
}

Site_2D::~Site_2D()
{
    for(unsigned int i = 0; i < _N_neighbour; ++i)
    {
        delete[] neighbour[i];
    }
    delete[] neighbour;
}

void Site_2D::Set_Spinor(unsigned int N_spin, unsigned int N_neighbour)
{
    _N_spin = N_spin;
    _N_neighbour = N_neighbour;

    //pointing up, with max spin
    spin = (_N_spin - 1);

    neighbour = new int*[_N_neighbour];
    for(unsigned int i = 0; i < _N_neighbour; ++i)
    {
        //2 means we only need 2 index to determine the neighbour
        //In fact one is always enough
        neighbour[i] = new int[2];
    }
}


Ising_2D::Ising_2D(unsigned int width, unsigned int height, double kBT_over_J, Deng::RandNumGen::RNG_Base* RNG_pt, double J)
{
    _width = width;
    _height = height;
    _kBT_over_J = kBT_over_J;
    _J = J;

    RandNum_pt = RNG_pt;

    sites = new Site_2D*[height];
    for(int i = 0; i < _height; ++i)
    {
        sites[i] = new Site_2D[_width];
    }
}
Ising_2D::~Ising_2D()
{
    //we do not release RandNum_pt here because the lifetime of the RNG could be independent with our Ising model
    for(int i = 0; i < _height; ++i)
    {
        delete[] sites[i];
    }
    delete[] sites;
}
void Ising_2D::Set_Spinors_Neighbours()
{
    for(int i = 0; i < _height; ++i)
    {
        for(int j = 0; j < _width; ++j)
        {
            sites[i][j].Set_Spinor(2, 4);

            //nearest sites
            //neighbour on top
            sites[i][j].neighbour[0][0] = ((i-1)<0) ? (i-1+_height):(i-1);
            sites[i][j].neighbour[0][1] = j;

            //bottom
            sites[i][j].neighbour[1][0] = ((i+1)>=_height) ? (i+1-_height):(i+1);
            sites[i][j].neighbour[1][1] = j;

            //left
            sites[i][j].neighbour[2][0] = i;
            sites[i][j].neighbour[2][1] = ((j-1)<0) ? (j-1+_width):(j-1);

            //right
            sites[i][j].neighbour[3][0] = i;
            sites[i][j].neighbour[3][1] = ((j+1)>=_width) ? (j+1-_width):(j+1);
        }
    }
}
void Ising_2D::Calc_Total_Energy()
{
    _energy = 0.0;

    for(int i = 0; i < _height; ++i)
    {
        for(int j = 0; j < _width; ++j)
        {
            _energy += Calc_One_Site_Single_Site_Energy(i, j);
            _energy += 0.5*Calc_One_Site_2_Sites_Interaction(i,j);
        }
    }
}
double Ising_2D::Calc_One_Site_2_Sites_Interaction(int i, int j)
{
    double E = 0.0;

    assert( i >= 0 && i <= _height && j >= 0 && j <= _width);

    const unsigned int N_neighbour = sites[i][j].N_neighbour();

    for(unsigned int count_neighbour = 0; count_neighbour < N_neighbour; ++count_neighbour)
    {
        E += -_J * sites[i][j].spin * sites[sites[i][j].neighbour[count_neighbour][0]][sites[i][j].neighbour[count_neighbour][1]].spin;
        //E += -sites[i][j].spin * sites[sites[i][j].neighbour[count_neighbour][0]][sites[i][j].neighbour[count_neighbour][1]].spin;
    }

    return E;
}
double Ising_2D::Calc_One_Site_Single_Site_Energy(int i, int j)
{
    return 0.0;
}
void Ising_2D::Flip()
{
    //randomly choose the site to be flipped
    const unsigned int i = _height*RandNum_pt->drand64();
    const unsigned int j = _width*RandNum_pt->drand64();


    const unsigned int N_spin = sites[i][j].N_spin();
    const int spin = sites[i][j].spin;
    unsigned int n_th_level = spin + N_spin - 1;
    //double check whether the spin is allowed
    //such function should be included in Site_2D!!!!!!!! will figure out later
    assert((n_th_level%2)==0);

    n_th_level = n_th_level/2;
    //uniformly randomly flip the spin to another spin state (in general not necessarily -spin)
    n_th_level = n_th_level + (int)(RandNum_pt->drand64()*(N_spin - 1)) + 1;
    n_th_level = n_th_level%N_spin;

    //calculate the energy difference
    double energy_before_flip = Calc_One_Site_2_Sites_Interaction(i, j) + Calc_One_Site_Single_Site_Energy(i, j);

    //flip the spin!
    sites[i][j].spin = 2*n_th_level - (N_spin - 1);
    double energy_after_flip = Calc_One_Site_2_Sites_Interaction(i, j) + Calc_One_Site_Single_Site_Energy(i, j);

    double dE = energy_after_flip - energy_before_flip;

    if(RandNum_pt->drand64()<exp(-dE/_kBT_over_J))//if dE<0, or drand64<probability,
    {
        //spin is already flipped, therefore we only need to update the energy
        _energy = _energy + dE;//update the total energy
    }
    else
    {
        sites[i][j].spin = spin;
    }


}
