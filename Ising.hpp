#ifndef DENG_ISING_HPP
#define DENG_ISING_HPP

#include <iostream>
#include <cmath>

#include "rand64.hpp"

namespace Deng
{
    namespace Monte_Carlo
    {
        namespace Ising
        {
            class Site_2D
            {
            protected:
                unsigned int _N_spin;
                unsigned int _N_neighbour;
            public:
                int spin;//spin direction
                int** neighbour;//nearest neighbors
                Site_2D();
                virtual ~Site_2D();
                unsigned int N_spin() const { return _N_spin; };
                unsigned int N_neighbour() const { return _N_neighbour; };
                virtual void Set_Spinor(unsigned int N_spin, unsigned int N_neighbour);
            };

            class Ising_2D
            {
            protected:
                unsigned int _width;//size of the network
                unsigned int _height;//size of the network
                double _kBT_over_J;//temperature scaled by interaction
                double _J = 1.0;//This quantity is useless since we are doing dimensionless calculation
                //But we still put it here in case for further modification
                double _energy;
                Deng::RandNumGen::RNG_Base* RandNum_pt;
            public:
                Site_2D** sites;
                Ising_2D(unsigned int width, unsigned int height, double kBT_over_J, Deng::RandNumGen::RNG_Base* RNG_pt, double J = 1.0);
                virtual ~Ising_2D();
                double Energy() const { return _energy; };
                void Modify_kBT (double kBT_over_J) { _kBT_over_J = kBT_over_J; };
                //set the type of spinors
                //this version is for spin half square lattice, should be overrode for other modeling
                virtual void Set_Spinors_Neighbours();


                //calculate total energy. In principle should be called only once
                //this version is for nearest neighbourhood, should be overrode for long range interaction, e.g. Heisenberg model
                virtual void Calc_Total_Energy();
            protected:
                //kernel of above function
                //I choose to hide it since we don't really have the chance to use it explicitly. could also be public
                virtual double Calc_Single_Site_Energy(int i, int j);

            public:
                //Flip returns the energy change
                //this version is for single flip, should be overrode, for example cluster flip
                virtual void Flip();

            };
        }
    }
}



#endif // DENG_ISING_HPP
