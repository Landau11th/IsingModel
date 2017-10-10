#include <iostream>
#include <ctime>

#include <omp.h>

#include "rand64.hpp"
#include "Ising.hpp"

int main()
{
//    time_t now = time(nullptr);
//    int now_int = time(nullptr);
//
//    std::cout << now << std::endl;
//    std::cout << now_int << std::endl;

    Deng::RandNumGen::LCG64 my_rand(time(nullptr));

    Deng::Monte_Carlo::Ising::Ising_2D Ising(10, 10, 1.0, &my_rand);

    Ising.Set_Spinors_Neighbours();
    Ising.Calc_Total_Energy();

    Ising.Flip();

    return 0;
}
