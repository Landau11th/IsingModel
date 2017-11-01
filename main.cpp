//io and other general headers
#include <iostream>
#include <ctime>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <fstream>

//to process arguments from command line
#include <vector>
#include <algorithm>
#include <string>
#include <map>

//OpenMP
#include <omp.h>

//My Modelings
#include "rand64.hpp"
#include "Ising.hpp"


class Comd_Line_Args
{
public:
	//for OpenMP
	unsigned int num_threads;
	//Monte Carlo sampling size
	unsigned long long int num_samples;//long long int is definitely enough
	//Specs of Ising model
	double kBT_over_J_min = 2;
	double kBT_over_J_max = 2.5;
	unsigned int num_kBT = 100;
	unsigned int grid_width;
	unsigned int grid_height;
	//record all the options and status
	std::map<std::string, bool> opt_status;
	//Constructor and public function
	Comd_Line_Args();
	void Input_Arguments(int argc_main, char *argv_main[]);
	bool Start_Rest();
private:
	void Treat_single_option(const std::string option_flag, const std::string option_value);

};
Comd_Line_Args::Comd_Line_Args()
{
	opt_status["--gridsize="] = false;
	opt_status["--kBT="] = false;
	opt_status["--samplingsize="] = false;
	opt_status["--threads="] = false;

	//-h is special
	opt_status["-h"] = true;
}
void Comd_Line_Args::Input_Arguments(int argc_main, char *argv_main[])
{
	std::string option;
	auto iter = opt_status.begin();

	std::string option_flag;
	std::string option_value;
	for (int i = 1; i < argc_main; ++i)
	{
		option = argv_main[i];
		for (iter = opt_status.begin(); iter != opt_status.end(); ++iter)
		{
			if (option.find(iter->first) == 0)
			{
				//split the string to option flag and value parts, and pass to other member function
				option_flag = iter->first;
				option_value = option.substr(option_flag.size());
				//std::cout << option_value << std::endl;
				//guarantees
				Treat_single_option(option_flag, option_value);
				break;
			}
		}
		if (iter == opt_status.end())
		{
			std::cout << "Unknown option" << option << std::endl;
		}
	}
}
void Comd_Line_Args::Treat_single_option(const std::string option_flag, const std::string option_value)
{

	if (option_flag == "--gridsize=")
	{
		size_t pos = option_value.find('*');
		std::string temp = option_value.substr(0, pos);
		grid_width = std::stoi(temp);
		temp = option_value.substr(pos+1);
		grid_height = std::stoi(temp);
	}
	else if (option_flag == "--samplingsize=")
	{
		double temp = std::stod(option_value);
		num_samples = (int)(temp + 0.5);
	}
	else if (option_flag == "--threads=")
    {
        num_threads = std::stoi(option_value);
    }
    else if (option_flag == "--kBT=")
    {
        size_t pos1 = option_value.find('~');
        size_t pos2 = option_value.find(',');
        kBT_over_J_min = std::stod(option_value.substr(0, pos1));
        kBT_over_J_max = std::stod(option_value.substr(pos1+1, pos2-pos1));
        num_kBT = std::stoi(option_value.substr(pos2+1));
    }

	opt_status[option_flag] = true;

	//if someone need help, we do not start the rest of the code
	//should be realized in exception and catch. Will revise later.
	if (option_flag == "-h")
	{
		std::cout << "Must give following parameters" << std::endl;
		std::cout << "--gridsize=X*Y" << std::endl;
		std::cout << "--samplingsize=X" << std::endl;
		std::cout << "--threads=X" << std::endl;
		std::cout << "--kBT=A~B,X" << std::endl;
		std::cout << "where X, Y are positive integer, A, B are positive double" << std::endl;
		std::cout << "Example of 5 by 5 lattice, temperature varies from 1.5 to 2.5 divided to 1000 parts, using 4 threads and 100000 samples:" << std::endl;
		std::cout << "--gridsize=5*5 --kBT=1.5~2.5,1000 --threads=4 --samplingsize=100000" << std::endl;
		//terminate();
		opt_status["-h"] = false;
	}
}
bool Comd_Line_Args::Start_Rest()
{
	if(!opt_status["-h"])
    {
        return false;
    }
	for (auto iter = opt_status.begin(); iter != opt_status.end(); ++iter)
	{
		//if some status is false
		if (!(iter->second))
		{
			std::cout << "Lack argument for " << iter->first << std::endl;
			return false;
		}
	}
	return true;
}


int main(int argc, char *argv[])
{
	//parse the input command line arguments to get necessary variables
	Comd_Line_Args args;
    args.Input_Arguments(argc, argv);
    //check if all required parameters are received
	if (!args.Start_Rest())
	{
		std::cout << "Calculation didn't start" << std::endl;
		return 0;
	}
    else
    {
        std::cout << "Calculation started" << std::endl;
    }


    assert(args.num_threads>0 && args.num_threads<=omp_get_max_threads() && "# of Threads out of range!\n" );

    //create output file name with time stamp, to avoid being covered
    std::string filename = "stats_";
    {
        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        char buffer[80];
        strftime(buffer,80,"%Y%m%d-%H%M%S",now);
        filename = filename+buffer+".dat";
    }
    //write important parameters at the beginning of the file
    std::ofstream outputfile;
    outputfile.open(filename);
    assert(outputfile.is_open() && "File not open!");
    outputfile << "# of threads used: " << args.num_threads << std::endl;
    outputfile << "# of samplings: " << args.num_samples << std::endl;
    outputfile << "gridsize: " << args.grid_width << "*" << args.grid_height << std::endl;
    outputfile << "Scaled T varies from " << args.kBT_over_J_min << " to " << args.kBT_over_J_max << ", divided to " << args.num_kBT << " parts" << std::endl << std::endl;




    //random seed, will be modified later for all threads
    int seed = time(nullptr);

    //for elapsed time
    double start_time = omp_get_wtime();


    double **moment_1st_all_threads = new double*[args.num_threads];
    double **moment_2nd_all_threads = new double*[args.num_threads];
    double **moment_4th_all_threads = new double*[args.num_threads];

    //double *kBTs = new double[args.num_kBT + 1];
    for(unsigned int i = 0; i < args.num_threads; ++i)
    {
        moment_1st_all_threads[i] = new double[args.num_kBT + 1];
        moment_2nd_all_threads[i] = new double[args.num_kBT + 1];
        moment_4th_all_threads[i] = new double[args.num_kBT + 1];
    }


    #pragma omp parallel num_threads(args.num_threads)
    {
        //may help improve the false sharing?
        //int avoid_cache_problem[32];

        //get omp configurations and choose a random number
        const unsigned int id = omp_get_thread_num();
        const unsigned int num_thrds = omp_get_num_threads();

        Deng::RandNumGen::LCG64 rand_this_thrd(seed + id*218459);
        //set up Ising model
        //notice Ising_2D has an pointer member linking to the RNG
        Deng::Monte_Carlo::Ising::Ising_2D Ising(args.grid_width, args.grid_height, 0.0, &rand_this_thrd);
        Ising.Set_Spinors_Neighbours();
        Ising.Calc_Total_Energy();
        //acquire necessary constants
        const unsigned long long int num_samples = args.num_samples;
        const double kBT_over_J_min = args.kBT_over_J_min;
        const double kBT_over_J_max = args.kBT_over_J_max;
        const unsigned int num_kBT = args.num_kBT;
        const double incre_kBT = (kBT_over_J_max - kBT_over_J_min)/num_kBT;

        double current_kBT_over_J = 0;
        //for each kBT do the same operation
        for(unsigned int count_kBT = 0; count_kBT <= num_kBT; ++count_kBT)
        {
            current_kBT_over_J = count_kBT*incre_kBT + kBT_over_J_min;

            //check if kBT is illegal
            if(current_kBT_over_J < 1E-4)
            {
                if(id == 0)
                {

                    std::cout << "Scaled T (" << current_kBT_over_J << ") is too low, skipped current iteration" << std::endl;
                }
                continue;
            }



            Ising.Modify_kBT(current_kBT_over_J);
//            //before doing stats, make the configuration more messy
//            for(int i = 0; i < 100; ++i)
//            {
//                Ising.Flip();
//            }

            double temp_energy = 0;
            double moment_1st = 0.0;
            double moment_2nd = 0.0;
            double moment_4th = 0.0;

            for(unsigned long long int sample_count = id; sample_count < num_samples; sample_count += num_thrds)
            {
                //Flip
                Ising.Flip();
                //record current energy and do statistics
                temp_energy = Ising.Energy();
                moment_1st += temp_energy;
                moment_2nd += temp_energy*temp_energy;
                moment_4th += pow(temp_energy, 4);
            }

            moment_1st_all_threads[id][count_kBT] = moment_1st;
            moment_2nd_all_threads[id][count_kBT] = moment_2nd;
            moment_4th_all_threads[id][count_kBT] = moment_4th;


            if(id == 0)
                std::cout << "Finished kBT=" << current_kBT_over_J << std::endl;
        }


    }
    std::cout << "Simulation costs " << omp_get_wtime() - start_time << "s" << std::endl;
    outputfile << "Simulation costs " << omp_get_wtime() - start_time << "s" << std::endl << std::endl;

    //output to file
    outputfile << "scaled T, " << "average Energy, " << "specific heat, "<< "error of spec heat, " << std::endl;
    double moment_1st, moment_2nd, moment_4th;
    const double total_sites = args.num_samples * args.grid_width * args.grid_height;
    for(unsigned int count_kBT = 0; count_kBT <= args.num_kBT; ++count_kBT)
    {
        moment_1st = moment_2nd = moment_4th = 0.0;
        for(int i = 0; i < args.num_threads; ++i)
        {
            moment_1st += moment_1st_all_threads[i][count_kBT];
            moment_2nd += moment_2nd_all_threads[i][count_kBT];
            moment_4th += moment_4th_all_threads[i][count_kBT];
        }

        moment_1st = moment_1st/args.num_samples;
        moment_2nd = moment_2nd/args.num_samples;
        moment_4th = moment_4th/args.num_samples;
        moment_4th = sqrt(moment_4th - moment_2nd*moment_2nd)/(args.grid_width * args.grid_height);
        moment_2nd = (moment_2nd - moment_1st*moment_1st)/(args.grid_width * args.grid_height);



        outputfile << args.kBT_over_J_min + count_kBT*(args.kBT_over_J_max - args.kBT_over_J_min)/args.num_kBT << ", " << moment_1st << ", " << moment_2nd << std::endl;
    }
    //2.269185314213 for square lattice

    std::cout << "With statistics costs " << omp_get_wtime() - start_time << "s" << std::endl;


    return 0;
}
