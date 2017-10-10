#ifndef DENG_RAND64_HPP
#define DENG_RAND64_HPP

#include <assert.h>
//#include <stdio.h>

//These codes are modified from Prof Wang Jian-sheng's random number generators

namespace Deng
{
    namespace RandNumGen
    {
        class RNG_Base
        {
        public:
            //as a function object
            virtual double operator()() { return drand64(); };
            //generate random number
            virtual double drand64() = 0;
            //input random seed
            virtual void srand64(int seed) = 0;
            //virtual void srand64(int seed, std::ofstream & output_file) {};//in case we need to output the seed to file
            virtual ~RNG_Base() = default;
        };

        class LCG64 : public RNG_Base
        {
        protected:
            //all generator of type LCG64 share the same x!!!
            static unsigned long long int x;
        public:
            //require the user to call constructor explicitly to initialize the seed
            LCG64(int seed) { srand64(seed); };
            virtual double drand64() override;
            virtual void srand64(int seed) override;
        };
    }
}

#endif // DENG_RAND64_HPP
