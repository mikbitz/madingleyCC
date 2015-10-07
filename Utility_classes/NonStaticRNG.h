#ifndef NONSTATICRNG_H
#define NONSTATICRNG_H
/** \file NonStaticRNG.h
 * \brief the NonStaticRNG header file
 */
      
        
/** \brief SimpleRNG is a simple random number generator */

     class NonStaticRNG
    {

    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;

/** \brief Constructor the random number generator
*/
         NonStaticRNG()
        {

        }

        // The random generator seed can be set two ways:
        // 1) specifying one non-zero unsigned integer
        // 2) setting the seed from the system time

/** \brief Set the seed of the random number generator using the specified integer
*/
@param u An integer 
         void SetSeed(unsigned u)
        {
            RandomNumberGenerator.seed(u);
        }

/** \brief
Set the seed of the random number generator based on the system time
*/
         void SetSeedFromSystemTime()
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            RandomNumberGenerator.seed(seed);
        }

/** \brief A random draw from a uniform distribution between 0 and 1

@return A random draw from a uniform distribution between 0 and 1
@remarksThis will not return either 0 or 1 */
         double GetUniform()
        {
          std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
          return randomNumber(RandomNumberGenerator);

        }

/** \brief Get a random unsigned integer
@return A random unsigned integer*/
        unsigned GetUint()
        {
            return std::default_random_engine(seed);
        }

/** \brief A random draw from a normal distribution with mean 0 and standard deviation 1
@return A random draw from a normal distribution with mean 0 and standard deviation 1*/
         double GetNormal()
        {
            // Use Box-Muller algorithm
            double u1 = GetUniform();
            double u2 = GetUniform();
            double r = sqrt(-2.0 * log(u1));
            double theta = 2.0 * acos(-1.) * u2;
            return r * sin(theta);

            
        }

/** \brief A random draw from a normal distribution

@param mean The mean of the normal distribution 
@param standardDeviation The standard deviation of the normal distribution 
@return A random draw from a normal distribution*/
         double GetNormal(double mean, double standardDeviation)
        {
            if (standardDeviation <= 0.0)
            {
                cout<<"Shape must be positive. Received "<< standardDeviation<<endl;
                exit(1);
            }
            return mean + standardDeviation * GetNormal();
        }

/** \brief A random draw from an exponential distribution with mean 1
@return A random draw from an exponential distribution with mean 1*/
         double GetExponential()
        {
            return -log(GetUniform());
        }

/** \brief
A random draw from the exponential distribution
*/
@param mean The mean of the exponential distribution 
@return A random draw from the exponential distribution*/
         double GetExponential(double mean)
        {
            if (mean <= 0.0)
            {
                cout<<"Mean must be positive. Received "<< mean<<endl;
                exit(1);
            }
            return mean * GetExponential();
        }


/** \brief
A random draw from a lognormal distribution
*/
@param mu Mean of the lognormal distribution 
@param sigma Standard deviation of the lognormal distribution 
@return A random draw from a lognormal distribution*/
         double GetLogNormal(double mu, double sigma)
        {
            return exp(GetNormal(mu, sigma));
        }

    }
}
     };
#endif
