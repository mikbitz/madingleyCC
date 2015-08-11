#ifndef NONSTATICRNG_H
#define NONSTATICRNG_H
/** \file NonStaticRNG.h
 * \brief the NonStaticRNG header file
 */





//
//namespace Madingley
//{
/** \brief
//    /// SimpleRNG is a simple random number generator 
//    /// </summary>
//     class NonStaticRNG
//    {
//
//        private Random rand;
//
/** \brief
Constructor the random number generator
*/
//         NonStaticRNG()
//        {
//
//            rand = new Random();
//
//            
//        }
//
//        // The random generator seed can be set two ways:
//        // 1) specifying one non-zero unsigned integer
//        // 2) setting the seed from the system time
//
/** \brief
Set the seed of the random number generator using the specified integer
*/
@param u An integer 
//         void SetSeed(unsigned u)
//        {
//            rand = new Random((int)u);
//        }
//
/** \brief
Set the seed of the random number generator based on the system time
*/
//         void SetSeedFromSystemTime()
//        {
//            System.DateTime dt = System.DateTime.Now;
//            long x = dt.ToFileTime();
//            
//            rand = new Random((int)x);
//        }
//
/** \brief
A random draw from a uniform distribution between 0 and 1
*/
@return A random draw from a uniform distribution between 0 and 1*/
<remarks>This will not return either 0 or 1</remarks>
//         double GetUniform()
//        {
//            return (rand.NextDouble());
//
//        }
//
/** \brief
Get a random unsigned integer
*/
@return A random unsigned integer*/
//        private unsigned GetUint()
//        {
//            return ((unsigned)rand.Next());
//        }
//
/** \brief
A random draw from a normal distribution with mean 0 and standard deviation 1
*/
@return A random draw from a normal distribution with mean 0 and standard deviation 1*/
//         double GetNormal()
//        {
//            // Use Box-Muller algorithm
//            double u1 = GetUniform();
//            double u2 = GetUniform();
//            double r = Math.Sqrt(-2.0 * Math.Log(u1));
//            double theta = 2.0 * Math.PI * u2;
//            return r * Math.Sin(theta);
//
//            
//        }
//
/** \brief
A random draw from a normal distribution
*/
@param mean The mean of the normal distribution 
@param standardDeviation The standard deviation of the normal distribution 
@return A random draw from a normal distribution*/
//         double GetNormal(double mean, double standardDeviation)
//        {
//            if (standardDeviation <= 0.0)
//            {
//                string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
//                throw new ArgumentOutOfRangeException(msg);
//            }
//            return mean + standardDeviation * GetNormal();
//        }
//
/** \brief
A random draw from an exponential distribution with mean 1
*/
@return A random draw from an exponential distribution with mean 1*/
//         double GetExponential()
//        {
//            return -Math.Log(GetUniform());
//        }
//
/** \brief
A random draw from the exponential distribution
*/
@param mean The mean of the exponential distribution 
@return A random draw from the exponential distribution*/
//         double GetExponential(double mean)
//        {
//            if (mean <= 0.0)
//            {
//                string msg = string.Format("Mean must be positive. Received {0}.", mean);
//                throw new ArgumentOutOfRangeException(msg);
//            }
//            return mean * GetExponential();
//        }
//
//
/** \brief
A random draw from a lognormal distribution
*/
@param mu Mean of the lognormal distribution 
@param sigma Standard deviation of the lognormal distribution 
@return A random draw from a lognormal distribution*/
//         double GetLogNormal(double mu, double sigma)
//        {
//            return Math.Exp(GetNormal(mu, sigma));
//        }
//
//    }
//}
#endif
