/**
 * FWaveVec.h
 *
 ****
 **** This is a vectorizable C++ implementation of the F-Wave solver (FWave.hpp).
 ****
 *
 * Created on: Nov 13, 2012
 * Last Update: Dec 28, 2013
 *
 ****
 *
 *  Author: Sebastian Rettenberger
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.
 *    E-Mail: rettenbs AT in.tum.de
 *  Some optimzations: Michael Bader
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Michael_Bader
 *    E-Mail: bader AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @article{bale2002wave,
 *            title={A wave propagation method for conservation laws and balance laws with spatially varying flux functions},
 *            author={Bale, D.S. and LeVeque, R.J. and Mitran, S. and Rossmanith, J.A.},
 *            journal={SIAM Journal on Scientific Computing},
 *            volume={24},
 *            number={3},
 *            pages={955--978},
 *            year={2002}}
 *
 *   @book{leveque2002finite,
 *         Author = {LeVeque, R. J.},
 *         Publisher = {Cambridge University Press},
 *         Title = {Finite Volume Methods for Hyperbolic Problems},
 *         Volume = {31},
 *         Year = {2002}}
 *
 *   @webpage{levequeclawpack,
 *            Author = {LeVeque, R. J.},
 *            Lastchecked = {January, 05, 2011},
 *            Title = {Clawpack Sofware},
 *            Url = {https://github.com/clawpack/clawpack-4.x/blob/master/geoclaw/2d/lib}}
 *
 ****
 */

#ifndef FWAVEVEC_HPP_
#define FWAVEVEC_HPP_

#include <cmath>

namespace solver
{

/**
 *
 */
template<typename T>
class FWaveVec
{
private:
	const T dryTol;
	const T half_gravity; // 0.5 * gravity constant
	const T sqrt_gravity; // square root of the gravity constant 
	const T zeroTol;

public:
	/**
	 * FWaveVec Constructor, takes three problem parameters
	 * @param dryTol "dry tolerance": if the water height falls below dryTol, wall boundary conditions are applied (default value is 100)
	 * @param gravity takes the value of the gravity constant (default value is 9.81 m/s^2)
	 * @param zeroTol computed f-waves with an absolute value < zeroTol are treated as static waves (default value is 10^{-7})
	 */
	FWaveVec(T i_dryTol = (T) 1.0,
			 T i_gravity = (T) 9.81,
			 T i_zeroTol = (T) 0.0000001)
		: dryTol(i_dryTol),
		  half_gravity( (T).5 * i_gravity ),
		  sqrt_gravity( std::sqrt(i_gravity) ),
		  zeroTol(i_zeroTol)
	{
	}

	/**
	 * takes the water height, discharge and bathymatry in the left and right cell
	 * and computes net updates (left and right going waves) according to the f-wave approach.
	 * It also returns the maximum wave speed.
	 */
#ifdef VECTORIZE
    #pragma omp declare simd
#endif
	void computeNetUpdates ( T i_hLeft,  T i_hRight,
                             T i_huLeft, T i_huRight,
                             T i_bLeft,  T i_bRight,

                             T &o_hUpdateLeft,
                             T &o_hUpdateRight,
                             T &o_huUpdateLeft,
                             T &o_huUpdateRight,
                             T &o_maxWaveSpeed ) const
	{
		  // determine the wet dry state and corr. values, if necessary.
		  if( i_hLeft >= dryTol ) {
		     if ( i_hRight < dryTol ) {
		      // Dry/Wet case
		      // Set values according to wall boundary condition
		      i_hRight = i_hLeft;
		      i_huRight = -i_huLeft;
		      i_bRight = i_bLeft;
		     } 
		  } else if ( i_hRight >= dryTol ) {
		      // Wet/Dry case
		      // Set values according to wall boundary condition
		      i_hLeft = i_hRight;
		      i_huLeft = -i_huRight;
		      i_bLeft = i_bRight;
		  } else {
		      // Dry/Dry case
		      // Set dummy values such that the result is zero
		      i_hLeft = dryTol;
		      i_huLeft = 0.; i_bLeft = 0.;
		      i_hRight = dryTol;
		      i_huRight = 0.; i_bRight = 0.;
		  };

		  //! velocity on the left side of the edge
		  T uLeft = i_huLeft / i_hLeft;                               // 1 FLOP (div)
		  //! velocity on the right side of the edge
		  T uRight = i_huRight / i_hRight;                            // 1 FLOP (div)

		  //! wave speeds of the f-waves
		  T waveSpeeds0 = 0., waveSpeeds1 = 0.;

		  //compute the wave speeds
		  fWaveComputeWaveSpeeds( i_hLeft, i_hRight,
		                          i_huLeft, i_huRight,
		                          uLeft, uRight,
		                          i_bLeft, i_bRight,

		                          waveSpeeds0, waveSpeeds1 );         // 20 FLOPs (incl. 3 sqrt, 1 div, 2 min/max)

		  //! variables to store the two f-waves
		  T fWaves0 = 0., fWaves1 = 0.;

		  //compute the decomposition into f-waves
		  fWaveComputeWaveDecomposition( i_hLeft, i_hRight,
		                                 i_huLeft, i_huRight,
		                                 uLeft, uRight,
		                                 i_bLeft, i_bRight,

		                                 waveSpeeds0, waveSpeeds1,
		                                 fWaves0, fWaves1);           // 23 FLOPs (incl. 1 div)

		  //compute the net-updates
		  o_hUpdateLeft = 0.;
		  o_hUpdateRight = 0.;
		  o_huUpdateLeft = 0.;
		  o_huUpdateRight = 0.;

		  //1st wave family
		  if(waveSpeeds0 < -zeroTol) { //left going
		    o_hUpdateLeft +=  fWaves0;
		    o_huUpdateLeft += fWaves0 * waveSpeeds0;                    // 3 FLOPs (assume left going wave ...)
		  }
		  else if(waveSpeeds0 > zeroTol) { //right going
		    o_hUpdateRight +=  fWaves0;
		    o_huUpdateRight += fWaves0 * waveSpeeds0;
		  }
		  else { //split waves, if waveSpeeds0 close to 0
		    o_hUpdateLeft +=   (T).5*fWaves0;
		    o_huUpdateLeft +=  (T).5*fWaves0 * waveSpeeds0;
		    o_hUpdateRight +=  (T).5*fWaves0;
		    o_huUpdateRight += (T).5*fWaves0 * waveSpeeds0;
		  }

		  //2nd wave family
		  if(waveSpeeds1 > zeroTol) { //right going
		    o_hUpdateRight +=  fWaves1;
		    o_huUpdateRight += fWaves1 * waveSpeeds1;                   // 3 FLOPs (assume right going wave ...)
		  }
		  else if(waveSpeeds1 < -zeroTol) { //left going
			o_hUpdateLeft +=  fWaves1;
			o_huUpdateLeft += fWaves1 * waveSpeeds1;
		  }
		  else { //split waves
		    o_hUpdateLeft +=   (T).5*fWaves1;
		    o_huUpdateLeft +=  (T).5*fWaves1 * waveSpeeds1;
		    o_hUpdateRight +=  (T).5*fWaves1;
		    o_huUpdateRight += (T).5*fWaves1 * waveSpeeds1;
		  }

		  //compute maximum wave speed (-> CFL-condition)
		  o_maxWaveSpeed = std::max( std::abs(waveSpeeds0) , std::abs(waveSpeeds1) ); 
		                                                              // 3 FLOPs (2 abs, 1 max)
		                                                              //========================
		                                                              // 54 FLOPs (3 sqrt, 4 div, 2 abs, 3 min/max)
    }

#ifdef VECTORIZE
    #pragma omp declare simd
#endif
	inline
	void fWaveComputeWaveSpeeds(
            const T i_hLeft,  const T i_hRight,
            const T i_huLeft, const T i_huRight,
            const T i_uLeft,  const T i_uRight,
            const T i_bLeft,  const T i_bRight,

            T &o_waveSpeed0, T &o_waveSpeed1 ) const
	{
		// helper variables for sqrt of h:
		T sqrt_hLeft = std::sqrt(i_hLeft);                                // 1 FLOP (sqrt)
		T sqrt_hRight = std::sqrt(i_hRight);                              // 1 FLOP (sqrt)

		// compute eigenvalues of the jacobian matrices 
		// in states Q_{i-1} and Q_{i}
		T characteristicSpeed0 = i_uLeft - sqrt_gravity * sqrt_hLeft;     // 2 FLOPs
		T characteristicSpeed1 = i_uRight + sqrt_gravity * sqrt_hRight;   // 2 FLOPs

		// compute "Roe averages"
		T hRoe = (T).5 * (i_hRight + i_hLeft);                            // 2 FLOPs
		T sqrt_hRoe = std::sqrt(hRoe);                                    // 1 FLOP (sqrt)
		T uRoe = i_uLeft * sqrt_hLeft + i_uRight * sqrt_hRight;           // 3 FLOPs
		uRoe /= sqrt_hLeft + sqrt_hRight;                                 // 2 FLOPs (1 div)

		// compute "Roe speeds" from Roe averages
		T roeSpeed0 = uRoe - sqrt_gravity * sqrt_hRoe;                    // 2 FLOPs
		T roeSpeed1 = uRoe + sqrt_gravity * sqrt_hRoe;                    // 2 FLOPs

		// compute Eindfeldt speeds (returned as output parameters)
		o_waveSpeed0 = std::min(characteristicSpeed0, roeSpeed0);         // 1 FLOP (min)
		o_waveSpeed1 = std::max(characteristicSpeed1, roeSpeed1);         // 1 FLOP (max)
		                                                                  //==============
		                                                                  //20 FLOPs (incl. 3 sqrt, 1 div, 2 min/max)
	}

#ifdef VECTORIZE
    #pragma omp declare simd
#endif
	inline
	void fWaveComputeWaveDecomposition(
			const T i_hLeft,  const T i_hRight,
			const T i_huLeft, const T i_huRight,
			const T i_uLeft,  const T i_uRight,
			const T i_bLeft,  const T i_bRight,
			const T i_waveSpeed0, const T i_waveSpeed1,

			T &o_fWave0, T &o_fWave1 ) const
	{
	  //calculate modified (bathymetry!) flux difference
	  // f(Q_i) - f(Q_{i-1}) -> serve as right hand sides
	  T fDif0 = i_huRight - i_huLeft;                                        // 1 FLOP
	  T fDif1 = i_huRight * i_uRight + half_gravity * i_hRight * i_hRight
	          -(i_huLeft  * i_uLeft  + half_gravity * i_hLeft  * i_hLeft);   // 9 FLOPs

	  // \delta x \Psi[2]
	  fDif1 += half_gravity * (i_hRight + i_hLeft)*(i_bRight - i_bLeft);     // 5 FLOPs

	  // solve linear system of equations to obtain f-waves:
	  // (       1            1      ) ( o_fWave0 ) = ( fDif0 )
	  // ( i_waveSpeed0 i_waveSpeed1 ) ( o_fWave1 )   ( fDif1 )
	  
	  // compute the inverse of the wave speed difference:
	  T inverseSpeedDiff = (T)1. / ( i_waveSpeed1 - i_waveSpeed0 );          // 2 FLOPs (1 div)
	  // compute f-waves:
	  o_fWave0 = (  i_waveSpeed1 * fDif0 - fDif1 ) * inverseSpeedDiff;       // 3 FLOPs
	  o_fWave1 = ( -i_waveSpeed0 * fDif0 + fDif1 ) * inverseSpeedDiff;       // 3 FLOPs
	                                                                         //=========
	                                                                         //23 FLOPs in total (incl. 1 div)
	}
};

}

#endif // FWAVEVEC_HPP_
