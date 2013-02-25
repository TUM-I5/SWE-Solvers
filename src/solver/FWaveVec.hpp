/**
 * FWaveCuda.h
 *
 ****
 **** This is a C++ wrapper for the Cuda implementation of the F-Wave solver (FWave.hpp).
 ****
 *
 * Created on: Nov 13, 2012
 * Last Update: Nov 16, 2012
 *
 ****
 *
 *  Author: Sebastian Rettenberger
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.
 *    E-Mail: rettenbs AT in.tum.de
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
 *            year={2002},
 *            publisher={Citeseer}}
 *
 *   @book{leveque2002finite,
 *         Author = {LeVeque, R. J.},
 *         Date-Added = {2011-09-13 14:09:31 +0000},
 *         Date-Modified = {2011-10-31 09:46:40 +0000},
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
	const T gravity;
	const T zeroTol;

public:
	FWaveVec(T i_dryTol = (T) 100,
			 T i_gravity = (T) 9.81,
			 T i_zeroTol = (T) 0.0000001)
		: dryTol(i_dryTol),
		  gravity(i_gravity),
		  zeroTol(i_zeroTol)
	{
	}

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
		  if( i_hLeft < dryTol && i_hRight < dryTol ) {
		      // Dry/Dry case
		      // Set dummy values such that the result is zero
		      i_hLeft = dryTol;
		      i_huLeft = 0.; i_bLeft = 0.;
		      i_hRight = dryTol;
		      i_huRight = 0.; i_bRight = 0.;
		  } else if ( i_hLeft < dryTol ) {
		      i_hLeft = i_hRight;
		      i_huLeft = -i_huRight;
		      i_bLeft = i_bRight;
		  } else if ( i_hRight < dryTol ) {
		      i_hRight = i_hLeft;
		      i_huRight = -i_huLeft;
		      i_bRight = i_bLeft;
		  }

		  //! velocity on the left side of the edge
		  T uLeft = i_huLeft / i_hLeft;
		  //! velocity on the right side of the edge
		  T uRight = i_huRight / i_hRight;

		  //! wave speeds of the f-waves
		  T waveSpeeds0 = 0., waveSpeeds1 = 0.;

		  //compute the wave speeds
		  fWaveComputeWaveSpeeds( i_hLeft, i_hRight,
		                          i_huLeft, i_huRight,
		                          uLeft, uRight,
		                          i_bLeft, i_bRight,

		                          waveSpeeds0, waveSpeeds1 );

		  //! where to store the two f-waves
		  T fWaves0 = 0., fWaves1 = 0.;

		  //compute the decomposition into f-waves
		  fWaveComputeWaveDecomposition( i_hLeft, i_hRight,
		                                 i_huLeft, i_huRight,
		                                 uLeft, uRight,
		                                 i_bLeft, i_bRight,

		                                 waveSpeeds0, waveSpeeds1,
		                                 fWaves0, fWaves1);

		  //compute the net-updates
		  T hUpdateLeft = 0.;
		  T hUpdateRight = 0.;
		  T huUpdateLeft = 0.;
		  T huUpdateRight = 0.;

		  //1st wave family
		  if(waveSpeeds0 < -zeroTol) { //left going
		    hUpdateLeft +=  fWaves0;
		    huUpdateLeft += fWaves0 * waveSpeeds0;
		  }
		  else if(waveSpeeds0 > zeroTol) { //right going
		    hUpdateRight +=  fWaves0;
		    huUpdateRight += fWaves0 * waveSpeeds0;
		  }
		  else { //split waves
		    hUpdateLeft +=   (T).5*fWaves0;
		    huUpdateLeft +=  (T).5*fWaves0 * waveSpeeds0;
		    hUpdateRight +=  (T).5*fWaves0;
		    huUpdateRight += (T).5*fWaves0 * waveSpeeds0;
		  }

		  //2nd wave family
		  if(waveSpeeds1 > zeroTol) { //right going
		    hUpdateRight +=  fWaves1;
		    huUpdateRight += fWaves1 * waveSpeeds1;
		  }
		  else if(waveSpeeds1 < -zeroTol) { //left going
			hUpdateLeft +=  fWaves1;
			huUpdateLeft += fWaves1 * waveSpeeds1;
		  }
		  else { //split waves
		    hUpdateLeft +=   (T).5*fWaves1;
		    huUpdateLeft +=  (T).5*fWaves1 * waveSpeeds1;
		    hUpdateRight +=  (T).5*fWaves1;
		    huUpdateRight += (T).5*fWaves1 * waveSpeeds1;
		  }

		  // Set output variables
		  o_hUpdateLeft = hUpdateLeft;
		  o_hUpdateRight = hUpdateRight;
		  o_huUpdateLeft = huUpdateLeft;
		  o_huUpdateRight = huUpdateRight;

		  //compute maximum wave speed (-> CFL-condition)
		  o_maxWaveSpeed = std::max( std::abs(waveSpeeds0) , std::abs(waveSpeeds1) );
    }

	inline
	void fWaveComputeWaveSpeeds(
            const T i_hLeft,  const T i_hRight,
            const T i_huLeft, const T i_huRight,
            const T i_uLeft,  const T i_uRight,
            const T i_bLeft,  const T i_bRight,

            T &o_waveSpeed0, T &o_waveSpeed1 ) const
	{
		//compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i}
		T characteristicSpeed0 = 0., characteristicSpeed1 = 0.;
		characteristicSpeed0 = i_uLeft - std::sqrt(gravity*i_hLeft);
		characteristicSpeed1 = i_uRight + std::sqrt(gravity*i_hRight);

		//compute "Roe speeds"
		T hRoe = (T).5 * (i_hRight + i_hLeft);
		T uRoe = i_uLeft * std::sqrt(i_hLeft) + i_uRight * std::sqrt(i_hRight);
		uRoe /= std::sqrt(i_hLeft) + std::sqrt(i_hRight);

		T roeSpeed0 = 0., roeSpeed1 = 0.;
		roeSpeed0 = uRoe - std::sqrt(gravity*hRoe);
		roeSpeed1 = uRoe + std::sqrt(gravity*hRoe);

		//computer eindfeldt speeds
		o_waveSpeed0 = std::min(characteristicSpeed0, roeSpeed0);
		o_waveSpeed1 = std::max(characteristicSpeed1, roeSpeed1);
	}

	inline
	void fWaveComputeWaveDecomposition(
			const T i_hLeft,  const T i_hRight,
			const T i_huLeft, const T i_huRight,
			const T i_uLeft,  const T i_uRight,
			const T i_bLeft,  const T i_bRight,
			const T i_waveSpeed0, const T i_waveSpeed1,

			T &o_fWave0, T &o_fWave1 ) const
	{
	  T lambdaDif = i_waveSpeed1 - i_waveSpeed0;

	  //compute the inverse matrix R^{-1}
	  T Rinv00 = 0., Rinv01 = 0., Rinv10 = 0., Rinv11 = 0.;

	  T oneDivLambdaDif = (T)1. /  lambdaDif;
	  Rinv00 = oneDivLambdaDif *  i_waveSpeed1;
	  Rinv01 = -oneDivLambdaDif;

	  Rinv10 = oneDivLambdaDif * -i_waveSpeed0;
	  Rinv11 = oneDivLambdaDif;

	  //right hand side
	  T fDif0 = 0., fDif1 = 0.;

	  //calculate modified (bathymetry!) flux difference
	  // f(Q_i) - f(Q_{i-1})
	  fDif0 = i_huRight - i_huLeft;
	  fDif1 = i_huRight * i_uRight + (T).5 * gravity * i_hRight * i_hRight
	          -(i_huLeft  * i_uLeft  + (T).5 * gravity * i_hLeft  * i_hLeft);

	  // \delta x \Psi[2]
	  T psi = -gravity * (T).5 * (i_hRight + i_hLeft)*(i_bRight - i_bLeft);
	  fDif1 -= psi;

	  //solve linear equations
	  o_fWave0 = Rinv00 * fDif0 + Rinv01 * fDif1;
	  o_fWave1 = Rinv10 * fDif0 + Rinv11 * fDif1;
	}
};

}

#endif // FWAVEVEC_HPP_
