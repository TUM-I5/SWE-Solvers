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
		  T waveSpeeds[2];

		  //compute the wave speeds
		  fWaveComputeWaveSpeeds( i_hLeft, i_hRight,
		                          i_huLeft, i_huRight,
		                          uLeft, uRight,
		                          i_bLeft, i_bRight,

		                          waveSpeeds[0], waveSpeeds[1] );

		  //! where to store the two f-waves
		  T fWaves[2];

		  //compute the decomposition into f-waves
		  fWaveComputeWaveDecomposition( i_hLeft, i_hRight,
		                                 i_huLeft, i_huRight,
		                                 uLeft, uRight,
		                                 i_bLeft, i_bRight,

		                                 waveSpeeds[0], waveSpeeds[1],
		                                 fWaves[0], fWaves[1]);

		  //compute the net-updates
		  T hUpdateLeft = 0.;
		  T hUpdateRight = 0.;
		  T huUpdateLeft = 0.;
		  T huUpdateRight = 0.;

		  //1st wave family
		  if(waveSpeeds[0] < -zeroTol) { //left going
		    hUpdateLeft +=  fWaves[0];
		    huUpdateLeft += fWaves[0] * waveSpeeds[0];
		  }
		  else if(waveSpeeds[0] > zeroTol) { //right going
		    hUpdateRight +=  fWaves[0];
		    huUpdateRight += fWaves[0] * waveSpeeds[0];
		  }
		  else { //split waves
		    hUpdateLeft +=   (T).5*fWaves[0];
		    huUpdateLeft +=  (T).5*fWaves[0] * waveSpeeds[0];
		    hUpdateRight +=  (T).5*fWaves[0];
		    huUpdateRight += (T).5*fWaves[0] * waveSpeeds[0];
		  }

		  //2nd wave family
		  if(waveSpeeds[1] < -zeroTol) { //left going
		    hUpdateLeft +=  fWaves[1];
		    huUpdateLeft += fWaves[1] * waveSpeeds[1];
		  }
		  else if(waveSpeeds[1] > zeroTol) {
		    hUpdateRight +=  fWaves[1];
		    huUpdateRight += fWaves[1] * waveSpeeds[1];
		  }
		  else { //split waves
		    hUpdateLeft +=   (T).5*fWaves[1];
		    huUpdateLeft +=  (T).5*fWaves[1] * waveSpeeds[1];
		    hUpdateRight +=  (T).5*fWaves[1];
		    huUpdateRight += (T).5*fWaves[1] * waveSpeeds[1];
		  }

		  // Set output variables
		  o_hUpdateLeft = hUpdateLeft;
		  o_hUpdateRight = hUpdateRight;
		  o_huUpdateLeft = huUpdateLeft;
		  o_huUpdateRight = huUpdateRight;

		  //compute maximum wave speed (-> CFL-condition)
		  o_maxWaveSpeed = std::max( std::abs(waveSpeeds[0]) , std::abs(waveSpeeds[1]) );
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
		T characteristicSpeed[2];
		characteristicSpeed[0] = i_uLeft - std::sqrt(gravity*i_hLeft);
		characteristicSpeed[1] = i_uRight + std::sqrt(gravity*i_hRight);

		//compute "Roe speeds"
		T hRoe = (T).5 * (i_hRight + i_hLeft);
		T uRoe = i_uLeft * std::sqrt(i_hLeft) + i_uRight * std::sqrt(i_hRight);
		uRoe /= std::sqrt(i_hLeft) + std::sqrt(i_hRight);

		T roeSpeed[2];
		roeSpeed[0] = uRoe - std::sqrt(gravity*hRoe);
		roeSpeed[1] = uRoe + std::sqrt(gravity*hRoe);

		//computer eindfeldt speeds
		o_waveSpeed0 = std::min(characteristicSpeed[0], roeSpeed[0]);
		o_waveSpeed1 = std::max(characteristicSpeed[1], roeSpeed[1]);
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
	  T Rinv[2][2];

	  T oneDivLambdaDif = (T)1. /  lambdaDif;
	  Rinv[0][0] = oneDivLambdaDif *  i_waveSpeed1;
	  Rinv[0][1] = -oneDivLambdaDif;

	  Rinv[1][0] = oneDivLambdaDif * -i_waveSpeed0;
	  Rinv[1][1] = oneDivLambdaDif;

	  //right hand side
	  T fDif[2];

	  //calculate modified (bathymetry!) flux difference
	  // f(Q_i) - f(Q_{i-1})
	  fDif[0] = i_huRight - i_huLeft;
	  fDif[1] = i_huRight * i_uRight + (T).5 * gravity * i_hRight * i_hRight
	          -(i_huLeft  * i_uLeft  + (T).5 * gravity * i_hLeft  * i_hLeft);

	  // \delta x \Psi[2]
	  T psi = -gravity * (T).5 * (i_hRight + i_hLeft)*(i_bRight - i_bLeft);
	  fDif[1] -= psi;

	  //solve linear equations
	  o_fWave0 = Rinv[0][0] * fDif[0] + Rinv[0][1] * fDif[1];
	  o_fWave1 = Rinv[1][0] * fDif[0] + Rinv[1][1] * fDif[1];
	}
};

}

#endif // FWAVEVEC_HPP_
