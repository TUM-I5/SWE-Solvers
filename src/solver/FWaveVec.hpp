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
		  // (voodoo magic to work around if-else statements)
		  bool hDryTol = i_hLeft >= dryTol;
		  i_hLeft = hDryTol*i_hLeft + (!hDryTol)*i_hRight;
		  i_huLeft = hDryTol*i_huLeft + (!hDryTol)*(-i_huRight);
		  i_bLeft = hDryTol*i_bLeft + (!hDryTol)*i_bRight;

		  hDryTol = i_hRight >= dryTol;
		  i_hRight = hDryTol*i_hRight + (!hDryTol)*i_hLeft;
		  i_huRight = hDryTol*i_huRight + (!hDryTol)*(-i_huLeft);
		  i_bRight = hDryTol*i_bRight + (!hDryTol)*i_bLeft;

		  // If both cell are dry, make sure we get no waves
		  hDryTol = i_hLeft >= dryTol || i_hRight >= dryTol;
		  i_hLeft = hDryTol*i_hLeft + (!hDryTol)*dryTol;
		  i_hRight = hDryTol*i_hRight + (!hDryTol)*dryTol;

		  //! velocity on the left side of the edge
		  T uLeft = i_huLeft / i_hLeft;
		  //! velocity on the right side of the edge
		  T uRight = i_huRight / i_hRight;

		  //! wave speeds of the f-waves
		  T waveSpeed[2];

		  //compute the wave speeds
		  fWaveComputeWaveSpeeds( i_hLeft, i_hRight,
		                          i_huLeft, i_huRight,
		                          uLeft, uRight,
		                          i_bLeft, i_bRight,

		                          waveSpeed[0], waveSpeed[1] );

		  //! where to store the two f-waves
		  T fWave[2];

		  //compute the decomposition into f-waves
		  fWaveComputeWaveDecomposition( i_hLeft, i_hRight,
		                                 i_huLeft, i_huRight,
		                                 uLeft, uRight,
		                                 i_bLeft, i_bRight,

		                                 waveSpeed[0], waveSpeed[1],
		                                 fWave[0], fWave[1]);

		  //compute the net-updates (some more voodoo magic)
		  //1st wave family
		  bool waveSpeedZeroTol = std::abs(waveSpeed[0]) > zeroTol;
		  bool waveSpeedPositive = waveSpeed[0] > 0;

		  o_hUpdateLeft =   ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*(!waveSpeedPositive))
				* fWave[0];
		  o_huUpdateLeft =  ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*(!waveSpeedPositive))
				* fWave[0] * waveSpeed[0];
		  o_hUpdateRight =  ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*waveSpeedPositive)
				* fWave[0];
		  o_huUpdateRight = ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*waveSpeedPositive)
				* fWave[0] * waveSpeed[0];

		  //2nd wave family
		  waveSpeedZeroTol = std::abs(waveSpeed[1]) > zeroTol;
		  waveSpeedPositive = waveSpeed[1] > 0;

		  o_hUpdateLeft +=   ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*(!waveSpeedPositive))
				 * fWave[1];
		  o_huUpdateLeft +=  ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*(!waveSpeedPositive))
				 * fWave[1] * waveSpeed[1];
		  o_hUpdateRight +=  ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*waveSpeedPositive)
				 * fWave[1];
		  o_huUpdateRight += ((!waveSpeedZeroTol)*(T).5 + waveSpeedZeroTol*waveSpeedPositive)
				 * fWave[1] * waveSpeed[1];

		  //compute maximum wave speed (-> CFL-condition)
		  o_maxWaveSpeed = std::max( std::abs(waveSpeed[0]) , std::abs(waveSpeed[1]) );
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
