/**
 * AugRie.hpp
 *
 ****
 **** Approximate Augmented Riemann Solver for the Shallow Water Equations
 ****
 *
 *  Created on: Sep 12, 2011
 *  Last Update: Feb 18, 2012
 *
 ****
 *
 *  Author: Alexander Breuer
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer
 *    E-Mail: breuera AT in.tum.de
 *
 *  Some optimizations: Martin Schreiber
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Martin_Schreiber
 *    E-Mail: schreibm AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @phdthesis{george2006finite,
 *              Author = {George, D.L.},
 *              Title = {Finite volume methods and adaptive refinement for tsunami propagation and inundation},
 *              Year = {2006}}
 *
 *   @article{george2008augmented,
 *            Author = {George, D.L.},
 *            Journal = {Journal of Computational Physics},
 *            Number = {6},
 *            Pages = {3089--3113},
 *            Publisher = {Elsevier},
 *            Title = {Augmented Riemann solvers for the shallow water equations over variable topography with steady states and inundation},
 *            Volume = {227},
 *            Year = {2008}}
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
 *
 * Acknowledgments:
 *   Special thanks go to R.J. LeVeque and D.L. George for publishing their code
 *   and the corresponding documentation (-> Literature).
 */


/*
 * TODO: store nLow/nHigh variables in array[2]
 */


#ifndef AUGRIE_HPP_
#define AUGRIE_HPP_

#include "WavePropagation.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

/** switch features of the solver on/off
 *
 *  The solver is not strictly positivity preserving with correctrarefactions.
 *
 *  Solver seems to fail in the "Single wave on a simple beach"-benchmark
 *  if correctrarefactions or complexsteadystatewave is used.
 *
 *  TODO: Further investigation is recommended.
 */
//#define correctrarefactions
//#define complexsteadystatewave

namespace solver {
  template <typename T>	class AugRie;
}

/**
 * Approximate Augmented Riemann Solver for the Shallow Water Equations.
 *
 * T should be double or float.
 */
template <typename T> class solver::AugRie: public WavePropagation<T> {
  private: //explicit for unit tests
    //use nondependent names (template base class)
    using solver::WavePropagation<T>::zeroTol;
    using solver::WavePropagation<T>::g;
    using solver::WavePropagation<T>::dryTol;

    using solver::WavePropagation<T>::hLeft;
    using solver::WavePropagation<T>::hRight;
    using solver::WavePropagation<T>::huLeft;
    using solver::WavePropagation<T>::huRight;
    using solver::WavePropagation<T>::bLeft;
    using solver::WavePropagation<T>::bRight;
    using solver::WavePropagation<T>::uLeft;
    using solver::WavePropagation<T>::uRight;

    using solver::WavePropagation<T>::wetDryState;
    using solver::WavePropagation<T>::DryDry;
    using solver::WavePropagation<T>::WetWet;
    using solver::WavePropagation<T>::WetDryInundation;
    using solver::WavePropagation<T>::WetDryWall;
    using solver::WavePropagation<T>::WetDryWallInundation;
    using solver::WavePropagation<T>::DryWetInundation;
    using solver::WavePropagation<T>::DryWetWall;
    using solver::WavePropagation<T>::DryWetWallInundation;

    using solver::WavePropagation<T>::storeParameters;


    //! Newton-tolerance (exit Newton-Raphson-method, if we are close enough to the root)
    const T newtonTol;
    //! maximum number of Newton-Raphson-Iterations
    const int maxNumberOfNewtonIterations;

    //! height of our homogeneous Riemann-problem at middle state (computed by determineMiddleState)
    T hMiddle;

    //! shock or inner rarefaction speeds of our homogeneous Riemann-problem (computed by computeMiddleState)
    T middleStateSpeeds[2];


    /**
     * the Riemann-struture of the homogeneous Riemann-problem.
     */
    enum RiemannStructure {
      DrySingleRarefaction,  /**< 1st wave family: contact discontinuity; 2nd wave family: rarefaction. */
      SingleRarefactionDry,  /**< 1st wave family: rarefaction; 2nd wave family: contact discontinuity. */
      ShockShock,            /**< 1st wave family: shock; 2nd wave family: shock. */
      ShockRarefaction,      /**< 1st wave family: shock; 2nd wave family: rarefaction. */
      RarefactionShock,      /**< 1st wave family: rarefaction; 2nd wave family: shock. */
      RarefactionRarefaction /**< 1st wave family: rarefaction; 2nd wave family: rarefaction. */
    };

    //! Riemann-structure of our homogeneous Riemann-problem (determined by determineRiemannStructure)
    RiemannStructure riemannStructure;

  public:
    /**
     * Constructor of the Augmented Riemann solver with optional parameters.
     * 
     * @param i_dryTolerance numerical definition of "dry".
     * @param i_gravity gravity constant.
     * @param i_newtonTolerance numerical definition of "convergence" (used in the AugRie solver only).
     * @param i_maxNumberOfNewtonIterations maximum steps for the Newton-Raphson method (used in the AugRie solver only).
     * @param i_zeroTolerance numerical definition of zero.
     */
    AugRie(	T i_dryTolerance =                  (T) 0.01,
            T i_gravity =                       (T) 9.81,
            T i_newtonTolerance =               (T) 0.000001,
            int i_maxNumberOfNewtonIterations =     10,
            T i_zeroTolerance =                   (T) 0.00001 ):
              WavePropagation<T>( i_dryTolerance, i_gravity, i_zeroTolerance ),
              newtonTol(i_newtonTolerance),
              maxNumberOfNewtonIterations(i_maxNumberOfNewtonIterations) {

#ifndef NDEBUG
#ifndef SUPPRESS_SOLVER_DEBUG_OUTPUT
      //print some information about the used solver
      std::cout << "  *** solver::AugRie created" << std::endl
                << "    zeroTolerance=" << zeroTol << std::endl
                << "    gravity=" << g << std::endl
                << "    dryTolerance=" << dryTol << std::endl
                << "    newtonTolerance=" << newtonTol << std::endl
                << "    maxNumberNumberOfNewtonIterations=" << maxNumberOfNewtonIterations << std::endl
                << "\n  optional features:"
                #ifdef correctrarefactions
                  << "   correctrarefactions"
                #endif
                #ifdef complexsteadystatewave
                  << "   complexsteadystatewave"
                #endif
                << "\n  ***\n\n";
#endif
#endif

    }
    //declare variables which are used over and over again
#if 0
    T sqrt_g_h[2];
    T sqrt_h[2];

#define sqrt_g_hLeft  (sqrt_g_h[0])
#define sqrt_g_hRight (sqrt_g_h[1])

#define sqrt_hLeft    (sqrt_h[0])
#define sqrt_hRight   (sqrt_h[1])

#else
    T sqrt_g_hLeft;
    T sqrt_g_hRight;

    T sqrt_hLeft;
    T sqrt_hRight;
#endif


    T sqrt_g;


    /**
     * Compute net updates for the left/right cell of the edge.
     *
     * maxWaveSpeed will be set to the maximum (linearized) wave speed => CFL
     */
    void computeNetUpdates ( const T &i_hLeft,  const T &i_hRight,
                             const T &i_huLeft, const T &i_huRight,
                             const T &i_bLeft,  const T &i_bRight,

                             T &o_hUpdateLeft,
                             T &o_hUpdateRight,
                             T &o_huUpdateLeft,
                             T &o_huUpdateRight,
                             T &o_maxWaveSpeed
#if CONFIG_TSUNAMI_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
                            ,T  o_eigenCoefficients[3]
#endif
                            ) {
      // store parameters to member variables
      storeParameters( i_hLeft, i_hRight,
                       i_huLeft, i_huRight,
                       i_bLeft, i_bRight );

      //set speeds to zero (will be determined later)
      uLeft = uRight = 0.;

      //reset net updates and the maximum wave speed
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight  = (T)0.;
      o_maxWaveSpeed = (T)0.;
      
#if CONFIG_TSUNAMI_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
      // reset eigen coefficients
      o_eigenCoefficients[0] = o_eigenCoefficients[1] = o_eigenCoefficients[2] = 0;
#endif

      //determine the wet/dry state and compute local variables correspondingly
      determineWetDryState();

      if (wetDryState == DryDry) //nothing to do in a dry region
        return;

      //precompute some terms which are fixed during
      //the computation after some specific point
      sqrt_g = std::sqrt(g);
      sqrt_hLeft = std::sqrt(hLeft);
      sqrt_hRight = std::sqrt(hRight);

      sqrt_g_hLeft = sqrt_g*sqrt_hLeft;
      sqrt_g_hRight = sqrt_g*sqrt_hRight;


      //where to store the three waves
      T fWaves[3][2];
      //and their speeds
      T waveSpeeds[3];

      //compute the augmented decomposition
      //  (thats the place where the computational work is done..)
      computeWaveDecomposition( fWaves,
                                waveSpeeds
#if CONFIG_TSUNAMI_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
                               ,o_eigenCoefficients
#endif
                               );


      //compute the updates from the three propagating waves
      //A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      //A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
        if (waveSpeeds[waveNumber] < -zeroTol) { //left going
          o_hUpdateLeft +=  fWaves[waveNumber][0];
          o_huUpdateLeft += fWaves[waveNumber][1];
        }

        else if (waveSpeeds[waveNumber] > zeroTol) { //right going
          o_hUpdateRight +=  fWaves[waveNumber][0];
          o_huUpdateRight += fWaves[waveNumber][1];
        }
        else { //TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
          o_hUpdateLeft +=  (T).5*fWaves[waveNumber][0];
          o_huUpdateLeft += (T).5*fWaves[waveNumber][1];

          o_hUpdateRight +=  (T).5*fWaves[waveNumber][0];
          o_huUpdateRight += (T).5*fWaves[waveNumber][1];
        }
        
        //no wave speeds => zero strength fWaves
//        assert( std::fabs(fWaves[waveNumber][0]) < zeroTol );
//        assert( std::fabs(fWaves[waveNumber][1]) < zeroTol );
      }

#ifndef NDEBUG
      if(wetDryState == DryWetWall)
        assert( std::fabs(o_hUpdateLeft) < zeroTol && std::fabs(o_huUpdateLeft) < zeroTol );
      else if(wetDryState == WetDryWall)
        assert( std::fabs(o_hUpdateRight) < zeroTol && std::fabs(o_huUpdateRight) < zeroTol );
#endif

      //compute maximum wave speed (-> CFL-condition)
      waveSpeeds[0] = std::fabs(waveSpeeds[0]);
      waveSpeeds[1] = std::fabs(waveSpeeds[1]);
      waveSpeeds[2] = std::fabs(waveSpeeds[2]);

      o_maxWaveSpeed = std::max(waveSpeeds[0], waveSpeeds[1]);
      o_maxWaveSpeed = std::max(o_maxWaveSpeed, waveSpeeds[2]);
    }

  private:
    /**
     * Determine the wet/dry state and set member variables accordingly.
     */
    void determineWetDryState() {
      //compute speeds or set them to zero (dry cells)
      if (hLeft > dryTol) {
        uLeft = huLeft / hLeft;
      }
      else {
        bLeft += hLeft;
        hLeft = huLeft = uLeft = 0;
      }

      if (hRight > dryTol) {
        uRight = huRight / hRight;
      }
      else {
        bRight += hRight;
        hRight = huRight = uRight = 0;
      }

      //test for simple wet/wet case since this is most probably the
      //most frequently executed branch
      if (hLeft >= dryTol && hRight >= dryTol) {
        wetDryState =  WetWet;
      }

      //check for the dry/dry-case
      else if (hLeft < dryTol && hRight < dryTol)
        wetDryState =  DryDry;

      //we have a shoreline: one cell dry, one cell wet

      //check for simple inundation problems
      // (=> dry cell lies lower than the wet cell)
      else if (hLeft < dryTol && hRight + bRight > bLeft)
        wetDryState =  DryWetInundation;

      else if (hRight < dryTol && hLeft + bLeft > bRight)
        wetDryState =  WetDryInundation;

      //dry cell lies higher than the wet cell
      else if (hLeft < dryTol) {
        //lets check if the momentum is able to overcome the difference in height
        //  => solve homogeneous Riemann-problem to determine the middle state height
        //     which would arise if there is a wall (wall-boundary-condition)
        //       \cite[ch. 6.8.2]{george2006finite})
        //       \cite[ch. 5.2]{george2008augmented}
        computeMiddleState( hRight,  hRight,
                           -uRight, uRight,
                           -huRight, huRight,
                           maxNumberOfNewtonIterations );

        if (hMiddle + bRight > bLeft) {
          //momentum is large enough, continue with the original values
//          bLeft = hMiddle + bRight;
          wetDryState =  DryWetWallInundation;
        }
        else {
          //momentum is not large enough, use wall-boundary-values
          hLeft = hRight;
          uLeft = -uRight;
          huLeft = -huRight;
          bLeft = bRight = (T)0.;
          wetDryState =  DryWetWall;
        }
      }

      else if (hRight < dryTol) {
        //lets check if the momentum is able to overcome the difference in height
        //  => solve homogeneous Riemann-problem to determine the middle state height
        //     which would arise if there is a wall (wall-boundary-condition)
        //       \cite[ch. 6.8.2]{george2006finite})
        //       \cite[ch. 5.2]{george2008augmented}
        computeMiddleState( hLeft,  hLeft,
                            uLeft, -uLeft,
                            huLeft, -huLeft,
                            maxNumberOfNewtonIterations );

        if (hMiddle + bLeft > bRight) {
          //momentum is large enough, continue with the original values
//          bRight = hMiddle + bLeft;
          wetDryState =  WetDryWallInundation;
        }
        else {
          hRight = hLeft;
          uRight = -uLeft;
          huRight = -huLeft;
          bRight = bLeft = (T)0.;
          wetDryState =  WetDryWall;
        }
      }
      //done with all cases
      else assert(false);

      //limit the effect of the source term if there is a "wall"
      //\cite[end of ch. 5.2?]{george2008augmented}
      //\cite[rpn2ez_fast_geo.f][levequeclawpack]
      if ( wetDryState == DryWetWallInundation )
        bLeft = hRight + bRight;

      else if( wetDryState == WetDryWallInundation )
        bRight = hLeft + bLeft;
    }

    /**
     * Compute the augmented wave decomposition.
     *
     * @param o_fWaves will be set to: Decomposition into f-Waves.
     * @param o_waveSpeeds will be set to: speeds of the linearized waves (eigenvalues).
     */
    inline void computeWaveDecomposition( T o_fWaves[3][2],
                                          T o_waveSpeeds[3]
#if CONFIG_TSUNAMI_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
                                         ,T o_eigenCoefficients[3]
#endif
                                         )	{
      //compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
      T characteristicSpeeds[2];
      characteristicSpeeds[0] = uLeft - sqrt_g_hLeft;
      characteristicSpeeds[1] = uRight + sqrt_g_hRight;

      //compute "Roe speeds"
      T hRoe = (T).5 * (hRight + hLeft);
      T uRoe = uLeft * sqrt_hLeft + uRight * sqrt_hRight;
      uRoe /= sqrt_hLeft + sqrt_hRight;

      T roeSpeeds[2];
      //optimization for dumb compilers
      T sqrt_g_hRoe = std::sqrt(g*hRoe);
      roeSpeeds[0] = uRoe - sqrt_g_hRoe;
      roeSpeeds[1] = uRoe + sqrt_g_hRoe;

      //compute the middle state of the homogeneous Riemann-Problem
      if(wetDryState != WetDryWall && wetDryState != DryWetWall) { //case WDW and DWW was computed in determineWetDryState already
        computeMiddleState( hLeft, hRight,
                            uLeft, uRight,
                            huLeft, huRight
        );
      }

      //compute extended eindfeldt speeds (einfeldt speeds + middle state speeds)
      //  \cite[ch. 5.2]{george2008augmented}, \cite[ch. 6.8]{george2006finite}
      T extEinfeldtSpeeds[2]= {(T)0, (T)0};
      if( wetDryState == WetWet ||
          wetDryState == WetDryWall ||
          wetDryState == DryWetWall ) {
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[0] = std::min(extEinfeldtSpeeds[0], middleStateSpeeds[1]);

        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);
        extEinfeldtSpeeds[1] = std::max(extEinfeldtSpeeds[1], middleStateSpeeds[0]);
      }
      else if (hLeft < dryTol) { //ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(roeSpeeds[0], middleStateSpeeds[1]);
        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);

        assert(middleStateSpeeds[0] < extEinfeldtSpeeds[1]);
      }
      else if (hRight < dryTol) { //ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[1] = std::max(roeSpeeds[1], middleStateSpeeds[0]);

        assert(middleStateSpeeds[1] > extEinfeldtSpeeds[0]);
      }
      else {
        assert(false);
      }

      //HLL middle state
      //  \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
      T hLLMiddleHeight = huLeft - huRight + extEinfeldtSpeeds[1] * hRight - extEinfeldtSpeeds[0] * hLeft;
      hLLMiddleHeight /= extEinfeldtSpeeds[1] - extEinfeldtSpeeds[0];
      hLLMiddleHeight = std::max(hLLMiddleHeight, (T)0.);

      //define eigenvalues
      T eigenValues[3];
      eigenValues[0] = extEinfeldtSpeeds[0];
      //eigenValues[1] --> corrector wave
      eigenValues[2] = extEinfeldtSpeeds[1];

      //define eigenvectors
      T eigenVectors[3][3];

      //set first and third eigenvector
      eigenVectors[0][0] = (T)1;
      eigenVectors[0][2] = (T)1;

      eigenVectors[1][0] = eigenValues[0];
      eigenVectors[1][2] = eigenValues[2];

      eigenVectors[2][0] = eigenValues[0]*eigenValues[0];
      eigenVectors[2][2] = eigenValues[2]*eigenValues[2];


      //compute rarefaction corrector wave
      //  \cite[ch. 6.7.2]{george2006finite}, \cite[ch. 5.1]{george2008augmented}
      bool strongRarefaction = false;
#ifdef correctrarefactions
      if ( (riemannStructure == ShockRarefaction ||
           riemannStructure == RarefactionShock ||
           riemannStructure == RarefactionRarefaction) &&
           hMiddle > dryTol ) { //limit to ensure non-negative depth

        //TODO: GeoClaw, riemann_aug_JCP; No explicit boundaries for "strong rarefaction" in literature?
        T rareBarrier[2] = {.5, .9};

        //lower rare barrier in the case of a transsonic rarefaction
        if ( (riemannStructure == RarefactionShock || riemannStructure == RarefactionRarefaction) &&
             eigenValues[0]*middleStateSpeeds[0] < (T).0 ) { //transsonic rarefaction, first wave family
          rareBarrier[0] = (T).2;
        }
        else if ( (riemannStructure == ShockRarefaction || riemannStructure == RarefactionRarefaction) &&
                  eigenValues[2]*middleStateSpeeds[1] < (T)0. ) { //transsonic rarefaction, second wave family
          rareBarrier[0] = (T).2;
        }


        T sqrt_g_hMiddle = std::sqrt(g*hMiddle);
        //determine the max. rarefaction size (distance between head and tail)
        T rareFactionSize[2];
        rareFactionSize[0] = (T)3. * (sqrt_g_hLeft - sqrt_g_hMiddle);
        rareFactionSize[1] = (T)3. * (sqrt_g_hRight - sqrt_g_hMiddle);

        T maxRareFactionSize = std::max(rareFactionSize[0], rareFactionSize[1]);

        //set the eigenvalue of the corrector wave in the case of a "strong rarefaction"
        if ( maxRareFactionSize > rareBarrier[0] * (eigenValues[2] - eigenValues[0]) &&
             maxRareFactionSize < rareBarrier[1] * (eigenValues[2] - eigenValues[0])
        ) {
          strongRarefaction = true;
          if (rareFactionSize[0] > rareFactionSize[1])
            eigenValues[1] = middleStateSpeeds[0];
          else
            eigenValues[1] = middleStateSpeeds[1];
        }

        //TODO: implemented in clawpack, why?
//        if ( hMiddle < std::min(hLeft, hRight)/5.) { //middle state in an HLL solve
//          strongRarefaction = false;
//        }
      }
#endif

      //set 2nd eigenvector
      if(strongRarefaction == false) {
        eigenValues[1] = (T).5*(eigenValues[0]+eigenValues[2]);
        eigenVectors[0][1] = 0.;
        eigenVectors[1][1] = 0.;
        eigenVectors[2][1] = 1.;
      }
      else {
        eigenVectors[0][1] = (T)1.;
        eigenVectors[1][1] = eigenValues[1];
        eigenVectors[2][1] = eigenValues[1]*eigenValues[1];
      }

      //compute the jump in state
      T rightHandSide[3];
      rightHandSide[0] = hRight - hLeft;
      rightHandSide[1] = huRight - huLeft;
      rightHandSide[2] = huRight * uRight + (T).5 * g * hRight * hRight
                      -( huLeft  * uLeft  + (T).5 * g * hLeft  * hLeft );

      //compute steady state wave
      //  \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
      T steadyStateWave[2];
      T hBar = (hLeft + hRight)*(T).5;
#ifdef complexsteadystatewave
      T lLambdaBar = (uLeft + uRight)*(uLeft + uRight)*.25 - g*hBar;

      T lLambdaTilde = std::max((T)0., uLeft*uRight) - g*hBar;

      //near sonic as defined in geoclaw (TODO: literature?)
      if( ( std::fabs(lLambdaBar) < zeroTol )                                             ||
          ( lLambdaBar*lLambdaTilde < zeroTol )                                           ||
          ( lLambdaBar*eigenValues[0]*eigenValues[1] < zeroTol )                          ||
          ( std::min( std::fabs( eigenValues[0]), std::fabs(eigenValues[2]) ) < zeroTol ) ||
          ( eigenValues[0] < (T)0. && middleStateSpeeds[0] > (T)0. )                      ||
          ( eigenValues[2] > (T)0. && middleStateSpeeds[1] < (T)0. )                      ||
          ( (uLeft + sqrt(g*hLeft)) * (uRight+sqrt(g*hRight)) < (T)0. )                   ||
          ( (uLeft - sqrt(g*hLeft)) * (uRight-sqrt(g*hRight)) < (T)0. ) ) {
#endif
        steadyStateWave[0] = -(bRight - bLeft);
        steadyStateWave[1] = -g * hBar * (bRight - bLeft);
#ifdef complexsteadystatewave
      }
      else {
        steadyStateWave[0] = g * (hBar / lLambdaBar) * (bRight - bLeft);
        T hTilde = hBar * lLambdaTilde / lLambdaBar;

        //set bounds for problems far from steady state
        hTilde = std::max(std::min(hLeft, hRight), hTilde);
        hTilde = std::min(std::max(hLeft, hRight), hTilde);
        steadyStateWave[1] = -g*hTilde*(bRight - bLeft);
      }
#endif

      //preserve depth-positivity
      //  \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
      if(eigenValues[0] < -zeroTol && eigenValues[2] > zeroTol) { //subsonic
        steadyStateWave[0] = std::max(steadyStateWave[0], hLLMiddleHeight*(eigenValues[2]-eigenValues[0])/eigenValues[0]);
        steadyStateWave[0] = std::min(steadyStateWave[0], hLLMiddleHeight*(eigenValues[2]-eigenValues[0])/eigenValues[2]);
      }
      else if(eigenValues[0] > zeroTol) { //supersonic right TODO: motivation?
        steadyStateWave[0] = std::max(steadyStateWave[0], -hLeft);
        steadyStateWave[0] = std::min(steadyStateWave[0], hLLMiddleHeight*(eigenValues[2]-eigenValues[0])/eigenValues[0]);
      }
      else if(eigenValues[2] < -zeroTol) { //supersonic left TODO: motivation?
        steadyStateWave[0] = std::max(steadyStateWave[0], hLLMiddleHeight*(eigenValues[2]-eigenValues[0])/eigenValues[2]);
        steadyStateWave[0] = std::min(steadyStateWave[0], hRight);
      }

      //Limit the effect of the source term
      //  \cite[ch. 6.4.2]{george2006finite}
      steadyStateWave[1] = std::min(steadyStateWave[1], g*std::max(-hLeft*(bRight-bLeft) , -hRight*(bRight-bLeft)));
      steadyStateWave[1] = std::max(steadyStateWave[1], g*std::min(-hLeft*(bRight-bLeft) , -hRight*(bRight-bLeft)));

      //no source term in the case of a wall
      if( wetDryState == WetDryWall || wetDryState == DryWetWall ) {
        assert(std::fabs(steadyStateWave[0])<zeroTol);
        assert(std::fabs(steadyStateWave[1])<zeroTol);
      }

      rightHandSide[0] -= steadyStateWave[0];
      //rightHandSide[1]: no source term
      rightHandSide[2] -= steadyStateWave[1];

      //everything is ready, solve the equations!
      T beta[3];
      solveLinearEquation(eigenVectors, rightHandSide, beta);

      //compute f-waves and wave-speeds
      if (wetDryState == WetDryWall) { //zero ghost updates (wall boundary)
        //care about the left going wave (0) only
        o_fWaves[0][0] = beta[0] * eigenVectors[1][0];
        o_fWaves[0][1] = beta[0] * eigenVectors[2][0];

        //set the rest to zero
        o_fWaves[1][0] = o_fWaves[1][1] = (T)0.;
        o_fWaves[2][0] = o_fWaves[2][1] = (T)0.;

        o_waveSpeeds[0] = eigenValues[0];
        o_waveSpeeds[1] = o_waveSpeeds[2] = (T)0.;

        assert(eigenValues[0] < zeroTol);
      }
      else if (wetDryState == DryWetWall) { //zero ghost updates (wall boundary)
        //care about the right going wave (2) only
        o_fWaves[2][0] = beta[2] * eigenVectors[1][2];
        o_fWaves[2][1] = beta[2] * eigenVectors[2][2];

        //set the rest to zero
        o_fWaves[0][0] = o_fWaves[0][1] = (T)0.;
        o_fWaves[1][0] = o_fWaves[1][1] = (T)0.;

        o_waveSpeeds[2] = eigenValues[2];
        o_waveSpeeds[0] = o_waveSpeeds[1] = 0.;

        assert(eigenValues[2] > -zeroTol);
      }
      else {
        //compute f-waves (default)
        for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
          o_fWaves[waveNumber][0] = beta[waveNumber] * eigenVectors[1][waveNumber];  //select 2nd and
          o_fWaves[waveNumber][1] = beta[waveNumber] * eigenVectors[2][waveNumber];  //3rd component of the augmented decomposition
#if CONFIG_TSUNAMI_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
          // store eigen-coefficients
          o_eigenCoefficients[waveNumber] = beta[waveNumber];
#endif
        }

        o_waveSpeeds[0] = eigenValues[0];
        o_waveSpeeds[1] = eigenValues[1];
        o_waveSpeeds[2] = eigenValues[2];
      }
    }

    /**
     * Computes the middle state of the homogeneous Riemann-problem.
     *   -> (\cite[ch. 13]{leveque2002finite})
     *
     * @param i_hLeft height on the left side of the edge.
     * @param i_hRight height on the right side of the edge.
     * @param i_uLeft velocity on the left side of the edge.
     * @param i_uRight velocity on the right side of the edge.
     * @param i_huLeft momentum on the left side of the edge.
     * @param i_huRight momentum on the right side of the edge.
     * @param i_maxNumberOfNewtonIterations maximum number of Newton iterations.
     */
    inline void computeMiddleState( const T &i_hLeft, const T &i_hRight,
                                    const T &i_uLeft, const T &i_uRight,
                                    const T &i_huLeft, const T &i_huRight,
                                    const int &i_maxNumberOfNewtonIterations = 1 ) {

      //set everything to zero
      hMiddle = (T)0.;
      middleStateSpeeds[0] = (T)0.;
      middleStateSpeeds[1] = (T)0.;

      //compute local square roots
      //(not necessarily the same ones as in computeNetUpdates!)
      T l_sqrt_g_hRight = std::sqrt(g*i_hRight);
      T l_sqrt_g_hLeft = std::sqrt(g*i_hLeft);

      //single rarefaction in the case of a wet/dry interface
      if (i_hLeft < dryTol) {
        middleStateSpeeds[1] = middleStateSpeeds[0] = i_uRight - (T)2*l_sqrt_g_hRight;
        riemannStructure = DrySingleRarefaction;
        return;
      }
      else if (i_hRight < dryTol) {
        middleStateSpeeds[0] = middleStateSpeeds[1] = i_uLeft + (T)2*l_sqrt_g_hLeft;
        riemannStructure = SingleRarefactionDry;
        return;
      }

      //determine the wave structure of the Riemann-problem
      riemannStructure = determineRiemannStructure( i_hLeft, i_hRight,
                                                    i_uLeft, i_uRight );

      //will be computed later
      T sqrt_g_hMiddle = 0.;

      if (riemannStructure == ShockShock) {
        /*compute All-Shock Riemann Solution
         *  \cite[ch. 13.7]{leveque2002finite}
         *
         *compute middle state h_m
         * => Solve non-linear scalar equation by Newton's method
         *    u_l - (h_m - h_l) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_l} \right) }
         * =  u_r + (h_m - h_r) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_r} \right) }
         *
         * Therefore determine the root of phi_{ss}:
         *\begin{equation}
         *  \begin{matrix}
         *    \phi_{ss}(h) =&&  u_r - u_l +
         *                  (h-h_l) \sqrt{
         *                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right)
         *                          } \\
         *             &&  +(h-h_r) \sqrt{
         *                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right)
         *                          }
         *  \end{matrix}
         *\end{equation}
         *
         *\begin{equation}
         *  \begin{matrix}
         *    \phi_{ss}'(h) = &&\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right) }
         *                      +\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right) }\\
         *                    && -\frac{g}{4}
         *                        \frac{h-h_l}
         *                        {
         *                         h^2\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right) }
         *                        }-
         *                        \frac{g}{4}
         *                        \frac{h-h_r}
         *                        {
         *                         h^2\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right) }
         *                        }
         *  \end{matrix}
         *\end{equation}
         */

//        hMiddle = (hLeft + hRight)*(T)0.618; //first estimate
        hMiddle = std::min(i_hLeft, i_hRight);

        T l_sqrtTermH[2] = {0,0};

        for (int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
//          sqrtTermHLow = std::sqrt((T).5 * g * ((T)1/hMiddle + (T)1/hLeft));
          l_sqrtTermH[0] = std::sqrt((T).5 * g * ((hMiddle+i_hLeft)/(hMiddle*i_hLeft)));
//          sqrtTermHHigh = std::sqrt((T).5 * g * ((T)1/hMiddle + (T)1/hRight));
          l_sqrtTermH[1] = std::sqrt((T).5 * g * ((hMiddle+i_hRight)/(hMiddle*i_hRight)));

          T phi = i_uRight - i_uLeft + (hMiddle - i_hLeft) * l_sqrtTermH[0] + (hMiddle - i_hRight) * l_sqrtTermH[1];

          if (std::fabs(phi) < newtonTol)
            break;

          T derivativePhi = l_sqrtTermH[0] + l_sqrtTermH[1]
                          - (T).25*g*(hMiddle - i_hLeft)  / (l_sqrtTermH[0]*hMiddle*hMiddle)
                          - (T).25*g*(hMiddle - i_hRight) / (l_sqrtTermH[1]*hMiddle*hMiddle);

          hMiddle = hMiddle - phi/derivativePhi; //Newton step
          assert(hMiddle >= dryTol);

//          if(i == i_maxNumberOfNewtonIterations-1) {
//            std::cerr << "Newton-Method did not converge" << std::endl;
//            std::cerr << "std::fabs(phi): " << std::fabs(phi) << std::endl;
//            std::cerr << "hMiddle: " << hMiddle << std::endl;
//            assert(false);
//          }
        }

        sqrt_g_hMiddle = std::sqrt(g*hMiddle);

        //compute middle speed u_m
        //\begin{equation}
        //  \label{eq:hmshock1stfamsimp}
        //  u_m = u_l - (h_m - h_l) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_l} \right) }
        //\end{equation}
        //
        //\begin{equation}
        //  \label{eq:hmshock2ndfamsimp}
        //  u_m = u_r + (h_m - h_r) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_r} \right) }
        //\end{equation}

//        T uMiddleEstimates[2];
//        uMiddleEstimates[0] = uLeft - (hMiddle - hLeft) * l_sqrtTermH[0];
//        uMiddleEstimates[1] = i_uRight + (hMiddle - hRight) * l_sqrtTermH[1];
//        uMiddle = (T).5 * (uMiddleEstimates[0] + uMiddleEstimates[1]);

        //middle state speeds as they are implemented in clawpack, TODO: why?
//        middleStateSpeeds[0] = uMiddleEstimates[0] - sqrt_g_hMiddle;
//        middleStateSpeeds[1] = uMiddleEstimates[1] + sqrt_g_hMiddle;
      }

      if (riemannStructure == RarefactionRarefaction) {
        //Compute All-Rarefaction Riemann Solution
        //  \cite[ch. 13.8.6]{leveque2002finite}

        //compute middle state height h_m
        hMiddle = std::max((T)0., i_uLeft - i_uRight + (T)2.*(l_sqrt_g_hLeft + l_sqrt_g_hRight)); //std::max -> Text after formula (13.56), page 279
        hMiddle = (T)1/((T)16*g)* hMiddle * hMiddle;

        sqrt_g_hMiddle = std::sqrt(g*hMiddle);

        //middle state speeds as they are implemented in clawpack, why?
//        middleStateSpeeds[0] = uLeft + (T)2.*l_sqrt_g_hLow - (T)3.*sqrt_g_hMiddle;
//        middleStateSpeeds[1] = i_uRight - (T)2.*l_sqrt_g_hHigh + (T)3.*sqrt_g_hMiddle;
      }

      if (riemannStructure == ShockRarefaction || riemannStructure == RarefactionShock) {
        //compute dam-break Riemann-solution
        //TODO: reference
        T hMin, hMax;
        if (riemannStructure == ShockRarefaction) {
          hMin = i_hLeft;
          hMax = i_hRight;
        }
        else {
          hMin = i_hRight;
          hMax = i_hLeft;
        }

        //hMiddle = (hLeft + hRight)*(T)0.618; //first estimate
        hMiddle = hMin;

        sqrt_g_hMiddle = std::sqrt(g*hMiddle);
        T sqrt_g_hMax = std::sqrt(g*hMax);
        for (int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
          /*
           *compute middle state h_m
           * => Solve non-linear scalar equation by Newton's method
           *
           * Therefore determine the root of phi_{sr}/phi_{rs}:
           *\begin{equation}
           *  \begin{matrix}
           *    h_{min} = \min(h_l, h_r)\\
           *    h_{max} = \max(h_l, h_r)
           *  \end{matrix}
           *\end{equation}
           *\begin{equation}
           *  \begin{matrix}
           *    \phi_{sr}(h) = \phi_{rs}(h) =&& u_r - u_l + (h - h_{min}) \sqrt{
           *                                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_{min}} \right)
           *                                          } \\
           *               && + 2 (\sqrt{gh} - \sqrt{gh_{max}})
           *  \end{matrix}
           *\end{equation}
           *\begin{equation}
           *    \phi'_{sr}(h) = \phi'_{rs}(h) = \sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_{min}} \right)  }
           *                     -\frac{g}{4} \frac{h-h_{min}}
           *                        { h^2
           *                          \sqrt{
           *                            \frac{g}{2}  \left(\frac{1}{h} + \frac{1}{h_{min}} \right)
           *                          }
           *                        }
           *                     + \sqrt{\frac{g}{h}}
           *\end{equation}
           */

          T sqrtTermHMin = std::sqrt( (T).5*g* ((hMiddle+hMin)/(hMiddle*hMin)) );

          T phi = i_uRight - i_uLeft + (hMiddle - hMin) * sqrtTermHMin + (T)2*(sqrt_g_hMiddle - sqrt_g_hMax);

          if (std::fabs(phi) < newtonTol)
            break;

          T derivativePhi = sqrtTermHMin - (T).25*g * (hMiddle - hMin) / (hMiddle * hMiddle * sqrtTermHMin) + sqrt_g/sqrt_g_hMiddle;

          hMiddle = hMiddle - phi/derivativePhi; //Newton step

          #ifndef NDEBUG
          if (hMiddle < hMin) {
            std::cout << phi << std::endl;
            std::cout << derivativePhi << std::endl;
            std::cerr << "hMiddle(" << hMiddle << ") < hMin(" << hMin << ")" << std::endl;
            assert(false);
          }
          #endif

//          if(i == i_maxNumberOfNewtonIterations-1) {
//            std::cerr << "Newton-Method did not converge" << std::endl;
//            std::cerr << "std::fabs(phi): " << std::fabs(phi) << std::endl;
//            std::cerr << "hMiddle: " << hMiddle << std::endl;
//            assert(false);
//          }

          sqrt_g_hMiddle = std::sqrt(g*hMiddle);
        }

//        //middle state speeds as they are implemented in clawpack, TODO: why?
//        if (riemannStructure == ShockRarefaction) {
//          uMiddle = i_uRight - (T)2.*( l_sqrt_g_hHigh - sqrt_g_hMiddle );
//          middleStateSpeeds[0] = i_uRight - (T)2.*l_sqrt_g_hHigh + sqrt_g_hMiddle;
//          middleStateSpeeds[1] = i_uRight - (T)2.*l_sqrt_g_hHigh + (T)3.*sqrt_g_hMiddle;
//        }
//        else {
//          uMiddle = uLeft + (T)2.*( l_sqrt_g_hLow - sqrt_g_hMiddle );
//          middleStateSpeeds[0] = uLeft + (T)2.*l_sqrt_g_hLow - (T)3.*sqrt_g_hMiddle;
//          middleStateSpeeds[1] = uLeft + (T)2.*l_sqrt_g_hLow - sqrt_g_hMiddle;
//        }
      }

      middleStateSpeeds[0] = i_uLeft + (T)2.*l_sqrt_g_hLeft - (T)3.*sqrt_g_hMiddle;
      middleStateSpeeds[1] = i_uRight - (T)2.*l_sqrt_g_hRight + (T)3.*sqrt_g_hMiddle;

      assert(hMiddle >= 0);
    }

    /**
     * Determine the Riemann-structure of a given problem.
     *   -> \cite[theorem 4.2]{george2006finite}, \cite[appendix B]{george2008augmented}
     *
     * @param i_hLeft height on the left side of the edge.
     * @param i_hRight height on the right side of the edge.
     * @param i_uLeft velocity on the left side of the edge.
     * @param i_uRight velocity on the right side of the edge.
     *
     * @return Riemann-structure of a given problem.
     */
    inline RiemannStructure determineRiemannStructure( const T &i_hLeft, const T &i_hRight,
                                                       const T &i_uLeft, const T &i_uRight ) const {
      T hMin = std::min(i_hLeft, i_hRight);
      T hMax = std::max(i_hLeft, i_hRight);

      T uDif = i_uRight - i_uLeft;

      if( 0 <= (T)2.0*(std::sqrt(g*hMin) - std::sqrt(g*hMax)) + uDif)
        return RarefactionRarefaction;

      if( (hMax - hMin) * std::sqrt(g*(T)0.5*( 1/hMax + 1/hMin )) + uDif <= 0 )
        return ShockShock;

      if( i_hLeft < i_hRight )
        return ShockRarefaction;

      return RarefactionShock;
    }

    /**
     * Solve the linear equation:
     * A * x = b with A \in \mathbb{R}^{3\times3}, x,b \in \mathbb{R}^3
     *
     * @param i_matrix the matrix
     * @param i_b right hand side
     * @param o_x solution
     */
    static inline void solveLinearEquation( const T i_matrix[3][3],
                                            const T i_b[3],
                                            T o_x[3] ) {
//#if 1
      // compute inverse of 3x3 matrix
      const T m[3][3] =
        {	{	 (i_matrix[1][1]*i_matrix[2][2] - i_matrix[1][2]*i_matrix[2][1]), -(i_matrix[0][1]*i_matrix[2][2] - i_matrix[0][2]*i_matrix[2][1]),  (i_matrix[0][1]*i_matrix[1][2] - i_matrix[0][2]*i_matrix[1][1]) },
        {   -(i_matrix[1][0]*i_matrix[2][2] - i_matrix[1][2]*i_matrix[2][0]),  (i_matrix[0][0]*i_matrix[2][2] - i_matrix[0][2]*i_matrix[2][0]), -(i_matrix[0][0]*i_matrix[1][2] - i_matrix[0][2]*i_matrix[1][0]) },
        {    (i_matrix[1][0]*i_matrix[2][1] - i_matrix[1][1]*i_matrix[2][0]), -(i_matrix[0][0]*i_matrix[2][1] - i_matrix[0][1]*i_matrix[2][0]),  (i_matrix[0][0]*i_matrix[1][1] - i_matrix[0][1]*i_matrix[1][0]) }
      };
      T d = (i_matrix[0][0]*m[0][0] + i_matrix[0][1]*m[1][0] + i_matrix[0][2]*m[2][0]);

#ifndef NDEBUG
      if (std::fabs(d) < 0.000000001) {
        std::cerr << "Division close to zero!" << std::endl;
        std::cout << "Matrix: " << std::endl;
        std::cout << i_matrix[0][0] << ", " << i_matrix[0][1] << ", " << i_matrix[0][2] << std::endl;
        std::cout << i_matrix[1][0] << ", " << i_matrix[1][1] << ", " << i_matrix[1][2] << std::endl;
        std::cout << i_matrix[2][0] << ", " << i_matrix[2][1] << ", " << i_matrix[2][2] << std::endl;
        std::cout << "b: " << std::endl;
        std::cout << i_b[0] << ", " << i_b[1] << ", " << i_b[2] << std::endl;
        std::cout << "d: " << d << std::endl;
        assert(false);
      }
#endif

      // m stores not really the inverse matrix, but the inverse multiplied by d
      T s = (T)1/d;

      // compute m*i_b
      o_x[0] = (m[0][0]*i_b[0] + m[0][1]*i_b[1] + m[0][2]*i_b[2])*s;
      o_x[1] = (m[1][0]*i_b[0] + m[1][1]*i_b[1] + m[1][2]*i_b[2])*s;
      o_x[2] = (m[2][0]*i_b[0] + m[2][1]*i_b[1] + m[2][2]*i_b[2])*s;

//      return;
//#else
//      T origDet = computeDeterminant(i_matrix);
//
//      if (std::fabs(origDet) > zeroTol) {
//        T modifiedMatrix[3][3];
//
//        for (int column = 0; column < 3; column++)
//        {
//          memcpy(modifiedMatrix, i_matrix, sizeof(T)*3*3);
//
//          //set a column of the matrix to i_b
//          modifiedMatrix[0][column] = i_b[0];
//          modifiedMatrix[1][column] = i_b[1];
//          modifiedMatrix[2][column] = i_b[2];
//
//          o_x[column] = computeDeterminant(modifiedMatrix) / origDet;
//        }
//      }
//      else {
//        std::cerr << "Warning: Linear dependent eigenvectors! (using Jacobi Solver)" << std::endl;
//        std::cerr << "Determinant: " << origDet << std::endl;
//        std::cerr << i_matrix[0][0] << "\t" << i_matrix[0][1] << "\t" << i_matrix[0][2] << std::endl;
//        std::cerr << i_matrix[1][0] << "\t" << i_matrix[1][1] << "\t" << i_matrix[1][2] << std::endl;
//        std::cerr << i_matrix[2][0] << "\t" << i_matrix[2][1] << "\t" << i_matrix[2][2] << std::endl;
//#if 0
//        T xTemp[3];
//        xTemp[0] = xTemp[1] = xTemp[2] = 0.;
//        for (int m = 0; m < 20; m++) {
//          for (int row = 0; row < 3; row++) {
//            o_x[row] = 0.;
//            for (int col = 0; col < 3; col++) {
//              if (col != row)
//                o_x[row] += i_matrix[row][col]*xTemp[col];
//            }
//            o_x[row] = (i_b[row] - o_x[row])/i_matrix[row][row];
//          }
//
//          if (fabs(o_x[0]-xTemp[0]) + fabs(o_x[1]-xTemp[1]) + fabs(o_x[2]-xTemp[2]) < zeroTol*10.)
//            break;
//          else {
//            std::cout << "Error: " << fabs(o_x[0]-xTemp[0]) + fabs(o_x[1]-xTemp[1]) + fabs(o_x[2]-xTemp[2]) << std::endl;
//            xTemp[0] = o_x[0];
//            xTemp[1] = o_x[1];
//            xTemp[2] = o_x[2];
//          }
//        }
//        std::cout << "Solution:" << std::endl;
//        std::cout << "\t" << o_x[0] << std::endl;
//        std::cout << "x=\t" << o_x[1] << std::endl;
//        std::cout << "\t" << o_x[2] << std::endl;
//        std::cout << "*********" << std::endl;
//        std::cout << "\t" << i_b[0] << std::endl;
//        std::cout << "b=\t" << i_b[1] << std::endl;
//        std::cout << "\t" << i_b[2] << std::endl;
//        std::cout << "***A*x****" << std::endl;
//        for (int row = 0; row < 3; row++) {
//          std::cout << "\t" <<
//              i_matrix[row][0] * o_x[0] +
//              i_matrix[row][1] * o_x[1] +
//              i_matrix[row][2] * o_x[2] << std::endl;
//        }
//#endif
//        assert(false);
//      }
//#endif
    }
//#if 0
//    /**
//     * Compute the determinant of a 3x3 matrix
//     * using formula: |A|=a11 * (a22a33-a32a23) + a12 * (a23a31-a33a21) + a13 * (a21a32-a31a22)
//     *                |A|=a00 * (a11a22-a21a12) + a01 * (a12a20-a22a10) + a02 * (a10a21-a20a11)
//     */
//    T computeDeterminant(T matrix[3][3]) const
//    {
//
//      T determinant =     matrix[0][0] * ( matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2] )
//                        		    + matrix[0][1] * ( matrix[1][2] * matrix[2][0] - matrix[2][2] * matrix[1][0] )
//                        		    + matrix[0][2] * ( matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1] );
//
//      return determinant;
//    }
//#endif


#undef sqrt_g_hLeft
#undef sqrt_g_hRight
#undef sqrt_hLeft
#undef sqrt_hRight

};

#undef complexsteadystatewave
#undef correctrarefactions


#endif /* AUGRIE_HPP_ */
