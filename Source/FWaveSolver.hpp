/**
 * FWaveSolver.hpp
 *
 ****
 **** F-Wave Riemann Solver for the Shallow Water Equation
 ****
 *
 *  Created on: Aug 25, 2011
 *  Last Update: Feb 18, 2012
 *
 ****
 *
 *  Author: Alexander Breuer
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer
 *    E-Mail: breuera AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @article{bale2002wave,
 *            title={A wave propagation method for conservation laws and balance laws with spatially varying flux
 *functions}, author={Bale, D.S. and LeVeque, R.J. and Mitran, S. and Rossmanith, J.A.}, journal={SIAM Journal on
 *Scientific Computing}, volume={24}, number={3}, pages={955--978}, year={2002}, publisher={Citeseer}}
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

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "WavePropagationSolver.hpp"

namespace Solvers {

  /**
   * FWave Riemann Solver for the Shallow Water Equations.
   *
   * T should be double or float.
   */
  template <class T>
  class FWaveSolver: public WavePropagationSolver<T> {
  private:
    // Use nondependent names (template base class)
    using WavePropagationSolver<T>::dryTol_;
    using WavePropagationSolver<T>::gravity_;
    using WavePropagationSolver<T>::zeroTol_;

    using WavePropagationSolver<T>::hLeft_;
    using WavePropagationSolver<T>::hRight_;
    using WavePropagationSolver<T>::huLeft_;
    using WavePropagationSolver<T>::huRight_;
    using WavePropagationSolver<T>::bLeft_;
    using WavePropagationSolver<T>::bRight_;
    using WavePropagationSolver<T>::uLeft_;
    using WavePropagationSolver<T>::uRight_;

    using WavePropagationSolver<T>::wetDryState_;

    /**
     * Compute the edge local eigenvalues.
     *
     * @param o_waveSpeeds will be set to: speeds of the linearized waves (eigenvalues).
     */
    void computeWaveSpeeds(T o_waveSpeeds[2]) const {
      // Compute eigenvalues of the Jacobian matrices in states Q_{i-1} and Q_{i}
      T characteristicSpeeds[2]{};
      characteristicSpeeds[0] = uLeft_ - std::sqrt(gravity_ * hLeft_);
      characteristicSpeeds[1] = uRight_ + std::sqrt(gravity_ * hRight_);

      // Compute "Roe speeds"
      T hRoe = T(0.5) * (hRight_ + hLeft_);
      T uRoe = uLeft_ * std::sqrt(hLeft_) + uRight_ * std::sqrt(hRight_);
      uRoe /= std::sqrt(hLeft_) + std::sqrt(hRight_);

      T roeSpeeds[2]{};
      roeSpeeds[0] = uRoe - std::sqrt(gravity_ * hRoe);
      roeSpeeds[1] = uRoe + std::sqrt(gravity_ * hRoe);

      // Compute eindfeldt speeds
      T einfeldtSpeeds[2]{};
      einfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
      einfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);

      // Set wave speeds
      o_waveSpeeds[0] = einfeldtSpeeds[0];
      o_waveSpeeds[1] = einfeldtSpeeds[1];
    }

    /**
     * Compute the decomposition into f-Waves.
     *
     * @param waveSpeeds speeds of the linearized waves (eigenvalues).
     * @param o_fWaves will be set to: Decomposition into f-Waves.
     */
    void computeWaveDecomposition(const T waveSpeeds[2], T o_fWaves[2][2]) const {
      // Eigenvalues***********************************************************************************************
      // Computed somewhere before.
      // An option would be to use the char. Speeds:
      //
      // lambda^1 = u_{i-1} - sqrt(g * h_{i-1})
      // lambda^2 = u_i     + sqrt(g * h_i)
      // Matrix of right eigenvectors******************************************************************************
      //     1                              1
      // R =
      //     u_{i-1} - sqrt(g * h_{i-1})    u_i + sqrt(g * h_i)
      // **********************************************************************************************************
      //                                                                      u_i + sqrt(g * h_i)              -1
      // R^{-1} = 1 / (u_i - sqrt(g * h_i) - u_{i-1} + sqrt(g * h_{i-1}) *
      //                                                                   -( u_{i-1} - sqrt(g * h_{i-1}) )     1
      // **********************************************************************************************************
      //                hu
      // f(q) =
      //         hu^2 + 1/2 g * h^2
      // **********************************************************************************************************
      //                                    0
      // \delta x \Psi =
      //                  -g * 1/2 * (h_i + h_{i-1}) * (b_i - b_{i+1})
      // **********************************************************************************************************
      // beta = R^{-1} * (f(Q_i) - f(Q_{i-1}) - \delta x \Psi)
      // **********************************************************************************************************

      // assert: wave speed of the 1st wave family should be less than the speed of the 2nd wave family.
      assert(waveSpeeds[0] < waveSpeeds[1]);

      T lambdaDif = waveSpeeds[1] - waveSpeeds[0];

      // assert: no division by zero
      assert(std::abs(lambdaDif) > zeroTol_);

      // Compute the inverse matrix R^{-1}
      T Rinv[2][2]{};

      T oneDivLambdaDif = T(1.0) / lambdaDif;
      Rinv[0][0]        = oneDivLambdaDif * waveSpeeds[1];
      Rinv[0][1]        = -oneDivLambdaDif;

      Rinv[1][0] = oneDivLambdaDif * -waveSpeeds[0];
      Rinv[1][1] = oneDivLambdaDif;

      // Right hand side
      T fDif[2]{};

      // Calculate modified (bathymetry!) flux difference
      // f(Q_i) - f(Q_{i-1})
      fDif[0] = huRight_ - huLeft_;
      fDif[1] = huRight_ * uRight_ + T(0.5) * gravity_ * hRight_ * hRight_
                - (huLeft_ * uLeft_ + T(0.5) * gravity_ * hLeft_ * hLeft_);

      // \delta x \Psi[2]
      T psi = -gravity_ * T(0.5) * (hRight_ + hLeft_) * (bRight_ - bLeft_);
      fDif[1] -= psi;

      // Solve linear equations
      T beta[2]{};
      beta[0] = Rinv[0][0] * fDif[0] + Rinv[0][1] * fDif[1];
      beta[1] = Rinv[1][0] * fDif[0] + Rinv[1][1] * fDif[1];

      // Return f-waves
      o_fWaves[0][0] = beta[0];
      o_fWaves[0][1] = beta[0] * waveSpeeds[0];

      o_fWaves[1][0] = beta[1];
      o_fWaves[1][1] = beta[1] * waveSpeeds[1];
    }

    /**
     * Compute net updates for the cell on the left/right side of the edge.
     * Its assumed that the member variables are set already.
     *
     * @param waveSpeeds speeds of the linearized waves (eigenvalues).
     *
     * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
     * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
     * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
     * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
     * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
     */
    void computeNetUpdatesWithWaveSpeeds(
      const T waveSpeeds[2],
      T&      o_hUpdateLeft,
      T&      o_hUpdateRight,
      T&      o_huUpdateLeft,
      T&      o_huUpdateRight,
      T&      o_maxWaveSpeed
    ) {
      // Reset net updates
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = T(0.0);

      //! Where to store the two f-waves
      T fWaves[2][2];

      // Compute the decomposition into f-waves
      computeWaveDecomposition(waveSpeeds, fWaves);

      // Compute the net-updates
      // 1st wave family
      if (waveSpeeds[0] < -zeroTol_) { // Left going
        o_hUpdateLeft += fWaves[0][0];
        o_huUpdateLeft += fWaves[0][1];
      } else if (waveSpeeds[0] > zeroTol_) { // Right going
        o_hUpdateRight += fWaves[0][0];
        o_huUpdateRight += fWaves[0][1];
      } else { // Split waves
        o_hUpdateLeft += T(0.5) * fWaves[0][0];
        o_huUpdateLeft += T(0.5) * fWaves[0][1];
        o_hUpdateRight += T(0.5) * fWaves[0][0];
        o_huUpdateRight += T(0.5) * fWaves[0][1];
      }

      // 2nd wave family
      if (waveSpeeds[1] < -zeroTol_) { // Left going
        o_hUpdateLeft += fWaves[1][0];
        o_huUpdateLeft += fWaves[1][1];
      } else if (waveSpeeds[1] > zeroTol_) { // Right going
        o_hUpdateRight += fWaves[1][0];
        o_huUpdateRight += fWaves[1][1];
      } else { // Split waves
        o_hUpdateLeft += T(0.5) * fWaves[1][0];
        o_huUpdateLeft += T(0.5) * fWaves[1][1];
        o_hUpdateRight += T(0.5) * fWaves[1][0];
        o_huUpdateRight += T(0.5) * fWaves[1][1];
      }

      // Compute maximum wave speed (-> CFL-condition)
      o_maxWaveSpeed = std::max(std::fabs(waveSpeeds[0]), std::fabs(waveSpeeds[1]));
    }

  protected:
    /**
     * Determine the wet/dry state and set member variables accordingly.
     */
    void determineWetDryState() override {
      // Determine the wet/dry state
      if (hLeft_ < dryTol_ && hRight_ < dryTol_) { // Both cells are dry
        wetDryState_ = WavePropagationSolver<T>::WetDryState::DryDry;
      } else if (hLeft_ < dryTol_) { // Left cell dry, right cell wet
        uRight_ = huRight_ / hRight_;

        // Set wall boundary conditions.
        // This is not correct in the case of inundation problems.
        hLeft_       = hRight_;
        bLeft_       = bRight_;
        huLeft_      = -huRight_;
        uLeft_       = -uRight_;
        wetDryState_ = WavePropagationSolver<T>::WetDryState::DryWetWall;
      } else if (hRight_ < dryTol_) { // Left cell wet, right cell dry
        uLeft_ = huLeft_ / hLeft_;

        // Set wall boundary conditions.
        // This is not correct in the case of inundation problems.
        hRight_      = hLeft_;
        bRight_      = bLeft_;
        huRight_     = -huLeft_;
        uRight_      = -uLeft_;
        wetDryState_ = WavePropagationSolver<T>::WetDryState::WetDryWall;
      } else { // Both cells wet
        uLeft_  = huLeft_ / hLeft_;
        uRight_ = huRight_ / hRight_;

        wetDryState_ = WavePropagationSolver<T>::WetDryState::WetWet;
      }
    }

  public:
    /**
     * Constructor of the f-Wave solver with optional parameters.
     *
     * @param dryTolerance numerical definition of "dry".
     * @param gravity gravity constant.
     * @param zeroTolerance numerical definition of zero.
     */
    FWaveSolver(
      T dryTolerance  = static_cast<T>(0.01),
      T gravity       = static_cast<T>(9.81),
      T zeroTolerance = static_cast<T>(0.000000001)
    ):
      WavePropagationSolver<T>(dryTolerance, gravity, zeroTolerance) {}

    ~FWaveSolver() override = default;

    /**
     * Compute net updates for the cell on the left/right side of the edge.
     * This is the default method of a standalone f-Wave solver.
     *
     * Please note:
     *   In the case of a Dry/Wet- or Wet/Dry-boundary, wall boundary conditions will be set.
     *   The f-Wave solver is not positivity preserving.
     *   -> You as the programmer have to take care about "negative water heights"!
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param huLeft momentum on the left side of the edge.
     * @param huRight momentum on the right side of the edge.
     * @param bLeft bathymetry on the left side of the edge.
     * @param bRight bathymetry on the right side of the edge.
     *
     * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
     * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
     * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
     * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
     * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
     */
    void computeNetUpdates(
      const T& hLeft,
      const T& hRight,
      const T& huLeft,
      const T& huRight,
      const T& bLeft,
      const T& bRight,
      T&       o_hUpdateLeft,
      T&       o_hUpdateRight,
      T&       o_huUpdateLeft,
      T&       o_huUpdateRight,
      T&       o_maxWaveSpeed
    ) override {
      // Set speeds to zero (will be determined later)
      uLeft_ = uRight_ = 0;

      // Reset the maximum wave speed
      o_maxWaveSpeed = 0;

      //! Wave speeds of the f-waves
      T waveSpeeds[2];

      // Store parameters to member variables
      WavePropagationSolver<T>::storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

      // Determine the wet/dry state and compute local variables correspondingly
      determineWetDryState();

      // Zero updates and return in the case of dry cells
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryDry) {
        o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = T(0.0);
        return;
      }

      // Compute the wave speeds
      computeWaveSpeeds(waveSpeeds);

      // Use the wave speeds to compute the net-updates
      computeNetUpdatesWithWaveSpeeds(
        waveSpeeds, o_hUpdateLeft, o_hUpdateRight, o_huUpdateLeft, o_huUpdateRight, o_maxWaveSpeed
      );

      // Zero ghost updates (wall boundary)
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::WetDryWall) {
        o_hUpdateRight  = 0;
        o_huUpdateRight = 0;
      } else if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWall) {
        o_hUpdateLeft  = 0;
        o_huUpdateLeft = 0;
      }
    }

    /**
     * Compute net updates for the cell on the left/right side of the edge.
     * This is an expert method, because a lot of (numerical-)knowledge about the problem is assumed/has to be provided.
     * It is the f-Wave entry point for the hybrid solver,  which combines the "simple" F-Wave approach with the more
     * complex Augmented Riemann Solver.
     *
     * wetDryState is assumed to be WetWet.
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param huLeft momentum on the left side of the edge.
     * @param huRight momentum on the right side of the edge.
     * @param bLeft bathymetry on the left side of the edge.
     * @param bRight bathymetry on the right side of the edge.
     * @param uLeft velocity on the left side of the edge.
     * @param uRight velocity on the right side of the edge.
     * @param waveSpeeds speeds of the linearized waves (eigenvalues).
     *                   A hybrid solver will typically provide its own values.
     *
     * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
     * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
     * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
     * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
     * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
     */
    void computeNetUpdatesHybrid(
      const T& hLeft,
      const T& hRight,
      const T& huLeft,
      const T& huRight,
      const T& bLeft,
      const T& bRight,
      const T& uLeft,
      const T& uRight,
      const T  waveSpeeds[2],
      T&       o_hUpdateLeft,
      T&       o_hUpdateRight,
      T&       o_huUpdateLeft,
      T&       o_huUpdateRight,
      T&       o_maxWaveSpeed
    ) {
      // Store parameters to member variables
      storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight, uLeft, uRight);

      computeNetUpdatesWithWaveSpeeds(
        waveSpeeds, o_hUpdateLeft, o_hUpdateRight, o_huUpdateLeft, o_huUpdateRight, o_maxWaveSpeed
      );
    }
  };

} // namespace Solvers
