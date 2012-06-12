/**
 * AugRieGeoClaw.hpp
 *
 ****
 **** Augmented Riemann solver which uses the underlying Fortran routines of GeoClaw directly.
 ****
 *
 *  Created on: Mar 13, 2012
 *  Last Update: Mar 14, 2012
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
 *   @webpage{levequeclawpack,
 *            Author = {LeVeque, R. J.},
 *            Lastchecked = {Mar, 14, 2011},
 *            Title = {Clawpack Sofware},
 *            Url = {https://github.com/clawpack/geoclaw}
 *
 *
 ****
 *
 * Acknowledgments:
 *   Special thanks go to R.J. LeVeque and D.L. George for publishing their code.
 */

#ifndef AUGRIEGEOCLAW_HPP_
#define AUGRIEGEOCLAW_HPP_

#include <cassert>
#include <cmath>
#include "WavePropagation.hpp"

/**
 * Extern declaration of the c_bing_geoclaw_riemann_aug_JCP routine.
 *
 * Remark: Not constant variables might change during execution.
 *
 * @param i_maxNumberOfRiemannIterations maximum number Riemann iterations (solver might iterate over the Riemann problem one day, currently fixed to 1)
 * @param i_numberOfFWaves number of fWaves (might change to 2 for standard f-Wave solver for example, currently fixed to 3)
 * @param i_hLeft height on the left side of the edge.
 * @param i_hRight height on the right side of the edge.
 * @param i_huLeft momentum on the left side of the edge in normal direction.
 * @param i_huRight momentum on the right side of the edge in normal direction.
 * @param i_hvLeft momentum on the left side of the edge in parallel direction.
 * @param i_hvRight momentum on the right side of the edge in parallel direction.
 * @param i_bLeft bathymetry on the left side of the edge.
 * @param i_bRight bathymetry on the right side of the edge.
 * @param i_dryTol dry tolerance (definition of wet/dry cells).
 * @param i_g gravity constant (typically 9.81).
 * @param o_netUpdatesLeft will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the cell on the left side of the edge.
 * @param o_netUpdatesRight will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the cell on the right side of the edge.
 * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
 */
extern "C" void c_bind_geoclaw_riemann_aug_JCP( const int &i_maxNumberOfRiemannIterations, const int &i_numberOfFWaves,
                                                double &i_hLeft,   double &i_hRight,
                                                double &i_huLeft,  double &i_huRight,
                                                double &i_hvLeft,  double &i_hvRight,
                                                double &i_bLeft,   double &i_bRight,
                                                const double &i_dryTol, const double &i_g,
                                                double o_netUpdatesLeft[3], double o_netUpdatesRight[3],
                                                double &o_maxWaveSpeed);

namespace solver {
  template <typename T> class AugRieGeoClaw;
}

/**
 *  Augmented Riemann solver which uses the underlying Fortran routines of GeoClaw directly.
 *
 *  T should be double or float. (currently fixed to double precision)
 */
template <typename T> class solver::AugRieGeoClaw: public WavePropagation<T> {
  //private:
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

    using solver::WavePropagation<T>::storeParameters;

    //! momentum on the left side of the edge (could change during execution) in parallel direction.
    T hvLeft;
    //! momentum on the right side of the edge (could change during execution) in parallel direction.
    T hvRight;

    //! number of iterations over the Riemann problem (not used)
    const int maxNumberOfRiemannIterations;

    /**
     * Implemented within the GeoClaw-Fortran routines. Not used here.
     */
    void determineWetDryState() {assert(false);}

    /**
     * Store parameters to member variables.
     *
     * @param i_hLeft height on the left side of the edge.
     * @param i_hRight height on the right side of the edge.
     * @param i_huLeft momentum on the left side of the edge in normal direction.
     * @param i_huRight momentum on the right side of the edge in normal direction.
     * @param i_hvLeft momentum on the left side of the edge in parallel direction.
     * @param i_hvRight momentum on the right side of the edge in parallel direction.
     * @param i_bLeft bathymetry on the left side of the edge.
     * @param i_bRight bathymetry on the right side of the edge.
     */
    inline void storeParameters( const T &i_hLeft,  const T &i_hRight,
                                 const T &i_huLeft, const T &i_huRight,
                                 const T &i_hvLeft, const T &i_hvRight,
                                 const T &i_bLeft,  const T &i_bRight ) {
      //store common parameters to member variables
      storeParameters( i_hLeft,  i_hRight,
                       i_huLeft, i_huRight,
                       i_bLeft,  i_bRight );

      hvLeft = i_hvLeft;
      hvRight = i_hvRight;
    }

  public:
    /**
     * Constructor of the GeoClaw Augmented Riemann solver with optional parameters.
     *
     * @param i_dryTolerance numerical definition of "dry".
     * @param i_gravity gravity constant.
     * @param i_zeroTolerance numerical definition of zero (currently not used).
     * @param i_maxNumberOfRiemannIterations maximum number Riemann iterations (solver might iterate over the Riemann problem one day, currently fixed to 1)
     */
    AugRieGeoClaw( T i_dryTolerance =  (T) 0.01,
                   T i_gravity =       (T) 9.81,
                   T i_zeroTolerance = (T) 0.000000001,
                   int i_maxNumberOfRiemannIterations = 1 ):
                     WavePropagation<T>( i_dryTolerance, i_gravity, i_zeroTolerance ),
                     maxNumberOfRiemannIterations(i_maxNumberOfRiemannIterations){};

    /**
     * Not implemented.
     */
    void computeNetUpdates ( const T &i_hLeft,  const T &i_hRight,
                             const T &i_huLeft, const T &i_huRight,
                             const T &i_bLeft,  const T &i_bRight,

                             T &o_hUpdateLeft,
                             T &o_hUpdateRight,
                             T &o_huUpdateLeft,
                             T &o_huUpdateRight,
                             T &o_maxWaveSpeed ) {
      //not implemented
      assert(false);
    }

    /**
     * Compute net updates for the cell on the left/right side of the edge.
     *
     * @param i_hLeft height on the left side of the edge.
     * @param i_hRight height on the right side of the edge.
     * @param i_huLeft momentum on the left side of the edge in normal direction.
     * @param i_huRight momentum on the right side of the edge in normal direction.
     * @param i_hvLeft momentum on the left side of the edge in parallel direction.
     * @param i_hvRight momentum on the right side of the edge in parallel direction.
     * @param i_bLeft bathymetry on the left side of the edge.
     * @param i_bRight bathymetry on the right side of the edge.
     *
     * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
     * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
     * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge in normal direction.
     * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge in normal direction.
     * @param o_hvUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge in parallel direction.
     * @param o_hvUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge in parallel direction.
     * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
     */
    void computeNetUpdates ( const T &i_hLeft,  const T &i_hRight,
                             const T &i_huLeft, const T &i_huRight,
                             const T &i_hvLeft, const T &i_hvRight,
                             const T &i_bLeft,  const T &i_bRight,

                             T &o_hUpdateLeft,
                             T &o_hUpdateRight,
                             T &o_huUpdateLeft,
                             T &o_huUpdateRight,
                             T &o_hvUpdateLeft,
                             T &o_hvUpdateRight,
                             T &o_maxWaveSpeed ) {
      // store parameters to member variables
      storeParameters( i_hLeft, i_hRight,
                       i_huLeft, i_huRight,
                       i_hvLeft, i_hvRight,
                       i_bLeft, i_bRight );

      double netUpdatesLeft[3], netUpdatesRight[3];

      c_bind_geoclaw_riemann_aug_JCP( maxNumberOfRiemannIterations, 3,
                                      hLeft,   hRight,
                                      huLeft,  huRight,
                                      hvLeft,  hvRight,
                                      bLeft,   bRight,
                                      dryTol, g,
                                      netUpdatesLeft, netUpdatesRight,
                                      o_maxWaveSpeed );

      //copy net-updates
      o_hUpdateLeft = netUpdatesLeft[0];
      o_hUpdateRight = netUpdatesRight[0];

      o_huUpdateLeft = netUpdatesLeft[1];
      o_huUpdateRight = netUpdatesRight[1];

      o_hvUpdateLeft = netUpdatesLeft[2];
      o_hvUpdateRight = netUpdatesRight[2];
    }
};

#endif /* AUGRIEGEOCLAW_HPP_ */
