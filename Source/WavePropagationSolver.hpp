#pragma once

namespace Solvers {

  /**
   * Abstract wave propagation solver for the Shallow Water Equations.
   *
   * T should be double or float.
   */
  template <class T>
  class WavePropagationSolver {
  protected:
    T       dryTol_;  //! Numerical definition of "dry".
    const T gravity_; //! Gravity constant.
    const T zeroTol_; //! Numerical definition of zero.

#if 0
    //! Parameters for computeNetUpdates.
    T h_[2];
    T hu_[2];
    T b_[2];
    T u_[2];

#define hLeft_ (h_[0])
#define hRight_ (h_[1])

#define huLeft_ (hu_[0])
#define huRight_ (hu_[1])

#define bLeft_ (b_[0])
#define bRight_ (b_[1])

#define uLeft_ (u_[0])
#define uRight_ (u_[1])
#else
    //! Edge-local variables.
    T hLeft_;   //! Height on the left side of the edge (could change during execution).
    T hRight_;  //! Height on the right side of the edge (could change during execution).
    T huLeft_;  //! Momentum on the left side of the edge (could change during execution).
    T huRight_; //! Momentum on the right side of the edge (could change during execution).
    T bLeft_;   //! Bathymetry on the left side of the edge (could change during execution).
    T bRight_;  //! Bathymetry on the right side of the edge (could change during execution).
    T uLeft_;   //! Velocity on the left side of the edge (computed by determineWetDryState).
    T uRight_;  //! Velocity on the right side of the edge (computed by determineWetDryState).
#endif

    /**
     * The wet/dry state of the Riemann-problem.
     */
    enum WetDryState {
      DryDry,           /**< Both cells are dry. */
      WetWet,           /**< Both cells are wet. */
      WetDryInundation, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies higher than the 2nd one. */
      WetDryWall, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies lower than the 2nd one. Momentum is not large enough
                     to overcome the difference. */
      WetDryWallInundation, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies lower than the 2nd one. Momentum is large
                               enough to overcome the difference. */
      DryWetInundation,     /**< 1st cell: dry, 2nd cell: wet. 1st cell lies lower than the 2nd one. */
      DryWetWall, /**< 1st cell: dry, 2nd cell: wet. 1st cell lies higher than the 2nd one. Momentum is not large enough
                     to overcome the difference. */
      DryWetWallInundation /**< 1st cell: dry, 2nd cell: wet. 1st cell lies higher than the 2nd one. Momentum is large
                              enough to overcome the difference. */
    };

    WetDryState wetDryState_; //! wet/dry state of our Riemann-problem (determined by determineWetDryState).

    //! Determine the wet/dry-state and set local values if we have to.
    virtual void determineWetDryState() = 0;

    /**
     * Constructor of a wave propagation solver.
     *
     * @param gravity gravity constant.
     * @param dryTolerance numerical definition of "dry".
     * @param zeroTolerance numerical definition of zero.
     */
    WavePropagationSolver(T dryTolerance, T gravity, T zeroTolerance):
      dryTol_(dryTolerance),
      gravity_(gravity),
      zeroTol_(zeroTolerance) {}

    /**
     * Store parameters to member variables.
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param huLeft momentum on the left side of the edge.
     * @param huRight momentum on the right side of the edge.
     * @param bLeft bathymetry on the left side of the edge.
     * @param bRight bathymetry on the right side of the edge.
     */
    void storeParameters(
      const T& hLeft, const T& hRight, const T& huLeft, const T& huRight, const T& bLeft, const T& bRight
    ) {
      hLeft_  = hLeft;
      hRight_ = hRight;

      huLeft_  = huLeft;
      huRight_ = huRight;

      bLeft_  = bLeft;
      bRight_ = bRight;
    }

    /**
     * Store parameters to member variables.
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param huLeft momentum on the left side of the edge.
     * @param huRight momentum on the right side of the edge.
     * @param bLeft bathymetry on the left side of the edge.
     * @param bRight bathymetry on the right side of the edge.
     * @param uLeft velocity on the left side of the edge.
     * @param uRight velocity on the right side of the edge.
     */
    void storeParameters(
      const T& hLeft,
      const T& hRight,
      const T& huLeft,
      const T& huRight,
      const T& bLeft,
      const T& bRight,
      const T& uLeft,
      const T& uRight
    ) {
      storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

      uLeft_  = uLeft;
      uRight_ = uRight;
    }

  public:
    virtual ~WavePropagationSolver() = default;

    /**
     * Compute net updates for the cell on the left/right side of the edge.
     * This is the default method every standalone wave propagation solver should provide.
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
    virtual void computeNetUpdates(
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
#ifdef ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
      ,
      T o_eigenCoefficients[3]
#endif
    ) = 0;

    /**
     * Sets the dry tolerance of the solver.
     *
     * @param dryTolerance dry tolerance.
     */
    void setDryTolerance(const T dryTolerance) { dryTol_ = dryTolerance; }

#undef hLeft
#undef hRight

#undef huLeft
#undef huRight

#undef bLeft
#undef bRight

#undef uLeft
#undef uRight
  };

} // namespace Solvers
