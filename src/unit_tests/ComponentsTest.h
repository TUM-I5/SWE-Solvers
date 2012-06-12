/**
 * ComponentsTest.h
 *
 ****
 **** Collection of elementary unit tests for the Augmented Riemann solver.
 ****
 *
 *  Created on: Nov 30, 2011
 *  Last Update: June 12, 2012
 *
 ****
 *
 *  Author: Alexander Breuer
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer
 *    E-Mail: breuera AT in.tum.de
 *
 ****
 */

#ifndef COMPONENTS_HPP_
#define COMPONENTS_HPP_

#define private public
#define protected public
#include "../solver/AugRie.hpp"
#undef private
#undef protected

#include <string>

class ComponentsTest {
  //private:
  //! gravity constant.
  const double gravity;
  //! numerical definition of "dry".
  const double dryTol;
  //! numerical defition of zero.
  const double zeroTol;
  //! maximum water height within the tests.
  const double maxWaterHeight;
  //! maximum wave speed within the tests.
  //850km/h = 14166.66 m/s
  const double maxWaveSpeed;
  //! maximum variation of the speed within the tests.
  //50km/h = 833.33 m/s
  const double maxWaveSpeedVariation;
  //! how many of the tests with random numbers should run.
  const int numberOfRandomTests;
  //! number of precomputed values available for the tests of the function phi.
  const int testPhiSize;
  //! number of phi-values, which correspond to "walls", focus: Inundation.
  const int testPhiWallSamples;
  //! number of phi-values, which correspond to random values.
  const int testPhiRandomSamples;
  //! maximum number of newton iterations in a regular setup.
  const int maxNumberOfNewtonIterationsDefault;
  //! maximum number of newton iterations in an expensive setup, focus: test for convergency.
  const int maxNumberOfNewtonIterationsExpensive;
  //! newton tolerance (when to stop iterating).
  const double newtonTol;

  //! path of the test file, which contains precomputed values.
  const std::string testFileName;

  //! Augmented Riemann solver.
  solver::AugRie<double> mySolver;

  //! water heights and velocities of the Riemann problem.
  double hLeft, hRight, uLeft, uRight;

  //! array of precomputed test values.
  double** testPhi;

  //functions
  double createRandomNumber(const double min, const double max);
  void readTestPhi();
  std::string printSolverVariables();
  std::string printSolverVariablesWithResults( const double hUpdateLeft, const double hUpdateRight,
                                               const double huUpdateLeft, const double huUpdateRight,
                                               const double maxWaveSpeed );


  public:
  ComponentsTest();
  ~ComponentsTest();

  //functions
  void operator() ();

  void checkDetermineRiemannStructure();
  void checkComputeMiddleState();
  void checkDetermineWetDryStateSimple();
  void checkInundationLimitsInOneDirection(const bool i_leftCellDry);
  void checkInundationLimits();
  void checkComputeNetUpdatesForWallsInOneDirection(const bool i_leftCellDry);
  void checkComputeNetUpdatesForWalls();
  void checkComputeNetUpdatesForPositivity();
};


#endif /* COMPONENTS_HPP_ */
