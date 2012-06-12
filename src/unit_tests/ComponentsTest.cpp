/**
 * ComponentsTest.cpp
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
#include "ComponentsTest.h"

//Access private members
#define private public
#include "../solver/AugRie.hpp"
#undef private

#include "cute.h"

#include <cmath>
#include <string>
#include <netcdfcpp.h>

/**
 * Contructor.
 *
 * Sets all variables to reasonable default values and reads the pre-computed values of the function phi.
 */
ComponentsTest::ComponentsTest():
  gravity(9.81),
  dryTol(0.001),
  zeroTol(0.00001),
  maxWaterHeight(8500),
  maxWaveSpeed(14166.66),
  maxWaveSpeedVariation(833.33),
  numberOfRandomTests(500000),
  testPhiSize(1000000),
  testPhiWallSamples(500000),
  testPhiRandomSamples(500000),
  maxNumberOfNewtonIterationsDefault(10),
  maxNumberOfNewtonIterationsExpensive(50000),
  newtonTol(0.0001),
  testFileName("unit_tests/files/testphi_1000000.nc"),
  mySolver(dryTol, gravity, newtonTol, maxNumberOfNewtonIterationsDefault, zeroTol) {

  // read pre-computed values of phi
  testPhi = new double* [5];
  for(int i = 0; i < 5; i++) {
    testPhi[i] = new double[testPhiSize];
  }
  readTestPhi();
}

/**
 * Destructor.
 *
 */
ComponentsTest::~ComponentsTest() {
//  for(int i = 0; i < 5; i++) {
//    delete[] testPhi[i];
//  }
//  delete[] testPhi;
  testPhi = NULL;
}

/**
 * Operator with no arguments, which is called in the main function.
 */
void ComponentsTest::operator() (){
  checkDetermineRiemannStructure();
  checkComputeMiddleState();
  checkDetermineWetDryStateSimple();
  checkInundationLimits();
  checkComputeNetUpdatesForWalls();
  checkComputeNetUpdatesForPositivity();
};

/**
 * Create a random number between l_min and l_max.
 *
 * @param i_min minimum value.
 * @param i_max maximum value.
 * @return random number.
 */
double ComponentsTest::createRandomNumber( const double i_min, const double i_max) {
    double randomNumber = ((double)rand()/(double)RAND_MAX); //between 0 and 1
    return i_min + (i_max-i_min) * randomNumber;
}

/**
 * Read a test file with pre-computed solutions for
 * the non-linear algebraic funtion phi.
 */
void ComponentsTest::readTestPhi() {
  NcFile dataFile(testFileName.c_str(), NcFile::ReadOnly);
  if (!dataFile.is_valid()) {
    std::cout << "Couldn't open file: " << testFileName << std::endl;
    assert(false);
  }
  ASSERT( dataFile.get_dim(0)->size() == testPhiSize );

  NcVar* data;
  //Retrieve the variables
  data = dataFile.get_var("hLow");
  data->get(testPhi[0], testPhiSize, 0);

  data = dataFile.get_var("hHigh");
  data->get(testPhi[1], testPhiSize, 0);

  data = dataFile.get_var("huLow");
  data->get(testPhi[2], testPhiSize, 0);

  data = dataFile.get_var("huHigh");
  data->get(testPhi[3], testPhiSize, 0);

  data = dataFile.get_var("hStar");
  data->get(testPhi[4], testPhiSize, 0);
}

/**
 * Generates a string which contains the variables defined in the solver.
 *
 * @return string containing solver variables.
 */
std::string ComponentsTest::printSolverVariables() {
  std::stringstream errMsg;
  errMsg << std::endl
         << "  dryTol: " << mySolver.dryTol << std::endl
         << "  g: " << mySolver.g << std::endl
         << "  zeroTol: " << mySolver.zeroTol << std::endl
         << "  maxNumberOfNewtonIterations: " << mySolver.maxNumberOfNewtonIterations << std::endl
         << "  newtonTol: " << mySolver.newtonTol << std::endl
         << "  middleStateSpeeds[0]: " << mySolver.middleStateSpeeds[0] << std::endl
         << "  middleStateSpeeds[1]: " << mySolver.middleStateSpeeds[1] << std::endl
         << "  hLeft: " << mySolver.hLeft << ", hRight: " << mySolver.hRight << std::endl
         << "  uLeft: " << mySolver.uLeft << ", uRight: " << mySolver.uRight << std::endl
         << "  huLeft: " << mySolver.huLeft << ", huRight: " << mySolver.huRight << std::endl
         << "  bLeft: " << mySolver.bLeft << ", bLeft: " << mySolver.bRight << std::endl
         << "  hMiddle: " << mySolver.hMiddle << std::endl
         << "  wetDryState: " << mySolver.wetDryState << std::endl
         << "  riemannStructure: " << mySolver.riemannStructure << std::endl
         << "\n  *** WetDryState" << std::endl
         << "  DryDry: " << mySolver.DryDry
         << ", WetWet: " << mySolver.WetWet
         << ", WetDryInundation: " << mySolver.WetDryInundation
         << ", WetDryWall: " << mySolver.WetDryWall
         << ", WetDryWallInundation: " << mySolver.WetDryWallInundation
         << ", DryWetInundation: " << mySolver.DryWetInundation
         << ", DryWetWall: " << mySolver.DryWetWall
         << ", DryWetWallInundation: " << mySolver.DryWetWallInundation << std::endl
         << "\n  *** riemannStructure" << std::endl
         << "  DrySingleRarefaction: " << mySolver.DrySingleRarefaction
         << ", SingleRarefactionDry: " << mySolver.SingleRarefactionDry
         << ", ShockShock: " << mySolver.ShockShock
         << ", ShockRarefaction: " << mySolver.ShockRarefaction
         << ", RarefactionShock: " << mySolver.RarefactionShock
         << ", RarefactionRarefaction: " << mySolver.RarefactionRarefaction << std::endl;
  errMsg.flush();

  return errMsg.str();
}

/**
 * Generates a string, which contains the variables of the solver and the given Riemann solution in the net-update formulation.
 *
 * @param i_hUpdateLeft net-update for the height of the cell on the left side of the edge.
 * @param i_hUpdateRight net-update for the height of the cell on the right side of the edge.
 * @param i_huUpdateLeft net-update for the momentum of the cell on the left side of the edge.
 * @param i_huUpdateRight net-update for the momentum of the cell on the right side of the edge.
 * @param i_maxWaveSpeed maximum wave speed.
 * @return string conataining solver variables and Riemann solution.
 */
std::string ComponentsTest::printSolverVariablesWithResults( const double i_hUpdateLeft,
                                                             const double i_hUpdateRight,
                                                             const double i_huUpdateLeft,
                                                             const double i_huUpdateRight,
                                                             const double i_maxWaveSpeed ) {
  std::stringstream errMsg;
  errMsg << printSolverVariables()
         << "\n  *** return values of computeNetUpdates" << std::endl
         << "  hUpdateLeft: " << i_hUpdateLeft << std::endl
         << "  hUpdateRight: " << i_hUpdateRight << std::endl
         << "  huUpdateLeft: " << i_huUpdateLeft << std::endl
         << "  huUpdateRight: " << i_huUpdateRight << std::endl
         << "  maxWaveSpeed: " << i_maxWaveSpeed << std::endl;
  errMsg.flush();

  return errMsg.str();
}

/**
 * Tests the function determineRiemannStructure(...) of the Augmented Riemann solver.
 */
void ComponentsTest::checkDetermineRiemannStructure() {
  solver::AugRie<double>::RiemannStructure myRiemannStructure;

  for(int i = 0; i < numberOfRandomTests; i++) {

    //test for shock/shock problems
    hLeft = hRight = createRandomNumber(zeroTol, maxWaterHeight);


    uLeft = createRandomNumber(14166.6, 0);
    uRight = createRandomNumber(0, -maxWaveSpeed);

    myRiemannStructure = mySolver.determineRiemannStructure(hLeft, hRight, uLeft, uRight);
    ASSERT(myRiemannStructure == solver::AugRie<double>::ShockShock);

    //test for rare/rare problems
    uLeft = createRandomNumber(-14166.6, 0);
    uRight = createRandomNumber(0, +maxWaveSpeed);

    myRiemannStructure = mySolver.determineRiemannStructure(hLeft, hRight, uLeft, uRight);
    ASSERT(myRiemannStructure == solver::AugRie<double>::RarefactionRarefaction);

    //test for shock/rare problems
    hLeft = createRandomNumber(zeroTol, 8000);
    uLeft = uRight = 0.;
    myRiemannStructure = mySolver.determineRiemannStructure(hLeft, hRight, uLeft, uRight);
    if(hLeft < hRight) {
      ASSERT(myRiemannStructure == solver::AugRie<double>::ShockRarefaction);
    }
    else if(hLeft > hRight) {
      ASSERT(myRiemannStructure == solver::AugRie<double>::RarefactionShock);
    }
  }
}

/**
 * Test the function determineMiddleState(...) of the Augmented Riemann solver with pre-computed solutions.
 */
void ComponentsTest::checkComputeMiddleState() {
  //allow lots of Newton Iterations
  solver::AugRie<double> mySolver(dryTol, gravity, newtonTol, maxNumberOfNewtonIterationsExpensive, zeroTol);

  //do the tests
  for(int i = 0; i < testPhiSize; i++) {
    double hLeft = testPhi[0][i];
    double hRight = testPhi[1][i];

    double huLeft = testPhi[2][i];
    double huRight = testPhi[3][i];

    double uLeft = huLeft / hLeft;
    double uRight = huRight / hRight;

    mySolver.computeMiddleState( hLeft, hRight,
                                 uLeft, uRight,
                                 huLeft, huRight,
                                 maxNumberOfNewtonIterationsExpensive
                               );
    if(testPhi[4][i] < 25.)
      ASSERT_EQUAL_DELTA(mySolver.hMiddle, testPhi[4][i], .1);
    else if(testPhi[4][i] < 100.)
      ASSERT_EQUAL_DELTA(mySolver.hMiddle, testPhi[4][i], 1.);
    else
      ASSERT_EQUAL_DELTA(mySolver.hMiddle, testPhi[4][i], 2.);
  }

}

/**
 * Test the function determinWetDryState(...) of the Augmented Riemann solver in a simple way
 */
void ComponentsTest::checkDetermineWetDryStateSimple() {
  for(int i = 0; i < numberOfRandomTests; i++) {
    mySolver.bLeft = 0.;
    mySolver.bRight = 0.;

    mySolver.uLeft = createRandomNumber(1., maxWaveSpeed);
    mySolver.uRight = createRandomNumber(-maxWaveSpeed, -1);

    //test case dry/dry
    mySolver.hLeft = createRandomNumber(0., dryTol-zeroTol);
    mySolver.hRight = createRandomNumber(0.,dryTol-zeroTol);

    mySolver.determineWetDryState();
    ASSERTM(printSolverVariables() ,mySolver.wetDryState == solver::AugRie<double>::DryDry);

    //test case wet/wet
    mySolver.hLeft = createRandomNumber(dryTol+zeroTol, maxWaterHeight);
    mySolver.hRight = createRandomNumber(dryTol+zeroTol, maxWaterHeight);

    mySolver.determineWetDryState();
    ASSERTM(printSolverVariables(), mySolver.wetDryState == solver::AugRie<double>::WetWet);

    //test case WetDryInundation
    mySolver.bLeft = createRandomNumber(-maxWaterHeight, maxWaterHeight);
    mySolver.bRight = mySolver.bLeft+createRandomNumber(0., 100.);

    mySolver.hLeft = mySolver.bRight - mySolver.bLeft + createRandomNumber(dryTol+zeroTol, maxWaterHeight);
    mySolver.hRight = createRandomNumber(0.,dryTol-zeroTol);

    mySolver.determineWetDryState();
    ASSERTM(printSolverVariables(), mySolver.wetDryState == solver::AugRie<double>::WetDryInundation);

    //test case DryWetInundation
    mySolver.bRight = createRandomNumber(-maxWaterHeight, maxWaterHeight);
    mySolver.bLeft = mySolver.bRight+createRandomNumber(0., 100.);

    mySolver.hLeft = createRandomNumber(0.,dryTol-zeroTol);
    mySolver.hRight =  mySolver.bLeft - mySolver.bRight + createRandomNumber(dryTol+zeroTol, maxWaterHeight);

    mySolver.determineWetDryState();
    ASSERTM(printSolverVariables(), mySolver.wetDryState == solver::AugRie<double>::DryWetInundation);
  }
}

/**
 * Test the function determineWetDrystate(...) of the Augmented Riemann solver for the right inundation limits in one direction.
 *
 * @param i_leftCellDry true, if the left cell is dry. else: false
 */
void ComponentsTest::checkInundationLimitsInOneDirection(const bool i_leftCellDry) {
  //do the tests
  for(int i = 0; i < testPhiWallSamples; i++) {
    double hLeft = testPhi[0][i];
    double hRight = testPhi[1][i];

    double huLeft = testPhi[2][i];
    double huRight = testPhi[3][i];

    double hMiddle = testPhi[4][i];


    double bLeft = 0.;
    double bRight = 0.;
    if(i_leftCellDry == true) {
      bLeft = createRandomNumber(-50, 50);
      bRight = bLeft - hRight - createRandomNumber(0, 10);
    }
    else {
      bRight = createRandomNumber(-50,50);
      bLeft = bRight - hLeft - createRandomNumber(0, 10);
    }

    solver::AugRie<double>::WetDryState preComputedWetDryState;

    if(i_leftCellDry == false) {
     if(bLeft + hMiddle < bRight)
       preComputedWetDryState = solver::AugRie<double>::WetDryWall;
     else
       preComputedWetDryState = solver::AugRie<double>::WetDryWallInundation;
    }

    if(i_leftCellDry == true) {
      if(bRight + hMiddle < bLeft)
        preComputedWetDryState = solver::AugRie<double>::DryWetWall;
      else
        preComputedWetDryState = solver::AugRie<double>::DryWetWallInundation;
    }

    if(i_leftCellDry == true) {
      mySolver.hLeft = 0.;
      mySolver.huLeft = 0.;

      mySolver.hRight = hRight;
      mySolver.huRight = huRight;
    }
    else {
      mySolver.hLeft = hLeft;
      mySolver.huLeft = huLeft;

      mySolver.hRight = 0.;
      mySolver.huRight = 0.;
    }

    mySolver.bLeft = bLeft;
    mySolver.bRight = bRight;

    mySolver.determineWetDryState();
    ASSERT(preComputedWetDryState == mySolver.wetDryState);
  }
}

/**
 * Test the function determineWetDrystate(...) of the Augmented Riemann solver for inundation with much more effort.
 */
void ComponentsTest::checkInundationLimits() {
  checkInundationLimitsInOneDirection(true);
  checkInundationLimitsInOneDirection(false);
}

/**
 * Check if the function computeNetUpdate(...) of the Augmented Riemann solver computes updates for walls
 * in one direction. And Check the positivity of the update.
 * @param i_leftCellDry true, if the left cell is dry. else: false.
 */
void ComponentsTest::checkComputeNetUpdatesForWallsInOneDirection(bool i_leftCellDry) {
  //do the tests
  for(int i = 0; i < testPhiWallSamples; i++) {

    double hLeft, hRight, huLeft, huRight, bLeft, bRight;
    hLeft = hRight = huLeft = huRight = bLeft = bRight = 0.;

    if(i_leftCellDry == true) {
      hRight = testPhi[1][i];
      huRight = testPhi[3][i];
    }
    else {
      hLeft = testPhi[0][i];
      huLeft = testPhi[2][i];
    }

    if(i_leftCellDry == true) {
      bLeft = createRandomNumber(-100, 100);
      bRight = bLeft - hRight - createRandomNumber(0, 20);
    }
    else {
      bRight = createRandomNumber(-50,50);
      bLeft = bRight - hLeft - createRandomNumber(0, 10);
    }

    double hUpdateLeft, hUpdateRight, huUpdateLeft, huUpdateRight, maxWaveSpeed = 0.;

    mySolver.computeNetUpdates( hLeft, hRight,
                                huLeft, huRight,
                                bLeft, bRight,
                                hUpdateLeft, hUpdateRight,
                                huUpdateLeft, huUpdateRight,
                                maxWaveSpeed);

    if(mySolver.wetDryState == solver::AugRie<double>::WetDryWall) {
      ASSERT_EQUAL_DELTA(hUpdateRight, 0., zeroTol);
      ASSERT_EQUAL_DELTA(huUpdateRight, 0., zeroTol);
    }
    else if(mySolver.wetDryState == solver::AugRie<double>::DryWetWall) {
      ASSERT_EQUAL_DELTA(hUpdateLeft, 0., zeroTol);
      ASSERT_EQUAL_DELTA(huUpdateLeft, 0., zeroTol);
    }

    //check if we are positivity preserving
    //Mini-Update with CFL
    ASSERTM( printSolverVariablesWithResults( hUpdateLeft, hUpdateRight,
                                              huUpdateLeft, huUpdateRight,
                                              maxWaveSpeed ),
             hLeft - (.4 / maxWaveSpeed)*hUpdateLeft > -zeroTol );

    //Mini-Update with CFL
    ASSERTM( printSolverVariablesWithResults( hUpdateLeft, hUpdateRight,
                                              huUpdateLeft, huUpdateRight,
                                              maxWaveSpeed ),
             hRight - (.4 / maxWaveSpeed)*hUpdateRight > -zeroTol);
  }
}

/**
 * Check if the function computeNetUpdate(...) of the
 * augmented Riemann solver computes the right updates for walls.
 */
void ComponentsTest::checkComputeNetUpdatesForWalls() {
  checkComputeNetUpdatesForWallsInOneDirection(true);
  checkComputeNetUpdatesForWallsInOneDirection(false);
}

void ComponentsTest::checkComputeNetUpdatesForPositivity() {
  for(int i = 0; i < numberOfRandomTests; i++) {
    double hUpdateLeft, hUpdateRight, huUpdateLeft, huUpdateRight, maxWaveSpeed = 0.;
    double hLeft = createRandomNumber(0, 200);
    double hRight = hLeft - createRandomNumber(-20, 20);
    hRight = std::max(0., hRight);

    double uLeft = createRandomNumber(-maxWaveSpeed, maxWaveSpeed);
    double uRight = maxWaveSpeed + createRandomNumber(-maxWaveSpeedVariation, maxWaveSpeedVariation);

    double huLeft = hLeft * uLeft;
    double huRight = hRight * uRight;

    double bLeft = createRandomNumber(-100, 100);
    double bRight = bLeft + createRandomNumber(-2, 2)*(hLeft-hRight);

    mySolver.computeNetUpdates( hLeft, hRight,
                                huLeft, huRight,
                                bLeft, bRight,
                                hUpdateLeft, hUpdateRight,
                                huUpdateLeft, huUpdateRight,
                                maxWaveSpeed);

    if(maxWaveSpeed < zeroTol)
      continue;

    ASSERTM( printSolverVariablesWithResults( hUpdateLeft, hUpdateRight,
                                              huUpdateLeft, huUpdateRight,
                                              maxWaveSpeed ),
             hLeft - (.45 / maxWaveSpeed)*hUpdateLeft > -zeroTol);
    ASSERTM( printSolverVariablesWithResults( hUpdateLeft, hUpdateRight,
                                              huUpdateLeft, huUpdateRight,
                                              maxWaveSpeed ),
             hRight - (.45 / maxWaveSpeed)*hUpdateRight > -zeroTol);
  }
}





