/**
 * UnitTests.cpp
 *
 ****
 **** Collection of unit tests for the solvers of the shallow water equations.
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

#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"

#include "ComponentsTest.h"

#include "../solver/AugRieGeoClaw.hpp"

#include <iostream>

/**
 * Runs the test suite of unit tests.
 */
void runSuite(){
	cute::suite testSuite;

	testSuite.push_back(ComponentsTest());
	cute::ide_listener lis;

	cute::makeRunner(lis)(testSuite, "running the test suite");
}

/**
 * Main function.
 *
 * @return 0 if completed, 1 else.
 */
int main(){
  runSuite();



  /**
   * TODO: Benchmarks will be available with the new 1D-code again (end of 2012?)
   */
//  benchmarks::SingleWaveOnASimpleBeach( benchmarks::OneDimensional::AugRie,
//                                        2000,
//                                        480,
//                                        "/work/breuera/workspace/tsunami-benchmarks/analytical/single_wave_on_simple_beach/output/SingleWaveOnASimpleBeachAugRie.nc"
//                                      );
//
//  benchmarks::SingleWaveOnASimpleBeach( benchmarks::OneDimensional::Hybrid,
//                                        2000,
//                                        480,
//                                        "/work/breuera/workspace/tsunami-benchmarks/analytical/single_wave_on_simple_beach/output/SingleWaveOnASimpleBeachHybrid.nc"
//                                       );
//  benchmarks::DamBreakDry damBreakDry("/work/breuera/workspace/tsunami-benchmarks/analytical/single_wave_on_simple_beach/output/DamBreakDry");

  std::cout << "\nUnit tests finished" << std::endl;

  return 0;
}



