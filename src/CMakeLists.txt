set(SIM_SOURCES main.cxx DetectorConstruction.cxx GasGapSensitiveDetector.cxx GasGapHit.cxx GasGapRunAction.cxx GasGapEventAction.cxx GasGapEventActionMessenger.cxx GasGapSteepingAction.cxx GasGapTrajectory.cxx GasGapTrackInformation.cxx PrimaryGenerator.cxx PrimaryGeneratorMessenger.cc PhysicsList.cc PhysicsListMessenger.cc StepMax.cc StepMaxMessenger.cc TrGEMAnalysis.cxx DetectorMessenger.cc GasGapElectricFieldSetup.cc GasGapFieldMessenger.cc SteppingAction.cc)


add_executable(sim ${SIM_SOURCES})
target_link_libraries(sim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
install(TARGETS sim DESTINATION ${CMAKE_INSTALL_PREFIX})
