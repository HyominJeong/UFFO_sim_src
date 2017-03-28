#include "OTDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

OTDetectorConstruction::OTDetectorConstruction()
: G4VUserDetectorConstruction()
{
}

OTDetectorConstruction::~OTDetectorConstruction()
{
}

G4VPhysicalVolume* OTDetectorConstruction::Construct()
{  
  G4NistManager* nist = G4NistManager::Instance();


  // -----------------------------------------------------
  // World

  // G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");
  // define material, Vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* Vacuum =
    new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);

G4Material* world_mat = Vacuum;


  G4double world_size = 1200*mm;

  G4Box* solidWorld =    
    new G4Box("World",                       // its name
              0.5*world_size,                // half x
              0.5*world_size,                // half y
              0.5*world_size);               // half z
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      true);                 //overlaps checking


  // -----------------------------------------------------
  // Detector

//  G4Material* detector_mat = nist -> FindOrBuildMaterial("G4_Pb");
  // define material, Detector
  G4double detectorAtomicNumber = 1.;
  G4double detectorMassOfMole = 1.008*g/mole;
  G4double detectorDensity = 1.e-20*g/cm3;
  G4double detectorTemperature = 2.73*kelvin;
  G4double detectorPressure = 3.e-18*pascal;
  G4Material* detector_mat =
    new G4Material("interGalactic", detectorAtomicNumber, detectorMassOfMole, detectorDensity, kStateSolid, detectorTemperature, detectorPressure);

  
  G4double detector_size = 20*mm;
  G4double detector_offset_z = 520*mm;

  G4Box* solidDetector =    
    new G4Box("Detector",
              2.5*detector_size,
              2.5*detector_size,
              0.01*detector_size);
      
  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,
                        detector_mat,
                        "Detector");
                                   
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,detector_offset_z),
                      logicDetector,
                      "Detector",
                      logicWorld,
                      false,
                      1,
                      true);

  // Lead block

  G4Material* Lead_block_mat = nist -> FindOrBuildMaterial("G4_Pb");
  G4double Lead_block_size = 10.*mm;
  G4double Lead_block_offset_z = 110.*mm;

  G4Box* solidLead_block =    
    new G4Box("Lead_block",
              20.0*Lead_block_size,
              20.0*Lead_block_size,
              10.0*Lead_block_size);
      
  G4LogicalVolume* logicLead_block =                         
    new G4LogicalVolume(solidLead_block,
                        Lead_block_mat,
                        "Lead_block");
                                   
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,Lead_block_offset_z),
                      logicLead_block,
                      "Lead_block",
                      logicWorld,
                      false,
                      1,
                      true);

  return physWorld;
}



