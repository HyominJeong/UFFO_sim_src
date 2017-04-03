//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include <G4SDKineticEnergyFilter.hh>
#include <G4SDParticleFilter.hh>
#include <G4VSDFilter.hh>

#include "B1DetectorConstruction.hh"
#include "B1ShieldSD.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4Trap.hh>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <CLHEP/Vector/Rotation.h>
#include "G4RotationMatrix.hh"
#include <G4PVReplica.hh>
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include <G4PVParameterised.hh>
#include "B1CellParameterisation.hh"
#include "B1EmCalorimeterSD.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpticalSurface.hh" // for build mirror


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 200*cm, env_sizeZ = 300*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
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
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
/*
  //     
  // Aluminum Hopper structure
  //

  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Al");

  // Trapezoid shape       
  G4double trapAl_dxa_outer = 17.648*cm, trapAl_dxb_outer = 40.086*cm;
  G4double trapAl_dya_outer = 17.648*cm, trapAl_dyb_outer = 40.086*cm;
  
  G4double trapAl_dxa_inner = 17.002*cm, trapAl_dxb_inner = 39.44*cm;
  G4double trapAl_dya_inner = 17.002*cm, trapAl_dyb_inner = 39.44*cm;

  G4double  pX_Al_cor = 170.02*0.5*mm, pY_Al_cor = 170.02*0.5*mm, pZ_Al_cor = 4*0.5*mm;

  G4double trapAl_dz  = 28*cm;      

  G4Trd* solidOuterTrapAl =    
    new G4Trd("Outer trapezoid Al",                      //its name
              0.5*trapAl_dxa_outer, 0.5*trapAl_dxb_outer, 
              0.5*trapAl_dya_outer, 0.5*trapAl_dyb_outer, 0.5*trapAl_dz); //its size

  G4Trd* solidInnerTrapAl =    
    new G4Trd("TrapAlShapeInner",                      //its name
              0.5*trapAl_dxa_inner, 0.5*trapAl_dxb_inner, 
              0.5*trapAl_dya_inner, 0.5*trapAl_dyb_inner, 0.5*trapAl_dz); //its size

  G4Box* solidCorrection_Al = 
    new G4Box("correction_Al" ,pX_Al_cor, pY_Al_cor, pZ_Al_cor);

    // This box is purely for cosmetic purposes. Without this box, the trapezoid would appear
    // as having a bottom plate in the visualization, although for physics purposes the plate is non-existent.

  G4SubtractionSolid *TrapAl_noCor = 
    new G4SubtractionSolid("Hollow Trap Al no correction",solidOuterTrapAl,solidInnerTrapAl);

  G4SubtractionSolid *solidTrapAl = 
    new G4SubtractionSolid("Hollow Trap Al",TrapAl_noCor,solidCorrection_Al, 0, G4ThreeVector(0, 0, -14*cm));

   //
   // Lower box Al
   //

  G4double  pX_Al = 199.6*0.5*mm, pY_Al = 199.6*0.5*mm, pZ_Al = 89.6*0.5*mm;
  G4double  pX_Air_big_Al = 193.6*0.5*mm, pY_Air_big_Al = 193.6*0.5*mm, pZ_Air_big_Al = 83.6*0.5*mm;
  G4double  pX_Air_small_Al = 170.02*0.5*mm, pY_Air_small_Al = 170.02*0.5*mm, pZ_Air_small_Al = 4*0.5*mm;
  G4ThreeVector box_pos = G4ThreeVector(0, 0*cm, -(140+45)*mm);
  G4ThreeVector box_small_shift_Al = G4ThreeVector(0, 0*mm, 43.3*mm);

    G4Box* solidLowBox_Al = 
    new G4Box("lowBox_Al" ,pX_Al, pY_Al, pZ_Al);

    G4Box* solidLowBox_Air_big_Al = 
    new G4Box("lowBox_Air_big" ,pX_Air_big_Al, pY_Air_big_Al, pZ_Air_big_Al);

  G4Box* solidLowBox_Air_small_Al = 
    new G4Box("lowBox_Air_small" ,pX_Air_small_Al, pY_Air_small_Al, pZ_Air_small_Al);

  G4UnionSolid* union_Air_Al =
    new G4UnionSolid("Big_air+Small_air", solidLowBox_Air_big_Al, solidLowBox_Air_small_Al,0, box_small_shift_Al); 

  G4SubtractionSolid *solidHollowBox_Al = 
    new G4SubtractionSolid("lowBox",solidLowBox_Al,union_Air_Al,0, G4ThreeVector(0, 0*cm, 0*mm));

    //
    // UnionSolid of lower box and cone of aluminium
    

    G4ThreeVector lowerBoxShift_Al = G4ThreeVector(0, 0*mm, -(140+44.8)*mm);

     G4UnionSolid* solidTrapBoxAl =
    new G4UnionSolid("Union of trapezoid and box Al", solidTrapAl, solidHollowBox_Al,0, lowerBoxShift_Al); 

                
  G4LogicalVolume* logicTrapBoxAl =                         
    new G4LogicalVolume(solidTrapBoxAl,         //its solid
                        shape2_mat,          //its material
                        "GapLV");           //its name
               
    new G4PVPlacement(0,                       //no rotation
                   G4ThreeVector(),                    //at position
                    logicTrapBoxAl,             //its logical volume
                    "Union of trapezoid and box Al",                //its name
                 logicEnv,                //its mother  volume
                  false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

 // G4Box* solidBox_big_Al = 
 //   new G4Box("Al box big" ,10*cm, 10*cm, 10*cm);

 // G4Box* solidBox_small_Al = 
  //  new G4Box("Al box small" ,100-6*mm, 100-6*mm, 100-6*mm);

 // G4SubtractionSolid *solidHollowBox_Al = 
 //   new G4SubtractionSolid("Al box solid",solidBox_big_Al,solidBox_small_Al,0, G4ThreeVector());

//  G4LogicalVolume* logicBoxAl =                         
 //   new G4LogicalVolume(solidHollowBox_Al,         //its solid
//                        shape2_mat,          //its material
//                        "GapLV");           //its name
               
//    new G4PVPlacement(0,                       //no rotation
//                    G4ThreeVector(),                    //at position
//                    logicBoxAl,             //its logical volume
//                    "Al box physical",                //its name
//                    logicEnv,                //its mother  volume
//                    false,                   //no boolean operation
  //                  0,                       //copy number
    //                checkOverlaps);          //overlaps checking




    // Tungsten HOPPER Structure
    // Trapezoid shape
    //

  G4Material* side_mat = nist->FindOrBuildMaterial("G4_W");

  G4double trapW_dxa_outer = 17.6912*cm, trapW_dxb_outer = 40.1292*cm;
  G4double trapW_dya_outer = 17.6912*cm, trapW_dyb_outer = 40.1292*cm;
  
  G4double trapW_dxa_inner = 17.6482*cm, trapW_dxb_inner = 40.0862*cm;
  G4double trapW_dya_inner = 17.6482*cm, trapW_dyb_inner = 40.0862*cm;

  G4double  pX_W_cor = 176.6482*0.5*mm, pY_W_cor = 176.6482*0.5*mm, pZ_W_cor = 4*0.5*mm;

  G4double trapW_dz  = 28*cm;      

  G4Trd* solidOuterTrapW =    
    new G4Trd("Outer trapezoid W",                      //its name
              0.5*trapW_dxa_outer, 0.5*trapW_dxb_outer, 
              0.5*trapW_dya_outer, 0.5*trapW_dyb_outer, 0.5*trapW_dz); //its size

  G4Trd* solidInnerTrapW =    
    new G4Trd("Inner trapezoid W",                      //its name
              0.5*trapW_dxa_inner, 0.5*trapW_dxb_inner, 
              0.5*trapW_dya_inner, 0.5*trapW_dyb_inner, 0.5*trapW_dz); //its size

  G4Box* solidCorrection_W = 
    new G4Box("correction box W" ,pX_W_cor, pY_W_cor, pZ_W_cor);

    // This box is purely for cosmetic purposes. Without this box, the trapezoid would appear
    // as having a bottom plate in the visualization, although for physics purposes the plate is non-existent.

  G4SubtractionSolid *TrapW_noCor = 
    new G4SubtractionSolid("Hollow trapezoid W no correction",solidOuterTrapW,solidInnerTrapW);

  G4SubtractionSolid *solidTrapW = 
    new G4SubtractionSolid("Hollow trapezoid W",TrapW_noCor,solidCorrection_W, 0, G4ThreeVector(0, 0, -14*cm));

   //
   // Lower box W
   //

  G4double  pX_W = 200*0.5*mm, pY_W = 200*0.5*mm, pZ_W = 90*0.5*mm;
  G4double  pX_Air_big_W = 199.602*0.5*mm, pY_Air_big_W = 199.602*0.5*mm, pZ_Air_big_W = 89.602*0.5*mm;
  G4double  pX_Air_small_W = 176.900*0.5*mm, pY_Air_small_W = 176.900*0.5*mm, pZ_Air_small_W = 1*0.5*mm;
  //G4ThreeVector box_pos = G4ThreeVector(0, 0*cm, -(140+45)*mm);
  G4ThreeVector box_small_shift_W = G4ThreeVector(0, 0*mm, 44.9*mm);

    G4Box* solidLowBox_W = 
    new G4Box("Lower box W" ,pX_W, pY_W, pZ_W);

  G4Box* solidLowBox_Air_big_W = 
    new G4Box("lowBox_Air_big" ,pX_Air_big_W, pY_Air_big_W, pZ_Air_big_W);

  G4Box* solidLowBox_Air_small_W = 
    new G4Box("lowBox_Air_small" ,pX_Air_small_W, pY_Air_small_W, pZ_Air_small_W);

  G4UnionSolid* union_Air_W =
    new G4UnionSolid("Big_air+Small_air", solidLowBox_Air_big_W, solidLowBox_Air_small_W,0, box_small_shift_W); 

  G4SubtractionSolid *solidHollowBox_W = 
    new G4SubtractionSolid("lowBox",solidLowBox_W,union_Air_W,0, G4ThreeVector(0, 0*cm, 0*mm));

    //
    // UnionSolid of lower box and cone of tungsten

    //

    G4ThreeVector lowerBoxShift_W = G4ThreeVector(0, 0*mm, -(140+45-0.2)*mm);

     G4UnionSolid* solidTrapBoxW =
    new G4UnionSolid("Union of trapezoid and box W", solidTrapW, solidHollowBox_W,0, lowerBoxShift_W); 


	G4Box* solidMask_W = 
    	new G4Box("W mask" ,394.4*0.5*mm, 394.4*0.5*mm, 1*mm*0.5);

     G4UnionSolid* solidTrapBoxMaskW =
    new G4UnionSolid("Union of trapezoid, andmask  box W", solidTrapBoxW, solidMask_W,0, G4ThreeVector(0, 0,(140 + 0.5)*mm)); 

                
  G4LogicalVolume* logicTrapBoxW =                         
    new G4LogicalVolume(solidTrapBoxMaskW,         //its solid
                        side_mat,          //its material
                        "AbsoLV");           //its name
               
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0, 0),                    //at position
                    logicTrapBoxW,             //its logical volume
                    "Union of trapezoid and box W",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

 //G4Box* solidBox_big_W = 
 //   new G4Box("W box big" ,15*cm, 15*cm, 15*cm);

 // G4Box* solidBox_small_W = 
 //   new G4Box("W box small" ,150-0.4*mm, 150-0.4*mm, 150-0.4*mm);

  //G4SubtractionSolid *solidHollowBox_W = 
  //  new G4SubtractionSolid("W box solid",solidBox_big_W,solidBox_small_W,0, G4ThreeVector());

  //G4LogicalVolume* logicBoxW =                         
  //  new G4LogicalVolume(solidHollowBox_W,         //its solid
   //                     side_mat,          //its material
  //                      "AbsoLV");           //its name
               
  //  new G4PVPlacement(0,                       //no rotation
  //                  G4ThreeVector(),                    //at position
  //                  logicBoxW,             //its logical volume
  //                  "W box physical",                //its name
  //                  logicEnv,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking
*/


    //
    // 
    //
    // 
    //
    // Time to build the actual YSO detector 
    //
    // Mother volume of parameterisation

  G4Box* solidCellMother =    
    new G4Box("solidCellMother",                       		//its name
       0.5*2.88*8*mm, 0.5*2.88*8*mm, 0.5*3*mm);   		//its size
      
  G4LogicalVolume* logicCellMother =                         
    new G4LogicalVolume(solidCellMother,         			//its solid
                        world_mat,           				//its material
                        "logicCellMother");            		//its name
                                   
  G4VPhysicalVolume* physCellMother = 
    new G4PVPlacement(0,                     				//no rotation
                      G4ThreeVector(0.,0.,-141.50*mm),       //at (0,0,0)
                      logicCellMother,            			//its logical volume
                      "physCellMother",               		//its name
                      logicEnv,                  		    //its mother  volume
                      false,       					        //no boolean operation
                      0,                     				//copy number
                      checkOverlaps);        				//overlaps checking


    //
    //
    //

  // YSO material

  G4double a, density, z;
  G4String name, symbol;
  G4int ncomponents, natoms; 
     	
	a = 88.90585*g/mole;
	G4Element* elY  = new G4Element(name="Yttrium",symbol="Y" , z= 39., a);

 	a = 28.085*g/mole;
	G4Element* elSi  = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

	a = 16.00*g/mole;
	G4Element* elO  = new G4Element(name="Oxygen",symbol="O" , z= 8., a);

  // Ce dopping
  a = 140.116*g/mole;
  G4Element* elCe  = new G4Element(name="Cerium",symbol="Ce" , z= 58., a);

 	density = 4.45*g/cm3;
 	G4Material* YSO = new G4Material(name="YSO",density,ncomponents=4);
	YSO->AddElement(elY, natoms=2);
	YSO->AddElement(elSi, natoms=1);
	YSO->AddElement(elO, natoms=5);
  YSO->AddElement(elCe, natoms=0.4);

  // Add Material Property table

  G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  // Assumed refrection index
  G4double refractiveIndex1[] =
            { 1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800,  1.800,  1.800,  1.800,
              1.800,  1.800};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));
  
  
/*  G4double absorption[] =
            {  3.448*m,  4.082*m,  6.329*m,  9.174*m,
              12.346*m, 13.889*m, 15.152*m, 17.241*m,
              18.868*m, 20.000*m, 26.316*m, 35.714*m,
              45.455*m, 47.619*m, 52.632*m, 52.632*m,
              55.556*m, 52.632*m, 52.632*m, 47.619*m,
              45.455*m, 41.667*m, 37.037*m, 33.333*m,
              30.000*m, 28.500*m, 27.000*m, 24.500*m,
              22.000*m, 19.500*m, 17.500*m, 14.500*m };
*/

  G4double absorption[] =
            {  1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm,
               1.16*cm,  1.16*cm,  1.16*cm,  1.16*cm };

  assert(sizeof(absorption) == sizeof(photonEnergy));
  
  /*
  G4double scintilFast[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));
  */
  G4double scintilFast[] = //scintilSlow[] =
            { 0.01, 0.01, 0.10, 0.20,
              0.21, 0.22, 0.23, 0.24,
              0.25, 0.26, 0.27, 0.28,
              0.29, 0.30, 0.35, 0.40,
              0.45, 0.50, 0.60, 0.70,
              0.40, 0.10, 0.01, 0.01,
              0.01, 0.01, 0.01, 0.01,
              0.01, 0.01, 0.01, 0.01 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));
  

  G4MaterialPropertiesTable* YSOPMT = new G4MaterialPropertiesTable();
  
  YSOPMT->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
  
  YSOPMT->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
        ->SetSpline(true);
  YSOPMT->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
        ->SetSpline(true);
  //YSOPMT->AddProperty("SLOWCOMPONENT",photonEnergy, scintilFast,     nEntries)
  //      ->SetSpline(true);
  //YSOPMT->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
  //      ->SetSpline(true);

  YSOPMT->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  //YSOPMT->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  YSOPMT->AddConstProperty("RESOLUTIONSCALE",1.0);
  YSOPMT->AddConstProperty("FASTTIMECONSTANT", 35.*ns);
  //YSOPMT->AddConstProperty("SLOWTIMECONSTANT",35.*ns);
  YSOPMT->AddConstProperty("YIELDRATIO",0.8);
  
  G4cout << "YSO G4MaterialPropertiesTable" << G4endl;
  YSOPMT->DumpTable();

  YSO->SetMaterialPropertiesTable(YSOPMT);

  // Sensitive detector setup

  //G4double detectSize = 169.92*mm;

  G4double CrystalX = (2.88-0.002)*0.5*mm;
  G4double CrystalY = CrystalX;
  G4double CrystalZ = 3*0.5*mm;


  G4Box *CrystalSolid = 
      new G4Box("CrystalSolid", CrystalX, CrystalY, CrystalZ);
 
  G4LogicalVolume* fCrystalLog = 
      new G4LogicalVolume(CrystalSolid,YSO,
                                            "logCrystal");


  G4String tName1("Crystal");        // Allow all target physicals to share
 
  // Surface properties of YSO - mirror

  //G4OpticalSurface* OpYSOSurface = new G4OpticalSurface()

  ////////////////////////////
  // Set regions for SetCut //
  ////////////////////////////
  G4Region* regCrystal = new G4Region("regCrystal");
  logicCellMother -> SetRegion(regCrystal);
  regCrystal -> AddRootLogicalVolume(logicCellMother);


  G4VPVParameterisation* cellParam = new B1CellParameterisation();
      new G4PVParameterised("cellPhysical",fCrystalLog,logicCellMother,
                          kXAxis,64,cellParam);


  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;

  G4VSensitiveDetector* emCalorimeter 
    = new B1EmCalorimeterSD(SDname="/EMcalorimeter");
  G4SDParticleFilter* electronFilter
    = new G4SDParticleFilter("electronFilter","e-");
  G4SDParticleFilter* gammaFilter
    = new G4SDParticleFilter("gammaFilter","gamma");
  G4SDParticleFilter* positronFilter
    = new G4SDParticleFilter("positronFilter","e+");
  // G4SDKineticEnergyFilter* energyFilter = 
  // new G4SDKineticEnergyFilter("energyFilter",5*keV, 200*keV);
  emCalorimeter->SetFilter(positronFilter);
  emCalorimeter->SetFilter(electronFilter);
  emCalorimeter->SetFilter(gammaFilter);
  
  // emCalorimeter->SetFilter(energyFilter);
  SDman->AddNewDetector(emCalorimeter);
  fCrystalLog->SetSensitiveDetector(emCalorimeter); // 


/*
  B1ShieldSD* absoSD 
    = new B1ShieldSD("AbsorberSD", "AbsorberHitsCollection", 1);
  SDman->AddNewDetector(absoSD);
  logicTrapBoxAl->SetSensitiveDetector(absoSD);

  B1ShieldSD* gapSD 
    = new B1ShieldSD("GapSD", "GapHitsCollection", 1);
  SDman->AddNewDetector(gapSD);
  logicTrapBoxW->SetSensitiveDetector(gapSD);
*/
  /////////////////
  // Set Scorers //
  /////////////////
  /*
  G4MultiFunctionalDetector* myScorer = new G4MultiFunctionalDetector("myCellScorer");
  G4VPrimitiveSensitivity* totalSurfFlux = new G4PSFlatSurfaceFlux("TotalSurfFlux");
  myScorer->Register(totalSurfFlux)
  G4VPrimitiveSensitivity* protonSurfFlux = new G4PSFlatSurfaceFlux("protonSurfFlux");
  protonFilter->Add("proton");
  protonSurfFlux->SetFilter(protonFilter);
  myScorer->Resiger(protonSurfFlux);

  SetSensitiveDetector("myLogVol", myScorer);
  */
  /////////////////
  // Set Scorers //
  /////////////////
  // 
  //

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
