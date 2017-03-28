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
// $Id: B1RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1Analysis.hh"
#include "G4RunManager.hh"


#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(B1EventAction* eventAction)
 : G4UserRunAction(),
   fEventAction(eventAction)
{ 
  G4RunManager::GetRunManager()->SetPrintProgress(1); 
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B5Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("B1");

  // Book histograms, ntuple
  //
  
  // Creating 1D histograms
  analysisManager
    ->CreateH1("Time","Time", 900000, 0., 900); // h1 Id = 0
  //analysisManager
  //  ->CreateH1("Chamber2","Drift Chamber 2 # Hits", 50, 0., 50); // h1 Id = 1
  
  // Creating 2D histograms
   analysisManager                                                
    ->CreateH2("Cell XY","Detector cells X vs Y",           // h2 Id = 0
               4, -20., 20, 4, -20., 20.); 


  //analysisManager                                                
  //  ->CreateH2("Chamber2 XY","Drift Chamber 2 X vs Y",           // h2 Id = 1
  //             50, -1500., 1500, 50, -300., 300.);

  // Creating ntuple
  //
  if ( fEventAction ) {
    analysisManager->CreateNtuple("B1", "Hits");
    analysisManager->CreateNtupleDColumn("Time");
    //analysisManager->CreateNtupleIColumn("DetectorTriggerMinusOne");
   // analysisManager->CreateNtupleIColumn("Dc1Hits");  // column Id = 0
   // analysisManager->CreateNtupleIColumn("Dc2Hits");  // column Id = 1
    analysisManager->CreateNtupleDColumn("ECEnergy"); // column Id = 2
   // analysisManager->CreateNtupleDColumn("HCEnergy"); // column Id = 3
   // analysisManager->CreateNtupleDColumn("Time1");    // column Id = 4
   // analysisManager->CreateNtupleDColumn("Time2");    // column Id = 5
    analysisManager                                   // column Id = 6
      ->CreateNtupleDColumn("ECEnergyVector", fEventAction->GetEmCalEdep());
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Egap"); 
    analysisManager->CreateNtupleDColumn("Etotal");  
   // analysisManager                                   // column Id = 7
   //   ->CreateNtupleDColumn("HCEnergyVector", fEventAction->GetHadCalEdep());
    analysisManager->FinishNtuple();
  }

 // Creating 3D histograms
   analysisManager                                                
    ->CreateH3("Cell XY and eDep","Detector cells X vs Y vs eDep",           // h3 Id = 0
               4, -20., 20, 4, -20., 20., 5, 0, 50); 

  analysisManager                                                
    ->CreateH2("Cell ID and eDep","Detector cells ID vs eDep",           // h2 Id = 0
               8, 0., 64, 5, 0., 50.); 


 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file 
  // The default file name is set in B5RunAction::B5RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
