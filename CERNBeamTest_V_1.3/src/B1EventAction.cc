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
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1EmCalorimeterHit.hh"
#include "B1Analysis.hh"

#include "B1ShieldSD.hh"
#include "B1ShieldHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction()
: G4UserEventAction(), 
  fECHCID(-1),
  fEmCalEdep(),
  HitVec(),
  TimeVec1(0),
  TimeVec2(0),
  fAbsHCID(-1),
  fGapHCID(-1)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // initialize the vectors
  fEmCalEdep.resize(64, 0.);
  HitVec.resize(64, 0.);
  TimeVec1.resize(64, 0.);
  TimeVec2.resize(64, 0.);
      for (G4int i=0;i<64;i++)
	{
		TimeVec2[i] = 0;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}


B1ShieldHitsCollection* 
B1EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  B1ShieldHitsCollection* hitsCollection 
    = static_cast<B1ShieldHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B1EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{
    if (fECHCID==-1) {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      fECHCID = sdManager->GetCollectionID("EMcalorimeter/EMcalorimeterColl");
    }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) 
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found." << G4endl; 
        G4Exception("B1EventAction::EndOfEventAction()",
                    "B1Code001", JustWarning, msg);
        return;
    }   


    // Get hits collections 

    B1EmCalorimeterHitsCollection* ecHC 
      = static_cast<B1EmCalorimeterHitsCollection*>(hce->GetHC(fECHCID));
      
      
    if (  (!ecHC)) 
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found." << G4endl; 
        G4Exception("B1EventAction::EndOfEventAction()",
                    "B1Code001", JustWarning, msg);
        return;
    }   
    
    //
    // Fill histograms & ntuple
    // 
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
    // Fill histograms


 


//	for (G4int i=0;i<j;i++)
//    {
//       B1EmCalorimeterHit* hit = (*ecHC)[HitVec[i]];
//       G4ThreeVector localPos = hit->GetPos();
//       analysisManager->FillH2(0, localPos.x(), localPos.y());
//    }
 

    // Fill ntuple
    
    // Dc1Hits
    //analysisManager->FillNtupleIColumn(0, ecHC->entries());
    // Dc2Hits
    //analysisManager->FillNtupleIColumn(1, ecHC->entries());
    
    // ECEnergy
    G4int totalEmHit = 0;
    G4double totalEmE = 0.;
    G4int j = 0;
    for (G4int i=0;i<64;i++)
    {
        B1EmCalorimeterHit* hit = (*ecHC)[i];
        G4double eDep = hit->GetEdep()/keV;
        G4double time = hit->GetGTime()/s;
        analysisManager->FillNtupleDColumn(0, time);
	TimeVec1[i] = time;
        if (eDep>0.)
        {
            totalEmHit++;
            totalEmE += eDep;
            
	    j++;
            if ( (TimeVec2[i]==0) || ((TimeVec1[i] - TimeVec2[i]) > 0.000005) )
            { 
              analysisManager->FillH1(0, time);
	      TimeVec2[i] = time;
            }
        }
        fEmCalEdep[i] = eDep;
    }

    analysisManager->FillNtupleDColumn(1, totalEmE);

	for (G4int i=0;i<j;i++)
    {
       B1EmCalorimeterHit* hit = (*ecHC)[HitVec[i]];
       G4ThreeVector localPos = hit->GetPos()/cm;
       analysisManager->FillH3(0, localPos.x(), localPos.y(), fEmCalEdep[HitVec[i]]);

       G4int CellIDNo = hit->GetCellID();
       analysisManager->FillH2(1, CellIDNo, fEmCalEdep[HitVec[i]]);
    }

// Get hits collections IDs (only once)
/*
  if ( fAbsHCID == -1 ) {
    fAbsHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }
*/
  // Get hits collections
  //B1ShieldHitsCollection* absoHC = GetHitsCollection(fAbsHCID, event);
  //B1ShieldHitsCollection* gapHC = GetHitsCollection(fGapHCID, event);

  // Get hit with total values
  //B1ShieldHit* absoHit = (*absoHC)[absoHC->entries()-1];
  //B1ShieldHit* gapHit = (*gapHC)[gapHC->entries()-1];
 
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    /*
    PrintEventStatistics(
      absoHit->GetEdep(), absoHit->GetTrackLength(),
      gapHit->GetEdep(), gapHit->GetTrackLength());
    */
  }  

  //analysisManager->FillNtupleDColumn(3, absoHit->GetEdep()/keV);
  //analysisManager->FillNtupleDColumn(4, gapHit->GetEdep()/keV);

 G4double totalEdep = 0;
 //totalEdep = totalEmE + absoHit->GetEdep()/keV + gapHit->GetEdep()/keV;

  analysisManager->FillNtupleDColumn(5, totalEdep);

    
    analysisManager->AddNtupleRow();  
    
    //
    // Print diagnostics
    // 
    

    if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
    
    G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout << G4endl
           << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
           << primary->GetG4code()->GetParticleName()
           << " " << primary->GetMomentum() << G4endl;
    
    

    // EM calorimeter
    G4cout << "EM Calorimeter has " << totalEmHit << " hits. Total Edep is "
    << totalEmE/keV << " (keV)" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

