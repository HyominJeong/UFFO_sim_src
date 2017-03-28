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
// $Id: B1TrackingAction.cc 66587 2012-12-21 11:06:44Z ihrivnac $
//
/// \file src/B1TrackingAction.cc
/// \brief Implementation of the B1TrackingAction class
//


#include "B1TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"


void B1TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	// print out lots for non-primary particle
	// modified >1 to >=0 for check input particle
	if(aTrack->GetParentID()>=1)  //not the input particle (assumed to be just one particle)
	{
	/*
	    G4cout << "B1TrackingAction: " << G4endl;
	    G4cout << " Track ID:          " << aTrack->GetTrackID() << G4endl;
	    G4cout << " particle:          " << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() << G4endl;
	    G4cout << " Parent ID:         " << aTrack->GetParentID() << G4endl;
	    G4cout << " created by:        " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
	    G4cout << " kin. energy (keV): " << aTrack->GetKineticEnergy() / keV << G4endl;
	    G4cout << " volume:            " << aTrack->GetVolume()->GetName() << G4endl;
	    G4cout << " global time:       " << aTrack->GetGlobalTime() << G4endl;
	*/
	}
}
