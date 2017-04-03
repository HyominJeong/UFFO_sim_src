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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4PhysicalConstants.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0),
  i(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="proton");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(1000.*MeV);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //


	// GS: Generate random position on the 4PIsphere to create a unif. distrib.
	// GS: on the sphere

//      G4ThreeVector pos0;
//      G4ThreeVector dir0;
//      G4ThreeVector vertex0; // = G4ThreeVector(x0,y0,z0);
      
//      dir0 = G4ThreeVector(0.,0.,-1.);
      
//      G4double theta, phi, y, f;

//	phi = G4UniformRand() * twopi;
//	do {
//	  y = G4UniformRand()*1.0;
//	  theta = G4UniformRand() * pi;//
//	  f = std::sin(theta);
//	} while (y > f);
//	vertex0 = G4ThreeVector(1.,0.,0.);
//	vertex0.setMag(1*m);
//	vertex0.setTheta(theta);
//	vertex0.setPhi(phi);
//	fParticleGun->SetParticlePosition(vertex0);
//	
//	dir0 = G4ThreeVector(1.,0.,0.);
//	do {
//	  phi = G4UniformRand() * twopi;
//	  do {
//	    y = G4UniformRand()*1.0;
//	    theta = G4UniformRand() * pi;
//	    f = std::sin(theta);
//	  } while (y > f);
//	  dir0.setPhi(phi);
//	  dir0.setTheta(theta);
//	} while (vertex0.dot(dir0) >= -0.7 * vertex0.mag());
//	fParticleGun->SetParticleMomentumDirection((G4ParticleMomentum)dir0);


  //fParticleGun->SetParticlePosition(G4ThreeVector(1.77*mm,1.77*mm,10*cm));

  G4double alphaMin =  0*deg;      //alpha in [0,pi]
  G4double alphaMax = 180*deg;
  G4double fCosAlphaMin = std::cos(alphaMin);
  G4double fCosAlphaMax = std::cos(alphaMax);
  
 G4double fPsiMin = 0*deg;       //psi in [0, 2*pi]
 G4double fPsiMax = 360*deg;

 //G4double x0 = 0*mm, y0 = 0*mm, z0 = 0*mm;
  //direction uniform in solid angle
  //
  G4double cosAlpha = fCosAlphaMin-G4UniformRand()*(fCosAlphaMin-fCosAlphaMax);
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi = fPsiMin + G4UniformRand()*(fPsiMax - fPsiMin);

  G4double ux = 1*m*sinAlpha*std::cos(psi),
           uy = 1*m*sinAlpha*std::sin(psi),
           uz = 1*m*cosAlpha;

  fParticleGun->SetParticlePosition(G4ThreeVector(-ux,-uy,-uz));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz-141.5*mm));
	

        fParticleGun->SetParticleTime(i/841.25*s); 



      G4double pEnergy, y, f;

	do {
	  pEnergy = G4UniformRand() * 9.95 * GeV + 0.05*GeV;
	  if (pEnergy < 0.5*GeV)
		{
		f = std::pow(pEnergy * (1/GeV), 0.6);
		y = G4UniformRand()*0.66;
		}
	  else if (pEnergy >= 0.5*GeV)
		{
		f = std::pow(pEnergy * (1/GeV), -1.1);
	 	y = G4UniformRand()*2.15;
		}

	} while (y > f);

	fParticleGun->SetParticleEnergy(pEnergy);
	


  // Cancel randomization
  fParticleGun->SetParticlePosition(G4ThreeVector(2.88*G4UniformRand()*mm,2.88*G4UniformRand()*mm,0.01*mm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));

  fParticleGun->GeneratePrimaryVertex(anEvent);
  i++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

