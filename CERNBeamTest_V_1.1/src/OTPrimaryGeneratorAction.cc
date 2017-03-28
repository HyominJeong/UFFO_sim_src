#include "OTPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

OTPrimaryGeneratorAction::OTPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable -> FindParticle(particleName = "alpha");

  fParticleGun -> SetParticleDefinition(particle);
  fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fParticleGun -> SetParticleEnergy(63.2*GeV);
  G4double positionDeltaX;
  G4double positionDeltaY;
  G4double energyDelta;
  positionDeltaX = G4UniformRand();
  positionDeltaY = G4UniformRand();
  energyDelta    = G4UniformRand();
  //fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(velosityX,velosityY,velosityZ));
  fParticleGun -> SetParticlePosition(G4ThreeVector((positionDeltaX - 0.5)*cm, 
                                                    (positionDeltaY - 0.5)*cm,
                                                    (0.                  )));
  fParticleGun -> SetParticleEnergy((63.2 + 10. * (energyDelta - 0.5)) * GeV);
}

OTPrimaryGeneratorAction::~OTPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void OTPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  // Modify Energy, initial position of particle
  G4double positionDeltaX;
  G4double positionDeltaY;
  G4double energyDelta;
  positionDeltaX = G4UniformRand();
  positionDeltaY = G4UniformRand();
  energyDelta    = G4UniformRand();
  //fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(velosityX,velosityY,velosityZ));
  fParticleGun -> SetParticlePosition(G4ThreeVector((positionDeltaX - 0.5)*cm, 
						                                        (positionDeltaY - 0.5)*cm,
						                                        (0.                  )));
  fParticleGun -> SetParticleEnergy((63.2 - 0.5 + energyDelta) * GeV);


  fParticleGun -> GeneratePrimaryVertex(anEvent);
}
