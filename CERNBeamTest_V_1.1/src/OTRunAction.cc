#include "OTRunAction.hh"

OTRunAction::OTRunAction()
: G4UserRunAction()
{

}

OTRunAction::~OTRunAction()
{
  delete G4AnalysisManager::Instance();
}

void OTRunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile("data");
  // Creating histograms



  // Creating ntuple
  analysisManager -> CreateNtuple("0", "step");
  analysisManager -> CreateNtupleIColumn("eventID");
  analysisManager -> CreateNtupleIColumn("volumeID");
  analysisManager -> CreateNtupleDColumn("edep");
  analysisManager -> FinishNtuple();



}

void OTRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> Write();
  analysisManager -> CloseFile();
}
