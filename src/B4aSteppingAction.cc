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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  G4VPhysicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

//  std::cout << "Step Length / mm     : " << stepLength << std::endl;

//  std::cout << "Step Z position : " << step->GetPreStepPoint()->GetPosition().z() << std::endl;
//  std::cout << "Step Energy     : " << edep << std::endl;

  G4double z = step->GetPreStepPoint()->GetPosition().z();
  G4int layer = this->getLayerNumber(z,volume);

  fEventAction->AddEZ(edep, z, layer);

  if ( volume == fDetConstruction->GetAbsorberPV() or volume == fDetConstruction->GetAbsorberPV2() ) {
    fEventAction->AddAbsE(edep,layer);
  }
  else if ( volume == fDetConstruction->GetGapPV() or volume == fDetConstruction->GetGapPV2() ) {
    fEventAction->AddGapE(edep,layer);
  }

//  std::cout << "Z : " << z << std::endl;
//  std::cout << "Layer number : " << this->getLayerNumber(z,volume) << std::endl;

  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    fEventAction->AddAbs(edep,stepLength);
    //std::cout << "Step Length / mm     : " << stepLength/mm << std::endl;
  }

  if ( volume == fDetConstruction->GetAbsorberPV2() ) {
    fEventAction->AddAbs2(edep,stepLength);
    //std::cout << "Step Length / mm     : " << stepLength/mm << std::endl;
  }
  
  if ( volume == fDetConstruction->GetGapPV() ) {
    fEventAction->AddGap(edep,stepLength);
    //std::cout << "Step Length / mm     : " << stepLength/mm << std::endl;
  }

  if ( volume == fDetConstruction->GetGapPV2() ) {
    fEventAction->AddGap2(edep,stepLength);
    //std::cout << "Step Length / mm     : " << stepLength/mm << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int B4aSteppingAction::getLayerNumber(float z, G4VPhysicalVolume* volume)
{
    const G4double calorThickness = fDetConstruction->getCalorThickness();
    const G4double layerThickness = fDetConstruction->getLayerThickness();
    const G4double layerThickness2 = fDetConstruction->getLayerThickness2();
    const G4int nofLayers = fDetConstruction->getNumberOfLayers();
    G4int layer(-2);

    if (volume == fDetConstruction->GetAbsorberPV() or volume == fDetConstruction->GetGapPV())
    {
        layer = floor((z + calorThickness/2) / layerThickness);
    }

    else if ((volume == fDetConstruction->GetAbsorberPV2() or volume == fDetConstruction->GetGapPV2()))
    {
        layer = floor((z - calorThickness/2) / layerThickness2) + nofLayers;
    }

    return layer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
