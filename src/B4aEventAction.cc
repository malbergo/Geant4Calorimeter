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
// $Id: B4aEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class

#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  m_Z.clear();
  m_Energy.clear();
  m_Layer.clear();
  m_AbsEnergy.clear();
  m_GapEnergy.clear();

  for (unsigned int i = 0; i < 30; i++)
  {
    m_AbsEnergy.push_back(0);
    m_GapEnergy.push_back(0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  analysisManager->FillH1(1, fEnergyAbs);
  analysisManager->FillH1(2, fEnergyGap);
  analysisManager->FillH1(3, fTrackLAbs);
  analysisManager->FillH1(4, fTrackLGap);

  for (unsigned int i = 0; i < m_Z.size(); i++)
  {
    analysisManager->FillH1(5, m_Z.at(i));
    analysisManager->FillH1(6, m_Layer.at(i));
    analysisManager->FillH1(7, m_Energy.at(i));
  }
  
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  analysisManager->FillNtupleDColumn(1, fEnergyGap);
  analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  analysisManager->FillNtupleDColumn(3, fTrackLGap);

  for (unsigned int i = 0; i < 30; i++)
  {
    analysisManager->FillNtupleDColumn((2*i)+4, m_AbsEnergy.at(i));
    analysisManager->FillNtupleDColumn((2*i)+5, m_GapEnergy.at(i));
  }

  analysisManager->AddNtupleRow();  

  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLGap,"Length")
       << G4endl;
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
