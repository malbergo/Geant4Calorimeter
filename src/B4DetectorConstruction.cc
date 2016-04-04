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
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   m_nofLayers(799),
   m_nofLayers2(1),
   m_absoThickness(0.21*mm),
   m_absoThickness2(0.21*mm),
   m_gapThickness(0.05*mm),
   m_gapThickness2(0.05*mm),
   m_calorSizeXY(100.*cm),
   fAbsorberPV(0),
   fAbsorberPV2(0),
   fGapPV(0),
   fGapPV2(0),
   fCheckOverlaps(true)
{
    m_layerThickness = m_absoThickness + m_gapThickness;
    m_layerThickness2 = m_absoThickness2 + m_gapThickness2;
    m_calorThickness = m_nofLayers * m_layerThickness;
    m_calorThickness2 = m_nofLayers2 * m_layerThickness2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

// Tungsten Atomic Weight : 183.84 g/mol
// Tungsten Density : 19.3 grams per cubic centimeter
  G4double tungsten_atomic_weight = 183.84*g/mole;
  G4double tungsten_density = 19.3*g/cm3;
  new G4Material("tungsten_19.3gccm", z=74, tungsten_atomic_weight, tungsten_density);

// Silicon Atomic Weight : 183.84 g/mol
// Silicon Density : 20.0855 grams per cubic centimeter
  G4double silicon_atomic_weight = 20.0855*g/mole;
  G4double silicon_density = 2.33*g/cm3;
  new G4Material("silicon_2.33gccm", z=14, silicon_atomic_weight, silicon_density);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double worldSizeXY = 1.2 * m_calorSizeXY;
  G4double worldSizeZ  = 3.2 * m_calorThickness; 
  
  G4double caloOffset = m_calorThickness/2+m_calorThickness2/2;

  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
//  G4Material* absorberMaterial = G4Material::GetMaterial("G4_Pb"); <- Original
  G4Material* absorberMaterial = G4Material::GetMaterial("tungsten_19.3gccm");
//  G4Material* gapMaterial = G4Material::GetMaterial("liquidArgon"); <- Original
  G4Material* gapMaterial = G4Material::GetMaterial("silicon_2.33gccm");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
 
  //                               
  // Calorimeter
  //  
  G4VSolid* calorimeterS
    = new G4Box("Calorimeter",     // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_calorThickness/2); // its size
                         
  G4LogicalVolume* calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  //  Calorimeter 2
  //  
  G4VSolid* calorimeterS2
    = new G4Box("Calorimeter2",     // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_calorThickness2/2); // its size

  G4LogicalVolume* calorLV2
    = new G4LogicalVolume(
                 calorimeterS2,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter2");   // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,0,caloOffset),  // at (0,0,0)
                 calorLV2,         // its logical volume
                 "Calorimeter2",   // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

 
  //                                 
  // Layer
  //
  G4VSolid* layerS 
    = new G4Box("Layer",           // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_layerThickness/2); // its size
                         
  G4LogicalVolume* layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 m_nofLayers,        // number of replica
                 m_layerThickness);  // witdth of replica

  //
  // Layer
  //
  G4VSolid* layerS2
    = new G4Box("Layer2",           // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_layerThickness2/2); // its size

  G4LogicalVolume* layerLV2
    = new G4LogicalVolume(
                 layerS2,          // its solid
                 defaultMaterial,  // its material
                 "Layer2");        // its name

  new G4PVReplica(
                 "Layer2",         // its name
                 layerLV2,         // its logical volume
                 calorLV2,         // its mother
                 kZAxis,           // axis of replication
                 m_nofLayers2,       // number of replica
                 m_layerThickness2); // witdth of replica

 
  //                               
  // Absorber
  //
  G4VSolid* absorberS 
    = new G4Box("Abso",            // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_absoThickness/2); // its size
                         
  G4LogicalVolume* absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Abso");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -m_gapThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Absorber
  //
  G4VSolid* absorberS2
    = new G4Box("Abso2",            // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_absoThickness2/2); // its size

  G4LogicalVolume* absorberLV2
    = new G4LogicalVolume(
                 absorberS2,       // its solid
                 absorberMaterial, // its material
                 "Abso2");         // its name

  fAbsorberPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -m_gapThickness2/2), // its position, relative to mother volume
                 absorberLV2,      // its logical volume
                 "Abso2",          // its name
                 layerLV2,         // its mother volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //                               
  // Gap
  //
  G4VSolid* gapS 
    = new G4Box("Gap",             // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_gapThickness/2); // its size
                         
  G4LogicalVolume* gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "Gap");           // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., m_absoThickness/2), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Gap
  //
  G4VSolid* gapS2
    = new G4Box("Gap2",             // its name
                 m_calorSizeXY/2, m_calorSizeXY/2, m_gapThickness2/2); // its size

  G4LogicalVolume* gapLV2
    = new G4LogicalVolume(
                 gapS2,             // its solid
                 gapMaterial,      // its material
                 "Gap2");          // its name

  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., m_absoThickness2/2), // its position relative to mother volume
                 gapLV2,           // its logical volume
                 "Gap2",           // its name
                 layerLV2,         // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.001*mm;
  G4UserLimits* fStepLimit = new G4UserLimits(maxStep);
  worldLV->SetUserLimits(fStepLimit);
  calorLV->SetUserLimits(fStepLimit);
  calorLV2->SetUserLimits(fStepLimit);
  layerLV->SetUserLimits(fStepLimit);
  layerLV2->SetUserLimits(fStepLimit);
  absorberLV->SetUserLimits(fStepLimit);
  absorberLV2->SetUserLimits(fStepLimit);
  gapLV->SetUserLimits(fStepLimit);
  gapLV2->SetUserLimits(fStepLimit);

  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  //                                           maxLength,
  //                                           maxTime,
  //                                           minEkin));

  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << m_nofLayers << " layers of: [ "
    << m_absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << m_gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);
  calorLV2->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
