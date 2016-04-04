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
// $Id: B4DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4DetectorConstruction.hh
/// \brief Definition of the B4DetectorConstruction class

#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <math.h>

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4UserLimits;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class B4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B4DetectorConstruction();
    virtual ~B4DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // get methods
    //
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetAbsorberPV2() const;
    const G4VPhysicalVolume* GetGapPV() const;
    const G4VPhysicalVolume* GetGapPV2() const;
 
    G4double getCalorThickness() const;
    G4double getLayerThickness() const;
    G4double getLayerThickness2() const;
    G4int getNumberOfLayers() const;

  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger

    const G4int m_nofLayers;
    const G4int m_nofLayers2;
    const G4double m_absoThickness;
    const G4double m_absoThickness2;
    const G4double m_gapThickness;
    const G4double m_gapThickness2;
    const G4double m_calorSizeXY;

    G4double m_layerThickness;
    G4double m_layerThickness2;
    G4double m_calorThickness;
    G4double m_calorThickness2;
     
    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4VPhysicalVolume*   fAbsorberPV2;// the absorber physical volume 2
    G4VPhysicalVolume*   fGapPV;      // the gap physical volume
    G4VPhysicalVolume*   fGapPV2;     // the gap physical volume 2
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV() const { 
  return fAbsorberPV; 
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV2() const {
  return fAbsorberPV2;
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV() const  { 
  return fGapPV; 
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV2() const  {
  return fGapPV2;
}

inline G4double B4DetectorConstruction::getCalorThickness() const {
  return m_calorThickness;
}

inline G4double B4DetectorConstruction::getLayerThickness() const {
  return m_layerThickness;
}

inline G4double B4DetectorConstruction::getLayerThickness2() const {
  return m_layerThickness2;
}

inline G4int B4DetectorConstruction::getNumberOfLayers() const {
  return m_nofLayers;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

