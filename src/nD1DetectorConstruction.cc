//// based on template provided by G4 example B2a.
// adapted and developed by Vitor Jose Shen
// at department of physics, Tsinghua University.
// This simulation program "nD1" stands for:
// Simulation of "1(Single) neutron (liquid scintillation) Detector"
// 28 Dec 2020
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
//
/// \file nD1DetectorConstruction.cc
/// \brief Implementation of the nD1DetectorConstruction class
 
#include "nD1DetectorConstruction.hh"
#include "nD1DetectorMessenger.hh"
#include "nD1TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* nD1DetectorConstruction::fMagFieldMessenger = 0;

nD1DetectorConstruction::nD1DetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new nD1DetectorMessenger(this);

//  fNbOfChambers = 1;									// Altering num#########
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
nD1DetectorConstruction::~nD1DetectorConstruction()
{
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* nD1DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1DetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Xe");
  
  //Adding what we needed #########################
  // Hydrogen defined using NIST Manager
  G4Material* H = nistManager->FindOrBuildMaterial("G4_H");
  // Carbon defined using NIST Manager
  G4Material* C = nistManager->FindOrBuildMaterial("G4_C");
  // Nitrogen
  G4Material* N = nistManager->FindOrBuildMaterial("G4_N");
  // Oxygen
  G4Material* O = nistManager->FindOrBuildMaterial("G4_O");
  // Aluminum
  G4Material* Al = nistManager->FindOrBuildMaterial("G4_Al"); 
  // Concrete
  G4Material* Concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE"); 
  //#####
  //how it works!
/*  G4NistManager* man = G4NistManager::Instance();
  G4Material* H2O = man->FindOrBuildMaterial("G4_WATER"); 
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");   */
  //G4Material* air  = G4Material::GetMaterial("G4_AIR");
  
  //Elements
  //G4double z, a;
  G4double z, a;
  G4String name, symbol;
  a = 1.0079*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z=1, a);
  a=12.0107*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6, a);
//  a=14.0067*g/mole;
//  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7, a);
//  a=15.9994*g/mole;
//  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8, a);

  
  //G4double z, a;
  G4double fractionmass, density;
  G4int ncomponents;
  
  G4int iz, n;
  G4Isotope* H2 = new G4Isotope(name="H2Isotope",iz=1,n=2,a=2.0141*g/mole);
  G4Element* D  = new G4Element(name="deuterium",symbol="D",ncomponents=1);
  D->AddIsotope(H2,100. *perCent);

  //***********************************************************************
  //Liquid Scintillator Material C8H10 & C6D6 Below
  //***********************************************************************
  //NE-213 or BC-501A  (C8H10)
  density = 0.874*g/cm3;
  G4Material* NE213 = new G4Material(name="NE213",density,ncomponents=2);
  // H:C \approx 4.82 : 3.98 or 1.213 (in different measurement data)
  // H \approx 54.812%
  // C \approx 45.188%
  NE213->AddElement(elH,fractionmass=54.812*perCent);
  NE213->AddElement(elC, fractionmass=45.188*perCent);

  //NE-230 or BC-537  (C6D6)
  density = 0.945*g/cm3;
  G4Material* NE230 = new G4Material(name="NE230",density,ncomponents=3);
  //C : D : H = 4.10 : 4.06 : 0.0287
  // H \approx 0.35%
  // D \approx 49.58%
  // C \approx 50.07%
  NE230->AddElement(elH,fractionmass=0.35*perCent);
  NE230->AddElement(D,fractionmass=49.58*perCent);
  NE230->AddElement(elC, fractionmass=50.07*perCent);
  //***********************************************************************

  // define a vacuum with a restgas pressure  typical for accelerators
  G4double const Torr  = atmosphere/760.;         // 1 Torr
  G4double pressure = 10e-9*Torr, temperature = 296.150*kelvin;    // 23  Celsius
  G4Material* Vacuum = new G4Material("Vacuum", z=7., a=14.01*g/mole, density= 1.516784e-11*kg/m3,
                 kStateGas, temperature, pressure);

  // Testing for Hydrogen and Cabron
  nistManager->FindOrBuildMaterial("G4_H");
  nistManager->FindOrBuildMaterial("G4_C");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* nD1DetectorConstruction::DefineVolumes()
{
  G4Material* Air  = G4Material::GetMaterial("G4_AIR");
  G4Material* Al  = G4Material::GetMaterial("G4_Al");
  G4Material* Concrete  = G4Material::GetMaterial("G4_CONCRETE");
  // Sizes of the principal geometrical components (solids)
  
  //G4double chamberSpacing = 80*cm; // from chamber center to center!

  //G4double chamberWidth = 20.0*cm; // width of the chambers  ######
  //G4double targetLength =  5.0*cm; // full length of Target  ######
  
//  G4double trackerLength = (fNbOfChambers)*chamberSpacing;	//fNbOfChamber+1 #############
	
//  G4double worldLength = 1.2 * (2*targetLength + trackerLength);   //####
//	G4double worldLength = 1080*cm;

  	G4double worldLength = 2*106.35*cm;
  	//G4double worldLength = 2*112.7*cm;
  	//G4double worldLength = 2*125.4*cm;
  	
//	G4double InnerL = 4.6*cm;
//	G4double OutterL = 6.6*cm;

  //G4double targetRadius  = 0.5*targetLength;   // Radius of Target		   ##########
  //targetLength = 0.5*targetLength;             // Half length of the Target  ##########
  //G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker  which is 200cm=2m
  
  // NE213: length=diameter || 2*length=diameter
  
  G4double trackerDiam = 12.7*cm;
  //G4double trackerDiam = 25.4*cm;

  //G4double trackerLength = 12.7*cm;
  G4double trackerLength = 25.4*cm;
  //G4double trackerLength = 50.8*cm;

  G4double trackerRadius = 0.5*trackerDiam;

  G4double trackerSize = 0.5*trackerLength;

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  /*G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;*/												//#####

 
  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldLength/2,worldLength/2,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 Air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps           

  // Target
  
 /* G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetLength+trackerSize));

  G4Tubs* targetS
    = new G4Tubs("target",0.,targetRadius,targetLength,0.*deg,360.*deg);
  fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
  new G4PVPlacement(0,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << 2*targetLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;							*///#########
  
  //Container
  
  
  //to be continued

  // Tracker
 																		//what about only tracker
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  G4Tubs* trackerS
    = new G4Tubs("tracker",			//pName: it's name	
    			 0,					//pRmin: inner radius
    			 trackerRadius,		//pRmax: outer radius
    			 trackerSize,		//pDz: half length in z
    			 0.*deg, 			//pSphi: phi angle in radius
    			 360.*deg);			//pDphi: angle of segment of radians
    			 
  G4LogicalVolume* trackerLV
  = new G4LogicalVolume(trackerS, G4Material::GetMaterial("NE213"), "Tracker",0,0,0);  
  //= new G4LogicalVolume(trackerS, G4Material::GetMaterial("NE230"), "Tracker",0,0,0);  
  //= new G4LogicalVolume(trackerS, G4Material::GetMaterial("G4_WATER"), "Tracker",0,0,0);  
  
  new G4PVPlacement(0,               // no rotation
                    positionTracker, // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
//  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));  //#####

  worldLV      ->SetVisAttributes(boxVisAtt);
//  fLogicTarget ->SetVisAttributes(boxVisAtt);  //####
  trackerLV    ->SetVisAttributes(boxVisAtt);

  // Tracker segments

/*  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << G4endl
         << "The chambers are " << chamberWidth/cm << " cm of "
         << fChamberMaterial->GetName() << G4endl
         << "The distance between chamber is " << chamberSpacing/cm << " cm" 
         << G4endl;																*/ //#####
  
  /*G4double firstPosition = -trackerSize + chamberSpacing;
  G4double firstLength   = trackerLength/10;
  G4double lastLength    = trackerLength;

  G4double halfWidth = 0.5*chamberWidth;
  G4double rmaxFirst = 0.5 * firstLength;

  G4double rmaxIncr = 0.0;
  if( fNbOfChambers > 0 ){
    rmaxIncr =  0.5 * (lastLength-firstLength)/(fNbOfChambers-1);
    if (chamberSpacing  < chamberWidth) {
       G4Exception("nD1DetectorConstruction::DefineVolumes()",
                   "InvalidSetup", FatalException,
                   "Width>Spacing");
    }
  }																			*/ //########

//  for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {			//#####################

      //G4double Zposition = firstPosition + copyNo * chamberSpacing; //###################
      //G4double Zposition = firstPosition + chamberSpacing;
      //G4double rmax =  rmaxFirst + copyNo * rmaxIncr;				  //###################
//      G4double rmax =  rmaxFirst + rmaxIncr;

//      G4Tubs* chamberS
//        = new G4Tubs("Chamber_solid", 0, rmax, halfWidth, 0.*deg, 360.*deg);

      //fLogicChamber[copyNo] =											//#################
      //        new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber_LV",0,0,0);

    //  fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);			//#################

//      new G4PVPlacement(0,                            // no rotation
//                        G4ThreeVector(0,0,Zposition), // at (x,y,z)
//                        fLogicChamber
                        //fLogicChamber[copyNo],        // its logical volume
//                        "Chamber_PV",                 // its name
//                        trackerLV,                    // its mother  volume
//                        false,                        // no boolean operations
//                        0,
                        //copyNo,                       // copy number   ###################
//                        fCheckOverlaps);              // checking overlaps 

//  }																	#####################

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  //G4double maxStep = 0.5*chamberWidth;	//###
  //fStepLimit = new G4UserLimits(maxStep);	//###
  //trackerLV->SetUserLimits(fStepLimit);	//###
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void nD1DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "nD1/TrackerChamberSD";
  nD1TrackerSD* aTrackerSD = new nD1TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
//  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);		//##############################
  SetSensitiveDetector("Tracker", aTrackerSD, true);

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
 
void nD1DetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1DetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
