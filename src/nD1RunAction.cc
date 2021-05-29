// based on template provided by G4 example B2a.
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
/// \file nD1RunAction.cc
/// \brief Implementation of the nD1RunAction class

#include "nD1RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include <g4root.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nD1RunAction::nD1RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);  
  // Create analysis manager
	auto anaMgr = G4AnalysisManager::Instance();
	anaMgr->SetFileName("neutron");
	anaMgr->SetVerboseLevel(1);
	anaMgr->SetFirstHistoId(1);
	anaMgr->SetNtupleMerging(true);
	anaMgr->CreateH1("Edep", "Total Energy Deposit of Neutron [MeV]", 100, 0., 120.);
	anaMgr->CreateNtuple("neutron", "Detected Neutron");
	anaMgr->CreateNtupleDColumn("Edep");//0
	anaMgr->CreateNtupleDColumn("PosX");//1
	anaMgr->CreateNtupleDColumn("PosY");//2
	anaMgr->CreateNtupleDColumn("PosZ");//3
	anaMgr->CreateNtupleIColumn("copyNb");//4
	anaMgr->CreateNtupleIColumn("Escaped");//5
	anaMgr->CreateNtupleDColumn("EscapedTheta");//6
	anaMgr->CreateNtupleDColumn("EscapedPhi");//7
	anaMgr->CreateNtupleDColumn("EscapedPosX");//8
	anaMgr->CreateNtupleDColumn("EscapedPosY");//9
	anaMgr->CreateNtupleDColumn("EscapedPosZ");//10
	anaMgr->CreateNtupleDColumn("Eesc");//11
	anaMgr->CreateNtupleDColumn("EscMomentumX");//12
    anaMgr->CreateNtupleDColumn("EscMomentumY");//13
	anaMgr->CreateNtupleDColumn("EscMomentumZ");//14
	anaMgr->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nD1RunAction::~nD1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  // Initialize Analysis Manager
  auto anaMgr = G4AnalysisManager::Instance();
  anaMgr->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1RunAction::EndOfRunAction(const G4Run* )
{
	auto anaMgr = G4AnalysisManager::Instance();
	anaMgr->Write();
	anaMgr->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
