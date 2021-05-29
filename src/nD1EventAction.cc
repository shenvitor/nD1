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
/// \file nD1EventAction.cc
/// \brief Implementation of the nD1EventAction class

#include "nD1EventAction.hh"
#include "nD1TrackerHit.hh"
#include "nD1PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4Track.hh"

#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nD1EventAction::nD1EventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nD1EventAction::~nD1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nD1EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  G4int eventID = event->GetEventID();
  G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
    }
    G4cout << "    "  
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }
  
  double Edep = 0., PosX=0., PosY=0., PosZ=0.;
  int copyNb=0;
  if(hc->GetSize()) {
     auto hit = (nD1TrackerHit*)(hc->GetHit(0));
     PosX = hit->GetPos().x(); 
     PosY = hit->GetPos().y(); 
     PosZ = hit->GetPos().z(); 
     copyNb = hit->GetChamberNb();
     }

  //double Energy = (nD1PrimaryGeneratorAction*)(hc->GetHit(0))->GetParticleEnergy();
  //Eesc = Energy;

  for(auto i=0; i<hc->GetSize(); ++i) 
  Edep+=((nD1TrackerHit*)(hc->GetHit(i)))->GetEdep();

  //for(auto i=0; i<hc->GetSize(); ++i) 
  //Eesc+=((nD1TrackerHit*)(hc->GetHit(i)))->GetEesc();
  
  auto anaMgr=G4AnalysisManager::Instance();

  anaMgr->FillH1(1, Edep);
  anaMgr->FillNtupleDColumn(0, Edep);
  anaMgr->FillNtupleDColumn(1, PosX);
  anaMgr->FillNtupleDColumn(2, PosY);
  anaMgr->FillNtupleDColumn(3, PosZ);
  anaMgr->FillNtupleIColumn(4, copyNb);

  for(auto i=0; i<hc->GetSize(); ++i)
  {
    auto hit=(nD1TrackerHit*)(hc->GetHit(i));
    auto trackID=hit->GetTrackID();
    if(trackID==1&&hit->GetIsBoundary())
    {
      anaMgr->FillNtupleIColumn(5, 1);
      anaMgr->FillNtupleDColumn(6, hit->GetDirection().getTheta()/CLHEP::degree);
      anaMgr->FillNtupleDColumn(7, hit->GetDirection().getPhi()/CLHEP::degree);
      anaMgr->FillNtupleDColumn(8, hit->GetPos().x());
      anaMgr->FillNtupleDColumn(9, hit->GetPos().y()); 
      anaMgr->FillNtupleDColumn(10,hit->GetPos().z());
      anaMgr->FillNtupleDColumn(11,hit->GetEkin());
      anaMgr->FillNtupleDColumn(12,hit->GetMomentum().x());
      anaMgr->FillNtupleDColumn(13,hit->GetMomentum().y());
      anaMgr->FillNtupleDColumn(14,hit->GetMomentum().z());
      break;
    }
    else
    {
      anaMgr->FillNtupleIColumn(5, 0);
      continue;
    }
  }
  anaMgr->AddNtupleRow();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
