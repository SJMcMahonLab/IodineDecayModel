// Scorer for NikjooRadicals
//

// ******************************************************************
// 
// This software is made freely available in accordance with the simplifed BSD
// license:
// 
// Copyright (c) <2021>, <Stephen McMahon>
// All rights reserved
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Contacts: Stephen McMahon,	stephen.mcmahon@qub.ac.uk
// 
// ******************************************************************
//                                                                   
//  This scorer implements a specialised radical scorer based on the 
//  Tuple scorer in TOPAS-nBio, to follow approach of                
//  Nikjoo et al, Acta Oncologica, 1996							  
//  https://doi.org/10.3109/02841869609104036						  
//                                                                   
//  This tracks first radical interactions within a backbone, and    
//  then scavenges them from system to prevent subsequent interaction
//                                                                   
// ******************************************************************
//

#include "NikjooTuple.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Molecule.hh"
#include "G4MoleculeFinder.hh"

NikjooTuple::NikjooTuple(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fIncludeKineticEnergy(false), fIncludeEventID(false),  fIncludeTrackID(false), fIncludeParentID(false), fIncludeStepNumber(false),
fIncludeParticleName(false), fIncludeProcessName(false), fIncludeVolumeName(false), fIncludeVolumeCopyNumber(false),
fIncludeGlobalTime(false), fIncludeEnergyDeposited(false), fIncludeVertex(false)
{
	SetUnit("");

	fNtuple->RegisterColumnI(&fMoleculeID, "MoleculeID or ParticlePDG");
	fNtuple->RegisterColumnF(&fPosX, "Position X", "um");
	fNtuple->RegisterColumnF(&fPosY, "Position Y", "um");
	fNtuple->RegisterColumnF(&fPosZ, "Position Z", "um");
	
	fTimeCut = 1.0 * us;
	if ( fPm->ParameterExists(GetFullParmName("TimeCut") ) )
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"),"Time");
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeKineticEnergy")) )
		fIncludeKineticEnergy = fPm->GetBooleanParameter(GetFullParmName("IncludeKineticEnergy"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeEventID")) )
		fIncludeEventID = fPm->GetBooleanParameter(GetFullParmName("IncludeEventID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeTrackID")) )
		fIncludeTrackID = fPm->GetBooleanParameter(GetFullParmName("IncludeTrackID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeParentID")) )
		fIncludeParentID = fPm->GetBooleanParameter(GetFullParmName("IncludeParentID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeStepNumber")) )
		fIncludeStepNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeStepNumber"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeParticleName")) )
		fIncludeParticleName = fPm->GetBooleanParameter(GetFullParmName("IncludeParticleName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludePhysicalProcessName")) )
		fIncludeProcessName = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalProcessName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVolumeName")) )
		fIncludeVolumeName = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVolumeCopyNumber")) )
		fIncludeVolumeCopyNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeCopyNumber"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeGlobalTime")) )
		fIncludeGlobalTime = fPm->GetBooleanParameter(GetFullParmName("IncludeGlobalTime"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeEnergyDeposited")) )
		fIncludeEnergyDeposited = fPm->GetBooleanParameter(GetFullParmName("IncludeEnergyDeposited"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVertexPosition")) )
		fIncludeVertex = fPm->GetBooleanParameter(GetFullParmName("IncludeVertexPosition"));
	
	if ( fIncludeEventID )
		fNtuple->RegisterColumnI(&fEvt,      "EventID");
	
	if ( fIncludeTrackID )
		fNtuple->RegisterColumnI(&fTrackID,  "Track ID");
	
	if ( fIncludeStepNumber)
		fNtuple->RegisterColumnI(&fStepNumber, "Step number");
	
	if ( fIncludeParticleName )
		fNtuple->RegisterColumnS(&fParticleName, "Particle name");
	
	if ( fIncludeProcessName )
		fNtuple->RegisterColumnS(&fProcessName, "Process name");
		
	if ( fIncludeVolumeName )
		fNtuple->RegisterColumnS(&fVolumeName, "Volume name");
	
	if ( fIncludeVolumeCopyNumber )
		fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume copy number");
	
	if ( fIncludeParentID ) {
		fNtuple->RegisterColumnI(&fParentAID, "ParentA ID");
		fNtuple->RegisterColumnI(&fParentBID, "ParentB ID");
	}
	
	if ( fIncludeVertex ) {
		fNtuple->RegisterColumnF(&fVertexPositionX, "Vertex position x", "um");
		fNtuple->RegisterColumnF(&fVertexPositionY, "Vertex position y", "um");
		fNtuple->RegisterColumnF(&fVertexPositionZ, "Vertex position z", "um");
	}
	
	if ( fIncludeGlobalTime )
		fNtuple->RegisterColumnF(&fTime, "Global time", "ps");
	
	if ( fIncludeEnergyDeposited )
		fNtuple->RegisterColumnF(&fEnergyDeposited, "Energy deposited", "keV");
	
	if ( fIncludeKineticEnergy )
		fNtuple->RegisterColumnF(&fKineticEnergy, "Kinetice energy", "keV");

	G4String strandMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("StrandMaterial")))
		strandMaterialName = fPm->GetStringParameter(GetFullParmName("StrandMaterial"));
	fStrandMaterial = GetMaterial(strandMaterialName);
	
}


NikjooTuple::~NikjooTuple() {;}


G4bool NikjooTuple::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	G4Track* aTrack = aStep->GetTrack();
	// Only want radicals, so return if we're a physical particle
	if(aTrack->GetTrackID()>=0) {
		return false;
	}

	G4bool status = false;

	G4StepPoint* preStep = aStep->GetPreStepPoint();
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	fTime = aStep->GetPreStepPoint()->GetGlobalTime();

	
	// Filter out anything created within DNA volume
	G4ThreeVector vpos = aTrack->GetVertexPosition();
	G4double vrSq = vpos.x()*vpos.x()+vpos.y()*vpos.y();

	if ( vrSq <(1.15*nm*1.15*nm)) {
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	} else {
		if ( material == fStrandMaterial ) {
			G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
			fVolumeName = touchable->GetVolume()->GetName();
			fVolumeCopyNumber  = touchable->GetVolume()->GetCopyNo();
			
			G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
			fEvt = GetEventID();
			fKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
			fPosX = pos.x();
			fPosY = pos.y();
			fPosZ = pos.z();
			if ( fIncludeVertex ) {
				G4ThreeVector vpos = aTrack->GetVertexPosition();
				fVertexPositionX = vpos.x();
				fVertexPositionY = vpos.y();
				fVertexPositionZ = vpos.z();
			}
			
			fParentAID = -1;
			fParentBID = -1;
			fTrackID = aTrack->GetTrackID();
			
			fParticleName = GetMolecule(aTrack)->GetName();
			fMoleculeID = GetMolecule(aTrack)->GetMoleculeID();
			
			GetMolecule(aTrack)->GetParentID(fParentAID, fParentBID);
			fEnergyDeposited = 0.0;
			fProcessName = "none";
			fStepNumber = aTrack->GetCurrentStepNumber();
		
			
			if(fParticleName == "OH^0" || fParticleName == "OH^-1") {
				fNtuple->Fill();
				status=true;
			} else {
				status=false;
			}
			aTrack->SetTrackStatus(fStopAndKill);
		}
	}

	G4MoleculeFinder::Instance()->UpdatePositionMap();
	return status;
}
