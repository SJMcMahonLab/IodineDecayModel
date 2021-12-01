
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

#ifndef NikjooTuple_hh
#define NikjooTuple_hh

#include "TsVNtupleScorer.hh"

class NikjooTuple : public TsVNtupleScorer
{
public:
	NikjooTuple(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	
	virtual ~NikjooTuple();
	
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	
private:
    G4Material* fStrandMaterial;

	G4int    fEvt;
	G4String fParticleName;
	G4String fProcessName;
	G4String fVolumeName;
	G4int    fVolumeCopyNumber;
	G4float  fPosX;
	G4float  fPosY;
	G4float  fPosZ;
	G4float  fKineticEnergy;
	G4float  fEnergyDeposited;
	G4float  fTime;
	G4int    fTrackID;
	G4int    fParentAID;
	G4int    fParentBID;
	G4int    fMoleculeID;
	G4int    fStepNumber;
	G4float  fVertexPositionX;
	G4float  fVertexPositionY;
	G4float  fVertexPositionZ;
	
	G4double fTimeCut;
	
private:
	G4bool fIncludeKineticEnergy;
	G4bool fIncludeEventID;
	G4bool fIncludeTrackID;
	G4bool fIncludeParentID;
	G4bool fIncludeStepNumber;
	G4bool fIncludeParticleName;
	G4bool fIncludeProcessName;
	G4bool fIncludeVolumeName;
	G4bool fIncludeVolumeCopyNumber;
	G4bool fIncludeGlobalTime;
	G4bool fIncludeEnergyDeposited;
	G4bool fIncludeVertex;

};
#endif

