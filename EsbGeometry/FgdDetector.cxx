/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *              GNU Lesser General Public Licence (LGPL) version 3,             *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#include "EsbGeometry/FgdDetector.h"
#include "EsbData/FgdDetectorPoint.h" 

#include "FairVolume.h"
#include "FairRootManager.h"
#include "FairGenericStack.h"
#include "FairStack.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoBuilder.h"
#include "FairGeoMedia.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGraph.h"

#include <iostream>
using std::cout;
using std::endl;

namespace esbroot {

namespace geometry {
// PC: work around
static const Int_t kFgdDetector = 2;

//___________________________________________________________________
FgdDetector::FgdDetector()
  : FairDetector("FgdDetector", kTRUE, kFgdDetector),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fFgdDetectorPointCollection(new TClonesArray(data::FgdDetectorPoint::Class()))
{
}

//___________________________________________________________________
  FgdDetector::FgdDetector(const char* name, Bool_t active)
  : FairDetector(name, active, kFgdDetector),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fFgdDetectorPointCollection(new TClonesArray(data::FgdDetectorPoint::Class())) 
{
}

//___________________________________________________________________
FgdDetector::~FgdDetector()
{
  if (fFgdDetectorPointCollection) {
    fFgdDetectorPointCollection->Delete();
    delete fFgdDetectorPointCollection;
  }
}

//___________________________________________________________________
void FgdDetector::Initialize()
{
  FairDetector::Initialize();
}

//___________________________________________________________________
Bool_t  FgdDetector::ProcessHits(FairVolume* vol)
{
	cout << __PRETTY_FUNCTION__ << endl;
  /** This method is called from the MC stepping */

  if ( TVirtualMC::GetMC()->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
    fLength = TVirtualMC::GetMC()->TrackLength();
    TVirtualMC::GetMC()->TrackPosition(fPos);
    TVirtualMC::GetMC()->TrackMomentum(fMom);
  }

  // Sum energy loss for all steps in the active volume
  fELoss += TVirtualMC::GetMC()->Edep();

  // Create FairTutorialDet1Point at exit of active volume
  if ( TVirtualMC::GetMC()->IsTrackExiting()    ||
       TVirtualMC::GetMC()->IsTrackStop()       ||
       TVirtualMC::GetMC()->IsTrackDisappeared()   ) {

    fTrackID  = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
    fVolumeID = vol->getMCid();
    //~ if (fELoss == 0. ) { return kFALSE; }
    AddHit(fTrackID, fVolumeID, TVector3(fPos.X(),  fPos.Y(),  fPos.Z()),
           TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime);
  }

  return kTRUE;
}

void FgdDetector::EndOfEvent()
{
  fFgdDetectorPointCollection->Clear();
}



void FgdDetector::Register()
{

  /** This will create a branch in the output tree called
      EsbFgdDetectorPoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("EsbFgdDetectorPoint", "FgdDetector", 
                                        fFgdDetectorPointCollection, kTRUE);

}


TClonesArray* FgdDetector::GetCollection(Int_t iColl) const
{
  if (iColl == 0) { return fFgdDetectorPointCollection; }
  else { return NULL; }
}

void FgdDetector::Reset()
{
  fFgdDetectorPointCollection->Clear();
}

void FgdDetector::ConstructGeometry()
{
	//TODO3: Create the real Fgd geometry
	FairGeoLoader *geoLoad = FairGeoLoader::Instance();
	FairGeoInterface *geoFace = geoLoad->getGeoInterface();
	
	FairGeoMedia *geoMedia = geoFace->getMedia();
	FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();
	
	FairGeoMedium* mWC = geoMedia->getMedium("H2O_ESSnuSB");
	geoBuild->createMedium(mWC);
  TGeoMedium *WC_med = gGeoManager->GetMedium("H2O_ESSnuSB");
  
	//TODO: Change this to use media.geo file
  TGeoMaterial *Al_mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *Al_med = new TGeoMedium("Al", 101, Al_mat);
  
  // Create water cylinder
  TGeoVolume *wc = gGeoManager->MakeTube("wc", WC_med, 0.0, 1, 1);
  
  // Create thin wall around the water cylinder, 1 cm thick, to act as sensitive volume
  TGeoVolume *wall = gGeoManager->MakeTube("wall", Al_med, 1, 1+0.1, 1);
  TGeoVolume *endwall = gGeoManager->MakeTube("endwall", Al_med, 0, 1, 0.05);

  TList* media = gGeoManager->GetListOfMedia();
  for(TObject *obj : *media) {
		obj->Print();
	}	

  AddSensitiveVolume(wall); //From FairModule
  AddSensitiveVolume(endwall); //From FairModule

  //TODO: Top volume should be a parameter in the constructor
  TGeoVolume *top = gGeoManager->GetTopVolume();
  top->AddNode(wc, 1);
  top->AddNode(wall, 1); 
  top->AddNode(endwall, 1, new TGeoTranslation(0.0, 0.0, 1+0.05));
  top->AddNode(endwall, 2, new TGeoTranslation(0.0, 0.0, -1-0.05));

  wc->SetLineColor(kRed);

}

//___________________________________________________________________
data::FgdDetectorPoint* FgdDetector::AddHit(Int_t trackID, Int_t detID, 
					  TVector3 pos, TVector3 mom,
					  Double_t time)
{
  TClonesArray& clref = *fFgdDetectorPointCollection;
  Int_t size = clref.GetEntriesFast();

  return new(clref[size]) data::FgdDetectorPoint(trackID, detID, pos, mom, 
					     time);
}

}
}

