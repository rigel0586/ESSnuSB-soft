/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *              GNU Lesser General Public Licence (LGPL) version 3,             *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef ESBWCDETECTORPOINT_H
#define ESBWCDETECTORPOINT_H 1


#include "FairMCPoint.h"

#include "TObject.h"
#include "TVector3.h"

namespace esbroot {

class EsbWCDetectorPoint : public FairMCPoint
{

  public:

    /** Default constructor **/
    EsbWCDetectorPoint();


    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param pos      Ccoordinates at entrance to active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     **/
    EsbWCDetectorPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
		       Double_t tof);
    
    /** Destructor **/
    virtual ~EsbWCDetectorPoint();


    Bool_t IsSortable() const { return kTRUE; }

    /// To be able to sort the photon hits according to their z position
    Int_t Compare(const TObject *obj) const;
    
    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

  private:
    /** Copy constructor **/
    EsbWCDetectorPoint(const EsbWCDetectorPoint& point);
    EsbWCDetectorPoint operator=(const EsbWCDetectorPoint& point);

    ClassDef(EsbWCDetectorPoint,2)
};

}

#endif
