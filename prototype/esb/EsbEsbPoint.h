/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *              GNU Lesser General Public Licence (LGPL) version 3,             *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef NEWDETECTORPOINT_H
#define NEWDETECTORPOINT_H 1


#include "FairMCPoint.h"

#include "TObject.h"
#include "TVector3.h"

class EsbEsbPoint : public FairMCPoint
{

  public:

    /** Default constructor **/
    EsbEsbPoint();


    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param pos      Ccoordinates at entrance to active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/
    EsbEsbPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
                     Double_t tof, Double_t length, Double_t eLoss);




    /** Destructor **/
    virtual ~EsbEsbPoint();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

  private:
    /** Copy constructor **/
    EsbEsbPoint(const EsbEsbPoint& point);
    EsbEsbPoint operator=(const EsbEsbPoint& point);

    ClassDef(EsbEsbPoint,1)

};

#endif
