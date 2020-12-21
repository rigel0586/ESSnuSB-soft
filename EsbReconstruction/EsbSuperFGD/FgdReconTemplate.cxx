#include "EsbReconstruction/EsbSuperFGD/FgdReconTemplate.h"

// FairRoot
#include "FairLogger.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <algorithm>

// Genfit headers
#include "AbsBField.h"
#include "AbsMeasurement.h"
#include "ConstField.h"
#include <Exception.h>
#include <EventDisplay.h>
#include <FieldManager.h>
#include "FitStatus.h"
#include <KalmanFitterRefTrack.h>
#include <KalmanFitter.h>
#include "MaterialEffects.h"
#include "MeasuredStateOnPlane.h"
#include <PlanarMeasurement.h>
#include <RKTrackRep.h>
#include "SpacepointMeasurement.h"
#include <StateOnPlane.h>
#include "TDatabasePDG.h"
#include <TGeoMaterialInterface.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>
#include <TRandom.h>
#include <TVector3.h>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

using namespace std;

#define PERMUTATIONS 23

FgdReconTemplate::FgdReconTemplate() : fgdConstructor("")
{
    LoadTemplates();
}

FgdReconTemplate::FgdReconTemplate(const char* geoConfigFile, const char* graphTrackConfig, bool validateTrackGrad)
    : fvalidateTrackGrad(validateTrackGrad)
    , fgdConstructor(geoConfigFile)
    , fminGenFitInterations(2)
    , fmaxGenFitIterations(4)
    , fDebuglvl_genfit(0)
{
    LoadTemplates();
    fParams.LoadPartParams(geoConfigFile);

    // Get dimentions from geometry file
    flunit = fParams.GetLenghtUnit(); // [cm]

    f_step_X  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_X) * flunit;
    f_step_Y  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_Y) * flunit;
    f_step_Z  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_Z) * flunit;

    f_bin_X = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_X);
    f_bin_Y = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_Y);
    f_bin_Z = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_Z);

    f_total_X = f_step_X * f_bin_X;
    f_total_Y = f_step_Y * f_bin_Y;
    f_total_Z = f_step_Z * f_bin_Z;

    fsmoothDepth = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_SMOOTH_GRAPH_DEPTH);
    fsmoothErrLimit = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::FGD_SMOOTH_GRAPH_ERR_LIMIT);

    initGraphRecord(graphTrackConfig);
}

FgdReconTemplate::~FgdReconTemplate()
{
}

Bool_t FgdReconTemplate::IsLeaf(ReconHit* hit)
{
    Bool_t isHitLeaf(false);

    if(hit->fAllHits.size()==1)
    {
        LOG(debug)<< "Single leaf " << hit->fLocalId;
        isHitLeaf = true;
    }
    else
    {
        Int_t permutation(0);
        std::vector<TVector3> vecs;
        GetHitVectors(hit, vecs);
        for(size_t temp=0; !isHitLeaf && temp < fLeafVectors.size(); ++temp)
        {
            if(fLeafVectors[temp].hitVectors.size() == hit->fAllHits.size())
            {
                std::vector<TVector3>& tempVecs = fLeafVectors[temp].hitVectors;
                isHitLeaf = AreVectorsEqual(tempVecs, vecs, permutation);
            }
        }

        if(isHitLeaf)
        {
            LOG(debug) << "Leaf  " << hit->fLocalId << " with local hits " << hit->fAllHits.size() << endl;
        }
    }
    
    return isHitLeaf;
}



Bool_t FgdReconTemplate::GetNextHit(ReconHit* previous, ReconHit* current, ReconHit*& next)
{
    Bool_t nextFound(false);
    next = nullptr;

    if(current==nullptr)
    {
        throw "Invalid current hit";
    }

    if(current->fAllHits.empty())
    {
        LOG(info) <<"Empty cube ";
        return nextFound;
    }

    Int_t previousId = (previous==nullptr) ? -1 : previous->fLocalId;
    size_t nearestId(0);
    Double_t min_dist = std::numeric_limits<Int_t>::max();

    for(size_t nid = 0; nid< current->fAllHits.size(); ++nid)
    {
        ReconHit* candidateHit = current->fAllHits[nid];

        if(!candidateHit->fIsVisited
            && candidateHit->fLocalId != previousId
            && !candidateHit->fIsLeaf)
        {
            TVector3 vecPosition = current->fmppcLoc - candidateHit->fmppcLoc;
            Double_t dist = vecPosition.X()*vecPosition.X() + vecPosition.Y()*vecPosition.Y() + vecPosition.Z()*vecPosition.Z();

            if(dist < min_dist)
            {
                min_dist = dist;
                nearestId = nid;
            }
            
        }
    }

    next = current->fAllHits[nearestId];
    nextFound = (min_dist!=std::numeric_limits<Int_t>::max());
    if(!nextFound)
    {
        LOG(debug) << "X " << current->fmppcLoc.X() << " Y " << current->fmppcLoc.Y() << " Z " << current->fmppcLoc.Z();
        LOG(debug) <<"No next cubes found " << current->fAllHits.size();
    }
        

    return nextFound;
}

Bool_t FgdReconTemplate::CheckAllVisited(ReconHit* hit)
{
    Bool_t allVisited(false);

    if(hit==nullptr)
    {
        throw "Invalid hit!";
    }

    if(hit->fAllHits.empty())
    {
        return allVisited;
    }

    uint visitHits(0);

    for(size_t i = 0; i< hit->fAllHits.size(); ++i)
    {
        ReconHit* neighbourHit = hit->fAllHits[i];
        if(neighbourHit->fIsVisited)
        {
            ++visitHits;
        }
    }

    allVisited = (visitHits == hit->fAllHits.size());

    return allVisited;
}

void FgdReconTemplate::LoadTemplates()
{
    // ==============================
    // 1. Single leaf templates
    // ==============================
    //     Leaf 
    // --- --- ---
    // --- -C- -X-
    // --- --- ---
    TVector3 leaf1_1(0,0,-1);
    FgdReconTemplate::HitTemplate leaf1_1_temp;
    leaf1_1_temp.hitVectors.emplace_back(leaf1_1);

    fLeafVectors.push_back(leaf1_1_temp);

    // Leaf 
    // --- --- -X-
    // --- -C- ---
    // --- --- ---
    TVector3 leaf1_2(0,-1,-1);
    FgdReconTemplate::HitTemplate leaf1_2_temp;
    leaf1_2_temp.hitVectors.emplace_back(leaf1_2);

    fLeafVectors.push_back(leaf1_2_temp);

    // Leaf 
    // --- --- --X
    // --- -C- ---
    // --- --- ---
    TVector3 leaf1_3(1,-1,-1);
    FgdReconTemplate::HitTemplate leaf1_3_temp;
    leaf1_3_temp.hitVectors.emplace_back(leaf1_3);

    fLeafVectors.push_back(leaf1_3_temp);


    // ==============================
    // 2. Two hits leaf templates
    // ==============================

    // Leaf 
    // --- --- -X-
    // --- -C- -X-
    // --- --- ---
    TVector3 leaf2_1_1(0,-1,-1);
    TVector3 leaf2_1_2(0,0,-1);
    FgdReconTemplate::HitTemplate leaf2_1_temp;
    leaf2_1_temp.hitVectors.emplace_back(leaf2_1_1);
    leaf2_1_temp.hitVectors.emplace_back(leaf2_1_2);

    fLeafVectors.push_back(leaf2_1_temp);

    // Leaf 
    // --- --- --X
    // --- -C- -X-
    // --- --- ---
    TVector3 leaf2_2_1(1,-1,-1);
    TVector3 leaf2_2_2(0,0,-1);
    FgdReconTemplate::HitTemplate leaf2_2_temp;
    leaf2_2_temp.hitVectors.emplace_back(leaf2_2_1);
    leaf2_2_temp.hitVectors.emplace_back(leaf2_2_2);

    fLeafVectors.push_back(leaf2_2_temp);

    // Leaf 
    // --- --- X--
    // --- -C- X--
    // --- --- ---
    TVector3 leaf2_3_1(-1,-1,-1);
    TVector3 leaf2_3_2(-1,0,-1);
    FgdReconTemplate::HitTemplate leaf2_3_temp;
    leaf2_3_temp.hitVectors.emplace_back(leaf2_3_1);
    leaf2_3_temp.hitVectors.emplace_back(leaf2_3_2);

    fLeafVectors.push_back(leaf2_3_temp);

    // Leaf 
    // --- --- --X
    // --- -C- --X
    // --- --- ---
    TVector3 leaf2_4_1(1,-1,-1);
    TVector3 leaf2_4_2(1,0,-1);
    FgdReconTemplate::HitTemplate leaf2_4_temp;
    leaf2_4_temp.hitVectors.emplace_back(leaf2_4_1);
    leaf2_4_temp.hitVectors.emplace_back(leaf2_4_2);

    fLeafVectors.push_back(leaf2_4_temp);

    // ==============================
    // 3. Three hits leaf templates
    // ==============================

    // Leaf 
    // --- --- -XX
    // --- -C- -X-
    // --- --- ---
    TVector3 leaf3_1_1(0,-1,-1);
    TVector3 leaf3_1_2(1,-1,-1);
    TVector3 leaf3_1_3(0,0,-1);
    FgdReconTemplate::HitTemplate leaf3_1_temp;
    leaf3_1_temp.hitVectors.emplace_back(leaf3_1_1);
    leaf3_1_temp.hitVectors.emplace_back(leaf3_1_2);
    leaf3_1_temp.hitVectors.emplace_back(leaf3_1_3);

    fLeafVectors.push_back(leaf3_1_temp);


    // Leaf 
    // --- --- XX-
    // --- -C- -X-
    // --- --- ---
    TVector3 leaf3_2_1(-1,-1,-1);
    TVector3 leaf3_2_2(0,-1,-1);
    TVector3 leaf3_2_3(0,0,-1);
    FgdReconTemplate::HitTemplate leaf3_2_temp;
    leaf3_2_temp.hitVectors.emplace_back(leaf3_2_1);
    leaf3_2_temp.hitVectors.emplace_back(leaf3_2_2);
    leaf3_2_temp.hitVectors.emplace_back(leaf3_2_3);

    fLeafVectors.push_back(leaf3_2_temp);

    LOG(debug) << " Leaf templates found " << fLeafVectors.size();
}

void FgdReconTemplate::GetHitVectors(ReconHit* hit , std::vector<TVector3>& vecs)
{
    for(size_t i=0; i< hit->fAllHits.size(); ++i)
    {
        ReconHit& neightbourHit = *hit->fAllHits[i];
        TVector3 result = hit->fmppcLoc - neightbourHit.fmppcLoc;
        vecs.emplace_back(result);
    }
}

// Compares if the two vectors are equal
// This also includes rotational symmetry
// Make a copy of the template since it will be modified
Bool_t FgdReconTemplate::AreVectorsEqual(const std::vector<TVector3>& tempVecs, const std::vector<TVector3>& vecs,  Int_t& foundPermutation )
{
    Bool_t areEqual(false);

    if(tempVecs.size()!=vecs.size())
    {
        return areEqual;
    }

    std::vector<TVector3> tempVecPermut = tempVecs;

    Int_t permutation(0);
    Int_t limitPermutations = PERMUTATIONS;

    while(!areEqual && permutation<=limitPermutations)
    {
        for(size_t i=0; i<tempVecs.size(); ++i)
        {
            TVector3 tmp =  GetPermutation(tempVecs[i],permutation);
            tempVecPermut[i].SetX(tmp.X());
            tempVecPermut[i].SetY(tmp.Y());
            tempVecPermut[i].SetZ(tmp.Z());
        }

        Bool_t allVecsAreEqual(true);
        for(size_t i=0; allVecsAreEqual && i<vecs.size(); ++i)
        {
            allVecsAreEqual = std::find(tempVecPermut.begin(), tempVecPermut.end(), vecs[i])     !=  tempVecPermut.end();
        }

        areEqual = allVecsAreEqual;

        if(areEqual)
        {
            foundPermutation = permutation;
            break;
        }

        ++permutation;
    }

    return areEqual;
}

TVector3 FgdReconTemplate::GetPermutation(TVector3 vec, Int_t numPermutation)
{
    Double_t rot90deg = TMath::Pi()/2;
    Double_t rot180deg = TMath::Pi();

    if(numPermutation<0 || numPermutation>PERMUTATIONS)
    {
        throw "Invalid permutation";
    }

    switch(numPermutation)
    {
        case 0:
                // No rotation
                break;
        case 1:
                vec.RotateZ(rot90deg);
                break;
        case 2:
                vec.RotateZ(2*rot90deg);
                break;
        case 3:
                vec.RotateZ(3*rot90deg);
                break;
        case 4:
                vec.RotateY(rot90deg);
                break;
        case 5:
                vec.RotateY(rot90deg);
                vec.RotateX(rot90deg);
                break;
        case 6:
                vec.RotateY(rot90deg);
                vec.RotateX(2*rot90deg);
                break;
        case 7:
                vec.RotateY(rot90deg);
                vec.RotateX(3*rot90deg);
                break;
        case 8:
                vec.RotateY(2*rot90deg);
                break;
        case 9:
                vec.RotateY(2*rot90deg);
                vec.RotateZ(rot90deg);
                break;
        case 10:
                vec.RotateY(2*rot90deg);
                vec.RotateZ(2*rot90deg);
                break;
        case 11:
                vec.RotateY(2*rot90deg);
                vec.RotateZ(3*rot90deg);
                break;
        case 12:
                vec.RotateY(3*rot90deg);
                break;
        case 13:
                vec.RotateY(3*rot90deg);
                vec.RotateX(rot90deg);
                break;
        case 14:
                vec.RotateY(3*rot90deg);
                vec.RotateX(2*rot90deg);
                break;
        case 15:
                vec.RotateY(3*rot90deg);
                vec.RotateX(3*rot90deg);
                break;
        case 16:
                vec.RotateX(rot90deg);
                break;
        case 17:
                vec.RotateX(rot90deg);
                vec.RotateY(rot90deg);
                break;
        case 18:
                vec.RotateX(rot90deg);
                vec.RotateY(2*rot90deg);
                break;
        case 19:
                vec.RotateX(rot90deg);
                vec.RotateY(3*rot90deg);
                break;
        case 20:
                vec.RotateX(-rot90deg);
                break;
        case 21:
                vec.RotateX(-rot90deg);
                vec.RotateY(rot90deg);
                break;
        case 22:
                vec.RotateX(-rot90deg);
                vec.RotateY(2*rot90deg);
                break;
        case 23:
                vec.RotateX(-rot90deg);
                vec.RotateY(3*rot90deg);
                break;
        default: 
                break;
    }

    Double_t&& x = std::round(vec.X());
    Double_t&& y = std::round(vec.Y());
    Double_t&& z = std::round(vec.Z());

    return TVector3(x,y,z);
}

void FgdReconTemplate::SmoothGraph(std::vector<ReconHit>& hits)
{
    for(size_t i =0; i< hits.size(); ++i)
    {
        ReconHit& hit = hits[i];
        TVector3 smoth = hit.fsmoothco;
        std::set<Long_t> visited;
        visited.insert(hitId(hit));

        SmoothCoordinate(&hit,smoth, visited);
        hit.fsmoothco.SetX(smoth.X()/visited.size());
        hit.fsmoothco.SetY(smoth.Y()/visited.size());
        hit.fsmoothco.SetZ(smoth.Z()/visited.size());


        // if(hit.getSmoothErr().Mag() > fsmoothErrLimit)
        // {
        //     LOG(info) << "before "<< " X " << hit.fmppcLoc.X() << " Y " << hit.fmppcLoc.Y()<< " Z " << hit.fmppcLoc.Z();
        //     LOG(info) << "after "<< " X " << hit.fsmoothco.X() << " Y " << hit.fsmoothco.Y()<< " Z " << hit.fsmoothco.Z();
        //     LOG(info) << "visited "<< visited.size() << " => " << hit.getSmoothErr().Mag() << " > " << fsmoothErrLimit;
        //     LOG(info) << "===================================================";
        // }
        // LOG(info) << "before "<< " X " << hit.fmppcLoc.X() << " Y " << hit.fmppcLoc.Y()<< " Z " << hit.fmppcLoc.Z();
        // LOG(info) << "after "<< " X " << hit.fsmoothco.X() << " Y " << hit.fsmoothco.Y()<< " Z " << hit.fsmoothco.Z();

        // if(hit.getSmoothErr().Mag() < fsmoothErrLimit) LOG(info) << "Err "<< hit.getSmoothErr().Mag();
        // if(hit.getSmoothErr().Mag() > fsmoothErrLimit) LOG(warning) << "Err "<< hit.getSmoothErr().Mag();
        
        // LOG(info) << "visited "<< visited.size() << " " << fsmoothErrLimit;
        // LOG(info) << "===================================================";
    }
}


void FgdReconTemplate::SmoothCoordinate(ReconHit* hit, TVector3& cord, std::set<Long_t>& visited , size_t depth)
{
    if(depth>fsmoothDepth)
    {
        return;
    }

    
    for(size_t i =0; i< hit->fAllHits.size(); ++i)
    {  
        ReconHit& child = *hit->fAllHits[i];
        const bool isVisited = visited.find(hitId(child)) != visited.end();
        if(isVisited) 
        {
            //LOG(info) << "Continue "<< " X " << child.fmppcLoc.X() << " Y " << child.fmppcLoc.Y()<< " Z " << child.fmppcLoc.Z();
            continue;
        }

        visited.insert(hitId(child));
        Double_t&& x = child.fmppcLoc.X() + cord.X();
        cord.SetX(x);
        Double_t&& y = child.fmppcLoc.Y() + cord.Y();
        cord.SetY(y);
        Double_t&& z = child.fmppcLoc.Z() + cord.Z();
        cord.SetZ(z);
        //LOG(info) << "Visited "<< " X " << child.fmppcLoc.X() << " Y " << child.fmppcLoc.Y()<< " Z " << child.fmppcLoc.Z()
        //            << " => depth "<< depth << " size " << hit->fAllHits.size();
        SmoothCoordinate(&child, cord, visited, depth + 1);
    }
}


void FgdReconTemplate::BuildGraph(std::vector<ReconHit>& hits)
{
    // Create the position to which index in the vector it is pointing
    std::map<Long_t, Int_t> positionToId;
    for(Int_t i=0; i<hits.size(); ++i)
    {
      Int_t&& x = hits[i].fmppcLoc.X();
      Int_t&& y = hits[i].fmppcLoc.Y();
      Int_t&& z = hits[i].fmppcLoc.Z();

      Int_t&& ind = ArrInd(x,y,z);

      // GUARD agains double or more hits in the same cube
      if(positionToId.find(ind)==positionToId.end())
      {
        positionToId[ind] = i;
      }
      

      hits[i].fAllHits.clear(); // Clear previous index positions
      hits[i].fLocalId = i;
    }


    auto checkNext = [&](Int_t x_pos, Int_t y_pos, Int_t z_pos, Int_t ind){
                                                                  Long_t&& key = ArrInd(x_pos,y_pos,z_pos);
                                                                  if(positionToId.find(key)!=positionToId.end())
                                                                  {
                                                                    ReconHit* toAdd = &hits[positionToId[key]];
                                                                    hits[ind].fAllHits.push_back(toAdd);
                                                                  }
                                                                };

    for(Int_t i=0; i<hits.size(); ++i)
    {
      Int_t&& x = hits[i].fmppcLoc.X();
      Int_t&& y = hits[i].fmppcLoc.Y();
      Int_t&& z = hits[i].fmppcLoc.Z();

      // Check in X axis
      checkNext(x+1,y,z, i);
      checkNext(x-1,y,z, i);

      // Check in Y axis
      checkNext(x,y+1,z, i);
      checkNext(x,y-1,z, i);

      // Check in Z axis
      checkNext(x,y,z+1, i);
      checkNext(x,y,z-1, i);

      // Check in X,Y corners
      checkNext(x+1,y+1,z, i);
      checkNext(x+1,y-1,z, i);
      checkNext(x-1,y+1,z, i);
      checkNext(x-1,y-1,z, i);

      // Check in X,Z corners
      checkNext(x+1,y,z+1, i);
      checkNext(x+1,y,z-1, i);
      checkNext(x-1,y,z+1, i);
      checkNext(x-1,y,z-1, i);

      // Check in Y,Z corners
      checkNext(x,y+1,z+1, i);
      checkNext(x,y+1,z-1, i);
      checkNext(x,y-1,z+1, i);
      checkNext(x,y-1,z-1, i);

      // Check in X,Y,Z corners
      checkNext(x+1,y+1,z+1, i);
      checkNext(x+1,y+1,z-1, i);
      checkNext(x+1,y-1,z+1, i);
      checkNext(x+1,y-1,z-1, i);

      checkNext(x-1,y+1,z+1, i);
      checkNext(x-1,y+1,z-1, i);
      checkNext(x-1,y-1,z+1, i);
      checkNext(x-1,y-1,z-1, i);
    }
}

void FgdReconTemplate::FindLeaves(std::vector<ReconHit>& hits)
{
    for(Int_t i=0; i<hits.size(); ++i)
    {
        ReconHit& hit = hits[i];
        hit.fIsLeaf = IsLeaf(&hit);
    }
}


void FgdReconTemplate::FindTracks(std::vector<ReconHit>& hits, std::vector<std::vector<ReconHit*>>& tracks)
{
  ReconHit* currentHit = nullptr;
  ReconHit* nextHit = nullptr;
  ReconHit* previousHit = nullptr;
  for(size_t i=0; i<hits.size(); ++i)
  {
    if(hits[i].fIsVisited)
    {
      continue;
    }

    if(hits[i].fIsLeaf)
    {
      std::vector<ReconHit*> track;
      track.push_back(&hits[i]);

      hits[i].fIsVisited = true;
      currentHit = &hits[i];

      uint returnTries(0);
      bool hasNext = GetNextHit(previousHit, currentHit, nextHit);

      while(hasNext) 
      {
        if(nextHit->fIsLeaf || nextHit->fIsVisited)
        {
          break;
        }
        track.push_back(nextHit);
        if(fvalidateTrackGrad && !ValidateGradTrack(track))
        {
            nextHit->fIsVisited=false;
            nextHit->fIsLeaf = true;
            i = 0; // reset the counter
            break;
        }

        currentHit->fIsVisited = true;
        previousHit = currentHit;
        currentHit = nextHit;


        hasNext = GetNextHit(previousHit, currentHit, nextHit);
        // If we fail, check if all neightbour cubes are visited.
        // if yes, then we may have gotten inside a loop.
        // try to return 2 cubes back and check if we can go around
        // the loop by taking another not visited cube.
        if( !hasNext   
              && RETURN_TRY_LIMIT>returnTries
              && track.size()>=RETURN_PREVIOUS
              && CheckAllVisited(currentHit))
        {
          ++returnTries;
          currentHit->fIsVisited = true;
          currentHit = track[track.size()-2];
          previousHit = track[track.size()-3];
          hasNext = GetNextHit(previousHit, currentHit, nextHit);
        }
      }

      tracks.push_back(track);
    }
  }
}

void FgdReconTemplate::CalculateGrad(std::vector<std::vector<ReconHit*>>& tracks, bool useSmoothCoor)
{
  const Int_t distToCalc = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_GRAD_DIST);
  const Int_t intervalToCal = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_GRAD_INTERVAL_DIST);

  static const Double_t radToDeg = 180/TMath::Pi();
  static const TVector3 xAxisVec(1,0,0);
  static const TVector3 yAxisVec(0,1,0);
  static const TVector3 zAxisVec(0,0,1);

  for(Int_t i=0; i<tracks.size(); ++i)
  {
    std::vector<ReconHit*>& track = tracks[i];
    LOG(debug2) << "=================================== ";
    LOG(debug2) << "Track " << i;

    for(Int_t j=0; j<track.size(); ++j)
    {

      ReconHit* currentHit = track[j];

      // 1. Calculate direction from the previous hit
      TVector3 diffVec(0,0,0);
      if(j>=1)
      {
        ReconHit* previous = track[j-1];
        diffVec = useSmoothCoor ? 
                  currentHit->fsmoothco - previous->fsmoothco:
                  currentHit->fmppcLoc - previous->fmppcLoc;
      }

      // 2.Calculate cosine change beween 2 consecutive vectors - gradient
      Double_t diffAngle(0);    // Difference angle between two vectors - the vectors are the difference in position (grad (position))
                                // between two hits distToCalc cubes apart. The angles measures the angle change of the track
                                // from hit to hit. distToCalc distance is selected to smooth out the change in angle.

      Double_t xAxisAngle(0);
      Double_t yAxisAngle(0);
      Double_t zAxisAngle(0);   // Vector angle between the z axis and the vector difference between the current hit and the hit distToCalc 
                                // poisitions back.
      

      TVector3 diffVec1(0,0,0);
      TVector3 diffVec2(0,0,0);

      Int_t indOne = j - distToCalc + 1;
      Int_t indTwo = j; // current hit

      if(j>= (distToCalc -1) )
      {
        ReconHit* one = track[indOne];
        ReconHit* two = track[indTwo]; 

        diffVec1 = useSmoothCoor ? 
                  two->fsmoothco - one->fsmoothco:
                  two->fmppcLoc - one->fmppcLoc;
        
        xAxisAngle = radToDeg * xAxisVec.Angle(diffVec1);
        yAxisAngle = radToDeg * yAxisVec.Angle(diffVec1);
        zAxisAngle = radToDeg * zAxisVec.Angle(diffVec1);

        // The z angle is calculated for the current hit
        currentHit->fXaxisAngle = xAxisAngle;
        currentHit->fYaxisAngle = yAxisAngle;
        currentHit->fZaxisAngle = zAxisAngle;
      }
      
      if(j>= (distToCalc + intervalToCal -1) )
      {
        Int_t indOneP = indOne - intervalToCal;
        Int_t indTwoP = indTwo - intervalToCal;
        ReconHit* oneP = track[indOneP];
        ReconHit* twoP = track[indTwoP];
        diffVec2 = useSmoothCoor ? 
                  twoP->fsmoothco - oneP->fsmoothco:
                  twoP->fmppcLoc - oneP->fmppcLoc;

        diffAngle = radToDeg * diffVec1.Angle(diffVec2);

        // The angle is calculated for the current hit
        currentHit->fChangeAngle = diffAngle;
      }

      LOG(debug2) << "Id " << currentHit->fLocalId
                  << " \tX " << currentHit->fmppcLoc.X()
                  << " \tY " << currentHit->fmppcLoc.Y()
                  << " \tZ " << currentHit->fmppcLoc.Z() 
                  << " \tPhotons  " << "(" << currentHit->fphotons.X() << "," << currentHit->fphotons.Y() << "," << currentHit->fphotons.Z() << ")"
                  << " \tChange  " << "(" << diffVec.X() << "," << diffVec.Y() << "," << diffVec.Z() << ")"
                  // << " \t(" << diffVec1.X() << "," << diffVec1.Y() << "," << diffVec1.Z() << ") "
                  // << " \t(" << diffVec2.X() << "," << diffVec2.Y() << "," << diffVec2.Z() << ") "
                  << " \tAngle  " << diffAngle
                  << " \tX axis Angle  " << xAxisAngle
                  << " \tY axis Angle  " << yAxisAngle
                  << " \tZ axis Angle  " << zAxisAngle;
    }

    LOG(debug2) << "=================================== ";
  }
}

Bool_t FgdReconTemplate::ValidateGradTrack(std::vector<ReconHit*>& track)
{
    static const Double_t radToDeg = 180/TMath::Pi();

    Bool_t rc(true);
    for(size_t i = 0; rc && i < fgraphRecords.size(); ++i)
    {
        FgdGraphRecord& graphRec = fgraphRecords[i];
        Int_t minLength = graphRec.fvecDist + graphRec.fdisDiff;

        if(minLength < track.size())
        {
            int lastHitInd = track.size() - 1;
            ReconHit* oneCurrent = track[lastHitInd];
            ReconHit* onePrevious = track[lastHitInd - graphRec.fvecDist];

            TVector3 one = oneCurrent->fmppcLoc - onePrevious->fmppcLoc;

            ReconHit* twoCurrentHit = track[lastHitInd - graphRec.fdisDiff];
            ReconHit* twoPrevious = track[lastHitInd - graphRec.fdisDiff - graphRec.fvecDist];

            TVector3 two = twoCurrentHit->fmppcLoc - twoPrevious->fmppcLoc;

            Double_t diffAngle = radToDeg * one.Angle(two);

            rc = (diffAngle <= graphRec.fgradDiffAlowed);
        }
    }

    return rc;
}

void FgdReconTemplate::SplitTrack(std::vector<std::vector<ReconHit*>>& originalTracks, std::vector<std::vector<ReconHit*>>& splitTracks)
{
  const double gradAllowedDiff = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::FGD_GRAD_ALLOWABLE_DIFF);
  const int previousHitToCompare = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_GRAD_INTERVAL_DIST);
  
  LOG(debug) << "Initial tracks size " << originalTracks.size();

  for(Int_t i=0; i<originalTracks.size(); ++i)
  {
    std::vector<ReconHit*>& track = originalTracks[i];

    if(track.empty())
    {
      continue;
    }

    std::vector<ReconHit*> trackToAdd;
    trackToAdd.push_back(track[0]); // Add the first hit to start the track

    for(size_t j=previousHitToCompare; j<track.size(); ++j)
    {
      ReconHit* currentHit = track[j];
      ReconHit* previousHit = track[j-previousHitToCompare];

      if(currentHit->fChangeAngle == 0)
      {
        // If there is no change add the track
        // Also '0' indicates that the track hit is too short to calculate the 
        // the change.
        // Only compare if both the current and previous angle is not zero.
        trackToAdd.push_back(currentHit);
      }
      else 
      {
        Double_t differenceInAngle = currentHit->fZaxisAngle - previousHit->fZaxisAngle;
        Double_t absDiff = std::fabs(differenceInAngle);
        if(absDiff >= gradAllowedDiff)
        {
          // If the difference is not allowed add the current track hit points
          // and reset it to start a new track.
          splitTracks.emplace_back(trackToAdd);
          trackToAdd.clear();
        }
        trackToAdd.push_back(currentHit);
      }
    }

    splitTracks.emplace_back(trackToAdd);
  }

  LOG(debug) << "Split tracks size " << splitTracks.size();
}


Bool_t FgdReconTemplate::FindUsingHough(std::vector<ReconHit>& hits
                                , std::vector<std::vector<ReconHit*>>& foundTracks
                                , std::vector<std::vector<TVector3>>& foundTracksPoints
                                , HoughType houghType
                                , TVector3 trackVertex
                                , bool useVertex)
{
    double maxdistxy = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_MAXDISTXY);
    double maxdistsz = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_MAXDISTSZ);
    double maxdistxyfit = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_MAXDISTXYFIT);
    double maxdistszfit = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_MAXDISTSZFIT);
    unsigned int minhitnumber = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_MINHITNUMBER);
    unsigned int xythetabins = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_XYTHETABINS);
    unsigned int xyd0bins = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_XYD0BINS);
    unsigned int xyomegabins = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_XYOMEGABINS);
    unsigned int szthetabins = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_SZTHETABINS);
    unsigned int szd0bins = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_SZD0BINS);
    double maxdxy = f_total_X + f_total_Y;
    double maxdsz = f_total_Z;
    unsigned int searchneighborhood = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::PF_SEACHINTHENEIGHBORHOOD);

    pathfinder::FinderParameter* newFinderParameter = nullptr;
    switch(houghType)
    {
        case HoughType::HELIX:
            newFinderParameter= new pathfinder::FinderParameter(false, true); 
            newFinderParameter -> setFindCurler(false);
            break;
        case HoughType::CURL:
            newFinderParameter= new pathfinder::FinderParameter(false, true); 
            newFinderParameter -> setFindCurler(true);
            break;
        case HoughType::STRAIGHT_LINE:
        default:
            newFinderParameter= new pathfinder::FinderParameter(true, false); 
            newFinderParameter -> setFindCurler(false);
            break;
    }


    if(useVertex)
    {
        std::pair<double, double> vertex(trackVertex.X(), trackVertex.Y());
        newFinderParameter -> setVertex(vertex);
    }

    newFinderParameter -> setMaxXYDistance(maxdistxy);
    newFinderParameter -> setMaxSZDistance(maxdistsz);
    newFinderParameter -> setMaxXYDistanceFit(maxdistxyfit);
    newFinderParameter -> setMaxSZDistanceFit(maxdistszfit);
    newFinderParameter -> setMinimumHitNumber(minhitnumber);
    newFinderParameter -> setNumberXYThetaBins(xythetabins);
    newFinderParameter -> setNumberXYDzeroBins(xyd0bins);
    newFinderParameter -> setNumberXYOmegaBins(xyomegabins);
    newFinderParameter -> setNumberSZThetaBins(szthetabins);
    newFinderParameter -> setNumberSZDzeroBins(szd0bins);
    newFinderParameter -> setMaxDxy(maxdxy);
    newFinderParameter -> setMaxDsz(maxdsz);
    
    if(searchneighborhood == 0)
    {
        newFinderParameter -> setSearchNeighborhood(false);
    }
    else
    {
        newFinderParameter -> setSearchNeighborhood(true);
    }
    newFinderParameter -> setSaveRootFile(false);

    std::vector<pathfinder::basicHit> allHits;
    for(size_t i=0; i< hits.size(); ++i)
    {
        TVector3 pointLoc = hits[i].fHitPos;
        allHits.emplace_back(pathfinder::basicHit(  pointLoc.X()
                                                    , pointLoc.Y()
                                                    , pointLoc.Z()
                                                    )
                                );
    }
    
    pathfinder::HoughTrafoTrackFinder newTrackFinder;

    //setting steering parameter
    newTrackFinder.setFinderParameter(*newFinderParameter);

    // If there is not time interval include all hits
    newTrackFinder.setInitialHits(allHits);

    //do the actual track finding
    Bool_t found = newTrackFinder.find();
    if(found)
    {
        ReconHit tempCopyHit;
        std::vector<pathfinder::TrackFinderTrack> pfTracks = newTrackFinder.getTracks();
        for(Int_t i =0; i <  pfTracks.size() ; i++)
        {
        LOG(debug2) << "Track " << i;
        pathfinder::TrackFinderTrack& trFinder = pfTracks[i];
        const std::vector<pathfinder::basicHit>& hitsOnTrack = trFinder.getHitsOnTrack();
        std::vector<ReconHit*> track;
        std::vector<TVector3> trackPoints;
        for(size_t j =0; j< hitsOnTrack.size(); ++j)
        {
            trackPoints.emplace_back(hitsOnTrack[j].getX(),hitsOnTrack[j].getY(),hitsOnTrack[j].getZ());

            tempCopyHit.fHitPos.SetXYZ(hitsOnTrack[j].getX(),hitsOnTrack[j].getY(),hitsOnTrack[j].getZ());
            auto it = std::find(hits.begin(), hits.end(), tempCopyHit);
            if(it != hits.end())
            {
                track.push_back(&(*it));
                ReconHit& argHit = (*it);
                LOG(debug2) << "X " << argHit.fmppcLoc.X() << " Y " << argHit.fmppcLoc.Y()  << " Z " << argHit.fmppcLoc.Z();
            }

            LOG(debug2) << "X " << hitsOnTrack[j].getX() << " Y " << hitsOnTrack[j].getY() << " Z " << hitsOnTrack[j].getZ();
        }
        foundTracks.push_back(track);
        foundTracksPoints.push_back(trackPoints);
        LOG(debug2) << "============================ ";
        }
    }

    return found;
}

Long_t FgdReconTemplate::hitId(ReconHit& hit)
{
    return ArrInd(hit.fmppcLoc.X(), hit.fmppcLoc.Y(), hit.fmppcLoc.Z());
}

Long_t FgdReconTemplate::ArrInd(int x, int y, int z)
{
  return (x*f_bin_Y*f_bin_Z + y*f_bin_Z+z);
}

void FgdReconTemplate::initGraphRecord(const char* graphTrackConfig)
{
    std::ifstream graphFileStream;
    const char spaceChar(' ');
    std::vector<std::string> values;
    try
    {        
        graphFileStream.open(graphTrackConfig, std::ios::in);

        if(graphFileStream.is_open())
        {
            values.clear();
            std::string line;
            while(std::getline(graphFileStream,line))
            {
                if(line.empty() || line[0]=='#') continue;

                std::istringstream ss(line);
                std::string token;
                while(std::getline(ss, token, spaceChar))
                {
                    if(!token.empty())
                    {
                        values.push_back(token);
                    }
                }
                Int_t segLenght = std::stoi(values[Data::SEGMENT_LENGTH]);
                Int_t distToSegment = std::stoi(values[Data::DISTANCE_TO_SEGMENT]);
                Double_t gradAllowed = std::stod(values[Data::GRAD_ALLOWED_DEGREES]);
                fgraphRecords.push_back(FgdGraphRecord(segLenght,distToSegment,gradAllowed));
                values.clear();
            }
        }
    }
    catch(const std::exception& e)
    {
        LOG(fatal) << e.what();
    }
}

bool FgdReconTemplate::FitTrack(
       std::vector<ReconHit>& track
    ,  Bool_t useRefFitter
    ,  Int_t pdg
    ,  TVector3& posM
    ,  TVector3& momM
    ,  Bool_t smearMomentum
    ,  Bool_t useSmoothPos
    ,  Bool_t usePhotonResolution
    ,  TVector3& momentum)
{

    std::shared_ptr<genfit::AbsKalmanFitter> fitter{nullptr};
    if(useRefFitter)
    {
        fitter = make_shared<genfit::KalmanFitterRefTrack>();
    }
    else
    {
        fitter = make_shared<genfit::KalmanFitter>();
    }
    fitter->setMinIterations(fminGenFitInterations);
    fitter->setMaxIterations(fmaxGenFitIterations);
    //fitter->setDeltaPval(1e-7);
    fitter->setDebugLvl(fDebuglvl_genfit);


    TVector3 magField = fgdConstructor.GetMagneticField(); // values are in kGauss

    static int detId(1); // Detector id, it is the same, we only have one detector
    Bool_t rc(false);

    // genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(false);
    // genfit::MaterialEffects::getInstance()->drawdEdx(pdg);
    // genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(true);
    // genfit::MaterialEffects::getInstance()->setNoiseCoulomb(true);
    // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
    // genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
    // genfit::MaterialEffects::getInstance()->ignoreBoundariesBetweenEqualMaterials(true);

    // approximate covariance
    static const double poisson_resolution = 0.1;  // Default in example is 0.1;
                                                               // cm; resolution of generated measurements

    Int_t maxPhotonsPerCube = 1;
    for(Int_t ph = 0; ph < track.size(); ++ph)
    {
        Int_t sumPhotons = track[ph].fphotons.X() + track[ph].fphotons.Y() + track[ph].fphotons.Z();
        if(maxPhotonsPerCube < sumPhotons) {maxPhotonsPerCube = sumPhotons;}
    }

    TMatrixDSym hitCov(3);
    hitCov(0,0) = poisson_resolution*poisson_resolution;
    hitCov(1,1) = poisson_resolution*poisson_resolution;
    hitCov(2,2) = poisson_resolution*poisson_resolution;

    TMatrixDSym covM(6);
    for (int ci = 0; ci < 3; ++ci)
        covM(ci,ci) = poisson_resolution*poisson_resolution;
    for (int ci = 3; ci < 6; ++ci)
        covM(ci,ci) = covM(ci,ci) = pow(  ((poisson_resolution / track.size()) / sqrt(3)), 2); 

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);

    stateSmeared.setPosMomCov(posM, momM, covM);
    //stateSmeared.setPosMomCov(posM, calMom, covM);

    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);

    genfit::Track* toFitTrack = nullptr;
    if(smearMomentum) { toFitTrack = new genfit::Track(rep, seedState, seedCov);}
    else {toFitTrack = new genfit::Track(rep, posM, momM); }
    
    // genfit::Track* toFitTrack = new genfit::Track();
    // toFitTrack->setStateSeed(seedState);
    // toFitTrack->setCovSeed(seedCov);
    // toFitTrack->addTrackRep(rep);
    
    //std::vector<genfit::AbsMeasurement*> measurements;
    for(Int_t bh = 0; bh < track.size(); ++bh)
    {
        static const double cube_diag = f_step_X*f_step_X + f_step_Y*f_step_Y + f_step_Z*f_step_Z;
        static const double res_length = (sqrt(cube_diag)/2);
        TVectorD hitPos(3);
        hitPos(0) = useSmoothPos? track[bh].fsmoothco.X() : track[bh].fHitPos.X();
        hitPos(1) = useSmoothPos? track[bh].fsmoothco.Y() : track[bh].fHitPos.Y();
        hitPos(2) = useSmoothPos? track[bh].fsmoothco.Z() : track[bh].fHitPos.Z();

        // Use monte carlo coordiantes
        // hitPos(0) = hitsOnTrack[bh].fMCPos.X();
        // hitPos(1) = hitsOnTrack[bh].fMCPos.Y();
        // hitPos(2) = hitsOnTrack[bh].fMCPos.Z();

        Int_t sumPhotons = track[bh].fphotons.X() + track[bh].fphotons.Y() + track[bh].fphotons.Z();
        double normCoeff = res_length * (1 - sumPhotons/maxPhotonsPerCube);     //std::sqrt(sumPhotons);
        normCoeff = normCoeff == 0 ? 1 : normCoeff;
        TMatrixDSym phhitCov(3);
        // The more photons, the more likely the particle went through the center
        phhitCov(0,0) = poisson_resolution * poisson_resolution * normCoeff;
        phhitCov(1,1) = poisson_resolution * poisson_resolution * normCoeff;
        phhitCov(2,2) = poisson_resolution * poisson_resolution * normCoeff;

        genfit::AbsMeasurement* measurement = nullptr;
        if(usePhotonResolution) {measurement = new genfit::SpacepointMeasurement(hitPos, hitCov, detId, bh, nullptr); }
        else {measurement = new genfit::SpacepointMeasurement(hitPos, phhitCov, detId, bh, nullptr); }

        toFitTrack->insertMeasurement(measurement, bh);
        //measurements.emplace_back(measurement);
    }
    //toFitTrack->insertPoint(new genfit::TrackPoint(measurements, toFitTrack));

    try
    {
      //check
      toFitTrack->checkConsistency();

      // do the fit
      LOG(info) << "Somthing ";
      fitter->processTrack(toFitTrack, true);
      toFitTrack->checkConsistency();

      genfit::FitStatus* fiStatuStatus = toFitTrack->getFitStatus(rep);
      rc = fiStatuStatus->isFitted() && (fiStatuStatus->isFitConverged() || fiStatuStatus->isFitConvergedPartially());
      
      const genfit::MeasuredStateOnPlane& me = (*toFitTrack).getFittedState();
      momentum = me.getMom();
    }
    catch(genfit::Exception& e)
    {
        LOG(error) <<"Exception, when tryng to fit track";
        LOG(error) << "e.what() " << e.what();
        LOG(error) << "e.getExcString() " << e.getExcString();
    }

    return rc;
}

} //superfgd
} //reconstruction
} //esbroot