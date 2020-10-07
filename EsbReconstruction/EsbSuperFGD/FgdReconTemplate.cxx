#include "EsbReconstruction/EsbSuperFGD/FgdReconTemplate.h"

#include "FairLogger.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

using namespace std;

#define PERMUTATIONS 23

FgdReconTemplate::FgdReconTemplate() 
{
    LoadTemplates();
}

FgdReconTemplate::FgdReconTemplate(const char* geoConfigFile)
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

    fsmoothDepth = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_SMOOTH_GRAPH_DEPTH);
    fsmoothErrLimit = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::FGD_SMOOTH_GRAPH_ERR_LIMIT);
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

Long_t FgdReconTemplate::hitId(ReconHit& hit)
{
    return ArrInd(hit.fmppcLoc.X(), hit.fmppcLoc.Y(), hit.fmppcLoc.Z());
}

Long_t FgdReconTemplate::ArrInd(int x, int y, int z)
{
  return (x*f_bin_Y*f_bin_Z + y*f_bin_Z+z);
}

} //superfgd
} //reconstruction
} //esbroot