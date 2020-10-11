#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_RECON_TEMPLATE_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_RECON_TEMPLATE_H

// Essnusb
#include "EsbReconstruction/EsbSuperFGD/FgdReconHit.h"
#include "EsbGeometry/EsbSuperFGD/EsbFgdDetectorParameters.h"
#include "EsbReconstruction/EsbSuperFGD/FgdGraphRecord.h"

// Root
#include "TObject.h"
#include <TVector3.h>

// PathFinder headers
#include "FinderParameter.h"
#include "HoughTrafoTrackFinder.h"
#include "TrackParameterFull.h"

#include <algorithm>
#include <vector>
#include <set>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

#define RETURN_TRY_LIMIT 2
#define RETURN_CURRENT 2
#define RETURN_PREVIOUS 3

class FgdReconTemplate : public TObject
{

public:

    enum HoughType : size_t
    {
        STRAIGHT_LINE,
        HELIX,
        CURL
    };

    /** Default constructor **/  
    FgdReconTemplate();

    /** Constructor with argument
        *@param geoConfigFile  - Configuration file detector
        *@param graphTrackConfig -  Configuration for the graph traversal algorithm
        *@param validateTrackGrad -  validate the track using the grahConfig file
    **/  
    FgdReconTemplate(const char* geoConfigFile, const char* graphTrackConfig, bool validateTrackGrad = false);

    ~FgdReconTemplate();

    Bool_t IsLeaf(ReconHit* hit);

    Bool_t GetNextHit(ReconHit* previous, ReconHit* current, ReconHit*& next);
    Bool_t CheckAllVisited(ReconHit* hit);

    void BuildGraph(std::vector<ReconHit>& hits);
    void FindLeaves(std::vector<ReconHit>& hits);
    void SmoothGraph(std::vector<ReconHit>& hits);
    
    void FindTracks(std::vector<ReconHit>& hits, std::vector<std::vector<ReconHit*>>& tracks);
    void CalculateGrad(std::vector<std::vector<ReconHit*>>& tracks, bool useSmoothCoor = false);
    Bool_t ValidateGradTrack(std::vector<ReconHit*>& track);
    void SplitTrack(std::vector<std::vector<ReconHit*>>& originalTracks, std::vector<std::vector<ReconHit*>>& splitTracks);


    Bool_t FindUsingHough(std::vector<ReconHit>& hits
                                , std::vector<std::vector<ReconHit*>>& foundTracks
                                , std::vector<std::vector<TVector3>>& foundTracksPoints
                                , HoughType houghType
                                , TVector3 trackVertex
                                , bool useVertex);

private:

    /** Class Containing the vectors for each template found**/ 
    class HitTemplate
    {
    public:
        HitTemplate(){}

        ~HitTemplate(){}

        HitTemplate(const HitTemplate& c)
        {
            this->previousHit = c.previousHit;
            this->nextHit = c.nextHit;
            this->hitVectors = c.hitVectors;
        }

        HitTemplate& operator=(const HitTemplate& c)
        {
            this->previousHit = c.previousHit;
            this->nextHit = c.nextHit;
            this->hitVectors = c.hitVectors;
            return *this;
        }

        size_t Length(){return hitVectors.size();}
    
        TVector3 nextHit;
        TVector3 previousHit;
        std::vector<TVector3> hitVectors;
    };

    enum Data : int   // Position in the tokens which value represents which data 
    {
        SEGMENT_LENGTH = 0,
        DISTANCE_TO_SEGMENT = 1,
        GRAD_ALLOWED_DEGREES = 2
    };

    void LoadTemplates();
    void GetHitVectors(ReconHit* hit, std::vector<TVector3>& vecs);

    Bool_t AreVectorsEqual(const std::vector<TVector3>& tempVecs, const std::vector<TVector3>& vecs, Int_t& foundPermutation );
    TVector3 GetPermutation(TVector3 vec, Int_t numPermutation);
    void SmoothCoordinate(ReconHit* hit, TVector3& cord, std::set<Long_t>& visited, size_t depth = 1);
    Long_t hitId(ReconHit& hit);
    Long_t ArrInd(int i, int j, int k);
    void initGraphRecord(const char* graphTrackConfig);

    std::vector<FgdReconTemplate::HitTemplate> fLeafVectors;//!<!  
    std::vector<FgdGraphRecord> fgraphRecords;//!<! 

    /** Class to hold the Detector parameters read from external file **/
    esbroot::geometry::superfgd::FgdDetectorParameters fParams;//!<! 

    /** Detector dimentions **/
    Double_t flunit;
    double f_step_X;
    double f_step_Y;
    double f_step_Z;
    int f_bin_X;
    int f_bin_Y;
    int f_bin_Z;
    double f_total_X;
    double f_total_Y;
    double f_total_Z;

    int fsmoothDepth;
    Double_t fsmoothErrLimit;
    bool fvalidateTrackGrad;

    ClassDef(FgdReconTemplate, 2);
};

} //superfgd
} //reconstruction
} //esbroot

#endif