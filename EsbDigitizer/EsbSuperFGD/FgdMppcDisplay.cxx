#include "EsbDigitizer/EsbSuperFGD/FgdMppcDisplay.h"
#include "EsbData/EsbSuperFGD/FgdDetectorPoint.h"

#include <TClonesArray.h>
#include <TSpline.h>
#include <TRandom.h>
#include <TH2F.h>
#include <TH1F.h>

#include <FairRootManager.h>
#include "FairLogger.h"

#include "EsbData/WCDetectorPoint.h"
#include "EsbGeometry/PMTube.h"
#include "EsbData/PMTubeHit.h"

#include <iostream>
#include <sstream>
using std::cout;
using std::endl;

namespace esbroot {
namespace digitizer {
namespace superfgd {

// -----   Default constructor   -------------------------------------------
FgdMppcDisplay::FgdMppcDisplay() :
  FairTask(), fX(0), fY(0), fZ(0),fevNum(0),
  f_xy_hist(nullptr), f_yz_hist(nullptr), f_xz_hist(nullptr), fHitArray(nullptr)
  , f_photo_ave(0.),f_photo_count(0)
{ 
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdMppcDisplay::FgdMppcDisplay(const char* name
                          ,const char* geoConfigFile
                          ,double x, double y, double z
                          , Int_t verbose) :
  FairTask(name, verbose), fX(x), fY(y), fZ(z),fevNum(0),
  f_xy_hist(nullptr), f_yz_hist(nullptr), f_xz_hist(nullptr), fHitArray(nullptr)
  , f_photo_ave(0.), f_photo_count(0)
{ 
  fParams.LoadPartParams(geoConfigFile);
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdMppcDisplay::~FgdMppcDisplay() 
{
  if(fHitArray) {
    fHitArray->Delete();
    delete fHitArray;
  }
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------
InitStatus FgdMppcDisplay::Init() 
{   
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

  // Get RootManager
  FairRootManager* manager = FairRootManager::Instance();
  if ( !manager ) {
    LOG(error) << "-E- FgdDigitizer::Init: " << "FairRootManager not instantised!";
    return kFATAL;
  }

  fHitArray = (TClonesArray*) manager->GetObject(geometry::superfgd::DP::FGD_HIT.c_str());
  if (!fHitArray) 
  {
      LOG(fatal) << "Exec", "No fgd hits array";
      return kFATAL;
  }

  fcanvas = new TCanvas();
  f_xy_hist = new TH2F("hist_xy","XY histogram",f_bin_X,0,f_bin_X,f_bin_Y,0,f_bin_Y);
  f_yz_hist = new TH2F("hist_yz","YZ histogram",f_bin_Y,0,f_bin_Y,f_bin_Z,0,f_bin_Z);
  f_xz_hist = new TH2F("hist_zx","XZ histogram",f_bin_X,0,f_bin_X,f_bin_Z,0,f_bin_Z);

  f_ph_hist = new TH1F("hist_ph","Photons per cube histogram",50,0,300);
  f_ph_hist_total = new TH1F("hist_ph_total","Total Photons histogram",100,0,30000);
  f_time_hist  = new TH1F("hist_time","Total Photons histogram",100,0,2);

  return kSUCCESS;
}


// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------
void FgdMppcDisplay::FinishEvent()
{
  if(f_xy_hist)
  {
    f_xy_hist->Draw("colz");
    WriteCanvas("XY");
    f_xy_hist->Reset();
  } 

  if(f_yz_hist)
  {
    f_yz_hist->Draw("colz");
    WriteCanvas("YZ");
    f_yz_hist->Reset();
  } 

  if(f_xz_hist)
  {
    f_xz_hist->Draw("colz");
    WriteCanvas("XZ");
    f_xz_hist->Reset();
  } 

  if(f_ph_hist)
  {
    Double_t mean = f_ph_hist->GetMean();
    Double_t meanErr = f_ph_hist->GetMeanError();
    Double_t stdDev = f_ph_hist->GetStdDev();
    Double_t stdDevErr = f_ph_hist->GetStdDevError();

    Double_t entries = f_ph_hist->GetEntries();
    std::cout << "Mean " << mean << " StdDev " << stdDev << std::endl;
    std::cout << "Mean Err " << meanErr << " StdDev Err " << stdDevErr << std::endl;
    std::cout << "Emtries " << entries << std::endl;
    f_ph_hist->Draw("colz");
    WriteCanvas("Photons");
    /* 
    if(f_ph_hist)
    {
      std::cout << " Event " << fevNum << std::endl;
      for(size_t i =0; i < 300 ; ++i)
      {
        std::cout << i << " -> " << f_ph_hist->GetBinContent(i) << std::endl;
      }
    }
    */
    f_ph_hist->Reset();
  } 

  if(f_time_hist)
  {
    f_time_hist->Draw("colz");
    WriteCanvas("Total_time");
    f_time_hist->Reset();
  }

  if(fcanvas)
  {
    fcanvas->ResetDrawn();
  }
  fevNum++;
}

void FgdMppcDisplay::FinishTask()
{
  if(f_ph_hist_total)
  {
    f_ph_hist_total->Draw("colz");
    WriteCanvas("Total_ph");
    f_ph_hist_total->Reset();
  }
}


void FgdMppcDisplay::Exec(Option_t* opt) 
{
  const Int_t hits = fHitArray->GetEntries();
  float max(0);
  float sum(0);
  for(Int_t i =0; i < hits; i++)
  {
      data::superfgd::FgdHit* hit = (data::superfgd::FgdHit*)fHitArray->At(i);
      TVector3  photoE = std::move(hit->GetPhotoE()); //returns a temp obj
      TVector3  mppcLoc = std::move(hit->GetMppcLoc()); //returns a temp obj

      if(f_xy_hist) f_xy_hist->Fill(mppcLoc.X(), mppcLoc.Y(), photoE.Z());
      if(f_yz_hist) f_yz_hist->Fill(mppcLoc.Y(), mppcLoc.Z(), photoE.X());
      if(f_xz_hist) f_xz_hist->Fill(mppcLoc.X(), mppcLoc.Z(), photoE.Y());

      float temp = mppcLoc.X() + mppcLoc.Y() + mppcLoc.Z();
      max = max < temp? temp : max;
      sum += temp;
      //cout  << "Time " << hit->GetTime() << endl;
      f_ph_hist->Fill(temp);
      f_time_hist->Fill(hit->GetTime());
  }
  f_ph_hist_total->Fill(sum);
  cout  << "Max photons per cube are " << max << endl;
  cout  << "sum photons are " << sum << endl;
  cout  << endl;
  cout  << endl;
}
// -------------------------------------------------------------------------
  

void FgdMppcDisplay::WriteCanvas(string hist)
{
    fcanvas->Update();
    fcanvas->Draw();
    fcanvas->Write();
    std::ostringstream strb;
    strb<< "hist_" << hist << fevNum << ".gif";
    fcanvas->SaveAs((strb.str()).c_str());
    fcanvas->ResetDrawn();
}

}// namespace superfgd
}// namespace digitizer
}// namespace esbroot
