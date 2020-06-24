#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_CALORIMETRIC_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_CALORIMETRIC_H

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "EsbReconstruction/EsbSuperFGD/FgdReconHit.h"

static const double EdepToPhotConv_FGD = 70.8 / CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
static const double EdepToPhotConv_SuperFGD = EdepToPhotConv_FGD * 1.3;
static const double DistMPPCscint_FGD = 41*CLHEP::mm;
static const double LongCompFrac_FGD = 0.816;
static const double LongAtt_FGD = 11926.*CLHEP::mm;
static const double ShortAtt_FGD = 312.*CLHEP::mm;
static const double DecayLength_FGD = 0.0858 / CLHEP::mm;
static const double Lbar_FGD = 1864.3 * CLHEP::mm;
static const double TransTimeInFiber = 1./28.;  //  1cm/2.8e10[cm/s] * 1e9 [ns]



static const double a_coeff = 1.52375;
static const double b_coeff = -0.0825768;

static const double coeff_sep = 5; // MeV / cm

static const double a_coeff_2 = 0.980461;
static const double b_coeff_2 = 0.00125502;

static const TVector3 x_axis(1,0,0);
static const TVector3 y_axis(0,1,0);
static const TVector3 z_axis(0,0,1);





double RevertScintiResponse(double edep, double trackLength, double charge, double pe);
double RevertFiberAttenuation(double Nphot0,double x);
double RevertFiberTime(double &time, double x);
double RevertyMPPCResponse(double npe);

void RevertFiberResponse(double &numPhotons, double &time, double distance);

double RevertToDeDx(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track);
double RevertToDeDxMC(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track);
double RevertEdepToDeDx(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track);

double RevertEdep(esbroot::reconstruction::superfgd::ReconHit& hit);
double Calcl(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track, size_t i);

double dedxToP(double dedx);

TVector3 getCalorimetricMomentum(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track);

#endif // ESBROOT_ESBDRECONSTRUCTION_FGD_CALORIMETRIC_H