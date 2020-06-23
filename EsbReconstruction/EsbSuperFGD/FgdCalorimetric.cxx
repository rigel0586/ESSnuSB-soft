#include "EsbReconstruction/EsbSuperFGD/FgdCalorimetric.h"
#include <math.h>


double RevertScintiResponse(double edep, double trackLength, double charge, double pe)
{
  return pe/EdepToPhotConv_SuperFGD;
}

void RevertFiberResponse(double &numPhotons, double &time, double distance)
{
    double temp_np = numPhotons;
    numPhotons = RevertFiberAttenuation(temp_np, distance);
}

double RevertFiberAttenuation(double Nphot0,double x)
{
    double a=0.;        // long attenuation component fraction
    double d=0.;        // distance MPPC-scint outside the bar
    double LongAtt=0.;  // long attenuation length
    double ShortAtt=0.; // short attenuation length
    double Ldecay=0.;   // decay length
    double Lbar=0.;     // bar length

    a = LongCompFrac_FGD;
    d = DistMPPCscint_FGD; 
    LongAtt = LongAtt_FGD;
    ShortAtt = ShortAtt_FGD;  
    Ldecay= DecayLength_FGD;
    Lbar = Lbar_FGD;
  
    double Nphot = Nphot0;
    Nphot /= ( a*exp((-x-d)/LongAtt) + (1-a)*exp((-x-d)/ShortAtt) );

    return Nphot;
}

double RevertFiberTime(double &time, double x)
{
  double fiberTime = time - TransTimeInFiber*x;
  return fiberTime;
}

double RevertyMPPCResponse(double npe)
{
  static const double Inverse_MPPCEff_SuperFGD = 2.15; // From trial an error the best match is chosen for the revert value
  return (npe * Inverse_MPPCEff_SuperFGD);
}



double RevertEdep(esbroot::reconstruction::superfgd::ReconHit& hit)
{
    TVector3& pe_1_dir = hit.fph1;
    TVector3& dist_1 = hit.fmppc1;
    TVector3& pe_2_dir = hit.fph2;
    TVector3& dist_2 = hit.fmppc2;

    double& time = hit.ftime;
    double charge = 1.0; // not used, but is a parameter from legacy


    double pe_1_x = RevertyMPPCResponse(pe_1_dir.X());
    double pe_1_y = RevertyMPPCResponse(pe_1_dir.Y());
    double pe_1_Z = RevertyMPPCResponse(pe_1_dir.Z());

    double pe_2_x = RevertyMPPCResponse(pe_2_dir.X());
    double pe_2_y = RevertyMPPCResponse(pe_2_dir.Y());
    double pe_2_z = RevertyMPPCResponse(pe_2_dir.Z());

    RevertFiberResponse(pe_1_x, time, dist_1.X());
    RevertFiberResponse(pe_1_y, time, dist_1.Y());
    RevertFiberResponse(pe_1_Z, time, dist_1.Z());

    RevertFiberResponse(pe_2_x, time, dist_2.X());
    RevertFiberResponse(pe_2_y, time, dist_2.Y());
    RevertFiberResponse(pe_2_z, time, dist_2.Z());

    // Deposited energy hit.fEdep and tracklength hit.ftrackLength
    // are used to calculate the CBIRKS coefficients together with dedx (energy losses)
    // to be able to revert from the photons

    Double_t x_1 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_x);
    Double_t y_1 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_y);
    Double_t z_1 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_Z);

    Double_t x_2 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_x);
    Double_t y_2 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_y);
    Double_t z_2 = RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_z);

    Double_t x = x_1 + x_2;
    Double_t y = y_1 + y_2;
    Double_t z = z_1 + z_2;

    Double_t totalEdep =  x + y + z;
    return totalEdep;
}

double Calcl(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track, size_t i)
{
    size_t lastId = track.size() - 1;
    if(i == lastId)
    {
        return 0;
    }

    size_t nextId = i+1;
    double x = track[i].fmppcLoc.X() - track[nextId].fmppcLoc.X();
    double y = track[i].fmppcLoc.Y() - track[nextId].fmppcLoc.Y();
    double z = track[i].fmppcLoc.Z() - track[nextId].fmppcLoc.Z();

    double l = x*x + y*y + z*z;
    return std::sqrt(l);
}

double RevertToDeDx(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track)
{
    double l(0);
    double E(0);

    for(size_t i = 0; i < track.size(); ++i)
    {
        E += RevertEdep(const_cast<esbroot::reconstruction::superfgd::ReconHit&>(track[i]));
        l += Calcl(track, i);
    }

    return E/l;
}



double RevertToDeDxMC(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track)
{
    double l(0);
    double E(0);

    Double_t firstL = track[0].ftrackLengthOrigin;
    Double_t lastL = track[track.size() - 1].ftrackLengthOrigin;

    l = firstL < lastL ? lastL : firstL;

    for(size_t i = 0; i < track.size(); ++i)
    {
        E += RevertEdep(const_cast<esbroot::reconstruction::superfgd::ReconHit&>(track[i]));
    }

    return E/l;
}

double RevertEdepToDeDx(const std::vector<esbroot::reconstruction::superfgd::ReconHit>& track)
{
    double l(0);
    double E(0);

    Double_t firstL = track[0].ftrackLengthOrigin;
    Double_t lastL = track[track.size() - 1].ftrackLengthOrigin;

    l = firstL < lastL ? lastL : firstL;

    for(size_t i = 0; i < track.size(); ++i)
    {
        const esbroot::reconstruction::superfgd::ReconHit& hit = track[i];
        E += hit.fEdep;
    }

    return E/l;
}

double dedxToP(double dedx)
{
    double a = dedx < coeff_sep ? a_coeff : a_coeff_2;
    double b = dedx < coeff_sep ? b_coeff : b_coeff_2;
    return ( 
                        ( a /  dedx ) + b
                    );
}