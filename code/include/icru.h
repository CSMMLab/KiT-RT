#ifndef ICRU_H
#define ICRU_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "common/globalconstants.h"
#include "common/typedef.h"
#include "interpolation.h"
#include "toolboxes/errormessages.h"

enum ParticleType { ELECTRON, POSITRON };

class ICRU
{
  private:
    const unsigned materialID   = 278;
    const ParticleType particle = ParticleType::ELECTRON;
    const double PI             = 3.1415926535897932e0;
    const double R4PI           = 1.0e0 / ( 4.0e0 * PI );
    // const double FOURPI         = 4.0e0 * PI;
    // const double AVOG           = 6.02214199e23;
    const double EMC2  = 510.998902e3;
    const double A0B   = 5.291772083e-9;
    const double HREV  = 27.2113834e0;
    const double TEMC2 = EMC2 + EMC2;
    const double CONS  = 2.0e0 * PI * ( HREV * A0B ) * ( HREV * A0B ) / EMC2;
    const std::vector<double> EGRD{
        1.00e0, 1.25e0, 1.50e0, 1.75e0, 2.00e0, 2.50e0, 3.00e0, 3.50e0, 4.00e0, 4.50e0, 5.00e0, 6.00e0, 7.00e0, 8.00e0, 9.00e0, 1.00e1 };
    const std::vector<double> ELAW{
        1.007900e0, 4.002600e0, 6.941000e0, 9.012200e0, 1.081100e1, 1.201070e1, 1.400670e1, 1.599940e1, 1.899840e1, 2.017970e1, 2.298980e1,
        2.430500e1, 2.698150e1, 2.808550e1, 3.097380e1, 3.206600e1, 3.545270e1, 3.994800e1, 3.909830e1, 4.007800e1, 4.495590e1, 4.786700e1,
        5.094150e1, 5.199610e1, 5.493800e1, 5.584500e1, 5.893320e1, 5.869340e1, 6.354600e1, 6.539000e1, 6.972300e1, 7.261000e1, 7.492160e1,
        7.896000e1, 7.990400e1, 8.380000e1, 8.546780e1, 8.762000e1, 8.890590e1, 9.122400e1, 9.290640e1, 9.594000e1, 9.890630e1, 1.010700e2,
        1.029055e2, 1.064200e2, 1.078682e2, 1.124110e2, 1.148180e2, 1.187100e2, 1.217600e2, 1.276000e2, 1.269045e2, 1.312900e2, 1.329054e2,
        1.373270e2, 1.389055e2, 1.401160e2, 1.409076e2, 1.442400e2, 1.449127e2, 1.503600e2, 1.519640e2, 1.572500e2, 1.589253e2, 1.625000e2,
        1.649303e2, 1.672600e2, 1.689342e2, 1.730400e2, 1.749670e2, 1.784900e2, 1.809479e2, 1.838400e2, 1.862070e2, 1.902300e2, 1.922170e2,
        1.950780e2, 1.969666e2, 2.005900e2, 2.043833e2, 2.072000e2, 2.089804e2, 2.089824e2, 2.099871e2, 2.220176e2, 2.230197e2, 2.260254e2,
        2.270277e2, 2.320381e2, 2.310359e2, 2.380289e2, 2.370482e2, 2.440642e2, 2.430614e2, 2.470000e2, 2.470000e2, 2.510000e2, 2.520000e2,
        2.570000e2, 2.580000e2, 2.590000e2, 2.620000e2 };

    unsigned _NE;
    unsigned _NA;
    unsigned _NELEM;
    unsigned _NOSC;

    double _CSI;
    double _TCS1I, _TCS2I;
    double _RMU1, _RMU2;
    double _B1, _B2, _B3, _B4;
    double _C0, _C1, _C2, _C3, _C4;
    double _BETA2;
    double _R1;

    Vector _E, _QMU;
    std::vector<double> _ET, _ETL;
    std::vector<double> _XMU, _XMUL;
    std::vector<double> _ECS, _PCS;
    std::vector<double> _ETCS1, _ETCS2;
    std::vector<double> _PTCS1, _PTCS2;
    std::vector<double> _DCSI, _DCSIL;
    std::vector<double> _STF;
    std::vector<double> _IZ;
    std::vector<double> _F, _WRI;

    Matrix _EDCS, _PDCS;

    void ELINIT( std::vector<double> IZ, std::vector<double> STF );
    void DCSEL0( double E );
    void SPLINE( const std::vector<double>& X,
                 const std::vector<double>& Y,
                 std::vector<double>& A,
                 std::vector<double>& B,
                 std::vector<double>& C,
                 std::vector<double>& D,
                 double S1,
                 double SN,
                 unsigned N );
    void FINDI( std::vector<double> X, double XC, unsigned N, unsigned& I );
    double DCSIN3( double E, double RMU );
    double DCSEL( double RMU );
    double DODCSD( double RMU );
    double DODCSC( double RMU, double EL );
    void DOICSS( double WR, double E, double RMU, double& DCSD, double& DCSC );
    void angdcs( unsigned ncor, double e, std::vector<double>& dxse, std::vector<double>& dxsi, std::vector<double>& dxs );

    ICRU() = delete;

  public:
    ICRU( const Vector& mu, const Vector& energy ) : _NE( 96 ), _NA( 606 ), _E( 1e6 * energy ), _QMU( mu ) {
        std::vector<double> TH( _NA ), THR( _NA );
        _XMU.resize( _NA );
        for( unsigned i = 0; i < _NA; ++i ) {
            if( i == 0 ) {
                TH[i]   = 0.0e0;
                THR[i]  = TH[i] * M_PI / 180.0e0;
                _XMU[i] = ( 1.0e0 - std::cos( THR[i] ) ) / 2.0e0;
            }
            else if( i == 1 ) {
                TH[i]   = 1.0e-4;
                THR[i]  = TH[i] * M_PI / 180.0e0;
                _XMU[i] = ( 1.0e0 - std::cos( THR[i] ) ) / 2.0e0;
            }
            else {
                if( TH[i - 1] < 0.9999e-3 )
                    TH[i] = TH[i - 1] + 2.5e-5;
                else if( TH[i - 1] < 0.9999e-2 )
                    TH[i] = TH[i - 1] + 2.5e-4;
                else if( TH[i - 1] < 0.9999e-1 )
                    TH[i] = TH[i - 1] + 2.5e-3;
                else if( TH[i - 1] < 0.9999e+0 )
                    TH[i] = TH[i - 1] + 2.5e-2;
                else if( TH[i - 1] < 0.9999e+1 )
                    TH[i] = TH[i - 1] + 1.0e-1;
                else if( TH[i - 1] < 2.4999e+1 )
                    TH[i] = TH[i - 1] + 2.5e-1;
                else
                    TH[i] = TH[i - 1] + 5.0e-1;

                THR[i]  = TH[i] * M_PI / 180.0e0;
                _XMU[i] = std::max( ( 1.0e0 - std::cos( THR[i] ) ) / 2.0e0, 1.0e-35 );
            }
        }

        _ET.resize( _NE );
        _ETL.resize( _NE );
        _XMUL.resize( _NA );
        _ECS.resize( _NE );
        _ETCS1.resize( _NE );
        _ETCS2.resize( _NE );
        _PCS.resize( _NE );
        _PTCS1.resize( _NE );
        _PTCS2.resize( _NE );
        _DCSI.resize( _NA );
        _DCSIL.resize( _NA );

        _EDCS.resize( _NE, _NA );
        _PDCS.resize( _NE, _NA );
    }
    ~ICRU(){};

    void GetAngularScatteringXS( Matrix& angularXS, Vector& integratedXS );
    void GetStoppingPower( Vector& stoppingPower );
};

#endif    // ICRU_H
