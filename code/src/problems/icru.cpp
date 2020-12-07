#include "problems/icru.h"
#include "common/config.h"

ICRU::ICRU( const Vector& mu, const Vector& energy, Config* settings ) : _NE( 96 ), _NA( 606 ), _E( 1e6 * energy ), _QMU( mu ) {
    _settings = settings;

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

void ICRU::angdcs( unsigned ncor, double e, std::vector<double>& dxse, std::vector<double>& dxsi, std::vector<double>& dxs ) {

    unsigned NM = materialID;
    if( NM < 1 || NM > 280 ) {
        ErrorMessages::Error( "The material number " + std::to_string( NM ) + " is not defined.", CURRENT_FUNCTION );
    }

    std::string PATHMT = _settings->GetDataDir() + "material/";
    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << materialID;

    std::string FILE1 = PATHMT + "mater" + ss.str() + ".stp";

    std::ifstream f( FILE1 );
    if( !f.good() ) {
        ErrorMessages::FileNotFoundError( FILE1, CURRENT_FUNCTION );
    }

    std::string line;
    for( unsigned i = 0; i < 5; ++i ) std::getline( f, line );
    // double RHO = std::stod( line.substr( line.find_first_of( "=" ) + 1, 16 ) );
    std::getline( f, line );
    _NELEM     = std::stoul( line.substr( line.find_first_of( "=" ) + 1, 3 ) );
    double ATW = 0.0, ZTOT = 0.0;
    _IZ.resize( _NELEM );
    _STF.resize( _NELEM );
    for( unsigned IEL = 0; IEL < _NELEM; ++IEL ) {
        std::getline( f, line );
        _IZ[IEL]  = std::stoul( line.substr( line.find_first_of( "=" ) + 1, 3 ) );
        _STF[IEL] = std::stod( line.substr( line.find_last_of( "=" ) + 1, 3 ) );
        ZTOT += _STF[IEL] * _IZ[IEL];
        ATW += ELAW[_IZ[IEL]] * _STF[IEL];
    }
    std::getline( f, line );
    std::getline( f, line );
    double EXPOT = std::stod( line.substr( line.find_first_of( "=" ) + 1, 16 ) );
    std::getline( f, line );
    _NOSC = std::stoul( line.substr( line.find_first_of( "=" ) + 1, line.back() ) );
    std::getline( f, line );
    std::getline( f, line );
    std::getline( f, line );
    double EXPOTR = 0.0;
    _F.resize( _NOSC );
    _WRI.resize( _NOSC );
    double buffer;
    for( unsigned I = 0; I < _NOSC; ++I ) {
        f >> buffer >> _F[I] >> buffer >> _WRI[I];
        EXPOTR += _F[I] * std::log( _WRI[I] );
    }
    f.close();
    EXPOTR = std::exp( EXPOTR / ZTOT );
    if( std::fabs( EXPOT - EXPOTR ) > 1.0e-5 * EXPOT ) {
        ErrorMessages::Error( "Inconsistent oscillator data", CURRENT_FUNCTION );
    }

    unsigned NDAT = _NA;
    ELINIT( _IZ, _STF );
    DCSEL0( e );
    dxse.resize( NDAT );
    dxsi.resize( NDAT );
    dxs.resize( NDAT );

    for( unsigned IA = 0; IA < NDAT; ++IA ) {
        dxse[IA] = DCSEL( _XMU[IA] );
    }

    if( ncor == 1u || ncor == 2u ) {    // 1: Inelastic DCS = elastic DCS/Z | 2: Morse's incoherent scattering function.
        ErrorMessages::Error( "Unsupported scattering model", CURRENT_FUNCTION );
    }
    else if( ncor == 3u ) {    // Sternheimer -Liljequist GOS model.
        for( unsigned IA = 0; IA < NDAT; ++IA ) {
            dxsi[IA] = DCSIN3( e, _XMU[IA] );
        }
    }
    else {    // No inelastic scattering correction
        for( unsigned IA = 0; IA < NDAT; ++IA ) {
            dxsi[IA] = 0.0e0;
        }
    }

    for( unsigned IA = 0; IA < _NA; ++IA ) {
        dxs[IA] = dxse[IA] + dxsi[IA];
    }
}

void ICRU::angdcs( unsigned ncor, std::vector<std::vector<double>>& dxs ) {

    unsigned NM = materialID;
    if( NM < 1 || NM > 280 ) {
        ErrorMessages::Error( "The material number " + std::to_string( NM ) + " is not defined.", CURRENT_FUNCTION );
    }

    std::string PATHMT = _settings->GetDataDir() + "material/";
    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << materialID;

    std::string FILE1 = PATHMT + "mater" + ss.str() + ".stp";

    std::ifstream f( FILE1 );
    if( !f.good() ) {
        ErrorMessages::FileNotFoundError( FILE1, CURRENT_FUNCTION );
    }

    std::string line;
    for( unsigned i = 0; i < 5; ++i ) std::getline( f, line );
    // double RHO = std::stod( line.substr( line.find_first_of( "=" ) + 1, 16 ) );
    std::getline( f, line );
    _NELEM     = std::stoul( line.substr( line.find_first_of( "=" ) + 1, 3 ) );
    double ATW = 0.0, ZTOT = 0.0;
    _IZ.resize( _NELEM );
    _STF.resize( _NELEM );
    for( unsigned IEL = 0; IEL < _NELEM; ++IEL ) {
        std::getline( f, line );
        _IZ[IEL]  = std::stoul( line.substr( line.find_first_of( "=" ) + 1, 3 ) );
        _STF[IEL] = std::stod( line.substr( line.find_last_of( "=" ) + 1, 3 ) );
        ZTOT += _STF[IEL] * _IZ[IEL];
        ATW += ELAW[_IZ[IEL]] * _STF[IEL];
    }
    std::getline( f, line );
    std::getline( f, line );
    double EXPOT = std::stod( line.substr( line.find_first_of( "=" ) + 1, 16 ) );
    std::getline( f, line );
    _NOSC = std::stoul( line.substr( line.find_first_of( "=" ) + 1, line.back() ) );
    std::getline( f, line );
    std::getline( f, line );
    std::getline( f, line );
    double EXPOTR = 0.0;
    _F.resize( _NOSC );
    _WRI.resize( _NOSC );
    double buffer;
    for( unsigned I = 0; I < _NOSC; ++I ) {
        f >> buffer >> _F[I] >> buffer >> _WRI[I];
        EXPOTR += _F[I] * std::log( _WRI[I] );
    }
    f.close();
    EXPOTR = std::exp( EXPOTR / ZTOT );
    if( std::fabs( EXPOT - EXPOTR ) > 1.0e-5 * EXPOT ) {
        ErrorMessages::Error( "Inconsistent oscillator data", CURRENT_FUNCTION );
    }

    unsigned NDAT = _NA;
    ELINIT( _IZ, _STF );

    std::vector<double> dxse;
    std::vector<double> dxsi;
    dxs.resize( _E.size() );

    for( unsigned n = 0; n < _E.size(); ++n ) {
        DCSEL0( _E[n] );
        dxse.resize( NDAT );
        dxsi.resize( NDAT );
        dxs[n].resize( NDAT );

        for( unsigned IA = 0; IA < NDAT; ++IA ) {
            dxse[IA] = DCSEL( _XMU[IA] );
        }

        if( ncor == 1u || ncor == 2u ) {    // 1: Inelastic DCS = elastic DCS/Z | 2: Morse's incoherent scattering function.
            ErrorMessages::Error( "Unsupported scattering model", CURRENT_FUNCTION );
        }
        else if( ncor == 3u ) {    // Sternheimer -Liljequist GOS model.
            for( unsigned IA = 0; IA < NDAT; ++IA ) {
                dxsi[IA] = DCSIN3( _E[n], _XMU[IA] );
            }
        }
        else {    // No inelastic scattering correction
            for( unsigned IA = 0; IA < NDAT; ++IA ) {
                dxsi[IA] = 0.0e0;
            }
        }

        for( unsigned IA = 0; IA < _NA; ++IA ) {
            dxs[n][IA] = dxse[IA] + dxsi[IA];
        }
    }
}

void ICRU::ELINIT( std::vector<double> IZ, std::vector<double> STF ) {
    unsigned IE    = 0;
    double FGRID   = 10.0e0;
    unsigned IGRID = 10;
    while( IE < _NE ) {
        double EV = EGRD[IGRID] * FGRID;
        if( IGRID == EGRD.size() - 1 ) {
            IGRID = 0;
            FGRID *= 10.0e0;
        }
        _ET[IE]  = EV;
        _ETL[IE] = std::log( _ET[IE] );
        IE++;
        IGRID++;
    }
    for( unsigned i = 0; i < _NA; ++i ) {
        _XMUL[i] = std::log( std::max( _XMU[i], 1.0e-35 ) );
    }
    for( unsigned IE = 0; IE < _NE; ++IE ) {
        _ECS[IE]   = 0.0e0;
        _ETCS1[IE] = 0.0e0;
        _ETCS2[IE] = 0.0e0;
        _PCS[IE]   = 0.0e0;
        _PTCS1[IE] = 0.0e0;
        _PTCS2[IE] = 0.0e0;
        for( unsigned IA = 0; IA < _NA; ++IA ) {
            _EDCS( IE, IA ) = 0.0e0;
            _PDCS( IE, IA ) = 0.0e0;
        }
    }

    for( unsigned IEL = 0; IEL < _NELEM; ++IEL ) {
        int IZZ     = IZ[IEL];
        double STFF = STF[IEL];
        unsigned NS = IZ[IEL];
        if( NS > 999 ) NS = 999;

        //  Electrons.
        std::string PATHEL = _settings->GetDataDir() + "elsep-database/";
        std::stringstream ss;
        ss << std::setw( 3 ) << std::setfill( '0' ) << NS;

        std::string FILE1 = PATHEL + "eeldx" + ss.str() + ".tab";

        _DCSI.resize( _NA );
        std::ifstream file( FILE1 );
        if( !file.good() ) {
            ErrorMessages::FileNotFoundError( FILE1, CURRENT_FUNCTION );
        }

        for( unsigned IE = 0; IE < _NE; ++IE ) {
            int IELEC, IZR;
            double ENR, CSE, TCS1E, TCS2E;
            file >> IELEC >> IZR >> ENR >> CSE >> TCS1E >> TCS2E;
            if( particle != ParticleType::ELECTRON || IZR != IZZ || std::fabs( ENR - _ET[IE] ) > 1.0e-3 ) {
                ErrorMessages::Error( "Corrupt data file for elsep-database.", CURRENT_FUNCTION );
            }
            for( unsigned IA = 0; IA < _NA; ++IA ) file >> _DCSI[IA];
            _ECS[IE]   = _ECS[IE] + STFF * CSE;
            _ETCS1[IE] = _ETCS1[IE] + STFF * TCS1E;
            _ETCS2[IE] = _ETCS2[IE] + STFF * TCS2E;
            for( unsigned IA = 0; IA < _NA; ++IA ) {
                _EDCS( IE, IA ) = _EDCS( IE, IA ) + STFF * _DCSI[IA];
            }
        }
        file.close();

        FILE1 = PATHEL + "peldx" + ss.str() + ".tab";
        file.open( FILE1 );
        for( unsigned IE = 0; IE < _NE; ++IE ) {
            int IELEC, IZR;
            double ENR, CSP, TCS1P, TCS2P;
            file >> IELEC >> IZR >> ENR >> CSP >> TCS1P >> TCS2P;
            if( IELEC != +1 || IZR != IZZ || std::fabs( ENR - _ET[IE] ) > 1.0e-3 ) {
                ErrorMessages::Error( "Corrupt data file", CURRENT_FUNCTION );
            }

            for( unsigned IA = 0; IA < _NA; ++IA ) file >> _DCSI[IA];
            _PCS[IE]   = _PCS[IE] + STFF * CSP;
            _PTCS1[IE] = _PTCS1[IE] + STFF * TCS1P;
            _PTCS2[IE] = _PTCS2[IE] + STFF * TCS2P;
            for( unsigned IA = 0; IA < _NA; ++IA ) {
                _PDCS( IE, IA ) += STFF * _DCSI[IA];
            }
        }

        file.close();
    }
}

void ICRU::DCSEL0( double E ) {
    if( E < 49.9999e0 || E > 1.0001e8 ) {
        ErrorMessages::Error( "Outside of supported energy range", CURRENT_FUNCTION );
    }
    double EL = std::log( E );
    unsigned JE;
    FINDI( _ETL, EL, _NE, JE );
    std::vector<double> Y( _NE );
    for( unsigned IA = 0; IA < _NA; ++IA ) {
        for( unsigned IE = 0; IE < _NE; ++IE ) {
            if( particle == ParticleType::ELECTRON )
                Y[IE] = std::log( _EDCS( IE, IA ) );
            else if( particle == ParticleType::POSITRON )
                Y[IE] = std::log( _PDCS( IE, IA ) );
            else {
                ErrorMessages::Error( "Unsupported particle type", CURRENT_FUNCTION );
            }
        }
        std::vector<double> A, B, C, D;
        SPLINE( _ETL, Y, A, B, C, D, 0.0e0, 0.0e0, _NE );
        _DCSIL[IA] = A[JE] + EL * ( B[JE] + EL * ( C[JE] + EL * D[JE] ) );
        _DCSI[IA]  = std::exp( _DCSIL[IA] );
    }
    std::vector<double> AA, BB, CC, DD;
    SPLINE( _XMUL, _DCSIL, AA, BB, CC, DD, 0.0e0, 0.0e0, _NA );
    for( unsigned IE = 0; IE < _NE; ++IE ) {
        if( particle == ParticleType::ELECTRON )
            Y[IE] = std::log( _ECS[IE] );
        else if( particle == ParticleType::POSITRON )
            Y[IE] = std::log( _PCS[IE] );
        else {
            ErrorMessages::Error( "Unsupported particle type", CURRENT_FUNCTION );
        }
    }
    std::vector<double> A, B, C, D;
    SPLINE( _ETL, Y, A, B, C, D, 0.0e0, 0.0e0, _NE );
    _CSI = std::exp( A[JE] + EL * ( B[JE] + EL * ( C[JE] + EL * D[JE] ) ) );

    for( unsigned IE = 0; IE < _NE; ++IE ) {
        if( particle == ParticleType::ELECTRON )
            Y[IE] = std::log( _ETCS1[IE] );
        else if( particle == ParticleType::POSITRON )
            Y[IE] = std::log( _PTCS1[IE] );
        else {
            std::cerr << "Unsupported particle type" << std::endl;
            exit( EXIT_FAILURE );
        }
    }
    SPLINE( _ETL, Y, A, B, C, D, 0.0e0, 0.0e0, _NE );
    _TCS1I = std::exp( A[JE] + EL * ( B[JE] + EL * ( C[JE] + EL * D[JE] ) ) );

    for( unsigned IE = 0; IE < _NE; ++IE ) {
        if( particle == ParticleType::ELECTRON )
            Y[IE] = std::log( _ETCS2[IE] );
        else if( particle == ParticleType::POSITRON )
            Y[IE] = std::log( _PTCS2[IE] );
        else {
            ErrorMessages::Error( "Unsupported particle type", CURRENT_FUNCTION );
        }
        SPLINE( _ETL, Y, A, B, C, D, 0.0e0, 0.0e0, _NE );
        _TCS2I = std::exp( A[JE] + EL * ( B[JE] + EL * ( C[JE] + EL * D[JE] ) ) );
    }
}

double ICRU::DCSEL( double RMU ) {
    double RMUL = std::log( std::max( RMU, 1.0e-35 ) );
    unsigned I;
    FINDI( _XMUL, RMUL, _NA, I );
    return std::exp( _DCSIL[I] + ( _DCSIL[I + 1] - _DCSIL[I] ) * ( ( RMUL - _XMUL[I] ) / ( _XMUL[I + 1] - _XMUL[I] ) ) );
}

void ICRU::SPLINE( const std::vector<double>& X,
                   const std::vector<double>& Y,
                   std::vector<double>& A,
                   std::vector<double>& B,
                   std::vector<double>& C,
                   std::vector<double>& D,
                   double S1,
                   double SN,
                   unsigned N ) {
    if( N < 4 ) {
        ErrorMessages::Error( "Too few points", CURRENT_FUNCTION );
    }
    unsigned N1 = N - 2;
    unsigned N2 = N - 3;

    A.resize( N );
    B.resize( N );
    C.resize( N );
    D.resize( N );

    for( unsigned I = 0; I < N1; ++I ) {
        if( X[I + 1] - X[I] < 1.0e-13 ) {
            std::cerr << "x values not increasing" << std::endl;
            exit( EXIT_FAILURE );
        }
        A[I] = X[I + 1] - X[I];
        D[I] = ( Y[I + 1] - Y[I] ) / A[I];
    }

    for( unsigned I = 0; I < N2; ++I ) {
        B[I]       = 2.0e0 * ( A[I] + A[I + 1] );
        unsigned K = N1 - I;
        D[K]       = 6.0e0 * ( D[K] - D[K - 1] );
    }
    D[1] -= A[0] * S1;
    D[N1] = D[N1] - A[N1] * SN;

    for( unsigned I = 1; I < N2 + 1; ++I ) {
        double R = A[I] / B[I - 1];
        B[I] -= R * A[I];
        D[I + 1] -= R * D[I];
    }

    D[N1] = D[N1] / B[N2];
    for( unsigned I = 1; I < N2; ++I ) {
        unsigned K = N1 - I;
        D[K]       = ( D[K] - A[K] * D[K + 1] ) / B[K - 1];
    }
    D[N - 1] = SN;

    double SI1 = S1;
    for( unsigned I = 0; I < N1; ++I ) {
        double SI = SI1;
        SI1       = D[I + 1];
        double H  = A[I];
        double HI = 1.0e0 / H;
        A[I]      = ( HI / 6.0e0 ) * ( SI * X[I + 1] * X[I + 1] * X[I + 1] - SI1 * X[I] * X[I] * X[I] ) + HI * ( Y[I] * X[I + 1] - Y[I + 1] * X[I] ) +
               ( H / 6.0e0 ) * ( SI1 * X[I] - SI * X[I + 1] );
        B[I] = ( HI / 2.0e0 ) * ( SI1 * X[I] * X[I] - SI * X[I + 1] * X[I + 1] ) + HI * ( Y[I + 1] - Y[I] ) + ( H / 6.0e0 ) * ( SI - SI1 );
        C[I] = ( HI / 2.0e0 ) * ( SI * X[I + 1] - SI1 * X[I] );
        D[I] = ( HI / 6.0e0 ) * ( SI1 - SI );
    }

    double FN  = Y[N - 1];
    double FNP = B[N1] + X[N - 1] * ( 2.0e0 * C[N1] + X[N - 1] * 3.0e0 * D[N1] );
    A[N - 1]   = FN - X[N - 1] * FNP;
    B[N - 1]   = FNP;
    C[N - 1]   = 0.0e0;
    D[N - 1]   = 0.0e0;
}

void ICRU::FINDI( std::vector<double> X, double XC, unsigned N, unsigned& I ) {
    if( XC > X[N - 1] ) {
        I = N - 1;
        return;
    }
    if( XC < X[0] ) {
        I = 0;
        return;
    }
    I           = 0;
    unsigned I1 = N;
    while( I1 - I > 1 ) {
        double IT = ( I + I1 ) / 2;
        if( XC > X[IT] )
            I = IT;
        else
            I1 = IT;
    }
}

double ICRU::DCSIN3( double E, double RMU ) {
    double res = 0.0e0;
    for( unsigned K = 0; K < _NOSC; ++K ) {
        if( E > _WRI[K] ) {
            double WR  = _WRI[K];
            double CSD = 0.0, CSC = 0.0;
            DOICSS( WR, E, RMU, CSD, CSC );
            res += _F[K] * ( CSC + CSD );
        }
    }
    res *= R4PI;
    return res;
}

void ICRU::DOICSS( double WR, double E, double RMU, double& DCSD, double& DCSC ) {
    double EL = E;
    double WU;
    if( particle == ParticleType::ELECTRON )
        WU = 0.5e0 * EL;
    else if( particle == ParticleType::POSITRON )
        WU = EL;
    else {
        ErrorMessages::Error( "Unsupported particle type", CURRENT_FUNCTION );
    }

    DCSD = 0.0e0;
    DCSC = 0.0e0;
    if( EL < WR + 1.0e-3 ) return;

    double CP2   = E * ( E + TEMC2 );
    double CP    = std::sqrt( CP2 );
    _BETA2       = CP2 / ( ( E + EMC2 ) * ( E + EMC2 ) );
    double GAMM1 = E / EMC2;
    double GAMMA = 1.0e0 + GAMM1;
    double CONS1 = CONS / _BETA2;

    double CPP2 = ( E - WR ) * ( E - WR + TEMC2 );
    double CPP  = std::sqrt( CPP2 );
    double XPP  = ( CP - CPP ) / EMC2;
    if( XPP > 1.0e-4 )
        _C1 = XPP * XPP;
    else {
        double WREG = WR / ( E * ( GAMMA + 1.0e0 ) );
        double XC   = 1.0e0 + 0.5e0 * WREG * ( ( 1.0e0 / GAMMA ) + WREG );
        _C1         = ( ( WR * XC ) * ( WR * XC ) ) / ( _BETA2 * EMC2 * EMC2 );
    }

    _R1 = ( GAMM1 / GAMMA ) * ( GAMM1 / GAMMA );

    if( particle == ParticleType::POSITRON ) {
        double G12 = ( GAMMA + 1.0e0 ) * ( GAMMA + 1.0e0 );
        _B1        = _R1 * ( 2.0e0 * G12 - 1.0e0 ) / ( GAMMA * GAMMA - 1.0e0 );
        _B2        = _R1 * ( 3.0e0 + 1.0e0 / G12 );
        _B3        = _R1 * 2.0e0 * GAMMA * ( GAMMA - 1.0e0 ) / G12;
        _B4        = _R1 * GAMM1 * GAMM1 / G12;
    }

    _C2   = 4.0e0 * CP * CPP / ( EMC2 * EMC2 );
    _C3   = 0.5e0 * ( WR / EMC2 ) * ( WR / EMC2 );
    _C0   = CONS1 * 0.5e0 * _C2 / WR;
    _RMU1 = ( ( WR / EMC2 ) * ( ( WR / EMC2 ) + 2.0e0 ) - _C1 ) / _C2;

    double CPPP2 = ( E - WU ) * ( E - WU + TEMC2 );
    double CPPP  = std::sqrt( CPPP2 );
    if( CPPP > 1.0e-9 )
        _RMU2 = ( WU * ( WU + TEMC2 ) - ( CP - CPPP ) * ( CP - CPPP ) ) / ( 4.0e0 * CP * CPPP );
    else
        _RMU2 = 0.5e0;
    _C4 = CONS1 * CP2 * TEMC2;

    DCSD = std::max( DODCSD( RMU ), 0.0e0 );
    DCSC = std::max( DODCSC( RMU, EL ), 0.0e0 );
}

double ICRU::DODCSD( double RMU ) {
    double res;

    if( RMU < 0.0e0 || RMU >= _RMU1 ) {
        res = 0.0e0;
        return res;
    }
    double QQM = _C1 + _C2 * RMU;
    double QR;
    if( QQM > 1.0e-4 ) {
        QR  = std::sqrt( QQM + 1.0e0 ) - 1.0e0;
        QQM = 0.5e0 * QQM;
    }
    else {
        QQM = 0.5e0 * QQM;
        QR  = QQM * ( 1.0e0 - 0.5e0 * QQM * ( 1.0e0 - QQM ) );
    }
    double FACT = 1.0e0 + _C3 * ( _BETA2 * QQM - _C3 ) / ( ( QQM - _C3 ) * ( QQM - _C3 ) );
    res         = ( _C0 / ( QQM * ( 1.0e0 + QR ) ) ) * FACT;
    return res;
}

double ICRU::DODCSC( double RMU, double EL ) {
    double res;
    if( RMU <= _RMU1 || RMU >= _RMU2 ) {
        res = 0.0e0;
        return res;
    }
    double EMU = EL * 2.0e0 * RMU * ( 1.0e0 - RMU );
    double WN  = ( EL + TEMC2 ) * EMU;
    double WE  = ( WN / ( EMU + EMC2 ) ) / EL;
    double R   = 0.0;
    if( particle == ParticleType::ELECTRON ) {
        R = 1.0e0 + ( WE / ( 1.0e0 - WE ) ) * ( WE / ( 1.0e0 - WE ) ) - ( 1.0e0 - _R1 ) * ( WE / ( 1.0e0 - WE ) ) + _R1 * WE * WE;
    }
    else if( particle == ParticleType::POSITRON ) {
        R = 1.0e0 - WE * ( _B1 - WE * ( _B2 - WE * ( _B3 - WE * _B4 ) ) );
    }
    else {
        ErrorMessages::Error( "Unsupported particle type", CURRENT_FUNCTION );
    }
    res = _C4 * R * ( 1.0e0 - 2.0e0 * RMU ) / ( WN * WN );
    return res;
}

void ICRU::GetAngularScatteringXS( Matrix& angularXS, Vector& integratedXS ) {
    angularXS.resize( _QMU.size(), _E.size() );

    Vector QMU_trafo( _QMU.size() );
    // for( unsigned n = 0; n < mu.size(); ++n ) mu[n] = 1.0 - 2.0 * _XMU[n];    // trafo

    // Ruecktrafo für _QM
    for( unsigned n = 0; n < _QMU.size(); ++n ) {
        QMU_trafo[n] = ( 1 - _QMU[n] ) / 2.0;    // trafo
    }
    // std::sort( mu.begin(), mu.end() );

    integratedXS.resize( _E.size() );
    for( unsigned i = 0; i < _E.size(); ++i ) {
        std::vector<double> dxse, dxsi, dxs;
        angdcs( 3u, _E[i], dxse, dxsi, dxs );    // speedup here possible
        Interpolation interpDXS( _XMU, dxs );
        blaze::column( angularXS, i ) = interpDXS( QMU_trafo );

        for( unsigned j = 0; j < _XMU.size() - 1; ++j ) {
            integratedXS[i] += 0.5 * ( _XMU[j + 1] - _XMU[j] ) * ( dxs[j + 1] + dxs[j] );
        }
    }
    angularXS *= H2OMolecularDensity;
    integratedXS *= H2OMolecularDensity;
}

void ICRU::GetTransportCoefficients( Matrix& xi ) {
    Vector mu( _XMU.size() );
    std::vector<std::vector<double>> dxs;
    // get scattering cross sections for all energies E
    angdcs( 3u, dxs );
    for( unsigned n = 0; n < mu.size(); ++n ) mu[n] = 1.0 - 2.0 * _XMU[n];    // mu starts at 1 and goes to 0
#pragma omp parallel for
    for( unsigned i = 0; i < _E.size(); ++i ) {
        // compute moments with trapezoidal rule
        for( unsigned n = 0; n < xi.rows(); ++n ) {
            xi( n, i ) = 0.0;
            // integration over mu
            for( unsigned k = 0; k < dxs[i].size() - 1; ++k ) {
                xi( n, i ) -= 0.5 * ( pow( 1.0 - mu[k + 1], n ) * dxs[i][k + 1] + pow( 1.0 - mu[k], n ) * dxs[i][k] ) * ( mu[k + 1] - mu[k] );
            }
        }
    }
    xi *= 2.0 * PI * H2OMolecularDensity;
}

void ICRU::GetStoppingPower( Vector& stoppingPower ) {
    std::string PATHMT = _settings->GetDataDir() + "material/";
    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << materialID;

    std::string FILE1 = PATHMT + "mater" + ss.str() + ".stp";

    std::vector<double> energy, stppwr;

    std::ifstream f( FILE1 );
    if( !f.good() ) {
        ErrorMessages::FileNotFoundError( FILE1, CURRENT_FUNCTION );
    }

    std::string line;
    for( unsigned i = 0; i < 23; ++i ) std::getline( f, line );
    double e, ecol, erad, pcol, prad;
    while( std::getline( f, line ) ) {
        std::stringstream stream( line );
        stream >> e >> ecol >> erad >> pcol >> prad;
        energy.push_back( e * 1e-6 );    // store as MeV
        if( particle == ELECTRON ) {
            stppwr.push_back( ecol + erad );
        }
        else if( particle == POSITRON ) {
            stppwr.push_back( pcol + prad );
        }
        else {
            ErrorMessages::Error( "Invalid particle type!", CURRENT_FUNCTION );
        }
    }
    Interpolation interp( energy, stppwr );
    stoppingPower.resize( _E.size() );
    for( unsigned i = 0; i < _E.size(); ++i ) {
        stoppingPower[i] = interp( _E[i] / 1e6 );
    }
}
