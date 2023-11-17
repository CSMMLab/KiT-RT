#include "problems/lattice.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"

// ---- Checkerboard Sn ----
// Constructor for Ckeckerboard case with Sn
Lattice_SN::Lattice_SN(Config *settings, Mesh *mesh) : ProblemBase(settings, mesh)
{

    // Initialise crosssections to 1
    _scatteringXS = Vector(_mesh->GetNumCells(), 1.0);
    _totalXS = Vector(_mesh->GetNumCells(), 1.0);

    // For absorption cells: set scattering XS to 0 and absorption to 10
    auto cellMids = _mesh->GetCellMidPoints();
    for (unsigned j = 0; j < cellMids.size(); ++j)
    {
        if (isAbsorption(cellMids[j]))
        {
            _scatteringXS[j] = 0.0;
            _totalXS[j] = 10.0;
        }
    }
    SetGhostCells();
}

Lattice_SN::~Lattice_SN() {}

VectorVector Lattice_SN::GetScatteringXS(const Vector &energies) { return VectorVector(energies.size(), _scatteringXS); }

VectorVector Lattice_SN::GetTotalXS(const Vector &energies) { return VectorVector(energies.size(), _totalXS); }

std::vector<VectorVector> Lattice_SN::GetExternalSource(const Vector & /*energies*/)
{
    VectorVector Q(_mesh->GetNumCells(), Vector(1u, 0.0));
    auto cellMids = _mesh->GetCellMidPoints();
    for (unsigned j = 0; j < cellMids.size(); ++j)
    {
        if (isSource(cellMids[j]))
            Q[j] = _settings->GetSourceMagnitude(); // isotropic source
    }
    return std::vector<VectorVector>(1u, Q);
}

VectorVector Lattice_SN::SetupIC()
{
    VectorVector psi(_mesh->GetNumCells(), Vector(_settings->GetNQuadPoints(), 1e-10));
    return psi;
}

bool Lattice_SN::isAbsorption(const Vector &pos) const
{
    // Check whether pos is inside absorbing squares
    double xy_corrector = -3.5;
    std::vector<double> lbounds{1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector};
    std::vector<double> ubounds{2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector};
    for (unsigned k = 0; k < lbounds.size(); ++k)
    {
        for (unsigned l = 0; l < lbounds.size(); ++l)
        {
            if ((l + k) % 2 == 1 || (k == 2 && l == 2) || (k == 2 && l == 4))
                continue;
            if (pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l])
            {
                return true;
            }
        }
    }
    return false;
}

bool Lattice_SN::isSource(const Vector &pos) const
{
    // Check whether pos is part of source region
    if (pos[0] >= 3 - 3.5 && pos[0] <= 4 - 3.5 && pos[1] >= 3 - 3.5 && pos[1] <= 4 - 3.5)
        return true;
    else
        return false;
}

void Lattice_SN::SetGhostCells()
{
    // Loop over all cells. If its a Dirichlet boundary, add cell to dict with {cell_idx, boundary_value}
    auto cellBoundaries = _mesh->GetBoundaryTypes();
    std::map<int, Vector> ghostCellMap;

    QuadratureBase *quad = QuadratureBase::Create(_settings);
    VectorVector vq = quad->GetPoints();
    unsigned nq = quad->GetNq();

    Vector void_ghostcell(nq, 0.0);

    for (unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++)
    {
        if (cellBoundaries[idx_cell] == BOUNDARY_TYPE::NEUMANN || cellBoundaries[idx_cell] == BOUNDARY_TYPE::DIRICHLET)
        {
            ghostCellMap.insert({idx_cell, void_ghostcell});
        }
    }
    _ghostCells = ghostCellMap;

    delete quad;
}

const Vector &Lattice_SN::GetGhostCellValue(int idx_cell, const Vector &cell_sol)
{ // re-write to use pointers
    return _ghostCells[idx_cell];
}

// ---- Checkerboard Moments ----

// Constructor for checkerboard case with Pn
Lattice_Moment::Lattice_Moment(Config *settings, Mesh *mesh) : ProblemBase(settings, mesh)
{

    // Initialise crosssections = 1 (scattering)
    _scatteringXS = Vector(_mesh->GetNumCells(), 1.0);
    _totalXS = Vector(_mesh->GetNumCells(), 1.0);

    // for absorption regions change crosssections to all absorption
    auto cellMids = _mesh->GetCellMidPoints();
    for (unsigned j = 0; j < cellMids.size(); ++j)
    {
        if (isAbsorption(cellMids[j]))
        {
            _scatteringXS[j] = 0.0;
            _totalXS[j] = 10.0;
        }
    }
}

Lattice_Moment::~Lattice_Moment() {}

VectorVector Lattice_Moment::GetScatteringXS(const Vector &energies) { return VectorVector(energies.size(), _scatteringXS); }

VectorVector Lattice_Moment::GetTotalXS(const Vector &energies) { return VectorVector(energies.size(), _totalXS); }

std::vector<VectorVector> Lattice_Moment::GetExternalSource(const Vector & /*energies*/)
{
    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS

    double integrationFactor = (4 * M_PI);
    if (_settings->GetDim() == 2)
    {
        integrationFactor = M_PI;
    }
    SphericalBase *tempBase = SphericalBase::Create(_settings);
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector Q(_mesh->GetNumCells(), Vector(ntotalEquations, 0.0)); // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector uIC(ntotalEquations, 0);

    if (_settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS)
    {
        QuadratureBase *quad = QuadratureBase::Create(_settings);
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w = quad->GetWeights();

        double my, phi;
        VectorVector moments = VectorVector(quad->GetNq(), Vector(tempBase->GetBasisSize(), 0.0));

        for (unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++)
        {
            my = quadPointsSphere[idx_quad][0];
            phi = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis(my, phi);
        }
        // Integrate <1*m> to get factors for monomial basis in isotropic scattering
        for (unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++)
        {
            uIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }
    double kinetic_density = _settings->GetSourceMagnitude();
    for (unsigned j = 0; j < cellMids.size(); ++j)
    {
        if (isSource(cellMids[j]))
        {
            if (_settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS)
            {
                Q[j] = kinetic_density * uIC / uIC[0] / integrationFactor; // Remember scaling
            }
            if (_settings->GetSphericalBasisName() == SPHERICAL_HARMONICS)
            {
                Q[j][0] = kinetic_density / integrationFactor; // first bassis function is 1/ ( 4 * M_PI )
            }
        }
    }
    delete tempBase; // Only temporally needed
    return std::vector<VectorVector>(1u, Q);
}

VectorVector Lattice_Moment::SetupIC()
{
    double integrationFactor = (4 * M_PI);
    if (_settings->GetDim() == 2)
    {
        integrationFactor = M_PI;
    }
    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase *tempBase = SphericalBase::Create(_settings);
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector initialSolution(_mesh->GetNumCells(), Vector(ntotalEquations, 0)); // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector tempIC(ntotalEquations, 0);

    if (_settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS)
    {
        QuadratureBase *quad = QuadratureBase::Create(_settings);
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w = quad->GetWeights();

        double my, phi;
        VectorVector moments = VectorVector(quad->GetNq(), Vector(tempBase->GetBasisSize(), 0.0));

        for (unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++)
        {
            my = quadPointsSphere[idx_quad][0];
            phi = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis(my, phi);
        }
        // Integrate <1*m> to get factors for monomial basis in isotropic scattering
        for (unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++)
        {
            tempIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }
    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
    double kinetic_density = 1e-4;
    for (unsigned j = 0; j < cellMids.size(); ++j)
    {
        if (_settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS)
        {
            initialSolution[j] = kinetic_density * tempIC / tempIC[0] / integrationFactor; // Remember scaling
        }
        if (_settings->GetSphericalBasisName() == SPHERICAL_HARMONICS)
        {
            initialSolution[j][0] = kinetic_density / integrationFactor; // first bassis function is 1/ ( 4 * M_PI )
        }
    }
    delete tempBase; // Only temporally needed
    return initialSolution;
}

bool Lattice_Moment::isAbsorption(const Vector &pos) const
{
    // Check whether pos is in absorption region
    double xy_corrector = -3.5;
    std::vector<double> lbounds{1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector};
    std::vector<double> ubounds{2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector};
    for (unsigned k = 0; k < lbounds.size(); ++k)
    {
        for (unsigned l = 0; l < lbounds.size(); ++l)
        {
            if ((l + k) % 2 == 1 || (k == 2 && l == 2) || (k == 2 && l == 4))
                continue;
            if (pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l])
            {
                return true;
            }
        }
    }
    return false;
}

bool Lattice_Moment::isSource(const Vector &pos) const
{
    // Check whether pos is in source region
    if (pos[0] >= 3 - 3.5 && pos[0] <= 4 - 3.5 && pos[1] >= 3 - 3.5 && pos[1] <= 4 - 3.5)
        return true;
    else
        return false;
}
