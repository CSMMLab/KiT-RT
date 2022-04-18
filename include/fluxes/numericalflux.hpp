#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H

#include "common/typedef.hpp"

class NumericalFluxBase
{
  public:
    NumericalFluxBase();
    virtual ~NumericalFluxBase() {}

    static NumericalFluxBase* Create();

    /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return numerical flux value
     */
    virtual double Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const = 0;

     /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge (in XZ plane)
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return numerical flux value
     */
    virtual double FluxXZ( const Vector& Omega, double psiL, double psiR, const Vector& n ) const = 0;

    /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge for pseudo 1D case
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return numerical flux value
     */
    virtual double Flux1D( const Vector& Omega, double psiL, double psiR, const Vector& n ) const = 0;

    /**
     * @brief Flux       Computes "Steger Warming" upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
     *                   resultFlux
     * @param AxPlus     Positive part of the flux jacobian in x direction
     * @param AxMinus    Negative part of the flux jacobian in x direction
     * @param AyPlus     Positive part of the flux jacobian in y direction
     * @param AyMinus    Negative part of the flux jacobian in y direction
     * @param AzPlus     Positive part of the flux jacobian in z direction
     * @param AzMinus    Negative part of the flux jacobian in z direction
     * @param psiL       Solution state of left hand side control volume
     * @param psiR       Solution state of right hand side control volume
     * @param n          Normal vector at the edge between left and right control volume
     * @return           Vector with resulting flux.
     */
    virtual Vector Flux1D( const Matrix AxPlus, const Matrix AxMinus, const Vector psiL, const Vector psiR, const Vector n ) const = 0;

    /**
     * @brief Flux       Computes "Steger Warming" upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
     *                   resultFlux
     * @param AxPlus     Positive part of the flux jacobian in x direction
     * @param AxMinus    Negative part of the flux jacobian in x direction
     * @param AyPlus     Positive part of the flux jacobian in y direction
     * @param AyMinus    Negative part of the flux jacobian in y direction
     * @param AzPlus     Positive part of the flux jacobian in z direction
     * @param AzMinus    Negative part of the flux jacobian in z direction
     * @param psiL       Solution state of left hand side control volume
     * @param psiR       Solution state of right hand side control volume
     * @param n          Normal vector at the edge between left and right control volume
     * @return           Vector with resulting flux.
     */
    virtual Vector Flux( const Matrix AxPlus,
                         const Matrix AxMinus,
                         const Matrix AyPlus,
                         const Matrix AyMinus,
                         const Matrix AzPlus,
                         const Matrix AzMinus,
                         const Vector psiL,
                         const Vector psiR,
                         const Vector n ) const = 0;

    /**
     * @brief Flux       Computes "Steger Warming" upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
     *                   resultFlux in 2D, X and Z direction
     * @param AxPlus     Positive part of the flux jacobian in x direction
     * @param AxMinus    Negative part of the flux jacobian in x direction
     * @param AyPlus     Positive part of the flux jacobian in y direction
     * @param AyMinus    Negative part of the flux jacobian in y direction
     * @param AzPlus     Positive part of the flux jacobian in z direction
     * @param AzMinus    Negative part of the flux jacobian in z direction
     * @param psiL       Solution state of left hand side control volume
     * @param psiR       Solution state of right hand side control volume
     * @param n          Normal vector at the edge between left and right control volume
     * @return           Vector with resulting flux.
     */
    virtual Vector FluxXZ( const Matrix AxPlus,
                           const Matrix AxMinus,
                           const Matrix AyPlus,
                           const Matrix AyMinus,
                           const Matrix AzPlus,
                           const Matrix AzMinus,
                           const Vector psiL,
                           const Vector psiR,
                           const Vector n ) const = 0;

    virtual void FluxVanLeer( const Matrix& Ax,
                              const Matrix& AxAbs,
                              const Matrix& Ay,
                              const Matrix& AyAbs,
                              const Matrix& Az,
                              const Matrix& AzAbs,
                              const Vector& psiL,
                              const Vector& psiR,
                              const Vector& n,
                              Vector& resultFlux ) const = 0;
};

#endif    // NUMERICALFLUX_H
