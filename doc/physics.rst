================
Transport Theory
================

The kinetic theory is dedicated to describe the dynamical behavior of a many-particle system through ensemble averaging.
Particles, e.g. molecules, photons, neutrons, electrons and plasmas, travel along the trajectories and undergo occasional collisions that change their directions and energies.
Such dynamics can be formulated via operator splitting approach in the Boltzmann equation, i.e.

.. math::
    :label: boltzmann

    \frac{\partial \psi}{\partial t}+v \cdot \nabla_x \psi = Q(\psi)

where :math:`\psi` is the one-particle probability distribution function, :math:`v` is velocity and :math:`Q(\psi)` is the scattering term.
In the Kit-RT solver, we consider the liner Boltzmann equation

.. math::
    :label: linearbz

    \partial_{t} f(v)+v \cdot \nabla_{x} f(v)=\sigma \int k\left(v, v^{\prime}\right)\left(f\left(v^{\prime}\right)-f(v)\right) d v^{\prime}-\tau f(v)

where the particles don't interact with one another but scatter with the background material.
For convenience, we reformulate the particle velocity into polar coordinates.
