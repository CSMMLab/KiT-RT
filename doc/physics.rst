================
Kinetic Theory
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

The physical world shows a diverse set of behaviors on different
characteristic scales. Consider the molecular motion of gases as an
example. Down to the finest scale of a many-particle system, the
Newton’s second law depicts particle motions via

.. math::

   \mathbf{F} = m \mathbf{a}

As a first order system it reads

.. math::

   \frac{d \mathbf x}{dt} = \mathbf v, \ \frac{d \mathbf v}{dt} = \frac{\mathbf F}{m}

An intuitive numerical algorithm is to get the numerous particles on
board and track the trajectories of them. A typical example is the
`Molecular Dynamics`_. This is not going to be efficient since there are
more than ``2e25`` molecules per cubic meter in normal atmosphere, and
things get even more complicated when you count on the N-body
interactions all the time. Some methods have been proposed to simplify
the computation. As an example, the `Direct simulation Monte Carlo`_
employs certain molecular models and conduct the intermolecular
collisions in a stochastic manner. It significantly reduces the
computational cost, while the trade-off is the artificial fluctuations.
Many realizations must be simulated successively to average the
solutions and reduce the errors.

An alternative strategy is made from ensemble averaging, where the
coarse-grained modeling is used to provide a bottom-up view. At the mean
free path and collision time scale of molecules, particles travel freely
during most of time with mild intermolecular collisions. Such dynamics
can be described with an operator splitting approach, i.e. the kinetic
transport equation

.. math::

   \frac{\partial f}{\partial t}+ \mathbf v \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf v f = Q(f)

where the left and right hand sides model particle transports and
collisions correspondingly. Different collision models can be inserted
into such equation. If the particles only collide with a background
material one obtains linear Boltzmann collision operator

.. math::

   Q(f)=\int_{\mathbb R^3} \mathcal B(\mathbf v_*, \mathbf v) \left[ f(\mathbf v_*)-f(\mathbf v)\right] d\mathbf v_*

where the collision kernel ``\mathcal B`` models the strength of
collisions at different velocities. If the interactions among particles
are considered, the collision operator becomes nonlinear. For example,
the two-body collision results in nonlinear Boltzmann equation

.. math::

   Q(f)=\int_{\mathbb R^3} \int_{\mathcal S^2} \mathcal B(\cos \beta, |\mathbf{v}-\mathbf{v_*}|) \left[ f(\mathbf v')f(\mathbf v_*')-f(\mathbf v)f(\mathbf v_*)\right] d\mathbf \Omega d\mathbf v_*

