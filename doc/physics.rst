================
Theory
================

The Boltzmann equation
----------------------

The particle transport phenomena enjoy rich academic research value and application prospects.
A many-particle system can exhibit different behaviors at characteristic different scales.
Down to the finest scale of a many-particle system, the Newton’s second law depicts particle motions via

.. math::

   F = m a,

which leads

.. math::

   \frac{d x}{dt} = v, \ \frac{d v}{dt} = \frac{F}{m}.

An intuitive numerical solution algorithm is to get the numerous particles on board and track the trajectories of them. 
A typical example is the molecular dynamics (MD) method.
This is not going to be efficient since there are more than :math:`2\times 10^{25}` molecules per cubic meter in normal atmosphere, 
and things get extremely complicated if the N-body interactions are counted all the time. 

Simplifications can be conducted to accelerate the numerical computation.
As an example, the Monte Carlo method employs certain particle models and conduct the interactions in a stochastic manner. 
It significantly reduces the computational cost, while the trade-off is the artificial fluctuations.
Many realizations must be simulated successively to average the solutions and reduce the errors.

An alternative strategy can be made from ensemble averaging, where the
coarse-grained modeling is used to provide a bottom-up view. 
At the mean free path and collision time scale of particles. Such dynamics can be described with kinetic theory.
The Boltzmann equation can be formulated via an operator splitting approach, i.e.

.. math::

   \partial_{t} f(v)+v \cdot \nabla_{x} f(v)=\int_{\mathcal R^3} \int_{\mathcal S^2} k\left(v, v^{\prime}\right) \left(f\left(v^{\prime}\right)f\left(v_*^{\prime}\right)-f(v)f(v_*)\right) d\Omega d v_*,

where the left and right hand sides model particle transports and collisions correspondingly. 
The distribution function :math:`f` is the probability of finding a particle with certain location, and :math:`\{v, v_*\}` denotes the velocities of two classes of colliding particles. 
The collision kernel :math:`k` models the strength of collisions at different velocities.

Different collision models can be inserted into the Boltzmann equation.
In the KiT-RT solver, we are interested in the linear Boltzmann equation, where the particles don't interact with one another but scatter with the background material.
Therefore, the Boltzmann can be simplified as the linear equation with respect to :math:`f`

.. math::

    \partial_{t} f(v)+v \cdot \nabla_{x} f(v)=\int k\left(v, v^{\prime}\right)\left(f\left(v^{\prime}\right)-f(v)\right) d v^{\prime}-\tau f(v).

For convenience, it is often reformulated with polar coordinates :math:`\{r, \phi, \theta \}`,

.. math::

    &\left[\frac{1}{v(E)} \partial_{t} +\Omega \cdot \nabla+\Sigma_t (t, r, E)\right] \psi(t, r, \Omega, E) \\
    &=\int_{0}^{\infty} d E^{\prime} \int_{\mathcal R^2} d \Omega^{\prime} \Sigma_{s}\left(r, \Omega^{\prime} \bullet \Omega, E^{\prime} \rightarrow E\right) \psi\left(t, r, \Omega^{\prime}, E^{\prime}\right) + Q(t, r, \Omega, E).

The particle distribution :math:`\psi(r, \Omega, E, t)` here is often named as angular flux, :math:`\{\Sigma_s, \Sigma_t \}` are the scattering and total cross sections correspondingly, and :math:`Q` denotes a source term.


The spherical harmonics moment equations
----------------------------------------

The spherical harmonics (:math:`P_N`) method (cf. Brunner and Holloway [2005]) is one of several ways to solve the equation of radiative transfer. 
It serves as an approximate method, i.e. the method of moments, to reduce the high dimensionality when the original kinetic equation of radiative transfer, which is formulated on a seven-dimensional domain.
Let us consider the radiative transfer equation in one-dimensional physical space with only one energy group, i.e.

.. math::

   \partial_{t} \psi(t, z, \mu) &+\mu \nabla_{z} \psi(t, z, \mu)+\Sigma_{t}(t, z) \psi(t, z, \mu) \\
   &=\int_{\mathcal S^{2}} \Sigma_{s}\left(t, z, \mu \cdot \Omega^{\prime}\right) \psi\left(t, z, \mu^{\prime}\right) d \mu^{\prime}+q(t, z, \mu),

where :math:`\mu` is the projected angular variable on :math:`z` axis.

To obtain the :math:`P_N` equations, we express the angular dependence of the distribution function in terms of a Fourier series,

.. math::

   \psi(t, z, \mu)=\sum_{\ell=0}^{\infty} \psi_{\ell}(t, z) \frac{2 \ell+1}{2} P_{\ell}(\mu),

where :math:`P_{\ell}` are the Legendre polynomials.
These form an orthogonal basis of the space
of polynomials with respect to the standard scalar product on :math:`[−1, 1]`.
We can then obtain

.. math::

   \partial_{t} \psi_{\ell}+\partial_{z} \int_{=1}^{1} \mu P_{\ell} \psi \mathrm{d} \mu+\Sigma_{t \ell} \psi_{\ell}=q_{\ell},

where 

.. math::

   \Sigma_{t \ell}=\Sigma_{t}-\Sigma_{s \ell}=\Sigma_{a}+\Sigma_{s 0}-\Sigma_{s \ell},  \quad \Sigma_{s \ell}=2 \pi \int_{-1}^{1} P_{\ell}(\mu) \Sigma_{s}(\mu) \mathrm{d} \mu.

Two properties of the spherical harmonics are crucial for our method. These appear here as properties of the Legendre polynomials. First, we observe that, by this
procedure, we have diagonalized the scattering operator on the right-hand side (the
Legendre polynomials are the eigenfunctions of scattering operator). 
Second, a general property of orthogonal polynomials is that they satisfy a recursion relation. In
particular, the Legendre polynomials satisfy

.. math::

   \mu P_{\ell}(\mu)=\frac{\ell}{2 \ell+1} P_{\ell-1}(\mu)+\frac{\ell+1}{2 \ell+1} P_{\ell+1}(\mu).

Using this fact and truncating the expansion at :math:`\ell = N`, we arrive at the slab-geometry
:math:`P_N` equations,

.. math::

   \partial_{t} \psi_{\ell}+\partial_{z}\left(\frac{\ell+1}{2 \ell+1} \psi_{\ell+1}+\frac{\ell}{2 \ell+1} \psi_{\ell-1}\right)+\Sigma_{t \ell} \psi_{\ell}=q_{\ell}.

The above method can be extended to multi-dimensional case with the help of spherical harmonics, which are defined as

.. math::

   Y_{\ell}^{m}(\mu, \phi)=(-1)^{m} \sqrt{\frac{2 \ell+1}{4 \pi} \frac{(\ell-m) !}{(\ell+m) !}} e^{i m \phi} P_{\ell}^{m}(\mu),

where :math:`\ell \leq 0` and :math:`\ell \leq m \leq -\ell`.


The entropy closure moment equations
------------------------------------

Another method of moments comes from the minimal principle of a convex entropy to close the moment system.
Derivation of such moment system begins with the choice of a vector-valued function
:math:`m: \mathbb{S}^{2} \rightarrow \mathbb{R}^{n}, \Omega \mapsto\left[m_{0}(\Omega), \ldots, m_{n-1}(\Omega)\right]^{T}`,
whose n components are linearly independent functions of :math:`\Omega`.
Evolution equations for the moments u(x, t) :=
hmψ(x, ·, t)i are found by multiplying the transport equation by m and integrating
over all angles to give

.. math::

   \frac{1}{v} \partial_{t} u+\nabla_{x} \cdot\langle\Omega m \psi\rangle=\langle m \mathcal{C}(\psi)\rangle.

The system above is not closed; a recipe, or closure, must be prescribed to express
unknown quantities in terms of the given moments. Often this is done via an
approximation for :math:`\psi` that depends on :math:`u`,

.. math::

   \psi(x, \Omega, t) \simeq \mathcal{E}(u(x, t))(\Omega).

A general strategy for prescribing a closure is to
use the solution of a constrained optimization problem

.. math::
   :label: closure

   \min_{g \in \operatorname{Dom}(\mathcal{H})} & \mathcal{H}(g) \\
   \quad \text { s.t. } & \langle\mathbf{m} g\rangle=\langle\mathbf{m} \psi\rangle=u,

where :math:`\mathcal H(g)=\langle \eta(g) \rangle` and $\eta: \mathbb R \rightarrow \mathbb R$
is a convex function that is related to
the entropy of the system. For photons, the physically relevant entropy comes from
Bose-Einstein statistics

.. math::

   \eta(g)=\frac{2 k \nu^{2}}{v^{3}}\left[n_{g} \log \left(n_{g}\right)-\left(n_{g}+1\right) \log \left(n_{g}+1\right)\right],

where :math:`n_g` is the occupation number associated with g,

.. math::

   n_{g}:=\frac{v^{2}}{2 h \nu^{3}} g.

The solution of :eq:`closure` is expressed in terms of the Legendre dual

.. math::

   \eta_{*}(f)=-\frac{2 k \nu^{2}}{v^{3}} \log \left(1-\exp \left(-\frac{h \nu c}{k} f\right)\right).

Let

.. math::

   \mathcal{B}(\boldsymbol{\alpha}):=\eta_{*}^{\prime}\left(\boldsymbol{\alpha}^{T} \mathbf{m}\right)=\frac{2 h \nu^{3}}{v^{2}} \frac{1}{\exp \left(-\frac{h \nu c}{k} \boldsymbol{\alpha}^{T} \mathbf{m}\right)-1},

then the solution of :eq:`closure` is given by :math:`\mathcal B(\hat \alpha)`, where :math:`\hat \alpha= \hat \alpha(u)` solves the
dual problem

.. math::

   \min _{\boldsymbol{\alpha} \in \mathbb{R}^{n}}\left\{\left\langle\eta_{*}\left(\boldsymbol{\alpha}^{T} \mathbf{m}\right)\right\rangle-\boldsymbol{\alpha}^{T} \mathbf{u}\right\}.


The continuous slowing down approximation
-----------------------------------------

For the radiation therapy, the main goal is to compute the radiation dose accurately, which is defined as

.. math::

    D(x)=\frac{1}{\rho(x)}\int_0^{\infty}\int_{\mathbb{S}^2}S(E,x)\psi(E,x,\Omega)\,d\Omega dE.

Here, :math:`E\in\mathbb{R}_+` is the energy, :math:`\mathbf{x}\in \mathbf{X}\subset\mathbb{R}^3` denotes the spatial domain, 
and :math:`\mathbf{\Omega}\in\mathbb{S}^2` is the flight direction of particles. 
Moreover, :math:`\psi:\mathbb{R}_+\times\mathbb{R}^3\times\mathbb{S}^2\rightarrow\mathbb{R}` denotes the radiation flux density and 
:math:`\rho:\mathbb{R}^3\rightarrow\mathbb{R}` is the patient tissue density. 
The stopping power :math:`S:\mathbb{R}_+\times\mathbb{R}^3 \rightarrow \mathbb{R}` models the continuous energy loss of particles due to scattering with tissue and is defined as

.. math::
    S(E,x) = \int_0^{\infty} E'\int_{-1}^1\Sigma(E,E',x,\mu)d\mu dE'.

with the scattering cross section :math:`\Sigma:\mathbb{R}_+\times \mathbb{R}_+\times \mathbb{R}^3\times[-1,1]\rightarrow \mathbb{R}`.
The radiation flux density, which describes the probability of finding a particle at a certain region in phase space, can be computed from the continuous slowing down (CSD) approximation of the energy dependent linear Boltzmann equation

.. math::
    &-\partial_E\left(S(E,x)\psi(E,x,\Omega)\right)+\Omega\cdot\nabla_x\psi(E,x,\Omega)+\Sigma_t(E,x)\psi(E,x,\Omega) \\
    &= \int_{\mathbb{S}^2}\Sigma_s(E,x,\Omega\cdot\Omega')\psi(E,x,\Omega')d\Omega'.

Here :math:`E\in\mathbb{R}_+` is energy, :math:`x\in D\subset \mathbb{R}^3` is the spatial domain and :math:`\Omega\in\mathbb{S}^2` is the direction of travel. 
This model assumes a continuous energy loss of particles traveling through a background material, which is modeled using the stopping power :math:`S`. 
The scattering cross section :math:`\Sigma_s(E,\mathbf x,\mathbf \Omega\cdot\mathbf \Omega')` denotes the probability of particles at position :math:`\mathbf x` with energy :math:`E` changing their flight direction from :math:`\mathbf \Omega'` 
to :math:`\mathbf\Omega` due to a collision with the patient tissue. The total cross section :math:`\Sigma_t` is given by

Let us assume there are no absorption effects, and thus the total cross section is given by

.. math::

    \Sigma_t(E,x) = \Sigma_{s,0}(E,x)=2\pi \int_{-1}^1\Sigma_s(E,x,\mu)d\mu.

To simplify the evaluation of material properties, we follow the common assumption that all materials are
water-equivalent and differ only in density, i.e.,

.. math::
    &S(E,x) = S^{H_2O}(E)\rho(x), \\
    &\Sigma_t(E,x) = \Sigma_t^{H_2O}(E)\rho(x), \\
    &\Sigma_s(E,x,\Omega\cdot\Omega') = \Sigma_s(E,\Omega\cdot\Omega')\rho(x).

Leaving out the superscript :math:`H_2O`, the CSD equation can be simplified as

.. math::
   :label: CSD2

    &-\partial_E\left(\rho(x)S(E)\psi(E,x,\Omega)\right)+\Omega\cdot\nabla_x\psi(E,x,\Omega)+\rho(x)\Sigma_t(E)\psi(E,x,\Omega) \\
    &= \int_{\mathbb{S}^2}\rho(x)\Sigma_s(E,\Omega\cdot\Omega')\psi(E,x,\Omega')d\Omega'.    

Now, we bring this system in a form which resembles the standard Boltzmann equation. 
Multiplying :eq:`CSD2` with :math:`S(E)` gives

.. math::
   :label: CSD3

   \begin{align}
      -S(E)\partial_E\left(S(E)\rho(x)\psi(E,x,\Omega)\right)+&\Omega\cdot\nabla_x S(E)\psi(E,x,\Omega)+\Sigma_t(E)S(E)\rho(x)\psi(E,x,\Omega)\\ 
      &= \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')S(E)\rho(x)\psi(E,x,\Omega')d\Omega'.    
   \end{align}

Then, we substitute

.. math::
    \widehat{\psi}(E,x,\Omega):= S(E)\rho(x)\psi(E,x,\Omega)

into :eq:`CSD3`, which yields

.. math::
   :label: CSD4
    
    & -S(E)\partial_E\widehat{\psi}(E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widehat{\psi}(E,x,\Omega)}{\rho}+\Sigma_t(E)\widehat{\psi}(E,x,\Omega) \\
    & = \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')\widehat{\psi}(E,x,\Omega')d\Omega'.    

Now, to get rid of the stopping power in front of the energy derivative, we make use of the transformation

.. math::
   :label: TildeE

    \widetilde{E}(E) = \int_0^E \frac{1}{S(E')}\,dE'.

Now let us change to

.. math::
    \widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) := \widehat{\psi}(E(\widetilde E),x,\Omega)

In this case, the energy derivative becomes

.. math::
    \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{E}\widetilde{\widehat{\psi}}( E,x,\Omega)\partial_{\widetilde E }E(\widetilde E(\widetilde E) = \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega){S(E(\widetilde E))}.

And by rearranging the terms, we finally get

.. math::
    \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)\frac{1}{S(E(\widetilde E))},

since :math:`S(E(\widetilde E))` is nonzero.
Therefore, substituting :math:`\widetilde E` in :eq:`CSD4` gives

.. math::
   :label: CSD5

    & -\partial_{\widetilde E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)}{\rho}+\widetilde\Sigma_t(\widetilde E)\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) \\
    & = \int_{\mathbb{S}^2}\widetilde\Sigma_s(\widetilde E,\Omega\cdot\Omega')\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega')d\Omega'.

Here, we define :math:`\widetilde\Sigma_{t}(\widetilde E):=\Sigma_t(E(\widetilde E))` and :math:`\widetilde\Sigma_{s}(\widetilde E,\Omega\cdot\Omega'):=\Sigma_s(E(\widetilde E),\Omega\cdot\Omega')`. Finally, to obtain a positive sign in front of the energy derivative, we transform to

.. math::
    \bar{E}(\widetilde{E}) = \widetilde{E}_{\text{max}}-\widetilde{E}.

Then, with :math:`\bar{\psi}(\bar{E},x,\Omega):=\widetilde{\widehat{\psi}}(\widetilde{E}(\bar{E}),x,\Omega)`, :math:`\bar\Sigma_{t}(\bar E):=\widetilde{\Sigma}_t(\widetilde{E}(\bar{E}))` as well as :math:`\bar\Sigma_{s}(\bar E,\Omega\cdot\Omega'):=\widetilde{\Sigma}_s(\widetilde{E}(\bar{E}),\Omega\cdot\Omega')` equation :eq:`CSD4` becomes

.. math::
   :label: CSD6

    \partial_{\bar{E}}\bar{\psi}(\bar{E},x,\Omega)+\Omega\cdot\nabla_x \frac{\bar{\psi}(\bar{E},x,\Omega)}{\rho}+\bar\Sigma_t(\bar E)\bar{\psi}(\bar{E},x,\Omega) = \int_{\mathbb{S}^2}\bar\Sigma_s(\bar{E},\Omega\cdot\Omega')\bar{\psi}(\bar{E},x,\Omega')d\Omega'.

Dropping the bar notation and treating :math:`\bar E` as a pseudo-time :math:`t` gives a slightly modified version of the Boltzmann equation

.. math::
   :label: CSDBoltzmann

    \partial_{t}\psi(t,x,\Omega)+&\Omega\cdot\nabla_x \frac{\psi(t,x,\Omega)}{\rho}+\Sigma_t(t)\psi(t,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(t,\Omega\cdot\Omega')\psi(t,x,\Omega')d\Omega'\\
    &\psi(t=0,x,\Omega) = S(E_{\text{max}})\rho(x)\psi(E_{\text{max}},x,\Omega).

We are interested in computing the dose, which (when again using the original energy :math:`E` and angular flux :math:`\psi`) reads

.. math::
    D(x) = \int_0^{\infty} \int_{\mathbb{S}^2} S(E)\psi(E,x,\Omega)\,d\Omega dE = \int_0^{\infty} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\widehat\psi(E,x,\Omega)\,d\Omega dE.

So let us check how we can compute the dose from our solution :math:`\bar \psi(\bar E,x,\Omega)`. For this, let us substitute

.. math::

    \bar E(E) = \tilde{E}(E_{max}) - \int_0^E \frac{1}{S(E')}dE'.

We have

.. math::

    \frac{d\bar E(E)}{dE} = -\frac{1}{S(E)}

which gives

.. math::
    D(x) =& -\int_{\infty}^{0} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\bar \psi(\bar E,x,\Omega)\frac{1}{S(E(\bar E))}\,d\Omega d\bar E\\
    =& \int_{0}^{\infty} \frac{1}{\rho(x)S(E(\bar E))}\int_{\mathbb{S}^2} \bar \psi(\bar E,x,\Omega)\,d\Omega d\bar E.
