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
    
    
Macroscopic Models
-----------------------------------------
This section discusses macroscopic models to :eq:`CSDBoltzmann`. These models are derived from nodal and modal discretizations of the directional domain. 

Modal discretizations
**************************
Modal discretizations of \eqref{eq:BoltzmannCSDTrafo} can be interpreted as a closure problem [Levermore1996Moment]_ , [Levermore1996Entropy]_. To present the derivation of different closures, we first formulate the moment equations which are not closed. Second, we close these equations with the $P_N$ closure and third, we derive the $M_N$ closure.

Moment equations
+++++++++++++++++++++++++
Let us derive an evolution equation to describe the moments of the radiation flux with respect to the real-valued spherical harmonics basis functions. These are defined as the real parts of the spherical harmonics

.. math::
Y_{\ell}^k(\mathbf{\Omega}) = \sqrt{\frac{2\ell +1}{4\pi}\frac{(\ell-k)!}{(\ell+k)!}}\ e^{ik\varphi}P_{\ell}^k(\mu) ,


where $P_{\ell}^k$ are the associated Legendre polynomials. Then, the real spherical harmonics are given as

.. math::
    m_{\ell}^k(\mathbf{\Omega}) = 
    \begin{cases}
        \frac{(-1)^k}{\sqrt{2}}\left( Y_{\ell}^k(\mathbf{\Omega}) + (-1)^k Y_{\ell}^{-k}(\mathbf{\Omega}) \right), & k > 0\;, \\
        Y_{\ell}^0(\mathbf{\Omega}) & k = 0 \;, \\
        -\frac{(-1)^k i}{\sqrt{2}}\left( Y_{\ell}^{-k}(\mathbf{\Omega}) - (-1)^k Y_{\ell}^{k}(\mathbf{\Omega}) \right), & k < 0\;,
    \end{cases}

where $i$ is the imaginary unit. Collecting all basis functions up to degree $N$ in a vector

.. math::
    \mathbf m = \left(m_0^0, m_1^{-1}, m_1^{0}, m_1^{1},\cdots, m_N^{N}\right)^T\in\mathbb{R}^{(N+1)^2}

yields the so-called moments

.. math::
    u_{\ell}^k(t,\mathbf x) := \int_{\mathbb{S}^2} \psi(t,\mathbf x,\mathbf\Omega) m_{\ell}^k(\mathbf \Omega) \intD \mathbf\Omega\;.

Evolution equations for the moments can be derived by testing :eq:`CSDBoltzmann` against :math:`\mathbf{m}_{\ell} = (m_{\ell}^{-\ell},\cdots,m_{\ell}^{\ell})`, which gives

.. math::
    \partial_{t}\mathbf u_{\ell}(t,\mathbf x)+&\nabla_x\cdot\int_{ \mathbb{S}^2}\mathbf\Omega\mathbf m_{\ell}(\mathbf\Omega)\frac{\psi(t,\mathbf x,\mathbf\Omega)}{\rho(\mathbf x)}\intD \mathbf{\Omega}+\Sigma_t(t)\mathbf u_{\ell}(t,\mathbf x)\nonumber\\
    &= \int_{\mathbb{S}^2}\int_{\mathbb{S}^2}\mathbf m_{\ell}(\mathbf\Omega)\Sigma_s(t,\mathbf\Omega\cdot\mathbf\Omega')\psi(t,\mathbf x,\mathbf\Omega')\intD \mathbf\Omega'\intD \mathbf\Omega\;.

To rewrite this equation, we use the spherical harmonics recursion relation
.. math::
    \Omega_i \mathbf{m}_{\ell} = \mathbf{a}_{\ell}^i\mathbf m_{\ell-1} + \mathbf{a}_{\ell+1}^i\mathbf m_{\ell+1} \enskip \text{ with } \mathbf{a}_{\ell}^i\in\mathbb{R}^{(2\ell-1)\times (2\ell+1)}

as well as the fact that there exists a diagonal matrix :math:`\bm{\Sigma}_{\ell}(t)$ with entries $\Sigma_{\ell,kk} = \Sigma_{\ell}^k := 2\pi\int_{[-1,1]}P_{\ell}(\mu)\Sigma_s(t,\mu)\intD \mu` such that
.. math::
    \int_{\mathbb{S}^2}\int_{\mathbb{S}^2}\mathbf m_{\ell}(\mathbf\Omega)\Sigma_s(t,\mathbf\Omega\cdot\mathbf\Omega')\psi(t,\mathbf x,\mathbf\Omega')\intD \mathbf\Omega'd\mathbf\Omega = \bm{\Sigma}_{\ell}(t) \mathbf u_{\ell}(t,\mathbf x)\;.

Then, the moment equations at degree $\ell$ become
.. math::
    \partial_{t}\mathbf u_{\ell}(t,\mathbf x)+&\sum_{i=1}^3\partial_{x_i}\left(\mathbf{a}_{\ell}^i\mathbf u_{\ell-1}(t,\mathbf x) + \mathbf{a}_{\ell+1}^i\mathbf u_{\ell+1}(t,\mathbf x)\right)/\rho(\mathbf{x})+\Sigma_t(t)\mathbf u_{\ell}(t,\mathbf x)\nonumber\\
    &= \bm{\Sigma}_{\ell}(t) \mathbf u_{\ell}(t,\mathbf x)\;.

Note that the equations for degree $\ell$ depend on the moments of degree $\ell+1$. Hence, to obtain a closed system of moments up to a fixed degree $N$, we need to define a closure relation 
.. math::
    \mathbf u_{N+1}(t,\mathbf x)\simeq \mathcal{U}(\mathbf u_{0}(t,\mathbf x),\cdots,\mathbf u_{N}(t,\mathbf x)) .
    

$P_N$ closure
++++++++++++++++++++++++++++

The most commonly used closure is the $P_N$ closure which leads to the spherical harmonics ($P_N$) method [Case1967linear]_. It expands the solution by spherical harmonics, i.e.,
.. math::
    \psi(t,\mathbf{x},\mathbf{\Omega}) \approx \psi_{\mathrm{P}_N}(t,\mathbf{x},\mathbf{\Omega}) := \mathbf{u}(t,\mathbf x)^T\mathbf{m}(\mathbf{\Omega}),

where :math:`\mathbf{u}\in\mathbb{R}^{(N+1)^2}` collects all moments according to :math:`\mathbf u = \left(u_0^0, u_1^{-1}, u_1^{0}, u_1^{1},\cdots, u_N^{N}\right)^T\in\mathbb{R}^{(N+1)^2}`.
Hence, the $P_N$ closure is simply given as $\mathcal{U}_{\mathrm{P}_N}\equiv \mathbf 0$. In this case, the moment equations read
.. math::
    \partial_t \mathbf u (t,\mathbf x) =-\mathbf A\cdot\nabla_{\mathbf{x}} \frac{\mathbf u(t,\mathbf x)}{\rho(\mathbf x)}-\Sigma_t(t)\mathbf u (t,\mathbf x)+\bm{\Sigma}\mathbf u (t,\mathbf x),

where :math:`\mathbf A\cdot\nabla_{\mathbf{x}} := \mathbf A_1\partial_{x} + \mathbf A_2\partial_y+ \mathbf A_3\partial_z` with :math:`\mathbf A_i := \int_{\mathbb{S}^2}\mathbf m\mathbf m^T \Omega_i \intD \mathbf{\Omega}` and :math:`\bm \Sigma = \mathrm{diag}\left(\Sigma_0^0, \Sigma_1^{-1}, \Sigma_1^{0}, \Sigma_1^{1},\cdots, \Sigma_N^{N}\right)`. While $P_N$ is a computationally efficient method (especially for scattering terms), it does not preserve positivity of the radiation flux approximation and can lead to spurious oscillations. A closure which mitigates oscillations and preserves positivity at significantly increased computational costs is the $M_N$ closure.

$M_N$ closure
+++++++++++++++++++++++++++++

The $M_N$ closure [Levermore1996Moment]_ , [Levermore1996Entropy]_ employs the principle of minimal mathematical, i.e., maximal physical entropy to close the moment system.
To this end, we define the twice differentiable, strictly convex kinetic entropy density :math:`\eta:\mathbb{R}_+\rightarrow\mathbb{R}`. Different, problem specific entropy densities can be defined, e.g. the Maxwell-Boltzmann entropy :math:`\eta(g)=g\log(g)-g`, or a quadratic entropy :math:`\eta(g)=g^2`, which recovers the $P_N$ method.
Thus, one can close the system by choosing the reconstructed radiation flux density :math:`\psi_{\mathbf u}` out of the set of all possible functions 
.. math::
    F_{\mathbf m}=\menge{g(t,x,\mathbf{\Omega})>0 : u=\int_{\mathbb{S}^2}{\mathbf m g}\intD \mathbf{\Omega}<\infty},

that fulfill the moment condition :math:`\mathbf u(t,\mathbf{x})=\inner{\mathbf m g}` as the one with minimal entropy. The modal basis $\mathbf m$ can be chosen arbitrarily. Common choices consist of spherical harmonics or other polynomial basis functions. The minimal entropy closure can be formulated as a constrained optimization problem for a given vector of moments $\mathbf u$,
.. math::
:label: EntropyOCP 
\min_{g\in F_{\mathbf m}} \int_{\mathbb{S}^2}\eta(g)\intD \mathbf{\Omega} \quad  \text{ s.t. }  \mathbf u=\int_{\mathbb{S}^2}{\mathbf m g}\intD \mathbf{\Omega}.

The minimal value of the objective function is denoted by :math:`h(u)=\inner{\eta(\psi_{\mathbf u})}` and describes the systems minimal entropy depending on time and space. :math:`\psi_{\mathbf u}` is the minimizer of :eq: `EntropyOCP`, which we use to close the moment system
.. math::
    \partial_{t}\mathbf u_{\ell}(t,\mathbf x)+&\nabla_x\cdot\int_{ \mathbb{S}^2}\mathbf\Omega\mathbf m_{\ell}(\mathbf\Omega)\frac{\psi_u(t,\mathbf x,\mathbf\Omega)}{\rho(\mathbf x)}\intD \mathbf{\Omega}+\Sigma_t(t)\mathbf u_{\ell}(t,\mathbf x)= \Sigma_{\ell}\mathbf u_{\ell} (t,\mathbf x);.

The minimal entropy method preserves important properties of the underlying equation  [Alldredge2018regularized]_ , [Levermore1996Moment]_ , i.e., positivity of the solution, hyperbolicity of the moment system, dissipation of mathematical entropy and the H-Theorem. The minimal entropy closed moment system is invariant under Galilean transforms. Lastly, if collision invariants of the Boltzmann equations are used as modal basis functions, then the moment system yields a local conservation law. 

The set of all moments corresponding to a radiation flux density :math:`\psi_{\mathbf u}>0` is called the realizable set 
.. math::
    \mathcal{R}=\menge{\mathbf u: \int_{\mathbb{S}^2}{\mathbf m g}\intD \mathbf{\Omega}=\mathbf u,\, g\in F_{\mathbf m}}.

Outside $\mathcal{R}$ the minimal entropy closure problem has no solution.  At the boundary $\partial \mathcal{R}$, the optimization problem becomes singular and :math:`\psi_{\mathbf u}` consists of a linear combination of dirac distributions. Near $\partial \mathcal{R}$ the entropy closure becomes ill conditioned and thus, a numerical optimizer requires a large amount of iterations to compute a solution.

To mitigate this issue, a regularized version of the entropy closure problem has been proposed by  [Alldredge2018regularized]_ ,
.. math::
:label: EntropyOCP_reg 
\inf_{g\in F_{\mathbf m}}  \int_{\mathbb{S}^2}\eta(g)\intD \mathbf{\Omega}+
\frac{1}{2\gamma}\norm{ \int_{\mathbb{S}^2}{\mathbf m g}\intD \mathbf{\Omega} - \mathbf u}^2_2,

where $\gamma$ is the regularization parameter. Generally, moments of the regularized reconstructed radiation flux density :math:`\int_{\mathbb{S}^2}\mathbf m\psi_{\mathbf u}\intD \mathbf{\Omega}` deviate from the non-regularized moments. 
For :math:`\gamma\rightarrow 0`, we recover the original entropy closure of :eq: `EntropyOCP` and the moments coincide again. The regularized entropy closure is solvable for any :math:`\mathbf u\in\mathbb{R}^{(N+1)^2}` and preserves all structural properties of the non-regularized entropy closure [Alldredge2018regularized]_. One can also choose to regularize only parts of the entropy closure, e.g. to preserve moments of specific interest. Then the partially regularized entropy closure reads
.. math::
:label: EntropyOCP_part_reg 
\inf_{g\in F_m}  \int_{\mathbb{S}^2}\eta(g)\intD \mathbf{\Omega} + \frac{1}{2\gamma}\norm{\int_{\mathbb{S}^2}{\mathbf m^r g} \intD \mathbf{\Omega} - u^r}^2_2\quad  \text{ s.t. }  \mathbf u^{nr}=\int_{\mathbb{S}^2}{\mathbf m^{nr}g}\intD \mathbf{\Omega},

where :math:`\mathbf u^{nr}` denotes non-regularized moment elements and :math:`\mathbf u^{r}` denotes regularized elements of the moment vector.

Recently, structure preserving deep neural networks have been successfully employed to approximate the entropy closure [Schotthoefer2021structurepreserving]_ to accelerate the $M_N$ method. The authors leverage convexity of the optimization problem and use corresponding input convex neural networks [AmosICNN]_ to preserve structural properties of the closure in the neural network based approximation.

Nodal discretizations
**********************************
The $S_N$ method [Lewis1984computational]_ employs a nodal discretization for the directional domain. To facilitate the computation of integral terms that arise due to scattering, the nodal point sets are commonly chosen according to a quadrature rule.

In the application case of radiative transport, the directional domain is assumed to be the unit sphere :math:`\mathbb{S}^2\subset\mathbb{R}^3`, thus a suitable parametrization is given by spherical coordinates
.. math::
    \mathbb{S}^2 =  \menge{ \begin{pmatrix}
           \sqrt{1-\mu^2}\sin(\varphi) \\
           \sqrt{1-\mu^2}\cos(\varphi) \\
           \mu
         \end{pmatrix}
     : \mu\in\left[-1,1\right],\, \varphi\in\left[0,2\pi\right)}.

Note, that we can allow different particle velocities by scaling the unit sphere with a given maximum velocity.
The implementation assumes a slab geometry setting, i.e., lower dimensional velocity spaces are described by a projection of $\mathbb{S}^2$ onto $\mathbb{R}^2$ and $\mathbb{R}$, respectively. Thus, the parametrization of the two-dimensional slab geometry is given by
.. math::
    P_{\mathbb{R}^2}(\mathbb{S}^2) =  \menge{ \begin{pmatrix}
           \sqrt{1-\mu^2}\sin(\varphi) \\
           \sqrt{1-\mu^2}\cos(\varphi)
         \end{pmatrix}
     : \mu\in\left[-1,1\right],\, \varphi\in\left[0,2\pi\right)}

and the one dimensional case is described by
.. math::
    P_{\mathbb{R}}(\mathbb{S}^2) =  \menge{ \mu     : \mu\in\left[-1,1\right]}

Hence the task is to derive a quadrature formula for the direction of travel. The perhaps most common approach is the product quadrature rule. Here. a Gauss quadrature is used for 
$\mu$ and equally weighted and spaced points for $\varphi$, i.e., when using $N_q$ points, we have
.. math::
\varphi_i = i \Delta\varphi \quad \text{for} \quad i=1,\ldots,N_q \quad \text{and} 
\quad \Delta\varphi = \frac{2\pi}{N_q}\;.

If the Gauss quadrature for :math:`\mu` uses :math:`N_q` points, then we obtain a total of :math:`Q = N_q^2` possible directions. Denoting the Gauss weights as :math:`w^G_k` with :math:`k = 1,\cdots,N_q`, we obtain the product quadrature weights 
.. math::
    w_{k\cdot N_q +\ell} = \frac{2\pi w^G_k}{N_q}

and points
.. math:: 
:label: SphericalCoordinatesProduct
 \mathbf \Omega_{k\cdot N_q +\ell}  = \begin{pmatrix}
           \sqrt{1-\mu_k^2}\sin(\varphi_{\ell}) \\
           \sqrt{1-\mu_k^2}\cos(\varphi_{\ell})
         \end{pmatrix} \;.

Other implemented quadrature methods include spherical Monte Carlo, Levelsymmetric [Longoni2004PhDT]_, LEBEDEV [Lebedev1986numerical]_ and LDFESA [Jarrel2011discrete]_ . A comparison of different quadrature sets and their approximation behaviour for $S_N$ methods can be found in [Camminady2021highly]_.

The evolution equations for :math:`\psi_q(t,\mathbf x):= \psi(t,\mathbf x,\mathbf \Omega_q)` are then given by
.. math::
:label: SNEqns
    \partial_{t}\psi_q(t,\mathbf x)+&\mathbf \Omega_q\cdot\nabla_x \frac{\psi_q(t,\mathbf x)}{\rho(\mathbf{x})}+\Sigma_t(t)\psi_q(t,\mathbf x) = \sum_{p=1}^{Q}w_p\Sigma_s(t,\mathbf \Omega_q\cdot\mathbf \Omega_p)\psi_p(t,\mathbf x)\;.

A main disadvantage of $S_N$ methods are so called ray-effects [Lathrop1968ray]_ , [Morel2003analysis]_ , [Mathews1999propagation]_ , which are spurious artifacts that stem from the limited number of directions in which particles can travel. Moreover, radiation therapy applications exhibit forward-peaked scattering, 
which cannot be well-captured by classical quadrature rules. 

To allow for moderate computational costs when computing scattering terms and to efficiently treat forward-peaked scattering, we transform the nodal solution to a modal description and apply the more efficient $P_N$ methodology for scattering terms. For this, we define a truncation order $N$ and construct the matrices :math:`\mathbf{O}\in\mathbb{R}^{Q \times (N+1)^2}` which maps the modal onto its nodal representation and :math:`\mathbf{M}\in\mathbb{R}^{(N+1)^2\times Q}` which maps the nodal onto its modal representation. Such matrices can be constructed by
.. math::
    \mathbf O = \left(\mathbf{m}(\Omega_k)\right)_{k=1}^{Q}\, , \text{ and } \enskip \mathbf M = \left(w_k\mathbf{m}(\Omega_k)\right)_{k=1}^{Q}.

In this case, we can replace the scattering term on the right-hand side of :eq: SNEqns by its $P_N$ counterpart. Collecting the nodal solution in the vector $\bm\psi$ then gives
.. math:: 
:label: SNEqns2
    \partial_{t}\bm\psi(t,\mathbf x)+&\mathbf \Omega_q\cdot\nabla_x \frac{\bm\psi(t,\mathbf x)}{\rho(\mathbf{x})}+\Sigma_t(t)\bm{\psi}(t,\mathbf x) = \mathbf{O}\bm{\Sigma}\mathbf{M}\bm{\psi}\;.

References
--------------------------------------

.. [Alldredge2018regularized] Graham W. Alldredge, Martin Frank, and Cory D. Hauck. 2018. A regularized entropy-based moment method for kinetic equations.
arXiv:1804.05447 [math.NA]

.. [AmosICNN] Brandon Amos, Lei Xu, and J. Zico Kolter. 2016. Input Convex Neural Networks. CoRR abs/1609.07152 (2016). arXiv:1609.07152
http://arxiv.org/abs/1609.07152

.. [Camminady2021highly] Thomas Camminady, Martin Frank, and Jonas Kusch. 2021. Highly uniform quadrature sets for the discrete ordinates method. Nuclear
Science and Engineering (2021).

.. [Case1967linear] Kenneth M Case and Paul Frederick Zweifel. 1967. Linear transport theory. (1967)

.. [Jarrel2011discrete] J.J. Jarrel and M.L. Adams. 2011. Discrete-ordinates quadrature sets based on linear discontinuous finite elements. Proc. International
Conference on Mathematics and Computational Methods applied to Nuclear Science and Engineering (2011)

.. [Lathrop1968ray] Kaye D Lathrop. 1968. Ray effects in discrete ordinates equations. Nuclear Science and Engineering 32, 3 (1968), 357–369.

.. [Lebedev1986numerical] G Marchuk and V I Lebedev. 1986. Numerical methods in the theory of neutron transport. (1986). https://www.osti.gov/biblio/7057084

.. [Levermore1996Moment] Levermore. 1996. Moment closure hierarchies for kinetic theories. Journal of Statistical Physics 83 (1996), 1021–1065.

.. [Levermore1996Entropy] C. David Levermore. 1997. Entropy-based moment closures for kinetic equations. Transport Theory and Statistical Physics 26, 4-5 (1997),
591–606

.. [Lewis1984computational] Elmer Eugene Lewis and Warren F Miller. 1984. Computational methods of neutron transport. (1984).

.. [Longoni2004PhDT] Gianluca Longoni. 2004. Advanced quadrature sets and acceleration and preconditioning techniques for the discrete ordinates method in
parallel computing environments. Ph. D. Dissertation. University of Florida.

.. [Mathews1999propagation] Kirk A Mathews. 1999. On the propagation of rays in discrete ordinates. Nuclear science and engineering 132, 2 (1999), 155–180.

.. [Morel2003analysis] JE Morel, TA Wareing, RB Lowrie, and DK Parsons. 2003. Analysis of ray-effect mitigation techniques. Nuclear science and engineering
144, 1 (2003), 1–22.

.. [Schotthoefer2021structurepreserving] Steffen Schotthöfer, Tianbai Xiao, Martin Frank, and Cory D. Hauck. 2022. Neural network-based, structure-preserving entropy closures for the Boltzmann moment system. https://doi.org/10.48550/ARXIV.2201.10364
