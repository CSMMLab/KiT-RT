================
Theory
================

The Boltzmann equation
---------

The particle transport phenomena enjoy rich academic research value and application prospects.
A many-particle system can exhibit different behaviors at characteristic different scales.
Down to the finest scale of a many-particle system, the Newton’s second law depicts particle motions via

.. math::

   F = m a

which leads

.. math::

   \frac{d x}{dt} = v, \ \frac{d v}{dt} = \frac{F}{m}

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
The Boltzmann equation can be formulated via an operator splitting approach.

.. math::
   \partial_{t} f(v)+v \cdot \nabla_{x} f(v)=\int_{\mathcal R^3} \int_{\mathcal S^2} k\left(v, v^{\prime}\right) \left(f\left(v^{\prime}\right)f\left(v_*^{\prime}\right)-f(v)f(v_*)\right) d\Omega d v_*

where the left and right hand sides model particle transports and collisions correspondingly. 
The distribution function :math:`f` is the probability of finding a particle with certain location, and :math:`\{v, v_*\}` denotes the velocities of two classes of colliding particles. 
The collision kernel :math:`k` models the strength of collisions at different velocities.

Different collision models can be inserted into the Boltzmann equation.
In the KiT-RT solver, we are interested in the linear Boltzmann equation, where the particles don't interact with one another but scatter with the background material.
Therefore, the Boltzmann can be simplified as the linear equation with respect to :math:`f`.

.. math::
    :label: linbz
    
    \partial_{t} f(v)+v \cdot \nabla_{x} f(v)=\int k\left(v, v^{\prime}\right)\left(f\left(v^{\prime}\right)-f(v)\right) d v^{\prime}-\tau f(v)

For convenience, it is often reformulated with polar coordinates :math:`\{r, \phi, \theta \}`.

.. math::
    :label: porbz
   
    &\left[\frac{1}{v(E)} \frac{\partial}{\partial t}+\Omega \cdot \nabla+\Sigma_t (r, E, t)\right] \psi(r, \Omega, E, t) \\
    &=\int_{0}^{\infty} d E^{\prime} \int_{\mathcal R^2} d \Omega^{\prime} \Sigma_{s}\left(r, \Omega^{\prime} \bullet \Omega, E^{\prime} \rightarrow E\right) \psi\left(r, \Omega^{\prime}, E^{\prime}, t\right) + Q(r, \Omega, E, t)

The particle distribution :math:`\psi(r, \Omega, E, t)` here is often named as angular flux, :math:`\{\Sigma_s, \Sigma_t \}` are the scattering and total cross sections correspondingly, and :math:`Q` denotes a source term.


The continuous slowing down approximation
---------

For the radiation therapy, the main goal is to compute the radiation dose accurately, which is defined as

.. math::

    D(x)=\frac{1}{\rho(x)}\int_0^{\infty}\int_{\mathbb{S}^2}S(E,x)\psi(E,x,\Omega)\,d\Omega dE.

where :math:`\rho(x)` is the background material density.
The angular flux :math:`\psi` can be approximated by a further approximation equation, i.e. the continuous slowing down (CSD) equation which reads

.. math::
    &-\partial_E\left(S(E,x)\psi(E,x,\Omega)\right)+\Omega\cdot\nabla_x\psi(E,x,\Omega)+\Sigma_t(E,x)\psi(E,x,\Omega) \\
    &= \int_{\mathbb{S}^2}\Sigma_s(E,x,\Omega\cdot\Omega')\psi(E,x,\Omega')d\Omega'.

Here :math:`E\in\mathbb{R}_+` is energy, :math:`x\in D\subset \mathbb{R}^3` is the spatial domain and :math:`\Omega\in\mathbb{S}^2` is the direction of travel. 
The stopping power :math:`S` is given by

.. math::
    S(E,x) = \int_0^{\infty} E'\int_{-1}^1\Sigma(E,E',x,\mu)d\mu dE'.

Let us assume there are no absorption effects, and thus the total cross section is given by

.. math::

    \Sigma_t(E,x) = \Sigma_{s,0}(E,x)=2\pi \int_{-1}^1\Sigma_s(E,x,\mu)d\mu.

With a given :math:`\rho(x)`, we now make the following assumptions

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
Multiplying :ref:`CSD2` with :math:`S(E)` gives

.. math::
   :label: CSD3

   \begin{align}
      -S(E)\partial_E\left(S(E)\rho(x)\psi(E,x,\Omega)\right)+&\Omega\cdot\nabla_x S(E)\psi(E,x,\Omega)+\Sigma_t(E)S(E)\rho(x)\psi(E,x,\Omega)\\ 
      &= \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')S(E)\rho(x)\psi(E,x,\Omega')d\Omega'.    
   \end{align}

Then, we substitute

.. math::
    \widehat{\psi}(E,x,\Omega):= S(E)\rho(x)\psi(E,x,\Omega)

into \eqref{eq:CSD3}, which yields

.. math::
   :label: CSD4
    -S(E)\partial_E\widehat{\psi}(E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widehat{\psi}(E,x,\Omega)}{\rho}+\Sigma_t(E)\widehat{\psi}(E,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')\widehat{\psi}(E,x,\Omega')d\Omega'.    

Now, to get rid of the stopping power in front of the energy derivative, we make use of the transformation

.. math::
   :label: TildeE

    \widetilde{E}(E) = \int_0^E \frac{1}{S(E')}\,dE'.

Now let us change to

.. math::
    \widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) := \widehat{\psi}(E(\widetilde E),x,\Omega)

In this case, the energy derivative becomes

.. math::
    \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{E}\widetilde{\widehat{\psi}}( E,x,\Omega)\partial_{\widetlde E }E(\widetilde E(\widetilde E) = \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega){S(E(\widetilde E))}.

And by rearranging the terms, we finally get

.. math::
    \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)\frac{1}{S(E(\widetilde E))},

since :math:`S(E(\widetilde E))` is nonzero.
Therefore, substituting :math:`\widetilde E` in :ref:`CSD4` gives

.. math::
   :label: CSD5

    -\partial_{\widetilde E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)}{\rho}+\widetilde\Sigma_t(\widetilde E)\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \int_{\mathbb{S}^2}\widetilde\Sigma_s(\widetilde E,\Omega\cdot\Omega')\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega')d\Omega'.

Here, we define :math:`\widetilde\Sigma_{t}(\widetilde E):=\Sigma_t(E(\widetilde E))` and :math:`\widetilde\Sigma_{s}(\widetilde E,\Omega\cdot\Omega'):=\Sigma_s(E(\widetilde E),\Omega\cdot\Omega')`. Finally, to obtain a positive sign in front of the energy derivative, we transform to

.. math::
    \bar{E}(\widetilde{E}) = \widetilde{E}_{\text{max}}-\widetilde{E}.

Then, with $\bar{\psi}(\bar{E},x,\Omega):=\widetilde{\widehat{\psi}}(\widetilde{E}(\bar{E}),x,\Omega)$ and $\bar\Sigma_{t}(\bar E):=\widetilde{\Sigma}_t(\widetilde{E}(\bar{E}))$ as well as $\bar\Sigma_{s}(\bar E,\Omega\cdot\Omega'):=\widetilde{\Sigma}_s(\widetilde{E}(\bar{E}),\Omega\cdot\Omega')$ equation \eqref{eq:CSD4} becomes

.. math::
   :label: CSD6

    \partial_{\bar{E}}\bar{\psi}(\bar{E},x,\Omega)+\Omega\cdot\nabla_x \frac{\bar{\psi}(\bar{E},x,\Omega)}{\rho}+\bar\Sigma_t(\bar E)\bar{\psi}(\bar{E},x,\Omega) = \int_{\mathbb{S}^2}\bar\Sigma_s(\bar{E},\Omega\cdot\Omega')\bar{\psi}(\bar{E},x,\Omega')d\Omega'.

Dropping the bar notation and treating :math:`\bar E` as a pseudo-time :math:`t` gives a slightly modified version of the Boltzmann equation

.. math::
   :label: CSDBoltzmann

    \partial_{t}\psi(t,x,\Omega)+&\Omega\cdot\nabla_x \frac{\psi(t,x,\Omega)}{\rho}+\Sigma_t(t)\psi(t,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(t,\Omega\cdot\Omega')\psi(t,x,\Omega')d\Omega'\\
    &\psi(t=0,x,\Omega) = S(E_{\text{max}})\rho(x)\psi(E_{\text{max}},x,\Omega).

We are interested in computing the dose, which (when again using the original energy $E$ and angular flux $\psi$) reads

.. math::
    D(x) = \int_0^{\infty} \int_{\mathbb{S}^2} S(E)\psi(E,x,\Omega)\,d\Omega dE = \int_0^{\infty} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\widehat\psi(E,x,\Omega)\,d\Omega dE.

So let us check how we can compute the dose from our solution :math:`\bar \psi(\bar E,x,\Omega)`. For this, let us substitute

.. math::
   :label: BarE

    \bar E(E) = \tilde{E}(E_{max}) - \int_0^E \frac{1}{S(E')}dE'.

We have

.. math::

    \frac{d\bar E(E)}{dE} = -\frac{1}{S(E)}

which gives

.. math::
    D(x) =& -\int_{\infty}^{0} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\bar \psi(\bar E,x,\Omega)\frac{1}{S(E(\bar E))}\,d\Omega d\bar E\\
    =& \int_{0}^{\infty} \frac{1}{\rho(x)S(E(\bar E))}\int_{\mathbb{S}^2} \bar \psi(\bar E,x,\Omega)\,d\Omega d\bar E.
