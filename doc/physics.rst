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



The continuous slowing down approximation
---------

Our main goal is to compute the radiation dose
\begin{align}
    D(x)=\frac{1}{\rho(x)}\int_0^{\infty}\int_{\mathbb{S}^2}S(E,x)\psi(E,x,\Omega)\,d\Omega dE.
\end{align}
The angular flux $\psi$ can be approximated by the continuous slowing down (CSD) equation, which reads
\begin{align}
    -\partial_E\left(S(E,x)\psi(E,x,\Omega)\right)+\Omega\cdot\nabla_x\psi(E,x,\Omega)+\Sigma_t(E,x)\psi(E,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(E,x,\Omega\cdot\Omega')\psi(E,x,\Omega')d\Omega'.
\end{align}
Here $E\in\mathbb{R}_+$ is energy, $x\in D\subset \mathbb{R}^3$ is the spatial domain and $\Omega\in\mathbb{S}^2$ is the direction of travel. The stopping power $S$ is given by
\begin{align}
    S(E,x) = \int_0^{\infty} E'\int_{-1}^1\Sigma(E,E',x,\mu)d\mu dE'.
\end{align}
Since there are no absorption effects, the total cross section is given by
\begin{align}
    \Sigma_t(E,x) = \Sigma_{s,0}(E,x)=2\pi \int_{-1}^1\Sigma_s(E,x,\mu)d\mu.
\end{align}
\jknote{We can double check this in our code.}
\jknote{Check if we need to transform when intgrating w.r.t. $\mu$ instead of $\Omega'$.}
\txnote{Why the integral deduces from $\Omega$ to $\mu$ automatically?}
\ssnote{Fixed that issue}
With a given background material density $\rho(x)$ now make the following assumptions
\begin{align}
    S(E,x) = S^{H_2O}(E)\rho(x), \\
    \Sigma_t(E,x) = \Sigma_t^{H_2O}(E)\rho(x), \\
    \Sigma_s(E,x,\Omega\cdot\Omega') = \Sigma_s(E,\Omega\cdot\Omega')\rho(x).
\end{align}
Leaving out the superscript $H_2O$, the CSD equation simplifies to
\begin{align}\label{eq:CSD2}
    -\partial_E\left(\rho(x)S(E)\psi(E,x,\Omega)\right)+\Omega\cdot\nabla_x\psi(E,x,\Omega)+\rho(x)\Sigma_t(E)\psi(E,x,\Omega) = \int_{\mathbb{S}^2}\rho(x)\Sigma_s(E,\Omega\cdot\Omega')\psi(E,x,\Omega')d\Omega'.    
\end{align}
Now, we bring this system in a form which resembles the standard Boltzmann equation. Multiplying \eqref{eq:CSD2} with $S(E)$ gives
\begin{subequations}\label{eq:CSD3}
\begin{align}
    -S(E)\partial_E\left(S(E)\rho(x)\psi(E,x,\Omega)\right)+&\Omega\cdot\nabla_x S(E)\psi(E,x,\Omega)+\Sigma_t(E)S(E)\rho(x)\psi(E,x,\Omega)\\ 
    &= \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')S(E)\rho(x)\psi(E,x,\Omega')d\Omega'.    
\end{align}
\end{subequations}
Then, we substitute  
\begin{align}
    \widehat{\psi}(E,x,\Omega):= S(E)\rho(x)\psi(E,x,\Omega)
\end{align}
into \eqref{eq:CSD3}, which yields
\begin{align}\label{eq:CSD4}
    -S(E)\partial_E\widehat{\psi}(E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widehat{\psi}(E,x,\Omega)}{\rho}+\Sigma_t(E)\widehat{\psi}(E,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(E,\Omega\cdot\Omega')\widehat{\psi}(E,x,\Omega')d\Omega'.    
\end{align}
Now, to get rid of the stopping power in front of the energy derivative, we make use of the transformation
\begin{align}\label{eq:TildeE}
    \widetilde{E}(E) = \int_0^E \frac{1}{S(E')}\,dE'.
\end{align}
Now let us change to
\begin{align}
    \widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) := \widehat{\psi}(E(\widetilde E),x,\Omega)
\end{align}
In this case, the energy derivative becomes
\begin{align*}
    \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{E}\widetilde{\widehat{\psi}}( E,x,\Omega)\partial_{\widetlde E }E(\widetilde E(\widetilde E) = \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega){S(E(\widetilde E))}.
\end{align*}
And by rearranging the terms, we finally get
\begin{align*}
    \partial_{ E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \partial_{\widetilde{E}}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)\frac{1}{S(E(\widetilde E))},
\end{align*}
since $S(E(\widetilde E))$ is nonzero \ssnote{Is S always nonzero? Would make sense, physically.}.
Therefore, substituting $\widetilde E$ in \eqref{eq:CSD4} gives
\begin{align}\label{eq:CSD5}
    -\partial_{\widetilde E}\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)+\Omega\cdot\nabla_x \frac{\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)}{\rho}+\widetilde\Sigma_t(\widetilde E)\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega) = \int_{\mathbb{S}^2}\widetilde\Sigma_s(\widetilde E,\Omega\cdot\Omega')\widetilde{\widehat{\psi}}(\widetilde E,x,\Omega')d\Omega'.
\end{align}
Here, we define $\widetilde\Sigma_{t}(\widetilde E):=\Sigma_t(E(\widetilde E))$ and $\widetilde\Sigma_{s}(\widetilde E,\Omega\cdot\Omega'):=\Sigma_s(E(\widetilde E),\Omega\cdot\Omega')$. Finally, to obtain a positive sign in front of the energy derivative, we transform to
\begin{align}
    \bar{E}(\widetilde{E}) = \widetilde{E}_{\text{max}}-\widetilde{E}.
\end{align}
Then, with $\bar{\psi}(\bar{E},x,\Omega):=\widetilde{\widehat{\psi}}(\widetilde{E}(\bar{E}),x,\Omega)$ and $\bar\Sigma_{t}(\bar E):=\widetilde{\Sigma}_t(\widetilde{E}(\bar{E}))$ as well as $\bar\Sigma_{s}(\bar E,\Omega\cdot\Omega'):=\widetilde{\Sigma}_s(\widetilde{E}(\bar{E}),\Omega\cdot\Omega')$ equation \eqref{eq:CSD4} becomes
\begin{align}\label{eq:CSD5}
    \partial_{\bar{E}}\bar{\psi}(\bar{E},x,\Omega)+\Omega\cdot\nabla_x \frac{\bar{\psi}(\bar{E},x,\Omega)}{\rho}+\bar\Sigma_t(\bar E)\bar{\psi}(\bar{E},x,\Omega) = \int_{\mathbb{S}^2}\bar\Sigma_s(\bar{E},\Omega\cdot\Omega')\bar{\psi}(\bar{E},x,\Omega')d\Omega'.
\end{align}
Dropping the bar notation and treating $\bar E$ as a pseudo-time $t$ gives a slightly modified version of the Boltzmann equation
\begin{align}\label{eq:CSDBoltzmann}
    \partial_{t}\psi(t,x,\Omega)+&\Omega\cdot\nabla_x \frac{\psi(t,x,\Omega)}{\rho}+\Sigma_t(t)\psi(t,x,\Omega) = \int_{\mathbb{S}^2}\Sigma_s(t,\Omega\cdot\Omega')\psi(t,x,\Omega')d\Omega'\\
    &\psi(t=0,x,\Omega) = S(E_{\text{max}})\rho(x)\psi(E_{\text{max}},x,\Omega).
\end{align}
We are interested in computing the dose, which (when again using the original energy $E$ and angular flux $\psi$) reads
\begin{align}
    D(x) = \int_0^{\infty} \int_{\mathbb{S}^2} S(E)\psi(E,x,\Omega)\,d\Omega dE = \int_0^{\infty} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\widehat\psi(E,x,\Omega)\,d\Omega dE.
\end{align}
So let us check how we can compute the dose from our solution $\bar \psi(\bar E,x,\Omega)$. For this, let us substitute
\begin{align}\label{eq:BarE}
    \bar E(E) = \tilde{E}(E_{max}) - \int_0^E \frac{1}{S(E')}dE'.
\end{align}
We have
\begin{align}
    \frac{d\bar E(E)}{dE} = -\frac{1}{S(E)}
\end{align}
which gives
\begin{align}
    D(x) =& -\int_{\infty}^{0} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\bar \psi(\bar E,x,\Omega)\frac{1}{S(E(\bar E))}\,d\Omega d\bar E\\
    =& \int_{0}^{\infty} \frac{1}{\rho(x)S(E(\bar E))}\int_{\mathbb{S}^2} \bar \psi(\bar E,x,\Omega)\,d\Omega d\bar E.
\end{align}
\jknote{@all: Can you double check this formulation of the dose? Kerstin had a different version in her thesis.}
\txnote{Kerstin (1.27)}
\begin{align}
    &\widehat{\psi}(E,x,\Omega) := \widetilde{\widehat{\psi}}(\widetilde E,x,\Omega)  :=\bar{\psi}(\bar{E},x,\Omega),\\
    &dE = -S(E) d\bar E(E), \\
    &D(x) = -\int_{\infty}^{0} \int_{\mathbb{S}^2} \frac{1}{\rho(x)}\bar \psi(\bar E,x,\Omega)S(E(\bar E))d\Omega d\bar E
    = \int_{0}^{\infty} \frac{S(E(\bar E))}{\rho(x)}\int_{\mathbb{S}^2} \bar \psi(\bar E,x,\Omega)\,d\Omega d\bar E.
\end{align}


\section{Numerical Implementation}
\subsection{Quadrature Grid}
In the field of radiation therapy, the radiative transfer equation is often discretized by Monte-Carlo methods, which commonly yield an increased runtime and stochastic noise. Therefore, we choose a discrete ordinates (S$_N$) discretization (which has its own disadvantages such as ray-effects) instead. The first step to obtain a discretized system is to derive a finite dimensional, nodal representation of the angular variable $\bm\Omega\in\mathbb{S}^2$. Here, different strategies to obtain an accurate discretization exist and we are going to start with a straight forward product quadrature rule on the sphere. To obtain a finite set of directions in $\mathbb{S}^2$, we describe the angular variable
\begin{align*}
\bm\Omega = (\cos\phi \sin \theta, \sin\phi \sin \theta, \cos \theta)^T\;
\end{align*}
at finite positions of $\theta$ and $\phi$.

\subsection{Energy Grid}
Furthermore, we choose a modal discretization of the energy (or pseudo-time) variable. Dividing the original energy $E$ into an equidistant mesh $E_{\text{min}} = E_1 < \cdots < E_{N_E} = E_{\text{max}}$ lets us define the transformed energies $\widetilde E$ according to \eqref{eq:TildeE} by
\begin{align}
    \widetilde E(E_n)\approx\widetilde E_n := \sum_{\ell=1}^n \Delta E \frac{1}{S(E_{\ell})} \qquad\text{for }n = 1,\cdots,N_E.
\end{align}
\jknote{We could also use a more accurate quadrature here in the code!}
Then, we change to $\bar{E}$ according to \eqref{eq:BarE}, i.e. on a discrete level we get
\begin{align}
    \bar E_n := \widetilde E_{N_E} - \widetilde E_n \qquad\text{for }n = 1,\cdots,N_E.
\end{align}
The cross section $\Sigma_t$ evaluated at $\bar E_n$, i.e. $\bar \Sigma_t(\bar E_n) = \Sigma_t(E(\bar E_n))$ is then simply given by evaluating $\Sigma_t$ at $E_{N_E-n}$. The same holds for the scattering cross section $\Sigma_s$.