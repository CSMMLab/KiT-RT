.. _implementation:

================
Implementation
================

Finite volume method
------------------------

In the KiT-RT, we employ the finite volume method (FVM) to model and compute the particle evolutions.
It's a generic method for conservation laws.
Consider the following PDE,

.. math::
	
	\frac{\partial \mathbf{u}}{\partial t}+\nabla \cdot \mathbf{f}(\mathbf{u})=\mathbf{0}

Here, :math:`\mathbf{u}` represents any vector of states and 
:math:`\mathbf{f}` represents the corresponding flux tensor.
To solve the equation numerically, we can sub-divide the spatial domain into finite cells.
For a particular cell :math:`i`, we take the volume integral over the total volume of the cell, which gives,

.. math::
	
	\int_{v_{i}} \frac{\partial \mathbf{u}}{\partial t} d v+\int_{v_{i}} \nabla \cdot \mathbf{f}(\mathbf{u}) d v=\mathbf{0}.

On integrating the first term to get the volume average and applying the divergence theorem to the second, this yields

.. math::

	v_{i} \frac{d \overline{\mathbf{u}}_{i}}{d t}+\oint_{S_{i}} \mathbf{f}(\mathbf{u}) \cdot \mathbf{n} d S=\mathbf{0},

where :math:`S_i` represents the total surface area of the cell and :math:`\mathbf n` is a unit vector normal to the surface and pointing outward. 
The equivalent formulation results

.. math::

	\frac{d \overline{\mathbf{u}}_{i}}{d t}+\frac{1}{v_{i}} \oint_{S_{i}} \mathbf{f}(\mathbf{u}) \cdot \mathbf{n} d S=\mathbf{0}.