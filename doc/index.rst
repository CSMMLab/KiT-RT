===========================================================
KiT-RT: A kinetic transport solver for radiation therapy
===========================================================

The KiT-RT framework is an open source project for radiation transport written in C++ programming language.
It is interested in the evolution of many-particle systems, e.g. photons, neutrons, electrons, etc. 
Based on the finite volume method (FVM), it provides an efficient tool to solve the Boltzmann and related equations in multiple dimensions.
Special attention has been paid to the application of radiation therapy and treatment planning.
The framework provides rich deterministic solver types for different domain-specific problems.
A list of current supported models and equations is as follows.

- linear Boltzmann (:math:`S_N`) equation
- spherical harmonics (:math:`P_N`) moment equations
- entropy-closure (:math:`M_N`) moment equations
- continuous slowing down equation

The source code is publicly available on `Github <https://github.com/CSMMLab/KiT-RT>`_.

Design philosophy
------------------------
The code hierarchy is designed as intuitive and neat as possible.
It's dedicated to providing a friendly interface for educational usage in kinetic theory and rich functionality for scientific research. 
Benefiting from the brilliant expressiveness and low-overhead abstraction provided by the C++ programming language, we implement the easily extensible code structure,
which allow the users to focus on physics or easily extend the codes with their own methods.


What is new?
------------------------
Finite volume method is a proven approach for simulating conservation laws. 
Compared with the existing open-source softwares, e.g. OpenFOAM, SU2 and Clawpack, Kit-RT holds the novelty through the following points:

- Comprehensive support for kinetic theory and phase-space equations
- Special focus on radiation therapy
- Lightweight design to ensure the flexibility for secondary development


How to get help?
------------------------
The software is being developed by members of the group `CSMM <https://www.scc.kit.edu/en/aboutus/rg-csmm.php>`_ at the Karlsruhe Institute of Technology (KIT).
If you are interested in using KiT-RT or are trying to figure out how to use it, please feel free to get in touch with `us <authors.html>`_.
Do open an issue or pull request if you have questions, suggestions or solutions.


Table of contents
------------------------

.. toctree::
   :maxdepth: 1
  
   installation
   physics
   implement
   configFiles
   cpp_doc
   developer_guide
   authors
