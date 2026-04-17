---
bibliography:
  - readme.bib
---
# The `mtln` solver module

This module allows to solve networks of multiconductor transmission line bundles in the time domain. 
In the context of transmission line theory, a working definition of each of the terms above is:

* A **multiconductor tranmission line** is tranmission line composed of more than one conductor (and the reference conductor), i.e a tranmssion line composed of 3 or more conductors, where one of them is taken as the reference.
* A **bundle** is a series of shielded transmission lines grouped together. A bundle can contain other bundles. 
* A **network** is a series of interconnected transmission lines, multiconductor tranmission lines and/or bundles.

It also includes the following features, representing physical rather than topological properties of the transmission lines:

* The possibility to include dispersive elements, i.e. elements whose properties depend on the frequency of the excitation
* The effect of transfer impedances between a shield and the shielded conductores therein. Real-life shields are not perfect, and the currents circulating on the shields induce a voltage on the conductors inside them. These transfer impedances are often frequency-dependent, and their treatment is similar to that of dispersive elements. 
* The possibility to simulate non-homogeneous tranmission lines, i.e. transmission lines whose per-unit-length parameters (inductance, capacitance, resistance and conductance) depend on the position along the transmission line. This can be used to simulate lumped non-dispersive elements on the TL.
* The capability to simulate non-uniform transmission lines, i.e. where the spatial step can vary along the line. This allows to use the MTLN solver in non-uniform grids. 

# Multiconductor Transmission Lines #

The MTL method assumes a quasi-TEM propagation applicable in systems that can be decomposed on groups of conductors with constant cross-sections, i.e., translational symmetry. This assumption permits us to relate the voltages and currents on each conductor with per-unit-length parameters [@Paul2007]. For a
lossless line with (N + 1) conductors, these relationships can be expressed in the following form 

$$
\begin{align}
	\partial_z \{V\} &= -[\textrm{R}]\{I\}  - [\textrm{L}] \, \partial_t \{I\} \\
    \partial_z \{I\} &= -[\textrm{G}]\{V\}  - [\textrm{C}] \, \partial_t \{V\} 
\end{align} 
$$

with $\{V\}(z,t) = \left[ V_1 \, V_2 \, ... V_N \right]^T$ being a vector of the voltages with respect to a conductor taken as ground and $\{I\} (z,t) = \left[ I_1 \, I_2 \, ... I_N \right]^T$ is a vector of currents on each of these conductors. 
The matrices $[\textrm{L}]$ and $[\textrm{C}]$ represent the inductance and capacitance per unit length between all conductors.
$\textrm{[R]}$ and $\textrm{[G]}$ represent the losses in the transmission line: resistance (along conductors) and conductance (between conductors) per unit length, respectively.

(1) is discretized using the finite-difference time-domain technique (FDTD) described in [@Paul2007], leading to the following voltage and current update equations:

$$
	\begin{align}
		\{V\}^{n+1}_{k} &= \left(\frac{\Delta z}{\Delta t}[\textrm{C}] + \frac{\Delta z}{2}[\textrm{G}] \right)^{-1} \left[ \left(\frac{\Delta z}{\Delta t}[\textrm{C}] - \frac{\Delta z}{2}[\textrm{G}] \right) \{V\}^{n}_{k} -  \left(\{I\}^{n+1/2}_{k}-\{I\}^{n+1/2}_{k-1}\right) \right] \\
		\{I\}^{n+3/2}_{k} &= \left(\frac{\Delta z}{\Delta t}[\textrm{L}] + \frac{\Delta z}{2}[\textrm{R}] \right)^{-1} \left[ \left(\frac{\Delta z}{\Delta t}[\textrm{L}] - \frac{\Delta z}{2}[\textrm{R}] \right) \{I\}^{n+1/2}_{k}  -  \left(\{V\}^{n+1}_{k+1}-\{V\}^{n+1}_{k}\right) \right]
	\end{align} 
$$
where $\Delta t$ and $\Delta z$ are the time and space steps of the discretization, respectively, $k$ locates the voltage nodes and current segments along the discretized MTL, and $n$ locates the current time on the discretized solution time. 

# Networks #

In a general way we call *networks* the connections at the ends of a transmission line. These connections can be to other tranmission lines (_interconnection_ networks, or _junctions_) or to the reference conductor (_termination_ networks). In its simplest form, a network is a short between elements: a short from a transmission line end to the reference conductor, or a short between the ends of tranmission lines. 

![Left: MTLs connected by interconnection networks (in green) or terminated in termination networks (red). Right: example of possible resistive connections inside a network. ](fig/networks.png)
<figcaption>

Left: MTLs connected by interconnection networks (in green) or terminated in termination networks (red). Right: example of possible resistive connections inside a network. 

</figcaption>

_Solving_ a network means computing the voltages on all the external nodes of the network, which are the start or end nodes of the tranmission lines connected to the network.

Using the state variables analysis, the network can be analyzed solving a system with a number of variables equal to the sum of independent inductors and capacitors.
State equations can in principle solve any network. However, the complexity of the method rises noticeably with increasing number of inductors and capactors. The state equations are very different for different types of connections and the potential for generalization is limited. In the current version of the code, networks composed of either a short, a resistance, a capacitance, or a resistance and a capacitance in parallel, in series with an inductor (RCpLs) can be solved.

An alternative, is to use a circuit simulator to solve the networks. A call to the simulator would be made from within the MTLN solver. At each time step, the currents on the extreme of the bundle connected to the network nodes would be provided as an input, returning the voltages on said nodes.

# Bundles #

A bundle is either a shielded multiconductor tranmission line, or a 

A bundle is divided in domains and levels

A bundle can contain other bundles

# Features #

## Coupling to NgSpice ##

An arbitrary network might include electronic devices in addition to electronic components and sources. With increasing number of components, and specially with the addition of devices, the analytical solution of the network grows dramatically in complexity. Instead, a circuit simulator might be used. [NgSpice](https://ngspice.sourceforge.io/) is an open source spice simulator for electric and electronic circuits. It is an ongoing and maintained project, with a wide user base.

While semba-fdtd is written in Fortran, NgSpice is written in C. Fortran and C are interfactable (functions from the former can be called from the latter, and viceversa), and an interface has been written to allow calls to NgSpice from semba-fdtd.
(ngspice_interface.F90). The module circuit.F90 is a wrapper around the interface functions, allowing for simple calls to NgSpice and providing a series of methods to update currents and voltages between the bundles and the NgSpice circuits. 

NgSpice permits only one simultaneous instance. This means that all networks are part of a single circuit, composed of disconnected subcircuits[^1].



[^1]: Subcircuits understood as parts of the circuit defined in the same instance, not as the .SUBCKT elements.

## Dispersive elements ##

Dispersive elements are those whose properties depend on frequency. The solver works in the time-domain, thus an implementation of the dispersive properties in the time domain is necessary. The implementation proceeds with the following steps:

1. A rational approximation of the frequency-dependent impedance of the dispersive element is written using the [vector-fitting technique](https://www.sintef.no/en/software/vector-fitting/). In this approximation, the impedance can be expressed as a sum of rational functions (poles/residues pairs) and optionally a constant term and a term proportional to the frequency:

    $$Z(s) \simeq d\,+\,se\,+\,\sum\limits_{i=1}\frac{r_i}{s-p_i}\,+\,\sum\limits_{i=1}\frac{r^*_i}{s-p^*_i}$$

    where $s\,=\,j\omega$, and $r_i$ and $p_i$ are the _i_-th residue and pole


2. In the time domain, the products of frequency-dependent terms are expressed as convolutions:

   $$\frac{d\hat{V}(z, w)}{dz} = -\hat{Z}(w)\cdot\hat{I}(z, w) - l\frac{d\hat{I}(z, w)}{dI} \rightarrow 
    \frac{dV(z, t)}{dz} = -Z(t)*I(z, t) - l\frac{dI(z, t)}{dI} $$

3. If the rational approximation of $Z(w)$ is substituted in the convolution:

    $$Z(t)*I(z, t) = \int_o^t Z(\tau) I(z, t-\tau)d\tau = dI + e\frac{\partial I}{\partial t} + \sum\limits_{i} r_i \int_o^t e^{p_i t}I(t-\tau) d\tau$$

4. Following [Sarto](https://ieeexplore.ieee.org/document/809798) the integral in the last term is computed applying the piecewise linear recursive formulation proposed [here](https://ieeexplore.ieee.org/document/391136). Using the time discretization described before, the last term takes the following form:

    $$\sum\limits_{i} r_i \int_o^t e^{p_i t}I(t-\tau) d\tau = \sum\limits_{i} \varphi_i$$

    where

    $$
    \begin{align}
        \varphi^{n}_{i} &= I^{n+3/2}_{k}\,q_{1i} + I^{n+1/2}_{k}\,q_{2i} + \varphi^{n}_{i}\,q_{3i}\\
        q_{1i} &= \frac{\alpha_i}{\beta_i}\,(e^{\beta_i}\,-\,\beta_i\,-\,1)\\
        q_{2i} &= \frac{\alpha_i}{\beta_i}\,(1\,+\,e^{\beta_i}(\beta_i\,-\,1))\\
    	q_{3i} &= e^{\beta_{i}}\textrm{,} \hspace{10pt} \alpha_i\,=\,\frac{r_i}{p_i}\textrm{,} \hspace{10pt} \beta_i\,=\,p_i\,\delta t
    \end{align}
    $$

5. Using this expression for the convolution term, the discretized update relation for the current yields:

	$$I^{n+3/2}_{k} = F_{1}^{-1}[F_{2}I^{n+1/2}_{k} - (V^{n+1}_{k+1}-V^{n+1}_{k}) \Delta z \sum_{i}q_{3i}\varphi^{n-1}_{i}]$$

    where
    
    $$
    \begin{align}
        F_{1} &= ( \frac{\Delta z}{\Delta t}L + \frac{\Delta z}{2}d + \frac{\Delta z}{\Delta t}e + \Delta z \sum_{i}q_{1i}) \\
        F_{2} &= ( \frac{\Delta z}{\Delta t}L - \frac{\Delta z}{2}d + \frac{\Delta z}{\Delta t}e - \Delta z \sum_{i}q_{2i})
    \end{align}
    $$






