Quantum Electrodynamics
=======================

We saw in :ref:`sec_the_canonical_formalism` how Lagrangian (densities) may be guessed, passed to the interaction picture, and verified against the free fields derived from Lorentz invariance and causality in :ref:`sec_quantum_fields_and_antiparticles`. There were, however, cases we couldn't handle back in :ref:`sec_massless_fields`, which include, in particular, massless vector fields. It is this difficulty we will try to handle in this section, which, as it turns out, leads to the most successful quantum field theory to date, namely, the quantum electrodynamics (QED).

The plan to quantize electrodynamics is the following. First, we'll show that the failed attempt to write down a free massless helicity-:math:`1` vector field as :eq:`eq_massless_vector_field_a` can be remedied by the introduction of a gauge symmetry. The introduction of gauge symmetry turns out to also introduce constraints on the canonical variables, which we'll try to fix using Dirac's method discussed in :ref:`sec_constraints_and_dirac_brackets`. Once the constraints have been dealt with, we pass to the interaction picture and write down the Hamiltonian as before. With the propagators calculated using the appropriate commutation relations, the machinery of Feynman diagrams may be applied to calculate S-matrix elements. This is done first theoretically by writing down a recipe, and then applied explicitly to one example -- the famous `Compton scattering <https://en.wikipedia.org/wiki/Compton_scattering>`__.


Gauge Invariance and Action
---------------------------

Recall from :ref:`sec_the_failure_for_vector_fields` that we've tried to write down a quantum vector field :math:`a_{\mu}(x)` for massless helicity-:math:`1` particle in :eq:`eq_massless_vector_field_a`. The issue, which stems from the obviously non-Lorentz-invariant gauge fixing condition :eq:`eq_massless_vector_field_spinor_e_and_p_condition`, is that :math:`a_{\mu}(x)` fails to be covariant in light of :eq:`eq_massless_vector_field_a_conjugation_by_lorentz_transformation`.

From perspective of Lagrangian formalism, however, it's completely fine that the field itself is not covariant, as long as the action is Lorentz invariant. This consideration motivates the introduction of the co-called *mass action* :math:`I_M`, which is required to be invariant under the gauge transformation

.. math:: A_{\mu}(x) \to A_{\mu}(x) + \p_{\mu} \epsilon(x)
    :label: eq_qed_a_field_gauge_transformation

Here we switch to the uppercase :math:`A_{\mu}` to indicate that it is considered to be a Heisenberg-picture field, i.e., a field that knows everything about the interactions. Note also that an explicit formula for :math:`I_M` will only be worked out much later, after we formally apply all sorts of symmetry principles.

Now the variation of :math:`I_M` with respect to :math:`A_{\mu}(x)` can be written as follows

.. math:: \delta I_M = \int d^4 x \frac{\delta I_M}{\delta A_{\mu}(x)} \p_{\mu} \epsilon(x)

It follows from integration-by-parts that in order for :math:`\delta I_M = 0`, one must have

.. math:: \p_{\mu} \frac{\delta I_M}{\delta A_{\mu}(x)} = 0
    :label: eq_qed_4_derivative_of_action_vanishes

Note that this is trivially true if :math:`I_M` is a functional of

.. math:: F_{\mu \nu} = \p_{\mu} A_{\nu} - \p_{\nu} A_{\mu}
    :label: eq_qed_defn_f_field

and its derivatives.

To construct an :math:`I_M` that satisfies :eq:`eq_qed_4_derivative_of_action_vanishes`, we need to couple :math:`A_{\mu}(x)` to a conserved :math:`4`-current :math:`J^{\mu}(x)` in the sense that :math:`\p_{\mu} J^{\mu} = 0`. This will guarantee

.. math:: \frac{\delta I_M}{\delta A_{\mu}(x)} = J^{\mu}(x)
    :label: eq_qed_matter_lagrangian_variational_derivative_in_a_field

and hence :eq:`eq_qed_4_derivative_of_action_vanishes`. Here we've secretly absorbed a potential proportional constant into the definition of :math:`J^{\mu}`, which, in turn, will be absorbed into the definition of charges will be discussed below.

Now recall from :ref:`sec_global_symmetries` that such conserved :math:`4`-current can be obtained by assuming symmetries on :math:`I_M` in the form of :eq:`eq_functional_infinitesimal_variation_of_field`. In QED, we're specifically interested in variations of the following form

.. math:: \delta \Psi_{\ell}(x) = \ifrak \epsilon(x) q_{\ell} \Psi_{\ell}(x)

which is diagonalized into charges :math:`q_{\ell}`. Following discussions in :ref:`Lagrangian density preserving symmetry <list_lagrangian_density_preserving_symmetry>`, if we assume

.. math:: I_M = \int d^4 x~\Lscr_M(\Psi_{\ell}(x), \p_{\mu} \Psi_{\ell}(x))

such that :math:`\Lscr_M` is invariant under :math:`\Psi_{\ell} \to \Psi_{\ell} + \delta \Psi_{\ell}` with a constant :math:`\epsilon(x) = \epsilon`. Then the conserved current, according to :eq:`eq_lagrangian_density_preserving_symmetry_conserved_density`, takes the following form

.. math:: J^{\mu}(x) = -\ifrak \sum_{\ell} \frac{\p \Lscr_M}{\p (\p_{\mu} \Psi_{\ell}(x))} q_{\ell} \Psi_{\ell}(x)
    :label: eq_qed_conserved_current_j

Let

.. math:: Q \coloneqq \int d^3 x~J^0
    :label: eq_qed_defn_charge_operator

be the conserved charge operator. Then it follows from :eq:`eq_lagrangian_formalism_conserved_f_acts_as_symmetry_generator` (modulo an extremely confusing shuffle of notations)

.. math:: \left[ Q, \Psi_{\ell}(x) \right] = -q_{\ell} \Psi_{\ell}(x)

where :math:`\Psi_{\ell}` is now officially one of the canonical variables.

To summarize, we have concluded that the matter action :math:`I_M` is invariant under the joint (local) transformation

.. math::
    :label: eq_qed_gauge_symmetry

    \delta A_{\mu}(x) &= \p_{\mu} \epsilon(x) \\
    \delta \Psi_{\ell}(x) &= \ifrak \epsilon(x) q_{\ell} \Psi_{\ell}(x)

on :math:`A_{\mu}(x)` and :math:`\Psi_{\ell}(x)`. Symmetries like this are known as *gauge symmetries*.

It turns out that in QED, besides the matter action, we also need a light action which takes the following form

.. math:: I_{\gamma} = -\frac{1}{4} \int d^4 x~F_{\mu \nu} F^{\mu \nu}

It's hard to argue logically for why such a term is needed, but it turns out to work, and it appears also in the Lagrangian density :eq:`eq_spin_1_vector_field_lagrangian_density` of massive spin-:math:`1` vector fields.

Combining the mass and light actions together, we see that the field equation is given by the variational principle as follows

.. math:: 0 = \frac{\delta (I_M + I_{\gamma})}{\delta A_{\nu}} = \p_{\mu} F^{\mu \nu} + J^{\nu}

which is recognized as the inhomogeneous Maxwell equations. The homogeneous Maxwell equations

.. math:: \p_{\mu} F_{\nu \kappa} + \p_{\nu} F_{\kappa \mu} + \p_{\kappa} F_{\mu \nu} = 0

follows readily from the definition :eq:`eq_qed_defn_f_field`.

To recap our (or rather, Weinberg's) approach to QED, it starts with a partially successful construction of a massless helicity-:math:`1` vector field :math:`a_{\mu}(x)` from :eq:`eq_massless_vector_field_a`. Its transformation law :eq:`eq_massless_vector_field_a_conjugation_by_lorentz_transformation`, abstracted into :eq:`eq_qed_a_field_gauge_transformation`, necessitates a conserved coupling field :math:`J^{\mu}`, which then comes out of a variational principle. In the end, we've *derived* the gauge symmetry :eq:`eq_qed_gauge_symmetry` as a consequence of the need of a Lorentz invariant Lagrangian density.

It's very interesting that the above argument can be totally reverted, as is done in most textbooks. We outline this reversed argument in the dropdown below.

.. dropdown:: Gauge theoretic approach to the QED Lagrangian
    :animate: fade-in-slide-down
    :icon: unlock

    Suppose the particles of interest have some global internal symmetry (e.g. electric charge) of the following form

    .. math:: \delta \Psi_{\ell}(x) = \ifrak \epsilon q_{\ell} \Psi_{\ell}(x)

    with constant :math:`\epsilon`. The starting point then is the promise to promote these global symmetries to local ones, i.e., gauge symmetries, of the following form

    .. math:: \delta \Psi_{\ell}(x) = \ifrak \epsilon(x) q_{\ell} \Psi_{\ell}(x)

    with function :math:`\epsilon(x)`.

    Now if the Lagrangian involves any derivatives of the fields :math:`\Psi_{\ell}(x)`, which is almost certainly the case in reality, then such promotion cannot be done for free because of the following simple rule of differentiation

    .. math:: \delta \p_{\mu} \Psi_{\ell}(x) = \ifrak \epsilon(x) q_{\ell} \p_{\mu} \Psi_{\ell}(x) + \ifrak q_{\ell} \Psi_{\ell}(x) \p_{\mu} \epsilon(x)
        :label: eq_qed_variation_of_field_derivative

    The cure to this problem is to "invent" a vector field :math:`A_{\mu}(x)` with transformation rule

    .. math::  \delta A_{\mu}(x) = \p_{\mu} \epsilon(x)

    and observe that if we introduce the following `covariant derivative <https://en.wikipedia.org/wiki/Covariant_derivative>`__

    .. math:: D_{\mu} \Psi_{\ell} \coloneqq \p_{\mu} \Psi_{\ell} - \ifrak q_{\ell} A_{\mu} \Psi_{\ell}
        :label: eq_qed_defn_covariant_derivative

    then :eq:`eq_qed_variation_of_field_derivative` can be replaced by the following

    .. math::

        \delta D_{\mu} \Psi_{\ell}(x) &= \delta \p_{\mu} \Psi_{\ell}(x) - \ifrak q_{\ell} \left( A_{\mu}(x) \delta \Psi_{\ell}(x) + \Psi_{\ell}(x) \delta A_{\mu}(x) \right) \\
            &= \ifrak \epsilon(x) q_{\ell} D_{\mu} \Psi_{\ell}(x)

    It follows that the gauge symmetry can be restored if any dependence of the matter Lagrangian :math:`\Lscr_M` on :math:`\p_{\mu} \Psi_{\ell}` actually depends on the covariant derivative :eq:`eq_qed_defn_covariant_derivative`. In other words, we can write

    .. math:: \Lscr_M = \Lscr_M(\Psi, D_{\mu} \Psi)
        :label: eq_qed_matter_lagrangian_in_covariant_derivative

    We can also verify :eq:`eq_qed_matter_lagrangian_variational_derivative_in_a_field` and :eq:`eq_qed_conserved_current_j` in this setup as follows

    .. math::

        \frac{\delta I_M}{\delta A_{\mu}(x)} &= \int d^4 y~\frac{\p \Lscr_M(\Psi_{\ell}(y), D_{\mu} \Psi_{\ell}(y))}{\p A_{\mu}(x)} \\
            &= \sum_{\ell} \frac{\p \Lscr_M}{\p D_{\mu} \Psi_{\ell}} \left( -\ifrak q_{\ell} \Psi_{\ell}(x) \right) \\
            &= -\ifrak \sum_{\ell} \frac{\p \Lscr_M}{\p (\p_{\mu} \Psi_{\ell})} q_{\ell} \Psi_{\ell}(x)

    Comparing with the Lagrangian density of massive spin-:math:`1` vector field :eq:`eq_spin_1_vector_field_lagrangian_density`, we see that the following term

    .. math::  -\frac{1}{2} m^2 A_{\mu} A^{\mu}

    is missing, or :math:`m=0`. Indeed, the appearance of such quadratic term would obviously violate the postulated gauge symmetry. This is really fantastic since we've deduced from just the existence of gauge symmetry that the particle represented by :math:`A_{\mu}`, which turns out to be photon, is massless. This is in fact one the two greatest predictions of QED that photons are massless -- the other one being that they travel at the speed of light, and has been tested against experiment up to great precision.


Constraints and Gauge Conditions
--------------------------------

Recall from the previous section that the QED Lagrangian density splits up into the matter part and the light parts as follows

.. math:: \Lscr = \Lscr_M + \Lscr_{\gamma}

where

.. math:: \Lscr_{\gamma} = -\frac{1}{4} F_{\mu \nu} F^{\mu \nu}
    :label: eq_qed_light_lagrangian_density

and :math:`\Lscr_M` involves both :math:`\Psi`-field and :math:`A`-field. Moreover, in light of :eq:`eq_qed_defn_covariant_derivative` and :eq:`eq_qed_matter_lagrangian_in_covariant_derivative`, there is no involvement of derivatives of :math:`A`-field in :math:`\Lscr_M`. It follows that if we define the field conjugate to :math:`A_{\mu}` as follows

.. math:: \Pi^{\mu} \coloneqq \frac{\p \Lscr}{\p \dot{A}_{\mu}}
    :label: eq_qed_defn_conjugate_field_pi

then just as in the theory of (massive) spin-:math:`1` vector fields, we get the following primary constraint

.. math:: \Pi^0 = \frac{\p \Lscr_{\gamma}}{\p \dot{A}_{\mu}} = \frac{\p \Lscr_{\gamma}}{\p F_{00}} = 0
    :label: eq_qed_primary_constraint

Moreover we use the Euler-Lagrange equation to find the secondary constraint as follows

.. math:: \p_i \Pi^i = -\p_i \frac{\p \Lscr}{\p F_{i0}} = -\frac{\p \Lscr}{\p A_0} = -J^0
    :label: eq_qed_secondary_constraint

where the middle equality is one of the Euler-Lagrange equations :eq:`eq_euler_lagrange` with the time-derivative term dropped out due to the primary constraint, the equality on the left is the definition :eq:`eq_qed_defn_conjugate_field_pi`, and the equality on the right follows from :eq:`eq_qed_matter_lagrangian_variational_derivative_in_a_field`.

Using :eq:`eq_qed_conserved_current_j` and letting :math:`Q \coloneqq \Psi` (not to be confused with the charge operator :eq:`eq_qed_defn_charge_operator` defined in the previous section) and :math:`P \coloneqq \delta \Lscr / \delta \dot{\Psi}` be the conjugate canonical variables, we can write the charge density :math:`J^0` as follows

.. math:: J^0 = -\ifrak \sum_n \frac{\p \Lscr}{\p \dot{\Psi}_n} q_n \Psi_n = -\ifrak \sum_n P^n q_n Q_n
    :label: eq_qed_charge_density_in_canonical_q_and_p

Together with :eq:`eq_qed_primary_constraint` and :eq:`eq_qed_secondary_constraint`, we have obtained two constraints

.. math::

    \chi_1 &\coloneqq \Pi^0 = 0 \\
    \chi_2 &\coloneqq \p_i \Pi^i - \ifrak \sum_n P^n q_n Q_n = 0

Obviously the Poisson bracket :math:`[\chi_1, \chi_2]_P = 0`. Hence the constraints are of first class as discussed in :ref:`sec_first_class_constraints`. Moreover, one cannot hope to solve explicitly for :math:`A_0`, the conjugate of :math:`\Pi^0`, in terms of the other canonical variables, as we did in the massive case :eq:`eq_spin_1_vector_field_heisenberg_v0`, due to the presence of gauge symmetry :eq:`eq_qed_gauge_symmetry`.

One way, which is unfortunately not Lorentz invariant, to solve this problem is to "fix the gauge" by imposing an artificial condition on the :math:`A`-field so that :math:`A_0` may be solved explicitly. [#brst_quantization]_ It turns out that there are many such conditions that find their use cases under various circumstances. Some of the most common `gauge-fixing conditions <https://en.wikipedia.org/wiki/Gauge_fixing>`__ are listed below

+-------------------------+-------------------------------+
| Name                    | Condition                     |
+=========================+===============================+
| Lorenz gauge [#lorenz]_ | :math:`\p_{\mu} A^{\mu} = 0`  |
+-------------------------+-------------------------------+
| Coulomb gauge           | :math:`\nabla \cdot \Abf = 0` |
+-------------------------+-------------------------------+
| Temporal gauge          | :math:`A^0 = 0`               |
+-------------------------+-------------------------------+
| Axial gauge             | :math:`A^3 = 0`               |
+-------------------------+-------------------------------+

.. important::

    We will use the Coulomb gauge unless otherwise specified throughout this chapter. The reason, according to J. Schwinger, is that it is in this gauge that photons come out as helicity-:math:`1` particles.

Two things remain to be settled. One is to argue that the :math:`A`-field can always be put into the Coulomb gauge, and the other is to solve for :math:`A^0` in terms of the other fields.

To address the first point, it suffices, according to :eq:`eq_qed_gauge_symmetry`, to show that for any :math:`A = (A_0, \Abf)`, there exists a function :math:`\lambda(\xbf)` such that

.. math:: \nabla \cdot (\Abf + \nabla \cdot \lambda) = 0 \iff \nabla^2 \lambda = -\nabla \cdot \Abf

The equation on the right is known as `Poisson's equation <https://en.wikipedia.org/wiki/Poisson%27s_equation>`__ and can be solved using Green's function. In fact, the same equation shows up again when solving for :math:`A^0`, which we now explain.

Combining the secondary constraint :eq:`eq_qed_secondary_constraint` (cf. :eq:`eq_qed_light_lagrangian_density`) with the Coulomb gauge condition, we have

.. math:: J^0 = \p_i \frac{\delta \Lscr}{\p F_{i0}} = \p_i \frac{\delta \Lscr_{\gamma}}{\p F_{i0}} = -\p_i F^{i0} = -\nabla^2 A^0
    :label: eq_qed_poisson_equation_j_and_a

This equation can be solved explicitly [#solve_poisson_equation]_ by

.. math:: A^0(t, \xbf) = \int d^3 y~\frac{J^0(t, \ybf)}{4\pi |\xbf - \ybf|}
    :label: eq_qed_explicit_solution_of_a0

where the charge density :math:`J^0` can be further written in the form of :eq:`eq_qed_charge_density_in_canonical_q_and_p`.

Now :eq:`eq_qed_primary_constraint` and :eq:`eq_qed_explicit_solution_of_a0` allow us to get rid of the constrained canonical variables :math:`A_0` and :math:`\Pi^0`, which, as we will show in the next section, removes the first class constraints.


Quantization in Coulomb Gauge
-----------------------------

In the previous section, we imposed a Coulomb gauge condition :math:`\nabla \cdot \Abf = 0` and used it to eliminate the constrained variables :math:`A_0` and :math:`\Pi^0`. Now we're facing a new pair (of families parametrized by spatial coordinates :math:`\xbf`) of constraints, listed as follows

.. math::
    :label: eq_qed_constraints_in_coulomb_gauge

    \chi_{1 \xbf} &\coloneqq \p_i A^i(\xbf) = 0 \\
    \chi_{2 \xbf} &\coloneqq \p_i \Pi^i(\xbf) + J^0(\xbf) = 0

where we also remember that :math:`J^0` can be further expressed in terms of canonical variables as in :eq:`eq_qed_charge_density_in_canonical_q_and_p`.

The :math:`C` matrix as defined by :eq:`eq_constraints_c_matrix` is given in this case by

.. math::
    C = \begin{bmatrix*}
        0 & -\nabla^2 \delta^3(\xbf - \ybf) \\
        \nabla^2 \delta^3(\xbf - \ybf) & 0
    \end{bmatrix*}

where, for example, the upper-right entry :math:`C_{1\xbf, 2\ybf}` may be calculated as follows

.. math::

    C_{1\xbf, 2\ybf} &= [\chi_{1\xbf}, \chi_{2\ybf}]_P \\
        &= \int d^3 \zbf \left( \frac{\p (\p_i A^i(\xbf))}{\p A^k(\zbf)} \frac{\p (\p_j \Pi^j(\ybf) + J^0(\ybf))}{\p \Pi_k(\zbf)} - \xbf \leftrightarrow \ybf \right) \\
        &= \int d^3 \zbf~\p_i \left( \frac{\p A^i(\xbf)}{\p A^k(\zbf)} \right) \p_j \left( \frac{\p \Pi^j(\ybf)}{\p \Pi_k(\zbf)} \right) \\
        &= \int d^3 \zbf~\delta^i_k \delta^{jk} \p_i \delta^3(\xbf - \zbf) \p_j \delta^3(\ybf - \zbf) \\
        &= -\nabla^2 \delta^3(\xbf - \ybf)

Clearly :math:`C` is non-singular, and therefore the constraints :eq:`eq_qed_constraints_in_coulomb_gauge` are of second class. To apply Dirac's method, we need the inverse matrix :math:`C^{-1}` given by

.. math::

    C^{-1} = \begin{bmatrix*}
        0 & -\frac{1}{4\pi |\xbf - \ybf|} \\
        \frac{1}{4\pi |\xbf - \ybf|} & 0
    \end{bmatrix*}

Indeed, one easily verifies that :math:`C C^{-1} = 1` by, for example, the following calculation

.. math::

    \int d^3 \zbf~C_{1\xbf, 2\zbf} (C^{-1})_{2\zbf, 1\ybf}
        = -\int d^3 \zbf~\frac{\nabla^2 \delta^3(\xbf - \zbf)}{4\pi |\zbf - \ybf|}
        = \delta^3(\xbf - \ybf)

where the last equality follows again from the solution to Poisson's equation.

Now we apply Dirac's recipe :eq:`eq_canonical_bracket_as_dirac_bracket` and :eq:`eq_defn_dirac_bracket` to compute the commutators as follows

.. math::

    [A_i(\xbf), \Pi_j(\ybf)] &= \ifrak [A_i(\xbf), \Pi_j(\ybf)]_P - \ifrak \int d^3 \zbf \int d^3 \wbf [A_i(\xbf), \chi_{2\zbf}] (C^{-1})_{2\zbf, 1\wbf} [\chi_{1\wbf}, \Pi_j(\ybf)] \\
        &= \ifrak \delta_{ij} \delta^3(\xbf - \ybf) + \ifrak \int d^3 \zbf \int d^3 \wbf \left( \p_i \delta^3(\xbf - \zbf) \frac{1}{4\pi |\zbf - \wbf|} \p_j \delta^3(\wbf - \ybf) \right) \\
        &= \ifrak \delta_{ij} \delta^3(\xbf - \ybf) + \ifrak \frac{\p^2}{\p x_i \p x_j} \left( \frac{1}{4\pi |\xbf - \ybf|} \right) \\
    [A_i(\xbf), A_j(\ybf)] &= [\Pi_i(\xbf), \Pi_j(\ybf)] = 0

It's straightforward to check that they are indeed compatible with the constraints :eq:`eq_qed_constraints_in_coulomb_gauge`.

.. dropdown:: A formula for :math:`\Pi^i(\xbf)`
    :animate: fade-in-slide-down
    :icon: unlock

    It's a real concern that the Coulomb gauge :math:`\nabla \cdot \Abf = 0` may have spoiled the very definition :eq:`eq_qed_defn_conjugate_field_pi` of :math:`\Pi^i`. To settle this, let's go back to the root of Lagrangian formalism and redefine :math:`\Pi^i` using :eq:`eq_general_lagrangian_conjugate_pi` as follows

    .. math:: \Pi^i \coloneqq \frac{\delta L}{\delta \dot{A}_i}

    Now if this were well-defined, then for any infinitesimal variation :math:`\delta \dot{\Abf}`, there should exists a unique :math:`\bm{\Pscr}` such that

    .. math:: \delta L = \int d^3 \xbf~\bm{\Pscr} \cdot \delta \dot{\Abf}
        :label: eq_qed_lagrangian_variation_by_a_dot

    But this is not true since a substitution :math:`\bm{\Pscr} \to \bm{\Pscr} + \nabla f` will also satisfy :eq:`eq_qed_lagrangian_variation_by_a_dot` for any (scalar) function :math:`f(\xbf)` due to the Coulomb gauge condition.

    Baring this difficulty in mind, let's just evaluate the potentially ill-defined :math:`\Pi^i` anyway

    .. math:: \Pi^i = \frac{\delta L}{\delta \dot{A}_i} = \frac{\delta \Lscr_{\gamma}}{\delta \dot{A}_i} = \dot{A}^i(\xbf) + \frac{\p A^0(\xbf)}{\p x_i}
        :label: eq_qed_pi_in_terms_of_a

    The miracle is that this :math:`\Pi^i` satisfies exactly the second constraint in :eq:`eq_qed_constraints_in_coulomb_gauge` (cf. :eq:`eq_qed_poisson_equation_j_and_a`), and hence no ambiguity like :math:`\nabla f` can appear. In other words :eq:`eq_qed_pi_in_terms_of_a` is the correct formula for :math:`\Pi^i`.

Next let's bring in the matter fields introduced in :eq:`eq_qed_charge_density_in_canonical_q_and_p`. Due to the appearance of :math:`J^0` in the constraints :eq:`eq_qed_constraints_in_coulomb_gauge`, one can check that while :math:`\Abf` commutes with the matter fields, :math:`\bm{\Pi}` doesn't. Indeed, if :math:`F` is any function of the matter fields :math:`Q_n` and :math:`P^n`, then we can calculate the Dirac bracket as follows

.. math::

    [F, \bm{\Pi}(\zbf)]_D &= -\int d^3 \xbf \int d^3 \ybf~[F, \chi_{2\xbf}]_P \frac{1}{4\pi |\xbf - \ybf|} [\chi_{1\ybf}, \bm{\Pi}(\zbf)]_P \\
        &= -\int d^3 \xbf \int d^3 \ybf~[F, J^0(\xbf)]_P \frac{1}{4\pi |\xbf - \ybf|} \nabla \delta^3(\ybf - \zbf) \\
        &= -\int d^3 \ybf [F, A^0(\ybf)]_P \nabla \delta^3(\ybf - \zbf) \\
        &= [F, \nabla A^0(\zbf)]_P = [F, \nabla A^0(\zbf)]_D

where we've used :eq:`eq_qed_poisson_equation_j_and_a` in the third equality.

This calculation motivates the following definition a new field

.. math:: \bm{\Pi}_{\bot} \coloneqq \bm{\Pi} - \nabla A^0 = \dot{\Abf}
    :label: eq_qed_defn_pi_bot

where the last equality follows from :eq:`eq_qed_pi_in_terms_of_a`, such that

.. math:: [F, \bm{\Pi}_{\bot}] = 0

for any function :math:`F` of matter fields. Moreover, one can verify that :math:`\bm{\Pi}_{\bot}` satisfies exactly the same commutation relations with :math:`\Abf` just as :math:`\bm{\Pi}` which we repeat as follows

.. math::
    :label: eq_qed_a_pi_bot_commutation_relations

    [A_i(\xbf), (\Pi_{\bot})_j(\ybf)] &= \ifrak \delta_{ij} \delta^3(\xbf - \ybf) + \ifrak \frac{\p^2}{\p x_i \p x_j} \left( \frac{1}{4\pi |\xbf - \ybf|} \right) \\
    [A_i(\xbf), A_j(\ybf)] &= [(\Pi_{\bot})_i(\xbf), (\Pi_{\bot})_j(\ybf)] = 0

Finally, note that the second constraint in :eq:`eq_qed_constraints_in_coulomb_gauge` can now be written as follows

.. math:: \nabla \cdot \bm{\Pi}_{\bot} = 0
    :label: eq_qed_pi_bot_is_divergence_free

With all the preparations above, let's first write down the Hamiltonian in a rather general form

.. math:: H = \int d^3 x \left( (\Pi_{\bot})_i \dot{A}^i + P_n \dot{Q}^n - \Lscr \right)
    :label: eq_qed_hamiltonian_general_form

despite the existence of constraints. For the rest of this chapter, we'll be considering Lagrangian density of the following form

.. math:: \Lscr = -\frac{1}{4} F_{\mu\nu} F^{\mu\nu} + J^{\mu} A_{\mu} + \Lscr_{\text{matter}}
    :label: eq_qed_semi_concrete_lagrangian_with_matter

where :math:`J^{\mu}` is the charge density as before, but :math:`\Lscr_{\text{matter}}` is different from :math:`\Lscr_M` considered in :eq:`eq_qed_matter_lagrangian_in_covariant_derivative` since it contains only the terms that don't interact with :math:`A`, i.e., only the :math:`Q` and :math:`P` fields.

.. note::

    Both :math:`J^{\mu}` and :math:`\Lscr_{\text{matter}}` will be further specified later for a theory specific to spin-:math:`1/2` fermions, e.g., electrons. In fact, according to [Wei95]_ (page 349), the QED for spinless particles would require a more complicated Lagrangian density than :eq:`eq_qed_semi_concrete_lagrangian_with_matter`, but it wasn't mentioned in the book which physical particle(s) such theory would describe.

Plugging :eq:`eq_qed_semi_concrete_lagrangian_with_matter` into :eq:`eq_qed_hamiltonian_general_form`, and making use of :eq:`eq_qed_defn_pi_bot` and :eq:`eq_qed_pi_in_terms_of_a`, we can write the Hamiltonian, again with a separation between matter and light, as follows

.. math::
    :label: eq_qed_hamiltonian_vector_form_raw

    H = \int d^3 x \left( \bm{\Pi}_{\bot}^2 + \frac{1}{2} (\nabla \times \Abf)^2 - \frac{1}{2} (\bm{\Pi}_{\bot} + \nabla A^0)^2 - \Jbf \cdot \Abf + J^0 A^0 \right) + H_{\text{matter}}

where

.. math:: H_{\text{matter}} \coloneqq \int d^3 x \left( P_n \dot{Q}^n - \Lscr_{\text{matter}} \right)

Expanding out :math:`(\bm{\Pi}_{\bot} + \nabla A^0)^2` in :eq:`eq_qed_hamiltonian_vector_form_raw` and using integration-by-parts together with :eq:`eq_qed_pi_bot_is_divergence_free` and :eq:`eq_qed_poisson_equation_j_and_a`, we can further rewrite :eq:`eq_qed_hamiltonian_vector_form_raw` as follows

.. math::
    :label: eq_qed_hamiltonian_vector_form

    H = \int d^3 x \left( \frac{1}{2} \bm{\Pi}_{\bot}^2 + \frac{1}{2} (\nabla \times \Abf)^2 - \Jbf \cdot \Abf + \frac{1}{2} J^0 A^0 \right) + H_{\text{matter}}

.. note::

    The term :math:`\tfrac{1}{2} J^0 A^0` in :eq:`eq_qed_hamiltonian_vector_form` gives nothing but the Coulomb energy as the following calculation (cf. :eq:`eq_qed_explicit_solution_of_a0`) shows

    .. math::
        :label: eq_qed_defn_coulomb_energy

        V_{\text{Coul}} \coloneqq \frac{1}{2} \int d^3 x~J^0 A^0 = \frac{1}{2} \int d^3 x \int d^3 y~\frac{J^0(\xbf) J^0(\ybf)}{4\pi |\xbf - \ybf|}

    where we've suppressed the :math:`t`-dependence as usual in the Lagrangian formalism.


Electrodynamics in the Interaction Picture
------------------------------------------

As before, let's split the Hamiltonian :eq:`eq_qed_hamiltonian_vector_form` into the free part and the interaction part as follows

.. math::

    H &= H_0 + V \\
    H_0 &= \int d^3 x \left( \frac{1}{2} \bm{\Pi}_{\bot}^2 + \frac{1}{2} (\nabla \times \Abf)^2 \right) + H_{\text{matter}, 0} \\
    V &= -\int d^3 x~\Jbf \cdot \Abf + V_{\text{Coul}} + V_{\text{matter}}

where :math:`V_{\text{Coul}}` is defined by :eq:`eq_qed_defn_coulomb_energy` and :math:`H_{\text{matter}} = H_{\text{matter}, 0} + V_{\text{matter}}` is the splitting of the matter Hamiltonian into the free and interaction parts.

Before passing to the interaction picture, let's make one more (potentially confusing) change of notation:

.. warning::

    For the rest of this chapter, we'll simply write :math:`\bm{\Pi}` in place of :math:`\bm{\Pi}_{\bot}`, which is *not* the original :math:`\bm{\Pi}` as in :eq:`eq_qed_defn_pi_bot`. Since :math:`\bm{\Pi}` and :math:`\bm{\Pi}_{\bot}` satisfy the same commutation relations with the other fields such as :math:`\Abf` and the matter fields, the only thing we need to keep in mind is that we should use the constraint :eq:`eq_qed_pi_bot_is_divergence_free` rather than the second equation in :eq:`eq_qed_constraints_in_coulomb_gauge`.

Finally we can introduce the interaction-picture fields :math:`\abf, \bm{\pi}, q, p` corresponding to the Heisenberg-picture fields :math:`\Abf, \bm{\Pi}, Q, P`, respectively. We'll focus in this section on the :math:`\abf` and :math:`\bm{\pi}` fields. In fact, the letters :math:`q, p` will soon be occupied by something completely different, namely, the momentum-space coordinates. It turns out in the theory for spin-:math:`1/2` particles, which will be completely specified in the next section, the matters fields will be named by another letter, so we'll not run into conflicts of notations.

The interaction-picture free Hamiltonian :math:`H_0` takes the following form

.. math:: H_0 = \int d^3 x \left( \frac{1}{2} \bm{\pi}^2 + \frac{1}{2} (\nabla \times \abf)^2 \right) + H_{\text{matter}, 0}(t)
    :label: eq_qed_interaction_picture_free_hamiltonian

The commutation relations between :math:`\abf` and :math:`\bm{\pi}` follow directly from :eq:`eq_qed_a_pi_bot_commutation_relations` and the definition of interaction-picture fields :eq:`eq_defn_interaction_perturbation_term`

.. math::
    :label: eq_qed_interaction_picture_a_pi_commutation_relations

    [a_i(t, \xbf), \pi_j(t, \ybf)] &= \ifrak \delta_{ij} \delta^3(\xbf - \ybf) + \ifrak \frac{\p^2}{\p x_i \p x_j} \frac{1}{4\pi |\xbf - \ybf|} \\
    [a_i(t, \xbf), a_j(t, \ybf)] &= [\pi_i(t, \xbf), \pi_j(t, \ybf)] = 0

Moreover, they satisfy the following constraints

.. math::
    :label: eq_qed_a_and_pi_divergence_free

    \nabla \cdot \abf &= 0 \\
    \nabla \cdot \bm{\pi} &= 0

due to the Coulomb gauge condition and :eq:`eq_qed_pi_bot_is_divergence_free`, respectively.

For later reference, let's note that the Coulomb interaction :math:`V_{\text{Coul}}` defined by :eq:`eq_qed_defn_coulomb_energy` becomes

.. math:: V_{\text{Coul}}(t) = \frac{1}{2} \int d^3 x \int d^3 y~\frac{j^0(t, \xbf) j^0(t, \ybf)}{4\pi |\xbf - \ybf|}
    :label: eq_qed_interaction_picture_coulomb_interaction

in the interaction picture, and the interaction part of the full Hamiltonian becomes

.. math:: V(t) = -\int d^3 x~j^{\mu}(t, \xbf) a_{\mu}(t, \xbf) + V_{\text{Coul}}(t) + V_{\text{matter}}(t)
    :label: eq_qed_interaction_picture_interaction_hamiltonian

To derive the field equations, we resort to Hamilton's equations :eq:`eq_free_field_hamilton_equation_q_and_p_dot` as follows

.. math::

    \dot{a}_i(t,\xbf) &= \ifrak [H_0, a_i(t, \xbf)] \\
        &= \ifrak \int d^3 y~\left[ \pi_j(t, \ybf), a_i(t, \xbf) \right] \pi^j(t, \ybf) \\
        &= \int d^3 y \left( \delta_{ij} \delta^3(\xbf - \ybf) + \frac{\p^2}{\p x_i \p x_j} \frac{1}{4\pi |\xbf - \ybf|} \right) \pi^j(t, \ybf) \\
        &= \pi_i(t, \xbf) - \int d^3 x \left( \frac{\p^2}{\p x_i \p y_j} \frac{1}{4\pi |\xbf - \ybf|} \right) \pi^j(t, \ybf) \\
        &= \pi_i(t, \xbf) \\\\
    \dot{\pi}_i(t, \xbf) &= \ifrak [H_0, \pi_i(t, \xbf)] \\
        &= \ifrak \int d^3 y~[a_j(t, \ybf), \pi_i(t, \xbf)] (\nabla \times \nabla \times \abf(t, \ybf))^j \\
        &= -\int d^3 y \left( \delta_{ij} \delta^3(\xbf - \ybf) + \frac{\p^2}{\p x_i \p x_j} \frac{1}{4\pi |\xbf - \ybf|} \right)
            \left( \nabla (\nabla \cdot \abf(t, \ybf)) - \nabla^2 \abf(t, \ybf) \right)^j \\
        &= \nabla^2 a_i(t, \xbf)

These two equations combined together give the familiar wave equation

.. math:: \square \abf = 0
    :label: eq_qed_a_3_vector_satisfies_wave_equation

To upgrade :math:`\abf` to a :math:`4`-vector field, we must set

.. math:: a_0 = 0
    :label: eq_qed_a0_vanishes

due to the assumption that :math:`a_0` shouldn't depend on the charge density, and it must vanish when charge density vanishes according to :eq:`eq_qed_explicit_solution_of_a0`.

The general solutions to :eq:`eq_qed_a_3_vector_satisfies_wave_equation` and :eq:`eq_qed_a0_vanishes`, under the constraint :eq:`eq_qed_a_and_pi_divergence_free`, can be written as follows

.. math::
    :label: eq_qed_a_field_general_solution

    a_{\mu}(x) = (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \left(
        e^{\ifrak p \cdot x} e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) + e^{-\ifrak p \cdot x} e_{\mu}^{\ast}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma)
    \right)

where :math:`p_0 = |\pbf|` and :math:`e_{\mu}(\pbf, \pm 1)` are two independent "polarization vectors" satisfying the following familiar conditions (cf. :eq:`eq_massless_vector_field_spinor_e_and_p_condition`)

.. math::

    e_0(\pbf, \pm 1) &= 0 \\
    \pbf \cdot \ebf(\pbf, \pm 1) &= 0

It follows that :math:`\ebf(\pbf, \sigma)` can be normalized as follows

.. math:: \sum_{\sigma = \pm 1} e_i(\pbf, \sigma) e_j^{\ast}(\pbf, \sigma) = \delta_{ij} - \frac{p_i p_j}{|\pbf|^2}

Without working out the calculations, we claim, following [Wei95]_ page 352 --353, that the commutation relations :eq:`eq_qed_interaction_picture_a_pi_commutation_relations` are satisfied if the operators :math:`a(\pbf, \sigma)` and :math:`a^{\dagger}(\pbf, \sigma)` satisfy

.. math::
    :label: eq_qed_a_operator_commutation_relations

    \left[ a(\pbf, \sigma), a^{\dagger}(\pbf', \sigma') \right] &= \delta^3(\pbf - \pbf') \delta_{\sigma \sigma'} \\
    \left[ a(\pbf, \sigma), a(\pbf', \sigma') \right] &= \left[ a^{\dagger}(\pbf, \sigma), a^{\dagger}(\pbf', \sigma') \right] = 0

Moreover, the free-photon Hamiltonian, i.e., :eq:`eq_qed_interaction_picture_free_hamiltonian` without the matter term, can be written as follows

.. math::

    H_{\gamma, 0} &= \frac{1}{2} \int d^3 x \left( \bm{\pi}^2 + (\nabla \times \abf)^2 \right) \\
        &= \frac{1}{2} \int d^3 p~p_0 \sum_{\sigma = \pm 1} \left[ a(\pbf, \sigma), a^{\dagger}(\pbf, \sigma) \right]_+ \\
        &= \int d^3 p~p_0 \sum_{\sigma = \pm 1} \left( a^{\dagger}(\pbf, \sigma) a(\pbf, \sigma) + \frac{1}{2} \delta^3(0) \right)

where the last expression contains one of the infinities in QED.


The Photon Propagator
---------------------

As explained in :ref:`sec_the_feynman_rules`, to calculate the S-matrix using Feynman diagrams, one must calculate the propagators defined in :eq:`eq_feynman_rule_propagator`. Using :eq:`eq_qed_a_field_general_solution` and the commutation relations :eq:`eq_qed_a_operator_commutation_relations`, the photon propagator can be calculated as follows

.. math::

    \underbracket{a_{\mu}(x) a_{\nu}^{\dagger}(y)}
        &= \theta(x_0 - y_0) (2\pi)^{-3}
            \int \frac{d^3 p~d^3 p'}{2\sqrt{p_0 p'_0}}
            \sum_{\sigma, \sigma'} e^{\ifrak (p \cdot x - p' \cdot y)} e_{\mu} e_{\nu}^{\ast}
            \left[ a, a^{\dagger} \right]
            + x \leftrightarrow y \\
        &= (2\pi)^{-3} \int \frac{d^3 p}{2p_0}
            \left( \sum_{\sigma} e_{\mu}(\pbf, \sigma) e_{\nu}^{\ast}(\pbf, \sigma) \right)
            \left( e^{\ifrak p \cdot (x-y)} \theta(x_0 - y_0) + e^{\ifrak p \cdot (y-x) \theta(y_0 - x_0)} \right) \\
        &\eqqcolon -\ifrak \Delta_{\mu \nu}(x - y)

Since the spinor sum

.. math::
    :label: eq_qed_photon_spinor_sum

    P_{\mu \nu}(\pbf) \coloneqq \sum_{\sigma} e_{\mu}(\pbf, \sigma) e_{\nu}^{\ast}(\pbf, \sigma) = \begin{cases}
        \delta_{\mu \nu} - p_{\mu} p_{\nu} / |\pbf|^2 & \text{ if } \mu\nu \neq 0 \\
        0 & \text{ otherwise}
    \end{cases}

doesn't depend on :math:`p_0`, it follows from :eq:`eq_spinor_sum_momentum_space_linear_extension` and :eq:`eq_propagator_as_momentum_space_integral_linear` and the massless condition that

.. math:: \Delta_{\mu \nu}(x - y) = (2\pi)^{-4} \int d^4 q~\frac{P_{\mu \nu}(\qbf)}{q^2 - \ifrak \epsilon} e^{\ifrak q \cdot (x - y)}
    :label: eq_qed_photon_propagator_non_covariant

where :math:`q^2 = q_0^2 - |\qbf|^2` as usual but :math:`q` is not constrained to the mass shell.

Using the momentum-space Feynman rules derived in :ref:`sec_feynman_rules_in_momentum_space`, one gets a contribution of

.. math:: \frac{-\ifrak}{(2\pi)^4} \frac{P_{\mu \nu}(\qbf)}{q^2 - \ifrak \epsilon}
    :label: eq_qed_photon_line_contribution_non_covariant

for each internal photon line in a Feynman diagram.

Note that the photon line contribution given by :eq:`eq_qed_photon_line_contribution_non_covariant` is not Lorentz covariant, which is ultimately a consequence of the :math:`A`-field :math:`a_{\mu}` (cf. :eq:`eq_qed_a_field_general_solution`) not being Lorentz covariant. It turns out, rather miraculously, that such Lorentz non-covariance can be countered by yet another Lorentz non-covariant term in the interaction Hamiltonian, namely, the Coulomb interaction :eq:`eq_qed_interaction_picture_coulomb_interaction`

The heuristic for such cancellation starts by rewriting the spinor sum :eq:`eq_qed_photon_spinor_sum` as follows

.. math:: P_{\mu\nu}(\qbf) = \eta_{\mu\nu} - \frac{q^2 n_{\mu} n_{\nu} - q_0 q_{\mu} n_{\nu} - q_0 q_{\nu} n_{\mu} + q_{\mu} q_{\nu}}{|\qbf|^2}

where :math:`n = (1, 0, 0, 0)`. Here :math:`q_0` is completely arbitrary, but its value will be fixed in Feynman diagram evaluations, to be discussed in the next section, by momentum conservation at each vertex.

The point is that the last three terms in the nominator are proportional to :math:`q_{\mu}` or :math:`q_{\nu}`, which means that they will not contribute to the S-matrix due to the coupling :math:`j^{\mu} a_{\mu}` in the interaction :eq:`eq_qed_interaction_picture_interaction_hamiltonian` and the conservation law :math:`\p_{\mu} j^{\mu} = 0`. The remaining :math:`-q^2 n_{\mu} n_{\nu} / |\qbf|^2`, which is non-vanishing only if :math:`\mu = \nu = 0`, can be plugged into :eq:`eq_qed_photon_propagator_non_covariant`, and then follow the Feynman rules to produce the following contribution to the S-matrix

.. math::

    & \frac{\ifrak}{2} \int d^4 x \int d^4 y
            \left( -\ifrak j^0(x) \right)
            \left( -\ifrak j^0(y) \right)
            \frac{-\ifrak}{(2\pi)^4} \int \frac{d^4 q}{|\qbf|^2} e^{\ifrak q \cdot (x-y)} \\
        &\quad = -\frac{1}{2} \int d^4 x \int d^4 y~\frac{j^0(x) j^0(y)}{(2\pi)^3} \delta(x_0 - y_0) \int \frac{d^3 q}{|\qbf|^2} e^{\ifrak \qbf \cdot (\xbf - \ybf)} \\
        &\quad = -\frac{1}{2} \int d^3 x \int d^3 y~\frac{j^0(t, \xbf) j^0(t, \ybf)}{4\pi |\xbf - \ybf|}

which is countered by :math:`V_{\text{Coul}}` (cf. :eq:`eq_qed_interaction_picture_coulomb_interaction`).

It means that effectively, one can replace :math:`\Delta_{\mu \nu}(x - y)` with the following

.. math:: \Delta_{\mu\nu}^{\text{eff}}(x - y) = (2\pi)^{-4} \int d^4 q~\frac{\eta_{\mu \nu}}{q^2 - \ifrak \epsilon} e^{\ifrak q \cdot (x - y)}

as long as one forgets about the (instantaneous) Coulomb interaction term :math:`V_{\text{Coul}}(t)` in :eq:`eq_qed_interaction_picture_interaction_hamiltonian`. The corresponding photon line contribution :eq:`eq_qed_photon_line_contribution_non_covariant` then takes the following form

.. math:: \frac{-\ifrak}{(2\pi)^4} \frac{\eta_{\mu\nu}}{q^2 - \ifrak \epsilon}
    :label: eq_qed_photon_line_contribution_covariant


Feynman Rules for QED
---------------------

Having settle the photon field in the previous sections, we're now ready to specify the matter field in the Lagrangian density :eq:`eq_qed_semi_concrete_lagrangian_with_matter`. We'll be considering the theory describing the interaction between a massive spin-:math:`1/2` fermion with charge :math:`q = -e`, [#letter_e]_ e.g., an electron, and photon. In this case, the Lagrangian density takes the following form

.. math:: \Lscr = -\frac{1}{4} F^{\mu\nu} F_{\mu\nu} - \bar{\Psi} \left( \gamma^{\mu} (\p_{\mu} + \ifrak e A_{\mu}) + m \right) \Psi

It follows from :eq:`eq_qed_matter_lagrangian_variational_derivative_in_a_field` that the charge density :math:`J^{\mu}` takes the following form

.. math:: J^{\mu} = \frac{\p \Lscr}{\p A_{\mu}} = -\ifrak e \bar{\Psi} \gamma^{\mu} \Psi

The interaction :eq:`eq_qed_interaction_picture_interaction_hamiltonian` now takes the following form

.. math:: V(t) = \ifrak e \int d^3 x~\bar{\psi}(t, \xbf) \gamma^{\mu} \psi(t, \xbf) a_{\mu}(t, \xbf) + V_{\text{Coul}}(t)
    :label: eq_qed_concrete_interaction_massive_spin_half_fermion

where we remember that :math:`V_{\text{Coul}}` becomes irrelevant if we take the photon line contribution to the Feynman diagram to be :eq:`eq_qed_photon_line_contribution_covariant`.

Following the general recipe described in :ref:`sec_feynman_rules_in_momentum_space`, let's spell out the key points, while ignoring routines, in constructing and evaluating Feynman diagrams in QED as follows.

First of all, a Feynman diagram consists of electron lines, photon lines, and :math:`3`-valent vertices, each of which has one incoming electron line, one outgoing electron line, and one photon line attached. Electron lines that flow backward in time may also be called positron lines. To each internal line an off-mass-shell :math:`4`-momentum is attached. To each external line an on-mass-shell :math:`4`-momentum as well as a spin :math:`z`-component or helicity, depending on whether it's an electron line or a photon line, are attached.

Next, to evaluate a Feynman diagram, we associate factors to each vertex, external line, and internal line as follows.

Vertices
    To each vertex, we assign a :math:`\gamma`-matrix (cf. :eq:`eq_dirac_field_defn_gamma_matrices`) index :math:`\alpha` to the incoming electron line, and a similar index :math:`\beta` to the outgoing electron line, and a spacetime index :math:`\mu` to the photon line. Moreover, such a vertex contributes the following factor

    .. math:: (2\pi)^4 e (\gamma^{\mu})_{\beta \alpha} \delta^4(k - k' \pm q)

    where :math:`k, k'` are the :math:`4`-momentum of the incoming and outgoing electrons, respectively, and :math:`\pm q` is the :math:`4`-momentum of the incoming/outgoing photon.

External lines
    In what follows, we'll use :math:`\alpha` (for incoming) and :math:`\beta` (for outgoing) to index Dirac spinors as in :eq:`eq_dirac_field_psi_field`, and :math:`\mu` to index photon spinors as in :eq:`eq_qed_a_field_general_solution`.

    * For each in-state electron line, include a factor :math:`(2\pi)^{-3/2} u_{\alpha}(\pbf, \sigma)`.
    * For each in-state positron line, include a factor :math:`(2\pi)^{-3/2} \bar{v}_{\beta}(\pbf, \sigma)`.
    * For each out-state electron line, include a factor :math:`(2\pi)^{-3/2} \bar{u}_{\beta}(\pbf, \sigma)`.
    * For each out-state positron line, include a factor :math:`(2\pi)^{-3/2} v_{\alpha}(\pbf, \sigma)`.
    * For each in-state photon line, include a factor :math:`(2\pi)^{-3/2} (2p_0)^{-1/2} e_{\mu}(\pbf, \sigma)`.
    * For each out-state photon line, include a factor :math:`(2\pi)^{-3/2} (2p_0)^{-1/2} e_{\mu}^{\ast}(\pbf, \sigma)`.

Internal lines
    In what follows, the same indexing convention as for external lines is used. In fact, to index multiple photon spinors, we'll use :math:`\mu` and :math:`\nu`.

    * For each internal electron line from :math:`\beta` to :math:`\alpha`, include a factor (cf. :eq:`eq_p_polynomial_dirac` and :eq:`eq_propagator_as_momentum_space_integral`)

      .. math::

            \frac{-\ifrak}{(2\pi)^4} \frac{\left( -\ifrak \kslash + m \right)_{\alpha \beta}}{k^2 + m^2 - \ifrak \epsilon}

      where the :math:`\beta` matrix in :eq:`eq_p_polynomial_dirac` is dropped because of the use of :math:`\bar{\psi}`, rather than :math:`\psi^{\dagger}`, in :eq:`eq_qed_concrete_interaction_massive_spin_half_fermion`. Moreover, the `Feynman slash notation <https://en.wikipedia.org/wiki/Feynman_slash_notation>`__

      .. math:: \kslash \coloneqq \gamma^{\mu} k_{\mu}

      is adopted here.
    * For each internal photon line between :math:`\mu` and :math:`\nu`, include a factor (cf. :eq:`eq_qed_photon_line_contribution_covariant`)

      .. math::

            \frac{-\ifrak}{(2\pi)^4} \frac{\eta_{\mu \nu}}{q^2 - \ifrak \epsilon}

The rest of the calculation is the same as any general evaluation of scattering amplitudes using Feynman diagrams, and was discussed in :ref:`sec_spacetime_feynman_rules`.

Although the full S-matrix calculation can be done in principle by summing up all Feynman diagrams, it quickly becomes computationally infeasible as the number of vertices in the diagram increases. Therefore the calculation is only practically possible if the contribution of a Feynman diagram also diminishes as the complexity grows. Fortunately, this is indeed the case as we now demonstrate.

Consider a connected Feynman diagram [#connected_feynman_diagram]_ with :math:`V` vertices, :math:`I` internal lines, :math:`E` external lines, and :math:`L` (independent) loops. Then the following relations hold

.. math::
    :label: eq_qed_feynman_diagram_graph_identities

    L &= I - V + 1 \\
    2I + E &= 3V

Now the goal is to calculate the constant coefficient involving :math:`e` (electric charge) and :math:`\pi` (mathematical constant) for any given diagram. According to the explicit Feynman rules described above

* Each vertex contributes a factor :math:`(2\pi)^4 e`.
* Each internal line contributes a factor :math:`(2\pi)^{-4}`.
* Each loop contributes a momentum-space integral (cf. :eq:`eq_boson_boson_momentum_space_scattering`), which, in turn, contributes a factor :math:`\pi^2` (cf. `Volume of an n-ball <https://en.wikipedia.org/wiki/Volume_of_an_n-ball>`__).

Note that we don't include contributions from external lines since they are fixed by the in- and out-states. Multiplying the contributions above all together (and eliminating :math:`I` and :math:`V` using :eq:`eq_qed_feynman_diagram_graph_identities`), the constant coefficient of a Feynman diagram is given as follows

.. math::

    (2\pi)^{4V} e^V (2\pi)^{-4I} \pi^{2L} &= (2\pi)^{-4(L-1)} e^{2L+E-2} \pi^{2L} \\
        &= (2\pi)^4 e^{E-2} \left( \frac{e^2}{16\pi^2} \right)^L

It follows that the Feynman diagrams may be organized by an increasing number of (independent) loops so that the higher order terms are suppressed by a power of

.. math:: \frac{e^2}{16\pi^2} \approx 5.81 \times 10^{-4}


Compton Scattering
------------------


.. rubric:: Footnotes

.. [#brst_quantization] The more modern way to quantize a field theory with gauge symmetry is via the so-called `BRST quantization <https://en.wikipedia.org/wiki/BRST_quantization>`__, which is Lorentz invariant. We will come back to it much later when we discuss non-Abelian gauge symmetries.

.. [#lorenz] It's really unfortunate for L. Lorenz to work in the same field as H. Lorentz and be completely overshadowed. Apparently Weinberg thought this gauge condition was named after the more famous Nobel laureate.

.. [#solve_poisson_equation] My favorite solution to Poisson's equation is given by Feynman in his `lecture on electric field <https://www.feynmanlectures.caltech.edu/II_06.html>`__.

.. [#letter_e] If using the same letter :math:`e` for the spinor and the mathematical constant was not confusing enough, it should be now by using it also for the electric charge.

.. [#connected_feynman_diagram] Recall from :ref:`sec_cluster_decomposable_hamiltonians` that only connected diagrams contribute to the S-matrix due to the cluster decomposition principle.
