Quantum Electrodynamics
=======================

We saw in :ref:`sec_the_canonical_formalism` how Lagrangian (densities) may be guessed, passed to the interaction picture, and verified against the free fields derived from Lorentz invariance and causality in :ref:`sec_quantum_fields_and_antiparticles`. There were, however, cases we couldn't handle back in :ref:`sec_massless_fields`, which include, in particular, massless vector fields. It is this difficulty we will try to handle in this section, which, as it turns out, leads to the most successful quantum field theory to date, namely, the quantum electrodynamics (QED).

.. todo::
    Write down the general plan which I don't yet know.


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

where the middle equality is one of the Euler-Lagrange equations :eq:`eq_euler_lagrange` with the time-derivative term dropped out due to the primary constraint, the equality on the left is the definition :eq:`eq_qed_defn_conjugate_field_pi`, and the equality on the right is due to the coupling term :math:`J^{\mu} A_{\mu}` in :math:`\Lscr`.

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


.. rubric:: Footnotes

.. [#brst_quantization] The more modern way to quantize a field theory with gauge symmetry is via the so-called `BRST quantization <https://en.wikipedia.org/wiki/BRST_quantization>`__, which is Lorentz invariant. We will come back to it much later when we discuss non-Abelian gauge symmetries.

.. [#lorenz] It's really unfortunate for L. Lorenz to work in the same field as H. Lorentz and be completely overshadowed. Apparently Weinberg thought this gauge condition was named after the more famous Nobel laureate.

.. [#solve_poisson_equation] My favorite solution to Poisson's equation is given by Feynman in his `lecture on electric field <https://www.feynmanlectures.caltech.edu/II_06.html>`__.
