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

.. math:: J^{\mu}(x) = -\ifrak \sum_{\ell} \frac{\delta \Lscr_M}{\delta (\p_{\mu} \Psi_{\ell}(x))} q_{\ell} \Psi_{\ell}(x)
    :label: eq_qed_conserved_current_j

Let

.. math:: Q \coloneqq \int d^3 x~J^0

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

    We can also verify :eq:`eq_qed_matter_lagrangian_variational_derivative_in_a_field` and :eq:`eq_qed_conserved_current_j` in this setup as follows

    .. math::

        \frac{\delta I_M}{\delta A_{\mu}(x)} &= \int d^4 y~\frac{\delta \Lscr_M(\Psi_{\ell}(y), D_{\mu} \Psi_{\ell}(y))}{\delta A_{\mu}(x)} \\
            &= \sum_{\ell} \frac{\delta \Lscr_M}{\delta D_{\mu} \Psi_{\ell}} \left( -\ifrak q_{\ell} \Psi_{\ell}(x) \right) \\
            &= -\ifrak \sum_{\ell} \frac{\delta \Lscr_M}{\delta (\p_{\mu} \Psi_{\ell})} q_{\ell} \Psi_{\ell}(x)

    Comparing with the Lagrangian density of massive spin-:math:`1` vector field :eq:`eq_spin_1_vector_field_lagrangian_density`, we see that the following term

    .. math::  -\frac{1}{2} m^2 A_{\mu} A^{\mu}

    is missing, or :math:`m=0`. Indeed, the appearance of such quadratic term would obviously violate the postulated gauge symmetry. This is really fantastic since we've deduced from just the existence of gauge symmetry that the particle represented by :math:`A_{\mu}`, which turns out to be photon, is massless. This is in fact one the two greatest predictions of QED that photons are massless -- the other one being that they travel at the speed of light, and has been tested against experiment up to great precision.
