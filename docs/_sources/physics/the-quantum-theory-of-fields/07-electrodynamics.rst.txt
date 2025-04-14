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
