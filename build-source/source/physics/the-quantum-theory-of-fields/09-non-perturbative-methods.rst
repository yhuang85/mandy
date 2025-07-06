Non-Perturbative Methods
========================

TBD

Symmetry
--------

In this section, we'll revisit various symmetry principles encountered in perturbation theory from a non-perturbative perspective.

Translational symmetry
^^^^^^^^^^^^^^^^^^^^^^

The translational symmetry generates the :math:`4`-momentum operator :math:`P_{\mu}`. According to :eq:`eq_momenta_act_as_spatial_derivative` we have

.. math:: \left[P_{\mu}, O(x)\right] = \ifrak \frac{\p}{\p x^{\mu}} O(x)

for any operator :math:`O(x)` that doesn't explicitly involve spacetime coordinates. If the in- and out-states are also taken to have definite momenta in the following sense

.. math:: P^{\mu} \Psi_{\alpha}^- = p_{\alpha}^{\mu} \Psi_{\alpha}^-, \qquad P^{\mu} \Psi_{\beta}^+ = p_{\beta}^{\mu} \Psi_{\beta}^+

then the time-ordered product amplitude considered in Gell-Mann and Low's theorem :eq:`eq_gell_mann_low_theorem` satisfies the following differential equation

.. math::

    &(p_{\beta\mu} - p_{\alpha\mu})\left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^-\right) \\
        &\qquad= \left(\Psi_{\beta}^+, \left[P_{\mu}, T\{O_a(x_1), O_b(x_2), \cdots\}\right] \Psi_{\alpha}^-\right) \\
        &\qquad= \ifrak \left(\frac{\p}{\p x_1^{\mu}} + \frac{\p}{\p x_2^{\mu}} + \cdots\right) \left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^-\right)

The general solution to the equation is given as follows

.. math::
    :label: eq_npm_timed_order_amplitude

    \left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^-\right)
        = e^{\ifrak(p_{\alpha}-p_{\beta}) \cdot x} F_{ab\cdots}(x_1-x_2, x_1-x_3, \cdots)

where :math:`x = c_1x_1 + c_2x_2 + \cdots` as long as :math:`c_1 + c_2 + \cdots = 1`. Moreover :math:`F_{ab\cdots}` can be any function that depends only on the coordinate differences. It turns out that this is enough to deduce the overall conservation of momentum in the following sense

.. math::
    :label: eq_npm_4_momentum_conservation

    &\int d^4x_1~d^4x_2\cdots \left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^-\right) \exp(-\ifrak k_1 \cdot x_1 - \ifrak k_2 \cdot x_2 - \cdots)  \\
        &\qquad\propto \delta^4(p_{\alpha} - p_{\beta} - k_1 - k_2 - \cdots)

To see this, let's consider a special case where :math:`x=x_1` in :eq:`eq_npm_timed_order_amplitude`. In this case, the power in the exponential in :eq:`eq_npm_4_momentum_conservation` can be rewritten as follows

.. math::

    -\ifrak k_1 \cdot x_1 - \ifrak k_2 \cdot x_2 - \cdots = -\ifrak\left(k_1 + k_2 + \cdots\right) \cdot x_1 - \ifrak k_2 \cdot (x_2 - x_1) - \cdots

Then the delta function in :eq:`eq_npm_4_momentum_conservation` comes out as the integral over :math:`x_1`.

The conservation of momentum of obvious from perturbation theory since the momentum is conserved at every vertex of every contributing Feynman diagram.

Internal symmetry
^^^^^^^^^^^^^^^^^

Recall from :eq:`eq_charge_of_annihilation_and_creation_operator` that if :math:`O_a(x)` is a field or some operator that either annihilates a charge :math:`q_a` or creates a charge :math:`-q_a`, then

.. math:: \left[Q, O_a(x)\right] = -q_a O_a(x)

Now if the in- and out-states :math:`\Psi_{\alpha}^-` and :math:`\Psi_{\beta}^+` have definite charges :math:`q_{\alpha}` and :math:`q_{\beta}`, respectively, then the charge is conserved in a process that involves fields/operators like :math:`O_a(x)` in the following sense

.. math::

    &(q_{\beta}-q_{\alpha}) \left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\}\Psi_{\alpha}^-\right) \\
        &\qquad = \left(\Psi_{\beta}^+, \left[Q, T\{O_a(x_1), O_b(x_2), \cdots\}\right] \Psi_{\alpha}^-\right) \\
        &\qquad = -(q_a + q_b + \cdots) \left(\Psi_{\beta}^+, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^-\right)

which implies

.. math:: q_{\beta} = q_{\alpha} - q_a - q_b - \cdots

Charge conjugation
^^^^^^^^^^^^^^^^^^

We'll consider the spin-:math:`1/2` particle, whose charge conjugation formula is given by :eq:`eq_dirac_field_charge_conjugation`. Recall from :eq:`eq_dirac_field_bilinear_form_transform_under_spatial_inversion_and_charge_conjugation` that the vector :math:`\bar{\psi} \gamma^{\mu} \psi` transforms as follows

.. math::

    U(\Ccal) \left(\bar{\psi} \gamma^{\mu} \psi\right) U^{-1}(\Ccal) = -\bar{\psi} \gamma^{\mu} \psi

Given the Lagrangian :eq:`eq_qed_lagrangian_density` for QED, the charge conjugation is to be conserved if, in addition,

.. math:: U(\Ccal) a^{\mu} U^{-1}(\Ccal) = -a^{\mu}

for free photon field :math:`a^{\mu}`. The same transformation laws hold also in the Heisenberg picture. It follows, in particular, that the vacuum expectation value of the time-ordered product of an odd number of electromagnetic fields/currents vanishes. That is, the sum of Feynman diagrams with odd number of external lines, all of which are photonic, whether on or off mass-shell, vanishes.

This result is known as `Furry's theorem <https://en.wikipedia.org/wiki/Furry%27s_theorem>`__. To see why it holds in perturbation theory, note that each diagram with an odd number of external photon lines (and no other external lines), can be countered by another diagram whose electron lines are all reverted, by applying the charge conjugation.


Polology
--------

This is a made-up a word to describe the study of the structure and distribution of poles in amplitudes. To keep the discussion on a general ground, consider the following momentum-space vacuum expectation value

.. math::

    G(q_1, q_2, \cdots, q_n) \coloneqq \int d^4x_1 \cdots d^4x_n~e^{-\ifrak q_1 \cdot x_1} \cdots e^{-\ifrak q_n \cdot x_n} \braket{T\{A_1(x_1) \cdots A_n(x_n)\}}_{\VAC}

Here the :math:`A`\s are Heisenberg-picture operators of arbitrary Lorentz type.

So far the only constraint on the off-mass-shell :math:`4`-momenta :math:`q_1, q_2, \cdots, q_n` is that they sum up to zero according to total momentum conservation. We're going to impose a further constraint that :math:`G` is a function of :math:`q^2` where

.. math:: q \coloneqq q_1 + \cdots + q_r = -q_{r+1} - \cdots - q_n

In other words, the first :math:`r` momenta combined together is constrained to the mass shell. Then we'll argue that :math:`G` has a pole at :math:`q^2 = -m^2`, where :math:`m` is the mass of any one-particle state that has non-vanishing matrix elements with the states :math:`A_1^{\dagger} \cdots A_r^{\dagger} \Psi_{\VAC}` and :math:`A_{r+1} \cdots A_n \Psi_{\VAC}`. Moreover the residue at the pole is given by

.. math::
    :label: eq_vacuum_expectation_pole_formula

    G \to \frac{-2\ifrak\sqrt{\qbf^2 + m^2}}{q^2 + m^2 - \ifrak\epsilon} (2\pi)^7 \delta^4(q_1 + \cdots + q_n) \sum_{\sigma} M_{0 | \qbf, \sigma}(q_2, \cdots, q_r) M_{\qbf, \sigma | 0}(q_{r+2}, \cdots, q_n)

where the :math:`M`\s are defined by

.. math::

    &\int d^4x_1 \cdots d^4x_r~e^{-\ifrak q_1 \cdot x_1} \cdots e^{-\ifrak q_r \cdot x_r} \left(\Psi_{\VAC}, T\{A_1(x_1), \cdots, A_r(x_r)\} \Psi_{\pbf, \sigma}\right) \\
        &\qquad = (2\pi)^4 \delta^4(q_1 + \cdots + q_r - p) M_{0 | \pbf, \sigma}(q_2, \cdots, q_r) \\
    &\int d^4x_{r+1} \cdots d^4x_n~e^{-\ifrak q_{r+1} \cdot x_{r+1}} \cdots e^{-\ifrak q_n \cdot x_n} \left(\Psi_{\pbf, \sigma}, T\{A_{r+1}(x_{r+1}), \cdots, A_n(x_n)\} \Psi_{\VAC}\right) \\
        &\qquad = (2\pi)^4 \delta^4(q_{r+1} + \cdots + q_n + p) M_{\pbf, \sigma | 0}(q_{r+2}, \cdots, q_n)

Before diving into the proof, it's instructive to rewrite :eq:`eq_vacuum_expectation_pole_formula` in a way closer to a Feynman diagram evaluation as follows

.. math::

    &G(q_1, \cdots, q_n) \to \sum_{\sigma} \int d^4k \\
        &\qquad \times \left((2\pi)^4 \delta^4(q_1 + \cdots + q_r - k) (2\pi)^{3/2} \left(2\sqrt{\kbf^2 + m^2}\right)^{1/2} M_{0 | \kbf, \sigma}(q_2, \cdots, q_r)\right) \\
        &\qquad \times \frac{-\ifrak}{(2\pi)^4} \frac{1}{k^2 + m^2 - \ifrak\epsilon} \\
        &\qquad \times \left((2\pi)^4 \delta^4(q_{r+1} + \cdots + q_n + k) (2\pi)^{3/2} \left(2\sqrt{\kbf^2 + m^2}\right)^{1/2} M_{0 | \kbf, \sigma}(q_{r+1}, \cdots, q_n)\right)
