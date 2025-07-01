Non-Perturbative Methods
========================

TBD

Symmetry
--------

In this section, we'll revisit various symmetry principles encountered in perturbation theory from a non-perturbative perspective.

First let's consider the translation symmetry, which generates a :math:`4`-momentum operator :math:`P_{\mu}`. According to :eq:`eq_momenta_act_as_spatial_derivative` we have

.. math:: \left[P_{\mu}, O(x)\right] = \ifrak \frac{\p}{\p x^{\mu}} O(x)

for any operator :math:`O(x)` that doesn't explicitly involve spacetime coordinates. If the in- and out-states are also taken to have definite momenta in the following sense

.. math:: P^{\mu} \Psi_{\alpha}^+ = p_{\alpha}^{\mu} \Psi_{\alpha}^+, \qquad P^{\mu} \Psi_{\beta}^- = p_{\beta}^{\mu} \Psi_{\beta}^-

then the time-ordered product amplitude considered in Gell-Mann and Low's theorem :eq:`eq_gell_mann_low_theorem` satisfies the following differential equation

.. math::

    &(p_{\beta\mu} - p_{\alpha\mu})\left(\Psi_{\beta}^-, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^+\right) \\
        &\qquad= \left(\Psi_{\beta}^-, \left[P_{\mu}, T\{O_a(x_1), O_b(x_2), \cdots\}\right] \Psi_{\alpha}^+\right) \\
        &\qquad= \ifrak \left(\frac{\p}{\p x_1^{\mu}} + \frac{\p}{\p x_2^{\mu}} + \cdots\right) \left(\Psi_{\beta}^-, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^+\right)

The general solution to the equation is given as follows

.. math::

    \left(\Psi_{\beta}^-, T\{O_a(x_1), O_b(x_2), \cdots\} \Psi_{\alpha}^+\right)
        = e^{\ifrak(p_{\alpha}-p_{\beta}) \cdot x} F_{ab\cdots}(x_1-x_2, x_1-x_3, \cdots)

where :math:`x = c_1x_1 + c_2x_2 + \cdots` as long as :math:`c_1 + c_2 + \cdots = 1`. Moreover :math:`F_{ab\cdots}` can be any function that depends only on the coordinate differences.
