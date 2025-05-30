Path-Integral Methods
=====================

In this chapter we'll learn about a new way to build quantum field theory that complements the canonical quantization. This is done via a path-integral formalism invented by Feynman. The advantage of the path-integral formalism over the canonical formalism, as we'll see, is the manifest Lorentz invariance in deriving the propagators, as opposed to the miraculous cancellation of local interactions, such as the Coulomb interaction :eq:`eq_qed_interaction_picture_coulomb_interaction` discussed in the previous chapter. The disadvantage, on the other hand, is that the verification of the Lorentz invariance of the S-matrix becomes obscure.

.. note::

    The gear is shifted a bit in this chapter in the following two aspects. First, we'll adopt Dirac's bra-ket notation, and second, it's the Schrödinger picture rather than the interaction picture that will be considered alongside the Heisenberg picture.

Path Integrals for Bosons
-------------------------

Similar to the canonical formalism :eq:`eq_canonical_commutation_relations`, let's consider (Hermitian) operators :math:`Q_a` and their conjugates :math:`P_a` satisfying the following canonical commutation relations

.. math::
    :label: eq_pib_canonical_commutation_relations

    [Q_a, P_b] &= \ifrak \delta_{ab} \\
    [Q_a, Q_b] &= [P_a, P_b] = 0

We'll think of :math:`Q_a` as (spacetime) coordinates and :math:`P_a` as momenta. Moreover, they are considered to be time-independent Schrödinger-picture operators. Note that the commutators are used here since the particles under consideration are bosons. The corresponding theory for fermions will be taken up in the next section.

Since all :math:`Q_a`\s commute each other, they can be simultaneously diagonalized so that there exist eigenstates :math:`\ket{q}` satisfying

.. math:: Q_a \ket{q} = q_a \ket{q}

Moreover, the eigenstates satisfy the following orthogonality condition

.. math:: \braket{q' | q} = \prod_{a} \delta(q'_a - q_a) \eqqcolon \delta(q' - q)
    :label: eq_pif_q_basis_completeness

and the completeness condition

.. math:: 1 = \int \prod_a dq_a \ketbra{q}{q}

Similar eigenstates exist for :math:`P_a`\s as well.

It's important to note that :math:`\ket{q}` and :math:`\ket{p}` live in the same Hilbert space. Indeed, let's compute :math:`\braket{q | p}` as follows. First, it follows from :eq:`eq_pib_canonical_commutation_relations` that :math:`P_a` acts on wave functions in :math:`q`-basis as

.. math:: P_b = -\ifrak \frac{\p}{\p q_b}
    :label: eq_pim_p_acts_as_dq

Indeed we have

.. math::

    \left[ Q_a, P_b \right] \left( f(q) \ket{q} \right)
        &= Q_a \left( -\ifrak \frac{\p f}{\p q_b} \ket{q} \right) - P_b \left( f(q) q_a \ket{q} \right) \\
        &= -\ifrak \frac{\p f}{\p q_b} q_a \ket{q} + \ifrak \frac{\p f}{\p q_b} q_a \ket{q} + \ifrak \delta_{ab} f(q) \ket{q} \\
        &= \ifrak \delta_{ab} f(q) \ket{q}

in agreement with :eq:`eq_pib_canonical_commutation_relations`. Next, using the completeness condition, one can write

.. math::

    \ket{p} = \int \prod_a dq_a \braket{q | p} \ket{q}

It follows then from :eq:`eq_pim_p_acts_as_dq` that

.. math::

    \int \prod_a dq_a~p_b \braket{q | p} \ket{q} = p_b \ket{p} = P_b \ket{p} = -\ifrak \int \prod_a dq_a \frac{\p \braket{q | p}}{\p q_b} \ket{q}

for any :math:`b`. Using the fact that the :math:`\ket{q}`\s form a basis, we conclude that [#qp_product_2pi_factor]_

.. math:: \braket{q | p} = \prod_a e^{\ifrak q_a p_a}
    :label: eq_pib_schrodinger_picture_qp_formula


The general path integral formula
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To derive the general path integral formula, we need to pass to the Heisenberg picture as follows

.. math::
    :label: eq_pif_defn_heisenberg_q_and_p

    Q_a(t) &= e^{\ifrak Ht} Q_a e^{-\ifrak Ht} \\
    P_a(t) &= e^{\ifrak Ht} P_a e^{-\ifrak Ht}

where the Hamiltonian :math:`H` is given as a function of :math:`P` and :math:`Q`. Their eigenstates

.. math::
    :label: eq_pif_defn_heisenberg_q_and_p_eigenstates

    Q_a(t) \ket{q, t} &= q_a \ket{q, t} \\
    P_a(t) \ket{p, t} &= p_a \ket{p, t}

are obviously given by

.. math::
    :label: eq_pif_time_dependent_p_and_q

    \ket{q, t} &= e^{\ifrak Ht} \ket{q} \\
    \ket{p, t} &= e^{\ifrak Ht} \ket{p}

.. warning::

    The eigenstates :math:`\ket{q, t}` and :math:`\ket{p, t}` given by :eq:`eq_pif_time_dependent_p_and_q` are *not* time-:math:`t` evolutions of :math:`\ket{q}` and :math:`\ket{p}` which, according to Schrödinger's equation, would be :math:`e^{-\ifrak Ht} \ket{q}` and :math:`e^{-\ifrak Ht} \ket{p}`, respectively.

The time-independent eigenstates satisfy similar orthogonality and completeness conditions as follows

.. math::

    \braket{q', t | q, t} &= \delta(q' - q) \\
    \braket{p', t | p, t} &= \delta(p' - p) \\
    1 &= \int \prod_a dq_a \ketbra{q, t}{q, t} \\
    1 &= \int \prod_a dp_a \ketbra{p, t}{p, t}

Moreover :eq:`eq_pib_schrodinger_picture_qp_formula` also carries over

.. math:: \braket{q, t | p, t} = \prod_a e^{\ifrak q_a p_a}
    :label: eq_pif_heisenberg_picture_qp_formula

Now the key idea in deriving the path integral formula is to evaluate how the eigenstates evolve in infinitesimal time steps :math:`\tau \to \tau + d\tau` as follows

.. math:: \braket{q', \tau + d\tau | q, \tau} = \braket{q', \tau | e^{-\ifrak H d\tau} | q, \tau}
    :label: eq_pif_infinitesimal_q_progression

In light of :eq:`eq_pif_defn_heisenberg_q_and_p_eigenstates`, it'll be convenient to rewrite :math:`H = H(Q, P)` in terms of :math:`Q(t)` and :math:`P(t)` defined by :eq:`eq_pif_defn_heisenberg_q_and_p`. This is done by the following calculation

.. math:: H = H(Q, P) = e^{\ifrak Ht} H(Q, P) e^{-\ifrak Ht} = H(Q(t), P(t))
    :label: eq_pif_hamiltonian_schrodinger_equals_heisenberg

Using the canonical commutation relations :eq:`eq_pib_canonical_commutation_relations`, we can make the following assumption without losing any generality.

.. admonition:: Assumption

    All the :math:`Q` operators in :math:`H` lie to the left of the :math:`P` operators.

Under this assumption, one can expand :eq:`eq_pif_infinitesimal_q_progression` for infinitesimal :math:`d\tau` using :eq:`eq_pif_defn_heisenberg_q_and_p_eigenstates` and :eq:`eq_pif_heisenberg_picture_qp_formula` as follows

.. math::
    :label: eq_pif_infinitesimal_q_progression_expanded

    \braket{q', \tau + d\tau | q, \tau} &= \braket{q', \tau | \exp\left( -\ifrak H(Q(\tau), P(\tau)) d\tau \right) | q, \tau} \\
        &= \int \prod_a dp_a \braket{q', \tau | \exp(-\ifrak H(Q(\tau), P(\tau)) d\tau) | p, \tau} \braket{p, \tau | q, \tau} \\
        &= \int \prod_a dp_a \exp(-\ifrak H(q', p) d\tau) \braket{q', \tau | p, \tau} \braket{p, \tau | q, \tau} \\
        &= \int \prod_a dp_a \exp\left( -\ifrak H(q', p) d\tau + \ifrak \sum_a (q'_a - q_a) p_a \right)

Note that the third equality holds only for infinitesimal :math:`d\tau`, which allows us to pretend that :math:`e^{-\ifrak H d\tau}` is linear in :math:`H`.

.. important::

    The function :math:`H(q', p)` in the last expression, or written simply as :math:`H(q, p)`, is an ordinary function of scalars. In particular, it makes no difference however :math:`q` and :math:`p` are ordered. It should therefore be remembered that when this process is reversed, i.e., the quantization of a classical Hamiltonian, the quantized Hamiltonian must have all the :math:`Q` operators lying to the left of the :math:`P` operators.

Now given two time :math:`t < t'` with a finite separation, one can divide the time-interval into :math:`N` steps

.. math:: t < \tau_1 < \tau_2 < \cdots < \tau_N < t'
    :label: eq_pif_time_intervals

where

.. math:: \tau_i = \frac{t' - t}{N + 1}

As :math:`N \to \infty`, one can apply :eq:`eq_pif_infinitesimal_q_progression_expanded` to each sub-interval to calculate the transition amplitude

.. math::

    \braket{q', t' | q, t} &= \int \prod_{k=1}^N dq_k \braket{q', t' | q_N, t_N} \braket{q_{N-1}, t_{N-1} | q_{N-2}, t_{N-2}} \cdots \braket{q_1, t_1 | q, t} \\
        &= \int \left( \prod_{k=1}^N \prod_a dq_{k, a} \right) \left( \prod_{k=0}^N \prod_a dp_{k, a} \right) \\
        &\quad \times \exp\left(
            \ifrak \sum_{k=0}^N \left( -H(q_{k+1}, p_k) d\tau + \sum_a (q_{k+1, a} - q_{k, a}) p_{k, a} \right)
        \right) \\
        &= \int_{\substack{q_a(t) = q_a \\ q_a(t') = q'_a }} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} dp_a(\tau) \exp\left(
            \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right)
        \right)

with the understanding that :math:`q_0 = q` and :math:`q_{N+1} = q'`. It's in the last equality where the limit :math:`N \to \infty`, or equivalently :math:`d\tau \to 0`, is taken. The integral is taken over all paths from state :math:`\ket{q}` at time :math:`t` to state :math:`\ket{q'}` at time :math:`t'`, and hence the name -- path integral.

It turns out that the same recipe for deriving the general path integral formula above can also be applied to calculate matrix elements of an operator :math:`\Oscr(P(t), Q(t))`, or more generally a time-ordered product of such operators. Note that in contrast to the Hamiltonian (cf. :eq:`eq_pif_hamiltonian_schrodinger_equals_heisenberg`), we've swapped the order of arguments :math:`Q, P` in :math:`\Oscr`. This is, for reasons which will become clear momentarily, due to the following arrangement.

.. admonition:: Assumption

    All the :math:`P` operators in :math:`\Oscr` lie to the left of the :math:`Q` operators.

As before, let's first calculate the infinitesimal matrix element as follows

.. math::

    \braket{q', \tau + d\tau | \Oscr(P(\tau), Q(\tau)) | q, \tau}
        &= \int \prod_a dp_a \braket{q', \tau | \exp(-\ifrak H d\tau) | p, \tau} \braket{p, \tau | \Oscr | q, \tau} \\
        &= \int \prod_a dp_a \exp\left( -\ifrak H(q', p) d\tau \right) \Oscr(p, q) \braket{q', \tau | p, \tau} \braket{p, \tau | q, \tau} \\
        &= \int \prod_a dp_a \exp\left( -\ifrak H(q', p) d\tau + \ifrak \sum_a (q'_a - q_a) p_a \right) \Oscr(p, q)

Consider a time-ordered sequence of operators

.. math:: \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots

such that :math:`t_A > t_B > \cdots`. We can calculate the matrix element of the product of the operators at a finite time difference by dividing the time-interval in the same way as in :eq:`eq_pif_time_intervals` and pay attention to the sub-intervals that contains :math:`t_A, t_B, \cdots`, as follows

.. math::

    &\braket{q', t' | \Oscr_A(P(t_A), Q(t_A)) \Oscr_B(P(t_B), Q(t_B)) \cdots | q, t} \\
    &\quad = \int_{\substack{q_a(t)=q_a \\ q_a(t')=q'_a}} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} dp_a(\tau)
        \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right) \right)

Since the right-hand-side doesn't rely on the time-ordering, we may replace the product of operators on the left-hand-side with the timed-ordered product as follows

.. math::
    :label: eq_pif_time_ordered_product_matrix_element

    &\braket{q', t' | T\left\{ \Oscr_A(P(t_A), Q(t_A)) \Oscr_B(P(t_B), Q(t_B)) \cdots \right\} | q, t} \\
    &\quad = \int_{\substack{q_a(t)=q_a \\ q_a(t')=q'_a}} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} dp_a(\tau)
        \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right) \right)

as long as :math:`t_A, t_B, \cdots` are all distinct.


Transition to the S-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^

From now on, we will restrict the discussion to quantum field theories where the index :math:`a` from the previous section becomes :math:`(\xbf, m)`, where :math:`\xbf` is the spatial coordinates and :math:`m` denotes other quantum labels such as spin. In this case we rewrite :eq:`eq_pif_time_ordered_product_matrix_element` as follows

.. math::

    &\braket{q', t' | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | q, t} \\
    &\quad = \int_{\substack{q_m(t, \xbf)=q_m(\xbf) \\ q_m(t', \xbf')=q_m(\xbf')}} \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf) \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3 x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) \right) \right)

Recall that the S-matrix involves matrix elements between in- and out-states, which are states are time :math:`t = \mp\infty`, respectively. Hence if we write :math:`\ket{\alpha, \op{in}}` for the in-state and :math:`\ket{\beta, \op{out}}` for the out-state, then the S-matrix element can be written as follows

.. math::
    :label: eq_pi_to_s_matrix_timed_ordered_matrix_element

    &\braket{\beta, \op{out} | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | \alpha, \op{in}} \\
    &\quad = \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf) \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3 x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) \right) \right) \\
    &\qquad \times \braket{\beta, \op{out} | q(\infty), \infty} \braket{q(-\infty), -\infty | \alpha, \op{in}}

where the path integral now has essentially no boundary conditions.

The goal now is to calculate the wave functions :math:`\braket{\beta, \op{out} | q(\infty), \infty}` and :math:`\braket{q(-\infty), -\infty | \alpha, \op{in}}`, if we choose a specific basis for the in- and out-states. It turns out, following discussions in :ref:`sec_external_edges_off_the_mass_shell`, that it suffices to consider the vacuum state :math:`\ket{\VAC}`. Moreover we'll not distinguish between :math:`\ket{\VAC, \op{in}}` and :math:`\ket{\VAC, \op{out}}` since the calculations will mostly be the same.

The vacuum state, being a state with no particles, can be characterized by

.. math:: a(\pbf, \sigma, n) \ket{\VAC} = 0
    :label: eq_pi_to_s_matrix_a_annihilates_vacuum

where :math:`a(\pbf, \sigma, n)` is the operator that annihilates a particle with momentum :math:`\pbf`, spin :math:`z`-component :math:`\sigma`, and other quantum numbers :math:`n`.

For simplicity, we'll focus on the real scalar field given by :eq:`eq_scalar_field_psi_by_creation_and_annihilation_operators` and turned into canonical variables following :eq:`eq_defn_q_and_p_scalar_field_self_dual` as follows

.. math::

    \Phi(t, \xbf) &= (2\pi)^{-3/2} (2E)^{-1/2} \int d^3 p \left( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right) \\
    \Pi(t, \xbf) &= \dot{\Phi}(t, \xbf) = -\ifrak (2\pi)^{-3/2} (E/2)^{1/2} \int d^3 p \left( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)

where :math:`E = p_0 = \sqrt{\pbf^2 + m^2}` on the on-mass-shell energy. From these one can solve for :math:`a(\pbf)` as follows

.. math::

    a(\pbf) &= (2\pi)^{-3/2} \int d^3 x~e^{-\ifrak p \cdot x} \left( (E/2)^{1/2} \Phi(t, \xbf) + \ifrak (2E)^{-1/2} \Pi(t, \xbf) \right) \\
        &= (2\pi)^{-3/2} e^{\ifrak Et} \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \left( (E/2)^{1/2} \Phi(t, \xbf) + \ifrak (2E)^{-1/2} \Pi(t, \xbf) \right)

where we've pulled out the time-dependency since in order to apply it to in- and out-state, we need to take the limits :math:`t \to \mp\infty`, respectively. More explicitly, one can write

.. math:: a_{\op{in}}(\pbf) = \lim_{t \to -\infty} a(\pbf), \quad a_{\op{out}}(\pbf) = \lim_{t \to \infty} a(\pbf)

It turns out that the time limits are not really relevant in calculating the wave functions since :math:`e^{\ifrak Et}` is never zero. Hence we'll continue to just use :math:`a(\pbf)` in calculations. In the same vein, define Schrödinger-picture operators

.. math:: \phi(\mp\infty, \xbf) = \lim_{t \to \mp\infty} \Phi(t, \xbf), \quad \pi(\mp\infty, \xbf) = \lim_{t \to \mp\infty} \Pi(t, \xbf)

In places where specifying :math:`t = \mp\infty` doesn't matter, we'll also simply write :math:`\phi(\xbf)` and :math:`\pi(\xbf)`.

Using :eq:`eq_pi_to_s_matrix_a_annihilates_vacuum`, one finds a differential equation that the wave functions :math:`\braket{\phi(\mp\infty, \xbf), \mp\infty | \VAC}` must satisfy as follows

.. math::
    :label: eq_pi_to_s_matrix_differential_equation_for_wave_function

    & \braket{\phi(\mp\infty), \mp\infty | a(\pbf) | \VAC} = 0 \\
    \implies & \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \left( \frac{\delta}{\delta \phi(\xbf)} + E(\pbf)\phi(\xbf) \right) \braket{\phi(\mp\infty, \xbf), \mp\infty | \VAC} = 0

where we have also used the interpretation of :math:`\pi(\xbf)` as variational derivative :math:`-\ifrak \delta/\delta \phi(\xbf)` (cf. :eq:`eq_pim_p_acts_as_dq`). Based on the experience of solving an analogous ODE by exponential function, it's quite natural to postulate a Gaussian solution

.. math:: \braket{\phi(\mp\infty, \xbf), \mp\infty | \VAC} = \Nscr \exp\left( -\frac{1}{2} \int d^3 x~d^3 y~\Escr(\xbf, \ybf) \phi(\xbf) \phi(\ybf) \right)
    :label: eq_pi_to_s_matrix_wave_functions

where :math:`\Nscr` is a constant. Indeed :eq:`eq_pi_to_s_matrix_differential_equation_for_wave_function` becomes equivalent to

.. math::

    0 &= \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \left( \int d^3 y~\Escr(\xbf, \ybf) \phi(\ybf) - E(\pbf) \phi(\xbf) \right) \\
        &= \int d^3 x~d^3 y~e^{\ifrak \pbf \cdot \xbf} \Escr(\xbf, \ybf) \phi(\ybf) - \int d^3 y~e^{\ifrak \pbf \cdot \ybf} E(\pbf) \phi(\ybf) \\
        &= \int d^3 y~\phi(\ybf) \left( \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \Escr(\xbf, \ybf) - e^{\ifrak \pbf \cdot \ybf} E(\pbf) \right)

For the right-hand-side to vanish for any :math:`\phi`, the quantity in the parenthesis must vanish. An inverse Fourier transform then gives

.. math:: \Escr(\xbf, \ybf) =  (2\pi)^{-3} \int d^3 p~e^{\ifrak \pbf \cdot (\xbf - \ybf)} E(\pbf)

where we recall once again that :math:`E(\pbf) = \sqrt{\pbf^2 + m^2}`. This solves :eq:`eq_pi_to_s_matrix_wave_functions` up to an unknown field-independent constant :math:`\Nscr`, which turns out to be insignificant. Indeed, the same constant :math:`\Nscr` also appears in :math:`\braket{\VAC, \op{out} | \VAC, \op{in}}` and hence can be eliminated by normalization.

We can continue the calculation :eq:`eq_pi_to_s_matrix_timed_ordered_matrix_element` in the case of vacuum expectation values for real scalar fields as follows

.. math::

    & \braket{\VAC, \op{out} | \phi(\infty), \infty} \braket{\phi(-\infty), -\infty | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \exp\left( -\frac{1}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \left( \phi(\infty, \xbf) \phi(\infty, \ybf) + \phi(-\infty, \xbf) \phi(-\infty, \ybf) \right) \right) \\
    &\quad = |\Nscr|^2 \lim_{\epsilon \to 0+} \exp\left( -\frac{\epsilon}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \int_{-\infty}^{\infty} d\tau~\phi(\tau, \xbf) \phi(\tau, \ybf) e^{-\epsilon |\tau|} \right)

and therefore

.. math::

    & \braket{\VAC, \op{out} | T\left\{ \Oscr_A(\Pi(t_A), \Phi(t_A)), \Oscr(\Pi(t_B), \Phi(t_B)), \cdots \right\} | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \int \prod_{\tau, \xbf} d\phi(\tau, \xbf) \prod_{\tau, \xbf} d\pi(\tau, \xbf)~\Oscr_A(\Pi(t_A), \Phi(t_A)) \Oscr_B(\Pi(t_B), \Phi(t_B)) \cdots \\
    &\qquad \times \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(\phi(\tau), \pi(\tau)) + \int d^3x~\dot{\phi}(\tau, \xbf) \pi(\tau, \xbf) \right.\right. \\
        &\qquad \left.\left. + \frac{\ifrak\epsilon}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \phi(\tau, \xbf) \phi(\tau, \ybf) e^{-\epsilon |\tau|} \right)
    \right)

Without working out the details, we claim that the only difference in the calculation for general fields is the term after :math:`\ifrak\epsilon/2`, whose exact form turns out to be insignificant. For later references, the final result is recorded as follows

.. math::
    :label: eq_pi_to_s_matrix_general_vacuum_matrix_element

    & \braket{\VAC, \op{out} | T\left\{ \Oscr_A(\Pi(t_A), \Phi(t_A)), \Oscr(\Pi(t_B), \Phi(t_B)), \cdots \right\} | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf} dp_m(\tau, \xbf)~\Oscr_A(P(t_A), Q(t_A)) \Oscr_B(P(t_B), Q(t_B)) \cdots \\
    &\qquad \times \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) + \ifrak\epsilon \text{ terms} \right)
    \right)


Lagrangian version of the path integral
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far the path integral formalism has been developed using the Hamiltonian. Now we'll develop a version based on the Lagrangian. In fact, the integrand in the exponential power in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element`, leaving alone the :math:`\ifrak\epsilon` terms, looks just like the corresponding Lagrangian (cf. :eq:`eq_legendre_transformation_lagrangian_from_hamiltonian`). However, there is an important difference, namely, the :math:`q` and :math:`p` variables in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` are independent variables, while in the Lagrangian formalism, they are related by :eq:`eq_hamilton_equation_in_heisenberg_picture`. As we'll see, it turns out that when the Hamiltonian :math:`H` is quadratic in :math:`p` and the (timed-ordered) operators :math:`\Oscr_A, \Oscr_B, \cdots`, are independent of the :math:`P`\s, one can explicitly evaluate the integral in :math:`p` in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element`, which will then produce the Lagrangian version of the path integral.


.. rubric:: Footnotes

.. [#qp_product_2pi_factor] It's unclear to me why a factor of :math:`1/\sqrt{2\pi}`, which clearly diminishes the (infinite) product, is inserted in [Wei95]_ page 379 eq. (9.1.12). It might simply be a mistake considering its fermionic counterpart eq. (9.5.27) in page 403.
