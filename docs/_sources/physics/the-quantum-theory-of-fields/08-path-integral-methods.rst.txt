Path-Integral Methods
=====================

In this chapter we'll learn about a new way to build quantum field theory that complements the canonical quantization. This is done via a path-integral formalism invented by Feynman. The advantage of the path-integral formalism over the canonical formalism, as we'll see, is the manifest Lorentz invariance in deriving the propagators, as opposed to the miraculous cancellation of local interactions, such as the Coulomb interaction :eq:`eq_qed_interaction_picture_coulomb_interaction` discussed in the previous chapter. The disadvantage, on the other hand, is that the verification of the Lorentz invariance of the S-matrix becomes obscure.

.. note::

    The gear is shifted a bit in this chapter in the following two aspects. First, we'll adopt Dirac's bra-ket notation, and second, it's the Schrödinger picture rather than the interaction picture that will be considered alongside the Heisenberg picture.

Path Integrals for Bosons
-------------------------

Similar to the canonical formalism :eq:`eq_canonical_commutation_relations`, let's consider Hermitian operators :math:`Q_a` and their conjugates :math:`P_a` satisfying the following canonical commutation relations

.. math::
    :label: eq_pib_canonical_commutation_relations

    [Q_a, P_b] &= \ifrak \delta_{ab} \\
    [Q_a, Q_b] &= [P_a, P_b] = 0

We'll think of :math:`Q_a` as (spacetime) coordinates and :math:`P_a` as momenta. Moreover, they are considered to be time-independent Schrödinger-picture operators. Note that the commutators are used here since the particles under consideration are bosons. The corresponding theory for fermions will be taken up in the next section.

Since all :math:`Q_a`\s commute each other, they can be simultaneously diagonalized so that there exist eigenstates :math:`\ket{q}` satisfying

.. math:: Q_a \ket{q} = q_a \ket{q}
    :label: eq_path_integral_bosonic_simultaneous_q_eigenstate

Moreover, the eigenstates satisfy the following orthogonality condition

.. math:: \braket{q' | q} = \prod_{a} \delta(q'_a - q_a) \eqqcolon \delta(q' - q)
    :label: eq_path_integral_bosonic_q_basis_orthogonality

and the completeness condition

.. math:: 1 = \int \prod_a dq_a \ketbra{q}{q}
    :label: eq_path_integral_bosonic_q_basis_completeness

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

.. math:: \ket{p} = \int \prod_a dq_a \braket{q | p} \ket{q}

It follows then from :eq:`eq_pim_p_acts_as_dq` and the fact that the :math:`\ket{q}`\s form a basis that

.. math::

    &\int \prod_a dq_a~p_b \braket{q | p} \ket{q} = p_b \ket{p} = P_b \ket{p} = -\ifrak \int \prod_a dq_a \frac{\p \braket{q | p}}{\p q_b} \ket{q} \\
    \implies & \frac{\p \braket{q | p}}{\p q_b} = \ifrak p_b \braket{q | p}

for any :math:`b`. It follows that

.. math:: \braket{q | p} = \prod_a \frac{1}{\sqrt{2\pi}} e^{\ifrak q_a p_a}
    :label: eq_pib_schrodinger_picture_qp_formula

where the factor :math:`(2\pi)^{-1/2}` is determined by the normalizing condition :math:`\braket{p' | p} = \delta(p' - p)`.


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

    Q_a(t) \ket{t, q} &= q_a \ket{t, q} \\
    P_a(t) \ket{t, p} &= p_a \ket{t, p}

are obviously given by

.. math::
    :label: eq_pif_time_dependent_p_and_q

    \ket{t, q} &= e^{\ifrak Ht} \ket{q} \\
    \ket{t, p} &= e^{\ifrak Ht} \ket{p}

.. warning::

    The eigenstates :math:`\ket{t, q}` and :math:`\ket{t, p}` given by :eq:`eq_pif_time_dependent_p_and_q` are *not* time-:math:`t` evolutions of :math:`\ket{q}` and :math:`\ket{p}` which, according to Schrödinger's equation, would be :math:`e^{-\ifrak Ht} \ket{q}` and :math:`e^{-\ifrak Ht} \ket{p}`, respectively.

The time-independent eigenstates satisfy similar orthogonality and completeness conditions as follows

.. math::
    :label: eq_path_integral_bosonic_time_dependent_q_and_p_orthogonality_and_completeness

    \braket{t, q' | t, q} &= \delta(q' - q) \\
    \braket{t, p' | t, p} &= \delta(p' - p) \\
    1 &= \int \prod_a dq_a \ketbra{t, q}{t, q} \\
    1 &= \int \prod_a dp_a \ketbra{t, p}{t, p}

Moreover :eq:`eq_pib_schrodinger_picture_qp_formula` also carries over

.. math:: \braket{t, q | t, p} = \prod_a \frac{1}{\sqrt{2\pi}} e^{\ifrak q_a p_a}
    :label: eq_pif_heisenberg_picture_qp_formula

Now the key idea in deriving the path integral formula is to evaluate how the eigenstates evolve in infinitesimal time steps :math:`\tau \to \tau + d\tau` as follows

.. math:: \braket{\tau + d\tau, q' | \tau, q} = \braket{\tau, q' | e^{-\ifrak H d\tau} | \tau, q}
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

    \braket{\tau + d\tau, q' | \tau, q} &= \braket{\tau, q' | \exp\left( -\ifrak H(Q(\tau), P(\tau)) d\tau \right) | \tau, q} \\
        &= \int \prod_a dp_a \braket{\tau, q' | \exp(-\ifrak H(Q(\tau), P(\tau)) d\tau) | \tau, p} \braket{\tau, p | \tau, q} \\
        &= \int \prod_a dp_a \exp(-\ifrak H(q', p) d\tau) \braket{\tau, q' | \tau, p} \braket{\tau, p | \tau, q} \\
        &= \int \prod_a \frac{dp_a}{2\pi} \exp\left( -\ifrak H(q', p) d\tau + \ifrak \sum_a (q'_a - q_a) p_a \right)

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

    &\braket{t', q' | t, q} \\
        &= \int \prod_{k=1}^N dq_k \braket{t', q' | t_N, q_N} \braket{t_{N-1}, q_{N-1} | t_{N-2}, q_{N-2}} \cdots \braket{t_1, q_1 | t, q} \\
        &= \int \left( \prod_{k=1}^N \prod_a dq_{k, a} \right) \left( \prod_{k=0}^N \prod_a \frac{dp_{k, a}}{2\pi} \right) \exp\left(
            \ifrak \sum_{k=0}^N \left( -H(q_{k+1}, p_k) d\tau + \sum_a (q_{k+1, a} - q_{k, a}) p_{k, a} \right)
        \right) \\
        &= \int_{\substack{q_a(t) = q_a \\ q_a(t') = q'_a }} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} \frac{dp_a(\tau)}{2\pi} \exp\left(
            \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right)
        \right)

with the understanding that :math:`q_0 = q` and :math:`q_{N+1} = q'`. It's in the last equality where the limit :math:`N \to \infty`, or equivalently :math:`d\tau \to 0`, is taken. The integral is taken over all paths from state :math:`\ket{q}` at time :math:`t` to state :math:`\ket{q'}` at time :math:`t'`, and hence the name -- path integral.

It turns out that the same recipe for deriving the general path integral formula above can also be applied to calculate matrix elements of an operator :math:`\Oscr(P(t), Q(t))`, or more generally a time-ordered product of such operators. Note that in contrast to the Hamiltonian (cf. :eq:`eq_pif_hamiltonian_schrodinger_equals_heisenberg`), we've swapped the order of arguments :math:`Q, P` in :math:`\Oscr`. This is, for reasons which will become clear momentarily, due to the following arrangement.

.. admonition:: Assumption

    All the :math:`P` operators in :math:`\Oscr` lie to the left of the :math:`Q` operators.

As before, let's first calculate the infinitesimal matrix element as follows

.. math::

    \braket{\tau + d\tau, q' | \Oscr(P(\tau), Q(\tau)) | \tau, q}
        &= \int \prod_a dp_a \braket{\tau, q' | \exp(-\ifrak H d\tau) | \tau, p} \braket{\tau, p | \Oscr | \tau, q} \\
        &= \int \prod_a dp_a \exp\left( -\ifrak H(q', p) d\tau \right) \Oscr(p, q) \braket{\tau, q' | \tau, p} \braket{\tau, p | \tau, q} \\
        &= \int \prod_a \frac{dp_a}{2\pi} \exp\left( -\ifrak H(q', p) d\tau + \ifrak \sum_a (q'_a - q_a) p_a \right) \Oscr(p, q)

Consider a time-ordered sequence of operators

.. math:: \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots

such that :math:`t_A > t_B > \cdots`. We can calculate the matrix element of the product of the operators at a finite time difference by dividing the time-interval in the same way as in :eq:`eq_pif_time_intervals` and pay attention to the sub-intervals that contains :math:`t_A, t_B, \cdots`, as follows

.. math::

    &\braket{t', q' | \Oscr_A(P(t_A), Q(t_A)) \Oscr_B(P(t_B), Q(t_B)) \cdots | t, q} \\
    &\quad = \int_{\substack{q_a(t)=q_a \\ q_a(t')=q'_a}} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} \frac{dp_a(\tau)}{2\pi}
        \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right) \right)

Since the right-hand-side doesn't rely on the time-ordering, we may replace the product of operators on the left-hand-side with the timed-ordered product as follows

.. math::
    :label: eq_pif_time_ordered_product_matrix_element

    &\braket{t', q' | T\left\{ \Oscr_A(P(t_A), Q(t_A)) \Oscr_B(P(t_B), Q(t_B)) \cdots \right\} | t, q} \\
    &\quad = \int_{\substack{q_a(t)=q_a \\ q_a(t')=q'_a}} \prod_{\tau, a} dq_a(\tau) \prod_{\tau, a} \frac{dp_a(\tau)}{2\pi}
        \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right) \right)

as long as :math:`t_A, t_B, \cdots` are all distinct.


Transition to the S-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^

From now on, we will restrict the discussion to quantum field theories where the index :math:`a` from the previous section becomes :math:`(\xbf, m)`, where :math:`\xbf` is the spatial coordinates and :math:`m` denotes other quantum labels such as spin. In this case we rewrite :eq:`eq_pif_time_ordered_product_matrix_element` as follows

.. math::

    &\braket{t', q' | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | t, q} \\
    &\quad = \int_{\substack{q_m(t, \xbf)=q_m(\xbf) \\ q_m(t', \xbf')=q_m(\xbf')}} \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} \frac{dp_m(\tau, \xbf)}{2\pi} \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3 x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) \right) \right)

Recall that the S-matrix involves matrix elements between in- and out-states, which are states are time :math:`t = \mp\infty`, respectively. Hence if we write :math:`\ket{\alpha, \op{in}}` for the in-state and :math:`\ket{\beta, \op{out}}` for the out-state, then the S-matrix element can be written as follows

.. math::
    :label: eq_pi_to_s_matrix_timed_ordered_matrix_element

    &\braket{\beta, \op{out} | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | \alpha, \op{in}} \\
    &\quad = \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} \frac{dp_m(\tau, \xbf)}{2\pi} \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
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

Using :eq:`eq_pi_to_s_matrix_a_annihilates_vacuum`, one finds a differential equation that the wave functions :math:`\braket{\mp\infty, \phi(\mp\infty, \xbf) | \VAC}` must satisfy as follows

.. math::
    :label: eq_pi_to_s_matrix_differential_equation_for_wave_function

    & \braket{\mp\infty, \phi(\mp\infty) | a(\pbf) | \VAC} = 0 \\
    \implies & \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \left( \frac{\delta}{\delta \phi(\xbf)} + E(\pbf)\phi(\xbf) \right) \braket{\mp\infty, \phi(\mp\infty, \xbf) | \VAC} = 0

where we have also used the interpretation of :math:`\pi(\xbf)` as variational derivative :math:`-\ifrak \delta/\delta \phi(\xbf)` (cf. :eq:`eq_pim_p_acts_as_dq`). Based on the experience of solving an analogous ODE by exponential function, it's quite natural to postulate a Gaussian solution

.. math:: \braket{\mp\infty, \phi(\mp\infty, \xbf) | \VAC} = \Nscr \exp\left( -\frac{1}{2} \int d^3 x~d^3 y~\Escr(\xbf, \ybf) \phi(\xbf) \phi(\ybf) \right)
    :label: eq_pi_to_s_matrix_wave_functions

where :math:`\Nscr` is a constant. Indeed :eq:`eq_pi_to_s_matrix_differential_equation_for_wave_function` becomes equivalent to

.. math::

    0 &= \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \left( \int d^3 y~\Escr(\xbf, \ybf) \phi(\ybf) - E(\pbf) \phi(\xbf) \right) \\
        &= \int d^3 x~d^3 y~e^{\ifrak \pbf \cdot \xbf} \Escr(\xbf, \ybf) \phi(\ybf) - \int d^3 y~e^{\ifrak \pbf \cdot \ybf} E(\pbf) \phi(\ybf) \\
        &= \int d^3 y~\phi(\ybf) \left( \int d^3 x~e^{\ifrak \pbf \cdot \xbf} \Escr(\xbf, \ybf) - e^{\ifrak \pbf \cdot \ybf} E(\pbf) \right)

For the right-hand-side to vanish for any :math:`\phi`, the quantity in the parenthesis must vanish. An inverse Fourier transform then gives

.. math:: \Escr(\xbf, \ybf) =  (2\pi)^{-3} \int d^3 p~e^{\ifrak \pbf \cdot (\xbf - \ybf)} E(\pbf)
    :label: eq_path_integral_scalar_field_curly_e

where we recall once again that :math:`E(\pbf) = \sqrt{\pbf^2 + m^2}`. This solves :eq:`eq_pi_to_s_matrix_wave_functions` up to an unknown field-independent constant :math:`\Nscr`, which turns out to be insignificant. Indeed, the same constant :math:`\Nscr` also appears in :math:`\braket{\VAC, \op{out} | \VAC, \op{in}}` and hence can be eliminated by normalization. More details about this will be discussed in the next section.

We can continue the calculation :eq:`eq_pi_to_s_matrix_timed_ordered_matrix_element` in the case of vacuum expectation values for real scalar fields as follows

.. math::
    :label: eq_path_integral_scalar_field_vacuum_wave_function

    & \braket{\VAC, \op{out} | \infty, \phi(\infty)} \braket{-\infty, \phi(-\infty) | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \exp\left( -\frac{1}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \left( \phi(\infty, \xbf) \phi(\infty, \ybf) + \phi(-\infty, \xbf) \phi(-\infty, \ybf) \right) \right) \\
    &\quad = |\Nscr|^2 \lim_{\epsilon \to 0+} \exp\left( -\frac{\epsilon}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \int_{-\infty}^{\infty} d\tau~\phi(\tau, \xbf) \phi(\tau, \ybf) e^{-\epsilon |\tau|} \right)

and therefore

.. math::
    :label: eq_path_integral_vacuum_expectation_value_scalar_field

    & \braket{\VAC, \op{out} | T\left\{ \Oscr_A(\Pi(t_A), \Phi(t_A)), \Oscr(\Pi(t_B), \Phi(t_B)), \cdots \right\} | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \int \prod_{\tau, \xbf} d\phi(\tau, \xbf) \prod_{\tau, \xbf} \frac{d\pi(\tau, \xbf)}{2\pi} \Oscr_A(\Pi(t_A), \Phi(t_A)) \Oscr_B(\Pi(t_B), \Phi(t_B)) \cdots \\
    &\qquad \times \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(\phi(\tau), \pi(\tau)) + \int d^3x~\dot{\phi}(\tau, \xbf) \pi(\tau, \xbf) \right.\right. \\
        &\qquad \left.\left. + \frac{\ifrak\epsilon}{2} \int d^3x~d^3y~\Escr(\xbf, \ybf) \phi(\tau, \xbf) \phi(\tau, \ybf) e^{-\epsilon |\tau|} \right)
    \right)

Without working out the details, we claim that the only difference in the calculation for general fields is the term after :math:`\ifrak\epsilon/2`, whose exact form turns out to be insignificant. For later references, the final result is recorded as follows

.. math::
    :label: eq_pi_to_s_matrix_general_vacuum_matrix_element

    & \braket{\VAC, \op{out} | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr(P(t_B), Q(t_B)), \cdots \right\} | \VAC, \op{in}} \\
    &\quad = |\Nscr|^2 \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf} \frac{dp_m(\tau, \xbf)}{2\pi} \Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) + \ifrak\epsilon \text{ terms} \right)
    \right)

where the :math:`\ifrak\epsilon` terms depend only on :math:`q`\s.

.. _sec_lagrangian_version_of_the_path_integral:

Lagrangian version of the path integral
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far the path integral formalism has been developed using the Hamiltonian. Now we'll develop a version based on the Lagrangian. In fact, the integrand in the exponential power in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element`, leaving alone the :math:`\ifrak\epsilon` terms, looks just like the corresponding Lagrangian (cf. :eq:`eq_legendre_transformation_lagrangian_from_hamiltonian`). However, there is an important difference, namely, the :math:`q` and :math:`p` variables in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` are independent variables, while in the Lagrangian formalism, they are related by :eq:`eq_hamilton_equation_in_heisenberg_picture`. As we'll see, it turns out that when the Hamiltonian :math:`H` is quadratic in :math:`p` and the (timed-ordered) operators :math:`\Oscr_A, \Oscr_B, \cdots`, are independent of the :math:`P`\s, one can explicitly evaluate the integral in :math:`p` in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element`, which will then produce the Lagrangian version of the path integral.

To spell out the details, let's write down the (Heisenberg-picture) Hamiltonian in the most general form as follows

.. math::
    :label: eq_hamiltonian_quadratic_in_p

    H(Q, P) &= \frac{1}{2} \sum_{n, m} \int d^3x~d^3y~A_{\xbf n, \ybf m}(Q) P_n(\xbf) P_m(\ybf) \\
        &\quad + \sum_n \int d^3x~B_{\xbf n}(Q) P_n(\xbf) + C(Q)

where :math:`A` is a real, symmetric, positive matrix. Moreover :math:`H` is written in the way that all the :math:`Q` operators lie to the left of the :math:`P` operators.

Now we can write the power in the exponential in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` without the :math:`\ifrak\epsilon` terms as follows

.. math::
    :label: eq_path_integral_exp_power_quadratic

    &\int d\tau \left( -H(q(\tau), p(\tau)) + \int d^3x \sum_n \dot{q}_n(\tau, \xbf) p_n(\tau, \xbf) \right) \\
    &\quad = -\frac{1}{2} \sum_{n, m} \int d\tau~d\tau'~d^3x~d^3y~A_{\xbf n, \ybf m}(q(\tau)) \delta(\tau - \tau') p_n(\tau, \xbf) p_m(\tau', \ybf) \\
    &\qquad - \sum_n \int d\tau~d^3x \left( B_{\xbf n}(q(\tau)) - \dot{q}_n(\tau, \xbf) \right) p_n(\tau, \xbf) - \int d\tau~C(q(\tau))

where it's organized so that the first summand on the right-hand-side is quadratic in :math:`p`, the second is linear, and the third is independent of :math:`p`. The reason to arrange the power in this form is because of the following (finite-dimensional) Gaussian integral formula.

    **Gaussian Integral Formula**

    .. math::
        :label: eq_gaussian_integral_formula

        &\int_{-\infty}^{\infty} \prod_s d\xi_s \exp\left( -\ifrak \left( \frac{1}{2} \sum_{s, r} \Ascr_{sr} \xi_s \xi_r + \sum_s \Bscr_s \xi_s + \Cscr_s \right) \right) \\
        &\quad = \left( \det(\ifrak \Ascr / 2\pi) \right)^{-1/2} \exp\left( -\ifrak \left( \sum_{s, r} \Ascr_{sr} \bar{\xi}_s \bar{\xi}_r + \sum_s \Bscr_s \bar{\xi}_s + \Cscr_s \right) \right)

    where :math:`\bar{\xi}` is the (unique) stationary point of the quadratic power given explicitly by

    .. math:: \bar{\xi}_s = -\sum_r (\Ascr^{-1})_{sr} \Bscr_r

.. note::

    In more general cases where :math:`H` is not quadratic in :math:`P`, approximation techniques such as the `stationary phase approximation <https://en.wikipedia.org/wiki/Stationary_phase_approximation>`__ may be applied.

To figure out the stationary point of the power in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` with respect to :math:`p`, let's calculate the following variational derivative assuming the :math:`\ifrak\epsilon` terms are independent of the :math:`p`\s

.. math::

    &\frac{\delta}{\delta p_n(t, \xbf)} \int_{-\infty}^{\infty} d\tau \left(
        -H(q(\tau), p(\tau)) + \int d^3y \sum_m \dot{q}_m(\tau, \ybf) p_m(\tau, \ybf) + \ifrak\epsilon \text{ terms}
    \right) \\
    &\quad = - \frac{\delta H}{\delta p_n(t, \xbf)} + \dot{q}_n(t, \xbf)

It follows that :math:`\bar{p}` is stationary if it satisfies Hamilton's equation

.. math:: \dot{q}_n(t, \xbf) = \left. \frac{\delta H}{\delta p_n(t, \xbf)} \right|_{p=\bar{p}}
    :label: eq_path_integral_stationary_p_bar

Assuming, in addition, that the (timed-ordered) operators :math:`\Oscr_A, \Oscr_B, \cdots`, are independent of the :math:`P`\s, we can evaluate the :math:`p`-integral in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` using the (infinite-dimensional) Gaussian integral formula :eq:`eq_gaussian_integral_formula` as follows

.. math::

    &\int \prod_{\tau, \xbf} \frac{dp_m(\tau, \xbf)}{2\pi} \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(q(\tau), p(\tau)) + \int d^3x \sum_m \dot{q}_m(\tau, \xbf) p_m(\tau, \xbf) + \ifrak\epsilon \text{ terms} \right)
    \right) \\
    &\quad = \left( \det(2\pi\ifrak\Ascr(q)) \right)^{-1/2} \exp\left(
        \ifrak \int_{-\infty}^{\infty} d\tau \left( L(q(\tau), \dot{q}(\tau)) + \ifrak\epsilon \text{ terms} \right)
    \right)

where :math:`L` is the Lagrangian defined by

.. math:: L(q(\tau), \dot{q}(\tau)) \coloneqq -H(q(\tau), \bar{p}(\tau)) + \int d^3x \sum_m \dot{q}_m(\tau, \xbf) \bar{p}_m(\tau, \xbf)
    :label: eq_path_integral_defn_lagrangian

with :math:`\bar{p}` satisfying :eq:`eq_path_integral_stationary_p_bar` and

.. math:: \Ascr_{\tau \xbf n, \tau' \ybf m}(q) \coloneqq A_{\xbf n, \ybf m}(q(\tau)) \delta(\tau-\tau')
    :label: eq_path_integral_a_matrix

is given by :eq:`eq_path_integral_exp_power_quadratic`.

Finally, we can write down the Lagrangian version of :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` as follows

.. math::
    :label: eq_path_integral_operator_vacuum_matrix_element_lagrangian

    &\braket{\VAC, \op{out} | T\left\{ \Oscr_A(Q(t_A)), \Oscr_B(Q(t_B)), \cdots \right\} | \VAC, \op{in}} \\
        &\quad = |\Nscr|^2 \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \Oscr_A(Q(t_A)) \Oscr_B(Q(t_B)) \cdots \\
        &\qquad \times \left( \det(2\pi\ifrak\Ascr(q)) \right)^{-1/2} \exp\left(
            \ifrak \int_{-\infty}^{\infty} d\tau \left( L(q(\tau), \dot{q}(\tau)) + \ifrak\epsilon \text{ terms} \right)
        \right)

The rest of this section is devoted to the determination of :math:`\Ascr(q)` in various examples.

Scalar fields with non-derivative coupling
    Following :eq:`eq_canonical_to_interaction_scalar_field_with_derivative_coupling_lagrangian`, consider the following Lagrangian density of a set of (massless) scalar fields :math:`\Phi_n` that have only non-derivative interaction :math:`V` and are coupled to external currents :math:`J_n`

    .. math:: \Lscr = -\sum_n \left( \frac{1}{2} \p_{\mu} \Phi_n \p^{\mu} \Phi_n + J_n^{\mu} \p_{\mu} \Phi_n \right) - V(\Phi)

    The canonical adjoint :math:`\Pi_n` is, according to :eq:`eq_general_lagrangian_conjugate_pi`, given by

    .. math:: \Pi_n = \frac{\delta \Lscr}{\delta \dot{\Phi}_n} = \dot{\Phi}_n - J_n^0

    and hence the Hamiltonian is, according to :eq:`eq_legendre_transformation_hamiltonian_from_lagrangian`, given by

    .. math::

        H &= \int d^3x \left( \sum_n \Pi_n \dot{\Phi}_n - \Lscr \right) \\
            &= \int d^3x \sum_n \left( \Pi_n (\Pi_n + J_n^0) - \frac{1}{2} (\Pi_n + J_n^0)^2 + \frac{1}{2} (\nabla \Phi_n)^2 + J_n^0 (\Pi_n + J_n^0) + \Jbf_n \cdot \nabla \Phi_n \right) \\
            &\quad + \int d^3x~V(\Phi) \\
            &= \int d^3x \sum_n \left( \frac{1}{2} (\Pi_n + J^0_n)^2 + \frac{1}{2} (\nabla \Phi_n)^2 + \Jbf_n \cdot \nabla \Phi_n \right) + \int d^3x~V(\Phi)

    Comparing with :eq:`eq_hamiltonian_quadratic_in_p` and following :eq:`eq_path_integral_a_matrix`, we see that

    .. math:: \Ascr_{x n, x' n'} = \delta^4(x-x') \delta_{nn'}

    which is field independent, and therefore can be eliminated in the same way that :math:`\Nscr` can be eliminated (cf. :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian`).

Nonlinear :math:`\sigma`-model
    The so-called nonlinear :math:`\sigma`-model is described by the following Lagrangian density

    .. math:: \Lscr = -\frac{1}{2} \sum_{n, m} \p_{\mu} \Phi_n \p^{\mu} \Phi_m (\delta_{nm} + U_{nm}(\Phi)) - V(\Phi)

    where the nonlinearity is carried by :math:`U_{nm}(\Phi)`.

    In this case the canonical adjoint :math:`\Pi_n` is given by

    .. math:: \Pi_n = \frac{\delta \Lscr}{\delta \dot{\Phi}_n} = \sum_m \dot{\Phi}_m (\delta_{nm} + U_{nm}(\Phi))

    and can be solved in matrix notation as follows

    .. math:: \dot{\Phi}_n = \sum_m (1+U(\Phi))^{-1}_{nm} \Pi_m

    hence the Hamiltonian

    .. math::

        H &= \int d^3x \left( \sum_n \Pi_n \dot{\Phi}_n - \Lscr \right) \\
            &= \int d^3x \sum_{n, m} \left(
                \frac{1}{2} \Pi_n (1+U(\Phi))^{-1}_{nm} \Pi_m
                + \frac{1}{2} \nabla \Phi_n \cdot \nabla \Phi_m (1+U(\Phi))^{-1}_{nm}
            \right) + \int d^3x~V(\Phi)

    In follows that

    .. math:: \Ascr_{xn, x'n'} = (1+U(\Phi))^{-1}_{nn'} \delta^4(x-x')
        :label: eq_path_integral_nonlinear_sigma_model_a_matrix

    which obviously depend on :math:`\Phi`, and therefore cannot be eliminated by the division by the vacuum expectation value. The idea then is to absorb it into the Lagrangian (density) which we now explain.

    Looking at :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian`, we note the following general identity

    .. math:: \det\Ascr = \exp \Tr \ln \Ascr
        :label: eq_det_eq_exp_tr_ln

    for any real symmetric positive :math:`\Ascr`. To evaluate the logarithm, it's convenient to discretize the Dirac delta function in :eq:`eq_path_integral_nonlinear_sigma_model_a_matrix` as follows

    .. math:: \delta^4(x-x') = \Omega^{-1} \delta_{xx'}

    where :math:`\Omega` denotes an infinitesimal volume in spacetime. It follows that

    .. math:: (\ln \Ascr)_{xn, x'n'} = \delta_{xx'} \left( -\ln(1+U(\Phi)) - \ln\Omega \right)_{nn'}
        :label: eq_path_integral_ln_a

    where :math:`\ln\Omega` is understood as a constant multiple as the identity matrix. Next note that the trace of :math:`\delta_{xx'}` can be evaluated by

    .. math:: \Tr~\delta_{xx'} \cdots = \Omega^{-1} \int d^4x \cdots

    It follows that

    .. math:: \det\Ascr \propto \exp\left( -\Omega^{-1} \int d^4x~\Tr\ln(1+U(\Phi)) \right)

    where the proportionality constant, coming from the constant :math:`-\ln\Omega` in :eq:`eq_path_integral_ln_a`, is field-independent. Plugging into :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian`, we see that Lagrangian density receives a correction term

    .. math:: \Delta\Lscr = -\frac{\ifrak}{2} \Omega^{-1} \Tr\ln(1+U(\Phi))

    which unfortunately contains a diverging term :math:`\Omega^{-1}`. This is known as an ultraviolet divergence since it comes from the infinitesimal spacetime volume. We'll not address how it may be handled here.

Vector fields
    The two examples considered so far admit a Lagrangian without auxiliary fields (cf. :eq:`eq_general_quantum_lagrangian`). To cover this case, consider the following Lagrangian for a set of non-interacting vector fields (cf. :eq:`eq_spin_1_vector_field_lagrangian_density`)

    .. math::
        :label: eq_path_integral_many_vector_fields_lagrangian

        \Lscr = -\sum_n \left(
            \frac{1}{4} F_{n \mu\nu} F_n^{\mu\nu} + \frac{1}{2} M^2 A_{n \mu} A_n^{\mu} + J_n^{\mu} A_{n\mu}
        \right)

    According to :eq:`eq_spin_1_vector_field_hamiltonian`, the corresponding Hamiltonian is given by

    .. math::
        :label: eq_many_vector_fields_hamiltonian

        H &= \int d^3x \sum_n \left(
            \frac{1}{2} \bm{\Pi}_n^2 + \frac{1}{2M_n^2} (\nabla \cdot \bm{\Pi}_n)^2 + \frac{1}{M_n^2} J_n^0 \nabla \cdot \bm{\Pi}_n \right. \\
            &\qquad \left. + \frac{1}{2} (\nabla \times \Abf_n)^2 + \frac{1}{2} M_n^2 \Abf_n^2 + \frac{1}{2M_n^2} (J_n^0)^2 - \Jbf_n \cdot \Abf_n
        \right)

    where the terms are ordered in descending power of :math:`\bm{\Pi}`. Using the following calculation

    .. math::

        \int d^3x~d^3y~\nabla_i \nabla_j \delta^3(\xbf-\ybf) \bm{\Pi}_n^i(x) \bm{\Pi}_n^j(y)
            &= -\int d^3x~d^3y~\nabla_j \delta^3(\xbf-\ybf) \p_i \bm{\Pi}_n^i(x) \bm{\Pi}_n^j(y) \\
            &= -\int d^3x~d^3y~\delta^3(\xbf-\ybf) \p_i \bm{\Pi}_n^i(x) \p_j \bm{\Pi}_n^j(y) \\
            &= -\int d^3x~\left( \nabla \cdot \bm{\Pi}_n \right)^2

    we conclude that

    .. math:: \Ascr_{x i n, y j m} = \delta_{nm} \left( \delta_{ij}\delta^4(x-y) - \frac{1}{M_n^2} \nabla_i\nabla_j\delta^4(x-y) \right)

    which is field-independent. As before, it means that the term :math:`\det(2\pi\ifrak \Ascr(q))^{-1/2}` in :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian` plays no role. Nonetheless, the Lagrangian defined by :eq:`eq_path_integral_defn_lagrangian` cannot be the same the original :eq:`eq_path_integral_many_vector_fields_lagrangian` since the former doesn't involve the time-component :math:`A_0`. As a consequence, the Lorentz invariance of :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian` is far from obvious.

    To restore the manifest Lorentz invariance, let's introduce, according to :eq:`eq_spin_1_vector_field_heisenberg_v0`, a correction term to the Hamiltonian :math:`H \to H + \Delta H` where

    .. math:: \Delta H = -\frac{1}{2} \sum_n M_n^2 \int d^3x \left( A_n^0 - M_n^{-2} \nabla \cdot \bm{\Pi}_n - M_n^{-2} J_n^0 \right)^2
        :label: eq_many_vector_fields_hamiltonian_correction

    Moreover, in addition to the integration of :math:`\Abf_n` and :math:`\bm{\Pi}_n` in :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element`, we also integrate over :math:`A_n^0`. This addition doesn't really make a difference to the physics since the integrant in :math:`\Delta H` being a perfect square means that the integration over :math:`A_n^0` will only introduce an insignificant field-independent factor to the matrix element.

    Combining :eq:`eq_many_vector_fields_hamiltonian` and :eq:`eq_many_vector_fields_hamiltonian_correction` together, we have

    .. math::

        H + \Delta H = \int d^3x \sum_n \left(
            \frac{1}{2} \bm{\Pi}_n^2 + A_n^0 \nabla \cdot \bm{\Pi}_n + \frac{1}{2} (\nabla \times \Abf_n)^2 + \frac{1}{2} M_n^2 A_n^2 - J_n \cdot A_n
        \right)

    We see that the integrand is still quadratic in :math:`\bm{\Pi}`, whose integration, according to the Gaussian integral formula, can be done by replacing :math:`\bm{\Pi}_n` with the solution to :eq:`eq_path_integral_stationary_p_bar`, which reads

    .. math:: \dot{\Abf}_n = \bm{\Pi}_n - \nabla A_n^0 \iff \bm{\Pi}_n = \dot{\Abf}_n + \nabla A_n^0

    One can then verify that the Legendre transformed quantity

    .. math:: -H - \Delta H + \sum_n \int d^3x~\dot{\Abf}_n \cdot \bm{\Pi}_n

    indeed recovers the original Lagrangian density :eq:`eq_path_integral_many_vector_fields_lagrangian`.

We see from the above examples that it's far from obvious to choose the correct Lagrangian in :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian`. Moreover, the choice of canonical fields may not be the initial :math:`q`\s. Examples of this kind include the vector fields discussed above as well as QED which will be discussed later. Under these considerations, let's rewrite :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian` as follows

.. math::
    :label: eq_path_integral_operator_vacuum_matrix_element_lagrangian_final_form

    &\braket{\VAC, \op{out} | T\{\Oscr_A(\Psi_A(t_A)), \Oscr_B(\Psi_B(t_B)), \cdots\} | \VAC, \op{in}} \\
        &\quad \propto \int \prod_{\tau, \xbf, n} d\psi_n(\tau, \xbf)~\Oscr_A(\psi(t_A)) \Oscr_B(\psi(t_B)) \cdots \\
        &\qquad \times \exp\left(
            \ifrak \int_{-\infty}^{\infty} d\tau \left( L(\psi(\tau), \dot{\psi}(\tau)) + \ifrak\epsilon\text{ terms} \right)
        \right)

where the field-independent constants :math:`|\Nscr|^2` and the part of :math:`\det(2\pi\ifrak\Ascr)^{-1/2}` are suppressed into the proportionality, and the field-dependent part of :math:`\det(2\pi\ifrak\Ascr)^{-1/2}` is absorbed into the Lagrangian. In addition, the dependence of the :math:`\psi`-fields on the right-hand-side on :math:`A, B, \cdots`, is suppressed into the index :math:`n` in the product measure.


Path integral derivation of Feynman rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vacuum expectation value of a time-ordered product of operators given by :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian_final_form` can be evaluated by Feynman diagrams, assuming the propagators have been worked out. However these diagrams may not all be connected. In particular, there is a set of :math:`2`-component diagrams: one of them consists of vertices from only the timed-ordered operators, and the other consists of vertices from the interaction density. Such diagrams can be gotten rid of by considering the following normalized vacuum expectation value

.. math::
    :label: eq_path_integral_defn_m_general

    M_{\ell_A, \ell_B, \cdots}(x_A, x_B, \cdots) \coloneqq \frac{
        \braket{\VAC, \op{out} | T\{\Psi_{\ell_A}(x_A), \Psi_{\ell_B}(x_B), \cdots\} | \VAC, \op{in}}
    }{
        \braket{\VAC, \op{out} | \VAC, \op{in}}
    }

Here a few notations have changed from the previous sections. Firstly, we've used :math:`\ell_A, \ell_B, \cdots`, instead of :math:`A, B, \cdots`, to label the time-ordered operators, which will allow us to unify the labels by :math:`\ell`. Secondly, the argument of the fields has changed from time such as :math:`t_A` to spacetime coordinates :math:`x_A`, which of courses contain :math:`t_A` as its time-component.

Now if the Hamiltonian is quadratic in the :math:`P`-operators as discussed in :ref:`sec_lagrangian_version_of_the_path_integral`, then :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian_final_form` implies that :eq:`eq_path_integral_defn_m_general` can be rewritten as

.. math::
    :label: eq_path_integral_defn_m_quadratic

    M_{\ell_A, \ell_B, \cdots}(x_A, x_B, \cdots) = \frac{
        \int \prod_{x, \ell} d\psi_{\ell}(x)~\psi_{\ell_A}(x_A) \psi_{\ell_B}(t_B) \cdots e^{\ifrak I[\psi]}
    }{
        \int \prod_{x, \ell} d\psi_{\ell}(x)~e^{\ifrak I[\psi]}
    }

where

.. math:: I[\psi] = \int_{-\infty}^{\infty} d\tau \left( L(\psi(\tau), \dot{\psi}(\tau)) + \ifrak\epsilon \text{ terms} \right)

is the action.

Suppose, in the same vein as discussed in :ref:`sec_perturbation_theory_of_s_matrix` and specifically :eq:`eq_defn_v_by_density`, the Lagrangian is given by a density :math:`\Lscr`. Then following the philosophy of perturbation theory, let's write it as the sum of a free part :math:`\Lscr_0` and an interacting part :math:`\Lscr_1`. In other words

.. math:: L(\psi(\tau), \dot{\psi}(\tau)) = \int d^3x \left(
        \Lscr_0(\psi(\tau, \xbf), \p_{\mu} \psi(\tau, \xbf)) + \Lscr_1(\psi(\tau, \xbf), \p_{\mu} \psi(\tau, \xbf))
    \right)

which, in turn, implies that

.. math::
    :label: eq_path_integral_defn_free_and_interacting_actions

    I[\psi] &= I_0[\psi] + I_1[\psi] \\
    I_0[\psi] &= \int d^4x \left( \Lscr_0(\psi, \p_{\mu} \psi) + \ifrak\epsilon\text{ terms} \right) \\
    I_1[\psi] &= \int d^4x~\Lscr_1(\psi, \p_{\mu} \psi)

Such decomposition then allows us to write the exponential term in :eq:`eq_path_integral_defn_m_quadratic` in the following form

.. math::
    :label: eq_path_integral_expand_interaction_action

    \exp(\ifrak I[\psi]) &= \exp(\ifrak I_0[\psi]) \exp(\ifrak I_1[\psi]) \\
        &= \exp(\ifrak I_0[\psi]) \sum_{N=0}^{\infty} \frac{\ifrak^N}{N!} (I_1[\psi])^N

where we've also expanded the second exponential of the interaction action. The reason to do so, or rather, to keep the first exponential of the free action, is that :math:`I_0[\psi]` is typically, and will be assumed to be, quadratic. Indeed, an explicit example was worked out for scalar field in :eq:`eq_path_integral_vacuum_expectation_value_scalar_field`, as long as we ignore the term :math:`e^{-\epsilon |\tau|}` which spoils the quadraticity only in higher orders of :math:`\epsilon`.

If we write

.. math:: I_0[\psi] = -\frac{1}{2} \int d^4x~d^4x' \sum_{\ell, \ell'} \Dscr_{x \ell, x' \ell'} \psi_{\ell}(x) \psi_{\ell'}(x')
    :label: eq_path_integral_quadratic_free_action

then according to :eq:`eq_path_integral_expand_interaction_action`, both the denominator and the numerator of :eq:`eq_path_integral_defn_m_quadratic` are sums of integrals of the following form

.. math:: \Iscr_{\ell_1, \ell_2, \cdots}(x_1, x_2, \cdots) \coloneqq \int \prod_{x, \ell} d\psi_{\ell}(x)~e^{\ifrak I_0[\psi]} \psi_{\ell_1}(x_1) \psi_{\ell_2}(x_2) \cdots
    :label: eq_path_integral_feynman_rules_generic_integral

where :math:`I_0[\psi]` is quadratic. In the case of finite-dimensional integrals, this is a well-known extension of the Gaussian Integral Formula discussed above by integration-by-parts. More formally, this is known as `Wick's theorem <https://en.wikipedia.org/wiki/Isserlis%27s_theorem>`__ which we recall as follows

    **Wick's theorem**

    .. math::

        &\int \prod_r d\xi_s~\xi_{s_1} \xi_{s_2} \cdots \xi_{s_{2N}} \exp\left( -\frac{\ifrak}{2} \sum_{s, r} \Dscr_{sr} \xi_s \xi_r \right) \\
            &\quad = \left( \det(\ifrak\Dscr / 2\pi) \right)^{-1/2} \sum_{\substack{\text{pairings} \\ \text{of } s_1, \cdots, s_{2N}}}
                \prod_{\text{pairs}} \left(-\ifrak\Dscr^{-1}\right)_{\text{paired indices}}

    where :math:`\Dscr` is a real, symmetric, positive matrix.

Applying Wick's theorem to :eq:`eq_path_integral_feynman_rules_generic_integral` we get

.. math::

    \Iscr_{\ell_1, \ell_2, \cdots}(x_1, x_2, \cdots) = (\det(\ifrak\Dscr / 2\pi))^{-1/2} \sum_{\substack{\text{pairings} \\ \text{of fields}}}
        ~\prod_{\text{pairs}} \left( -\ifrak\Dscr^{-1} \right)_{\text{paired fields}}

where :math:`\Dscr` is given by :eq:`eq_path_integral_quadratic_free_action`. Observe that this evaluation, besides the unimportant field-independent factor :math:`(\det(\ifrak\Dscr / 2\pi))^{-1/2}`, can be thought of as a sum over Feynman diagrams where the edges are paired fields that come from either the expansion of :math:`e^{\ifrak I_1[\psi]}` or the timed-ordered operators in the denominator of :eq:`eq_path_integral_defn_m_quadratic`. Moreover, the "propagator" :math:`-\ifrak\Delta` can be defined as follows

.. math:: \Delta_{\ell_1, \ell_2}(x_1, x_2) \coloneqq \Dscr^{-1}_{x_1 \ell_1, x_2 \ell_2}
    :label: eq_path_integral_defn_propagator

To invert :math:`\Dscr` in spacetime coordinates, let's rewrite :eq:`eq_path_integral_defn_propagator` as an integral equation as follows

.. math:: \int d^4 x_2 \sum_{\ell_2} \Dscr_{x_1 \ell_1, x_2 \ell_2} \Delta_{\ell_2, \ell_3}(x_2, x_3) = \delta^4(x_1-x_3) \delta_{\ell_1 \ell_3}

Assuming translation-invariance of the theory, it follows that :math:`\Dscr` can be written as a Fourier transform as follows

.. math:: \Dscr_{x_1 \ell_1, x_2 \ell_2} \eqqcolon (2\pi)^{-4} \int d^4p~e^{\ifrak p \cdot (x_1-x_2)} \Dscr_{\ell_1 \ell_2}(p)
    :label: eq_path_integral_d_matrix_translation_invariance

which, in turn, implies

.. math:: \Delta_{\ell_1 \ell_2}(x_1, x_2) = (2\pi)^{-4} \int d^4p~e^{\ifrak p \cdot (x_1-x_2)} \Dscr^{-1}_{\ell_1 \ell_2}(p)
    :label: eq_path_integral_propagator_as_d_inverse

We conclude the discussion with an example.

Scalar field
    Recall from :eq:`eq_free_real_scalar_field_lagrangian` that the free Lagrangian density takes the following form

    .. math:: \Lscr_0 = -\frac{1}{2} \p_{\mu} \phi \p^{\mu} \phi - \frac{1}{2} m^2 \phi^2

    It follows then from :eq:`eq_path_integral_vacuum_expectation_value_scalar_field` and :eq:`eq_path_integral_defn_free_and_interacting_actions` that the free action :math:`I_0[\phi]`, up to the first order of :math:`\epsilon`, takes the following form

    .. math::

        &I_0[\phi] \\
            &= -\frac{1}{2} \int d^4x \left( \p_{\mu} \phi \p^{\mu} \phi + m^2 \phi^2 \right)
            + \frac{1}{2} \ifrak\epsilon \int dt \int d^3x~d^3x'~\Escr(\xbf, \xbf') \phi(t, \xbf) \phi(t, \xbf') \\
            &= -\frac{1}{2} \int d^4x~d^4x' \left(
                    \delta^4(x-x') (\p_{\mu} \phi \p^{\mu} \phi + m^2 \phi^2) - \ifrak\epsilon \delta(t-t') \Escr(\xbf, \xbf') \phi(x) \phi(x')
                \right) \\
            &= -\frac{1}{2} \int d^4x~d^4x' \left( \frac{\p^2}{\p x^{\mu} \p x'_{\mu}} \delta^4(x-x') + m^2 \delta^4(x-x') -\ifrak\epsilon \delta(t-t') \Escr(\xbf, \xbf') \right) \phi(x) \phi(x') \\
            &= -\frac{1}{2} \int d^4x~d^4x' \left( (2\pi)^{-4} \int d^4p~e^{\ifrak p \cdot (x-x')} \left( p^2 + m^2 - \ifrak\epsilon E(\pbf) \right) \right) \phi(x) \phi(x')

    where in the last equality we've also used :eq:`eq_path_integral_scalar_field_curly_e`.

    Comparing with :eq:`eq_path_integral_quadratic_free_action` and :eq:`eq_path_integral_d_matrix_translation_invariance`, we find

    .. math:: \Dscr(p) = p^2 + m^2 - \ifrak \epsilon E(\pbf)

    and therefore the propagator

    .. math:: \Delta(x, y) = (2\pi)^{-4} \int d^4p~e^{\ifrak p \cdot (x-y)} \left( p^2 + m^2 - \ifrak\epsilon E(\pbf) \right)^{-1}

    according to :eq:`eq_path_integral_propagator_as_d_inverse`. This recovers the Feynman propagator defined by :eq:`eq_defn_feynman_propagator` and evaluated in :eq:`eq_feynman_propagator_as_momentum_space_integral`.


Path Integrals for Fermions
---------------------------

We'll develop the path integral formalism for fermions in parallel to the theory for bosons. The starting point is the commutation relations between Schrödinger-picture canonical variables

.. math::
    :label: eq_path_integral_fermionic_commutation_relation

    \{ Q_a, P_b \} &= \ifrak \delta_{ab} \\
    \{ Q_a, Q_b \} &= \{ P_a, P_b \} = 0

where the curly bracket denotes the anti-commutator. This is to be compared with the bosonic commutation relations :eq:`eq_pib_canonical_commutation_relations`. As in the bosonic case, the indices :math:`a, b` will be replaced by spacetime coordinates as we transit specifically to quantum field theory. A key difference, which will be discussed in more detail in :ref:`sec_fermionic_zero_states`, is that, unlike the bosonic canonical variables, the fermionic :math:`Q` and :math:`P` operators are *not* Hermitian.

.. _sec_fermionic_zero_states:

Fermionic zero states
^^^^^^^^^^^^^^^^^^^^^

It follows from :eq:`eq_path_integral_fermionic_commutation_relation` that

.. math:: Q_a^2 = P_a^2 = 0

Hence there must exist a ket-state :math:`\ket{0}` and a bra-state :math:`\bra{0}` such that

.. math:: Q_a \ket{0} = \bra{0} P_a = 0
    :label: eq_path_integral_fermion_zero_states_annihilated_by_q_and_p

Indeed they can be explicitly constructed as follows

.. math::

    \ket{0} &\propto \left( \prod_a Q_a \right) \ket{f} \\
    \bra{0} &\propto \bra{g} \left( \prod_a P_a \right)

where :math:`\ket{f}` and :math:`\bra{g}` can be any states that makes the right-hand-sides nonzero. In particular the zero states are not in general unique. It turns out that in the absence of bosonic degrees of freedom, the zero states are unique up to a scalar, which can be chosen to satisfy the following normalization property

.. math:: \braket{0 | 0} = 1
    :label: eq_path_integral_zero_state_normalization

.. note::

    The condition :eq:`eq_path_integral_fermion_zero_states_annihilated_by_q_and_p` may seem a bit strange given that the fermionic :math:`Q` and :math:`P` operators are completely interchangeable in light of :eq:`eq_path_integral_fermionic_commutation_relation`. However, it cannot be the case that :math:`Q_a \ket{0} = \bra{0} Q_a = 0` since it would imply :math:`\braket{0 | \{Q_a, P_b\} | 0} = 0` in contradiction with :eq:`eq_path_integral_zero_state_normalization`.

    Indeed, the relationship between :math:`Q` and :math:`P` operators may vary. In Dirac's theory of spin-:math:`1/2` particles, we have :math:`Q_a^{\dagger} = -\ifrak P_a` in light of :eq:`eq_dirac_field_defn_conjugate_pi` and :eq:`eq_dirac_field_psi_field_bar` (cf. :eq:`eq_dirac_field_defn_gamma_matrices` and :eq:`eq_dirac_field_beta_matrix`). In the theory of `ghost field <https://en.wikipedia.org/wiki/Ghost_(physics)>`__,  on the other hand, the :math:`Q` and :math:`P` operators are not related at all.

Fermionic eigenstates
^^^^^^^^^^^^^^^^^^^^^

In light of :eq:`eq_path_integral_fermion_zero_states_annihilated_by_q_and_p`, we can think :math:`P` as the creation operators and :math:`Q` as the annihilation operators. A complete basis of the states can then be obtained from :math:`\ket{0}` by applying an arbitrary number of :math:`P` operators as follows

.. math:: \ket{a_1, a_2, \cdots, a_N} \coloneqq P_{a_1} P_{a_2} \cdots P_{a_N} \ket{0}
    :label: eq_path_integral_defn_fermionic_ket_state

Note that the basis state is anti-symmetric in the following sense

.. math:: \ket{a_1, \cdots, a_i, a_{i+1}, \cdots, a_N} = -\ket{a_1, \cdots, a_{i+1}, a_i, \cdots, a_N}
    :label: eq_path_integral_fermion_state_is_antisymmetric

It follows from :eq:`eq_path_integral_fermion_zero_states_annihilated_by_q_and_p` that

.. math::

    Q_a \ket{a_1, a_2, \cdots, a_N} = \begin{cases}
        (-1)^{k+1} \ifrak~\ket{a_1, a_2, \cdots, \hat{a}_k, \cdots, a_N} & \text{ if } a = a_k \\
        0 & \text{ if } a \notin \{a_1, a_2, \cdots, a_N\}
    \end{cases}

where :math:`\hat{a}_k` means that it's removed from the sequence.

Similarly, the dual basis can be obtained from :math:`\bra{0}` as follows

.. math:: \bra{a_1, a_2, \cdots, a_N} \coloneqq \bra{0} (-\ifrak Q_{a_N}) \cdots (-\ifrak Q_{a_2}) (-\ifrak Q_{a_1})
    :label: eq_path_integral_defn_fermionic_bra_state

The reason to define the dual vector this way is to realize the following normalization condition

.. math::
    :label: eq_path_integral_fermionic_states_normalization

    \braket{b_1, b_2, \cdots, b_M | a_1, a_2, \cdots, a_N} &= \braket{0 | (-\ifrak Q_{b_M}) \cdots (-\ifrak Q_{b_1}) P_{a_1} \cdots P_{a_N} | 0} \\
        &= \begin{cases}
            0 & \text{ if } \{ b_1, b_2, \cdots, b_M \} \neq \{ a_1, a_2, \cdots, a_N \} \text{ as sets} \\
            1 & \text{ if } M=N \text{ and } b_1 = a_1, b_2 = a_2, \cdots, b_M = a_N
        \end{cases}

The cases when :math:`\{ b_1, b_2, \cdots, b_M \}` is a permutation of :math:`\{ a_1, a_2, \cdots, a_N \}` can be covered using :eq:`eq_path_integral_fermion_state_is_antisymmetric`.

The issue with the ket and bra-states defined by :eq:`eq_path_integral_defn_fermionic_ket_state` and :eq:`eq_path_integral_defn_fermionic_bra_state`, respectively, is that they are not eigenstates of :math:`Q` or :math:`P`. In fact, in sharp contrast to the bosonic case (cf. :eq:`eq_path_integral_bosonic_simultaneous_q_eigenstate`), there cannot be *any* eigenstate of, say, all :math:`Q` operators with nonzero (numeric) eigenvalues in the following sense

.. math::
    :label: eq_path_integral_fermionic_q_eigenstate

    Q_a \ket{q} &= q_a \ket{q} \\
    \bra{q} Q_a &= \bra{q} q_a

Indeed, the fermionic commutation relation :eq:`eq_path_integral_fermionic_commutation_relation` would demand

.. math:: q_a q_b + q_b q_a = 0

which cannot be satisfied if :math:`q_a, q_b` are nonzero complex numbers. It turns out that the solution to this difficulty, which may seem to be artificial, is to introduce a new set of "numbers" :math:`q_a`, known as `Grassmann numbers <https://en.wikipedia.org/wiki/Grassmann_number>`__ which satisfy the following anti-commutation relations

.. math:: \{ q_a, q_b \} = \{ q_a, Q_b \} = \{ q_a, P_b \} = 0

Now the fermionic eigenstate equation :eq:`eq_path_integral_fermionic_q_eigenstate` as well as its dual can be solved by the following

.. math::
    :label: eq_path_integral_defn_fermionic_q_eigenstate

    \ket{q} &\coloneqq \exp\left( -\ifrak \sum_a P_a q_a \right) \ket{0} \\
    \bra{q} &\coloneqq \bra{0} \left( \prod_a Q_a \right) \exp\left( \ifrak \sum_a P_a q_a \right)

where the exponential is defined using its Taylor expansion.

.. warning::

   1. In the definition of :math:`\bra{q}` there is a sign ambiguity depending on the ordering of the product of the :math:`Q` operators.
   2. The ket-state :math:`\ket{q}` (e.g. :math:`\ket{0}`) is not necessarily the adjoint of the corresponding bra-state :math:`\bra{q}` (e.g. :math:`\bra{0}`) since :math:`Q` is not Hermitian.

.. dropdown:: Verification of the fermionic eigenstate and its dual
    :animate: fade-in-slide-down
    :icon: unlock

    Let's first verify :math:`\ket{q}` given by :eq:`eq_path_integral_defn_fermionic_q_eigenstate` indeed satisfies :eq:`eq_path_integral_fermionic_q_eigenstate` as follows

    .. math::

        (Q_a - q_a) \ket{q} &= (Q_a - q_a) \exp(-\ifrak P_a q_a) \exp\left( -\ifrak \sum_{b \neq a} P_b q_b \right) \ket{0} \\
            &= (Q_a - q_a) (1 - \ifrak P_a q_a) \exp\left( -\ifrak \sum_{b \neq a} P_b q_b \right) \ket{0} \\
            &= - \ifrak (Q_a P_a - \ifrak) q_a \exp\left( -\ifrak \sum_{b \neq a} P_b q_b \right) \ket{0} \\
            &= \ifrak P_a Q_a q_a \exp\left( -\ifrak \sum_{b \neq a} P_b q_b \right) \ket{0} = 0

    The dual eigenstate :math:`\bra{q}` can be verified as follows

    .. math::

        \bra{q} (Q_a - q_a) &= \bra{0} \left( \prod_a Q_a \right) \exp\left( \ifrak\sum_{b \neq a} P_b q_b \right) \exp(\ifrak P_a q_a)(Q_a - q_a) \\
            &= \bra{0} \left( \prod_a Q_a \right) \exp\left( \ifrak\sum_{b \neq a} P_b q_b \right) (1 + \ifrak P_a q_a)(Q_a - q_a) \\
            &= \bra{0} \left( \prod_a Q_a \right) \exp\left( \ifrak\sum_{b \neq a} P_b q_b \right) (-\ifrak P_a Q_a q_a - q_a) \\
            &= \bra{0} \left( \prod_a Q_a \right) \exp\left( \ifrak\sum_{b \neq a} P_b q_b \right) \ifrak Q_a P_a q_a = 0

Moreover the scalar product of the ket and bra :math:`Q`-eigenstates can be evaluated as follow

.. math::
    :label: eq_path_integral_q_scalar_product

    \braket{q' | q} &= \braket{0 | \left( \prod_a Q_a \right) \exp\left( \ifrak \sum_b P_b \left( q'_b - q_b \right) \right) | 0} \\
        &= \braket{0 | \left( \prod_a Q_a \right) \prod_b \left( 1 + \ifrak P_b (q'_b - q_b) \right) | 0} \\
        &= \prod_a \left( q_a - q'_a \right)

where in the last step, we've used :eq:`eq_path_integral_fermionic_commutation_relation` to move the :math:`Q` operators to the right of the :math:`P` operators. Though not obvious at the moment, the right-hand-side will work as a delta function in fermionic integrals (cf. the bosonic case :eq:`eq_path_integral_bosonic_q_basis_orthogonality`).

The eigenstate of the :math:`P` operators satisfying

.. math::

    P_a \ket{p} &= p_a \ket{p} \\
    \bra{p} P_a &= \bra{p} p_a

can be constructed in a way similar to :eq:`eq_path_integral_defn_fermionic_q_eigenstate` as follows

.. math::
    :label: eq_path_integral_defn_fermionic_p_eigenstate

    \ket{p} &= \exp\left( -\ifrak\sum_a Q_a p_a \right)\left( \prod_a P_a \right) \ket{0} \\
    \bra{p} &= \bra{0} \exp\left( \ifrak\sum_a Q_a p_a \right)

where the order of the product :math:`\prod_a P_a` is, by convention, the same as the one in :eq:`eq_path_integral_defn_fermionic_q_eigenstate`.

The scalar product between the :math:`P`-eigenstates can be similarly evaluated to the following

.. math:: \braket{p' | p} = \prod_a \left( p'_a - p_a \right)
    :label: eq_path_integral_p_scalar_product

In analogy to :eq:`eq_pib_schrodinger_picture_qp_formula`, let's calculate the scalar products between :math:`Q` and :math:`P`-eigenstates as follows

.. math::
    :label: eq_path_integral_fermionic_qp_scalar_product

    \braket{q | p} &= \braket{q | \exp\left( -\ifrak \sum_a Q_a p_a \right) \left( \prod_a P_a \right) | 0} \\
        &= \exp\left( -\ifrak \sum_a q_a p_a \right) \braket{q | \prod_a P_a | 0} \\
        &= \exp\left( -\ifrak \sum_a q_a p_a \right) \braket{0 | \left( \prod_a Q_a \right) \exp\left( \ifrak \sum_a P_a q_a \right) \left( \prod_a P_a \right) | 0} \\
        &= \exp\left( -\ifrak \sum_a q_a p_a \right) \braket{0 | \left( \prod_a Q_a \right) \left( \prod_a P_a \right) | 0} \\
        &= \ifrak^N (-1)^{N(N+1)/2} \exp\left( -\ifrak \sum_a q_a p_a \right)

where :eq:`eq_path_integral_defn_fermionic_bra_state`, :eq:`eq_path_integral_fermionic_states_normalization`, and :eq:`eq_path_integral_zero_state_normalization` are used in the last equality. Here :math:`N` is the number of :math:`Q_a`\s, which is the same as the number of :math:`P_a`\s. Similarly, but more simply, we have

.. math::
    :label: eq_path_integral_fermionic_pq_scalar_product

    \braket{p | q} &= \braket{p | \exp\left( -\ifrak\sum_a P_a q_a \right) | 0} \\
        &= \exp\left( -\ifrak\sum_a p_a q_a \right) \braket{p | 0} \\
        &= \exp\left( -\ifrak\sum_a p_a q_a \right) \braket{0 | \exp\left( \ifrak\sum_a Q_a p_a \right) | 0} \\
        &= \exp\left( -\ifrak\sum_a p_a q_a \right)

We end this section with the note that the states :math:`\ket{q}` are complete in the following sense. If we expand :math:`\ket{q}` in :eq:`eq_path_integral_defn_fermionic_q_eigenstate` as a power series in products of the :math:`q_a`\s, then the coefficients span the whole space of states defined by :eq:`eq_path_integral_defn_fermionic_ket_state`.

Fermionic calculus
^^^^^^^^^^^^^^^^^^

Since the eigenvalues of fermionic eigenstates are Grassmann numbers rather than ordinary (complex) numbers, we need a framework to do calculus, in particular integration, for functions of Grassmann variables. It turns out that the fermionic integration can be formalized as the so-called `Berezin integration <https://en.wikipedia.org/wiki/Berezin_integral>`__, which can be determined by just two rules. Writing :math:`\xi` for generic Grassmann variables, the first rule consists of the evaluation of the integral on a single monomial

.. math:: \int \left( d\xi_N \cdots d\xi_2 d\xi_1 \right) \xi_1 \xi_2 \cdots \xi_N \xi_{N+1} \cdots \xi_M = \xi_{N+1} \cdots \xi_M
    :label: eq_path_integral_berezin_integral_monomial

and the second rule states that the integral is linear in both summation and multiplication by ordinary numbers. Here the ordering of the "differentials" :math:`d\xi_i` in :eq:`eq_path_integral_berezin_integral_monomial` is made so that the integral can be evaluated in steps as follows

.. math:: \int d\xi_N \cdots d\xi_2 d\xi_1~f(\xi) = \int d\xi_N \cdots \int d\xi_2 \int d\xi_1~f(\xi)

It turns out to be very convenient to introduce yet another anti-commutativity relation as follows

.. math:: \{\xi_i, d\xi_j\} = 0
    :label: eq_path_integral_xi_dxi_anti_commute

so we can move the integrand to the left of the "volume element" :math:`\prod_n d\xi_n` at the cost of a sign. Under this convention, it's straightforward to show that given an arbitrary function :math:`g(\xi')` of Grassmann variables that are *not* integrated, the following two formulae hold

.. math::

    \int \left( \prod_n d\xi_n \right) \left( f(\xi) g(\xi') \right) &= \left( \int \left( \prod_n d\xi_n \right) f(\xi) \right) g(\xi') \\
    \int g(\xi') \left( \prod_{n=1}^N d\xi_n \right) f(\xi) &= \int \left( \prod_{n=1}^N d\xi_n \right) \left( g((-1)^N \xi') f(\xi) \right) \\
        &= g(\xi') \int \left( \prod_{n=1}^N d\xi_n \right) f(\xi)

Another important formula in Berezin integration is to describe how the integral transforms under a (linear) change of variables. Consider the following transformation

.. math:: \xi_n \to \xi'_n = \sum_m \Sscr_{nm} \xi_m

where :math:`\Sscr = \left( \Sscr_{nm} \right)` is a non-singular matrix of ordinary numbers. It follows that

.. math:: \prod_n \xi'_n = \left( \det\Sscr \right) \prod_n \xi_n

and henceforth

.. math:: \int \left( \prod_n d\xi'_n \right) f = \left( \det\Sscr \right)^{-1} \int \left( \prod_n d\xi_n \right) f
    :label: eq_path_integral_fermionic_change_of_variable_formula

As an application of this formalism, we'll establish a fermionic analog of the completeness condition :eq:`eq_path_integral_bosonic_q_basis_completeness`. First, note that any state :math:`\ket{f}` can be written as an integral

.. math:: \ket{f} = \int \left( \prod_a dq_a \right) \ket{q} f(q)
    :label: eq_path_integral_f_state_as_integral

where :math:`f(q)` is a polynomial in the Grassmann variables :math:`q`.

.. dropdown:: Verification of :eq:`eq_path_integral_f_state_as_integral`
    :animate: fade-in-slide-down
    :icon: lock

    Rewrite :eq:`eq_path_integral_fermionic_q_eigenstate` as follows

    .. math::

        \ket{q} = \exp\left( -\ifrak\sum_a P_a q_a \right) \ket{0}
            = \left( \prod_a e^{-\ifrak P_a q_a} \right) \ket{0}
            = \left( \prod_a \left( 1 - \ifrak P_a q_a \right) \right) \ket{0}

    so that the right-hand-side is a linear combination of basis states

    .. math:: \ket{a_1, a_2, \cdots, a_k} = P_{a_1} P_{a_2} \cdots P_{a_k} \ket{0}

    whose coefficient is :math:`q_{a_1} q_{a_2} \cdots q_{a_k}` up to a phase.

    Now if we write

    .. math:: \ket{f} = \sum f_{a_1 a_2 \cdots a_k} \ket{a_1, a_2, \cdots, a_k}

    then each summand proportional to :math:`\ket{a_1, a_2, \cdots, a_k}` can be picked up in the right-hand-side of :eq:`eq_path_integral_f_state_as_integral` by a summand in :math:`f(q)` proportional to

    .. math:: \prod_{a \notin \{a_1, a_2, \cdots, a_k\}} q_a

It follows from :eq:`eq_path_integral_q_scalar_product` and :eq:`eq_path_integral_xi_dxi_anti_commute` that

.. math::
    :label: eq_path_integral_scalar_product_f_and_q

    \braket{q' | f} &= \int \braket{q' | q} \left( \prod_{n=1}^N dq_n \right) f(q) \\
        &= \int \left( \prod_{n=1}^N \left( q_n - q'_n \right) \right) \left( \prod_{n=1}^N dq_n \right) f(q) \\
        &= (-1)^N \int \left( \prod_{n=1}^N dq_n \right) \left( \prod_{n=1}^N \left( q_n - q'_n \right) \right) f(q) \\
        &= (-1)^N \int \left( \prod_{n=1}^N dq_n \right) \left( \prod_{n=1}^N \left( q_n - q'_n \right) \right) f(q') \\
        &= (-1)^N f(q')

where the easiest way to justify the second-to-last equality is to write :math:`f(q) = f(q' + (q - q'))` and expand it in powers of :math:`q-q'`, so that only the zeroth order term survive due to the product :math:`\prod_n \left( q_n - q'_n \right)` to the left.

Plugging :eq:`eq_path_integral_scalar_product_f_and_q` into :eq:`eq_path_integral_f_state_as_integral` we have

.. math::
    :label: eq_path_integral_fermionic_q_orthogonality

    \ket{f} = (-1)^N \int \left( \prod_{n=1}^N dq_n \right) \ket{q} \braket{q | f}
        \implies 1 = \int \left( \prod_a -dq_a \right) \ketbra{q}{q}

which is the fermionic version of :eq:`eq_path_integral_bosonic_q_basis_completeness`. The same calculation can be done to the :math:`P`-eigenstates to get the following

.. math:: 1 = \int \left( \prod_a dp_a \right) \ketbra{p}{p}
    :label: eq_path_integral_fermionic_p_orthogonality


The general path integral formula and transition to S-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`H` be the (full) Hamiltonian. Then just as in the bosonic case (cf. :eq:`eq_pif_defn_heisenberg_q_and_p`), we define Heisenberg-picture operators

.. math::

    Q_a(t) &= e^{\ifrak Ht} Q_a e^{-\ifrak Ht} \\
    P_a(t) &= e^{\ifrak Ht} P_a e^{-\ifrak Ht}

with right and left-eigenstates defined as follows

.. math::

    \begin{alignat*}{3}
        \ket{t, q} &\coloneqq e^{\ifrak Ht} \ket{q}, \qquad &&\ket{t, p} &&\coloneqq e^{\ifrak Ht} \ket{p} \\
        \bra{t, q} &\coloneqq \bra{q} e^{-\ifrak Ht}, \qquad &&\bra{t, p} &&\coloneqq \bra{p} e^{-\ifrak Ht}
    \end{alignat*}

Here we recall that :math:`H` must contain an even number of fermionic fields and therefore commute with any Grassmann numbers. Just as in the bosonic case :eq:`eq_path_integral_bosonic_time_dependent_q_and_p_orthogonality_and_completeness`, the time-dependent eigenstates satisfy the obviously analogous time-independent scalar product formulae :eq:`eq_path_integral_q_scalar_product`, :eq:`eq_path_integral_p_scalar_product`, :eq:`eq_path_integral_fermionic_qp_scalar_product`, :eq:`eq_path_integral_fermionic_pq_scalar_product` and orthogonality conditions :eq:`eq_path_integral_fermionic_q_orthogonality`, :eq:`eq_path_integral_fermionic_p_orthogonality`.

Assuming :math:`H = H(P, Q)` is arranged so that all the :math:`P` operators lie to the left of the :math:`Q` operators, we can calculate the infinitesimal transition amplitude in parallel to the bosonic case :eq:`eq_pif_infinitesimal_q_progression_expanded` (except for the ordering :math:`Q` and :math:`P` which is merely a matter of convenience) as follows

.. math::

    \braket{\tau+d\tau, q' | \tau, q} &= \braket{\tau, q' | \exp(-\ifrak H(P, Q)) d\tau | \tau, q} \\
        &= \int \prod_a dp_a \braket{\tau, q' | \tau, p} \braket{\tau, p | \exp(-\ifrak H(P, Q)) d\tau | \tau, q} \\
        &= \int \prod_a dp_a \braket{\tau, q' | \tau, p} \braket{\tau, p | \tau, q} \exp(-\ifrak H(p, q) d\tau) \\
        &\propto \int \prod_a dp_a \exp\left( \ifrak\sum_a p_a(q'_a - q_a) - \ifrak H(p, q) d\tau \right)

where in the last quantity we've thrown away an insignificant field-independent phase factor (coming from :eq:`eq_path_integral_fermionic_qp_scalar_product`), and hence the proportionality is used instead of equality.

Now given a sequence of time-ordered operators, the matrix element analogous to the bosonic :eq:`eq_pif_time_ordered_product_matrix_element` is given by

.. math::

    &\braket{t', q' | T\left\{ \Oscr_A(P(t_A), Q(t_A)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | t, q} \\
    &\quad \propto \int_{\substack{q_a(t)=q_a} \\ q_a(t')=q'_a} \prod_{\tau, a} dq_a(\tau) dp_a(\tau)~\Oscr_A(p(t_A), q(t_A)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak\int_t^{t'} d\tau \left( -H(p(\tau), q(\tau)) + \sum_a p_a(\tau) \dot{q}_a(\tau) \right) \right)

Note that, unlike the bosonic case, an extra sign is added to each permutation of the fermionic operators demanded by the time-ordering operator :math:`T`.

Transitioning to quantum field theory, we get the fermionic analog of :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` as follows

.. math::
    :label: eq_path_integral_fermionic_timed_ordered_operators_vacuum_expectation

    &\braket{\VAC, \op{out} | T\left\{ \Oscr_A(P(t_A), Q(t_B)), \Oscr_B(P(t_B), Q(t_B)), \cdots \right\} | \VAC, \op{in}} \\
    &\quad \propto \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf)~\Oscr_A(p(t_A), q(t_B)) \Oscr_B(p(t_B), q(t_B)) \cdots \\
    &\qquad \times \exp\left( \ifrak \int_{-\infty}^{\infty} d\tau \left( -H(p(\tau), q(\tau)) + \int d^3x \sum_m p_m(\tau, \xbf) \dot{q}_m(\tau, \xbf) + \ifrak \epsilon \text{ terms} \right) \right)

where we've, once again, omitted the details of the :math:`\ifrak\epsilon` terms. As in the bosonic case, they come from the vacuum wave functions (cf. :eq:`eq_path_integral_vacuum_expectation_value_scalar_field`).

The next step in the bosonic case, as discussed in :ref:`sec_lagrangian_version_of_the_path_integral`, is to assume that the Hamiltonian :math:`H(p, q)` is quadratic in :math:`p` and hence the integral in :math:`p`-fields can be evaluated using the Gaussian integral formula (cf. :eq:`eq_path_integral_operator_vacuum_matrix_element_lagrangian`). This is *not* the case for fermions. Unlike the bosonic case (cf. :eq:`eq_path_integral_stationary_p_bar`), the canonical conjugate :math:`p` is  unrelated to :math:`\dot{q}`. In fact, for each fermion that carries a nonzero quantum number, there are equal number of :math:`p`\s and :math:`q`\s. In particular, the free Hamiltonian :math:`H_0` is bilinear in :math:`p` and :math:`q` so that

.. math::
    :label: eq_path_integral_fermionic_defn_d_matrix

    &\int_{-\infty}^{\infty} d\tau \left(
        -H_0(p(\tau), q(\tau)) + \int d^3x \sum_m p_m(\tau, \xbf) \dot{q}_m(\tau, \xbf) + \ifrak\epsilon\text{ terms}
    \right) \\
    &\quad = -\sum_{m, n} \int d^4x~d^4y~\Dscr_{x m, y n} p_m(x) q_n(y)

where :math:`\Dscr` is a matrix in ordinary numbers.

Expanding both the time-ordered operators and the interaction Hamiltonian :math:`V = H - H_0` in a power series in :math:`p`\s and :math:`q`\s, the right-hand-side of :eq:`eq_path_integral_fermionic_timed_ordered_operators_vacuum_expectation` becomes a sum of Gaussian-like integrals and can be evaluated as follows

.. math::
    :label: eq_path_integral_fermionic_gaussian_like_integral_evaluation

    &\Iscr_{m_1 n_1 m_2 n_2 \cdots m_N n_N}(x_1, y_1, x_2, y_2, \cdots, x_N, y_N) \\
    &\quad = \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf) p_{m_1}(x_1) q_{n_1}(y_1) \cdots p_{m_N}(x_N) q_{n_N}(y_N) \\
    &\qquad \times \exp\left( -\ifrak \sum_{m, n} \int d^4x~d^4y~\Dscr_{xm, yn} p_m(x) q_n(y) \right) \\
    &\quad \propto \sum_{\sigma \in \op{Perm}_N} (-1)^{\op{sign}(\sigma)} \prod_{k=1}^N \left( -\ifrak \Dscr^{-1} \right)_{x_k m_k,~y_{\sigma(k)} m_{\sigma(k)}}

where :math:`\op{Perm}_N` denotes the permutation group of :math:`N` symbols. This reproduces the Feynman rules if we regard :math:`\left(\Dscr^{-1}\right)_{xm, yn}` as the propagator between :math:`p_m(x)` and :math:`q_n(y)`.

.. dropdown:: Derivation of the integral evaluation in :eq:`eq_path_integral_fermionic_gaussian_like_integral_evaluation`
    :animate: fade-in-slide-down
    :icon: lock

    Consider a generating function defined as follows

    .. math::
        :label: eq_path_integral_fermionic_generating_function

        \Iscr(f, g) &\coloneqq \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf) \\
            &\quad \times \exp\left( -\ifrak\sum_{m, n} \int d^4x~d^4y~\Dscr_{xm, yn} p_m(x) q_n(y) \right. \\
            &\qquad \left. -\ifrak\sum_m \int d^4x~p_m(x) f_m(x) - \ifrak\sum_n \int d^4y~g_n(y) q_n(y) \right)

    where :math:`f, g` are generic Grassmann numbers. Next, consider the following (translational) change of variables

    .. math::

        p'_m(x) &= p_m(x) - \sum_n \int d^4y~g_n(y) \left( \Dscr^{-1} \right)_{yn, xm} \\
        q'_n(y) &= q_n(y) - \sum_m \int d^4x \left( \Dscr^{-1} \right)_{yn, xm} f_m(x)

    Obviously the integral doesn't change by this change of variables. Let's calculate the change of exponential power in :eq:`eq_path_integral_fermionic_generating_function` (without the common factor :math:`-\ifrak`) as follows

    .. math::

        &\sum_{m,n} \int d^4x~d^4y~\Dscr_{xm,yn} p'_m(x) q'_n(y) + \sum_m \int d^4x~p'_m(x) f_m(x) + \sum_n \int d^4y~g_n(y) q'_n(y) \\
        &= \sum_{m,n} \int d^4x~d^4y~\Dscr_{xm,yn} p_m(x) q_n(y) - \sum_{m,n,k} \int d^4x~d^4y~d^4w~\Dscr_{xm,yn} p_m(x) \left(\Dscr^{-1}\right)_{yn,wk} f_k(w) \\
        &\quad - \sum_{m,n,k} \int d^4x~d^4y~d^4w~\Dscr_{xm,yn} g_k(w) \left(\Dscr^{-1}\right)_{wk,xm} q_n(y) \\
        &\quad + \sum_{m,n,k,r} \int d^4x~d^4y~d^4w~d^4v~\Dscr_{xm,yn} g_k(w) \left(\Dscr^{-1}\right)_{wk,xm} \left(\Dscr^{-1}\right)_{yn,vr} f_r(v) \\
        &\quad + \sum_m \int d^4x~p_m(x) f_m(x) - \sum_{m,n} \int d^4x~d^4y~g_n(y) \left(\Dscr^{-1}\right)_{yn,xm} f_m(x) \\
        &\quad + \sum_n \int d^4y~g_n(y) q_n(y) - \sum_{m,n} \int d^4x~d^4y~g_n(y) \left(\Dscr^{-1}\right)_{yn,xm} f_m(x) \\
        &= \sum_{m,n} \int d^4x~d^4y~\Dscr_{xm,yn} p_m(x) q_n(y) - \sum_{m,n} \int d^4x~d^4y~\left(\Dscr^{-1}\right)_{yn,xm} g_n(y) f_m(x)

    It follows that :eq:`eq_path_integral_fermionic_generating_function` can be rewritten as follows

    .. math::
        :label: eq_path_integral_fermionic_generating_function_reevaluated

        \Iscr(f, g) &= \exp\left( \ifrak \sum_{m,n} \int d^4x~d^4y \left(\Dscr^{-1}\right)_{yn,xm} g_n(y) f_m(x) \right) \\
            &\times \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf) \\
            &\times \exp\left( -\ifrak \sum_{m,n} \int d^4x~d^4y~\Dscr_{xm,yn} p_m(x) q_n(y) \right)

    Now the integral :eq:`eq_path_integral_fermionic_generating_function` is called a generating function since if the exponential is expanded as a power series in :math:`f` and :math:`g`, then the coefficient of

    .. math:: f_{m_1} g_{n_1} f_{m_2} g_{n_2} \cdots f_{m_N} g_{n_N}
        :label: eq_path_integral_fermionic_fg_monomial

    is, up to a constant factor independent of :math:`x,y,m,n`, precisely the integral :eq:`eq_path_integral_fermionic_gaussian_like_integral_evaluation`. Expanding :eq:`eq_path_integral_fermionic_generating_function_reevaluated` also in a power series of :math:`f` and :math:`g`, the coefficient of :eq:`eq_path_integral_fermionic_fg_monomial` then gives the final result in :eq:`eq_path_integral_fermionic_gaussian_like_integral_evaluation`.

As an example, consider the free-particle action of spin-:math:`1/2` particles as follows (cf. :eq:`eq_canonical_formalism_dirac_field_lagrangian_density` and :eq:`eq_dirac_field_free_hamiltonian`)

.. math::

    &\int_{-\infty}^{\infty} d\tau \left( -H_0(p(\tau), q(\tau)) + \int d^3x \sum_m p_m(\tau, \xbf) \dot{q}_m(\tau, \xbf) \right) \\
        &\quad = -\int d^4x~\bar{\psi}(x) \left( \gamma^{\mu} \p_{\mu} + M \right) \psi(x)

The canonical variables are given by (cf. :eq:`eq_dirac_field_defn_conjugate_pi`)

.. math::

    q_m(x) &= \psi_m(x) \\
    p_m(x) &= -\left( \bar{\psi}(x) \gamma^0 \right)_m

Comparing with :eq:`eq_path_integral_fermionic_defn_d_matrix` we have

.. math::

    \Dscr_{xm, yn} &= \left( \gamma^0 \left( \gamma^{\mu} \p_{x_{\mu}} + M - \ifrak\epsilon \right) \right)_{mn} \delta^4(x-y) \\
        &= \int \frac{d^4 k}{(2\pi)^4} \left( \gamma^0 \left( \ifrak \gamma^{\mu} k_{\mu} + M -\ifrak\epsilon \right) \right)_{mn} e^{\ifrak k \cdot (x-y)}

where we've omitted the details about the :math:`\ifrak\epsilon` term since it's not important here and can be worked out as the vacuum wave function in the same way as for (bosonic) scalar field (cf. :eq:`eq_path_integral_scalar_field_vacuum_wave_function`). It follows that the propagator can be written as follows

.. math::

    \left(\Dscr^{-1}\right)_{xm, yn} &= \int \frac{d^4k}{(2\pi)^4} \left( \left( \ifrak\gamma^{\mu} k_{\mu} + M -\ifrak\epsilon \right)^{-1} (-\gamma^0) \right)_{mn} e^{\ifrak k \cdot (x-y)} \\
        &= \int \frac{d^4k}{(2\pi)^4} \frac{\left( \left( -\gamma^{\mu} k_{\mu} + M \right) (-\gamma^0) \right)_{mn}}{k^2 + M^2 - \ifrak\epsilon} e^{\ifrak k \cdot (x-y)}

which recovers :eq:`eq_propagator_as_momentum_space_integral_linear` and :eq:`eq_p_polynomial_dirac`, except for that we've calculated here the propagator between :math:`\psi` and :math:`-\bar{\psi} \gamma^0`, rather than :math:`\psi^{\dagger}`.

To illustrate how path-integral method works in this case, let's consider the following interacting Lagrangian density

.. math:: \Lscr = -\bar{\psi} \left( \gamma^{\mu} \p_{\mu} + M + \Gamma \right) \psi

where :math:`\Gamma(x)` represents the interaction of the fermion with an external field. Applying :eq:`eq_path_integral_fermionic_timed_ordered_operators_vacuum_expectation`, we can calculate the vacuum persistence amplitude (in matrix notation) as follows

.. math::
    :label: eq_path_integral_interacting_fermion_vacuum_expectation

    \braket{\VAC, \op{out} | \VAC, \op{in}}_{\Gamma} &\propto \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf, m) \\
        &\quad \times \exp\left( -\ifrak \int d^4x~p^T \gamma^0 (\gamma^{\mu} \p_{\mu} + M + \Gamma - \ifrak\epsilon) q \right) \\
        &= \int \prod_{\tau, \xbf, m} dq_m(\tau, \xbf) \prod_{\tau, \xbf, m} dp_m(\tau, \xbf, m) \\
        &\quad \times \exp\left( -\ifrak \sum_{m, n} \int d^4x~d^4y~p_m(x) q_n(y) \Kscr(\Gamma)_{xm, yn} \right)

where

.. math:: \Kscr(\Gamma)_{xm, yn} \coloneqq \left( \gamma^0 \left( \gamma^{\mu} \p_{x_{\mu}} + M + \Gamma(x) - \ifrak\epsilon \right) \right)_{mn} \delta^4(x-y)
    :label: eq_path_integral_fermionic_interacting_k_matrix

Observe that the following change of variables

.. math:: q'_m(x) \coloneqq \sum_n \int d^4y~\Kscr(\Gamma)_{xm, yn} q_n(y)

will make the integral in :eq:`eq_path_integral_interacting_fermion_vacuum_expectation` independent of :math:`\Gamma`. Hence it follows from :eq:`eq_path_integral_fermionic_change_of_variable_formula` that

.. math:: \braket{\VAC, \op{out} | \VAC, \op{in}}_{\Gamma} \propto \det\Kscr(\Gamma)
    :label: eq_path_integral_fermionic_vacuum_expectation_for_interacting_dirac_fields

To see that this calculation is consistent with the Feynman rules derived in :ref:`sec_momentum_space_feynman_rules`, rewrite :eq:`eq_path_integral_fermionic_interacting_k_matrix` as follows

.. math:: \Kscr(\Gamma) = \Dscr + \Gscr(\Gamma)

where

.. math:: \Gscr(\Gamma)_{xm, yn} \coloneqq \left( \gamma^0 \Gamma(x) \right)_{mn} \delta^4(x-y)

It follows then (cf. :eq:`eq_det_eq_exp_tr_ln`)

.. math::

    \braket{\VAC, \op{out} | \VAC, \op{in}}_{\Gamma} &\propto \det\left( \Dscr \left( 1 + \Dscr^{-1} \Gscr(\Gamma) \right) \right) \\
        &= \left(\det\Dscr\right) \exp\left( \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n} \Tr\left( \Dscr^{-1} \Gscr(\Gamma) \right) \right)

which is exactly what one would expect from a Feynman diagram calculation. Indeed, the Feynman rules demand that each vertex contributes a factor :math:`-\ifrak\Gscr(\Gamma)` and each (internal) edge contributes a factor :math:`-\ifrak\Dscr`. In this example all connected Feynman diagram are simply loops. A loop with :math:`n` vertices then contributes a factor :math:`1/n` accounting for the cyclic symmetry. Finally the extra sign in :math:`(-1)^{n+1}` accounts for the fermionic loop.

It turns out that path integral calculations like :eq:`eq_path_integral_fermionic_vacuum_expectation_for_interacting_dirac_fields` do much more than just recovering perturbative calculations. But we'll only come back to this much later.


Path-integral formulation of QED
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In :ref:`sec_quantum_electrodynamics`, we derived the photon propagator :eq:`eq_qed_photon_line_contribution_covariant` in a rather cumbersome way by working with non-Lorentz-invariant :math:`A`-fields (cf. :eq:`eq_qed_a0_vanishes` and :eq:`eq_qed_a_field_general_solution`). In this section, we will re-derive :eq:`eq_qed_photon_line_contribution_covariant` in a cleaner manner using path integrals.

Recall from :eq:`eq_qed_hamiltonian_vector_form` that the QED Hamiltonian can be written as follows

.. math::

    H = \int d^3x \left( \frac{1}{2} \bm{\Pi}^2_{\bot} + \frac{1}{2} (\nabla \times \Abf)^2 - \Jbf \cdot \Abf \right)
        + V_{\op{Coul}} + H_{\op{matter}}

where the Coulomb potential :math:`V_{\op{Coul}}` is given by :eq:`eq_qed_defn_coulomb_energy`, under the Coulomb gauge condition

.. math:: \nabla \cdot \Abf = 0
    :label: eq_path_integral_qed_a_constraint

and the constraint (cf. :eq:`eq_qed_pi_bot_is_divergence_free`)

.. math:: \nabla \cdot \bm{\Pi}_{\bot} = 0
    :label: eq_path_integral_qed_pi_constraint

The general path integral calculation of the vacuum expectation value of a timed-ordered product of operators, whether the bosonic :eq:`eq_pi_to_s_matrix_general_vacuum_matrix_element` or the fermionic :eq:`eq_path_integral_fermionic_timed_ordered_operators_vacuum_expectation` can be applied to the QED Hamiltonian to give the following

.. math::
    :label: eq_path_integral_qed_operator_vacuum_expectation

    &\braket{T\{\Oscr_A \Oscr_B \cdots\}}_{\VAC} = \int \prod_{x, i} da_i(x) \prod_{x, i} d\pi_i(x) \prod_{x, \ell} d\psi_{\ell}(x)~\Oscr_A \Oscr_B \cdots \\
        &\qquad \times \exp\left( \ifrak \int d^4x \left( \bm{\pi} \cdot \dot{\abf} - \frac{1}{2} \bm{\pi}^2 - \frac{1}{2}(\nabla \times \abf)^2 + \jbf \cdot \abf + \Lscr_M \right) - \ifrak \int dt~V_{\op{Coul}}(t) \right) \\
        &\qquad \times \left( \prod_x \delta(\nabla \cdot \abf(x)) \right) \left( \prod_x \delta(\nabla \cdot \bm{\pi}(x)) \right)

where the lowercase fields represent simply variables of integration, rather than interaction-picture fields as in :ref:`sec_quantum_electrodynamics`, the index :math:`i` runs through the spatial :math:`1,2,3` as usual, the matter field variables are represented by :math:`\psi_{\ell}`, and the last two delta functions are introduced to enforce the constraints :eq:`eq_path_integral_qed_a_constraint` and :eq:`eq_path_integral_qed_pi_constraint`.

First thing to observe is that since :math:`\Lscr_M` doesn't involve :math:`\bm{\pi}` (cf. :eq:`eq_qed_lagrangian_density` and :eq:`eq_canonical_formalism_dirac_field_lagrangian_density`) and the power of the exponential is quadratic in :math:`\bm{\pi}`, it can be integrated out by the Gaussian integral formula :eq:`eq_gaussian_integral_formula`. More explicitly, it amounts to substitute :math:`\bm{\pi}` with the solution to the stationary point of the quadratic power, which is :math:`\bm{\pi} = \dot{\abf}`, so we can rewrite :eq:`eq_path_integral_qed_operator_vacuum_expectation` as follows

.. math::

    &\braket{T\{\Oscr_A \Oscr_B \cdots\}}_{\VAC} = \int \prod_{x, i} da_i(x) \prod_{x, \ell} d\psi_{\ell}(x) \Oscr_A \Oscr_B \cdots \\
        &\qquad \times \exp\left( \ifrak \int d^4x \left( \frac{1}{2} \dot{\abf}^2 - \frac{1}{2} \bm{\pi}^2 - \frac{1}{2} (\nabla \times \abf)^2 + \jbf \cdot \abf + \Lscr_M \right) - \ifrak \int dt~V_{\op{Coul}}(t) \right) \\
        &\qquad \times \left( \prod_x \delta(\nabla \cdot \abf(x)) \right)

Next we'd like to restore manifest Lorentz invariance by integrating also over :math:`a_0`. The trick is to consider the following quantity

.. math:: \int d^4x \left( -a_0(x) j_0(x) + \frac{1}{2} \left(\nabla a_0(x)\right)^2 \right)

whose exponential's (path) integral over :math:`a_0(x)` can be done by setting :math:`a_0(x)` to the stationary point. In other words :math:`a_0(x)` should solve the following differential equation

.. math:: j_0(x) + \nabla^2 a_0(x) = 0

This is a rather familiar equation (cf. :eq:`eq_qed_poisson_equation_j_and_a`), whose solution, given by :eq:`eq_qed_explicit_solution_of_a0`, is reproduced as follows

.. math:: a_0(t, \xbf) = \int d^3y \frac{j_0(t, \ybf)}{4\pi|\xbf-\ybf|}
