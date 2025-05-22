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
            \int_t^{t'} d\tau \left( -H(q(\tau), p(\tau)) + \sum_a \dot{q}_a(\tau) p_a(\tau) \right)
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


.. rubric:: Footnotes

.. [#qp_product_2pi_factor] It's unclear to me why a factor of :math:`1/\sqrt{2\pi}`, which clearly diminishes the (infinite) product, is inserted in [Wei95]_ page 379 eq. (9.1.12). It might simply be a mistake considering its fermionic counterpart eq. (9.5.27) in page 403.
