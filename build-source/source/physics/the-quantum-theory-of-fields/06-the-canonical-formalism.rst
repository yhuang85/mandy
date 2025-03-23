The Canonical Formalism
-----------------------

The quantum theory we've been developing so far has been based almost solely on the symmetry principles, especially Lorentz symmetries. This is a very satisfying approach since it's logically clean and relies only on the most fundamental principles, however, this is not the way quantum theory historically had been developed. Not surprisingly, the original development of quantum theory is much messier and requires substantial experience in "classical" physics. It's largely based on the so-called *Lagrangian formalism*, which is a readily well-established principle in classical physics and can be "quantized". The main goal of this chapter is to go through this formalism, not for historical sake, but because it offers a particularly convenient way to construct Hamiltonians that generate Lorentz-invariant S-matrices, which has been difficult for us as can be seen in :ref:`sec_feynman_rules_in_momentum_space`.

Canonical variables
^^^^^^^^^^^^^^^^^^^

We've seen in :ref:`sec_quantum_fields_and_antiparticles` a few ways of constructing (Lorentz-invariant) interaction densities. However, we don't have a systematic way to do so. The so-called Lagrangian formalism will not provide a systematic solution either, but it'll allow us to construct more interesting interaction densities (from classical physics theories), to the extent that all known quantum field theories arise in this way! In addition, it'll shed light on the mysterious local terms as for example in :eq:`eq_vector_field_propagator_needs_local_term`, that are needed to compensate for a Lorentz-invariant momentum space propagator.

The offer from the Lagrangian formalism regarding constructing a quantum field theory is the following. Instead of using the  creation and annihilation fields defined by :eq:`eq_defn_annihilation_and_creation_field` to construct the Hamiltonians, we'll use the so-called *canonical variables*, which have particularly simple (equal time) commutation relations. More precisely, it consists of a collection of quantum operators :math:`q_n(t, \xbf)` and its canonical conjugates :math:`p_n(t, \xbf)`, which satisfy the following (anti-)commutation relations

.. math::
	:label: eq_canonical_commutation_relations

	\left[ q_n(t, \xbf), p_{n'}(t, \ybf) \right]_{\pm} &= \ifrak \delta^3(\xbf - \ybf) \delta_{n n'} \\
	\left[ q_n(t, \xbf), q_{n'}(t, \ybf) \right]_{\pm} &= 0 \\
	\left[ p_n(t, \xbf), p_{n'}(t, \ybf) \right]_{\pm} &= 0 \\

where :math:`\pm` correspond to when the particle under question is fermionic or bosonic, respectively.

To see how canonical variables may be constructed from fields considered in :ref:`sec_quantum_fields_and_antiparticles`, let's consider a few examples.

Scalar fields
	Let's start by considering scalar fields of particles that are their own antiparticles. Using notations from :ref:`sec_scalar_field`, it means that :math:`\psi(x) = \psi^{
	\dagger}(x)`, i.e., the field is Hermitian. It follows then from :eq:`eq_scalar_field_commutator` and :eq:`eq_defn_Delta` that

	.. math::

		\left[ \psi(x), \psi(y) \right] = \Delta(x-y) = \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0} \left( e^{\ifrak p \cdot (x-y)} - e^{-\ifrak p \cdot (x-y)} \right)

	where :math:`p_0 = \sqrt{\pbf^2 + m^2}`.

	We claim that the canonical commutation relations :eq:`eq_canonical_commutation_relations` are satisfied by

	.. math::
		:label: eq_defn_q_and_p_scalar_field_self_dual

		q(t, \xbf) &\coloneqq \psi(t, \xbf) \\
		p(t, \xbf) &\coloneqq \dot{\psi}(t, \xbf)

	Indeed, it follows from the following calculations

	.. math::
		:label: eq_canonical_commutators_scalar_field_self_dual

		\begin{alignat*}{3}
			\left[ q(t, \xbf), p(t, \ybf) \right]
				&= \left[ \psi(t, \xbf), \dot{\psi}(t, \ybf) \right]
				&&= -\dot{\Delta}(0, \xbf-\ybf)
				&&= \ifrak \delta^3(\xbf-\ybf) \\
			\left[ q(t, \xbf), q(t, \ybf) \right]
				&= \left[ \psi(t, \xbf), \psi(t, \ybf) \right]
				&&= \Delta(0, \xbf-\ybf)
				&&= 0 \\
			\left[ p(t, \xbf), p(t, \ybf) \right] &= \left[ \dot{\psi}(t, \xbf), \dot{\psi}(t, \ybf) \right] &&= -\ddot{\Delta}(0, \xbf-\ybf) = 0
		\end{alignat*}

	Now for particles that are different from their antiparticles, we must modify :eq:`eq_defn_q_and_p_scalar_field_self_dual` as follows

	.. math::

		q(t, \xbf) &= \psi(t, \xbf) \\
		p(t, \xbf) &= \dot{\psi}^{\dagger}(t, \xbf)

	and note that in this case :math:`\left[ \psi(t, \xbf), \psi(t', \ybf) \right] = 0`, in contrast to the second equation in :eq:`eq_canonical_commutators_scalar_field_self_dual`.

Spin-:math:`1` vector fields
	Consider once again particles that are self-charge-dual. Using notations from :ref:`sec_spin_1_vector_fields`, we recall the commutation relation :eq:`eq_vector_field_commutator` as follows

	.. math::

		\left[ \psi_{\mu}(x), \psi_{\nu}(y) \right] = \left( \eta_{\mu\nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta(x-y)

	The canonical variables in this case can be defined as follows

	.. math::
		:label: eq_defn_q_and_p_vector_field_self_dual

		q_i(t, \xbf) &= \psi_i(t, \xbf) \\
		p_i(t, \xbf) &= \dot{\psi}_i(t, \xbf) - \frac{\p \psi_0(t, \xbf)}{\p x_i}

	where :math:`i=1,2,3`. Indeed, let's calculate the equal-time commutators as follows

	.. math::

		\left[ q_i(t, \xbf), p_j(t, \ybf) \right] &= \left[ \psi_i(t, \xbf), \dot{\psi}_j(t, \ybf) \right] - \left[ \psi_i(t, \xbf), \frac{\p \psi_0(t, \ybf)}{\p y_j} \right] \\
			&= -\left( \eta_{ij} -\frac{\p_i \p_j}{m^2} \right) \dot{\Delta}(0, \xbf-\ybf) - \left. \frac{\p_i \p_0}{m^2} \right|_{t=0} \left( \p_j \Delta(t, \xbf-\ybf) \right) \\
			&= \ifrak \delta^3(\xbf-\ybf) \delta_{ij} \\
		\left[ q_i(t, \xbf), q_j(t, \ybf) \right] &= \left( \eta_{ij} - \frac{\p_i \p_j}{m^2} \right) \Delta(0, \xbf-\ybf) = 0 \\
		\left[ p_i(t, \xbf), p_j(t, \ybf) \right] &= \left[ \dot{\psi}_i(t, \xbf), \dot{\psi}_j(t, \ybf)\right] + \p_{x_i} \p_{y_j} \left[ \psi_0(t, \xbf), \psi_0(t, \ybf) \right] \\
		&\qquad - \p_{x_i} \left[ \psi_0(t, \xbf), \dot{\psi}_j(t, \ybf) \right] - \p_{y_j} \left[ \dot{\psi}_i(t, \xbf), \psi_0(t, \ybf) \right] = 0

	We've omitted some details about the vanishing of the last quantities -- it turns out that the the first and second terms cancel out, and the third and the fourth terms also cancel out.

	In any case, we've constructed three pairs of canonical variables, one for each spatial index. But what about the time index? It turns out that :math:`\psi_0` is *not* an independent variable. Indeed, we can derive from :eq:`eq_defn_q_and_p_vector_field_self_dual`, using :eq:`eq_vector_field_gauge_fixing_condition` and :eq:`eq_klein_gordon`, an expression of :math:`\psi_0` as follows

	.. math::

		p_i = \p_0 \psi_i - \p_i \psi_0
			& \implies \p_i p_i = \p_0 \p_i \psi_i - \p^2_i \psi_0 \\
			& \implies \nabla \cdot \pbf = \p_0 \sum_{i=1}^3 \p_i \psi_i - \sum_{i=1}^3 \p^2_i \psi_0 \\
			& \implies \nabla \cdot \pbf = \p_0^2 \psi_0 - \sum_{i=1}^3 \p_i^2 \psi_0 = -\square \psi_0 \\
			& \implies \psi_0 = -m^{-2} \nabla \cdot \pbf

Spin-:math:`1/2` Dirac fields
	Recall the anti-commutator of Dirac fields :eq:`eq_dirac_field_commutator` as follows

	.. math::

		\left[ \psi_{\ell}(x), \psi^{\dagger}_{\ell'}(y) \right]_+ = \left( (-\gamma^{\mu} \p_{\mu} + m) \beta \right)_{\ell \ell'} \Delta(x-y)

	where :math:`\ell, \ell'` are indexes corresponding to the two spin :math:`z`-component :math:`\pm 1/2`. Assuming that particle under question has distinct antiparticle, i.e., it's not a Majorana fermion, the following holds trivially

	.. math:: \left[ \psi_{\ell}(x), \psi_{\ell'}(y) \right]_+ = 0

	It follows that the canonical variables can be defined by

	.. math::

		q_{\ell}(x) &= \psi_{\ell}(x) \\
		p_{\ell}(x) &= \ifrak \psi^{\dagger}_{\ell}(x)

	Indeed, the only nontrivial (and non-vanishing) anti-commutator can be calculated as follows

	.. math::

		\left[ q_{\ell}(t, \xbf), p_{\ell'}(t, \ybf) \right]_+ &= \ifrak \left[ \psi_{\ell}(t, \xbf), \psi_{\ell'}^{\dagger}(t, \ybf) \right]_+ \\
			&= -\ifrak \left( \gamma^0 \beta \right)_{\ell \ell'} \dot{\Delta}(0, \xbf-\ybf) \\
			&= \ifrak \delta^3(\xbf-\ybf) \delta_{\ell \ell'}

Through these examples, we see that there is no particular pattern in how one may define canonical variables. In fact, one doesn't really define canonical variables in this way either -- they are simply given for granted in the Lagrangian formalism as we will see.

We begin by a general discussion on functionals :math:`F[q(t), p(t)]` of canonical variables, since both Hamiltonians and Lagrangians will be such functionals. A few notes are in order. First we've used a shorthand notation :math:`q(t)` and :math:`p(t)` to denote a collection of canonical variables. Moreover, in writing :math:`q(t)` (and similarly for :math:`p(t)`) we implicitly think of them as fields at a given time. Indeed, as we'll see, the time variable plays an exceptional role in the Lagrangian formalism, in contrast to our mindset so far that space and time are all mixed up in a Lorentz invariant theory. Finally, we've used square bracket to differentiate it from regular functions of spacetime or momentum variables.

At the heart of the Lagrangian formalism lies a variational principle. Hence it's crucial to be able to take infinitesimal variations on :math:`F[q(t), p(t)]`, which we write as follows

.. math::
	:label: eq_infinitesimal_variation_of_functional_of_canonical_variables

	\delta F[q(t), p(t)] = \int d^3 x \sum_n \left( \delta q_n(t, \xbf) \frac{\delta F[q(t), p(t)]}{\delta q_n(t, \xbf)} + \frac{\delta F[q(t), p(t)]}{\delta p_n(t, \xbf)} \delta p_n(t, \xbf) \right)

Here the infinitesimal fields :math:`\delta q_n` and :math:`\delta p_n` are assumed to (anti-)commute with all other fields. Now assuming :math:`F[q(t), p(t)]` is written so that all the :math:`q` fields lie to the left of all the :math:`p` fields, then :eq:`eq_infinitesimal_variation_of_functional_of_canonical_variables` can be realized by the following definition of variational derivatives

.. math::

	\frac{\delta F[q(t), p(t)]}{\delta q_n(t, \xbf)} \coloneqq \ifrak \big[ p_n(t, \xbf), F[q(t), p(t)] \big] \\
	\frac{\delta F[q(t), p(t)]}{\delta p_n(t, \xbf)} \coloneqq \ifrak \big[ F[q(t), p(t)], q_n(t, \xbf) \big]


Hamiltonian and Lagrangian for free fields
++++++++++++++++++++++++++++++++++++++++++

For free fields we have

.. math::
	:label: eq_free_field_q_and_p_time_evolution

	q_n(t, \xbf) &= e^{\ifrak H_0 t} q_n(0, \xbf) e^{-\ifrak H_0 t} \\
	p_n(t, \xbf) &= e^{\ifrak H_0 t} p_n(0, \xbf) e^{-\ifrak H_0 t}

where :math:`H_0` is the free field Hamiltonian, also known as the symmetry generator for the time translation, or the energy operator. However, rather than thinking of it as an abstract operator as we've done so far, we'll (momentarily) make it a functional of canonical variables. With this in mind, we can take the time derivative of :eq:`eq_free_field_q_and_p_time_evolution` as follows

.. math::
	:label: eq_free_field_hamilton_equation_q_and_p_dot

	\begin{alignat*}{2}
		\dot{q}_n(t, \xbf) &= \ifrak \left[ H_0, q_n(t, \xbf) \right] &&= \frac{\delta H_0}{\delta p_n(t, \xbf)} \\
		\dot{p}_n(t, \xbf) &= \ifrak \left[ H_0, p_n(t, \xbf) \right] &&= -\frac{\delta H_0}{\delta q_n(t, \xbf)}
	\end{alignat*}

We recognize these as the quantum analog of `Hamilton's equation of motion <https://en.wikipedia.org/wiki/Hamiltonian_mechanics>`__.

To turn :math:`H_0` into a functional of canonical variables, we first make it a functional of creation and annihilation operators. Remembering that :math:`H_0` is the energy operator, and :math:`p_0 = \sqrt{\pbf^2 + m^2}` is the energy in the :math:`4`-momentum, we can write :math:`H_0` as a diagonal operator as follows

.. math:: H_0 = \sum_{n, \sigma} \int d^3 p~a^{\dagger}(\pbf, \sigma, n) a(\pbf, \sigma, n) \sqrt{\pbf^2 + m^2}
	:label: eq_free_field_hamiltonian_diagonal

For simplicity, let's consider the case of a real scalar field :math:`\psi(x)` given by :eq:`eq_scalar_field_psi_by_creation_and_annihilation_operators` as follows

.. math::

	q(t, \xbf) = \psi(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} \left( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)

The canonical conjugate variable is

.. math::

	p(t, \xbf) = \dot{\psi}(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} (-\ifrak p_0) \left( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)

These look a bit far from :eq:`eq_free_field_hamiltonian_diagonal`. But since :math:`H_0` involves products like :math:`a^{\dagger}(\pbf, \sigma, n) a(\pbf, \sigma, n)`, let's try to square the canonical variables as follows

.. math::

	\int d^3 x~q^2(t, \xbf) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}}
			\Big( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big)
			\Big( e^{\ifrak p' \cdot x} a(\pbf') + e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
		&= \int \frac{d^3 p}{2p_0} \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } + \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)

and

.. math::

	\int d^3 x~p^2(t, \xbf) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}} (-p_0 p'_0)
			\Big( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big) \\
		&\qquad \times \Big( e^{\ifrak p' \cdot x} a(\pbf') - e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
		&= \int \frac{d^3 p}{2p_0} \left( -p_0^2 \right) \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } - \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)

and finally, inspired by the calculations above

.. math::

	\int d^3 x~\left( \nabla q(t, \xbf) \right)^2 &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}} \left( -\pbf \cdot \pbf' \right)
			\Big( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big) \\
		&\qquad \times \Big( e^{\ifrak p' \cdot x} a(\pbf') - e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
		&= \int \frac{d^3 p}{2p_0} \pbf^2 \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } + \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)

Putting these calculations together in a specific way, and using the identity :math:`p_0^2 - \pbf^2 = m^2`, we can eliminate the blue terms as follows

.. math::
	:label: eq_calculate_free_real_scalar_field_hamiltonian

	\frac{1}{2} \int d^3 x \left( p^2 + \left( \nabla q \right)^2 + m^2 q^2 \right)
		&= \frac{1}{2} \int d^3 p~p_0 \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \\
		&= \int d^3 p~p_0 \left( a^{\dagger}(\pbf) a(\pbf) + \frac{1}{2} \delta^3(\pbf-\pbf) \right) \\
		&= H_0 + \blue{ \frac{1}{2} \int d^3 p~p_0 \delta^3(0) }

Here we've encountered for the first time an infinite term (which we've marked in blue). As long as the Hamiltonian dynamics :eq:`eq_free_field_hamilton_equation_q_and_p_dot` is concerned, it makes no difference adding a constant to the Hamiltonian. Hence we can write the free Hamiltonian for real scalar fields as follows

.. math::
	:label: eq_free_scalar_field_hamiltonian

	H_0^{\text{RSF}} = \frac{1}{2} \int d^3 x \left( p^2 + \left( \nabla q \right)^2 + m^2 q^2 \right)

.. warning::
	Throwing away the infinite term in :eq:`eq_calculate_free_real_scalar_field_hamiltonian` is an instance of a well-known criticism in quantum field theory -- just because something is infinite doesn't mean it's zero. Indeed, Weinberg mentioned in page 297 [Wei95]_ that such "infinities" shouldn't be thrown away when, for example, the fields are constrained within a finite space, or there is an involvement of gravity.

Now it's time to introduce the rather mysterious Lagrangian, which can be derived from the Hamiltonian via the so-called `Legendre transformation <https://en.wikipedia.org/wiki/Legendre_transformation>`__ as follows

.. math::
	:label: eq_legendre_transformation_lagrangian_from_hamiltonian

	L_0\left[ q(t), \dot{q}(t) \right] \coloneqq \sum_n \int d^3 x~p_n(t, \xbf) \dot{q}_n(t, \xbf) - H_0

where each occurrence of :math:`p_n(t)` is replaced by its expression in :math:`q_n(t)` and :math:`\dot{q}_n(t)`.

As a concrete example, let's consider again the real scalar field, where :math:`p = \dot{q}`. It follows that

.. math::
	:label: eq_free_real_scalar_field_lagrangian

	L_0^{\text{RSF}} &= \int d^3 x \left( p\dot{q} - \frac{1}{2} p^2 - \frac{1}{2} \left( \nabla q \right)^2 - \frac{1}{2} m^2 q^2 \right) \\
		&= \frac{1}{2} \int d^3 x \left( \dot{q}^2 - \left( \nabla q \right)^2 - m^2 q^2 \right) \\
		&= -\frac{1}{2} \int d^3 x \left( \p_{\mu} \psi \p^{\mu} \psi + m^2 \psi^2 \right)

It should be noted that expressing :math:`p` in terms of :math:`q` and :math:`\dot{q}` isn't always easy. Indeed, it's far from obvious how the :math:`p_i` defined by :eq:`eq_defn_q_and_p_vector_field_self_dual` could be expressed in the corresponding :math:`q_i` and :math:`\dot{q}_i`. (Un)Fortunately, we'd never really need to do so -- writing down a Lagrangian turns out to be mostly a guess work.


.. _sec_hamiltonian_and_lagrangian_for_interacting_fields:

Hamiltonian and Lagrangian for interacting fields
+++++++++++++++++++++++++++++++++++++++++++++++++

Let :math:`H` be the full Hamiltonian. Then the Heisenberg picture canonical variables can be defined as follows

.. math::
	:label: eq_defn_heisenberg_canonical_q_and_p

	Q_n(t, \xbf) &\coloneqq e^{\ifrak Ht} q_n(0, \xbf) e^{-\ifrak Ht} \\
	P_n(t, \xbf) &\coloneqq e^{\ifrak Ht} p_n(0, \xbf) e^{-\ifrak Ht}

Then obviously these canonical variables also satisfy the canonical (anti-)commutation relations

.. math::

	\left[ Q_n(t, \xbf), P_{n'}(t, \ybf) \right]_{\pm} &= \ifrak \delta^3(\xbf-\ybf) \delta_{n n'} \\
	\left[ Q_n(t, \xbf), Q_{n'}(t, \ybf) \right]_{\pm} &= 0 \\
	\left[ P_n(t, \xbf), P_{n'}(t, \ybf) \right]_{\pm} &= 0

Moreover, the analog of :eq:`eq_free_field_hamilton_equation_q_and_p_dot` holds as follows

.. math::

	\begin{alignat*}{2}
		\dot{Q}_n(t, \xbf) &= \ifrak \left[ H, Q_n(t, \xbf) \right] &&= \frac{\delta H}{\delta P_n(t, \xbf)} \\
		\dot{P}_n(t, \xbf) &= \ifrak \left[ H, P_n(t, \xbf) \right] &&= -\frac{\delta H}{ \delta Q_n(t, \xbf)}
	\end{alignat*}

As an example, we note that, in light of :eq:`eq_free_scalar_field_hamiltonian`, the full Hamiltonian for real scalar fields may be written as

.. math::

	H^{RSF} = \int d^3 x \left( \tfrac{1}{2} P^2 + \tfrac{1}{2} \left( \nabla Q \right)^2 + \tfrac{1}{2} m^2 Q^2 + \Hscr(Q) \right)

where :math:`\Hscr(Q)` is the perturbation term giving rise to the interaction.


The Lagrangian formalism
^^^^^^^^^^^^^^^^^^^^^^^^

We'll leave aside the discussion of canonical variables for a bit to introduce the Lagrangian formalism in its most general form. After that we'll play the game backwards. Namely, instead of constructing canonical variables out of the free fields that we've been exclusively considering since :ref:`sec_quantum_fields_and_antiparticles`, we'll get canonically conjugate fields out of the (magically appearing) Lagrangians, and then *impose* the canonical commutation relations :eq:`eq_canonical_commutation_relations` on them -- a procedure generally known as "quantization".

In the classical physical theory of fields, a Lagrangian is a functional :math:`L[\Psi(t), \dot{\Psi}(t)]`, where :math:`\Psi(t)` is any field and :math:`\dot{\Psi}(t)` is its time derivative. Here we've capitalized the field variables to distinguish them from the free fields considered in the previous section. Define the conjugate fields as follows

.. math:: \Pi_n(t, \xbf) \coloneqq \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)}
	:label: eq_general_lagrangian_conjugate_pi

so that the field equations are given by

.. math:: \dot{\Pi}_n(t, \xbf) = \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi_n(t, \xbf)}
	:label: eq_equation_of_motion_for_fields

.. warning::
	Unlike the functional derivatives considered in :eq:`eq_infinitesimal_variation_of_functional_of_canonical_variables` for canonical variables, the functional derivative :eq:`eq_general_lagrangian_conjugate_pi`, interpreted quantum mechanically, is not really well-defined since :math:`\Psi(t)` and :math:`\dot{\Psi}(t)` don't in general satisfy a simple (same time) commutation relation. According to Weinberg (see footnote on page 299 in [Wei95]_), "no important issues hinge on the details here". So we'll pretend that it behaves just like usual derivatives.

Indeed, recall that in the classical Lagrangian formalism, the field equations are given by a variational principle applied to the so-called *action*, defined as follows

.. math:: I[\Psi] \coloneqq \int_{-\infty}^{\infty} dt~L[\Psi(t), \dot{\Psi}(t)]
	:label: eq_defn_action_of_fields

The infinitesimal variation of :math:`I[\Psi]` is given by

.. math::

	\delta I[\Psi] &= \sum_n \int_{-\infty}^{\infty} dt \int d^3 x \left(
			\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi_n(t, \xbf)} \delta \Psi_n(t, \xbf) +
			\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)} \delta \dot{\Psi}_n(t, \xbf) \right) \\
		&= \sum_n \int_{-\infty}^{\infty} dt \int d^3 x \left(
			\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi(t, \xbf)} - \frac{d}{dt} \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)} \right) \delta \Psi_n(t, \xbf)

where for the last equality, integration by parts is used under the assumption that the infinitesimal variation :math:`\delta \Psi_n(t, \xbf)` vanishes at :math:`t \to \pm\infty`. Obviously :math:`\delta I[\Psi]` vanishes for any :math:`\delta \Psi_n(t, \xbf)` if and only if :eq:`eq_equation_of_motion_for_fields` is satisfied.

Now we're interested in constructing Lorentz invariant theories, but an action defined by :eq:`eq_defn_action_of_fields` apparently distinguishes the time from space variables. This motivates the hypothesis that the Lagrangian itself is given by a spatial integral of a so-called *Lagrangian density* as follows

.. math:: L[\Psi(t), \dot{\Psi}(t)] = \int d^3 x~\Lscr(\Psi(t, \xbf), \nabla\Psi(t, \xbf), \dot{\Psi}(t, \xbf))
	:label: eq_defn_lagrangian_density

In terms of the Lagrangian density, we can rewrite the action :eq:`eq_defn_action_of_fields` as a :math:`4`-integral as follows

.. math:: I[\Psi] = \int d^4 x~\Lscr(\Psi(x), \p_{\mu} \Psi(x))

We'd also like to reexpress the field equations :eq:`eq_equation_of_motion_for_fields` in terms of the Lagrangian density. To this end, let's first calculate the variation of :eq:`eq_defn_lagrangian_density` by an amount :math:`\delta \Psi_n(t, \xbf)` as follows

.. math::

	\delta L &= \sum_n \int d^3 x \left( \frac{\delta\Lscr}{\delta\Psi_n} \delta\Psi_n + \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \cdot \nabla \delta\Psi_n + \frac{\delta\Lscr}{\delta\dot{\Psi}_n} \delta\dot{\Psi}_n \right) \\
		&= \sum_n \int d^3 x \left( \left( \frac{\delta\Lscr}{\delta\Psi_n} - \nabla \cdot \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \right) \delta\Psi_n + \frac{\delta\Lscr}{\delta\dot{\Psi}_n} \delta\dot{\Psi}_n \right)

It follows that

.. math::

	\frac{\delta L}{\delta\Psi_n} &= \frac{\delta\Lscr}{\delta\Psi_n} - \nabla \cdot \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \\
	\frac{\delta L}{\delta\dot{\Psi}_n} &= \frac{\delta \Lscr}{\delta \dot{\Psi}_n}

Combining these with :eq:`eq_equation_of_motion_for_fields` and :eq:`eq_equation_of_motion_for_fields`, we've derived the so-called Euler-Lagrange equations for the Lagrangian density

.. math:: \frac{\delta \Lscr}{\delta \Psi_n} = \p_{\mu} \frac{\delta \Lscr}{\delta(\p_{\mu} \Psi_n)}
	:label: eq_euler_lagrange

Note that the summing :math:`4`-index :math:`\mu` here represents :math:`x_{\mu}`. Most importantly, the field equations given by :eq:`eq_euler_lagrange` will be Lorentz invariant if :math:`\Lscr` is. Indeed, guessing such :math:`\Lscr` will be more or less the only way to construct Lorentz invariant (quantum) field theories.

.. note::
	The Lagrangian density :math:`\Lscr` is assumed to be real for two reasons. First, if :math:`\Lscr` were complex, then splitting it into the real and imaginary parts, :eq:`eq_euler_lagrange` would contain twice as many equations as there are fields, regardless whether real or complex. This is undesirable because generically there will be no solutions. The second reason has to wait until the next section, where symmetries will be discussed. It turns out that the reality of :math:`\Lscr` will guarantee that the symmetry generators are Hermitian.

Now recall from the previous section that the anchor of our knowledge is the Hamiltonian -- we know how it must look like, at least for free fields. To go from the Lagrangian to the Hamiltonian, we use again the Legendre transformation (cf. :eq:`eq_legendre_transformation_lagrangian_from_hamiltonian`) to *define* the Hamiltonian as follows

.. math:: H[\Psi, \Pi] \coloneqq \sum_n \int d^3 x~\Pi_n(t, \xbf) \dot{\Psi}_n(t, \xbf) - L[\Psi(t), \dot{\Psi}(t)]
	:label: eq_legendre_transformation_hamiltonian_from_lagrangian

.. warning::
	In order to realize :math:`H` as a functional of :math:`\Psi` and :math:`\Pi`, one must in principle be able to solve for :math:`\dot{\Psi}_n` in terms of :math:`\Psi_n` and :math:`\Pi_n` from :eq:`eq_general_lagrangian_conjugate_pi`. This isn't always easy, if at all possible, but it rarely pose serious difficulties in applications either.

As a double check, let's verify that the Hamiltonian defined by :eq:`eq_legendre_transformation_hamiltonian_from_lagrangian` also satisfies Hamilton's equations (cf. :eq:`eq_free_field_hamilton_equation_q_and_p_dot`). Indeed, the variational derivatives are calculated, using :eq:`eq_general_lagrangian_conjugate_pi` and :eq:`eq_equation_of_motion_for_fields`, as follows

.. math::

	\frac{\delta H}{\delta \Pi_n(t, \xbf)}
		&= \sum_m \int d^3 y \left( \frac{\delta \Pi_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \dot{\Psi}_m(t, \ybf) + \Pi_m(t, \ybf) \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \right) \\
		&\qquad - \sum_m \int d^3 y \frac{\delta L}{\delta \dot{\Psi}_m(t, \ybf)} \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \\
		&= \sum_m \int d^3 y~\delta_{m,n} \delta^3(\ybf-\xbf) \dot{\Psi}_m(t, \ybf) \\
		&= \dot{\Psi}_n(t, \xbf) \\
	\frac{\delta H}{\delta \Psi_n(t, \xbf)} &= \sum_m \int d^3 y~\Pi_m(t, \ybf) \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} \\
		&\qquad - \sum_m \int d^3 y \left( \frac{\delta L}{\delta \Psi_m(t, \ybf)} \frac{\delta \Psi_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} + \frac{\delta L}{\delta \dot{\Psi}_m(t, \ybf)} \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} \right) \\
		&= -\sum_m \int d^3 y~\delta_{m, n} \delta^3(\ybf-\xbf) \dot{\Pi}_m(t, \ybf) \\
		&= -\dot{\Pi}_n(t, \xbf)

It's therefore attempting to demand, in the Lagrangian formalism, that :math:`\Psi_n` and :math:`\Pi_n`, defined by :eq:`eq_general_lagrangian_conjugate_pi`, satisfy the canonical commutation relations. In other words, they are (Heisenberg picture) canonically conjugate fields. But this is not true in general, as it turns out.

The issue is that the Lagrangian :math:`L[\Psi(t), \dot{\Psi}(t)]` may contain certain field, but not its time derivative. One example is spin-:math:`1` vector fields, where we see from :eq:`eq_defn_q_and_p_vector_field_self_dual` that the spatial fields :math:`\psi_i` are part of the canonical variables, but not :math:`\psi_0`, which nonetheless should present in the Lagrangian by Lorentz invariance. It turns out that what's missing from the Lagrangian is :math:`\dot{\psi}_0`, which causes its conjugate variable defined by :eq:`eq_general_lagrangian_conjugate_pi` to vanish.

But instead of dealing with vector fields further, we'll turn back to the general ground to establish the fundamental principles. Inspired by above discussion, we can rewrite the Lagrangian as

.. math:: L[Q(t), \dot{Q}(t), C(t)]
	:label: eq_general_quantum_lagrangian

where each :math:`Q_n(t)` has a corresponding :math:`\dot{Q}_n(t)`, but not for :math:`C(t)`. It follows that one can define the canonical conjugates by

.. math:: P_n(t, \xbf) \coloneqq \frac{\delta L[Q(t), \dot{Q}(t), C(t)]}{\delta \dot{Q}_n(t, \xbf)}

and hence the Hamiltonian takes the following form

.. math:: H[Q, P] = \sum_n \int d^3 x~P_n \dot{Q}_n - L[Q(t), \dot{Q}(t), C(t)]

.. dropdown:: Quantization of free scalar fields
	:animate: fade-in-slide-down
	:icon: unlock

	We'll illustrate how "quantization" works in the simplest case of free scalar fields :math:`\phi(t, \xbf)`. Namely, we'll reverse our earlier approach to the quantum theory by starting from a Lagrangian, and then work out the field equations, solve them for the fields, impose canonical commutation relations, and finally arrive at the familiar commutation relations between creation and annihilation operators introduced in :ref:`sec_the_cluster_decomposition_principle`.

	Following :eq:`eq_free_real_scalar_field_lagrangian`, let's consider the following Lagrangian

	.. math:: L_0[\phi, \dot{\phi}] = -\frac{1}{2} \int d^3 x \left( \p_{\mu} \phi \p^{\mu} \phi + m^2 \phi^2 \right)

	The canonical conjugate of :math:`\phi(t, \xbf)` is then

	.. math:: \pi(t, \xbf) \coloneqq \frac{\delta L_0}{\delta \dot{\phi}(t, \xbf)} = \dot{\phi}(t, \xbf)
		:label: eq_free_scalar_field_pi_equals_dot_phi

	Hence the Hamiltonian takes the following form

	.. math::

		H_0[\phi, \pi]
			&= \int d^3 x~\pi(t, \xbf) \dot{\phi}(t, \xbf) - L_0 \\
			&= \frac{1}{2} \int d^3 x \left( \pi^2(t, \xbf) + \big( \nabla \phi(t, \xbf) \big)^2 + m^2 \phi^2 \right)

	The field equations :eq:`eq_free_field_hamilton_equation_q_and_p_dot` are then given as follows

	.. math::

		\begin{alignat*}{2}
			\dot{\phi}(t, \xbf) &= \frac{\delta H_0}{\delta \pi(t, \xbf)} &&= \pi(t, \xbf) \\
			\dot{\pi}(t, \xbf) &= -\frac{\delta H_0}{\delta \phi(t, \xbf)} &&= \nabla^2 \phi(t, \xbf) - m^2 \phi(t, \xbf)
		\end{alignat*}

	Together, it implies that the field equation is precisely the Klein-Gordon equation

	.. math:: \left( \square - m^2 \right) \phi(x) = 0

	Using Fourier transform, the general Hermitian solution, up to a scalar, can be written as follows

	.. math::

		\phi(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} \left( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)

	where :math:`p_0 = \sqrt{\pbf^2+m^2}` and :math:`a(\pbf)` is, at the moment, just any operator function of :math:`\pbf`. Using :eq:`eq_free_scalar_field_pi_equals_dot_phi` we have

	.. math::

		\pi(x) = \frac{-\ifrak}{(2\pi)^{3/2}} \int d^3 p \sqrt{\frac{p_0}{2}} \left( e^{\ifrak p \cdot x}a(\pbf) - e^{-\ifrak p \cdot x}a^{\dagger}(\pbf) \right)

	One can then verify that if we impose the canonical commutation relations on the conjugate fields :math:`\phi(t, \xbf)` and :math:`\pi(t, \xbf)` as follows

	.. math::

		\left[ \phi(t, \xbf), \pi(y, \ybf) \right] &= \ifrak \delta^3(\xbf-\ybf) \\
		\left[ \phi(t, \xbf), \phi(t, \ybf) \right] &= 0 \\
		\left[ \pi(t, \xbf), \pi(t, \ybf) \right] &= 0

	then the familiar commutation relations

	.. math::

		\left[ a(\pbf), a^{\dagger}(\pbf') \right] &= \delta^3(\pbf-\pbf') \\
		\left[ a(\pbf), a(\pbf') \right] &= 0

	must hold. In this way we've completely reversed the process of deriving a Lagrangian from free fields made up of creation and annihilation operators.


Global symmetries
^^^^^^^^^^^^^^^^^

Of course, the reason for introducing the Lagrangian formalism is not to reproduce the Hamiltonians and the fields that we already knew. The main motivation is that, as we'll see, the Lagrangian formalism provides a framework for studying symmetries. Recall from :ref:`sec_what_is_a_symmetry` that a symmetry was defined to be a(n anti-)unitary transformation on the Hilbert space of states, i.e., a transformation that preserves amplitudes. Now in the Lagrangian formalism, field equations come out of the stationary action condition. Therefore in this context, we'll redefine a symmetry as an infinitesimal variation of the fields that leaves the action invariant. As it turns out, symmetries in this sense lead to conserved currents, which are nothing but the symmetry operators considered earlier. Hence besides a slight abuse of terminology, the notion of symmetries will be consistent.

.. note::
	Throughout this section, repeated indexes like :math:`n`, which are used to index various fields, in an equation are not automatically summed up. On the other hand, repeated :math:`4`-indexes like :math:`\mu` do follow the Einstein summation convention.

Consider an infinitesimal variation

.. math:: \Psi_n(x) \to \Psi_n(x) + \ifrak \epsilon \Fscr_n(x)
	:label: eq_infinitesimal_variation_of_field

which leaves the action :math:`I[\Psi]`  invariant

.. math:: 0 = \delta I = \ifrak \epsilon \sum_n \int d^4 x~\frac{\delta I[\Psi]}{\delta \Psi_n(x)} \Fscr_n(x)
	:label: eq_vanishing_of_action_under_infinitesimal_variation

A few remarks are in order. First of all, if we think of :eq:`eq_infinitesimal_variation_of_field` as an infinitesimal (unitary) symmetry transformation, then the coefficient :math:`\ifrak` can be justified by then intention of making :math:`\Fscr_n(x)` Hermitian. Next, although :eq:`eq_vanishing_of_action_under_infinitesimal_variation` *always* holds when :math:`\Psi_n(x)` is stationary, the infinitesimal :math:`\Fscr_n(x)` being a symmetry demands that :eq:`eq_vanishing_of_action_under_infinitesimal_variation` holds true for *any* :math:`\Psi_n(x)`. Finally, we emphasize the fact that :math:`\epsilon` is an infinitesimal *constant*, rather than a function of :math:`x`, is the defining property for the symmetry to be called "global". Indeed, we'll be dealing with symmetries that are not global in the next chapter, namely, the gauge symmetries.


.. _sec_from_symmetries_to_conservation_laws:

From symmetries to conservation laws
++++++++++++++++++++++++++++++++++++

The general principle that "symmetries imply conservation laws" is mathematically known as `Noether's theorem <https://en.wikipedia.org/wiki/Noether%27s_theorem>`__, but we'll not bother with any mathematical formality here. To see how to derive conserved quantitites from an assumed symmetry, let's change :eq:`eq_infinitesimal_variation_of_field` as follows

.. math:: \Psi_n(x) \to \Psi_n(x) + \ifrak \epsilon(x) \Fscr_n(x)
	:label: eq_functional_infinitesimal_variation_of_field

where :math:`\epsilon(x)` now is an infinitesimal function of :math:`x`. Under this variation, the corresponding :math:`\delta I` may not vanish. But it must take the following form

.. math:: \delta I = -\int d^4 x J^{\mu}(x) \p_{\mu} \epsilon(x)
	:label: eq_variation_of_action_by_functional_deformation

because it must vanishe when :math:`\epsilon(x)` is constant. Here :math:`J^{\mu}(x)` is a function(al) to be determined in individual cases, and is usually known as *current*. Now if :math:`\Psi_n(x)` satisfies the field equations, i.e., it's a stationary point of the action, then :eq:`eq_variation_of_action_by_functional_deformation` must vanishes for any :math:`\epsilon(x)`. Applying integration by parts (and assuming :math:`\Fscr_n(x)` vanishes at infinity), we must have

.. math:: \p_{\mu} J^{\mu}(x) = 0
	:label: eq_general_conservation_of_current

which is the conservation law for :math:`J`, which then can be called a conserved current. One gets also a conserved quantity, i.e., a quantity that doesn't change by time, by integrating :eq:`eq_general_conservation_of_current` over the :math:`3`-space as follows

.. math::

	\dot{J}^0(x) = -\nabla \cdot \Jbf(x)
		& \implies \int d^3 x~\dot{J}^0(x) = -\int d^3 x~\nabla \cdot \Jbf(x) = 0 \\
		& \implies F \coloneqq \int d^3 x~J^0(x) \text{ is conserved.}

Unfortunately, not much more can be said about the conserved current :math:`J` at this level of generality. This is, however, not the case if one imposes stronger assumptions on the symmetry, as we now explain.

Lagrangian-preserving symmetry
	This is the first strengthening of the symmetry assumption. Namely, instead of assuming that the variation :eq:`eq_infinitesimal_variation_of_field` fixes the action, we assume that it fixes the Lagrangian itself. Namely,

	.. math::
		:label: eq_stationary_lagrangian

		\delta L = \ifrak \epsilon \sum_n \int d^3 x \left( \frac{\delta L}{\delta \Psi_n(t, \xbf)} \Fscr_n(t, \xbf) + \frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \dot{\Fscr}_n(t, \xbf) \right) = 0

	Now let :math:`\epsilon(t)` be a time-dependent infinitesimal in :eq:`eq_functional_infinitesimal_variation_of_field`. Then we can calculate :math:`\delta I` under such variation as follows

	.. math::

		\delta I &= \ifrak \sum_n \int dt \int d^3 x \left( \frac{\delta L}{\delta \Psi_n(t, \xbf)} \epsilon(t) \Fscr_n(t, \xbf) + \frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \frac{d}{dt} \big( \epsilon(t) \Fscr_n(t, \xbf) \big) \right) \\
			&= \ifrak \sum_n \int dt \int d^3 x~\frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \dot{\epsilon}(t) \Fscr_n(t, \xbf)

	Comparing with :eq:`eq_variation_of_action_by_functional_deformation`, we can derive an explicit formula for the conserved quantity as follows

	.. math:: F = -\ifrak \sum_n \int d^3 x~\frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \Fscr_n(t, \xbf)
		:label: eq_lagrangian_preserving_symmetry_conserved_quantity

	Indeed, one can verify directly that :math:`\dot{F}(t) = 0` using :eq:`eq_stationary_lagrangian` together with the field equations :eq:`eq_general_lagrangian_conjugate_pi` and :eq:`eq_equation_of_motion_for_fields`.

Lagrangian-density-preserving symmetry
	Taking the previous assumption further, let's impose the even stronger condition that the Lagrangian density is invariant under :eq:`eq_infinitesimal_variation_of_field`. It means that

	.. math::
		:label: eq_stationary_lagrangian_density

		\delta \Lscr = \ifrak \epsilon \sum_n \left( \frac{\delta \Lscr}{\delta \Psi_n(x)} \Fscr_n(x) + \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \p_{\mu} \Fscr_n(x) \right) = 0

	Now under :eq:`eq_functional_infinitesimal_variation_of_field`, we can calculate the variation of the action as follows

	.. math::

		\delta I &= \ifrak \sum_n \int d^4 x~\left( \frac{\delta \Lscr}{\delta \Psi_n(x)} \epsilon(x) \Fscr_n(x) + \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \p_{\mu} \big( \epsilon(x) \Fscr_n(x) \big) \right) \\
			&= \ifrak \sum_n \int d^4 x~\frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \Fscr_n(x) \p_{\mu}\epsilon(x)

	Comparing with :eq:`eq_variation_of_action_by_functional_deformation` as before, we can derive an explicit formula for the conserved current as follows

	.. math:: J^{\mu}(x) = -\ifrak \sum_n \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \Fscr_n(x)
		:label: eq_lagrangian_density_preserving_symmetry_conserved_density

	Once again, one can directly verify that :math:`\p_{\mu} J^{\mu}(x) = 0` using :eq:`eq_stationary_lagrangian_density` together with the Euler-Lagrange equation :eq:`eq_euler_lagrange`.

So far everything has been completely classical. To make it a quantum theory, we'll involve the canonical fields introduced in :ref:`sec_hamiltonian_and_lagrangian_for_interacting_fields`. More precisely, instead of any :math:`\Fscr_n(t, \xbf)`, we'll suppose that it takes the following form

.. math:: \Fscr_n(Q(t), \xbf)

where :math:`Q(t)` is defined by :eq:`eq_defn_heisenberg_canonical_q_and_p`. Next, recall from :eq:`eq_general_quantum_lagrangian` that the field :math:`\Psi_n` is either a :math:`Q_n`, in which case :math:`\delta L / \delta \dot{Q}_n = P_n`, or a :math:`C_n`, in which case the functional derivative vaninshes.

Now in the case of a Lagrangian-preserving symmetry, the conserved quantity :eq:`eq_lagrangian_preserving_symmetry_conserved_quantity` takes the following form

.. math:: F = -\ifrak \sum_n \int d^3 x~P_n(t, \xbf) \Fscr_n(Q(t), \xbf)
	:label: eq_lagrangian_preserving_symmetry_generator_formula

which of course is time-independent. Moreover, one can show that :math:`F` in fact generates the quantum symmetry in the following sense

.. math::

	\left[ F, Q_n(t, \xbf) \right] = -\ifrak \sum_m \int d^3 y~\left[ P_m(t, \ybf), Q_n(t, \xbf) \right] \Fscr_m(Q(t), \ybf) = -\Fscr_n(Q(t), \xbf)

where we've taken advantage of the time-independency of :math:`F` to arrange the same-time commutator.


.. _sec_spacetime_translations:

Spacetime translations
++++++++++++++++++++++

So far the symmetries have been rather abstract, to make it more explicit, and also to get warmed up for the general case, let's assume the Lagrangian is invariant under the (spacetime) translation transformation given as follows

.. math:: \Psi_n(x) \to \Psi_n(x + \epsilon) = \Psi_n(x) + \epsilon^{\mu} \p_{\mu} \Psi_n(x)

Comparing with :eq:`eq_infinitesimal_variation_of_field` we see that

.. math:: \Fscr_{\mu} = -\ifrak \p_{\mu} \Psi_n

It follows from :eq:`eq_variation_of_action_by_functional_deformation` and :eq:`eq_general_conservation_of_current` that there exists a conserved :math:`4`-current :math:`{T^{\nu}}_{\mu}`, which is known as the `energy-momentum tensor <https://en.wikipedia.org/wiki/Stress%E2%80%93energy_tensor>`__, such that

.. math:: \p_{\nu} {T^{\nu}}_{\mu} = 0

The corresponding conserved currents then take the form

.. math:: P_{\mu} \coloneqq \int d^3 x~{T^0}_{\mu}
	:label: eq_spacetime_translation_conserved_quantity_is_momentum

such that :math:`\dot{P}_{\mu} = 0`. Here it's important to not confuse :math:`P_{\mu}` with a canonical variable -- it's just a conserved quantity which turns out to be the :math:`4`-momentum.

Now recall from :eq:`eq_defn_lagrangian_density` that the Lagrangian is usually the spatial integral of a density functional. Hence it's not unreasonable to suppose that the Lagrangian is indeed invariant under spatial translations. Under this assumption, we can rewrite :eq:`eq_lagrangian_preserving_symmetry_generator_formula` as follows

.. math:: \Pbf \coloneqq -\sum_n \int d^3 x~P_n(t, \xbf) \nabla Q_n(t, \xbf)
	:label: eq_spatial_translation_conserved_quantity

with the understanding that :math:`\Psi_n = Q_n`.

To verify that :math:`\Pbf` indeed generates spatial translations, let's calculate using the fact that :math:`\Pbf` is time-independent as follows

.. math::

	\left[ \Pbf, Q_n(t, \xbf) \right] &= \ifrak \nabla Q_n(t, \xbf) \\
	\left[ \Pbf, P_n(t, \xbf) \right] &= \ifrak \nabla P_n(t, \xbf)

It follows that

.. math:: \left[ \Pbf, \Gscr \right] = \ifrak \nabla \Gscr
	:label: eq_momenta_act_as_spatial_derivative

for any functional :math:`\Gscr` that doesn't explicitly involve :math:`\xbf`. This verifies that :math:`\Pbf` indeed generates the spatial translation.

In contrast, one cannot hope that the Lagrangian to be invariant under time translation, if there should be any interaction. But we already know the operator that generates time translation, namely, the Hamiltonian. In other words, we define :math:`P_0 \coloneqq -H` such that

.. math:: \left[ H, \Gscr \right] = -\ifrak \dot{\Gscr}
	:label: eq_hamiltonian_acts_as_time_derivative

for any functional :math:`\Gscr` that doesn't explicitly involve :math:`t`.

In general, the Lagrangian density is not invariant under spacetime translations. However, it turns out that the conserved current, which in this case is :math:`{T^{\mu}}_{\nu}`, can nonetheless be calculated. To spell out the details, let's consider the following variation

.. math:: \Psi_n(x) \to \Psi_n(x + \epsilon(x)) = \Psi_n(x) + \epsilon^{\mu}(x) \p_{\mu} \Psi_n(x)

The corresponding variation of the action is given as follows

.. math::

	\delta I[\Psi] &= \sum_n \int d^4 x \left( \frac{\delta \Lscr}{\delta \Psi_n} \epsilon^{\mu} \p_{\mu} \Psi_n + \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\nu}(\epsilon^{\mu} \p_{\mu} \Psi_n) \right) \\
		&= \int d^4 x \left( \epsilon^{\mu} \p_{\mu} \Lscr + \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu}\Psi_n \p_{\nu} \epsilon^{\mu} \right) \\
		&= -\int d^4 x \left( \delta^{\nu}_{\mu} \Lscr - \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu} \Psi_n \right) \p_{\nu} \epsilon^{\mu}

where we've used the chain rule for derivatives in the second equality, and integration by parts in the third. Comparing with :eq:`eq_variation_of_action_by_functional_deformation`, we see that

.. math::
	:label: eq_energy_momentum_tensor_from_translation_invariance

	{T^{\nu}}_{\mu} = \delta^{\nu}_{\mu} \Lscr - \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu} \Psi_n


.. _note_energy_momentum_tensor_not_symmetric:

.. note::
	The energy-momentum tensor :math:`{T^{\nu}}_{\mu}` is not yet suitable for general relativity since it's not symmetric. As we'll see in :ref:`sec_lorentz_symmetry`, when taking homogeneous Lorentz transformation symmetry into account, one can supplement :math:`{T^{\nu}}_{\mu}` with some extra terms to make it both conserved and symmetric.

Indeed, this calculation recovers :eq:`eq_spatial_translation_conserved_quantity` by letting :math:`\nu = 0` and :math:`\mu \neq 0`. Moreover, it recovers the Hamiltonian by letting :math:`\mu = \nu = 0` as follows

.. math:: H = -P_0 = \int d^3 x \left( \sum_n P_n \dot{Q}_n - \Lscr \right)


Linear transformations
++++++++++++++++++++++

As another example, let's consider linear variations as follows

.. math::

	Q_n(x) &\to Q_n(x) + \ifrak \epsilon^a {(t_a)_n}^m Q_m(x) \\
	C_r(x) &\to C_r(x) + \ifrak \epsilon^a {(\tau_a)_r}^s C_s(x)

where we've adopted the Einstein summation convention for repeated upper and lower indexes because it'd otherwise be too tedious to write out the summations. Here :math:`(t_{\square})^{\square}_{\square}` should furnish a representation of the Lie algebra of the symmetry group.

As before, the invariance of action under such variations implies the existstence of conserved currents :math:`J^{\mu}_a` such that

.. math:: \p_{\mu} J^{\mu}_a = 0

as well as the conserved quantity

.. math:: T_a \coloneqq \int d^3 x~J^0_a

If, in addition, the Lagrangian is invariant under such variations, then :math:`T_a` takes the following form by :eq:`eq_lagrangian_preserving_symmetry_generator_formula`

.. math:: T_a = -\ifrak \int d^3 x~P_n(t, \xbf) {(t_a)^n}_m Q^m(t, \xbf)

It follows that

.. math::

	\left[ T_a, Q^n(x) \right] &= -{(t_a)^n}_m Q^m(x) \\
	\left[ T_a, P_n(x) \right] &= {(t_a)_n}^m P_m(x)

In particular, when :math:`t_a` is diagonal (e.g., in electrodynamics), the operators :math:`Q^n` and :math:`P_n` may be regarded as raising/lowering operators. In fact, we claim that :math:`T_a` form a Lie algebra by the following calculation

.. math::

	\left[ T_a, T_b \right] &= -\left[ \int d^3 x~P_n(t, \xbf) {(t_a)^n}_m Q^m(t, \xbf), \int d^3 y~P_r(t, \ybf) {(t_b)^r}_s Q^s(t, \ybf) \right] \\
		&= -\int d^3 x~d^3 y~{(t_a)^n}_m {(t_b)^r}_s \left[ P_n(t, \xbf) Q^m(t, \xbf), P_r(t, \ybf) Q^s(t, \ybf) \right] \\
		&= -\int d^3 x~d^3 y~{(t_a)^n}_m {(t_b)^r}_s \Big( P_n(t, \xbf) \left[ Q^m(t, \xbf), P_r(t, \ybf) \right] Q^s(t, \ybf) \\
			&\qquad - P_r(t, \ybf) \left[ Q^s(t, \ybf), P_n(t, \xbf) \right] Q^m(t, \xbf) \Big) \\
		&= -\ifrak \int d^3 x \Big( {(t_a)^n}_m {(t_b)^m}_s P_n(t, \xbf) Q^s(t, \xbf) - {(t_a)^n}_m {(t_b)^r}_n P_r(t, \xbf) Q^m(t, \xbf) \Big) \\
		&= -\ifrak \int d^3 x~{\left[ t_a, t_b \right]^n}_m P_n(t, \xbf) Q^m(t, \xbf)

Now if :math:`t_a` form a Lie algebra with structure constants :math:`{f_{ab}}^c` as follows

.. math:: \left[ t_a, t_b \right] = \ifrak {f_{ab}}^c t_c

then

.. math:: \left[ T_a, T_b \right] = \ifrak {f_{ab}}^c T_c

In other words, the conserved quantities also form the same Lie algebra.

Now if, in addition, the Lagrangian density is also invariant, then :eq:`eq_lagrangian_density_preserving_symmetry_conserved_density` takes the following form

.. math::
	:label: eq_lagrangian_density_invariant_linear_transformation_conserved_current

	J^{\mu}_a = -\ifrak \frac{\delta \Lscr}{\delta (\p_{\mu} Q_n)} {(t_a)_n}^m Q_m - \ifrak \frac{\delta \Lscr}{\delta (\p_{\mu} C_r)} {(\tau_a)_r}^s C_s

.. dropdown:: Interacting equal-mass real scalar fields
	:animate: fade-in-slide-down
	:icon: unlock

	As a specific example, let's consider the following Lagrangian density for two interacting equal-mass real scalar fields

	.. math::

		\Lscr = -\frac{1}{2} \p_{\mu} \Phi_1 \p^{\mu} \Phi_1 - \frac{1}{2} m^2 \Phi_1^2 - \frac{1}{2} \p_{\mu} \Phi_2 \p^{\mu} \Phi_2 - \frac{1}{2} m^2 \Phi_2^2 - \Hscr(\Phi_1^2 + \Phi_2^2)

	This density is invariant under the following linear transformation

	.. math::

		\Phi_1 &\to \Phi_1 - \epsilon \Phi_2 \\
		\Phi_2 &\to \Phi_2 + \epsilon \Phi_1

	In this case, we can evaluate :eq:`eq_lagrangian_density_invariant_linear_transformation_conserved_current` as follows

	.. math:: J^{\mu} = -\Psi_2 \p^{\mu} \Phi_1 + \Phi_1 \p^{\mu} \Phi_2

Note that since :math:`\Lscr` doesn't have :math:`\dot{C}_r` dependencies, we have the following by letting :math:`\mu = 0` in :eq:`eq_lagrangian_density_invariant_linear_transformation_conserved_current`

.. math:: J^0_a = -\ifrak P^n {(t_a)_n}^m Q_m

whose equal-time commutation relations with canonical variables :math:`P` and :math:`Q` can be easily calculated.


.. _sec_lorentz_invariance:

Lorentz invariance
^^^^^^^^^^^^^^^^^^

The goal of this section is to show that the Lorentz invariance of the Lagrangian density implies the Lorentz invariance of the S-matrix, which justifies our interest in the Lagrangian formalism in the first place.

Recall from :eq:`eq_expansion_of_Lambda` and :eq:`eq_lorentz_lie_algebra_is_antisymmetric` that

.. math::
	:label: eq_lorentz_omega_is_antisymmetric

	{\Lambda_{\mu}}^{\nu} &= {\delta_{\mu}}^{\nu} + {\omega_{\mu}}^{\nu} \\
	\omega_{\mu \nu} &= -\omega_{\nu \mu}

is a :math:`(\mu, \nu)`-parametrized anti-symmetric variation. It follows then from :eq:`eq_variation_of_action_by_functional_deformation` that there exist :math:`(\mu, \nu)`-parametrized anti-symmetric conserved currents as follows

.. math::
	:label: eq_lorentz_invariance_m_conservation_and_antisymmetry

	\p_{\rho} \Mscr^{\rho \mu \nu} &= 0 \\
	\Mscr^{\rho \mu \nu} &= -\Mscr^{\rho \nu \mu}

which, in turn, make conversed quantities

.. math:: J^{\mu \nu} \coloneqq \int d^3 x~\Mscr^{0 \mu \nu}
	:label: eq_lorentz_invariance_conserved_j

such that :math:`\dot{J}^{\mu \nu} = 0`. These, as we'll see, turn out to be rather familiar objects that we've encountered as early as in :eq:`eq_u_lorentz_expansion`.

In light of :eq:`eq_lagrangian_density_preserving_symmetry_conserved_density`, one can work out an explicit formula for :math:`\Mscr^{\rho \mu \nu}` if the Lagrangian density is invariant under the symmetry transformation. Now since the Lagrangian density is expressed in terms of quantum fields, one'd like to know how they transform under Lorentz transformations. Since the translation symmetry has already been dealt with in :ref:`sec_spacetime_translations`, we'll consider here homogeneous Lorentz transformations. Luckily this has been worked out already in :ref:`sec_quantum_fields_and_antiparticles`. More precisely, recall from :eq:`eq_dirac_field_linearize_representation` that the variation term can be written as follows

.. math:: \delta \Psi_n = \frac{\ifrak}{2} \omega^{\mu \nu} {(\Jscr_{\mu \nu})_n}^m \Psi_m

where :math:`\Jscr` are matrices satisfying :eq:`eq_bracket_repr_j`. The corresponding derivatives then have the following variation term

.. math::

	\delta (\p_{\kappa} \Psi_n) = \frac{\ifrak}{2} \omega^{\mu \nu} {(\Jscr_{\mu \nu})_n}^m \p_{\kappa} \Psi_m + {\omega_{\kappa}}^{\lambda} \p_{\lambda} \Psi_n

where the second summand on the right-hand-side corresponds to the fact the the Lorentz transformation also acts on the spacetime coordinates.

Now the invariance of the Lagrangian density under such variation can be written as follows

.. math::

	\frac{\delta \Lscr}{\delta \Psi_n} \frac{\ifrak}{2} \omega^{\mu \nu} {(\Jscr_{\mu \nu})_n}^m \Psi_m
		+ \frac{\delta \Lscr}{\delta (\p_{\kappa} \Psi_n)} \left( \frac{\ifrak}{2} \omega^{\mu \nu} ({\Jscr_{\mu \nu})_n}^m \p_{\kappa} \Psi_m + {\omega_{\kappa}}^{\lambda} \p_{\lambda} \Psi_n \right) = 0

Since :math:`\omega^{\mu \nu}` is not in general zero, its coefficient must be zero, which, taking :eq:`eq_lorentz_omega_is_antisymmetric` into account, implies the following

.. math::
	:label: eq_lorentz_invariance_current_raw_identity

	& \frac{\ifrak}{2} \frac{\delta\Lscr}{\delta\Psi_n} {(\Jscr_{\mu\nu})_n}^m \Psi_m
		+ \frac{\ifrak}{2} \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} ({\Jscr_{\mu\nu})_n}^m \p_{\kappa}\Psi_m \\
		& \qquad + \frac{1}{2} \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} \left( \eta_{\kappa \mu} \p_{\nu} - \eta_{\kappa \nu} \p_{\mu} \right) \Psi_n = 0

Using :eq:`eq_euler_lagrange`, we can get rid of the :math:`\delta\Lscr / \delta\Psi_n` term in :eq:`eq_lorentz_invariance_current_raw_identity` to arrive at the following

.. math::
	:label: eq_lorentz_invariance_current_identity

	\ifrak \p_{\kappa} \left( \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} {(\Jscr_{\mu\nu})_n}^m \Psi_m \right) - T_{\mu\nu} + T_{\nu\mu} = 0

where we've also used :eq:`eq_energy_momentum_tensor_from_translation_invariance`. Now we can address the issue of :ref:`energy-momentum tensor not being symmetric <note_energy_momentum_tensor_not_symmetric>` by introducing the following so-called `Belinfante tensor <https://en.wikipedia.org/wiki/Belinfante%E2%80%93Rosenfeld_stress%E2%80%93energy_tensor>`__

.. math::
	:label: eq_defn_belinfante_tensor

	\Theta_{\mu\nu}
		&\coloneqq T_{\mu\nu} - \frac{\ifrak}{2} \p_{\kappa} \Big(
			\frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} {(\Jscr_{\mu\nu})_n}^m \Psi_m -
			\frac{\delta\Lscr}{\delta(\p_{\mu} \Psi_n)} {(\Jscr_{\kappa\nu})_n}^m \Psi_m \\
		&\qquad - \frac{\delta\Lscr}{\delta(\p_{\nu} \Psi_n)} {(\Jscr_{\kappa\mu})_n}^m \Psi_m \Big)

which is both conserved in the sense that

.. math:: \p^{\mu} \Theta_{\mu\nu} = 0
	:label: eq_belinfante_tensor_is_conserved

and symmetric in the sense that

.. math:: \Theta_{\mu\nu} = \Theta_{\nu\mu}
	:label: eq_belinfante_tensor_is_symmetric

Indeed :eq:`eq_belinfante_tensor_is_conserved` follows from the observation that the term inside the parenthesis of :eq:`eq_defn_belinfante_tensor` is anti-symmetric in :math:`\mu` and :math:`\kappa`, and :eq:`eq_belinfante_tensor_is_symmetric` is a direct consequence of :eq:`eq_lorentz_invariance_current_identity`.

The conserved quantities corresponding to :math:`\Theta_{\mu\nu}`, according to :eq:`eq_spacetime_translation_conserved_quantity_is_momentum` are

.. math:: \int d^3 x~{\Theta^0}_\nu = \int d^3 x~{T^0}_\nu = P_{\nu}
	:label: eq_p_as_integral_of_belinfante_tensor

where the first equality holds because, again, the item in the parenthesis of :eq:`eq_defn_belinfante_tensor` is anti-symmetric is :math:`\mu` and :math:`\kappa`, and therefore :math:`\kappa \neq 0` given :math:`\mu = 0`. Hence it's at least equally legitimate to call :math:`\Theta_{\mu \nu}` the energy-momentum tensor. Indeed, the fact that :math:`\Theta_{\mu \nu}` is the symmetric makes it suitable for general relativity.

Unlike the other conserved currents, which are derived under the general principles explained in :ref:`sec_from_symmetries_to_conservation_laws`, we'll construct the anti-symmetric :math:`\Mscr^{\rho \mu \nu}` declared in :eq:`eq_lorentz_invariance_m_conservation_and_antisymmetry` by hand as follows

.. math:: \Mscr^{\rho\mu\nu} \coloneqq x^{\mu} \Theta^{\rho\nu} - x^{\nu} \Theta^{\rho\mu}

While :eq:`eq_lorentz_omega_is_antisymmetric` is automatically satisfied by definition, we can verify :eq:`eq_lorentz_invariance_m_conservation_and_antisymmetry`, using :eq:`eq_belinfante_tensor_is_conserved` and :eq:`eq_belinfante_tensor_is_symmetric`, as follows

.. math:: \p_{\rho} \Mscr^{\rho\mu\nu} = \Theta^{\mu\nu} - \Theta^{\nu\mu} = 0

Moreover :eq:`eq_lorentz_invariance_conserved_j` takes the following form

.. math:: J^{\mu\nu} = \int d^3 x \left( x^{\mu} \Theta^{0\nu} - x^{\nu} \Theta^{0\mu} \right)
	:label: eq_j_by_belinfante_tensor

Now if we consider the rotation generators defined by

.. math:: J_i \coloneqq \tfrac{1}{2} \epsilon_{ijk} J^{jk}

then it follows from :eq:`eq_hamiltonian_acts_as_time_derivative` that

.. math:: [H, \Jbf] = -\ifrak \dot{\Jbf} = 0

since :math:`\Jbf` doesn't implicitly involve :math:`t`. This recovers one of the commutation relations :eq:`eq_poincare_algebra` for the Poincaré algebra. Next, let's verify the commutation relation between :math:`\Pbf` and :math:`\Jbf`, using :eq:`eq_momenta_act_as_spatial_derivative` and :eq:`eq_p_as_integral_of_belinfante_tensor`, as follows

.. math::

	[P_i, J_j] &= \frac{1}{2} \epsilon_{jk\ell} \left[ P_i, J^{k\ell} \right] \\
		&= \frac{\ifrak}{2} \epsilon_{jk\ell} \int d^3x \left( x^k \p_i \Theta^{0\ell} - x^{\ell} \p_i \Theta^{0k} \right) \\
		&= \frac{\ifrak}{2} \epsilon_{jk\ell} \int d^3x \left( -\delta^k_i \Theta^{0\ell} + \delta^{\ell}_i \Theta^{0k} \right) \\
		&= \ifrak \epsilon_{ijk} \int d^3x~\Theta^{0k} \\
		&= \ifrak \epsilon_{ijk} P^k

What come next are the boost operators defined as follows [#k_convention]_

.. math:: K^i \coloneqq J^{0i} = \int d^3x \left( x^0 \Theta^{0i} - x^i \Theta^{00} \right)
	:label: eq_boost_k_by_belinfante_tensor

Bringing down the index, we can rewrite :eq:`eq_boost_k_by_belinfante_tensor` in vector form as follows

.. math:: \Kbf = t \Pbf - \int d^3 x~\xbf \Theta^{00}
	:label: eq_boost_k_by_belinfante_tensor_vector_form

Now it follows from :eq:`eq_hamiltonian_acts_as_time_derivative` that

.. math::

	[H, \Kbf] &= t[H, \Pbf] + \ifrak \int d^3 x~\xbf \dot{\Theta}^{00} \\
		&= \ifrak \int d^3 x~\xbf \dot{\Theta}^{00} \\
		&= \ifrak (\Pbf - \dot{\Kbf}) = \ifrak \Pbf

which is consistent with :eq:`eq_poincare_algebra`.

Finally, using :eq:`eq_momenta_act_as_spatial_derivative` and :eq:`eq_boost_k_by_belinfante_tensor_vector_form` together with the fact that :math:`\Pbf` commutes with itself, one can evaluate the commutator between :math:`\Pbf` and :math:`\Kbf` as follows

.. math::

	\left[ P_j, K_k \right] = -\ifrak \int d^3 x~x_k \p_j \Theta^{00} = \ifrak \delta_{kj} \int d^3 x \Theta^{00} = \ifrak \delta_{kj} P^0 = \ifrak \delta_{kj} H

which, again, is consistent with :eq:`eq_poincare_algebra`.

It turns out, following the lines of argument in :ref:`Lorentz symmetry of S-matrix <sec_s_matrix_lorentz_symmetry>`, these commutation relations are enough to show the Lorentz invariance of S-matrix under the same "smootheness" assumptions on the interaction terms. In addition, the other commutation relations between :math:`H, \Pbf, \Jbf, \Kbf` also follows.

Though not necessary, it's indeed possible to verify the other Poincaré algebra relations directly. In particular, the commutation relations between the rotation generators are verified as follows.

.. dropdown:: An explicit formula for rotation generators :math:`J^{ij}`
	:animate: fade-in-slide-down
	:icon: unlock
	:open:

	According to :eq:`eq_j_by_belinfante_tensor`, :eq:`eq_defn_belinfante_tensor`, and :eq:`eq_energy_momentum_tensor_from_translation_invariance`, the rotation generator :math:`J^{ij}` can be calculated as follows

	.. math::

		J^{ij} &= \int d^3 x \left( x^i \Theta^{0j} - x^j \Theta^{0i} \right) \\
			&= \int d^3 x \left( x^i T^{0j} - x^j T^{0i} \right) \\
				&\mkern-24mu - \frac{\ifrak}{2} \int d^3 x~x^i \p_k \left(
					\frac{\delta \Lscr}{\delta (\p_k \Psi_n)} {\left( \Jscr^{0j} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta \dot{\Psi}_n} {\left( \Jscr^{kj} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta (\p_j \Psi_n)} {\left( \Jscr^{k0} \right)_n}^m \Psi_m
				\right) \\
				&\mkern-24mu + \frac{\ifrak}{2} \int d^3 x~x^j \p_k \left(
					\frac{\delta \Lscr}{\delta (\p_k \Psi_n)} {\left( \Jscr^{0i} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta \dot{\Psi}_n} {\left( \Jscr^{ki} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta (\p_i \Psi_n)} {\left( \Jscr^{k0} \right)_n}^m \Psi_m
				\right) \\
			&= \int d^3 x \left( x^i T^{0j} - x^j T^{0i} \right) \\
				&\mkern-24mu + \frac{\ifrak}{2} \int d^3 x \left(
					\frac{\delta \Lscr}{\delta (\p_i \Psi_n)} {\left( \Jscr^{0j} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta \dot{\Psi}_n} {\left( \Jscr^{ij} \right)_n}^m \Psi_m
					- \frac{\delta \Lscr}{\delta (\p_j \Psi_n)} {\left( \Jscr^{i0} \right)_n}^m \Psi_m
				\right) \\
				&\mkern-24mu - \frac{\ifrak}{2} \int d^3 x \left(
						\frac{\delta \Lscr}{\delta (\p_j \Psi_n)} {\left( \Jscr^{0i} \right)_n}^m \Psi_m
						- \frac{\delta \Lscr}{\delta \dot{\Psi}_n} {\left( \Jscr^{ji} \right)_n}^m \Psi_m
						- \frac{\delta \Lscr}{\delta (\p_i \Psi_n)} {\left( \Jscr^{j0} \right)_n}^m \Psi_m
					\right) \\
			&= \int d^3 x \left( x^i T^{0j} - x^j T^{0i} \right) - \ifrak \int d^3 x \frac{\delta \Lscr}{\delta \dot{\Psi}_n} {\left( \Jscr^{ij} \right)_n}^m \Psi_m \\
			&= \int d^3 x \frac{\delta \Lscr}{\delta \dot{\Psi}_n} \left( -x^i \p^j \Psi_n + x^j \p^i \Psi_n - \ifrak {\left( \Jscr^{ij} \right)_n}^m \Psi_m \right)


.. rubric:: Footnotes

.. [#k_convention] The definition :math:`K_k \coloneqq J^{k0}` on page 317 in [Wei95]_ differs from the very same definition on page 61 by a sign, which leads to, among others, inconsistency in the Poincaré algebra relations :eq:`eq_poincare_algebra`. We will stick to our convention in :ref:`sec_quantum_lorentz_symmetry`, which is consistent with Weinberg's convention on page 61.
