.. _sec_scattering_theory:

Scattering Theory
=================

Physics would have been rather boring if nothing interacts, like the free particles that we have been studying so far. On the flip side, physics would have been impossible if we try to know exactly what happens in the interactions. The middle ground, where we assume that the particles are non-interacting long before and after the interaction, and something mysterious happened in between, is called scattering theory -- a place where theories meet experiments.

Non-Interacting Many-Particles State
------------------------------------

We shall, as always, start from the easiest part of the theory, which is clearly the non-interacting parts. Recall our grand formulae for the Lorentz transformation on one-particle state :eq:`eq_translation_formula_for_particle_state` and :eq:`eq_lorentz_transformation_formula_for_particle_state`. For a non-interacting many-particles system, it's conceivable to assume that the Lorentz transformation law is simply a direct product of the individual particles. Since :math:`U(\Lambda, a) = U(1, a) U(\Lambda, 0)`, we have the following

.. math::
	:label: eq_lorentz_transformation_formula_for_many_free_particles

	U(\Lambda, a) \Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} =
		&~ \exp(-\ifrak a^{\mu} ((\Lambda p_1)_{\mu} + (\Lambda p_2)_{\mu} + \cdots)) \\
		&\times \sqrt{\frac{(\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p_1)_0 (p_2)_0 \cdots}} \\
		&\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma'_1 \sigma_1}(W_1(\Lambda, p_1)) D_{\sigma'_2 \sigma_2}(W_2(\Lambda, p_2)) \cdots \\
		&\times \Psi_{\Lambda p_1, \sigma'_1, n_1; ~\Lambda p_2, \sigma'_2, n_2; ~\cdots}

where the first component is the translation transformation :eq:`eq_translation_formula_for_particle_state`, the second component is the normalization factor, and the third component is the little group representation (cf. :eq:`eq_d_repr_of_little_group`), and the :math:`\sigma`'s are either the spin :math:`z`-component for massive particles or the helicity for massless particles, and the :math:`n`'s are additional (discrete) labels such as mass, charge, spin, etc.

Notice that by writing a many-particles state as :math:`\Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}`, we have given the particles an order, which is rather random. Hence the normalization of these states must take permutations of the particles into account. More precisely, we have

.. math::
	:label: eq_many_particles_state_normalization_rough

	\left( \Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}, \Psi_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots} \right)
		&= \delta^3(\pbf_1 - \pbf'_1) \delta_{\sigma_1 \sigma'_1} \delta_{n_1 n'_1} \delta^3(\pbf_2 - \pbf'_2) \delta_{\sigma_2 \sigma'_2} \delta_{n_2 n'_2} \\
		&\quad \pm \text{permutations}

The sign in front of the permutations has to do with the species of the particles, which will be discussed later. Note that although there are many terms in :eq:`eq_many_particles_state_normalization_rough`, there is at most one nonzero term, which happens exactly when the two states differ by a permutation.

To suppress the annoyingly many sub-indexes in the states, we shall use letters such as :math:`\alpha, \beta, \cdots` to denote the compound index such as :math:`(p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots)`, so that, for example, :eq:`eq_many_particles_state_normalization_rough` can be simplified as

.. math:: \left( \Psi_{\alpha}, \Psi_{\alpha'} \right) = \delta(\alpha - \alpha')

where the integral volume element reads

.. math:: \int d\alpha \cdots = \sum_{\sigma_1, n_1; ~\sigma_2, n_2; ~\cdots} \int d^3 \pbf_1 d^3 \pbf_2 \cdots

We have postulated that the transformation law :eq:`eq_lorentz_transformation_formula_for_many_free_particles` works for non-interacting particles, but in fact, it's also only possible for non-interacting particles. One way to see this is through an energy calculation by letting :math:`\Lambda = 1` and :math:`a = (\tau, 0, 0, 0)` in :eq:`eq_lorentz_transformation_formula_for_many_free_particles` to see that

.. math::

	\exp(\ifrak \tau E_{\alpha}) \Psi_{\alpha} = \exp(\ifrak \tau H) \Psi_{\alpha} = \exp(\ifrak \tau (E_1 + E_2 + \cdots)) \Psi_{\alpha} \implies E_\alpha = E_1 + E_2 + \cdots

where :math:`E_i \coloneqq (p_i)_0` is the energy of the :math:`i`-th particle. There is obviously no energy left for any interaction.

In- and Out-states
------------------

As mentioned earlier, scattering theory is concerned with a scenario where interactions happen within a finite time period, long before and after which the system can be regarded as non-interacting. We can therefore define the in-state :math:`\Psi_{\alpha}^-` and the out-state :math:`\Psi_{\alpha}^+`, where :math:`\alpha` is the compound index as defined in the previous section, such that the states appear to be non-interacting with the prescribed particle states when *observed* at :math:`t \to \mp \infty`, respectively. [#in_out_state_sign_convention]_

Now it's time to bring forward an implicit assumption on the quantum states that we've been studying so far: they're defined in one chosen inertial frame. Indeed, the Lorentz transformation law :eq:`eq_lorentz_transformation_formula_for_many_free_particles` tells us exactly how to transform the state to any other frame. States of this sort are called `Heisenberg picture <https://en.wikipedia.org/wiki/Heisenberg_picture>`_ states: they contain the entire history/future of the system and are not dynamical in time as opposed to the so-called `Schrödinger picture <https://en.wikipedia.org/wiki/Schr%C3%B6dinger_picture>`_ states.

Back to the scattering scenario, let's imagine a reference observer :math:`\Ocal`, who at :math:`t = 0` observes that the system is in a state :math:`\Psi`. Then imagine another observer :math:`\Ocal'` at rest with respect to :math:`\Ocal`, who sets his clock :math:`t' = 0` when :math:`t = \tau`, in other words :math:`t' = t - \tau`. Then from the viewpoint of :math:`\Ocal'`, the time-:math:`0` state should look like :math:`\exp(-\ifrak \tau H) \Psi`. It follows that the state :math:`\Psi`, viewed long before and long after the reference :math:`t = 0`, should look like :math:`\exp(-\ifrak \tau H) \Psi` for :math:`\tau \to \mp\infty`, respectively.

It follows that energy eigenstates such as :math:`\Psi_{\alpha}` will look the same at all time since

.. math:: \exp(-\ifrak \tau H) \Psi_{\alpha} = \exp(-\ifrak \tau E_{\alpha}) \Psi_{\alpha}

creates merely an inconsequential phase factor. This is one form of the uncertainty principle: if the energy is definitely known, then the time is completely unknown. Therefore we must consider a localized packet (or superposition) of states as follows

.. math::
	:label: eq_psi_packet

	\int d\alpha ~g(\alpha) \Psi_{\alpha}

where :math:`g(\alpha)` is a reasonably smooth function (e.g. without poles) which is non-vanishing within a finite range of energies. We can then demand that the time limits

.. math::

	\exp(-\ifrak \tau H) \int d\alpha ~g(\alpha) \Psi_{\alpha}^{\pm}
		= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm}

as :math:`\tau \to \pm\infty`, respectively, approach the corresponding superpositions of non-interacting particle states.

To be more precise, let's split the Hamiltonian into the free part and the interaction part as follows

.. math::
	:label: eq_h_as_h0_plus_v

	H = H_0 + V

such that the energy eigenstates :math:`\Phi_{\alpha}` of :math:`H_0` (in the same frame as :math:`\Psi_{\alpha}^{\pm}`) transform according to :eq:`eq_lorentz_transformation_formula_for_many_free_particles`. Then the asymptotic freeness translates into the following conditions

.. math::
	:label: eq_in_out_states_asymptotic_by_energy

	\lim_{\tau \to \pm\infty} \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm} = \lim_{\tau \to \pm\infty} \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha}

or equivalently in terms of the Hamiltonians

.. math::
	:label: eq_in_out_states_asymptotic_by_hamiltonian

	\lim_{\tau \to \pm\infty} \exp(-\ifrak \tau H) \int d\alpha ~g(\alpha) \Psi_{\alpha}^{\pm}
		= \lim_{\tau \to \pm\infty} \exp(-\ifrak \tau H_0) \int d\alpha ~g(\alpha) \Phi_{\alpha}

This motivates the following definition

.. math::
	:label: eq_defn_of_Omega

	\Omega(\tau) \coloneqq \exp(\ifrak \tau H) \exp(-\ifrak \tau H_0)

so that :math:`\Psi_{\alpha}^{\pm} = \Omega(\pm\infty) \Phi_{\alpha}`, at least formally. Moreover, since :math:`\Omega` is unitary, the in- and out-states :math:`\Psi_{\alpha}^{\pm}` are normalized as long as :math:`\Phi_{\alpha}` are normalized.

In practice it will be assumed that the interaction term :math:`V` in :eq:`eq_h_as_h0_plus_v` is relatively small so that a formal solution as power series in :math:`V` may be meaningful. As the first step, let's try to apply :eq:`eq_h_as_h0_plus_v` to :math:`\Psi_{\alpha}^{\pm}` as follows

.. math::

	E_{\alpha} \Psi_{\alpha}^{\pm}
		= H \Psi_{\alpha}^{\pm}
		= (H_0 + V) \Psi_{\alpha}^{\pm}
		\implies (E_{\alpha} - H_0) \Psi_{\alpha}^{\pm} = V \Psi_{\alpha}^{\pm}

Note that :math:`\Phi_{\alpha}` is also annihilated by :math:`E_{\alpha} - H_0`. Considering the asymptotic :eq:`eq_in_out_states_asymptotic_by_energy` or :eq:`eq_in_out_states_asymptotic_by_hamiltonian`, it's reasonable to guess the following formal solution

.. math::
	:label: eq_lippmann_schwinger_mixed

	\Psi_{\alpha}^{\pm} = \Phi_{\alpha} + (E_{\alpha} - H_0 \mp \ifrak \epsilon)^{-1} V \Psi_{\alpha}^{\pm}

where the infinitesimal :math:`\mp \ifrak \epsilon` is a mathematical trick added to avoid division by zero, and the signs will be justified momentarily. One can obviously apply :eq:`eq_lippmann_schwinger_mixed` recursively to get an expansion of :math:`\Psi_{\alpha}^{\pm}` as a power series in :math:`V`, and we shall come back to this point later. In order to express :math:`\Psi_{\alpha}^{\pm}` in terms of :math:`\Phi_{\alpha}`, let's expand the right-hand-side of :eq:`eq_lippmann_schwinger_mixed` as follows

.. math::
	:label: eq_lippmann_schwinger_pure

	\Psi_{\alpha}^{\pm} = \Phi_{\alpha} + \int d\beta ~\frac{(\Phi_{\beta}, V \Psi_{\alpha}^{\pm}) \Phi_{\beta}}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon}

Both :eq:`eq_lippmann_schwinger_mixed` and :eq:`eq_lippmann_schwinger_pure` are known as the `Lippmann-Schwinger equation <https://en.wikipedia.org/wiki/Lippmann%E2%80%93Schwinger_equation>`_.

Now let's justify the term :math:`\pm \ifrak \epsilon` by showing that :eq:`eq_lippmann_schwinger_pure` indeed satisfies the asymptotic condition :eq:`eq_in_out_states_asymptotic_by_energy` as follows

.. math::
	:label: eq_packet_expansion_by_lippmann_schwinger

	\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm}
		&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha} \\
		&\quad + \int d\alpha d\beta ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm}) \Phi_{\beta}}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon} \\
		&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha} \\
		&\quad + \int d\beta ~\Phi_{\beta} \blue{\int d\alpha ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm})}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon}}

Now the integral colored in blue can be integrated over :math:`E_{\alpha}` by a contour that runs from :math:`-\infty` to :math:`+\infty`, followed by a semicircle at infinity, in the upper-half-plane in the case of :math:`\Psi_{\alpha}^-` and the lower-half-plane in the case of :math:`\Psi_{\alpha}^+`, back to :math:`-\infty`. In either case, the sign in :math:`\mp \ifrak \epsilon` is chosen so that the integrant has no poles with infinitesimally small imaginary part, though both :math:`g(\alpha)` and :math:`(\Phi_{\beta}, V \Psi_{\alpha}^{\pm})`, viewed as complex functions, may have poles with finite imaginary parts. It follows then from the residual theorem and the damping factor :math:`\exp(-\ifrak \tau E_{\alpha})` as :math:`\tau \to \pm\infty` that the integral in blue vanishes, as desired.


.. _sec_s_matrix_and_its_symmetry:

S-matrix and its Symmetry
-------------------------

The `S-matrix <https://en.wikipedia.org/wiki/S-matrix>`_ defined by

.. math::
	:label: eq_defn_s_matrix_by_in_and_out_states

	S_{\beta \alpha} \coloneqq \left( \Psi_{\beta}^+, \Psi_{\alpha}^- \right)

records the probability amplitude of finding the out-state :math:`\Psi_{\beta}^+` given the in-state :math:`\Psi_{\alpha}^-`. Note that since the in- and out-states both form an orthonormal basis of the same Hilbert space, the S-matrix is unitary. However, the way :math:`S` is defined in :eq:`eq_defn_s_matrix_by_in_and_out_states` disqualifies it as an operator on the Hilbert space. Therefore it'll be convenient to convert both in- and out-states to the free states and define the *S-operator* by

.. math::
	:label: eq_defn_s_operator

	(\Phi_{\beta}, S \Phi_{\alpha}) \coloneqq S_{\beta \alpha}

Using :eq:`eq_defn_of_Omega` we see that

.. math::
	:label: eq_s_operator_by_u

	& \phantom{\implies} S_{\beta \alpha}
		= \left( \Omega(\infty) \Phi_{\beta}, \Omega(-\infty) \Phi_{\alpha} \right)
		= \left( \Phi_{\beta}, \Omega^{\dagger}(\infty) \Omega(-\infty) \Phi_{\alpha} \right) \\
	& \implies S = \Omega^{\dagger}(\infty) \Omega(-\infty) \eqqcolon U(\infty, -\infty)

where

.. math::
	:label: eq_defn_u_operator

	U(\tau_1, \tau_0) = \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)

The most straightforward way to calculate :math:`S_{\beta \alpha}` is probably to use :eq:`eq_lippmann_schwinger_pure` directly. However. this turns out to be rather involved, and doesn't lead to a simple result. The issue is that we don't really want to convert both the in- and out-states to the non-interacting states, but rather to push, say, the in-states from the far past to the far future and compare with the out-states. To spell out the details, let's first calculate the asymptotic of the in-packet as :math:`\tau \to \infty` (but omitting the :math:`\lim_{\tau \to \infty}` symbol) using :eq:`eq_packet_expansion_by_lippmann_schwinger`

.. math::
	:label: eq_positive_limit_of_in_state_by_lippmann_schwinger

	& \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- \\
		&\quad = \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} + \int d\beta ~\Phi_{\beta} \int d\alpha \frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)}{E_{\alpha} - E_{\beta} + \ifrak \epsilon} \\
		&\quad = \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} \\
		&\qquad - 2\pi\ifrak \int d\beta ~\Phi_{\beta} \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) \exp(-\ifrak \tau E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \\
		&\quad = \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \left( g(\beta) - 2\pi\ifrak \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \right) \\
		&\quad = \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \left( \blue{\delta(\alpha - \beta) - 2\pi\ifrak \delta(E_{\alpha} - E_{\beta}) (\Phi_{\beta}, V \Psi_{\alpha}^-)} \right)

where we've used the residue theorem again in the second equality. Next expand the left-hand-side of the equation in terms of the out-states and then let :math:`\tau \to \infty`

.. math::
	:label: eq_positive_limit_of_in_state_by_expanding_out_states

	\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^-
		&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \int d\beta ~(\Psi_{\beta}^+, \Psi_{\alpha}^-) \Psi_{\beta}^+ \\
		&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Psi_{\beta}^+ \int d\alpha ~g(\alpha) S_{\beta \alpha} \\
		&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \blue{S_{\beta \alpha}}

where we've used the fact that the S-matrix contains a :math:`\delta(E_{\alpha} - E_{\beta})` factor by energy conservation in the second equality, and the defining property :eq:`eq_in_out_states_asymptotic_by_energy` of the out-state in the third equality.

Equating the blue terms from :eq:`eq_positive_limit_of_in_state_by_lippmann_schwinger` and :eq:`eq_positive_limit_of_in_state_by_expanding_out_states`, we've derived the following formula

.. math::
	:label: eq_s_matrix_pre_born_approx

	S_{\beta \alpha} = \delta(\beta - \alpha) - 2\pi\ifrak \delta(E_{\beta} - E_{\alpha}) (\Phi_\beta, V \Psi_{\alpha}^-)

Up to the first order in :math:`V`, one can replace :math:`\Psi_{\alpha}^-` on the right-hand-side by :math:`\Phi_{\alpha}` and arrive at the so-called `Born approximation <https://en.wikipedia.org/wiki/Born_approximation>`_ of the S-matrix.


.. _sec_s_matrix_lorentz_symmetry:

Lorentz symmetry
++++++++++++++++

Recall that in :eq:`eq_lorentz_transformation_formula_for_many_free_particles`, or really in :ref:`Lorentz symmetry of one-particle states <sec_lorentz_symmetry>`, we understood how Lorentz transformations act on particle states. Now we'd like to understand how they act on the S-matrix. Of course, since :math:`U(\Lambda, a)` is unitary, we always have

.. math::

	S_{\beta \alpha} = (\Psi_{\beta}^+, \Psi_{\alpha}^-) = \left( U(\Lambda, a) \Psi_{\beta}^+, U(\Lambda, a) \Psi_{\alpha}^- \right)

but this is *not* what we mean by Lorentz symmetry. What we do want to know is, just like in :eq:`eq_lorentz_transformation_formula_for_many_free_particles`, how Lorentz transformation acts on the particle states, i.e., the (compound) indexes :math:`\alpha` and :math:`\beta`. Now although :eq:`eq_lorentz_transformation_formula_for_many_free_particles` doesn't work for general (interacting) states, it does work for, say, :math:`\Psi_{\alpha}^-` in the :math:`\tau \to -\infty` limit because of the asymptotic freeness. By Lorentz we mean that :math:`U(\Lambda, a)` acts the same way on both in- and out-states. In other words, we'll be looking for some :math:`U(\Lambda, a)` such that the following general formula holds.

.. math::
	:label: eq_lorentz_transformation_formula_for_s_matrix

	& S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} \\
	&\quad = \blue{\exp\left( \ifrak a^{\mu} {\Lambda_{\mu}}^{\nu} \left( (p'_1)_{\nu} + (p'_2)_{\nu} + \cdots - (p_1)_{\nu} - (p_2)_{\nu} - \cdots \right) \right)} \\
	&\qquad \times \sqrt{\frac{(\Lambda p'_1)_0 (\Lambda p'_2)_0 \cdots (\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p'_1)_0 (p'_2)_0 \cdots (p_1)_0 (p_2)_0 \cdots}} \\
	&\qquad \times \sum_{\bar{\sigma}'_1 \bar{\sigma}'_2 \cdots} D^{\ast}_{\bar{\sigma}'_1 \sigma'_1} (W(\Lambda, p'_1)) D^{\ast}_{\bar{\sigma}'_2 \sigma'_2} (W(\Lambda, p'_2)) \cdots \\
	&\qquad \times \sum_{\bar{\sigma}_1 \bar{\sigma}_2 \cdots} D_{\bar{\sigma}_1 \sigma_1} (W(\Lambda, p_1)) D_{\bar{\sigma}_2 \sigma_2} (W(\Lambda, p_2)) \cdots \\
	&\qquad \times S_{\Lambda p'_1, \bar{\sigma}'_1, n'_1; ~\Lambda p'_2, \bar{\sigma}'_2, n'_2; ~\cdots, ~~\Lambda p_1, \bar{\sigma}_1, n_1; ~\Lambda p_2, \bar{\sigma}_2, n_2, ~\cdots}

where we've used primes to distinguish between labels from in- and out-states, and bars to distinguish between labels, specifically the spin-:math:`z` or helicity, before and after the Lorentz transformation.

Since the left-hand-side doesn't depend on the translation parameter :math:`a`, the blue term on the right-hand-side must be :math:`1`. In other words,

.. math:: p_1 + p_2 + \cdots = p'_1 + p'_2 + \cdots

which is nothing but the conservation of (total) momentum. Note that a special case, which is the energy conservation, has already been used in the derivation of :eq:`eq_positive_limit_of_in_state_by_expanding_out_states` from the previous section.

As a consequence, we can now extract a delta function from the S-matrix as follows

.. math::
	:label: eq_s_matrix_with_m

	S_{\beta \alpha} \eqqcolon \delta(\beta - \alpha) - 2\pi\ifrak M_{\beta \alpha} \delta^4 (p_{\beta} - p_{\alpha})

which should be compared with :eq:`eq_s_matrix_pre_born_approx`.

Back to the core question of this section, how in the world can one engineer a magic :math:`U(\Lambda, a)` to satisfy the monstrous :eq:`eq_lorentz_transformation_formula_for_s_matrix`? One cannot. But remember that :eq:`eq_lorentz_transformation_formula_for_s_matrix` is readily satisfied for non-interacting particles. It follows that if we consider instead the S-operator defined by :eq:`eq_defn_s_operator`, and let :math:`U_0(\Lambda, a)` be the Lorentz transformation on free particles defined by :eq:`eq_lorentz_transformation_formula_for_many_free_particles`, then :eq:`eq_lorentz_transformation_formula_for_s_matrix` would be satisfied if :math:`U_0(\Lambda, a)` commutes with :math:`S`. Indeed, using shorthand notations, we have

.. math::

	S_{\beta \alpha} = \left( \Phi_{\beta}, S \Phi_{\alpha} \right) \
		= \left( U_0 \Phi_{\beta}, U_0 S \Phi_{\alpha} \right) \
		= \left( U_0 \Phi_{\beta}, S U_0 \Phi_{\alpha} \right)

where the last quantity is nothing but the right-hand-side in :eq:`eq_lorentz_transformation_formula_for_s_matrix`, as desired.

Now in order for :math:`S` to commute with :math:`U_0(\Lambda, a)`, it suffices that it commutes with the infinitesimal generators of :math:`U_0(\Lambda, a)`, namely,

.. math::
	:label: eq_s_commutes_with_hpjk

	\begin{alignat*}{2}
		&[H_0, S] &&= 0 \\
		&[\Pbf_0, S] &&= 0 \\
		&[\Jbf_0, S] &&= 0 \\
		&[\Kbf_0, S] &&= 0
	\end{alignat*}

where :math:`H_0, \Pbf_0, \Jbf_0, \Kbf_0` are discussed in :ref:`sec_quantum_lorentz_symmetry` and satisfy the commutation relations :eq:`eq_poincare_algebra`.

This shall be done in three steps, where the commutation between :math:`S` and :math:`P_0, J_0` will be handled first, followed by :math:`K_0`, and finally :math:`H_0`.

Step 1.
	Recall from :eq:`eq_s_operator_by_u` and :eq:`eq_defn_u_operator` that the S-operator can be understood as a composition of time translations governed by :math:`H` and :math:`H_0`.  It's therefore necessary to understand how the free infinitesimal Lorentz transformations commute with :math:`H`. To this end, let's consider the in-states at :math:`\tau \to -\infty`, which is approximately free. There we can similarly define infinitesimal operators :math:`\Pbf, \Jbf, \Kbf` that together with :math:`H` satisfy the same commutation relations :eq:`eq_poincare_algebra`.

	Now comes the crucial part, which is to make assumptions about :math:`H` so that :eq:`eq_s_commutes_with_hpjk` are satisfied. Recall from :eq:`eq_h_as_h0_plus_v` that :math:`H = H_0 + V` where :math:`V` describes the interactions. The first assumption we'll make is the following

	.. admonition:: Assumption on :math:`H` for Lorentz invariance of S-matrix #1

		The interaction :math:`V` affects neither the momentum :math:`\Pbf` nor the angular momentum :math:`\Jbf`. In other words, we assume that

		.. math::
			:label: eq_s_matrix_lorentz_invariance_assump_1

			\Pbf = \Pbf_0, ~~\Jbf = \Jbf_0, ~\text{ and }~ [V, \Pbf_0] = [V, \Jbf_0] = 0

	It follows readily from this assumption that :math:`H`, and hence :math:`S`, commutes with :math:`\Pbf_0` and :math:`\Jbf_0`.

Step 2.
	Next we turn to :math:`\Kbf`. This time we cannot "cheat" by assuming that :math:`\Kbf = \Kbf_0` because it would led to the undesirable consequence :math:`H = H_0` by :eq:`eq_poincare_algebra`. So instead, let's write

	.. math:: \Kbf = \Kbf_0 + \Wbf
		:label: eq_k_as_k0_plus_w

	where :math:`\Wbf` denotes the perturbation term. Let's calculate

	.. math::

		[\Kbf_0, S] &= \lim_{\substack{\tau_0 \to -\infty \\ \tau_1 \to \infty\phantom{-}}} [\Kbf_0, U(\tau_1, \tau_0)] \\
			&= \lim_{\substack{\tau_0 \to -\infty \\ \tau_1 \to \infty\phantom{-}}} [\Kbf_0, \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)]

	as follow. First using :eq:`eq_poincare_algebra` again, we have

	.. math::

		\begin{alignat*}{2}
			[\Kbf_0, \exp(\ifrak \tau H_0)] &= [\Kbf_0, \ifrak \tau H_0] \exp(\ifrak \tau H_0) &&= \tau \Pbf_0 \exp(\ifrak \tau H_0) \\
			[\Kbf, \exp(\ifrak \tau H)] &= [\Kbf, \ifrak \tau H] \exp(\ifrak \tau H) &&= \tau \Pbf \exp(\ifrak \tau H) =\tau \Pbf_0 \exp(\ifrak \tau H)
		\end{alignat*}

	from which we can calculate

	.. math::
		:label: eq_k30_u_commutation

		[\Kbf_0, U(\tau_1, \tau_0)]
			&= [\Kbf_0, \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)] \\
			&= [\Kbf_0, \exp(\ifrak \tau_1 H_0)] \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0) \\
			&\quad + \exp(\ifrak \tau_1 H_0) [\Kbf - \Wbf, \exp(\ifrak (\tau_0 - \tau_1) H)] \exp(-\ifrak \tau_0 H_0) \\
			&\quad + \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) [\Kbf_0, \exp(-\ifrak \tau_0 H_0)] \\
			&= \blue{\tau_1 \Pbf_0 \exp(\ifrak \tau H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \\
			&\quad \blue{+ (\tau_0 - \tau_1) \Pbf_0 \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \\
			&\quad \blue{- \tau_0 \Pbf_0 \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \\
			&\quad - \exp(\ifrak \tau_1 H_0) [\Wbf, \exp(\ifrak (\tau_0 - \tau_1) H)] \exp(-\ifrak \tau_0 H_0) \\
			&= -\Wbf(\tau_1) U(\tau_1, \tau_0) + U(\tau_1, \tau_0) \Wbf(\tau_0)

	where :math:`\Wbf(\tau) \coloneqq \exp(\ifrak \tau H_0) \Wbf \exp(-\ifrak \tau H_0)`. Note that the three blue terms cancel out.

	We see that :eq:`eq_k30_u_commutation` would vanish as :math:`\tau \to \pm\infty` if :math:`W(\tau) \to 0`. The latter, in turn, would follow from the following assumption

	.. admonition:: Assumption on :math:`H` for Lorentz invariance of S-matrix #2

		The matrix elements of :math:`W` with respect to the eigenstates :math:`\Phi_{\alpha}` of :math:`H_0` is smooth, so that :math:`W(\tau)` vanishes on any local packet of :math:`\Phi_{\alpha}` as in :eq:`eq_psi_packet` as :math:`\tau \to \pm\infty`.

	This assumption is line with :eq:`eq_in_out_states_asymptotic_by_energy` and :eq:`eq_in_out_states_asymptotic_by_hamiltonian`, which are fundamental to the validity of S-matrix theory.

Step 3.
	Finally let's handle the commutation between :math:`H_0` and :math:`S`. Recall from :eq:`eq_s_operator_by_u` that :math:`S = \Omega^{\dagger}(\infty) \Omega(-\infty)`. Hence the idea is to work out how :math:`H` and :math:`H_0` intertwine with :math:`\Omega(\pm\infty)`. To this end, let's use :eq:`eq_k30_u_commutation` by setting :math:`\tau_1 = 0` and :math:`\tau_0 = \mp\infty` as follows

	.. math:: [\Kbf_0, \Omega(\mp \infty)] = -\Wbf \Omega(\mp \infty) \implies \Kbf \Omega(\mp\infty) = \Omega(\mp\infty) \Kbf_0

	Moreover, by :eq:`eq_s_matrix_lorentz_invariance_assump_1`, we have also :math:`\Pbf \Omega(\mp\infty) = \Omega(\mp\infty) \Pbf_0`. Finally, using the commutation relations :eq:`eq_poincare_algebra` again we conclude that

	.. math:: H \Omega(\mp\infty) = \Omega(\mp\infty) H_0

	which completes the proof of :eq:`eq_s_commutes_with_hpjk`.

.. note::
	Besides showing that :eq:`eq_s_commutes_with_hpjk` hold, our calculations actually establish the following intertwining identities

	.. math::

		H \Omega(\pm\infty) &= \Omega(\pm\infty) H_0 \\
		\Pbf \Omega(\pm\infty) &= \Omega(\pm\infty) \Pbf_0 \\
		\Jbf \Omega(\pm\infty) &= \Omega(\pm\infty) \Jbf_0

	which imply, in particular, that the standard commutation relations :eq:`eq_poincare_algebra` also hold in a frame where :math:`\tau \to \infty`, as expected.

Internal symmetry
+++++++++++++++++

As we noticed in :eq:`eq_lorentz_transformation_formula_for_many_free_particles`, the Lorentz symmetry doesn't act on the labels :math:`n`. An internal symmetry, on the other hand, is a symmetry that leaves :math:`p` and :math:`\sigma` invariant and acts on the other labels such as charge, spin, and so on. We can write the general form of an internal symmetry on in- and out-states as follows

.. math::
	:label: eq_internal_symmetry_transformation_for_in_and_out_states

	U(T) \Psi^{\pm}_{p_1, \sigma_1, n_1;~p_2, \sigma_2, n_2;~\cdots}
		= \sum_{n'_1, n'_2, \cdots} \Dscr_{n'_1 n_1}(T) \Dscr_{n'_2 n_2}(T) \cdots \Psi^{\pm}_{p_1, \sigma_1, n'_1;~p_2, \sigma_2, n'_2;~\cdots}

where :math:`U(T)` is the unitary operator associated with the symmetry transformation :math:`T`, and the :math:`\Dscr`'s are analogs of the little group representations from :eq:`eq_d_repr_of_little_group`.

Similar to :eq:`eq_lorentz_transformation_formula_for_s_matrix`, we can formulate the internal symmetry of S-matrix as follows

.. math::
	:label: eq_internal_symmetry_transformation_formula_for_s_matrix

	S_{n'_1, n'_2, \cdots, ~n_1, n_2, \cdots} = \sum_{\bar{n}'_1, \bar{n}'_2, \cdots, \bar{n}_1, \bar{n}_2, \cdots}
		& \Dscr^{\ast}_{\bar{n}'_1 n'_1}(T) \Dscr^{\ast}_{\bar{n}'_2 n'_2}(T) \cdots \\
		& \times \Dscr_{\bar{n}_1 n_1}(T) \Dscr_{\bar{n}_2 n_2}(T) \cdots
		S_{\bar{n}'_1, \bar{n}'_2, \cdots, ~\bar{n}_1, \bar{n}_2, \cdots}

where we have suppressed the irrelevant :math:`p` and :math:`\sigma` labels.

For what kind of Hamiltonian :math:`H` does there exist an internal symmetry :math:`U(T)` that acts like :eq:`eq_internal_symmetry_transformation_for_in_and_out_states`? The answer is similar to the case of Lorentz symmetry. Namely, if we can split :math:`H = H_0 + V` into the free and perturbation terms, such that the free symmetry transformation :math:`U_0(T)`, which satisfies :eq:`eq_internal_symmetry_transformation_for_in_and_out_states` with :math:`\Phi` in place of :math:`\Psi^{\pm}`, commutes with both :math:`H_0` and :math:`V`.


Similar to the translations in Lorentz symmetry, let's consider a symmetry :math:`T(\theta)` parametrized by a real number. It follows from :eq:`eq_additive_symmetry` that we can write

.. math:: U(T(\theta)) = \exp(\ifrak \theta Q)
	:label: eq_charge_internal_symmetry


where :math:`Q` is a Hermitian operator called the charge. Probably the best known example of it is the electric charge. In this case, we can also write

.. math:: \Dscr_{n n'}(T(\theta)) = \delta_{n n'} \exp(\ifrak \theta q_n)

The general formula :eq:`eq_internal_symmetry_transformation_formula_for_s_matrix` then translates into

.. math:: q_1 + q_2 + \cdots = q'_1 + q'_2 + \cdots

which is nothing about the conservation of charges. Besides the electric charge, there exist also other similar conserved, or at least approximately conserved, quantities, such as baryon number and lepton number.


.. _sec_parity_symmetry:

Parity symmetry
+++++++++++++++

Our considerations on Lorentz symmetry has so far been restricted to proper orthochronous Lorentz transformations. Let's consider the effect of the spatial inversion on the S-matrix now. Recall from :ref:`sec_space_inversion_for_massive_particles` that for non-interacting massive particles

.. math::
	:label: eq_space_inversion_acts_on_in_and_out_states

	U(\Pcal) \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}
		= \eta_{n_1} \eta_{n_2} \cdots \Psi^{\pm}_{U(\Pcal)p_1, \sigma_1, n_1; ~U(\Pcal)p_2, \sigma_2, n_2; ~\cdots}

where :math:`\eta_n` denotes the intrinsic parity of particle :math:`n`. The S-matrix version of the parity symmetry is as follows

.. math::
	:label: eq_space_inversion_formula_for_s_matrix

	& S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} \\
	& \quad = \eta^{\ast}_{n'_1} \eta^{\ast}_{n'_2} \cdots \eta_{n_1} \eta_{n_2} \cdots
		S_{U(\Pcal)p'_1, \sigma'_1, n'_1; ~U(\Pcal)p'_2, \sigma'_2, n'_2; ~\cdots, ~~U(\Pcal)p_1, \sigma_1, n_1; ~U(\Pcal)p_2, \sigma_2, n_2; ~\cdots}

While the space inversion operator :math:`\Pcal` is defined explicitly in :eq:`eq_space_inversion`, the parity operator :math:`U(\Pcal)` is only characterized by :eq:`eq_space_inversion_acts_on_in_and_out_states` and :eq:`eq_space_inversion_formula_for_s_matrix`. In particular, it's not uniquely determined if the particle species under question possesses internal symmetries as discussed in the previous section, because their composition with :math:`\Pcal` will also satisfy :eq:`eq_space_inversion_acts_on_in_and_out_states` and :eq:`eq_space_inversion_formula_for_s_matrix`, and therefore may equally well be called a parity operator.

Since :math:`\Pcal^2 = 1`, it's an obvious question to ask whether :math:`U(\Pcal)^2 = 1` necessarily. This would have been the case if :math:`U` furnishes a genuine representation, but it doesn't have to. In general, we have

.. math::

	U(\Pcal)^2 \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \eta_{n_1}^2 \eta_{n_2}^2 \cdots
		\Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}

which looks just like an internal symmetry. Now if :math:`U(\Pcal)^2` belongs to a continuous family of internal symmetries, then it may be redefined, by suitably composing with internal symmetries, so that all :math:`\eta^2 = 1`. Examples of this kind include notably protons and neutrons. On the other hand, non-examples, i.e., those whose intrinsic parity cannot be reduced to :math:`\pm 1`, include only those hypothetical `Majorana fermions <https://en.wikipedia.org/wiki/Majorana_fermion>`_.

.. todo::
	Revise this part after I've learned more...

.. _dropdown_parities_of_elementary_particles:

.. dropdown:: Parities of elementary particles
	:icon: unlock
	:animate: fade-in-slide-down

	We shall in this section assume familiarity with angular momentum as discussed in :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`.

	Can parity be other than :math:`\pm 1`?
		Aside from electric charges, there exist other quantities that are (at least approximately) conserved by internal symmetries, for example, `baryon numbers <https://en.wikipedia.org/wiki/Baryon_number>`_ :math:`B` and `lepton numbers <https://en.wikipedia.org/wiki/Lepton_number>`_ :math:`L`. Examples of baryons include `protons <https://en.wikipedia.org/wiki/Proton>`_ and `neutrons <https://en.wikipedia.org/wiki/Neutron>`_. Examples of leptons include `electrons <https://en.wikipedia.org/wiki/Electron>`_, `muons <https://en.wikipedia.org/wiki/Muon>`_ and `neutrinos <https://en.wikipedia.org/wiki/Neutrino>`_. The internal symmetry operator generalizes :eq:`eq_charge_internal_symmetry` in a straightforward way as follows

		.. math:: U(T(\alpha, \beta, \gamma)) = \exp(\ifrak(\alpha B + \beta L + \gamma Q))
			:label: eq_baryon_lepton_charge_internal_symmetry

		so that :math:`T` is isomorphic to :math:`\Rbb^3` instead of :math:`\Rbb`. This will be the most general internal symmetry that will be considered here.

		By the conservation of angular momentum, the parity of the number of half-integer spin particles, which we denote by :math:`(-1)^F`, is conserved. Here :math:`F` stands for fermion. For all known (to Weinberg at least) particles, the following equality of parities holds

		.. math:: (-1)^F = (-1)^{B + L}
			:label: eq_fermion_count_eq_baryon_and_lepton_mod_2

		In particular, the above mentioned protons, neutrons, electrons, neutrinos are all spin-:math:`1/2` particles.

		If, for whatever reason, the following holds

		.. math:: \orange{U(\Pcal)^2 = (-1)^F}

		and in addition :eq:`eq_fermion_count_eq_baryon_and_lepton_mod_2` holds, then :math:`U(\Pcal)^2` is part of a continuous symmetry :eq:`eq_baryon_lepton_charge_internal_symmetry` and hence can be set to one.  A hypothetical example that breaks :eq:`eq_fermion_count_eq_baryon_and_lepton_mod_2` is the so-called Majorana fermions that are their own antiparticles, which implies :math:`B = L = 0`. For these particles, we have :math:`U(\Pcal)^4 = 1`, and hence the intrinsic parity may be :math:`\pm 1` or :math:`\pm \ifrak`.

	Can parity be :math:`-1`?
		The following reaction is observed experimentally

		.. math:: \pi^- + d \to n + n
			:label: eq_pion_deuteron_to_two_neutrons

		where a negative pion is absorbed by a `deuteron <https://en.wikipedia.org/wiki/Deuterium>`_ to produce two neutrons. Moreover, the reaction assumes that the initial state, i.e., the left-hand-side of :eq:`eq_pion_deuteron_to_two_neutrons` has orbital angular momentum :math:`\ell = 0` and total angular momentum :math:`j = 1`. Note that the spin of pion and deuteron is :math:`0` and :math:`1`, respectively.

		The conservation of angular momentum demands that the total angular momentum of the right-hand-side must also be :math:`1`, and this can be achieved, a priori, in a number of ways. Since neutrons have spin :math:`1/2`, the total spin :math:`\sfrak` of :math:`n + n` may be either :math:`0` or :math:`1` by :eq:`eq_composite_total_angular_momentum_range`. But since neutrons are fermions and therefore the state :math:`n + n` must be anti-symmetric, we conclude that :math:`\sfrak = 0` by :eq:`eq_second_highest_weight_am_pair_two`. [#pion_deuteron_reaction_final_state]_ Then it follows again from :eq:`eq_composite_total_angular_momentum_range` that the orbital angular momentum of the right-hand-side of :eq:`eq_pion_deuteron_to_two_neutrons` must be :math:`1`. We are left with only one choice.

		Now since the orbital angular momentum changes from :math:`0` in the initial state to :math:`1` in the final state, the S-matrix elements flip sign by the action of :math:`U(\Pcal)` (:red:`WHY? I guess I'm missing knowledge about how orbital angular momentum enters the S-matrix.`). It follows from :eq:`eq_space_inversion_formula_for_s_matrix` that

		.. math:: \eta_{\pi^-} \eta_d = -\eta_n^2

		Deuteron is a nucleus consisting of a proton and a neutron. By the previous discussions about internal symmetries, one can arrange so that they have the same intrinsic parity and hence :math:`\eta_d = \eta_n^2`. It follows that :math:`\eta_{\pi^-} = -1` and the pion :math:`\pi^-` is what we set out to look for. Indeed, all its companions :math:`\pi^0` and :math:`\pi^+` also have parity :math:`-1` due to the isospin symmetry.

		The fact that :math:`\eta_{\pi} = -1` had led to a profound consequence because it was discovered through experiments that there are two spin-:math:`0` particles, now known as :math:`K`-mesons, one of which decays into two pions and the other into three pions. By rotational invariance one can exclude the effects of orbital angular momentum and conclude, assuming parity conservation, that they must have opposite intrinsic parities. However, as more experimental evidence pointing towards the fact that the two :math:`K`-mesons look alike, walk alike and quack alike, it was finally suggested by T. D. Lee and C. N. Yang that they're really the same particle and it's the parity conservation that fails to hold in these reactions, now known as the weak interactions. This suggestion was later verified more directly by an experiment of `C. S. Wu <https://en.wikipedia.org/wiki/Wu_experiment>`_.


Time inversion symmetry
+++++++++++++++++++++++

Recall from :ref:`sec_time_inversion_for_massive_particles` that for a single massive particle

.. math:: U(\Tcal) \Psi_{p, \sigma, n} = \zeta (-1)^{\jfrak - \sigma} \Psi_{\Pcal p, -\sigma, n}

To generalize this to the in- and out-states, we need to remember that the time inversion also interchanges the very frame with respect to which the in- and out-states are defined. The result is as follows

.. math::
	:label: eq_time_inversion_acts_on_in_and_out_states

	U(\Tcal) \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \
		\zeta_{n_1} (-1)^{\jfrak_1 - \sigma_1} \zeta_{n_2} (-1)^{\jfrak_2 - \sigma_2} \cdots \
		\Psi^{\mp}_{\Pcal p_1, -\sigma_1, n_1; ~\Pcal p_2, -\sigma_2, n_2; ~\cdots}

The invariance of S-matrix can then be formulated as follows

.. math::
	:label: eq_time_inversion_acts_on_s_matrix

	& S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} \\
	& \quad = \zeta_{n'_1} (-1)^{\jfrak'_1 - \sigma'_1} \zeta_{n'_2} (-1)^{\jfrak'_2 - \sigma'_2} \cdots \
			\zeta^{\ast}_{n_1} (-1)^{\jfrak_1 - \sigma_1} \zeta^{\ast}_{n_2} (-1)^{\jfrak_2 - \sigma_2}
		\\
	& \qquad \times S_{\Pcal p_1, -\sigma_1, n_1; ~\Pcal p_2, -\sigma_2, n_2; ~\cdots; ~\Pcal p'_1, -\sigma'_1, n'_1; ~\Pcal p'_2, -\sigma'_2, n'_2; ~\cdots}

Since we'll be mainly concerned with the rate of interactions in this section, the phase factors in front of :math:`\Psi` play little role. So let's simplify the notations in :eq:`eq_time_inversion_acts_on_in_and_out_states` and :eq:`eq_time_inversion_acts_on_s_matrix` using compound indexes as follows

.. math::
	:label: eq_time_inversion_formula_for_s_matrix

	U(\Tcal) \Psi^{\pm}_{\alpha} &= \Psi^{\mp}_{\Tcal \alpha} \\
	S_{\beta, \alpha} &= S_{\Tcal\alpha, \Tcal\beta}

where the phase factors have been "absorbed" in the right-hand-side by a re-definition of :math:`\Tcal \alpha` and :math:`\Tcal \beta`.

Unlike the space inversions discussed in the previous section, time inversions don't directly lead to implications on reaction rates because, after all, we cannot turn time around in any experiment. However, under certain circumstances, one can use a trick to draw experimentally verifiable conclusions, which we now present.

The main assumption here is that one can split the S-operator as follows

.. math:: S_{\beta \alpha} = S^{(0)}_{\beta \alpha} + S^{(1)}_{\beta \alpha}
	:label: eq_s_operator_first_order_expansion

such that :math:`S^{(0)}` is also unitary, i.e., :math:`{S^{(0)}}^{\dagger} S^{(0)} = 1`, and is much larger than :math:`S^{(1)}`. Now the unitarity of :math:`S` can be written as follows

.. math:: 1 = S^{\dagger} S = {S^{(0)}}^{\dagger} S^{(0)} + {S^{(0)}}^{\dagger} S^{(1)} + {S^{(1)}}^{\dagger} S^{(0)}

which then implies

.. math:: S^{(1)} = -S^{(0)} {S^{(1)}}^{\dagger} S^{(0)}

Assuming :math:`S^{(0)}` and :math:`S^{(1)}` both satisfy :eq:`eq_time_inversion_formula_for_s_matrix`, the above equality can be rewritten as follows

.. math::
	:label: eq_first_order_s_matrix_is_hermitian

	S^{(1)}_{\beta \alpha} = -\int d\gamma' \int d\gamma ~S^{(0)}_{\beta \gamma'} ~S^{(1)~\ast}_{\Tcal \gamma' ~\Tcal \gamma} ~S^{(0)}_{\gamma \alpha}

where we recall that the adjoint :math:`\dagger` is the composition of the (complex) conjugation and transpose. Together with the unitarity of :math:`S^{(0)}`, we see that the rate of reaction :math:`\left| S^{(1)}_{\beta \alpha} \right|^2`, when summed up against a complete set of :math:`S^{(0)}`-eigenstates, remains the same after applying :math:`\Tcal` to both initial and final states.

The simplest case where :eq:`eq_first_order_s_matrix_is_hermitian` becomes applicable is when both :math:`\alpha` and :math:`\beta` are eigenstates of :math:`S^{(0)}`, with eigenvalues, say, :math:`\exp(\ifrak \theta_{\alpha})` and :math:`\exp(\ifrak \theta_{\beta})`, respectively. In this case :eq:`eq_first_order_s_matrix_is_hermitian` becomes

.. math::

	S^{(1)}_{\beta \alpha} = -\exp(\ifrak(\theta_{\alpha} + \theta_{\beta})) S^{(1)~\ast}_{\Tcal \beta ~\Tcal \alpha}
		\implies \left| S^{(1)}_{\beta \alpha} \right|^2 = \left| S^{(1)}_{\Tcal \beta ~\Tcal \alpha} \right|^2

This is to say that under the assumption that :eq:`eq_s_operator_first_order_expansion` is valid, at least approximately, the rate of reaction :math:`S^{(1)}_{\beta \alpha}` should be invariant under a flip of the :math:`3`-momentum as well as the spin :math:`z`-component. This is *not* contradicted by Wu's experiment which disproved the parity conservation.


Rates and Cross-Sections
------------------------

As we already mentioned, the S-matrix entries :math:`S_{\beta \alpha}` can be interpreted as probability amplitudes of a reaction that turns an in-state :math:`\Psi^-_{\alpha}` into an out-state :math:`\Psi^+_{\beta}`. In other words, the probability :math:`P(\Psi^-_{\alpha} \to \Psi^+_{\beta}) = \left| S_{\beta \alpha} \right|^2`. It is, however, not completely straightforward to square S-matrix entries because, as we've seen in :eq:`eq_s_matrix_with_m`, they contain Dirac delta functions.

Derivation in box model
+++++++++++++++++++++++

One trick that is often used in physics to deal with integration over an infinite space is to restrict the space to a (large) box, often with additional periodic boundary conditions, and hope that the final results will not depend on the size of the box, as long as it's large enough. This is exactly what we shall do.

Consider a cubic box whose sides have length :math:`L` and has volume :math:`V = L^3`. Imposing the periodic boundary condition on the cube, the :math:`3`-momentum is discretized as follows

.. math:: \pbf = \frac{2\pi}{L} (n_1, n_2, n_3)
	:label: eq_momentum_by_wave_number

where :math:`n_1, n_2, n_3` are nonnegative integers. Of course, the higher the :math:`n`, the shorter the wave length if we interpret it as wave mechanics. By analogy with the continuous case, we can define the Dirac delta function as follows

.. math::
	:label: eq_3_momentum_delta_in_a_box

	\delta^3_V (\pbf - \pbf') \coloneqq \frac{1}{(2\pi)^3} \int_V d^3 \xbf ~\exp(\ifrak (\pbf - \pbf') \cdot \xbf) = \frac{V}{(2\pi)^3} \delta_{\pbf \pbf'}

where :math:`\delta_{\pbf \pbf'}` is the usual Kronecker delta. With this setup, the states inner product :eq:`eq_many_particles_state_normalization_rough` will produce, from the Dirac deltas, an overall factor of :math:`\left( V/(2\pi)^3 \right)^N` where :math:`N` denotes the number of particles in the box. In order for the amplitudes to be independent of the size of the box, let's normalize the states as follows

.. math:: \Psi_{\alpha}^{\square} \coloneqq \left( \frac{(2\pi)^3}{V} \right)^{N/2} \Psi_{\alpha}

such that :math:`\left( \Psi^{\square}_{\beta}, \Psi^{\square}_{\alpha} \right) = \delta_{\beta \alpha}` is properly normalized. Correspondingly, we can express the S-matrix with respect to the box-normalized states as follows

.. math:: S^{\square}_{\beta \alpha} = \left( \frac{(2\pi)^3}{V} \right)^{(N_{\alpha} + N_{\beta})/2} S_{\beta \alpha}

where :math:`N_{\alpha}, N_{\beta}` are the numbers of particles in the in- and out-states, respectively.

Now the transition probability in the box model takes the following form

.. math::
	P(\alpha \to \beta) = \left| S^{\square}_{\beta \alpha} \right|^2 = \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha} + N_{\beta}} \left| S_{\beta \alpha} \right|^2

which we can further turn to a differential form as follows

.. math::
	:label: eq_differential_form_of_s_matrix_probability

	dP(\alpha \to \beta)
		&= P(\alpha \to \beta) d\Nscr_{\beta} \\
		&= P(\alpha \to \beta) \left( \frac{V}{(2\pi)^3} \right)^{N_{\beta}} d\beta \\
		&= \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha}} \left| S_{\beta \alpha} \right|^2 d\beta

where :math:`d\beta` denotes an infinitesimal volume element around the state :math:`\beta`, or more precisely, a product of :math:`d^3 \pbf`, one for each particle. Then :math:`\Nscr_{\beta}` counts the number of states within the infinitesimal :math:`d\beta`, which can be readily calculated from :eq:`eq_momentum_by_wave_number`.

Back to our core problem, which is to define :math:`\left| S_{\beta \alpha} \right|^2` as calculated by :eq:`eq_s_matrix_with_m`. The first assumption we will make, at least for now, is a genericity condition

.. _assump_genericity_s_matrix:

.. admonition:: Genericity assumption on the S-matrix
	:class: Important

	No subset of particles in the state :math:`\beta` have exactly the same (total) :math:`4`-momentum as some subset in the state :math:`\alpha`.

Under this assumption, we can remove the term :math:`\delta(\beta - \alpha)` from :eq:`eq_s_matrix_with_m` and write

.. math:: S_{\beta \alpha} = -2 \pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha}
	:label: eq_generic_s_matrix_with_m

and moreover, ensure that :math:`M_{\beta \alpha}` contains no more delta functions. Now the question becomes how to define :math:`\left| \delta^4(p_{\beta} - p_{\alpha}) \right|^2`. In fact, to align with the main theme of using in- and out-states to calculate the S-matrix, the interaction must only be turned on for a finite period of time, say, :math:`T`. Hence the timed delta function can be written as

.. math::
	:label: eq_time_delta_in_a_period

	\delta_T(E_{\beta} - E_{\alpha}) \coloneqq \frac{1}{2 \pi} \int_{-T/2}^{T/2} dt ~\exp(\ifrak (E_{\beta} - E_{\alpha}) t)

We can then modify :eq:`eq_generic_s_matrix_with_m` in a "timed box" as follows

.. math::
	:label: eq_generic_s_matrix_in_time_box

	S_{\beta \alpha} = -2\pi\ifrak \delta^3_V (\pbf_{\beta} - \pbf_{\alpha}) \delta_T(E_{\beta} - E_{\alpha}) M_{\beta \alpha}

Now using :eq:`eq_3_momentum_delta_in_a_box` and :eq:`eq_time_delta_in_a_period`, we can calculate the squares as follows

.. math::

	\begin{alignat*}{2}
		\left( \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \right)^2 &= \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \delta^3_V(0) &&= \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) V/(2\pi)^3 \\
		\left( \delta_T(E_{\beta} - E_{\alpha}) \right)^2 &= \delta_T(E_{\beta} - E_{\alpha}) \delta_T(0) &&= \delta_T(E_{\beta} - E_{\alpha}) T/(2\pi)
	\end{alignat*}

All together, we can now rewrite :eq:`eq_differential_form_of_s_matrix_probability` as follows

.. math::

	dP(\alpha \to \beta) &= \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha}} \left| S_{\beta \alpha} \right|^2 d\beta \\
		&= (2\pi)^2 \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha} - 1} \frac{T}{2\pi}
			\delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \delta_T(E_{\beta} - E_{\alpha}) \left| M_{\beta \alpha} \right|^2 d\beta \\
		&= (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} T \delta^4(p_{\beta} - p_{\alpha}) \left| M_{\beta \alpha} \right|^2 d\beta

where we have restored :math:`\delta^4(p_{\beta} - p_{\alpha})` by taking the large :math:`V` and :math:`T` limits. Dividing by time :math:`T`, the differential rate of transition can be defined as follows

.. math::
	:label: eq_rate_of_reaction_master_formula

	d\Gamma(\alpha \to \beta) \coloneqq dP(\alpha \to \beta) / T = (2\pi)^{3N_{\alpha}-2} V^{1-N_{\alpha}} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta

As will be explained in more detail in the example of the decay of one particle below, it should be kept in mind that :eq:`eq_rate_of_reaction_master_formula` is valid only when large :math:`T` limit can be justified. Nonetheless, such rates are closely related to what actual experiments report.

Examples with few initial particles
+++++++++++++++++++++++++++++++++++

One special case of interest is when :math:`N_{\alpha} = 1`, or in other words, processes where one particle decays into multi-particles. In this case :eq:`eq_rate_of_reaction_master_formula` becomes

.. math::
	:label: eq_differential_reaction_rate_one_particle

	d\Gamma(\alpha \to \beta) = 2\pi \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta

which becomes independent of the volume of the box. This is reasonable because the decay rate of one particle shouldn't care about the size of the containing box. However, the :math:`T \to \infty` limit in :math:`\delta^4(p_{\beta} - p_{\alpha})` is no longer valid. In fact, it cannot be longer than the (mean) lifetime :math:`\tau_{\alpha}` of the particle :math:`\alpha`, because the interaction wouldn't make sense if the particle itself already disintegrates. In this case, in order for :eq:`eq_time_delta_in_a_period` to still approximate a delta function, we must assume that any characteristic energy of the interaction satisfies

.. math:: |E_{\beta} - E_{\alpha}| \ll 1/\tau_{\alpha}

where the right-hand-side is known as the total decay rate.

Another case of interest is when :math:`N_{\alpha} = 2`. In this case :eq:`eq_rate_of_reaction_master_formula` takes the following form

.. math:: d\Gamma(\alpha \to \beta) = (2\pi)^4 V^{-1} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta
	:label: eq_differential_reaction_rate_two_particles

It turns out that in the world of experimentalists, it's more common to use, instead of the transition rate, something called *cross-section*, or equivalently, rate per flux, where the flux is defined as [#abuse_of_phi_as_both_state_vector_and_flux]_

.. math::
	:nowrap:

	\begin{equation*}
		\Phi_{\alpha} \coloneqq u_{\alpha} / V
	\end{equation*}

and :math:`u_{\alpha}` is the (relativistic) relative velocity between the two particles, to be discussed in more detail in the next section by considering Lorentz symmetry. We can then rewrite :eq:`eq_differential_reaction_rate_two_particles` in terms of the cross-section as follows

.. math::
	:label: eq_cross_section_two_particles

	d\sigma(\alpha \to \beta) \coloneqq d\Gamma(\alpha \to \beta) / \Phi_{\alpha} = (2\pi)^4 u_{\alpha}^{-1} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta

Note that :math:`d\sigma` has the dimension of an area.

Lorentz symmetry of rates and cross-sections
++++++++++++++++++++++++++++++++++++++++++++

We can investigate the Lorentz symmetry on the rates and cross-sections as follows. Squaring :eq:`eq_lorentz_transformation_formula_for_s_matrix`, and using the fact that the little group representations are unitary, we see that the following quantity

.. math:: R_{\beta \alpha} \coloneqq \sum_{\text{spins}} |M_{\beta \alpha}|^2 \prod_{\beta} E \prod_{\alpha} E

is Lorentz invariant, where :math:`E = p_0 = \sqrt{\pbf^2 + m^2}` for each particle in :math:`\alpha` and :math:`\beta`, respectively.

It follows that in the one-particle case, :eq:`eq_differential_reaction_rate_one_particle` gives

.. math::

	\sum_{\text{spins}} d\Gamma(\alpha \to \beta) = 2\pi E_{\alpha}^{-1} R_{\beta \alpha} \delta^4(p_{\beta} - p_{\alpha}) \frac{d\beta}{\prod_{\beta} E}

In particular, we recognize :math:`d\beta / \prod_{\beta} E` as a product of the Lorentz invariant :math:`3`-momentum volume elements constructed in :eq:`eq_lorentz_invariant_3_momentum_volume_element`. Hence the only factor in the right-hand-side which is not Lorentz invariant is :math:`E_{\alpha}^{-1}`. It follows that the decay rate of a particle, summed up over all spins, is inverse proportional to its energy, or in other words, a faster moving particle decays slower, which is consistent with the special theory of relativity and experimentally observed slow decay rates of high energy particles coming from cosmic rays.

Next, let's turn to the two-particles case. In this case :eq:`eq_cross_section_two_particles` gives

.. math::

	\sum_{\text{spins}} d\sigma(\alpha \to \beta) = (2\pi)^4 u_{\alpha}^{-1} E_1^{-1} E_2^{-1} R_{\beta \alpha} \delta^4(p_{\beta} - p_{\alpha}) \frac{d\beta}{\prod_{\beta} E}

where :math:`E_1, E_2` are the energies of the two particles in state :math:`\alpha`. As in the one-particle case, in order for the cross-section to be Lorentz invariant, we must define the relative velocity :math:`u_{\alpha}` such that the product :math:`u_{\alpha} E_1 E_2` is Lorentz invariant. Indeed, such a quantity is uniquely determined by the requirement that when one of the particles stays still, then :math:`u_{\alpha}` should be the velocity of the other particle, and it takes the following form

.. math:: u_{\alpha} = \frac{\sqrt{(p_1 \cdot p_2)^2 - m_1^2 m_2^2}}{E_1 E_2}
	:label: eq_scattering_two_particles_relative_velocity

For later use, let's rewrite :math:`u_{\alpha}` in the center-of-mass frame as follows. In the center-of-mass frame, the total momentum vanishes, and therefore we can write :math:`p_1 = (E_1, \pbf)` and :math:`p_2 = (E_2, -\pbf)`. It follows that

.. math::
	:label: eq_two_particles_relative_velocity_in_center_of_mass_frame

	u_{\alpha} &= \frac{\sqrt{(E_1 E_2 + \pbf^2)^2 - m_1^2 m_2^2}}{E_1 E_2} \\
		&= \frac{\sqrt{(E_1 E_2 + \pbf^2)^2 - (E_1^2 - \pbf^2)(E_2^2 - \pbf^2)}}{E_1 E_2} \\
		&= \frac{|\pbf| (E_1 + E_2)}{E_1 E_2} \\
		&= \left| \frac{\pbf_1}{E_1} - \frac{\pbf_2}{E_2} \right|

which indeed looks more like a relative velocity. Note, however, that this is *not* a physical velocity because its value may approach :math:`2` (i.e., faster than the speed of light) in relativistic limit.

The phase-space factor
++++++++++++++++++++++

By phase-space factor we mean the factor :math:`\delta^4(p_{\beta} - p_{\alpha}) d\beta` that appears in transition probabilities, rates and cross-sections discussed above. The goal of this section is to calculate it, particularly in the scenario where the final state consists of two particles. We'll use the center-of-mass frame with respect to the initial state so that :math:`\pbf_{\alpha} = 0`. Then the phase-space factors can be written as follows

.. math::

	\delta^4(p_{\beta} - p_{\alpha}) d\beta = \delta^3(\pbf'_1 + \pbf'_2 + \cdots) \delta(E'_1 + E'_2 + \cdots - E_{\alpha}) d^3 \pbf'_1 d^3 \pbf'_2 \cdots

where we recall that the primes indicate that the quantities are taken from state :math:`\beta`, and :math:`E_{\alpha}` denotes the total energy of state :math:`\alpha`. In the case where the final state consists of exactly two particles, the phase-space factor can be further simplified as follows

.. math::
	:label: eq_simplified_two_final_particles_phase_space_factor

	\delta^4(p_{\beta} - p_{\alpha}) d\beta &= \delta(E'_1 + E'_2 - E_{\alpha}) d^3 \pbf'_1 \\
		&= \delta \left( \sqrt{|\pbf'_1|^2 + {m'_1}^2} + \sqrt{|\pbf'_1|^2 + {m'_2}^2} - E_{\alpha} \right) |\pbf'_1|^2 d|\pbf'_1| d\Omega

where :math:`\Omega` is the solid angle in :math:`\pbf'_1`-space, if in the integration we replace every occurrence of :math:`\pbf'_2` with :math:`-\pbf'_1`.

To further simply the delta function in :eq:`eq_simplified_two_final_particles_phase_space_factor`, we recall the following identity, which is an incarnation of integration by substitution,

.. math:: \delta(f(x)) = \delta(x - x_0) / f'(x_0)

where :math:`x_0` is a simple zero of :math:`f`. In the case of :eq:`eq_simplified_two_final_particles_phase_space_factor`, we let

.. math:: f(|\pbf'_1|) = \sqrt{|\pbf'_1|^2 + {m'_1}^2} + \sqrt{|\pbf'_1|^2 + {m'_2}^2} - E_{\alpha}

so that the simple zero is found at

.. math::
	:label: eq_defn_root_k_prime

	k' = \frac{1}{2E_{\alpha}} \sqrt{\left( E_{\alpha}^2 - {m'_1}^2 - {m'_2}^2 \right)^2 - 4 {m'_1}^2 {m'_2}^2}

Now differentiating :math:`f` at :math:`k'` we get

.. math:: f'(k') = \frac{k'}{E'_1} + \frac{k'}{E'_2} = \frac{k' E_{\alpha}}{E_1 E_2}

where

.. math::
	:label: eq_defn_of_e1_and_e2

	E'_1 &= \sqrt{{k'}^2 + {m'_1}^2} = \frac{E_{\alpha}^2 + {m'_1}^2 - {m'_2}^2}{2E_{\alpha}} \\
	E'_2 &= \sqrt{{k'}^2 + {m'_2}^2} = \frac{E_{\alpha}^2 - {m'_1}^2 + {m'_2}^2}{2E_{\alpha}}

Putting all together, we can further simplify :eq:`eq_simplified_two_final_particles_phase_space_factor` as follows

.. math:: \delta^4(p_{\beta} - p_{\alpha}) d\beta = \frac{k' E'_1 E'_2}{E_{\alpha}} d\Omega
	:label: eq_two_particles_final_state_phase_factor_formula


where :math:`k', E'_1` and :math:`E'_2` are defined by :eq:`eq_defn_root_k_prime` and :eq:`eq_defn_of_e1_and_e2`, respectively.

Substituting :eq:`eq_two_particles_final_state_phase_factor_formula` into :eq:`eq_differential_reaction_rate_one_particle`, we see that in the case of one particle decaying into two particles

.. math::
	:nowrap:

	\begin{equation*}
		\frac{d\Gamma(\alpha \to \beta)}{d\Omega} = \frac{2\pi k' E'_1 E'_2}{E_{\alpha}} |M_{\beta \alpha}|^2
	\end{equation*}

and in the case of a two-body scattering :math:`1~2 \to 1'~2'`, we have, using also :eq:`eq_cross_section_two_particles` and :eq:`eq_two_particles_relative_velocity_in_center_of_mass_frame` (and remembering :math:`E_{\alpha} = E_1 + E_2`), the following

.. math::
	:label: eq_two_body_scattering_cross_section_per_solid_angle

	\frac{d\sigma(\alpha \to \beta)}{d\Omega} = \frac{(2\pi)^4 k' E'_1 E'_2}{u_{\alpha} E_{\alpha}} |M_{\beta \alpha}|^2 \
		= \frac{(2\pi)^4 k' E'_1 E'_2 E_1 E_2}{k E_{\alpha}^2} |M_{\beta \alpha}|^2

where :math:`k \coloneqq |\pbf_1| = |\pbf_2|`. These calculations will be used in the next section to get some insights into the scattering process.


Implications of the unitarity of S-matrix
+++++++++++++++++++++++++++++++++++++++++

In this section we'll no longer assume the :ref:`Genericity of the S-matrix <assump_genericity_s_matrix>`. This means that we'll get back to use :eq:`eq_s_matrix_with_m`, instead of :eq:`eq_generic_s_matrix_with_m`, which we recall as follows

.. math::

	S_{\beta \alpha} = \delta(\beta - \alpha) - 2\pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha}

However, all the calculations from the previous sections can still be used here because we'll be caring about, for example, the *total* rates, which are integrations over all possible final states, and the degenerate ones will *not* contribute to such integrals.

First, let's spell out the consequence of the unitarity of the S-matrix, or more precisely :math:`S^{\dagger} S = 1`, as follows

.. math::
	:label: eq_s_matrix_unitarity_first_half

	\delta(\gamma - \alpha) &= \int d\beta ~S^{\ast}_{\beta \gamma} S_{\beta \alpha} \\
		&= \int d\beta \left( \delta(\beta - \gamma) + 2\pi \ifrak \delta^4(p_{\beta} - p_{\gamma}) M^{\ast}_{\beta \gamma} \right)
			\left( \delta(\beta - \alpha) - 2\pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha} \right) \\
		&= \delta(\gamma - \alpha) + 2\pi \ifrak \delta^4(p_{\alpha} - p_{\gamma}) M^{\ast}_{\alpha \gamma} -
			2\pi \ifrak \delta^4(p_{\gamma} - p_{\alpha}) M_{\gamma \alpha} \\
		&\quad + 4\pi^2 \delta^4(p_{\gamma} - p_{\alpha}) \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) M^{\ast}_{\beta \gamma} M_{\beta \alpha}

which implies

.. math::
	:label: eq_s_matrix_unitarity_implication_on_m_general

	\ifrak M^{\ast}_{\alpha \gamma} - \ifrak M_{\gamma \alpha} + 2\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) M^{\ast}_{\beta \gamma} M_{\beta \alpha} = 0

In the special case where :math:`\alpha = \gamma`, :eq:`eq_s_matrix_unitarity_implication_on_m_general` gives the following key identity, known as the *generalized optical theorem*

.. math::
	:label: eq_s_matrix_unitarity_implication_on_m_special

	\op{Im} M_{\alpha \alpha} = -\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2

As an application we can calculate the total rate of all transitions produced by the initial state :math:`\alpha` using :eq:`eq_rate_of_reaction_master_formula` as follows

.. math::

	\Gamma_{\alpha} &\coloneqq \int d\beta ~\frac{d\Gamma(\alpha \to \beta)}{d\beta} \\
		&~= (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 \\
		&~= -\frac{1}{\pi} (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} \op{Im} M_{\alpha \alpha}

.. dropdown:: An example of two-body scattering
	:icon: unlock
	:animate: fade-in-slide-down

	In the case where :math:`\alpha` is a two-particles state, we can use :eq:`eq_cross_section_two_particles` to calculate the total cross-section as follows

	.. math::
		:label: eq_two_body_total_cross_section_by_m

		\sigma_{\alpha} &\coloneqq \int d\beta ~\frac{d\sigma(\alpha \to \beta)}{d\beta} \\
			&~= (2\pi)^4 u_{\alpha}^{-1} \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 \\
			&~= -16\pi^3 u_{\alpha}^{-1} \op{Im} M_{\alpha \alpha}

	We then recall from :eq:`eq_two_body_scattering_cross_section_per_solid_angle` that :math:`M_{\beta \alpha}` may also be expressed in terms of the differential cross-section by solid angle. Motivated by :eq:`eq_two_body_scattering_cross_section_per_solid_angle`, let's define the *scattering amplitude* as follows

	.. math::

		f(\alpha \to \beta) \coloneqq -\frac{4\pi^2}{E_{\alpha}} \sqrt{\frac{k' E'_1 E'_2 E_1 E_2}{k}} ~M_{\beta \alpha}

	so that :math:`d\sigma(\alpha \to \beta) / d\Omega = |f(\alpha \to \beta)|^2`, and the sign is just a convention. A particularly simple case is when the scattering is *elastic*, which means that :math:`E_i = E'_i` and :math:`|\pbf_i| = |\pbf'_i|` for :math:`i = 1, 2`, and consequentially :math:`k = k'`. In this case we have

	.. math::

		f(\alpha \to \beta) &= -\frac{4\pi^2 E_1 E_2}{E_{\alpha}} M_{\beta \alpha} \\
		u_{\alpha} &= \frac{k E_{\alpha}}{E_1 E_2}

	These, together with :eq:`eq_two_body_total_cross_section_by_m`, imply the following so-called `optical theorem <https://en.wikipedia.org/wiki/Optical_theorem>`_ for elastic two-body scattering

	.. math:: \op{Im} f(\alpha \to \alpha) = \frac{k}{4\pi} \sigma_{\alpha}
		:label: eq_optical_theorem

	From the optical theorem we can derive an estimate on the forward diffraction angle as follows

	.. math::

		\sigma_{\alpha} = \int d\Omega ~|f(\alpha \to \beta)|^2 > \frac{1}{2} |f(\alpha \to \alpha)|^2 \Delta\Omega \
			\geq \frac{1}{2} \left| \op{Im} f(\alpha \to \alpha) \right|^2 \Delta\Omega

	where the second inequality follows from the assumption that :math:`f` is continuous, and the factor :math:`1/2` is a rather random choice and can be anything less than :math:`1`. Roughly speaking :math:`\Delta\Omega` measures the peak of the diffraction in the direction of :math:`\alpha`. Combining with :eq:`eq_optical_theorem`, we conclude that

	.. math:: \Delta\Omega < \frac{32 \pi^2}{k^2 \sigma_{\alpha}}

	Assuming that total decay rate :math:`\sigma_{\alpha}` doesn't grow very fast at high energies, the peak in the direction of :math:`\alpha` shrink at the scale of :math:`1 / k^2`.

Another application of the unitary of the S-matrix is along the lines of statistical mechanics. Applying the same calculation in :eq:`eq_s_matrix_unitarity_first_half` to :math:`S S^{\dagger} = 1`, we get the counterpart to :eq:`eq_s_matrix_unitarity_implication_on_m_special`

.. math:: \op{Im} M_{\alpha \alpha} = -\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\alpha \beta}|^2

Combining with the master equation :eq:`eq_rate_of_reaction_master_formula` we have

.. math::
	:label: eq_s_matrix_unitarity_induced_symmetry_on_reaction_rates

	\int d\beta ~c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} = \int d\beta ~c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha}

where :math:`c_{\alpha} \coloneqq \left( V / (2\pi)^3 \right)^{N_{\alpha}}`.

We shall carry out an equilibrium analysis for state :math:`\alpha`. To this end, let :math:`P_{\alpha} d\alpha` be the infinitesimal probability of finding the system in state :math:`\alpha`. Then we have

.. math::

	\frac{dP_{\alpha}}{dt} = \int d\beta ~P_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} - P_{\alpha} \int d\beta ~\frac{d\Gamma(\alpha \to \beta)}{d\beta}

where the first term calculates the total rate that other states transit into :math:`\alpha`, and the second term calculates the total rate that the state :math:`\alpha` transits into other states. Recall that the *entropy* of the system is defined to be

.. math:: \Scal \coloneqq -\int d\alpha ~P_{\alpha} \ln(P_{\alpha} / c_{\alpha})

Its rate of change can be estimated as follows

.. math::

	\frac{d\Scal}{dt} \
		&= -\int d\alpha ~(\ln(P_{\alpha} / c_{\alpha}) + 1) \frac{dP_{\alpha}}{dt} \\
		&= -\int d\alpha \int d\beta ~(\ln(P_{\alpha} / c_{\alpha}) + 1) \left( P_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \
				- P_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} \right) \\
		&= \int d\alpha \int d\beta ~P_{\beta} \ln\left( \frac{P_{\beta} c_{\alpha}}{P_{\alpha} c_{\beta}} \right) \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \\
		&\geq \int d\alpha \int d\beta ~\left( \frac{P_{\beta}}{c_{\beta}} - \frac{P_{\alpha}}{c_{\alpha}} \right) c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \\
		&= \int d\alpha \int d\beta ~\frac{P_{\beta}}{c_{\beta}} \left( c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} - c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} \right) \\
		&= \int d\alpha ~\frac{P_{\alpha}}{c_{\alpha}} \int d\beta \left( c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} - c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \right) = 0

where the fourth inequality follows from the general inequality :math:`\ln(x) \geq 1 - 1 / x` for any :math:`x > 0`, and the last follows from :eq:`eq_s_matrix_unitarity_induced_symmetry_on_reaction_rates`. This is nothing but the famous slogan: entropy never decreases! And we see that as a consequence of the unitarity of the S-matrix.


Perturbation Theory of S-matrix
-------------------------------

Rather than being the epilogue of :ref:`sec_scattering_theory`, this section is more like a prelude to what comes next. In particular, we will work out a candidate Hamiltonian that satisfies the Lorentz invariance condition discussed in :ref:`sec_s_matrix_and_its_symmetry`.

One possible starting point of the perturbation theory is :eq:`eq_s_matrix_pre_born_approx` together with the Lippmann-Schwinger formula :eq:`eq_lippmann_schwinger_pure` which we recollect as follows

.. math::
	:label: eq_s_matrix_pre_born_approx_repeated

	S_{\beta \alpha} = \delta(\beta - \alpha) - 2\pi \ifrak \delta(E_{\beta} - E_{\alpha}) (\Phi_{\beta}, V\Psi_{\alpha}^-)

where

.. math::
	:label: eq_lippmann_schwinger_repeated

	\Psi_{\alpha}^- = \Phi_{\alpha} + \int d\beta ~\frac{(\Phi_{\beta}, V\Psi_{\alpha}^-) \Phi_{\beta}}{E_{\alpha} - E_{\beta} + \ifrak \epsilon}

Applying :math:`V` to :eq:`eq_lippmann_schwinger_repeated` and taking scalar product with :math:`\Phi_{\beta}`, we get

.. math::
	:label: eq_base_iter_old_fashioned_s_matrix_perturbation

	\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right) = V_{\beta \alpha} + \int d\beta ~\frac{\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right) V_{\beta \alpha}}{E_{\alpha} - E_{\beta} + \ifrak \epsilon}

where :math:`V_{\beta \alpha} \coloneqq \left( \Phi_{\beta}, V\Phi_{\alpha} \right)`. One can now apply :eq:`eq_base_iter_old_fashioned_s_matrix_perturbation` iteratively to get the following power series expansion

.. math::
	:label: eq_s_matrix_power_series_expansion_old_fashioned

	\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right)
		&= V_{\beta \alpha} + \int d\gamma ~\frac{V_{\beta \gamma} V_{\gamma \alpha}}{E_{\alpha} - E_{\gamma} + \ifrak \epsilon} \\
		&\quad + \int d\gamma \int d\gamma' ~\frac{V_{\beta \gamma} V_{\gamma \gamma'} V_{\gamma' \alpha}}{(E_{\alpha} - E_{\gamma} + \ifrak \epsilon)(E_{\alpha} - E_{\gamma'} + \ifrak \epsilon)} + \cdots

and therefore a power series expansion in :math:`V` of :math:`S_{\beta \alpha}` in view of :eq:`eq_s_matrix_pre_born_approx_repeated`.

One obvious drawback of the expansion :eq:`eq_s_matrix_power_series_expansion_old_fashioned` is that it obscures the Lorentz symmetry of the S-matrix because the denominators consist of only the energy terms. To overcome this, we shall use instead the other interpretation of the S-matrix in terms of the Hamiltonians given by :eq:`eq_s_operator_by_u` and :eq:`eq_defn_u_operator`, which we recall as follows

.. math:: S = U(\infty, -\infty)

where

.. math:: U(\tau, \tau_0) = \exp(\ifrak H_0 \tau) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0)
	:label: eq_defn_u_operator_repeated

Now differentiating :eq:`eq_defn_u_operator_repeated` with respect to :math:`\tau` gives

.. math::
	:label: eq_evolution_equation_of_u_operator

	\ifrak \frac{d}{d\tau} U(\tau, \tau_0) &= -H_0 \exp(\ifrak H_0 \tau) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \\
		&\quad + \exp(\ifrak H_0 \tau) H \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \\
		&= \exp(\ifrak H_0 \tau) (H - H_0) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \\
		&= \exp(\ifrak H_0 \tau) V \exp(-\ifrak H_0 \tau) U(\tau, \tau_0) \\
		&\eqqcolon V(\tau) U(\tau, \tau_0)

Here

.. math::
	:label: eq_defn_interaction_perturbation_term

	V(\tau) \coloneqq \exp(\ifrak H_0 \tau) V \exp(-\ifrak H_0 \tau)

is a time-dependent operator in the so-called *interaction picture*, to be distinguished from the Heisenberg picture operator where the true Hamiltonian :math:`H` should be used in place of :math:`H_0`. The differential equation :eq:`eq_evolution_equation_of_u_operator` can be easily solved as follows

.. math::

	U(\tau, \tau_0) = 1 - \ifrak \int_{\tau_0}^{\tau} dt ~V(t) U(t, \tau_0)

which can then be iterated to give the following

.. math::

	U(\tau, \tau_0)
		&= 1 - \ifrak \int_{\tau_0}^{\tau} dt_1 ~V(t_1) + (-\ifrak)^2 \int_{\tau_0}^{\tau} dt_1 \int_{\tau_0}^{t_1} dt_2 ~V(t_1) V(t_2) \\
	  	&\quad + (-\ifrak)^3 \int_{\tau_0}^{\tau} dt_1 \int_{\tau_0}^{t_1} dt_2 \int_{\tau_0}^{t_2} dt_3 ~V(t_1) V(t_2) V(t_3) + \cdots

Letting :math:`\tau \to \infty` and :math:`\tau_0 \to -\infty` we get another power series expansion of :math:`S` in :math:`V` as follows

.. math::
	:label: eq_s_matrix_power_series_expansion_raw

	S &= 1 - \ifrak \int_{-\infty}^{\infty} dt_1 ~V(t_1) + (-\ifrak)^2 \int_{-\infty}^{\infty} dt_1 \int_{-\infty}^{t_1} dt_2 ~V(t_1) V(t_2) \\
		&\quad + (-\ifrak)^3 \int_{-\infty}^{\infty} dt_1 \int_{-\infty}^{t_1} dt_2 \int_{-\infty}^{t_2} dt_3 ~V(t_1) V(t_2) V(t_3) + \cdots

It's somewhat inconvenient that the integral limits in :eq:`eq_s_matrix_power_series_expansion_raw` ruins the permutation symmetry of the products of :math:`V`. But this can be fixed by introducing a *time-ordered product* as follows

.. math::
	:label: eq_defn_time_ordered_product

	T\{ V(t) \} &\coloneqq V(t) \\
	T\{ V(t_1) V(t_2) \} &\coloneqq \theta(t_1 - t_2) V(t_1) V(t_2) + \theta(t_2 - t_1) V(t_2) V(t_1) \\
	T\{ V(t_1) V(t_2) V(t_3) \} &\coloneqq \theta(t_1 - t_2) \theta(t_2 - t_3) V(t_1) V(t_2) V(t_3) + \cdots \\
	&\cdots

where :math:`\theta(\tau)` is the step function which equals :math:`1` for :math:`\tau > 0` and :math:`0` for :math:`\tau < 0`, and it doesn't matter what the value at :math:`\tau = 0` is because it doesn't contribute to the integrals in :eq:`eq_s_matrix_power_series_expansion_raw` anyway. With this definition, we can rewrite :eq:`eq_s_matrix_power_series_expansion_raw` as follows

.. math::
	:label: eq_s_matrix_power_series_expansion_time_ordered

	S = 1 + \sum_{n=1}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 dt_2 \cdots dt_n ~T\{ V(t_1) V(t_2) \cdots V(t_n) \}

where the division by :math:`n!` is to account for the duplicated integrals introduced by the time-ordered product. Note that this power series looks much like the Taylor series of an exponential function. Indeed, in the unlikely event where :math:`V(t)` at different times all commute, one can remove the time-ordering and write :eq:`eq_s_matrix_power_series_expansion_time_ordered` as an exponential function.

One great benefit of writing :math:`S` as in the form of :eq:`eq_s_matrix_power_series_expansion_time_ordered` is that we can reformulate the condition of :math:`S` being Lorentz symmetric in terms of some condition on :math:`V`. Recall from :ref:`sec_s_matrix_and_its_symmetry` that a sufficient condition for a Lorentz invariant S-matrix is that the S-operator commutes with :math:`U_0(\Lambda, a)`, or equivalently in infinitesimal terms :eq:`eq_s_commutes_with_hpjk` are satisfied. Now the main postulation is to express :math:`V` using a density function as follows

.. math:: V(t) = \int d^3 x ~\Hscr(t, \xbf)
	:label: eq_defn_v_by_density

such that :math:`\Hscr(x)` is a scalar in the sense that

.. math:: U_0(\Lambda, a) \Hscr(x) U^{-1}_0(\Lambda, a) = \Hscr(\Lambda x + a)
	:label: eq_h_density_is_scalar

Under these assumptions, we can further rewrite :eq:`eq_s_matrix_power_series_expansion_time_ordered` in terms of :math:`\Hscr(x)` as follows

.. math::
	:label: eq_s_matrix_power_series_expansion_time_ordered_density

	S = 1 + \sum_{n=1}^{\infty} \frac{(-\ifrak)^n}{n!} \int d^4 x_1 \cdots d^4 x_n ~T\{ \Hscr(x_1) \cdots \Hscr(x_n) \}

This expression of :math:`S` is manifestly Lorentz invariant, except for the time-ordering part. In fact, the time-ordering between two spacetime points :math:`x_1, x_2` are Lorentz invariant if and only if :math:`x_1 - x_2` is time-like, namely, :math:`(x_1 - x_2)^2 \geq 0`. This is consistent with intuition because events with time-like (or light-like) separations may be observed by one observer, who definitely should know which event happened first. Therefore we obtain a sufficient condition for the Lorentz invariance of :math:`S` as follows

.. math:: [\Hscr(x_1), \Hscr(x_2)] = 0, \quad\forall ~(x_1 - x_2)^2 \leq 0
	:label: eq_h_commutativity_for_space_like_separations

where we've also included the light-like case for technical reasons that will only become clear later. This condition is also referred to as the *causality* condition as it may be interpreted as saying that interactions happening at space-like separations should not be correlated.

.. dropdown:: A formal proof of the Lorentz invariance of the S-matrix
	:icon: unlock
	:animate: fade-in-slide-down

	We shall verify that the S-operator defined by :eq:`eq_s_matrix_power_series_expansion_time_ordered_density` and satisfying :eq:`eq_h_commutativity_for_space_like_separations` indeed satisfies :eq:`eq_s_commutes_with_hpjk`, or rather, the most nontrivial :math:`[S, \Kbf_0] = 0`. Using the definition of :math:`\Kbf_0` from :ref:`sec_quantum_lorentz_symmetry`, it follows from :eq:`eq_h_density_is_scalar` that

	.. math:: -\ifrak [\Kbf_0, \Hscr(t, \xbf)] = t \nabla \Hscr(t, \xbf) + \xbf \p_t \Hscr(t, \xbf)

	Integrating over :math:`\xbf` and setting :math:`t = 0`, we get

	.. math::

		[\Kbf_0, V] &= \left[ \Kbf_0, \int d^3 x ~\Hscr(0, \xbf) \right] \\
			&= \int d^3 x ~\xbf \left. \frac{\p}{\p t} \right\vert_{t=0} \Hscr(t, \xbf) \\
			&= \left[ H_0, \int d^3 x ~\xbf \Hscr(0, \xbf) \right] \\
			&\eqqcolon [H_0, \Wbf]

	where the third equality follows again from :eq:`eq_h_density_is_scalar` by letting :math:`\Lambda` to be the infinitesimal time translation. Indeed the quantity

	.. math:: \Wbf \coloneqq \int d^3 x ~\xbf \Hscr(0, \xbf)

	here is the same as the :math:`\Wbf` in :eq:`eq_k_as_k0_plus_w`.

	Now for :math:`\Kbf = \Kbf_0 + \Wbf`, together with :math:`\Pbf` and :math:`H`, to satisfy the commutation relation :eq:`eq_poincare_algebra`, which is readily satisfied by the free particle state Lorentz symmetry generators :math:`\Kbf_0, \Pbf_0` and :math:`H_0`, it suffices that

	.. math:: [\Wbf, V] = 0
		:label: eq_w_v_commute

	because of the following calculation

	.. math::

		[H, \Kbf] &= [H_0 + V, \Kbf_0 + \Wbf] \\
			&= \ifrak \Pbf_0 + [V, \Kbf_0] + [H_0, W] + [V, \Wbf] \\
			&= \ifrak \Pbf_0 \\
			&= \ifrak \Pbf

	Finally, we note that :eq:`eq_w_v_commute` is equivalent to the following

	.. math:: \int d^3 x \int d^3 y ~\xbf [\Hscr(0, \xbf), \Hscr(0, \ybf)] = 0

	which is clearly a weaker condition than :eq:`eq_h_commutativity_for_space_like_separations`.

At last we've finally climbed the highest peak in scattering theory, namely :eq:`eq_h_commutativity_for_space_like_separations`. It is specific to the relativistic theory because time-ordering is always preserved in Galilean symmetry. It is also this restriction that eventually leads us to a quantum field theory. In the words of the author

	"It is this condition that makes the combination of Lorentz invariance and quantum mechanics so restrictive." [#weinberg_quote_on_lorentz_invariance_and_quantum_mechanics]_

	-- S. Weinberg

.. rubric:: Footnotes

.. [#in_out_state_sign_convention] I really don't like the sign convention of the in- and out-states in [Wei95]_, even though Weinberg said that it has become a tradition. So our sign convention, as far as the definitions of in- and out-states are concerned, is the opposite to the one used in [Wei95]_.

.. [#pion_deuteron_reaction_final_state] The final state claimed in [Wei95]_ page 126 has total spin :math:`1` rather than :math:`0`. I suspect that the claim in [Wei95]_ is wrong, but otherwise it doesn't affect any subsequent arguments anyway.

.. [#abuse_of_phi_as_both_state_vector_and_flux] It's really unfortunate that one constantly runs out symbols to represent physical quantities. So we have to live with the fact that the meaning of the symbol will depend on the context.

.. [#weinberg_quote_on_lorentz_invariance_and_quantum_mechanics] See [Wei95]_ page 145.
