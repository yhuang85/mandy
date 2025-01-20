.. _sec_quantum_fields_and_antiparticles:

Quantum Fields and Antiparticles
================================

In this chapter we will construct the Hamiltonians in the form of :math:`H = H_0 + V`, where :math:`H_0` is the Hamiltonian of free particles, and :math:`V = \int d^3 x~\Hscr(0, \xbf)` is a (small) interaction term in the form of :eq:`eq_defn_v_by_density`, and the interaction density :math:`\Hscr(x)` is a Lorentz scalar in the sense of :eq:`eq_h_density_is_scalar` and satisfies the cluster decomposition principle :eq:`eq_h_commutativity_for_space_like_separations`. As a byproduct of the construction, we'll also demystify the so-called *antiparticles* which have been mentioned a number of times so far without definition.


Symmetries and Quantum Fields
-----------------------------

Following the discussions in :ref:`sec_the_cluster_decomposition_principle`, we'll construct :math:`\Hscr(x)` out of creation and annihilation operators. However, as we've seen in :ref:`sec_the_lorentz_and_cpt_transformation_laws`, the Lorentz transformation of both :math:`a^{\dagger}(q)` and :math:`a(q)` involve in the coefficients a matrix element :math:`D_{\sigma \sigma'}(W(\Lambda, p))` that depend on the momenta of the particles, and hence are not scalars. The idea, then, is to construct :math:`\Hscr(x)` out of the so-called  annihilation and creation fields defined by

.. math::
	:label: eq_defn_annihilation_and_creation_field

	\psi_{\ell}^+(x) &\coloneqq \sum_{\sigma, n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) a(\pbf, \sigma, n) \\
	\psi_{\ell}^-(x) &\coloneqq \sum_{\sigma, n} \int d^3 p~v_{\ell}(x;~\pbf, \sigma, n) a^{\dagger}(\pbf, \sigma, n)

where :math:`\ell` is reserved for labeling particles later. We see, in particular, that the creation field :math:`\psi_{\ell}^-(x)` is a superposition of creation operators. Applying it to the vacuum state and let :math:`x` wander around the whole spacetime then creates a (quantum) field.

.. note::
	Just like particles, fields may be either bosonic or fermionic, but not mixed. For example :math:`\psi^-_{\ell}(x)` is a bosonic/fermionic if and only if all the particles created by :math:`a^{\dagger}` are bosonic/fermionic, respectively.

Lorentz symmetries
++++++++++++++++++

Now the hope is that the creation and annihilation fields transform, under (proper orthochronous) Lorentz symmetry, by a matrix that is independent of the spacetime coordinate :math:`x`. More precisely, we'd like the following to hold

.. math::
	:label: eq_conjugate_annihilation_and_creation_field

	U_0(\Lambda, b) \psi_{\ell}^+(x) U_0^{-1}(\Lambda, b) &= \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}^+ (\Lambda x + b) \\
	U_0(\Lambda, b) \psi_{\ell}^-(x) U_0^{-1}(\Lambda, b) &= \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}^- (\Lambda x + b)

Note that we've put :math:`\Lambda^{-1}` inside :math:`D_{\ell \ell'}` so that :math:`D` furnishes a representation of the homogeneous Lorentz transformations in the sense that :math:`D(\Lambda_1) D(\Lambda_2) = D(\Lambda_1 \Lambda_2)`. [#two_ways_of_representation]_ There is a priori no reason to use the same representation :math:`D` for both :math:`\psi^{\pm}_{\ell}`, but this turns out to be possible just by calculation. Moreover, the representation :math:`D` is not assumed to be irreducible. Indeed, as we'll see later, it generally decomposes into blocks fixed by further labels.

Then we can try to construct :math:`\Hscr(x)` by a formula similar to :eq:`eq_general_expansion_of_hamiltonian` as follows

.. math::
	:label: eq_construct_interaction_density_by_fields

	\Hscr(x) = \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \psi_{\ell'_1}^-(x) \cdots \psi_{\ell'_N}^-(x) \psi_{\ell_1}^+(x) \cdots \psi_{\ell_M}^+(x)

It follows from :eq:`eq_conjugate_annihilation_and_creation_field` that the interaction density defined by :eq:`eq_construct_interaction_density_by_fields` is a scalar if the following holds

.. math::
	:label: eq_coefficient_g_transformation_law

	g_{\bar{\ell}'_1, \cdots, \bar{\ell}'_N;~\bar{\ell}_1, \cdots, \bar{\ell}_M}
		& = \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} D_{\ell'_1 \bar{\ell}'_1}(\Lambda^{-1}) \cdots D_{\ell'_N \bar{\ell}'_N}(\Lambda^{-1}) \\
		&\qquad \times D_{\ell_1 \bar{\ell}_1}(\Lambda^{-1}) \cdots D_{\ell_M \bar{\ell}_M}(\Lambda^{-1}) g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M}

The solution to the last problem relies on a classification of the representations of the Lorentz group, which has been discussed (in terms of the little group representations) in :ref:`sec_one_particle_states`, and shall be dealt with at a later point. The main goal of this section is to pin down the conditions :math:`u_{\ell}` and :math:`v_{\ell}` must satisfy so that they can be explicitly solved in the following sections.

For this section, we'll focus on the massive case. Recall the Lorentz transformation laws for the creation and annihilation of massive particles :eq:`eq_lorentz_transformation_formula_for_creation_operator`, :eq:`eq_lorentz_transformation_formula_for_annihilation_operator` as follows

.. math::
	:label: eq_lorentz_transformation_formula_for_annihilation_and_creation_operator_revisited

	U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) &= e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \\
	U_0(\Lambda, b) a^{\dagger}(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) &= e^{-\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n) \ast}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a^{\dagger}(\pbf_{\Lambda}, \sigma', n)

where we've used the fact that :math:`D` is unitary in the sense that :math:`D^{\dagger} = D^{-1}` to invert :math:`W(\Lambda, p)` (and flip the indexes) for later convenience -- mostly because of the use of :math:`\Lambda^{-1}` in :eq:`eq_conjugate_annihilation_and_creation_field`.

Using :eq:`eq_lorentz_transformation_formula_for_annihilation_and_creation_operator_revisited` and :eq:`eq_defn_annihilation_and_creation_field`, we can compute the left-hand-side of :eq:`eq_conjugate_annihilation_and_creation_field` as follows

.. math::

	& U_0(\Lambda, b) \psi_{\ell}^+(x) U_0^{-1}(\Lambda, b) \\
		&\quad = \sum_{\sigma, n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) \\
		&\quad = \sum_{\sigma, \sigma', n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n)}_{\sigma \sigma'} (W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \\
		&\quad = \sum_{\sigma, \sigma', n} \int d^3 (\Lambda p)~u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{p_0}{(\Lambda p)_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \\
		&\quad = \sum_{\sigma', n} \int d^3 (\Lambda p) \blue{\sum_{\sigma} u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{p_0}{(\Lambda p)_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p))} a(\pbf_{\Lambda}, \sigma', n)

where in the last equality we've used that the fact from :eq:`eq_lorentz_invariant_3_momentum_volume_element` that :math:`d^3 p / p_0` is Lorentz invariant.

Now the right-hand-side of :eq:`eq_conjugate_annihilation_and_creation_field` can be calculated as follows

.. math::

	\begin{align*}
		& \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi^+_{\ell'}(\Lambda x + b) \\
			&\quad = \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \sum_{\sigma, n} \int d^3 p~u_{\ell'}(\Lambda x + b;~\pbf, \sigma, n) a(\pbf, \sigma, n) \\
			&\quad = \sum_{\sigma', n} \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \int d^3 (\Lambda p)~u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) a(\pbf_{\Lambda}, \sigma', n) \\
			&\quad = \sum_{\sigma', n} \int d^3 (\Lambda p)~\blue{\sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n)} a(\pbf_{\Lambda}, \sigma', n)
	\end{align*}

Equating the blue parts the two calculations, and inverting :math:`D^{(j_n)}_{\sigma \sigma'} (W^{-1}(\Lambda, p))` and :math:`D_{\ell \ell'}(\Lambda^{-1})`, we get

.. math::
	:label: eq_annihilation_u_transformation

	\sqrt{\frac{p_0}{(\Lambda p)_0}} e^{\ifrak b \cdot \Lambda p} \sum_{\ell} D_{\ell' \ell}(\Lambda) u_{\ell}(x;~\pbf, \sigma, n)
		= \sum_{\sigma'} u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) D_{\sigma' \sigma}^{(j_n)} (W(\Lambda, p))

A parallel calculation for the creation field in :eq:`eq_conjugate_annihilation_and_creation_field`, which we'll omit, gives the following

.. math::
	:label: eq_creation_v_transformation

	\sqrt{\frac{p_0}{(\Lambda p)_0}} e^{-\ifrak b \cdot \Lambda p} \sum_{\ell} D_{\ell' \ell}(\Lambda) v_{\ell}(x;~\pbf, \sigma, n)
		= \sum_{\sigma'} v_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) D^{(j_n) \ast}_{\sigma' \sigma}(W(\Lambda, p))

The identities :eq:`eq_annihilation_u_transformation` and :eq:`eq_creation_v_transformation` pose the fundamental conditions on :math:`u_{\ell}` and :math:`v_{\ell}`, respectively, which we'll utilize to eventually solve for their solutions. Currently both :math:`u_{\ell}` and :math:`v_{\ell}` depend on :math:`x, \pbf, \sigma` and :math:`n`, and the goal is to use the Lorentz symmetry to reduce the dependencies. This will be carried out in three steps, corresponding to the three types of Lorentz symmetries: translations, boosts, and rotations, as follows.


Translations
	Taking :math:`\Lambda = 1` in :eq:`eq_creation_v_transformation` gives :math:`\exp(\ifrak b \cdot p) u_{\ell}(x;~\pbf, \sigma, n) = u_{\ell}(x + b;~\pbf, \sigma, n)`, which then implies

	.. math:: u_{\ell}(x;~\pbf, \sigma, n) = (2\pi)^{-3/2} e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n)
		:label: eq_redefine_u_after_translation

	and similarly

	.. math:: v_{\ell}(x;~\pbf, \sigma, n) = (2\pi)^{-3/2} e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n)
		:label: eq_redefine_v_after_translation

	where we've slightly abused notations by keeping the names of :math:`u_{\ell}` and :math:`v_{\ell}`, while changing their arguments. Here the seemingly redundant :math:`(2\pi)^{-3/2}` is inserted so that the fields

	.. math::
		:label: eq_annihilation_and_creation_field_simplified_by_translation

		\psi^+_{\ell}(x) &= \sum_{\sigma, n} (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n) a(\pbf, \sigma, n) \\
		\psi^-_{\ell}(x) &= \sum_{\sigma, n} (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n) a^{\dagger}(\pbf, \sigma, n)

	look like the usual Fourier transforms.

	Plugging :eq:`eq_redefine_u_after_translation` and :eq:`eq_redefine_v_after_translation` into :eq:`eq_annihilation_u_transformation` and :eq:`eq_creation_v_transformation`, respectively, they can be simplified as follows

	.. math::
		:label: eq_annihilation_and_creation_transformation_simplified_by_translation

		\sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) u_{\ell}(\pbf, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(\pbf_{\Lambda}, \sigma', n) D^{(j_n)}_{\sigma' \sigma}(W(\Lambda, p)) \\
		\sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) v_{\ell}(\pbf, \sigma, n) &= \sum_{\sigma'} v_{\ell'}(\pbf_{\Lambda}, \sigma', n) D^{(j_n) \ast}_{\sigma' \sigma}(W(\Lambda, p))

	for any homogeneous Lorentz transformation :math:`\Lambda`.

Boosts
	Taking :math:`\pbf = 0` and :math:`\Lambda = L(q)` which takes a particle at rest to one with (arbitrary) momentum :math:`q`, we see, using :eq:`eq_d_repr_of_little_group`, that

	.. math::

		W(\Lambda, p) = L(\Lambda p)^{-1} \Lambda L(p) = L(q)^{-1} L(q) = 1

	In this case :eq:`eq_annihilation_and_creation_transformation_simplified_by_translation` take the following form (with :math:`\qbf` substituted by :math:`\pbf`)

	.. math::
		:label: eq_annihilation_and_creation_transformation_simplified_by_boost

		\sqrt{\frac{m}{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) u_{\ell}(0, \sigma, n) &= u_{\ell'}(\pbf, \sigma, n) \\
		\sqrt{\frac{m}{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) v_{\ell}(0, \sigma, n) &= v_{\ell'}(\pbf, \sigma, n)

	It follows that one can calculate :math:`u_{\ell}(\pbf, \sigma, n)` for any :math:`\pbf` from the special case of :math:`\pbf = 0` given a representation :math:`D`.

Rotations
	Taking :math:`\pbf = 0` and :math:`\Lambda = \Rcal` a :math:`3`-rotation, and recalling from :eq:`eq_little_group_rotation` that :math:`W(\Lambda, p) = \Rcal`, we get special cases of :eq:`eq_annihilation_and_creation_transformation_simplified_by_translation` as follows

	.. math::

		\sum_{\ell} D_{\ell' \ell}(\Rcal) u_{\ell}(0, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(0, \sigma', n) D_{\sigma' \sigma}^{(j_n)}(\Rcal) \\
		\sum_{\ell} D_{\ell' \ell}(\Rcal) v_{\ell}(0, \sigma, n) &= \sum_{\sigma'} v_{\ell'}(0, \sigma', n) D_{\sigma' \sigma}^{(j_n) \ast}(\Rcal)

	Using :eq:`eq_representation_rotation_first_order` we can further reduce it to the first order as follows

	.. math::
		:label: eq_j_intertwines_u_and_v

		\sum_{\ell} \hat{\Jbf}_{\ell' \ell} u_{\ell}(0, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(0, \sigma', n) \Jbf^{(j_n)}_{\sigma' \sigma} \\
		\sum_{\ell} \hat{\Jbf}_{\ell' \ell} v_{\ell}(0, \sigma, n) &= -\sum_{\sigma'} v_{\ell'}(0, \sigma', n) \Jbf^{(j_n) \ast}_{\sigma' \sigma}

	where :math:`\hat{\Jbf}` denotes the angular momentum vector for the representation :math:`D_{\ell' \ell}(\Rcal)`, in analogy with the usual angular momentum :math:`\Jbf^{(\jfrak)}` for :math:`D^{(j)}(\Rcal)`.


The cluster decomposition principle
+++++++++++++++++++++++++++++++++++

Let's verify that the fields defined by :eq:`eq_annihilation_and_creation_field_simplified_by_translation`, when plugged into :eq:`eq_construct_interaction_density_by_fields`, indeed satisfy the cluster decomposition principle as discussed in :ref:`sec_the_cluster_decomposition_principle`. It's really just a straightforward but tedious calculation which we spell out as follows

.. math::

	V(t) &= \int d^3 x~\Hscr(x) \\
		&= \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \int d^3 x~\psi^-_{\ell'_1}(x) \cdots \psi^-_{\ell'_N}(x) \psi^+_{\ell_1}(x) \cdots \psi^+_{\ell_M}(x) \\
		&= \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \sum_{\sigma'_1, \cdots, \sigma'_N} \sum_{n'_1, \cdots, n'_N} \sum_{\sigma_1, \cdots, \sigma_M} \sum_{n_1, \cdots, n_M} (2\pi)^{-3(N+M)/2} \\
		&\quad \times \int d^3 p'_1 \cdots d^3 p'_N d^3 p_1 \cdots d^3 p_M~\exp(\ifrak (E_1 + \cdots E_M - E'_1 - \cdots - E'_N)t) \\
		&\quad \times \blue{\int d^3 x~\exp(\ifrak (\pbf_1 + \cdots \pbf_M - \pbf'_1 - \cdots - \pbf'_N) \cdot \xbf)} \\
		&\quad \times v_{\ell'_1}(\pbf'_1, \sigma'_1, n'_1) \cdots v_{\ell'_N}(\pbf'_N, \sigma'_N, n'_N) u_{\ell_1}(\pbf_1, \sigma_1, n_1) \cdots u_{\ell_M}(\pbf_M, \sigma_M, n_M) \\
		&\quad \times a^{\dagger}(\pbf'_1, \sigma'_1, n'_1) \cdots a^{\dagger}(\pbf'_N, \sigma'_N, n'_N) a(\pbf_1, \sigma_1, n_1) \cdots a(\pbf_M, \sigma_M, n_M) \\
		&= \sum_{N,M=0}^{\infty} \sum_{\sigma'_1, \cdots, \sigma'_N} \sum_{n'_1, \cdots, n'_N} \sum_{\sigma_1, \cdots, \sigma_M} \sum_{n_1, \cdots, n_M} \int d^3 p'_1 \cdots d^3 p'_N d^3 p_1 \cdots d^3 p_M \\
		&\quad \times (2\pi)^{3 - 3N/2 - 3M/2} \exp(\ifrak (E_1 + \cdots E_M - E'_1 - \cdots - E'_N)t) \\
		&\quad \times a^{\dagger}(\pbf'_1, \sigma'_1, n'_1) \cdots a^{\dagger}(\pbf'_N, \sigma'_N, n'_N) a(\pbf_1, \sigma_1, n_1) \cdots a(\pbf_M, \sigma_M, n_M) \\
		&\quad \times \blue{\delta^3(\pbf_1 + \cdots + \pbf_M - \pbf'_1 - \cdots - \pbf'_N)} \\
		&\quad \times \Big( \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} v_{\ell'_1}(\pbf'_1, \sigma'_1, n'_1) \cdots v_{\ell'_N}(\pbf'_N, \sigma'_N, n'_N) \phantom{)} \\
		&\qquad \phantom{(} \times u_{\ell_1}(\pbf_1, \sigma_1, n_1) \cdots u_{\ell_M}(\pbf_M, \sigma_M, n_M) \Big)

Besides re-ordering the terms, the only actual calculation is highlighted in the two blue terms, where the second one is the integral of the first. One can compare this calculation with :eq:`eq_general_expansion_of_hamiltonian` and see that the cluster decomposition principle is indeed satisfied because there is a unique momentum conservation delta function in each coefficient, as long as :math:`g, u, v` are reasonably smooth, i.e., it's ok to have poles and/or branching singularities but no delta functions.


.. _sec_causality_and_antiparticles:

Causality and antiparticles
+++++++++++++++++++++++++++

We now turn to the other crucial condition on the Hamiltonian, namely, the causality condition :eq:`eq_h_commutativity_for_space_like_separations`. Given the general formula :eq:`eq_construct_interaction_density_by_fields` of the interaction density, we are forced to require that :math:`[\psi^+_{\ell}(x), \psi^-_{\ell'}(y)] = 0` whenever :math:`x - y` is space-like. However, according to :eq:`eq_annihilation_and_creation_field_simplified_by_translation`, we have

.. math::

	& [\psi^+_{\ell}(x), \psi^-_{\ell'}(y)]_{\pm} \\
		&\quad = \sum_{\sigma, n} (2\pi)^{-3} \int d^3 p~d^3 p'~e^{\ifrak (p \cdot x - p' \cdot y)} u_{\ell}(\pbf, \sigma, n) v_{\ell'}(\pbf', \sigma, n) [a(\pbf, \sigma, n), a^{\dagger}(\pbf', \sigma, n)]_{\pm} \\
		&\quad = \sum_{\sigma, n} (2\pi)^{-3} \int d^3 p~e^{\ifrak p \cdot (x - y)} u_{\ell}(\pbf, \sigma, n) v_{\ell'}(\pbf, \sigma, n)

where the sign :math:`\pm` is positive if the field is fermionic, and negative otherwise. This quantity is not necessarily vanishing even if :math:`x - y` is space-like.

Therefore in order to construct :math:`\Hscr` in the form of :eq:`eq_construct_interaction_density_by_fields` that satisfies :eq:`eq_h_commutativity_for_space_like_separations`, we must not just use :math:`\psi^{\pm}(x)` as the building blocks. It turns out that one may consider a linear combination of the two as follows

.. math:: \psi_{\ell}(x) \coloneqq \kappa_{\ell} \psi^+_{\ell}(x) + \lambda_{\ell} \psi^-_{\ell}(x)
	:label: eq_defn_psi_field

as well as its adjoint :math:`\psi^{\dagger}_{\ell}(x)`, and hope that they satisfy

.. math:: [\psi_{\ell}(x), \psi_{\ell}(y)]_{\pm} = [\psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y)]_{\pm} = 0
	:label: eq_space_like_commutativity_for_combined_field

whenever :math:`x-y` is space-like, and replace :math:`\psi^{\pm}_{\ell}(x)` with :math:`\psi_{\ell}(x), \psi^{\dagger}_{\ell}(x)` in :eq:`eq_construct_interaction_density_by_fields`. Under these assumptions, we can then construct the interaction density :math:`\Hscr` as a polynomial in :math:`\psi_{\ell}(x), \psi_{\ell}^{\dagger}(x)` with an even number of fermionic fields (so that the sign in :eq:`eq_space_like_commutativity_for_combined_field` is negative).

There remains, however, one issue with field like :eq:`eq_space_like_commutativity_for_combined_field` that mixes creation and annihilation fields. Namely, special conditions must hold in order for such fields to play well with conserved quantum numbers. To be more specific, let :math:`Q` be a conserved quantum number, e.g., the electric charge. Then the following hold

.. math::
	:label: eq_charge_of_annihilation_and_creation_operator

	[Q, a(\pbf, \sigma, n)] &= -q(n) a(\pbf, \sigma, n) \\
	[Q, a^{\dagger}(\pbf, \sigma, n)] &= q(n) a^{\dagger}(\pbf, \sigma, n)

where :math:`q(n)` denotes the quantum number of the particle species :math:`n`. These identities can be verified by applying both sides to :math:`\Psi_{\pbf, \sigma, n}` and :math:`\Psi_{\VAC}`, respectively.

Now in order for :math:`Q` to commute with :math:`\Hscr`, which is constructed as a polynomial of :math:`\psi_{\ell}(x)` and :math:`\psi_{\ell}^{\dagger}(x)`, we better have

.. math:: [Q, \psi_{\ell}(x)] = -q_{\ell} \psi_{\ell}(x)
	:label: eq_charge_of_psi_field

so that each monomial (with coefficient neglected) :math:`\psi^{\dagger}_{\ell'_1}(x) \cdots \psi^{\dagger}_{\ell'_M}(x) \psi_{\ell_1}(x) \cdots \psi_{\ell_N}(x)` in :math:`\Hscr` will commute with :math:`Q` if

.. math::
	:nowrap:

	\begin{equation*}
		q_{\ell_1} + \cdots + q_{\ell_N} = q_{\ell'_1} + \cdots + q_{\ell'_M}
	\end{equation*}

Note that the negative sign in :eq:`eq_charge_of_psi_field` is a formal analogy to :eq:`eq_charge_of_annihilation_and_creation_operator`, where we think of :math:`\psi_{\ell}(x)` as an annihilation field even though it's really not. Since :math:`\psi_{\ell}(x)` is a linear combination of :math:`\psi^+_{\ell}(x)` and :math:`\psi^-_{\ell}(x)`, which in turn are superpositions of annihilation and creation operators, respectively, it follows from :eq:`eq_charge_of_annihilation_and_creation_operator` that in order for :eq:`eq_charge_of_psi_field` to hold, the following conditions must be satisfied

1. all particles annihilated by :math:`\psi^+_{\ell}(x)` must have the same charge :math:`q(n) = q_{\ell}`,
2. all particles created by :math:`\psi^-_{\ell}(x)` must have the same charge :math:`q(n) = -q_{\ell}`, and
3. for any particle of species :math:`n`, which is annihilated by :math:`\psi^+_{\ell}(x)`, there exists a particle of species :math:`\bar{n}`, which is created by :math:`\psi^-_{\ell}(x)`, such that :math:`q(n) = -q(\bar{n})`.

The particles of species :math:`n` and :math:`\bar{n}` are called *antiparticles* of each other -- they are exactly the same except for the charges which are opposite. It is the last condition that demands the existence of particle-antiparticle pairs so that one can formulate a consistent (relativistic) quantum field theory.

.. dropdown:: The Klein-Gordon equation
	:icon: unlock
	:animate: fade-in-slide-down

	It follows from the definition :eq:`eq_defn_psi_field`, together with :eq:`eq_annihilation_and_creation_field_simplified_by_translation`, that the field :math:`\psi_{\ell}(x)` satisfies the following so-called `Klein-Gordon equation <https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation>`_

	.. math:: \left( \square - m^2 \right) \psi_{\ell}(x) = 0
		:label: eq_klein_gordon

	where :math:`\square \coloneqq \eta^{\mu \nu} \p_{\mu} \p_{\nu}` is the `d'Alembert operator <https://en.wikipedia.org/wiki/D%27Alembert_operator>`_ and :math:`m` is the (definite) mass of the field. This equation is traditionally one of the starting points of quantum field theory, from which creation/annihilation operators can be derived through the so-called canonical quantization formalism. However, we've derived the equation here from the other way around, namely, the creation/annihilation operators, which in turn come from the first principles of quantum mechanics and Lorentz symmetry.


.. _sec_scalar_field:

Scalar Fields
-------------

We'll start, as always, with the simplest case of scalar fields, namely, when :math:`\psi^+(x) = \psi^+_{\ell}(x)` and :math:`\psi^-(x) = \psi^-_{\ell}(x)` are scalar functions. We argue first that such fields can only create/annihilate spinless particles. Indeed, since :math:`\hat{\Jbf}` necessarily vanishes, it follows from :eq:`eq_j_intertwines_u_and_v` that :math:`u` and :math:`v` may be nonzero if and only if :math:`j_n = 0`. If we, for the moment, are concerned with just one particle species, then we can write :math:`u(\pbf, \sigma, n) = u(\pbf)` and :math:`v(\pbf, \sigma, n) = v(\pbf)`. Lastly, we note that since :math:`D = 1` in this case, :eq:`eq_annihilation_and_creation_transformation_simplified_by_translation` become

.. math::

	\sqrt{p_0}~u(\pbf) &= \sqrt{(\Lambda p)_0}~u(\pbf_{\Lambda}) \\
	\sqrt{p_0}~v(\pbf) &= \sqrt{(\Lambda p)_0}~v(\pbf_{\Lambda})

It follows that

.. math:: u(\pbf)	= v(\pbf) = (2p_0)^{-1/2}
	:label: eq_scalar_u_and_v

where the factor :math:`2` is just conventional. In particular :math:`u(0) = v(0) = (2m)^{-1/2}`.

Plugging :eq:`eq_scalar_u_and_v` into :eq:`eq_annihilation_and_creation_field_simplified_by_translation`, we get

.. math::
	:label: eq_scalar_field_psi_plus_and_minus_are_adjoints

	\psi^+(x) &= \int d^3 p~(2\pi)^{-3/2} e^{\ifrak p \cdot x} (2p_0)^{-1/2} a(\pbf) \\
	\psi^-(x) &= \int d^3 p~(2\pi)^{-3/2} e^{-\ifrak p \cdot x} (2p_0)^{-1/2} a^{\dagger}(\pbf) = \psi^{+ \dagger}(x)

In this case the interaction density :math:`\Hscr`, defined by :eq:`eq_construct_interaction_density_by_fields`, may be constructed as any polynomial in :math:`\psi^{\pm}(x)` since :eq:`eq_coefficient_g_transformation_law` holds trivial :math:`D` and scalar :math:`g`.

Next let's consider the causality condition which demands that :math:`\left[ \psi^+(x), \psi^-(y) \right] = 0` whenever :math:`x - y` is space-like. Using the canonical commutation relation :eq:`eq_creation_annihilation_commutator` we calculate

.. math::
	:label: eq_scalar_field_commutator_as_Delta

	\left[ \psi^+(x), \psi^-(y) \right]_{\pm} &= \int d^3 p~d^3 q~(2\pi)^{-3} e^{\ifrak (p \cdot x - q \cdot y)} (4 p_0 q_0)^{-1/2} \left[ a(\pbf), a^{\dagger}(\qbf) \right]_{\pm} \\
		&= \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0}~e^{\ifrak p \cdot (x - y)} \eqqcolon \Delta_+(x - y)

where

.. math:: \Delta_+(x) \coloneqq \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0}~e^{\ifrak p \cdot x}
	:label: eq_defn_Delta_plus

We notice that :math:`\Delta_+(x)` is manifestly (proper orthochronous) Lorentz invariant -- the invariance of the volume element comes from :eq:`eq_lorentz_invariant_3_momentum_volume_element`. It is, however, not in general invariant under transformations like :math:`x \to -x`. But, as we'll see, such invariance holds assuming :math:`x` is space-like.

.. note::
	The plus subscript in :math:`\Delta_+(x)` is there to distinguish it from an anti-symmetrized version :math:`\Delta(x)` to be introduced later.

Now we'll restrict ourselves to the special case of a space-like :math:`x` which, up to a Lorentz transformation, can be assumed to take the form :math:`x = (0, \xbf)` with :math:`|\xbf| > 0`. In this case, we can then calculate :math:`\Delta_+(x)` as follows [#wrong_integration_of_Delta_function]_

.. math::

	\Delta_+(x) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2\sqrt{\pbf^2 + m^2}}~\exp(\ifrak \pbf \cdot \xbf) \\
		&= \frac{4\pi}{(2\pi)^3} \int_0^{\infty} \frac{\pbf^2 d|\pbf|}{2\sqrt{\pbf^2 + m^2}} \int_{S^2} d^2 \hat{\pbf}~\exp(\ifrak |\pbf| |\xbf| \hat{\pbf} \cdot \hat{\xbf}) \\
		&= \frac{1}{2\pi} \int_0^{\infty} \frac{\pbf^2 d|\pbf|}{2\sqrt{\pbf^2 + m^2}} \int_0^{\pi} d\theta~\exp(\ifrak |\pbf| |\xbf| \cos\theta)

The last integral cannot be easily evaluated, at least without some knowledge about special functions. Nonetheless, we observe that :math:`\Delta_+(x) \neq 0`, which means that :math:`\Hscr` cannot be just any polynomial in :math:`\psi^{\pm}(x)`. Moreover, we note that :math:`\Delta_+(x) = \Delta_+(-x)` as promised earlier.

As already mentioned in :eq:`eq_defn_psi_field`, let's try

.. math:: \psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^-(x)
	:label: eq_scalar_field_first_defn_of_psi

Using :eq:`eq_scalar_field_psi_plus_and_minus_are_adjoints` and :eq:`eq_scalar_field_commutator_as_Delta`, we can then try to make :eq:`eq_space_like_commutativity_for_combined_field` hold by the following calculations

.. math::

	\left[ \psi(x), \psi(y) \right]_{\pm}
		&= \kappa\lambda \left(\left[ \psi^+(x), \psi^-(y) \right]_{\pm} + \left[ \psi^-(x), \psi^+(y) \right]_{\pm} \right) \\
		&= \kappa\lambda (1 \pm 1) \Delta(x - y) \\
	\left[ \psi(x), \psi^{\dagger}(y) \right]_{\pm}
		&= \left[ \kappa \psi^+(x) + \lambda \psi^-(x), \kappa^{\ast} \psi^-(y) + \lambda^{\ast} \psi^+(y) \right]_{\pm} \\
		&= |\kappa|^2 \left[ \psi^+(x), \psi^-(y) \right]_{\pm} + |\lambda|^2 \left[ \psi^-(x), \psi^+(y) \right]_{\pm} \\
		&= \left( |\kappa|^2 \pm |\lambda|^2 \right) \Delta(x - y)

We see that :eq:`eq_space_like_commutativity_for_combined_field` holds for scalar fields if the fields are bosonic, i.e., the bottom sign in :math:`\pm` applies, and :math:`|\kappa| = |\lambda|`. By adjust the phase of :math:`a(\pbf)`, we can actually arrange so that :math:`\kappa = \lambda`, in which case we have

.. math::
	:label: eq_scalar_field_psi_fixed_phase

	\psi(x) = \psi^+(x) + \psi^-(x) = \psi^+(x) + \psi^{+ \dagger}(x) = \psi^{\dagger}(x)

.. note::
	Although the arrangement of phase so that :math:`\kappa = \lambda` is a mere convention, it's a convention that needs to be applied to *all* scalar fields appearing in :math:`\Hscr`. Namely, one cannot have both :math:`\psi(x)` as in :eq:`eq_scalar_field_psi_fixed_phase` and another

	.. math:: \psi'(x) = e^{\ifrak \theta} \psi^+(x) + e^{-\ifrak \theta} \psi^{+ \dagger}(x)

	for some :math:`\theta`, because :math:`\psi(x)` won't commute with :math:`\psi'(y)` even if :math:`x - y` is space-like.

Now if the particle created and annihilated by :math:`\psi(x)` carries a (non-vanishing) conserved quantum number :math:`Q`, then by the discussions on the charge conservation from the previous section, a density :math:`\Hscr` made up of :math:`\psi(x)` as defined by :eq:`eq_scalar_field_psi_fixed_phase` will not commute with :math:`Q`. Instead, one must assume the existence of a field :math:`\psi^{+ c}(x)` that creates and annihilates the corresponding antiparticle, in the sense that

.. math::

	\left[ Q, \psi^+(x) \right] &= -q \psi^+(x) \\
	\left[ Q, \psi^{+ c}(x) \right] &= q \psi^{+ c}(x)

Here the supscript :math:`c` stands for charge (conjugation). Now instead of :eq:`eq_scalar_field_first_defn_of_psi`, let's try

.. math:: \psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x)

so that :math:`[Q, \psi(x)] = -q \psi(x)`. We calculate the commutators, assuming the antiparticle is different from the particle, just as before as follows

.. math::

	\left[\psi(x), \psi(y) \right]_{\pm} &= \left[ \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x), \kappa \psi^+(y) + \lambda \psi^{+ c \dagger}(y) \right]_{\pm} = 0 \\
	\left[\psi(x), \psi^{\dagger}(y) \right]_{\pm} &= \left[ \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x), \kappa^{\ast} \psi^{+ \dagger}(y) + \lambda^{\ast} \psi^{+ c}(y) \right]_{\pm} \\
		&= |\kappa|^2 \left[ \psi^+(x), \psi^{+ \dagger}(y) \right]_{\pm} + |\lambda|^2 \left[ \psi^{+ c \dagger}(x), \psi^{+ c}(y) \right]_{\pm} \\
		&= (|\kappa|^2 \pm |\lambda|^2) \Delta(x - y)

where we've assumed, in particular that the particle and its particle share the same mass so that :eq:`eq_scalar_field_commutator_as_Delta` equally applies.

By the same argument as in the case where no quantum number is involved, we see that a scalar field can satisfy the causality condition if it describes a boson. Moreover, by adjusting the phase of :math:`a(\pbf)`, one can arrange so that :math:`\kappa = \lambda` so that

.. math:: \psi(x) = \psi^+(x) + \psi^{+ c \dagger}(x)
	:label: eq_scalar_field_psi_fixed_phase_with_antiparticle

Note that this is compatible with :eq:`eq_scalar_field_psi_fixed_phase` in the case where the particle is its own antiparticle.

Using :eq:`eq_scalar_field_psi_plus_and_minus_are_adjoints`, we can write :math:`\psi(x)` in terms of the creation and annihilation operators as follows

.. math::
	:label: eq_scalar_field_psi_by_creation_and_annihilation_operators

	\psi(x) = \int \frac{d^3 p}{(2\pi)^{3/2} (2p_0)^{1/2}}~\left[ e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{c \dagger}(\pbf) \right]

with the possibility of :math:`a^{c \dagger}(\pbf) = a^{\dagger}(\pbf)` in the case where the created particle is its own antiparticle.

For later use (e.g., the evaluation of Feynman diagrams), we note the following identity which holds for any, and not just space-like, :math:`x` and :math:`y`.

.. math::
	:label: eq_scalar_field_commutator

	\left[ \psi(x), \psi^{\dagger}(y) \right] = \Delta(x - y)

where :math:`\Delta(x)` is defined as follows

.. math::
	:label: eq_defn_Delta

	\Delta(x) \coloneqq \Delta_+(x) - \Delta_+(-x) = \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0} \left( e^{\ifrak p \cdot x} - e^{-\ifrak p \cdot x} \right)


The CPT symmetries
++++++++++++++++++

Let's investigate how a scalar field transforms under spatial inversion :math:`\Pcal`, time inversion :math:`\Tcal`, and charge conjugation :math:`\Ccal`. This follows essentially from :eq:`eq_scalar_field_psi_by_creation_and_annihilation_operators` together with our knowledge about how creation/annihilation operators transform under CPT transformations in :ref:`sec_the_lorentz_and_cpt_transformation_laws`. Recall that we consider the case of massive particles here, leaving the massless case to a later section.

We start with the spatial inversion :math:`\Pcal` by recalling the following transformation rules

.. math::

	U(\Pcal) a(\pbf) U^{-1}(\Pcal) &= \eta^{\ast} a(-\pbf) \\
	U(\Pcal) a^{c \dagger}(\pbf) U^{-1}(\Pcal) &= \eta^c a^{c \dagger}(-\pbf)

where :math:`\eta` and :math:`\eta^c` are the intrinsic parities of the particle and antiparticle, respectively. In order for the scalar field :eq:`eq_scalar_field_psi_by_creation_and_annihilation_operators` to transform nicely with :math:`\Pcal`, one must have :math:`\eta^{\ast} = \eta^c` (or :math:`\eta^{\ast} = \eta` in the case where the particle is its own antiparticle). As a result, we have

.. math::
	:label: eq_scalar_field_spatial_inversion_transformation_law

	U(\Pcal) \psi(x) U^{-1}(\Pcal) = \eta^{\ast} \psi(\Pcal x)

Next let's consider the time inversion :math:`\Tcal`. We recall the transformation rules as follows

.. math::

	U(\Tcal) a(\pbf) U^{-1}(\Tcal) &= \zeta^{\ast} a(-\pbf) \\
	U(\Tcal) a^{c \dagger}(\pbf) U^{-1}(\Tcal) &= \zeta^c a^{c \dagger}(-\pbf)

Similar to the case of spatial inversions, in order for :math:`\psi(x)` to transform nicely with :math:`U(\Tcal)`, one must have :math:`\zeta^{\ast} = \zeta^c`. Moreover, since :math:`U(\Tcal)` is anti-unitary, we have

.. math:: U(\Tcal) \psi(x) U^{-1}(\Tcal) = \zeta^{\ast} \psi(-\Tcal x)

Finally let's consider the charge conjugation :math:`\Ccal` with the following transformation laws

.. math::

	U(\Ccal) a(\pbf) U^{-1}(\Ccal) &= \xi^{\ast} a^c(\pbf) \\
	U(\Ccal) a^{c \dagger}(\pbf) U^{-1}(\Ccal) &= \xi^c a^{\dagger}(\pbf)

As before, we must have :math:`\xi^{\ast} = \xi^c` and therefore

.. math:: U(\Ccal) \psi(x) U^{-1}(\Ccal) = \xi^{\ast} \psi^{\dagger}(x)


Vector Fields
-------------

The next simplest scenario after scalar field is vector field, where the representation :math:`D(\Lambda) = \Lambda`. Once again, let's consider particles of one species so that we can drop the :math:`n` label from, for example, :math:`a(\pbf, \sigma, n)`. In this case, we can rewrite :eq:`eq_annihilation_and_creation_field_simplified_by_translation` as follows

.. math::
	:label: eq_vector_field_psi

	\psi^+_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} u_{\mu}(\pbf, \sigma) a(\pbf, \sigma) \\
	\psi^-_{\nu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} v_{\nu}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma)

where :math:`\mu, \nu` are the :math:`4`-indexes. Moreover, the boost transformation formulae :eq:`eq_annihilation_and_creation_transformation_simplified_by_boost` take the following form

.. math::
	:label: eq_vector_field_boost_u_and_v

	u_{\mu}(\pbf, \sigma) &= (m / p_0)^{1/2} {L(p)_{\mu}}^{\nu} u_{\nu}(0, \sigma) \\
	v_{\mu}(\pbf, \sigma) &= (m / p_0)^{1/2} {L(p)_{\mu}}^{\nu} v_{\nu}(0, \sigma)

Finally the (linearized) rotation transformation formulae :eq:`eq_j_intertwines_u_and_v` take the following form

.. math::
	:label: eq_vector_field_angular_momentum_intertwines_u_and_v

	\sum_{\sigma'} u_{\mu}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_{\nu} \hat{\Jbf}_{\mu \nu} u_{\nu}(0, \sigma) \\
	-\sum_{\sigma'} v_{\mu}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_{\nu} \hat{\Jbf}_{\mu \nu} v_{\nu}(0, \sigma)

where :math:`\hat{\Jbf}` is the angular momentum vector associated with the (tautological) representation :math:`\Lambda`. It follows from :eq:`eq_expansion_of_Lambda` and :eq:`eq_u_lorentz_expansion` that

.. math::
	:label: eq_vector_field_j_intertwines_u_and_v

	\left( \hat{\Jbf}_k \right)_{00} = \left( \hat{\Jbf}_k \right)_{0i} = \left( \hat{\Jbf}_k \right)_{i0} &= 0 \\
	\left( \hat{\Jbf}_k \right)_{ij} &= -\ifrak \epsilon_{ijk}

where :math:`\{i,j,k\} = \{1,2,3\}`. From this one can then calculate :math:`\hat{\Jbf}^2` as follows

.. math::

	\left( \hat{\Jbf}^2 \right)_{00} &= \left( \hat{\Jbf}^2 \right)_{0i} = \left( \hat{\Jbf}^2 \right)_{i0} = 0 \\
	\left( \hat{\Jbf}^2 \right)_{ij} &= \sum_{k,m=1}^3 \left( \hat{\Jbf}_k \right)_{im} \left( \hat{\Jbf}_k \right)_{mj} = \sum_{k,m=1}^3 -\epsilon_{imk} \epsilon_{mjk} = 2\delta_{ij}

It follows then from :eq:`eq_vector_field_angular_momentum_intertwines_u_and_v` that

.. math::
	:label: eq_vector_field_u_and_v_multiply_j_sqaured

	\sum_{\sigma'} u_0(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma}
		&= \sum_{\nu} \left( \hat{\Jbf}^2 \right)_{0 \nu} u_{\nu}(0, \sigma)
		= 0 \\
	\sum_{\sigma'} u_i(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma}
		&= \sum_{\nu} \left( \hat{\Jbf}^2 \right)_{i \nu} u_{\nu}(0, \sigma)
		= \sum_j 2 \delta_{i j} u_j(0, \sigma)
		= 2 u_i(0, \sigma) \\
	\sum_{\sigma'} v_0(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma}
		&= 0 \\
	\sum_{\sigma'} v_i(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma}
		&=2v_i(0, \sigma)

where we've worked out the details of the calculations for :math:`u`, but not :math:`v` because they are essentially the same.

Now recall from :eq:`eq_angular_momentum_squared_eigenvalue` that

.. math:: \left( \Jbf^{(\jfrak)} \right)^2_{\sigma \sigma'} = \jfrak (\jfrak + 1) \delta_{\sigma \sigma'}

It follows that in order for :eq:`eq_vector_field_u_and_v_multiply_j_sqaured` to have nonzero solutions, one must have either :math:`\jfrak = 0`, in which case only the time-components :math:`u_0(0)` and :math:`v_0(0)` may be nonzero, where we've also suppressed :math:`\sigma` because spin vanishes, or :math:`\jfrak = 1`, in which case only the space-components :math:`u_i(0, \sigma)` and :math:`v_i(0, \sigma)` may be nonzero. These two cases are discussed in more details as follows.

.. _sec_spin_zero_vector_field:

Spin-:math:`0` vector fields
++++++++++++++++++++++++++++

In this case :math:`\jfrak = 0`. For reasons that will become clear momentarily, let's fix the constants :math:`u_0(0), v_0(0)` as follows

.. math::
	:nowrap:

	\begin{align*}
		u_0(0) &= \ifrak (m / 2)^{1/2} \\
		v_0(0) &= -\ifrak (m / 2)^{1/2}
	\end{align*}

It follows from :eq:`eq_vector_field_boost_u_and_v` (see also :eq:`eq_L_transformation_for_massive`) that

.. math::

	u_{\mu}(\pbf) &= (m / p_0)^{1/2} {L(p)_{\mu}}^0 u_0(0) \\
		&= (m / p_0)^{1/2} (p_{\mu} / m) \ifrak (m / 2)^{1/2} \\
		&= \ifrak p_{\mu} (2p_0)^{-1/2} \\
	v_{\mu}(\pbf) &= -\ifrak p_{\mu} (2p_0)^{-1/2}

where we once again have omitted the details of the calculation of :math:`v` because it's similar to that of :math:`u`. Plugging into :eq:`eq_vector_field_psi`, we see that the field components take the following form

.. math::

	\psi^+_{\mu}(x) &= (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} \ifrak p_{\mu} (2p_0)^{-1/2} a(\pbf) \\
	\psi^-_{\mu}(x) &= (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} (-\ifrak p_{\mu}) (2p_0)^{-1/2} a^{\dagger}(\pbf)

Comparing these with :eq:`eq_scalar_field_psi_plus_and_minus_are_adjoints`, and thanks to the choices of :math:`u_0(0)` and :math:`v_0(0)` above, we see that

.. math:: \psi^{\pm}_{\mu}(x) = \p_{\mu} \psi^{\pm}(x)

It follows that in fact a spinless vector field defined by :math:`\psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^-_{\mu}(x)` as usual is nothing but the gradient vector field of a (spinless) scalar field. Hence we get nothing new from spinless vector fields.

.. _sec_spin_1_vector_fields:

Spin-:math:`1` vector fields
++++++++++++++++++++++++++++

In this case :math:`\jfrak = 1`. We start with the states whose spin :math:`z`-component vanishes, i.e., :math:`u_{\mu}(0,0)` and :math:`v_{\mu}(0,0)`. First we claim that they are both in the :math:`z`-direction, i.e., :math:`u_{\mu}(0,0) = v_{\mu}(0,0) = 0` unless :math:`\mu=3`. Indeed, taking the :math:`z`-components of both sides of :eq:`eq_vector_field_angular_momentum_intertwines_u_and_v` and recalling that :math:`\left( J_3^{(1)} \right)_{0 \sigma} = 0`, we have for :math:`\mu = 1`

.. math:: 0 = \sum_{\nu} \left( \hat{\Jbf}_3 \right)_{1 \nu} u_{\nu}(0, 0) = -\ifrak u_2(0, 0) \implies u_2(0, 0) = 0

and for :math:`\mu = 2`

.. math:: 0 = \sum_{\nu} \left( \hat{\Jbf}_3 \right)_{2 \nu} u_{\nu}(0, 0) = \ifrak u_1(0, 0) \implies u_1(0, 0) = 0

These, together with the fact that :math:`u_0(0, 0) = 0` for :math:`\jfrak = 1`, imply that only :math:`u_3(0, 0)` can be nonzero. The same conclusion can also be drawn for :math:`v_3(0, 0)`. Therefore up to a normalization factor, we can write

.. math:: u_{\mu}(0, 0) = v_{\mu}(0, 0) = (2m)^{-1/2} \begin{bmatrix*}[r] 0 \\ 0 \\ 0 \\ 1 \end{bmatrix*}
	:label: eq_vector_field_uv_spin_z_0

Now to calculate :math:`u` and :math:`v` for the other spin :math:`z`-components, we'll try to use :eq:`eq_representation_rotation_first_order` as follows. First, according to :eq:`eq_vector_field_angular_momentum_intertwines_u_and_v` we have the following general equality

.. math::

	\sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{\mu \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{\mu \nu} \right) u_{\nu}(0, \sigma)
		&= \sum_{\sigma'} u_{\mu}(0, \sigma') \left( J^{(1)}_1 + \ifrak J^{(1)}_2 \right)_{\sigma \sigma'} \\
		&= \sum_{\sigma'} u_{\mu}(0, \sigma') \delta_{\sigma+1, \sigma'} \sqrt{(1 - \sigma)(2 + \sigma)}

Then, letting :math:`\sigma=0` and :math:`\mu=1`, we have

.. math::

	\sqrt{2}~u_1(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{1 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{1 \nu} \right) u_{\nu}(0, 0) = \ifrak \left( \hat{\Jbf}_2 \right)_{13} u_3(0, 0) = -(2m)^{-1/2}

Then, changing to :math:`\mu=2`, we have

.. math::

	\sqrt{2}~u_2(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{2 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{2 \nu} \right) u_{\nu}(0, 0) = (\hat{\Jbf}_1)_{23} u_3(0, 0) = -\ifrak (2m)^{-1/2}

Finally taking :math:`\mu=3`, we have

.. math::

	\sqrt{2}~u_3(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{3 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{3 \nu} \right) u_{\nu}(0, 0) = \left( \hat{\Jbf}_1 \right)_{32} u_2(0, 0) + \ifrak \left( \hat{\Jbf}_2 \right)_{31} u_1(0, 0) = 0

Putting these all together, we have calculated :math:`u_{\mu}(0, 1)` as follows

.. math::

	u_{\mu}(0, 1) = -\frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*}

Calculations for :math:`\sigma = -1` as well as for :math:`v` are similar and hence omitted. The results are listed for future reference as follows

.. math::
	:label: eq_vector_field_u_and_v_stationary

	\begin{alignat*}{2}
		u_{\mu}(0, 1) &= -v_{\mu}(0, -1) &&= -\frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*} \\
		u_{\mu}(0, -1) &= -v_{\mu}(0, 1) &&= \frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ -\ifrak \\ 0 \end{bmatrix*}
	\end{alignat*}

Applying the boosting formulae :eq:`eq_vector_field_boost_u_and_v` to :eq:`eq_vector_field_uv_spin_z_0` and :eq:`eq_vector_field_u_and_v_stationary`, we obtain the formulae for :math:`u` and :math:`v` with arbitrary momentum as follows

.. math::
	:label: eq_vector_field_defn_e_vector_at_p

	u_{\mu}(\pbf, \sigma) = v_{\mu}^{\ast}(\pbf, \sigma) = (2p_0)^{-1/2} {L(p)_{\mu}}^{\nu} e_{\nu}(0, \sigma) \eqqcolon (2p_0)^{-1/2} e_{\mu}(\pbf, \sigma)

where

.. math::
	:label: eq_vector_field_defn_e_vector_at_rest

	e_{\mu}(0, 0) = \begin{bmatrix*}[r] 0 \\ 0 \\ 0 \\ 1 \end{bmatrix*}, \quad \
	e_{\mu}(0, 1) = -\frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*}, \quad \
	e_{\mu}(0, -1) = \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ -\ifrak \\ 0 \end{bmatrix*}

Now we can rewrite the general :eq:`eq_vector_field_psi` more specifically as follows

.. math::
	:label: eq_vector_field_psi_minus_adjoint_to_plus

	\psi^+_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}}~\exp(\ifrak p \cdot x) e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) \\
	\psi^-_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}}~\exp(-\ifrak p \cdot x) e_{\mu}^{\ast}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma) = \psi^{+ \dagger}_{\mu}(x)

Similar to the calculation :eq:`eq_scalar_field_commutator_as_Delta` for scalar field, the (anti-)commutator can be calculated as follows

.. math::
	:label: eq_vector_field_commutator_by_Pi

	\left[ \psi^+_{\mu}(x), \psi^-_{\nu}(y) \right]_{\pm} = \int \frac{d^3 p}{(2\pi)^3 2p_0}~\exp(\ifrak p \cdot (x - y)) \Pi_{\mu \nu}(\pbf)

where

.. math::
	:label: eq_vector_field_Pi_matrix

	\Pi_{\mu \nu}(\pbf) \coloneqq \sum_{\sigma} e_{\mu}(\pbf, \sigma) e^{\ast}_{\nu}(\pbf, \sigma)

To better understand the quantity :math:`\Pi_{\mu \nu}(\pbf)`, let's first evaluate it at :math:`\pbf = 0` as follows

.. math::

	\Pi_{\mu \nu}(0) = \begin{bmatrix*}[r]
			0 & 0 & 0 & 0 \\
			0 & 0 & 0 & 0 \\
			0 & 0 & 0 & 0 \\
			0 & 0 & 0 & 1
		\end{bmatrix*} + \frac{1}{2} \begin{bmatrix*}[r]
			0 & 0 & 0 & 0 \\
			0 & 1 & -\ifrak & 0 \\
			0 & \ifrak & 1 & 0 \\
			0 & 0 & 0 & 0
		\end{bmatrix*} + \frac{1}{2} \begin{bmatrix*}[r]
			0 & 0 & 0 & 0 \\
			0 & 1 & \ifrak & 0 \\
			0 & -\ifrak & 1 & 0 \\
			0 & 0 & 0 & 0
		\end{bmatrix*} = \begin{bmatrix*}[r]
			0 & 0 & 0 & 0 \\
			0 & 1 & 0 & 0 \\
			0 & 0 & 1 & 0 \\
			0 & 0 & 0 & 1
		\end{bmatrix*}

which is nothing but the projection to the spatial :math:`3`-space, or phrased more invariantly, the orthogonal complement of the time direction. Considering the definition :math:`e_{\mu}(\pbf, \sigma) \coloneqq {L(p)_{\mu}}^{\nu} e_{\nu}(0, \sigma)` as in :eq:`eq_vector_field_defn_e_vector_at_p`, we see that the general :math:`\Pi_{\mu \nu}(\pbf)` is really just a projection to the orthogonal complement of :math:`p`, and therefore can be written as

.. math:: \Pi_{\mu \nu}(\pbf) = \eta_{\mu \nu} + \frac{p_{\mu} p_{\nu}}{m^2}
	:label: eq_vector_field_defn_Pi

because of the mass-shell condition :math:`p^2 + m^2 = 0`.

In light of :eq:`eq_scalar_field_commutator_as_Delta`, we can rewrite :eq:`eq_vector_field_commutator_by_Pi` as follows

.. math::
	:label: eq_vector_field_commutator_by_Delta

	\left[ \psi^+_{\mu}(x), \psi^-_{\nu}(y) \right]_{\pm} = \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x - y)

where :math:`\Delta_+(x - y)` is defined by :eq:`eq_defn_Delta_plus`. As in the case of scalar fields, this (anti-)commutator doesn't vanish even for space-like :math:`x - y`. Nonetheless, it's still an even function for space-like separations. The trick, as usual, is to consider a linear combination of :math:`\psi^+_{\mu}(x)` and :math:`\psi^-_{\mu}(x)` as follows

.. math:: \psi_{\mu}(x) \coloneqq \kappa \psi^+_{\mu}(x) + \lambda \psi^-_{\mu}(x)

Now for space-separated :math:`x` and :math:`y`, we can calculate using :eq:`eq_vector_field_commutator_by_Delta` and :eq:`eq_vector_field_psi_minus_adjoint_to_plus` as follows

.. math::

	\left[ \psi_{\mu}(x), \psi_{\nu}(y) \right]_{\pm} &= \kappa\lambda(1 \pm 1) \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x-y) \\
	\left[ \psi_{\mu}(x), \psi^{\dagger}_{\nu}(y) \right]_{\pm} &= (|\kappa|^2 \pm |\lambda|^2) \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x-y)

For them to vanishes, we see that first of all, we must adopt the top sign, i.e., take the commutator, or in other words, the vector field of spin :math:`1` must be bosonic. In addition, we must have :math:`|\kappa| = |\lambda|`. In fact, by adjusting the phase of the creation/annihilation operators, we can arrange so that :math:`\kappa = \lambda = 1`. To summarize, we can write a general vector field in the following form

.. math:: \psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^-_{\mu}(x) = \psi^+_{\mu}(x) + \psi^{+ \dagger}_{\mu}(x)
	:label: eq_vector_field_psi_fixed_phase

just like :eq:`eq_scalar_field_psi_fixed_phase`. It's also obvious that :math:`\psi_{\mu}(x)` is Hermitian.

Now if the vector field carries a nonzero (conserved) quantum charge, then one must adjust :eq:`eq_vector_field_psi_fixed_phase` as follows

.. math:: \psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^{+ c \dagger}_{\mu}(x)

in analogy with :eq:`eq_scalar_field_psi_fixed_phase_with_antiparticle` for scalar fields. Finally, we can express the vector field in terms of creation and annihilation operators as follows

.. math::
	:label: eq_vector_field_psi_by_creation_and_annihilation_operators

	\psi_{\mu}(x) = \sum_{\sigma} \int \frac{d^3 p}{(2\pi)^{3/2} (2p_0)^{1/2}}~\left[ e^{\ifrak p \cdot x} e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) \
		+ e^{-\ifrak p \cdot x} e_{\mu}^{\ast}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right]

in analogy with :eq:`eq_scalar_field_psi_by_creation_and_annihilation_operators` for scalar fields. Finally, let's calculate the commutator (for general :math:`x` and :math:`y`) for later use as follows

.. math::
	:label: eq_vector_field_commutator

	\left[ \psi_{\mu}(x), \psi^{\dagger}_{\nu}(y) \right] = \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta(x-y)

where :math:`\Delta(x-y)` as defined by :eq:`eq_defn_Delta`.

So far, besides the introduction of the vectors :math:`e_{\mu}(\pbf, \sigma)` in :eq:`eq_vector_field_defn_e_vector_at_p` and :eq:`eq_vector_field_defn_e_vector_at_rest`, the discussion on vector fields looks very much like scalar fields. A key difference, however, stems from the following observation

.. math:: e^{\mu}(\pbf, \sigma) p_{\mu} = 0
	:label: eq_vector_field_spinor_orthogonal_to_momentum

which, in turn, implies that

.. math:: \p_{\mu} \psi^{\mu}(x) = 0
	:label: eq_vector_field_gauge_fixing_condition

This condition turns out to be coincide with a so-called "gauge fixing" condition for spin-:math:`1` photons in quantum electrodynamics. However, it's known that photons are massless particles. Therefore we may wonder if a vanishing mass limit :math:`m \to 0` may be applied. Now the simplest way to construct a (scalar) interaction density :math:`\Hscr(x)` using :math:`\psi_{\mu}(x)` is

.. math:: \Hscr(x) = J^{\mu}(x) \psi_{\mu}(x)
	:label: eq_vector_field_j_coupling

where :math:`J^{\mu}(x)` is a :math:`4`-vector current. Suppose we fix the in- and out-states in the interaction. Then according to :eq:`eq_vector_field_psi_by_creation_and_annihilation_operators`, the rate of (anti-)particle emission is proportional to

.. math::

	\sum_{\sigma} \left| \langle J^{\mu} \rangle e_{\mu}^{\ast}(\pbf, \sigma) \right|^2
		= \langle J^{\mu} \rangle \langle J^{\nu} \rangle^{\ast} \Pi_{\mu \nu}(\pbf)
		= \langle J^{\mu} \rangle \langle J^{\nu} \rangle^{\ast} \left( \eta_{\mu \nu} - p_{\mu} p_{\nu} / m^2 \right)

where :math:`\Pi_{\mu \nu}` is evaluated by :eq:`eq_vector_field_defn_Pi`, and :math:`\langle J^{\mu} \rangle` denotes the matrix element of the current between the fixed in- and out-states. Now this rate blows up at :math:`m \to 0` limit unless :math:`p_{\mu} \langle J^{\mu} \rangle = 0`. This last condition can be translated to spacetime coordinates as follows

.. math:: \p_{\mu} J^{\mu}(x) = 0
	:label: eq_vector_field_j_coupling_condition

or in other words :math:`J^{\mu}(x)` is a conserved current.

The CPT symmetries
++++++++++++++++++

Let's start with the spatial inversion. First recall from :ref:`sec_the_lorentz_and_cpt_transformation_laws`

.. math::

	U(\Pcal) a(\pbf, \sigma) U^{-1}(\Pcal) &= \eta^{\ast} a(-\pbf, \sigma) \\
	U(\Pcal) a^{c \dagger}(\pbf, \sigma) U^{-1}(\Pcal) &= \eta^c a^{c \dagger}(-\pbf, \sigma)

It follows that we need to express :math:`e_{\mu}(-\pbf, \sigma)` in terms of :math:`e_{\mu}(\pbf, \sigma)`. To this end, let's calculate

.. math::
	:label: eq_vector_field_revert_momentum_transformation

	e_{\mu}(-\pbf, \sigma)
		&= {L(-\pbf)_{\mu}}^{\nu} e_{\nu}(0, \sigma) \\
		&= {\Pcal_{\mu}}^{\rho} {L(\pbf)_{\rho}}^{\tau} {\Pcal_{\tau}}^{\nu} e_{\nu}(0, \sigma) \\
		&= -{\Pcal_{\mu}}^{\rho} {L(\pbf)_{\rho}}^{\tau} e_{\tau}(0, \sigma) \\
		&= -{\Pcal_{\mu}}^{\rho} e_{\rho}(\pbf, \sigma)

It follows that the spatial inversion transformation law is given as follows

.. math:: U(\Pcal) \psi_{\mu}(x) U^{-1}(\Pcal) = -\eta^{\ast} {\Pcal_{\mu}}^{\nu} \psi_{\nu}(\Pcal x)
	:label: eq_vector_field_spatial_inversion_transformation_law

under the following assumption

.. math:: \eta^c = \eta^{\ast}
	:label: eq_vector_field_eta_assumption

Omitting further details, the transformation laws for time inversion and charge conjugation are given by

.. math::

	U(\Tcal) \psi_{\mu}(x) U^{-1}(\Tcal) &= \zeta^{\ast} {\Pcal_{\mu}}^{\nu} \psi_{\nu}(-\Pcal x) \\
	U(\Ccal) \psi_{\mu}(x) U^{-1}(\Ccal) &= \xi^{\ast} \psi_{\mu}^{\dagger}(x)

under the assumptions

.. math::
	:label: eq_vector_field_zeta_and_xi_assumptions

	\zeta^c = \zeta^{\ast}, \quad \xi^c = \xi^{\ast}

respectively.


Dirac Fields
------------

Here we'll encounter the first nontrivial representation of the (homogeneous orthochronous) Lorentz group, first discovered by P. Dirac in a completely different (and more physical) context. Our treatment here will be purely mathematical, and will serve as a warm-up for the general representation theory.

.. _sec_dirac_representation_and_gamma_matrices:

Dirac representation and gamma matrices
+++++++++++++++++++++++++++++++++++++++

Let :math:`D` be a representation of the Lorentz group in the sense that :math:`D(\Lambda_1) D(\Lambda_2) = D(\Lambda_1 \Lambda_2)`. By the discussion in :ref:`sec_quantum_lorentz_symmetry` and ignoring the translation part, we can write :math:`D(\Lambda)` up to first order as follows

.. math:: D(\Lambda) = 1 + \frac{\ifrak}{2} \omega^{\mu \nu} \Jscr_{\mu \nu} + \cdots
	:label: eq_dirac_field_linearize_representation

where :math:`\Jscr_{\mu \nu} = -\Jscr_{\nu \mu}` are (Hermitian) matrices that, according to :eq:`eq_bracket_j4_j4`, satisfy in addition the following Lie-algebraic condition

.. math::
	:label: eq_bracket_repr_j

	\left[ \Jscr_{\mu \nu}, \Jscr_{\rho \kappa} \right] = \ifrak\left( \eta_{\mu \rho} \Jscr_{\nu \kappa} - \eta_{\nu \rho} \Jscr_{\mu \kappa} + \eta_{\kappa \mu} \Jscr_{\rho \nu} - \eta_{\kappa \nu} \Jscr_{\rho \mu} \right)

Putting it this way, it may appear hopeless to find any solution to the equation above. Surprisingly, there is in fact a systematic way to find *all* solutions to :eq:`eq_bracket_repr_j`, which will be shown in :ref:`sec_general_fields`. The aim of this section, however, is to explain a seemingly unmotivated, but rather ingenious, solution, which also bares a great deal of significance in the quantum theory of electromagnetism.

The trick here is to assume the existence of a set of matrices :math:`\gamma_{\mu}` such that

.. math:: \left\{ \gamma_{\mu}, \gamma_{\nu} \right\} = 2\eta_{\mu \nu}
	:label: eq_dirac_field_clifford_algebra

where the curly bracket denotes the anti-commutator, and is equivalent to the notation :math:`[~,~]_+` used in previous chapters. Here the right-hand-side, written as a number, should be interpreted as a multiple of the identity matrix. Such matrices :math:`\gamma_{\mu}` form a so-called `Clifford algebra <https://en.wikipedia.org/wiki/Clifford_algebra>`_ of the symmetric bilinear form :math:`\eta_{\mu \nu}`. Then we simply claim that the set of :math:`\Jscr_{\mu \nu}` defined by

.. math:: \Jscr_{\mu \nu} \coloneqq -\frac{\ifrak}{4} \left[ \gamma_{\mu}, \gamma_{\nu} \right]
	:label: eq_dirac_field_defn_j

solves :eq:`eq_bracket_repr_j`. To see this, let's first do a preparational calculation as follows

.. math::
	:label: eq_dirac_field_j_gamma_commutator

	\left[ \Jscr_{\mu \nu}, \gamma_{\rho} \right]
		&= -\frac{\ifrak}{4} \left[ \left[ \gamma_{\mu}, \gamma_{\nu} \right], \gamma_{\rho} \right] \\
		&= -\frac{\ifrak}{4} \left( \gamma_{\mu}\gamma_{\nu}\gamma_{\rho} - \gamma_{\nu}\gamma_{\mu}\gamma_{\rho} - \gamma_{\rho}\gamma_{\mu}\gamma_{\nu} + \gamma_{\rho}\gamma_{\nu}\gamma_{\mu} \right) \\
		&= -\frac{\ifrak}{4} \big(
			(\gamma_{\mu}\gamma_{\nu}\gamma_{\rho} + \gamma_{\mu}\gamma_{\rho}\gamma_{\nu})
			- (\gamma_{\mu}\gamma_{\rho}\gamma_{\nu} + \gamma_{\rho}\gamma_{\mu}\gamma_{\nu}) \phantom{)} \\
		&\phantom{(}\qquad\quad - (\gamma_{\nu}\gamma_{\mu}\gamma_{\rho} + \gamma_{\nu}\gamma_{\rho}\gamma_{\mu})
			+ (\gamma_{\nu}\gamma_{\rho}\gamma_{\mu} + \gamma_{\rho}\gamma_{\nu}\gamma_{\mu}) \big) \\
		&= -\frac{\ifrak}{4} \left( 2\gamma_{\mu}\eta_{\nu \rho} - 2\eta_{\mu \rho}\gamma_{\nu} - 2\gamma_{\nu}\eta_{\mu \rho} + 2\eta_{\nu \rho}\gamma_{\mu} \right) \\
		&= -\ifrak \eta_{\nu \rho}\gamma_{\mu} + \ifrak \eta_{\mu \rho}\gamma_{\nu}

Then we can verify :eq:`eq_bracket_repr_j`, starting from the left-hand-side, as follows

.. math::

	\left[ \Jscr_{\mu \nu}, \Jscr_{\rho \kappa} \right] &= -\frac{\ifrak}{4} \left[ \Jscr_{\mu \nu}, \left[ \gamma_{\rho}, \gamma_{\kappa} \right] \right] \\
		&= \frac{\ifrak}{4} \left[ \gamma_{\rho}, \left[ \gamma_{\kappa}, \Jscr_{\mu \nu} \right] \right] + \frac{\ifrak}{4} \left[ \gamma_{\kappa}, \left[ \Jscr_{\mu \nu}, \gamma_{\rho} \right] \right] \\
		&= \frac{\ifrak}{4} \left[ \gamma_{\rho}, \left( \ifrak \eta_{\nu \kappa} \gamma_{\mu} - \ifrak \eta_{\mu \kappa} \gamma_{\nu} \right) \right] + \frac{\ifrak}{4} \left[ \gamma_{\kappa}, \left( -\ifrak \eta_{\nu \rho} \gamma_{\mu} + \ifrak \eta_{\mu \rho} \gamma_{\nu} \right) \right] \\
		&= -\ifrak \eta_{\nu \kappa} \Jscr_{\rho \mu} + \ifrak \eta_{\mu \kappa} \Jscr_{\rho \nu} + \ifrak \eta_{\nu \rho} \Jscr_{\kappa \mu} - \ifrak \eta_{\mu \rho} \Jscr_{\kappa \nu}

The last expression is easily seen to be equal to the right-hand-side of :eq:`eq_bracket_repr_j` using the anti-symmetry of :math:`\Jscr_{\mu \nu}`.

In fact, the calculation :eq:`eq_dirac_field_j_gamma_commutator` may be rephrased more compactly as follows

.. math:: D(\Lambda) \gamma_{\mu} D^{-1}(\Lambda) = {\Lambda^{\nu}}_{\mu} \gamma_{\nu}
	:label: eq_dirac_field_gamma_is_vector

or in plain words, :math:`\gamma_{\mu}` is a vector.

.. dropdown:: Proof of :eq:`eq_dirac_field_gamma_is_vector`
	:animate: fade-in-slide-down
	:icon: unlock

	.. math::

		D(\Lambda) \gamma_{\mu} D^{-1}(\Lambda) &= \left( 1 + \tfrac{\ifrak}{2} \omega^{\nu \rho} \Jscr_{\nu \rho} \right) \gamma_{\mu} \left( 1 - \tfrac{\ifrak}{2} \omega^{\nu \rho} \Jscr_{\nu \rho} \right) \\
			&= \gamma_{\mu} + \tfrac{\ifrak}{2} \omega^{\nu \rho} \left[ \Jscr_{\nu \rho}, \gamma_{\mu} \right] \\
			&= \gamma_{\mu} + \tfrac{\ifrak}{2} \omega^{\nu \rho} \left( -\ifrak \eta_{\rho \mu} \gamma_{\nu} + \ifrak \eta_{\nu \mu} \gamma_{\rho} \right) \\
		&= \left( {\delta^{\nu}}_{\mu} + {\omega^{\nu}}_{\mu} \right) \gamma_{\nu} = {\Lambda^{\nu}}_{\mu} \gamma_{\nu}

Using the very definition :eq:`eq_dirac_field_defn_j`, one then sees that :math:`\Jscr_{\mu \nu}` is an anti-symmetric tensor in the sense that

.. math:: D(\Lambda) \Jscr_{\mu \nu} D^{-1}(\Lambda) = \Lambda^{\rho}_{\mu} \Lambda^{\kappa}_{\nu} \Jscr_{\rho \kappa}

Indeed, one can construct *all* anti-symmetric tensors as follows

.. math::

	\Ascr_{\mu \nu \rho} &\coloneqq \gamma_{\mu}\gamma_{\nu}\gamma_{\rho} \pm \text{ signed permutations} \\
	\Pscr_{\mu \nu \rho \kappa} &\coloneqq \gamma_{\mu}\gamma_{\nu}\gamma_{\rho}\gamma_{\kappa} \pm \text{ signed permutations}

There can be no more because we're constrained by the :math:`4`-dimensional spacetime. We note that these anti-symmetric tensors form a complete basis of all matrices that can be constructed out of the :math:`\gamma`-matrices. This is because using :eq:`eq_dirac_field_clifford_algebra`, any product of :math:`\gamma`-matrices can be written as a linear combination of the anti-symmetric tensors with coefficients the metric tensors.

Now we claim that the matrices :math:`1, \gamma_{\mu}, \Jscr_{\mu \nu}, \Ascr_{\mu \nu \rho}` and :math:`\Pscr_{\mu \nu \rho \kappa}` are all linearly independent.

.. dropdown:: Proof of the anti-symmetric tensors being linearly independent
	:animate: fade-in-slide-down
	:icon: unlock

	One way to see that :math:`1, \gamma_{\mu}, \Jscr_{\mu \nu}, \Ascr_{\mu \nu \rho}` and :math:`\Pscr_{\mu \nu \rho \kappa}` are linearly independent is to observe that they transform differently under conjugation by :math:`D(\Lambda)`. But more directly, an inner product on the matrices can be defined by taking the trace of the product matrix. We claim that the anti-symmetric tensors are orthogonal to each other with respect to this inner product, and hence linearly independent.

	Instead of working out all the details, let's take a look at a few prototypical cases.

	1. :math:`\op{tr}(\Jscr_{\mu \nu}) = 0` because the trace of a commutator vanishes.
	2. :math:`\op{tr}(\gamma_{\mu} \gamma_{\nu}) = 0` for :math:`\mu \neq \nu` by :eq:`eq_dirac_field_clifford_algebra`.
	3. It's slightly tricker to see that :math:`\gamma_{\mu}` itself is also traceless, but this is again a consequence of the Clifford algebra relations :eq:`eq_dirac_field_clifford_algebra`, which we demonstrate as follows

	   .. math::

			\op{tr}(\gamma_{\mu})
				= \op{tr}(\gamma_{\nu} \gamma_{\mu} \gamma^{-1}_{\nu})
				= -\op{tr}(\gamma_{\mu} \gamma_{\nu} \gamma^{-1}_{\nu})
				= -\op{tr}(\gamma_{\mu})
				\implies \op{tr}(\gamma_{\mu}) = 0

	   where :math:`\nu` is any index different from :math:`\mu` so that :math:`\{ \gamma_{\mu}, \gamma_{\nu} \} = 0`.

Counting these linearly independent matrices, we see that there are :math:`1 + \binom{4}{1} + \binom{4}{2} + \binom{4}{3} + \binom{4}{4} = 16` of them. It means that the size of the :math:`\gamma_{\mu}` matrices is at least :math:`4 \times 4`.

It turns out that there exists indeed a solution of :eq:`eq_dirac_field_clifford_algebra` in terms of :math:`4 \times 4` matrices, conveniently known as the `gamma matrices <https://en.wikipedia.org/wiki/Gamma_matrices>`_, which we define as follows

.. math::
	:label: eq_dirac_field_defn_gamma_matrices

	\gamma_0 \coloneqq -\ifrak \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad \
		\bm{\gamma} \coloneqq -\ifrak \begin{bmatrix*}[r] 0 & \bm{\sigma} \\ -\bm{\sigma} & 0 \end{bmatrix*}

where :math:`\bm{\sigma} = (\sigma_1, \sigma_2, \sigma_3)` is made up of the so-called `Pauli matrices <https://en.wikipedia.org/wiki/Pauli_matrices>`_ defined as follows

.. math::
	:label: eq_pauli_matrices

	\sigma_1 \coloneqq \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad \
		\sigma_2 \coloneqq \begin{bmatrix*}[r] 0 & -\ifrak \\ \ifrak & 0 \end{bmatrix*}, \quad \
		\sigma_3 \coloneqq \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}

Indeed the Pauli matrices make up a solution to not only a :math:`3`-dimensional Clifford algebra with respect to the Euclidean inner product, but also an angular momentum representation if multiplied by :math:`1/2`, or more precisely, a spin-:math:`1/2` representation. Note that this representation, in terms of Hermitian matrices, is different from the one given by :eq:`eq_rotation_j_matrix` with :math:`\jfrak = 1/2`, in terms of real matrices. One can verify that they differ by a change of basis. Note, however, that as far as the terminology is concerned, one often uses Hermitian and real interchangeably.

Now using the Clifford relations for both gamma and Pauli matrices, we can evaluate :eq:`eq_dirac_field_defn_j` as follows

.. math::
	:label: eq_dirac_field_jscr_matrix

	\Jscr_{ij}
		&= -\frac{\ifrak}{4} \left[ \gamma_i, \gamma_j \right]
		= -\frac{\ifrak}{2} \epsilon_{ij} \gamma_i \gamma_j
		= -\frac{\ifrak}{2} \epsilon_{ij} \begin{bmatrix} \sigma_i \sigma_j & 0 \\ 0 & \sigma_i \sigma_j \end{bmatrix}
		= \frac{1}{2} \epsilon_{ijk} \begin{bmatrix} \sigma_k & 0 \\ 0 & \sigma_k \end{bmatrix} \\
	\Jscr_{i0}
		&= -\frac{\ifrak}{4} \left[ \gamma_i, \gamma_0 \right]
		= -\frac{\ifrak}{2} \gamma_i \gamma_0
		= \frac{\ifrak}{2} \begin{bmatrix} \sigma_i & 0 \\ 0 & \sigma_i \end{bmatrix}

.. _paragraph_dirac_field_representation_not_unitary:

where :math:`i, j \in \{1,2,3\}` and :math:`\epsilon` is the totally anti-symmetric sign. We see that the representation :math:`\Jscr_{\mu \nu}` is in fact reducible. Moreover, we see that that the corresponding representation :math:`D` of the Lorentz group given by :eq:`eq_dirac_field_linearize_representation` is not unitary, since while :math:`\Jscr_{ij}` are Hermitian, :math:`\Jscr_{i0}` are anti-Hermitian. The fact that :math:`D` is not unitary will have consequences when we try to construct the interaction density as in :eq:`eq_construct_interaction_density_by_fields`, because products like :math:`\psi^{\dagger} \psi` will not be a scalar (see :ref:`Construction of the Interaction Density for Dirac Fields <sec_construction_of_the_interaction_density>`).

Next let's consider the parity transformation, i.e., the transformation under spatial inversion, in the context of gamma matrices. In comparison with the transformation laws :eq:`eq_hpjk_cojugated_by_space_and_time_inversions`, we can define

.. math:: \beta \coloneqq \ifrak \gamma_0 = \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}
	:label: eq_dirac_field_beta_matrix

as the parity transformation. Indeed, we clearly have :math:`\beta^2 = 1`. Moreover, it follows from the Clifford relations :eq:`eq_dirac_field_clifford_algebra` that

.. math::
	:label: eq_dirac_field_beta_conjugate_gamma

	\beta \gamma_i \beta^{-1} = -\gamma_i, \quad \beta \gamma_0 \beta^{-1} = \gamma_0

which, in turn, implies that

.. math::

	\beta \Jscr_{ij} \beta^{-1} = \Jscr_{ij}, \quad \beta \Jscr_{0i} \beta^{-1} = - \Jscr_{0i}

which is consistent with :eq:`eq_hpjk_cojugated_by_space_and_time_inversions` if we think of :math:`\Jscr_{ij}` as the angular momenta and :math:`\Jscr_{0i}` as the boosts.

In connection to the non-unitarity of :math:`D(\Lambda)`, let's note that since :math:`\beta \gamma_{\mu}^{\dagger} \beta^{-1} = -\gamma_{\mu}` (which can be verified by :eq:`eq_dirac_field_beta_conjugate_gamma` and :eq:`eq_dirac_field_defn_gamma_matrices`), we have :math:`\beta \Jscr_{\mu \nu}^{\dagger} \beta^{-1} = \Jscr_{\mu \nu}`, and therefore

.. math:: \beta D^{\dagger}(\Lambda) \beta^{-1} = D^{-1}(\Lambda)
	:label: eq_dirac_field_pseudo_unitarity_of_d_matrix

in light of :eq:`eq_dirac_field_linearize_representation`. This identity will be useful when we later construct the interaction density.

At last we'll introduce yet another special element to the family of gamma matrices, namely :math:`\gamma_5`, defined as follows

.. math::
	:label: eq_dirac_field_defn_gamma_5

	\gamma_5 \coloneqq -\ifrak \gamma_0 \gamma_1 \gamma_2 \gamma_3 = \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}

One nice thing about :math:`\gamma_5` is that it anti-commutes with all other :math:`\gamma` matrices, and in particular

.. math:: \beta \gamma_5 = -\gamma_5 \beta
	:label: eq_dirac_field_gamma_5_anti_commutes_beta

In fact, the collection :math:`\gamma_0, \gamma_1, \gamma_2, \gamma_3, \gamma_5` makes exactly a :math:`5`-dimensional spacetime Clifford algebra.


Construction of Dirac fields
++++++++++++++++++++++++++++

As in the case of scalar and vector fields, let's write the Dirac fields as follows

.. math::
	:label: eq_dirac_field_psi

	\psi^+_{\ell}(x) &= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma) a(\pbf, \sigma) \\
	\psi^{-c}_{\ell}(x) &= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma)

Moreover, using :math:`\Jscr_{ij}` as given by :eq:`eq_dirac_field_jscr_matrix`, we can write the :math:`\pbf = 0` conditions :eq:`eq_j_intertwines_u_and_v` as follows

.. math::
	:label: eq_dirac_field_sigma_intertwines_j_by_u_and_v

	\frac{1}{2} \sum_m \bm{\sigma}_{m' m} u_{m \pm}(0, \sigma) &= \sum_{\sigma'} u_{m' \pm}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} \\
	-\frac{1}{2} \sum_m \bm{\sigma}_{m' m} v_{m \pm}(0, \sigma) &= \sum_{\sigma'} v_{m' \pm}(0, \sigma') \Jbf^{(\jfrak) \ast}_{\sigma' \sigma}

Here the signs :math:`\pm` correspond to the two identical irreducible representations of :math:`\Jscr_{ij}`, which is obvious from :eq:`eq_dirac_field_jscr_matrix`, while :math:`m` and :math:`m'` index the Pauli matrices :eq:`eq_pauli_matrices`.

Now if we think of :math:`u_{m \pm}(0, \sigma)` as matrix elements of a matrix :math:`U_{\pm}` and similarly for :math:`v`, then :eq:`eq_dirac_field_sigma_intertwines_j_by_u_and_v` can be rewritten compactly in matrix notation as follows

.. math::
	:label: eq_dirac_field_sigma_intertwines_j_by_u_and_v_matrix_form

	\tfrac{1}{2} \bm{\sigma} U_{\pm} &= U_{\pm} \Jbf^{(\jfrak)} \\
	-\tfrac{1}{2} \bm{\sigma} V_{\pm} &= V_{\pm} \Jbf^{(\jfrak) \ast}

We recall that both :math:`\tfrac{1}{2} \bm{\sigma}` and :math:`\Jbf^{(\jfrak)}` (as well as :math:`-\Jbf^{(\jfrak) \ast}`) are irreducible representations of the rotation group (or rather, its Lie algebra). We first claim that :math:`U_{\pm}` must be isomorphism. Indeed, the kernel of :math:`U_{\pm}` is easily seen to be an invariant subspace under the action of :math:`\Jbf^{(\jfrak)}`, and hence must be null if :math:`U_{\pm} \neq 0`. On the other hand, the image of :math:`U_{\pm}` is an invariant subspace under the action of :math:`\bm{\sigma}`, and hence must be the whole space if :math:`U_{\pm} \neq 0`. It follows then that the rank of :math:`\Jbf^{(\jfrak)}` and :math:`U_{\pm}` must be the same as :math:`\bm{\sigma}`, which is :math:`2`. The same argument applies also to :math:`V_{\pm}`. In particular we must have :math:`\jfrak = \tfrac{1}{2}`, or in other words, the Dirac fields describe :math:`\tfrac{1}{2}`-spin particles.

.. note::
	The mathematical argument above is commonly known as `Schur's lemma <https://en.wikipedia.org/wiki/Schur%27s_lemma>`_.

The matrix form of :math:`\Jbf^{(1/2)}` according to :eq:`eq_rotation_j_matrix` is given By

.. math::

	J_1^{(1/2)} = \frac{1}{2} \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad
		J_2^{(1/2)} = \frac{\ifrak}{2} \begin{bmatrix*}[r] 0 & -1 \\ 1 & 0 \end{bmatrix*}, \quad
		J_3^{(1/2)} = \frac{1}{2} \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}

Comparing with the Pauli matrices :eq:`eq_pauli_matrices`, we see that

.. math::

	\Jbf^{(1/2)} = \tfrac{1}{2} \bm{\sigma}, \quad -\Jbf^{(1/2) \ast} = \tfrac{1}{2} \sigma_2 \bm{\sigma} \sigma_2

Hence the first equation in :eq:`eq_dirac_field_sigma_intertwines_j_by_u_and_v_matrix_form` may be rewritten as :math:`\bm{\sigma} U_{\pm} = U_{\pm} \bm{\sigma}`. One can apply Schur's lemma once again to conclude that :math:`U_{\pm}` must be a scalar.

.. dropdown:: Proof of :math:`U_{\pm}` being scalar
	:animate: fade-in-slide-down
	:icon: unlock

	Since :math:`\bm{\sigma}` commutes with any scalar it commutes in particular with :math:`U_{\pm} - \lambda` for any :math:`\lambda \in \Cbb`. It follows that :math:`U_{\pm} - \lambda` is either an isomorphism or zero. The later must be the case if :math:`\lambda` is an eigenvalue of :math:`U_{\pm}`. Hence :math:`U_{\pm}` must be a scalar.

A similar argument can be applied to :math:`V_{\pm}` by rewriting the second equation in :eq:`eq_dirac_field_sigma_intertwines_j_by_u_and_v_matrix_form` as :math:`\bm{\sigma} (V_{\pm} \sigma_2) = (V_{\pm} \sigma_2) \bm{\sigma}`. Hence :math:`V_{\pm}` must be proportional to :math:`\sigma_2`.

Going back to :math:`u_{m \pm}(0, \sigma)` and :math:`v_{m \pm}(0, \sigma)` from :math:`U_{\pm}` and :math:`V_{\pm}`, respectively, we have concluded that

.. math::

	u_{m \pm}(0, \sigma) &= c_{\pm} \delta_{m \sigma} \\
	v_{m \pm}(0, \sigma) &= -\ifrak d_{\pm} (\sigma_2)_{m \sigma}

where we've inserted the extra factor :math:`-\ifrak` in front of :math:`d_{\pm}` to make the final results look more uniform. We can further unwrap this result in matrix notations as follows

.. math::
	:label: eq_dirac_field_u_and_v_matrix_prelim

	\begin{alignat*}{2}
		u(0, 1/2) &= \begin{bmatrix*}[l] c_+ \\ 0 \\ c_- \\ 0 \end{bmatrix*}, \quad u(0, -1/2) &&= \begin{bmatrix*}[l] 0 \\ c_+ \\ 0 \\ c_- \end{bmatrix*} \\
		v(0, 1/2) &= \begin{bmatrix*}[l] 0 \\ d_+ \\ 0 \\ d_- \end{bmatrix*}, \quad v(0, -1/2) &&= -\begin{bmatrix*}[l] d_+ \\ 0 \\ d_- \\ 0 \end{bmatrix*}
	\end{alignat*}

Now spinors at finite momentum can be determined as usual by :eq:`eq_annihilation_and_creation_transformation_simplified_by_boost` as follows

.. math::
	:label: eq_dirac_field_spinor_u_and_v_at_finite_momentum

	u(\pbf, \sigma) &= \sqrt{m/p_0} D(L(p)) u(0, \sigma) \\
	v(\pbf, \sigma) &= \sqrt{m/p_0} D(L(p)) v(0, \sigma)

In general the constants :math:`c_{\pm}` and :math:`d_{\pm}` may be arbitrary. However, if we further assume the conservation of parity, or in other words, that the fields transform covariantly under the spatial inversion, then the constants, and henceforth the spinors, can be determined uniquely.

To spell out the details, let's apply the spatial inversion transformation laws from :eq:`eq_creation_operator_cpt_conjugation` to :eq:`eq_dirac_field_psi` as follows

.. math::
	:label: eq_dirac_field_spatial_inversion_acts_on_psi

	& U(\Pcal) \psi^+(x) U^{-1}(\Pcal) \\
		&\quad = (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot x} u(\pbf, \sigma) a(-\pbf, \sigma) \\
		&\quad = (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot \Pcal x} u(-\pbf, \sigma) a(\pbf, \sigma) \\
		&\quad = (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot \Pcal x} \sqrt{m/p_0}~\beta D(L(p)) \beta u(0, \sigma) a(\pbf, \sigma) \\
	& U(\Pcal) \psi^{- c}(x) U^{-1}(\Pcal) \\
		&\quad = (2\pi)^{-3/2} \eta^c \sum_{\sigma} \int d^3 p~e^{-\ifrak p \cdot \Pcal x} \sqrt{m/p_0}~\beta D(L(p))\beta v(0, \sigma) a^{c \dagger}(\pbf, \sigma)

Here we've evaluated :math:`u(-\pbf, \sigma)` and :math:`v(-\pbf, \sigma)` in a way similar to :eq:`eq_vector_field_revert_momentum_transformation`. Namely, we've used :eq:`eq_dirac_field_spinor_u_and_v_at_finite_momentum` together with the fact that :math:`D(\Pcal) = \beta` defined by :eq:`eq_dirac_field_beta_matrix`.

Now in order for :math:`U(\Pcal) \psi^+(x) U^{-1}(\Pcal)` to be proportional to :math:`\psi^+(x)`, we observe that it would be sufficient (and most likely also necessary) if :math:`u(0, \sigma)` is an eigenvector of :math:`\beta`. Similar argument applies to :math:`\psi^-(x)` as well. So let's suppose

.. math::
	:label: eq_dirac_field_beta_eigenvalue_of_u_and_v

	\beta u(0, \sigma) &= b_+ u(0, \sigma) \\
	\beta v(0, \sigma) &= b_- v(0, \sigma)

with :math:`b^2_{\pm} = 1` since :math:`\beta^2 = 1`. Given this, we can rewrite :eq:`eq_dirac_field_spatial_inversion_acts_on_psi` as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Pcal) \psi^+(x) U^{-1}(\Pcal) &= \eta^{\ast} b_+ \beta \psi^+(\Pcal x)
		\label{eq_dirac_field_psi_plus_conjugated_by_u} \\
		U(\Pcal) \psi^{- c}(x) U^{-1}(\Pcal) &= \eta^c b_- \beta \psi^{- c}(\Pcal x)
		\label{eq_dirac_field_psi_minus_conjugated_by_u}
	\end{align}

so that the parity is indeed conserved, or rather, transformed covariantly.

Now using :eq:`eq_dirac_field_beta_eigenvalue_of_u_and_v` (and appropriately rescaling), we can rewrite :eq:`eq_dirac_field_u_and_v_matrix_prelim` as follows

.. math::

	\begin{alignat*}{2}
		u(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 1 \\ 0 \\ b_+ \\ 0 \end{bmatrix*}, \quad u(0, -1/2) &&= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 0 \\ 1 \\ 0 \\ b_+ \end{bmatrix*} \\
		v(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 0 \\ 1 \\ 0 \\ b_- \end{bmatrix*}, \quad v(0, -1/2) &&= -\frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 1 \\ 0 \\ b_- \\ 0 \end{bmatrix*}
	\end{alignat*}

So far by assuming the parity conservation, we've managed to reduced the free parameters from :math:`c_{\pm}, d_{\pm}` to :math:`b_{\pm}`. What eventually pins down the spinors is, once again, the causality condition. To see this, let's try to construct the Dirac field following :eq:`eq_defn_psi_field` as follows

.. math::
	:label: eq_dirac_field_psi_field_raw

	\psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^{- c}(x)

As usual, the (anti-)commutator can be calculated using :eq:`eq_dirac_field_psi` as follows

.. math::
	:label: eq_dirac_field_commutator_raw

	& \left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_{\pm} \\
	&\quad = (2\pi)^{-3} \int d^3 p~\left( |\kappa|^2 N_{\ell \ell'}(\pbf) \exp(\ifrak p \cdot (x-y)) \pm |\lambda|^2 M_{\ell \ell'}(\pbf) \exp(-\ifrak p \cdot (x-y)) \right)

where :math:`N_{\ell \ell'}` and :math:`M_{\ell \ell'}` are the spin sums defined as follows

.. math::
	:label: eq_dirac_field_n_and_m_matrix_as_spinor_sum

	N_{\ell \ell'}(\pbf) &= \sum_{\sigma} u_{\ell}(\pbf, \sigma) u^{\ast}_{\ell'}(\pbf, \sigma) \\
	M_{\ell \ell'}(\pbf) &= \sum_{\sigma} v_{\ell}(\pbf, \sigma) v^{\ast}_{\ell'}(\pbf, \sigma)

To evaluate the spin sums, we first turn back to the zero-momentum case and use :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_u}` and :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_v}` to express the values in terms of :math:`\beta` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		N(0) &= \sum_{\sigma} u(0, \sigma) u^{\dagger}(0, \sigma) &&= \frac{1 + b_+ \beta}{2}
		\label{eq_dirac_field_n_matrix_at_zero} \\
		M(0) &= \sum_{\sigma} v(0, \sigma) v^{\dagger}(0, \sigma) &&= \frac{1 + b_- \beta}{2}
		\label{eq_dirac_field_m_matrix_at_zero}
	\end{alignat}

where :math:`\dagger` here means transpose conjugation. The easiest way to see this is probably to momentarily forget about the spin :math:`z`-component :math:`\sigma` so that :math:`\beta` as in :eq:`eq_dirac_field_beta_matrix` behaves like a :math:`2 \times 2` matrix with obvious eigenvectors :math:`[1, 1]^T` for :math:`b_{\pm} = 1` and :math:`[1, -1]^T` for :math:`b_{\pm} = -1`. Then :math:`\eqref{eq_dirac_field_n_matrix_at_zero}` and :math:`\eqref{eq_dirac_field_m_matrix_at_zero}` can be verified by a direct calculation. Here the superscript :math:`T` means taking transpose so we have column vectors.

Now we can evaluate the :math:`N` and :math:`M` matrices in terms of :math:`\beta` as follows

.. math::
	:nowrap:

	\begin{align}
		N(\pbf) &\xlongequal{\eqref{eq_dirac_field_n_matrix_as_spinor_sum}} \sum_{\sigma} u(\pbf, \sigma) u^{\dagger}(\pbf, \sigma)
			\label{eq_dirac_field_n_matrix_first_evaluation} \\
			&\xlongequal{\eqref{eq_dirac_field_spinor_u_at_finite_momentum}} \frac{m}{p_0} D(L(p)) \left( \sum_{\sigma} u(0, \sigma) u^{\dagger}(0, \sigma) \right) D^{\dagger}(L(p)) \nonumber \\
			&\xlongequal{\eqref{eq_dirac_field_n_matrix_at_zero}} \frac{m}{2p_0} D(L(p)) (1 + b_+ \beta) D^{\dagger}(L(p)) \nonumber \\
		M(\pbf) &\xlongequal{\phantom{\eqref{eq_dirac_field_n_matrix_at_zero}}} \frac{m}{2p_0} D(L(p)) (1 + b_- \beta) D^{\dagger}(L(p))
		\label{eq_dirac_field_m_matrix_first_evaluation}
	\end{align}

To go further, we need to invoke the the transformation laws of the gamma matrices and their relatives, e.g., :math:`\beta`, under :math:`D(\Lambda)` from :ref:`sec_dirac_representation_and_gamma_matrices`. More precisely, using the pseudo-unitarity :eq:`eq_dirac_field_pseudo_unitarity_of_d_matrix` we have the following

.. math::
	:nowrap:

	\begin{align}
		D(L(p)) \beta D^{\dagger}(L(p)) &\xlongequal{\phantom{\eqref{eq_dirac_field_beta_matrix}}} D(L(p)) D^{-1}(L(p)) \beta = \beta \nonumber \\
		D(L(p)) D^{\dagger}(L(p)) &\xlongequal{\phantom{\eqref{eq_dirac_field_beta_matrix}}} D(L(p)) \beta D^{-1}(L(p)) \beta
			\label{eq_dirac_field_d_d_dagger_product} \\
			&\xlongequal{\eqref{eq_dirac_field_beta_matrix}} \ifrak D(L(p)) \gamma_0 D(L(p)) \beta \nonumber \\
			&\xlongequal{\eqref{eq_dirac_field_gamma_is_vector}} \ifrak \left( L^{-1}(p) \right)_0^{\mu} \gamma_{\mu} \beta \nonumber \\
			&\xlongequal{\eqref{eq_L_transformation_for_massive_2}} -\ifrak p^{\mu} \gamma_{\mu} \beta / m \
			= -\ifrak p_{\mu} \gamma^{\mu} \beta / m \nonumber
	\end{align}

Plugging them into :math:`\eqref{eq_dirac_field_n_matrix_first_evaluation}` and :math:`\eqref{eq_dirac_field_m_matrix_first_evaluation}`, respectively, we can continue our evaluation as follows

.. math::
	:nowrap:

	\begin{align*}
		N(\pbf) &= \frac{1}{2p_0} \left( -\ifrak p_{\mu} \gamma^{\mu} + b_+ m \right) \beta \\
		M(\pbf) &= \frac{1}{2p_0} \left( -\ifrak p_{\mu} \gamma^{\mu} + b_- m \right) \beta
	\end{align*}

Now that we've finished evaluating the spin sums, we can plug them into :eq:`eq_dirac_field_commutator_raw` to get the following evaluation of the (anti-)commutator

.. math::
	:nowrap:

	\begin{align}
		\left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_{\pm} &= (2\pi)^{-3} \int \frac{d^3 p}{2p_0}~\big( |\kappa|^2 (-\ifrak p_{\mu} \gamma^{\mu} + b_+ m) e^{\ifrak p \cdot (x-y)} \beta \phantom{\big)}
		\label{eq_dirac_field_commutator_first_evaluation} \\
			&\phantom{= \phantom{\big(}} \pm |\lambda|^2 (-\ifrak p_{\mu} \gamma^{\mu} + b_- m) e^{-\ifrak p \cdot (x-y)} \beta \big)_{\ell \ell'} \nonumber \\
			&= \left( |\kappa|^2 \left( -\ifrak \gamma^{\mu} \p_{x_{\mu}} + b_+ m \right) \Delta_+(x-y) \beta \pm |\lambda|^2 \left( -\ifrak \gamma^{\mu} \p_{y_{\mu}} + b_- m \right) \Delta_+(y-x) \beta \right)_{\ell \ell'} \nonumber
	\end{align}

where :math:`\Delta_+` is defined by :eq:`eq_defn_Delta_plus`. Recall that :math:`\Delta_+(x) = \Delta_+(-x)` for space-like :math:`x`. Hence for space-like :math:`x-y`, the following holds

.. math::
	:nowrap:

	\begin{equation}
		\p_{x_{\mu}} \Delta_+(x-y) = \p_{x_{\mu}} \Delta_+(y-x) = -\p_{y_{\mu}} \Delta_+(y-x)
		\label{eq_delta_plus_derivative_is_odd}
	\end{equation}

.. error::
	:math:`\eqref{eq_delta_plus_derivative_is_odd}` is what is claimed in the bottom of page 223 in [Wei95]_, but it's not true! Indeed, the fact that :math:`\Delta_+(x)` is even for space-like :math:`x` only implies that :math:`\p_i \Delta_+(x)` is odd for spatial indexes :math:`i=1,2,3`. However :math:`\dot{\Delta}_+(x)` is not odd even for space-like :math:`x`. It has a serious consequence that the anti-commutator doesn't vanish, which in turn means that the Dirac fields, as constructed here, cannot be arbitrarily assembled into a causal interaction density. This error is obviously not because of Weinberg's ignorance, since he later also pointed out the non-vanishing of the anti-commutator on page 295.

It follows that in order for :math:`\eqref{eq_dirac_field_commutator_first_evaluation}` to vanish for space-separated :math:`x` and :math:`y`, the following must hold

.. math::
	:nowrap:

	\begin{align*}
		|\kappa|^2 \mp |\lambda|^2 &= 0 \\
		|\kappa|^2 b_+ \pm |\lambda|^2 b_- &= 0
	\end{align*}

We see that first of all, the top sign applies, which means in particular that we must be considering the anti-commutator in :math:`\eqref{eq_dirac_field_commutator_first_evaluation}`, or in other words, the Dirac fields must be fermionic. In addition, we must also have :math:`|\kappa| = |\lambda|` and :math:`b_+ + b_- = 0`.

By the usual phase adjustments on the creation and annihilation operators and rescaling, we can arrange so that :math:`\kappa = \lambda = 1`. Recalling :eq:`eq_dirac_field_gamma_5_anti_commutes_beta` and replacing :math:`\psi` with :math:`\gamma_5 \psi` if necessary, we can arrange so that :math:`b_{\pm} = \pm 1`. Putting these all together, we have evaluated the Dirac field :eq:`eq_dirac_field_psi_field_raw` as follows

.. math::
	:nowrap:

	\begin{align}
		\psi(x) &= \psi^+(x) + \psi^{- c}(x)
			\label{eq_dirac_field_psi_field} \\
			&= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~\left( e^{\ifrak p \cdot x} u(\pbf, \sigma) a(\pbf, \sigma) + e^{-\ifrak p \cdot x} v(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right) \nonumber
	\end{align}

where the zero-momentum spinors are

.. math::
	:nowrap:

	\begin{alignat}{2}
		u(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 1 \\ 0 \\ 1 \\ 0 \end{bmatrix*}, \quad &&u(0, -1/2) = \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ 0 \\ 1 \end{bmatrix*}
		\label{eq_dirac_field_u_spinor_zero_momentum} \\
		v(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ 0 \\ -1 \end{bmatrix*}, \quad &&v(0, -1/2) = -\frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 1 \\ 0 \\ -1 \\ 0 \end{bmatrix*}
		\label{eq_dirac_field_v_spinor_zero_momentum}
	\end{alignat}

and the spin sums are

.. math::
	:nowrap:

	\begin{align}
		N(\pbf) &= \frac{1}{2p_0} \left( -\ifrak p^{\mu} \gamma_{\mu} + m \right) \beta
		\label{eq_dirac_field_spin_sum_u} \\
		M(\pbf) &= \frac{1}{2p_0} \left( -\ifrak p^{\mu} \gamma_{\mu} - m\right) \beta
		\label{eq_dirac_field_spin_sum_v}
	\end{align}

and the anti-commutator, calculated by plugging the spin sums into :eq:`eq_dirac_field_commutator_raw`, is

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_+ = \left( \left( -\gamma^{\mu}\p_{\mu} + m \right) \beta \right)_{\ell \ell'} \Delta(x-y)
		\label{eq_dirac_field_commutator}
	\end{equation}

where :math:`\p_{\mu} = \p_{x_{\mu}}` and :math:`\Delta(x-y)` is defined by :eq:`eq_defn_Delta`.

For :math:`\psi(x)` defined by :math:`\eqref{eq_dirac_field_psi_field}` to transform covariantly under spatial inversion, we recall :math:`\eqref{eq_dirac_field_psi_plus_conjugated_by_u}` and :math:`\eqref{eq_dirac_field_psi_minus_conjugated_by_u}` to conclude that

.. math::
	:nowrap:

	\begin{equation}
		\eta^{\ast} + \eta^c = 0
		\label{eq_dirac_field_spatial_parity_relation}
	\end{equation}

It follows that :math:`\eta \eta^c = -1`, or in other words, the intrinsic parity of the state consisting of a spin :math:`1/2` particle and its antiparticle is odd in the sense of :ref:`sec_parity_symmetry`. The parity transformation law for Dirac fields is as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \psi(x) U^{-1}(\Pcal) = \eta^{\ast} \beta \psi(\Pcal x)
		\label{eq_dirac_field_spatial_inversion_transformation_law}
	\end{equation}

.. dropdown:: The Dirac equation
	:animate: fade-in-slide-down

	Recall from :math:`\eqref{eq_klein_gordon}` that a general quantum field satisfies Klein-Gordon equation, which is a second-order differential equation. It turns out that the Dirac fields satisfy a first-order equation, known as the `Dirac equation <https://en.wikipedia.org/wiki/Dirac_equation>`_.

	To see this, we recall that in the derivation of :math:`\eqref{eq_dirac_field_d_d_dagger_product}`, we've essentially proved the following identity

	.. math::
		:nowrap:

		\begin{equation*}
			D(L(p)) \beta D^{-1}(L(p)) = -\ifrak p_{\mu} \gamma^{\mu} / m
		\end{equation*}

	Applying the left-hand-side to :math:`u(\pbf, \sigma)`, we see that

	.. math::
		:nowrap:

		\begin{align*}
			D(L(p)) \beta D^{-1}(L(p)) u(\pbf, \sigma) &\xlongequal{\eqref{eq_dirac_field_spinor_u_at_finite_momentum}} \sqrt{m/p_0} D^{-1}(L(p)) \beta u(0, \sigma) \\
				&\xlongequal{\eqref{eq_dirac_field_beta_eigenvalue_of_u}} \sqrt{m/p_0} D(L(p)) u(0, \sigma) \\
				&\xlongequal{\eqref{eq_dirac_field_spinor_u_at_finite_momentum}} u(\pbf, \sigma)
		\end{align*}

	It follows that

	.. math::
		:nowrap:

		\begin{equation*}
			\left( \ifrak p_{\mu} \gamma^{\mu} + m \right) u(\pbf, \sigma) = 0
		\end{equation*}

	Similarly one can show that :math:`v(\pbf, \sigma)` satisfies the following

	.. math::
		:nowrap:

		\begin{equation*}
			(-\ifrak p_{\mu} \gamma^{\mu} + m) v(\pbf, \sigma) = 0
		\end{equation*}

	Combining these identities with :math:`\eqref{eq_dirac_field_psi_field}`, we conclude that the Dirac fields satisfy the following Dirac equation

	.. math::
		:nowrap:

		\begin{equation*}
			\left( \gamma^{\mu} \p_{\mu} + m \right) \psi(x) = 0
		\end{equation*}

	However, unlike the original derivation of Dirac, we've derived it here as a consequence of parity conservation. In fact, we note that the Dirac equation is not something completely different from the Klein-Gordon equation, because using the Clifford algebra relations :eq:`eq_dirac_field_clifford_algebra`, we see that :math:`\gamma^{\mu} \p_{\mu}` is actually a square root of the D'Alembert operator :math:`\square`. This is allegedly one of the motivations of Dirac to find these gamma matrices in the first place.

The CPT symmetries
++++++++++++++++++

The transformation law of Dirac fields under spatial inversion has already been worked out in :math:`\eqref{eq_dirac_field_spatial_inversion_transformation_law}`, so we're left to work out the transformation laws under time inversion and charge conjugation.

Recall from :ref:`sec_space_and_time_inversions` that the time-inversion operator :math:`U(\Tcal)` is complex anti-linear. Hence to work out the transformation law under time inversion, we'll need to work out the complex-conjugated spinors :math:`u^{\ast}(\pbf, \sigma)` and :math:`v^{\ast}(\pbf, \sigma)`. Now in light of :math:`\eqref{eq_dirac_field_spinor_u_at_finite_momentum}` and :math:`\eqref{eq_dirac_field_spinor_v_at_finite_momentum}`, and the fact that the spinors are real at zero-momentum, we just need to work out the complex conjugate :math:`D^{\ast}(L(p))` in terms of :math:`D(L(p))` and the gamma matrices. Now according to :eq:`eq_dirac_field_linearize_representation`, it suffices to work out :math:`\Jscr^{\ast}_{\mu \nu}`. Finally according to :eq:`eq_dirac_field_defn_j`, it suffices to work out :math:`\gamma^{\ast}_{\mu}`.

Inspecting the explicit forms of the gamma matrices given by :eq:`eq_dirac_field_defn_gamma_matrices` and :eq:`eq_pauli_matrices` we see that :math:`\gamma_0, \gamma_1, \gamma_3` are anti-Hermitian while :math:`\gamma_2` is Hermitian, or more explicitly

.. math::
	:nowrap:

	\begin{equation*}
		\gamma^{\ast}_0 = -\gamma_0, \quad \gamma^{\ast}_1 = -\gamma_1, \quad \gamma^{\ast}_2 = \gamma_2, \quad \gamma^{\ast}_3 = -\gamma_3
	\end{equation*}

Using the Clifford algebra relations, this can be written more concisely as follows

.. math::
	:nowrap:

	\begin{equation}
		\gamma^{\ast}_{\mu} = \gamma_2 \gamma_{\mu} \gamma_2
		\label{eq_dirac_field_gamma_conjugation_by_gamma_2}
	\end{equation}

While this result could've be satisfactory in its own right, as we'll see, it'll be more convenient to factor out a :math:`\beta` matrix. Hence we're motivated to introduce yet another special matrix

.. math::
	:nowrap:

	\begin{equation}
		\Cscr \coloneqq \gamma_2 \beta = -\ifrak \begin{bmatrix*}[r] \sigma_2 & 0 \\ 0 & -\sigma_2 \end{bmatrix*}
		\label{eq_dirac_field_defn_c_matrix}
	\end{equation}

and rewrite :math:`\eqref{eq_dirac_field_gamma_conjugation_by_gamma_2}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\gamma^{\ast}_{\mu} = \beta \Cscr \gamma_{\mu} \Cscr^{-1} \beta
	\end{equation*}

where we also note that :math:`(\Cscr^{-1} \beta)^{-1} = \beta \Cscr`. It follows from :eq:`eq_dirac_field_defn_j` that

.. math::
	:nowrap:

	\begin{equation*}
		\Jscr^{\ast}_{\mu \nu} = -\beta\Cscr \Jscr_{\mu \nu} \Cscr^{-1}\beta
	\end{equation*}

and hence from :eq:`eq_dirac_field_linearize_representation` that

.. math::
	:nowrap:

	\begin{equation}
		D^{\ast}(L(p)) = 1 - \frac{\ifrak}{2} \omega^{\mu \nu} \Jscr^{\ast}_{\mu \nu} \
			= 1 + \frac{\ifrak}{2} \omega^{\mu \nu} \beta\Cscr \Jscr_{\mu \nu} \Cscr^{-1}\beta \
			= \beta \Cscr D(L(p)) \Cscr^{-1} \beta
		\label{eq_dirac_field_d_matrix_conjugation_1}
	\end{equation}

Using the explicit formula :math:`\eqref{eq_dirac_field_defn_c_matrix}` for :math:`\Cscr` as well :math:`\eqref{eq_dirac_field_u_spinor_zero_momentum}` and :math:`\eqref{eq_dirac_field_v_spinor_zero_momentum}` for :math:`u(0, \sigma)` and :math:`v(0, \sigma)`, we get the following

.. math::
	:nowrap:

	\begin{alignat}{2}
		u^{\ast}(\pbf, \sigma) &= \sqrt{m/p_0} D^{\ast}(L(p)) u(0, \sigma) &&= -\beta\Cscr v(\pbf, \sigma)
		\label{eq_dirac_field_u_conjugate_to_v} \\
		v^{\ast}(\pbf, \sigma) &= \sqrt{m/p_0} D^{\ast}(L(p)) v(0, \sigma) &&= -\beta\Cscr u(\pbf, \sigma)
		\label{eq_dirac_field_v_conjugate_to_u}
	\end{alignat}

These relations turns out to be useful for the charge conjugation transformation, but not for the time inversion because the spinors :math:`u` and :math:`v` are swapped. To remedy this, we notice from :math:`\eqref{eq_dirac_field_u_spinor_zero_momentum}` and :math:`\eqref{eq_dirac_field_v_spinor_zero_momentum}` that the :math:`u` and :math:`v` spinors at zero-momentum are related by :math:`\gamma_5` defined by :eq:`eq_dirac_field_defn_gamma_5`. Moreover, in order to cancel the :math:`\beta` matrices at the two ends of right-hand-side of :math:`\eqref{eq_dirac_field_d_matrix_conjugation_1}`, we can replace :math:`\pbf` with :math:`-\pbf`, which is also desirable as far as the time inversion is concerned in light of :math:`\eqref{eq_creation_operator_time_inversion_conjugation_massive}`. Putting all these considerations together, let's try the following

.. math::
	:nowrap:

	\begin{align*}
		D^{\ast}(L(-\pbf)) &= D^{\ast}(\Pcal L(p) \Pcal) \\
			&= \beta D^{\ast}(L(p)) \beta \\
			&= \gamma_5 \beta D^{\ast}(L(p)) \beta \gamma_5 \\
			&= \gamma_5 \Cscr D(L(p)) \Cscr^{-1} \gamma_5
	\end{align*}

where the third equality holds because of the Clifford relations, namely, :math:`\gamma_5` commutes with :math:`\Jscr_{\mu \nu}`, and hence :math:`D^{\ast}(L(p))`, and anti-commutes with :math:`\beta`. Now instead of :math:`\eqref{eq_dirac_field_u_conjugate_to_v}`, we can calculate as follows

.. math::
	:nowrap:

	\begin{align}
		u^{\ast}(-\pbf, \sigma) &= \sqrt{m/p_0} D^{\ast}(L(-\pbf)) u(0, \sigma)
		\label{eq_dirac_field_u_conjugate_to_u} \\
			&= \sqrt{m/p_0} \gamma_5 \Cscr D^{\ast}(L(p)) \Cscr^{-1} \gamma_5 u(0, \sigma) \nonumber \\
			&= \sqrt{m/p_0} \gamma_5 \Cscr D^{\ast}(L(p)) (-1)^{1/2 + \sigma} u(0, -\sigma) \nonumber \\
			&= (-1)^{1/2 + \sigma} \gamma_5 \Cscr u(\pbf, -\sigma) \nonumber
	\end{align}

A similar calculation can be done to show that in fact :math:`v(-\pbf, \sigma)` satisfies exactly the same conjugation formula.

With all the preparations above, we can now calculate the spatial inversion transformation laws as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Tcal) \psi(x) U^{-1}(\Tcal) &\xlongequal{\eqref{eq_dirac_field_psi_field}} (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~\big( e^{-\ifrak p \cdot x} u^{\ast}(\pbf, \sigma) U(\Tcal) a(\pbf, \sigma) U^{-1}(\Tcal) \\
	 			&\phantom{\eqref{eq_dirac_field_psi_field}} + e^{\ifrak p \cdot x} v^{\ast}(\pbf, \sigma) U(\Tcal) a^{c \dagger}(\pbf, \sigma) U^{-1}(\Tcal) \big) \\
			&\xlongequal{\eqref{eq_creation_operator_time_inversion_conjugation_massive}} (2\pi)^{-3/2} \sum_{\sigma} (-1)^{1/2 - \sigma} \int d^3 p~\big( e^{-\ifrak p \cdot x} \zeta^{\ast} u^{\ast}(\pbf, \sigma) a(-\pbf, -\sigma) \phantom{\big)} \\
				&\phantom{\eqref{eq_dirac_field_psi_field} \big(} + e^{\ifrak p \cdot x} \zeta^c v^{\ast}(\pbf, \sigma) a^{c \dagger}(-\pbf, -\sigma) \big) \\
			&\xlongequal{\phantom{\eqref{eq_dirac_field_psi_field}}} -(2\pi)^{-3/2} \sum_{\sigma} (-1)^{1/2 + \sigma} \int d^3 p~\big( e^{-\ifrak p \cdot \Pcal x} \zeta^{\ast} u^{\ast}(-\pbf, -\sigma) a(\pbf, \sigma) \phantom{\big)} \\
			&\phantom{\eqref{eq_dirac_field_psi_field} \big(} + e^{\ifrak p \cdot \Pcal x} \zeta^c v^{\ast}(-\pbf, -\sigma) a^{c \dagger}(\pbf, \sigma) \big) \\
			&\xlongequal{\phantom{\eqref{eq_dirac_field_psi_field}}} (2\pi)^{-3/2} \gamma_5 \Cscr \sum_{\sigma} \int d^3 p~\left( \zeta^{\ast} e^{-\ifrak p \cdot \Pcal x} u(\pbf, \sigma) a(\pbf, \sigma) + \zeta^c e^{\ifrak p \cdot \Pcal x} v(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
	\end{align*}

In order for :math:`\psi(x)` to transform nicely under the time inversion, we're forced to make the following assumption

.. math::
	:nowrap:

	\begin{equation}
		\zeta^{\ast} = \zeta^c
		\label{eq_dirac_field_time_parity_relation}
	\end{equation}

Under this assumption we've finally worked out the time inversion transformation law for Dirac fields

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \psi(x) U^{-1}(\Tcal) = \zeta^{\ast} \gamma_5 \Cscr \psi(-\Pcal x)
	\end{equation*}

Next let's calculate the charge inversion transformation as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Ccal) \psi(x) U^{-1}(\Ccal) &\xlongequal{\eqref{eq_creation_operator_charge_inversion_conjugation}} (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~\left( e^{\ifrak p \cdot x} u(\pbf, \sigma) \xi^{\ast} a^c(\pbf, \sigma) + e^{-\ifrak p \cdot x} v(\pbf, \sigma) \xi^c a^{\dagger}(\pbf, \sigma) \right) \\
			&\xlongequal{\substack{\eqref{eq_dirac_field_u_conjugate_to_v} \\ \eqref{eq_dirac_field_v_conjugate_to_u}}} (2\pi)^{-3/2} \Cscr \beta \sum_{\sigma} \int d^3p~\left( \xi^{\ast} e^{\ifrak p \cdot x} v^{\ast}(\pbf, \sigma) a^c(\pbf, \sigma) + \xi^c e^{-\ifrak p \cdot x} u^{\ast}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma) \right)
	\end{align*}

Just as for the time inversion, we are forced to assuming the following condition on the charge conjugation parities

.. math::
	:nowrap:

	\begin{equation}
		\xi^{\ast} = \xi^c
		\label{eq_dirac_field_charge_parity_relation}
	\end{equation}

Under this assumption, we can work out the charge conjugation transformation law as follows [#charge_inversion_on_dirac_fields_sign]_

.. math::
	:nowrap:

	\begin{equation*}
		U(\Ccal) \psi(x) U^{-1}(\Ccal) = \xi^{\ast} \Cscr \beta \psi^{\ast}(x)
	\end{equation*}

Here we've used :math:`\psi^{\ast}(x)` instead of :math:`\psi^{\dagger}(x)` because we don't want to transpose the spinors, but it should be understood that the :math:`\ast` when applied to the creation/annihilation operators are the same as taking the adjoint.

Finally, let's consider the special case where the spin-:math:`1/2` particles are their own antiparticles. These particles are known as Majorana fermions, as already discussed in :ref:`Parities of elementary particles <dropdown_parities_of_elementary_particles>`. According to :math:`\eqref{eq_dirac_field_spatial_parity_relation}, \eqref{eq_dirac_field_time_parity_relation}` and :math:`\eqref{eq_dirac_field_charge_parity_relation}`, we see that the spatial parity of a Majorana fermion must be :math:`\pm \ifrak`, while the time and charge parity must be :math:`\pm 1`.

.. _sec_construction_of_the_interaction_density:

Construction of the interaction density
+++++++++++++++++++++++++++++++++++++++

As mentioned in :ref:`Dirac representation and gamma matrices <paragraph_dirac_field_representation_not_unitary>`, the fact that the Dirac representation is not unitary means that we cannot construct the interaction density using :math:`\psi^{\dagger} \psi` because it won't be a scalar. Indeed, let's work out how :math:`\psi^{\dagger}` transforms under a (homogeneous orthochronous) Lorentz transformation using :eq:`eq_dirac_field_pseudo_unitarity_of_d_matrix` as follows

.. math::
	:nowrap:

	\begin{equation*}
		U_0(\Lambda) \psi^{\dagger}(x) U_0^{-1}(\Lambda) = \left( U_0(\Lambda) \psi(x) U_0^{-1}(\Lambda) \right)^{\dagger}
			= \left( D^{-1}(\Lambda) \psi(\Lambda x) \right)^{\dagger}
			= \psi^{\dagger}(\Lambda x) \left( D^{-1}(\Lambda) \right)^{\dagger}
			= \psi^{\dagger}(\Lambda x) \beta D(\Lambda) \beta
	\end{equation*}

We see that if we define a new adjoint

.. math::
	:nowrap:

	\begin{equation*}
		\bar{\psi} \coloneqq \psi^{\dagger} \beta
	\end{equation*}

then :math:`\bar{\psi}` transforms nicely as follows

.. math::
	:nowrap:

	\begin{equation*}
		U_0(\Lambda) \bar{\psi}(x) U_0^{-1}(\Lambda) = \bar{\psi}(\Lambda x) D(\Lambda)
	\end{equation*}

It follows that we can construct a bilinear form as follows

.. math::
	:nowrap:

	\begin{equation*}
		\bar{\psi}(x) M \psi(x)
	\end{equation*}

where :math:`M` is a :math:`4 \times 4` matrix, so that

.. math::
	:nowrap:

	\begin{equation*}
		U_0(\Lambda) \bar{\psi}(x) M \psi(x) U_0^{-1}(\Lambda) = \bar{\psi}(\Lambda x) D(\Lambda) M D^{-1}(\Lambda) \psi(\Lambda x)
	\end{equation*}

Letting :math:`M` to be :math:`1, \gamma_{\mu}, \Jscr_{\mu \nu}, \gamma_5 \gamma_{\mu}` or :math:`\gamma_5` then produces a scalar, vector, tensor, axial vector or pseudo-scalar, respectively. Here the adjectives "axial" and "pseudo-" refer to the opposite to usual parities under spatial and/or time inversion.

An important example is Fermi's theory of beta-decay, which involves an interaction density of the following form

.. math::
	:nowrap:

	\begin{equation*}
		\bar{\psi}_p \gamma^{\mu} \psi_n \bar{\psi}_e \gamma_{\mu} \psi_{\nu}
	\end{equation*}

where :math:`p, n, e, \nu` stand for proton, neutron, electron and neutrino, respectively.


.. _sec_general_fields:

General Fields
--------------

We've now seen how scalar, vector, and Dirac fields can be constructed out of specific representations of the (homogeneous orthochronous) Lorentz group. These constructions can be generalized and unified by understanding the general representation theory of the Lorentz group.

General representation theory of the Lorentz group
++++++++++++++++++++++++++++++++++++++++++++++++++

The starting point, as in the case of Dirac fields, is the general commutation relation :eq:`eq_bracket_repr_j` that the :math:`\Jscr_{\mu \nu}` matrices must satisfy. As explained in :ref:`sec_quantum_lorentz_symmetry`, we can rename the :math:`\Jscr_{\mu \nu}` matrices as follows

.. math::
	:nowrap:

	\begin{alignat*}{3}
		\Jscr_1 &\coloneqq \Jscr_{23}, \quad \Jscr_2 &&\coloneqq \Jscr_{31}, \quad \Jscr_3 &&\coloneqq \Jscr_{12} \\
		\Kscr_1 &\coloneqq \Jscr_{10}, \quad \Kscr_2 &&\coloneqq \Jscr_{20}, \quad \Kscr_3 &&\coloneqq \Jscr_{30}
	\end{alignat*}

and rewrite :eq:`eq_bracket_repr_j` as a set of equations as follows

.. math::
	:nowrap:

	\begin{align}
		\left[ \Jscr_i, \Jscr_j \right] &= \ifrak \epsilon_{ijk} \Jscr_k
		\label{eq_jjj_commutation_general_repr} \\
		\left[ \Jscr_i, \Kscr_j \right] &= \ifrak \epsilon_{ijk} \Kscr_k
		\label{eq_jkk_commutation_general_repr} \\
		\left[ \Kscr_i, \Kscr_j \right] &= -\ifrak \epsilon_{ijk} \Jscr_k
		\label{eq_kkj_commutation_general_repr}
	\end{align}

which correspond to :math:`\eqref{eq_jjj_commutation}, \eqref{eq_jkk_commutation}` and :math:`\eqref{eq_kkj_commutation}`, respectively.

Let's write :math:`\bm{\Jscr} \coloneqq \left(\Jscr_1, \Jscr_2, \Jscr_3\right)` and :math:`\bm{\Kscr} \coloneqq \left(\Kscr_1, \Kscr_2, \Kscr_3\right)`. It turns out that this Lie algebra generated by :math:`\bm{\Jscr}` and :math:`\bm{\Kscr}` can be *complex* linearly transformed into one that splits. The transformation is defined as follows

.. math::
	:nowrap:

	\begin{align}
		\bm{\Ascr} &= \frac{1}{2} \left( \bm{\Jscr} + \ifrak \bm{\Kscr} \right)
		\label{eq_general_field_a_from_jk} \\
		\bm{\Bscr} &= \frac{1}{2} \left( \bm{\Jscr} - \ifrak \bm{\Kscr} \right)
		\label{eq_general_field_b_from_jk}
	\end{align}

so that :math:`\eqref{eq_jjj_commutation_general_repr}` -- :math:`\eqref{eq_kkj_commutation_general_repr}` take the following form

.. math::
	:nowrap:

	\begin{align}
		\left[ \bm{\Ascr}_i, \bm{\Ascr}_j \right] &= \ifrak \epsilon_{ijk} \bm{\Ascr}_k
		\label{eq_aaa_commutation} \\
		\left[ \bm{\Bscr}_i, \bm{\Bscr}_j \right] &= \ifrak \epsilon_{ijk} \bm{\Bscr}_k
		\label{eq_bbb_commutation} \\
		\left[ \bm{\Ascr}_i, \bm{\Bscr}_j \right] &= 0
		\label{eq_ab0_commutation}
	\end{align}

In other words, both :math:`\bm{\Ascr}` and :math:`\bm{\Bscr}` form a Lie algebra of the :math:`3`-dimensional rotation group and they commute each other. It follows then from :ref:`Representations of angular momentum <dropdown_repr_of_angular_momenta>` that representations of the Lie algebra defined by :math:`\eqref{eq_aaa_commutation}` -- :math:`\eqref{eq_ab0_commutation}` can be parametrized by two nonnegative (half-)integers :math:`A` and :math:`B` such that

.. math::
	:nowrap:

	\begin{align}
		\bm{\Ascr}_{a'b', ab} &= \delta_{b'b} \Jbf^{(A)}_{a'a}
		\label{eq_general_field_a_repr} \\
		\bm{\Bscr}_{a'b', ab} &= \delta_{a'a} \Jbf^{(B)}_{b'b}
		\label{eq_general_field_b_repr}
	\end{align}

where :math:`a, a' \in \{-A, -A+1, \cdots, A\}` and :math:`b, b' \in \{-B, -B+1, \cdots, B\}`, and :math:`\Jbf^{(A)}` and :math:`\Jbf^{(B)}` are matrices given by :math:`\eqref{eq_j1_j2_matrix}` and :math:`\eqref{eq_j3_matrix}`. In particular, all these representations have dimension :math:`(2A+1)(2B+1)` and are unitary.

Now each one of these representations gives rise to a representation of the Lorentz group, which will be referred to as the :math:`(A, B)` representation. As we've seen for Dirac fields, these representations are not unitary, because while :math:`\bm{\Jscr} = \bm{\Ascr} + \bm{\Bscr}` is Hermitian, :math:`\bm{\Kscr} = -\ifrak \left( \bm{\Ascr} - \bm{\Bscr} \right)` is anti-Hermitian. For the corresponding unitary representation of the :math:`3`-dimensional rotation group, we recall from :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>` that it may be split into irreducible components of spin :math:`\jfrak`, which takes values in the following range

.. math::
	:nowrap:

	\begin{equation}
		\jfrak = |A-B|, |A-B| + 1, \cdots, A+B
		\label{eq_general_field_j_range}
	\end{equation}

according to :math:`\eqref{eq_composite_total_angular_momentum_range}`. Under this setup, the scalar, vector, and Dirac fields discussed before correspond to the following three scenarios.

1. Scalar field: :math:`A = B = 0`. In this case :math:`\jfrak` must vanish, and hence no spin is possible.
2. Vector field: :math:`A = B = \tfrac{1}{2}`. In this case :math:`\jfrak` may be :math:`0` or :math:`1`, which correspond to the time and space components of the vector field, respectively.
3. Dirac field: :math:`A = \tfrac{1}{2}, B = 0` or :math:`A = 0, B = \tfrac{1}{2}`. In either case :math:`\jfrak` must be :math:`\tfrac{1}{2}`. Indeed, they correspond to the two irreducible components of (the angular momentum part of) the Dirac representation :eq:`eq_dirac_field_jscr_matrix`. Therefore the Dirac field may be written in short hand as :math:`\left( \tfrac{1}{2}, 0 \right) \oplus \left( 0, \tfrac{1}{2} \right)`.

It turns out that any general :math:`(A, B)` fields can be derived from the above basic ones by taking tensor products and irreducible components. For example :math:`(1, 0)` and :math:`(0, 1)` fields can be derived, using again :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`, by the following calculation

.. math::
	:nowrap:

	\begin{equation*}
		\left( \tfrac{1}{2}, \tfrac{1}{2} \right) \otimes \left( \tfrac{1}{2}, \tfrac{1}{2} \right) \
			= (0, 0) \oplus (0, 1) \oplus (1, 0) \oplus (1, 1)
	\end{equation*}

In fact, all :math:`(A, B)` fields with :math:`A + B` being an integer can be obtained in this way by tensoring copies of :math:`\left( \tfrac{1}{2}, \tfrac{1}{2} \right)`. To get those fields with :math:`A + B` being a half-integer, we can consider the following calculation

.. math::
	:nowrap:

	\begin{equation*}
		\left( \tfrac{1}{2}, \tfrac{1}{2} \right) \otimes \left( \left(\tfrac{1}{2}, 0\right) \oplus \left(0, \tfrac{1}{2}\right) \right) \
			= \left(0, \tfrac{1}{2}\right) \oplus \left(1, \tfrac{1}{2}\right) \oplus \left(\tfrac{1}{2}, 0\right) \oplus \left(\tfrac{1}{2}, 1\right)
	\end{equation*}


Construction of general fields
++++++++++++++++++++++++++++++

We've seen that general fields can be indexed by two (half-)integers :math:`a` and :math:`b`, and take the following general form

.. math::
	:nowrap:

	\begin{equation}
		\psi_{ab}(x) = (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~\left( \kappa e^{\ifrak p \cdot x} u_{ab}(\pbf, \sigma) a(\pbf, \sigma) + \lambda e^{-\ifrak p \cdot x} v_{ab}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
		\label{eq_general_field_defn_psi_field}
	\end{equation}

As usual, let's translate :eq:`eq_j_intertwines_u_and_v` in the context of :math:`(A, B)` representations as follows

.. math::
	:nowrap:

	\begin{align*}
		\sum_{\sigma'} u_{a'b'}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_{a, b} \bm{\Jscr}_{a'b', ab} u_{ab}(0, \sigma) \\
		-\sum_{\sigma'} v_{a'b'}(0, \sigma') \Jbf^{(\jfrak) \ast}_{\sigma' \sigma} &= \sum_{a, b} \bm{\Jscr}_{a'b', ab} v_{ab}(0, \sigma)
	\end{align*}

Using the fact that :math:`\bm{\Jscr} = \bm{\Ascr} + \bm{\Bscr}` and :math:`\eqref{eq_general_field_a_repr}` and :math:`\eqref{eq_general_field_b_repr}`, we can further rewrite these conditions as follows

.. math::
	:nowrap:

	\begin{align}
		\sum_{\sigma'} u_{a'b'}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_a \Jbf^{(A)}_{a'a} u_{ab'}(0, \sigma) + \sum_b \Jbf^{(B)}_{b'b} u_{a'b}(0, \sigma)
		\label{eq_general_field_uj_relation} \\
		-\sum_{\sigma'} v_{a'b'}(0, \sigma') \Jbf^{(\jfrak) \ast}_{\sigma' \sigma} &= \sum_a \Jbf^{(A)}_{a'a} v_{ab'}(0, \sigma) + \sum_b \Jbf^{(B)}_{b'b} v_{a'b}(0, \sigma)
		\label{eq_general_field_vj_relation}
	\end{align}

Looking at :math:`\eqref{eq_general_field_uj_relation}`, it's an identity that relates an angular momentum representation of spin :math:`\jfrak` on the left to the sum of two independent angular momentum representations on the right. Hence it's not unreasonable to guess that :math:`u` might have something to do with the Clebsch-Gordan coefficients defined by :math:`\eqref{eq_defn_clebsch_gordan_coefficients}`. This turns out to be indeed the case as we now demonstrate. Since :math:`\bigoplus_{\jfrak} \Jbf^{(\jfrak)} = \Jbf^{(A)} + \Jbf^{(B)}`, where the direct sum is taken over the range :math:`\eqref{eq_general_field_j_range}` and can be thought of as a block diagonal matrix, we can calculate as follows

.. math::
	:nowrap:

	\begin{align}
		\Jbf^{(\jfrak)}_{\sigma' \sigma} &=\left ( \Psi^{AB~\jfrak}_{\sigma'},~\bigoplus_{\jfrak'} \Jbf^{(\jfrak')} \Psi^{AB~\jfrak}_{\sigma} \right)
		\label{eq_general_field_calculate_j_in_clebsch_gordan_coefficients} \\
			&\xlongequal{\eqref{eq_defn_clebsch_gordan_coefficients}} \left( \sum_{a', b'} C^{AB}(\jfrak, \sigma'; a', b') \Psi^{AB}_{a'b'},~\sum_{a, b} C^{AB}(\jfrak, \sigma; a, b) \bigoplus_{\jfrak'} \Jbf^{(\jfrak')} \Psi^{AB}_{ab} \right)
			\nonumber \\
			&= \sum_{a', b', a, b} C^{AB}(\jfrak, \sigma'; a', b') C^{AB}(\jfrak, \sigma; a, b) \left( \Psi^{AB}_{a'b'},~\left( \Jbf^{(A)} + \Jbf^{(B)} \right) \Psi^{AB}_{ab} \right)
			\nonumber \\
			&= \sum_{a', b', a, b} C^{AB}(\jfrak, \sigma'; a', b') C^{AB}(\jfrak, \sigma; a, b) \left( \delta_{b'b} \Jbf^{(A)}_{a'a} + \delta_{a'a} \Jbf^{(B)}_{b'b} \right)
			\nonumber \\
			&= \sum_{a', b'} C^{AB}(\jfrak, \sigma'; a', b') \left( \blue{\sum_a C^{AB}(\jfrak, \sigma; a, b') \Jbf^{(A)}_{a'a} + \sum_b C^{AB}(\jfrak, \sigma; a', b) \Jbf^{(B)}_{b'b}} \right)
			\nonumber
	\end{align}

We are now one step away from being able to compare with :math:`\eqref{eq_general_field_uj_relation}`. Namely we need to bring the coefficients :math:`C^{AB}(\jfrak, \sigma'; a', b')` to the left side of the equation. To to this, we recall the following identity of Clebsch-Gordan coefficients

.. math::
	:nowrap:

	\begin{equation}
		\sum_{a', b'} C^{AB}(\jfrak', \sigma'; a', b') C^{AB}(\jfrak, \sigma; a', b') = \delta_{\jfrak'\jfrak} \delta_{\sigma'\sigma}
		\label{eq_clebsch_gordan_coefficients_orthonormal_relation}
	\end{equation}

which follows from the orthonormality of the states :math:`\Psi^{AB}_{ab}` and the reality of the Clebsch-Gordan coefficients as constructed in :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`. [#clebsch_gordan_coefficients_orthonormality]_ Now :math:`\eqref{eq_general_field_calculate_j_in_clebsch_gordan_coefficients}` will be satisfied if we set the blue terms to equal to the following quantity

.. math::
	:nowrap:

	\begin{equation}
		\sum_a C^{AB}(\jfrak, \sigma; a, b') \Jbf^{(A)}_{a'a} + \sum_b C^{AB}(\jfrak, \sigma; a', b) \Jbf^{(B)}_{b'b} \
			= \sum_{\sigma''} C^{AB}(\jfrak, \sigma; a', b') \Jbf^{(\jfrak)}_{\sigma'' \sigma}
		\label{eq_general_field_condition_on_j_by_clebsch_gordan_coefficients}
	\end{equation}

due to :math:`\eqref{eq_clebsch_gordan_coefficients_orthonormal_relation}`.

Now compare :math:`\eqref{eq_general_field_condition_on_j_by_clebsch_gordan_coefficients}` with :math:`\eqref{eq_general_field_uj_relation}`, we've solved the :math:`u`-fields of dimension :math:`(2A+1)(2B+1)` as follows

.. math::
	:nowrap:

	\begin{equation}
		u_{ab}(0, \sigma) = (2m)^{-1/2} C^{AB}(\jfrak, \sigma; a, b)
		\label{eq_general_field_u_at_zero_momentum}
	\end{equation}

where :math:`(2m)^{-1/2}` is a conventional coefficient add here to cancel the mass term in :math:`\eqref{eq_annihilation_u_transformation_simplified_by_boost}` later. Using the fact that

.. math::
	:nowrap:

	\begin{equation}
		-\Jbf^{(\jfrak) \ast}_{\sigma'\sigma} = (-1)^{\sigma'-\sigma} \Jbf^{(\jfrak)}_{-\sigma', -\sigma}
		\label{eq_angular_momentum_representation_conjugate_formula}
	\end{equation}

which can be verified directly using :math:`\eqref{eq_j1_j2_matrix}` and :math:`\eqref{eq_j3_matrix}`, we can express the :math:`v`-fields in terms of the :math:`u`-fields as follows

.. math::
	:nowrap:

	\begin{equation}
		v_{ab}(0, \sigma) = (-1)^{\jfrak + \sigma} u_{ab}(0, -\sigma)
		\label{eq_general_field_v_at_zero_momentum}
	\end{equation}

To get the :math:`u` and :math:`v` fields at finite momentum, we need to invoke the general boost formulae :math:`\eqref{eq_annihilation_u_transformation_simplified_by_boost}` -- :math:`\eqref{eq_creation_v_transformation_simplified_by_boost}`, as well as the :math:`L` transformation :math:`\eqref{eq_L_transformation_for_massive_1}` -- :math:`\eqref{eq_L_transformation_for_massive_3}`. Here we'll think of a boost as a :math:`1`-parameter transformation in a given direction :math:`\hat{\pbf} \coloneqq \pbf / |\pbf|`. It turns out to be neat to use a `hyperbolic angle <https://en.wikipedia.org/wiki/Hyperbolic_functions>`_ :math:`\theta`, rather than :math:`|\pbf|`, defined by

.. math::
	:nowrap:

	\begin{equation}
		\cosh\theta = \sqrt{\pbf^2 + m^2} / m, \quad \sinh\theta = |\pbf| / m
		\label{eq_general_field_defn_theta}
	\end{equation}

to parametrize the boost as follows

.. math::
	:nowrap:

	\begin{align*}
		L(\theta)^0_0 &= \cosh\theta \\
		L(\theta)^i_0 = L(\theta)_i^0 &= \hat{\pbf}_i \sinh\theta \\
		L(\theta)_{ij} &= \delta_{ij} + \hat{\pbf}_i \hat{\pbf}_j (\cosh\theta - 1)
	\end{align*}

The nice thing about this parametrization is that :math:`L(\theta)` becomes additive in :math:`\theta` in the following sense

.. math::
	:nowrap:

	\begin{equation*}
		L(\theta') L(\theta) = L(\theta' + \theta)
	\end{equation*}

Indeed, one can verify it by, for example, the following calculation

.. math::
	:nowrap:

	\begin{align*}
		L(\theta')^i_{\mu} L(\theta)^{\mu}_j &= L(\theta')^i_0 L(\theta)^0_j + \sum_{k=1}^3 L(\theta')^i_k L(\theta)^k_j \\
			&= \hat{\pbf}_i \hat{\pbf}_j \sinh\theta' \sinh\theta + \sum_{k=1}^3 \left( \delta_{ik} + \hat{\pbf}_i \hat{\pbf}_k (\cosh\theta' - 1) \right) \left( \delta_{kj} + \hat{\pbf}_k \hat{\pbf}_j (\cosh\theta - 1) \right) \\
			&= \hat{\pbf}_i \hat{\pbf}_j \sinh\theta' \sinh\theta + \delta_{ij} + \hat{\pbf}_i \hat{\pbf}_j \left( \cosh\theta' - 1 + \cosh\theta - 1 + (\cosh\theta' - 1)(\cosh\theta - 1) \right) \\
			&= \hat{\pbf}_i \hat{\pbf}_j \sinh\theta' \sinh\theta + \delta_{ij} + \hat{\pbf}_i \hat{\pbf}_j (\cosh\theta' \cosh\theta - 1) \\
			&= \delta_{ij} + \hat{\pbf}_i \hat{\pbf}_j (\cosh(\theta' + \theta) - 1) = L(\theta' + \theta)^i_j
	\end{align*}

In light of :eq:`eq_dirac_field_linearize_representation`, we can then write

.. math::
	:nowrap:

	\begin{equation*}
		D(L(p)) = \exp\left(-\ifrak \theta~\hat{\pbf} \cdot \bm{\Kscr}\right)
	\end{equation*}

at least for :math:`\theta` infinitesimal. Here the minus sign comes from the fact that we have to bring the upper-index :math:`\mu=0` in :math:`\omega^{\mu\nu}` down.

For an :math:`(A, B)` representation, we can further write :math:`\ifrak \bm{\Kscr} = \bm{\Ascr} - \bm{\Bscr}`, and henceforth

.. math::
	:nowrap:

	\begin{equation}
		D(L(p))_{a'b', ab} = \exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(A)}\right)_{a'a} \exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(B)}\right)_{b'b}
		\label{eq_general_field_d_transformation_in_ab_repr}
	\end{equation}

since the representation splits into a direct sum of :math:`\bm{\Ascr}` and :math:`\bm{\Bscr}`. Combining :math:`\eqref{eq_general_field_d_transformation_in_ab_repr}` with :math:`\eqref{eq_general_field_u_at_zero_momentum}` and :math:`\eqref{eq_annihilation_u_transformation_simplified_by_boost}`, we obtain the following formula for the :math:`u`-field at finite momentum for an :math:`(A, B)` representation

.. math::
	:nowrap:

	\begin{equation}
		u_{ab}(\pbf, \sigma) = \frac{1}{\sqrt{2p_0}} \sum_{a'b'} \left(\exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(A)}\right)\right)_{a'a} \left(\exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(B)}\right)\right)_{b'b} C^{AB}(\jfrak \sigma; a'b')
		\label{eq_general_field_u_at_finite_momentum}
	\end{equation}

where we assume implicitly that the spin :math:`\jfrak` is within the range :math:`\eqref{eq_general_field_j_range}`, and :math:`\sigma` is the corresponding spin :math:`z`-component.

Parallel to :math:`\eqref{eq_general_field_v_at_zero_momentum}`, we can express the :math:`v`-field at finite momentum in terms of the :math:`u`-field as follows

.. math::
	:nowrap:

	\begin{equation}
		v_{ab}(\pbf, \sigma) = (-1)^{\jfrak + \sigma} u_{ab}(\pbf, -\sigma)
		\label{eq_general_field_v_at_finite_momentum}
	\end{equation}

The construction of interaction densities for general :math:`(A, B)` fields relies on Clebsch-Gordan coefficients, and is discussed in some detail in the following dropdown block.

.. dropdown:: Construction of interaction densities
	:icon: unlock
	:animate: fade-in-slide-down

	According to :eq:`eq_construct_interaction_density_by_fields` and :eq:`eq_coefficient_g_transformation_law`, a general interaction density can be constructed as follows

	.. math::
		:nowrap:

		\begin{equation}
			\Hscr(x) = \sum_{a_1 a_2 \cdots a_n} \sum_{b_1 b_2 \cdots b_n} g_{a_1 a_2 \cdots a_n;~b_1 b_2 \cdots b_n} \psi^{(1)}_{a_1 b_1}(x) \psi^{(2)}_{a_2 b_2}(x) \cdots \psi^{(n)}_{a_n b_n}(x)
			\label{eq_general_field_interaction_density}
		\end{equation}

	where :math:`\psi^{(i)}_{a_i b_i}(x)` is an :math:`(A_i, B_i)` field, and the coefficients :math:`g_{\underline{a}, \underline{b}} \coloneqq g_{a_1 a_2 \cdots a_n;~b_1 b_2 \cdots b_n}` are covariant under the product of the :math:`(A_i, B_i)` representations. Looking at :math:`\eqref{eq_general_field_u_at_finite_momentum}` and :math:`\eqref{eq_general_field_v_at_finite_momentum}`, we see that the :math:`D` matrices, as in :math:`\eqref{eq_conjugate_annihilation_field}` -- :math:`\eqref{eq_conjugate_creation_field}`, act on general fields as a product of the angular momentum representations associated with the :math:`A`'s and :math:`B`'s. Therefore we may also split the coefficients :math:`g_{\underline{a}, \underline{b}}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			g_{a_1 a_2 \cdots a_n;~b_1 b_2 \cdots b_n} = g_{a_1 a_2 \cdots a_n} g_{b_1 b_2 \cdots b_n}
		\end{equation*}

	such that :math:`g_{\underline{a}}` and :math:`g_{\underline{b}}` are covariant under the angular momentum representations in the sense that for any :math:`3`-vector :math:`\bm{\theta}` and :math:`J \coloneqq \bm{\theta} \cdot \Jbf`, the following holds

	.. math::
		:nowrap:

		\begin{equation}
			J^{(A_1)}_{a'_1 a_1} J^{(A_2)}_{a'_2 a_2} \cdots J^{(A_n)}_{a'_n a_n} g_{a_1 a_2 \cdots a_n} = g_{a'_1 a'_2 \cdots a'_n}
			\label{eq_general_field_g_coefficients_covariance}
		\end{equation}

	for :math:`g_{\underline{a}}`, and a similar relation holds for :math:`g_{\underline{b}}`.

	Now a particularly neat set of solutions to :math:`\eqref{eq_general_field_g_coefficients_covariance}` is given by identifying :math:`g` with the Clebsch-Gordan coefficients, which we think of as the coefficients of expressing a state with definite total angular momentum in terms of states with definite individual angular momenta. In this setting, the solutions to :math:`\eqref{eq_general_field_g_coefficients_covariance}` correspond to states with zero total angular momenta.

	A particularly interesting example of this kind is Wigner's `3j-symbol <https://en.wikipedia.org/wiki/3-j_symbol>`__ defined as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\begin{pmatrix} A_1 & A_2 & A_3 \\ a_1 & a_2 & a_3 \end{pmatrix} \coloneqq \sum_{a'_3} C^{A_3 A_3}(0 0; a'_3 a_3) C^{A_1 A_3}(A_3 a'_3; a_1 a_2)
		\end{equation*}

	We read the definition as a two-step process. Namely, we first (linearly) combine states :math:`\Psi^{A_1}_{a_1}` and :math:`\Psi^{A_2}_{a_2}` to a state with total angular momentum :math:`A_3` (and spin :math:`z`-component :math:`a'_3`), and then (linearly) combine with the state :math:`\Psi^{A_3}_{a_3}` to end up in a spinless final state. [#clebsch_gordan_coefficient_zero_total_angular_momentum]_

We will now turn to the arguably most interesting causality condition :math:`\eqref{eq_h_commutativity_for_space_like_separations}`. Indeed, it is this condition that clarifies the correlation between the spin and whether a particle/field is bosonic or fermionic. As before, we need to evaluate the (anti-)commutator between the fields using :math:`\eqref{eq_general_field_defn_psi_field}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi_{ab}(x), \psi^{\prime\, \dagger}_{a'b'}(y) \right]_{\pm} = (2\pi)^{-3} \int~\frac{d^3 p}{2p_0}~\pi_{ab,a'b'}(\pbf) \left( \kappa \kappa'^{~\ast} e^{\ifrak p \cdot (x-y}) \pm \lambda \lambda'^{~\ast} e^{-\ifrak p \cdot (x-y)} \right)
		\label{eq_general_field_psi_commutator}
	\end{equation}

where :math:`\pi(\pbf)` is the (rescaled) spin sum defined by

.. math::
	:nowrap:

	\begin{equation}
		(2p_0)^{-1} \pi_{ab,a'b'}(\pbf) \coloneqq \sum_{\sigma} u_{ab}(\pbf, \sigma) u_{a'b'}^{\prime~\ast}(\pbf, \sigma) = \sum_{\sigma} v_{ab}(\pbf, \sigma) v_{a'b'}^{\prime~\ast}(\pbf, \sigma)
		\label{eq_general_field_spin_sums_as_pi}
	\end{equation}

Here the second equality can be mostly easily seen using :math:`\eqref{eq_general_field_v_at_finite_momentum}`. Note also that we are considering the general scenario where :math:`\psi(x)` is an :math:`(A, B)` field, while :math:`\psi'(x)` is a possibly different :math:`(A', B')` field.

Using :math:`\eqref{eq_general_field_u_at_finite_momentum}`, we can spell out more details of the spin sum as follows

.. math::
	:nowrap:

	\begin{align}
		\pi_{ab, a'b'}(\pbf) &= \sum_{\bar{a}~\bar{b}} \sum_{\bar{a}'~\bar{b}'} \sum_{\sigma} C^{AB}(\jfrak \sigma; \bar{a} \bar{b}) C^{A'B'}(\jfrak \sigma; \bar{a}' \bar{b}')
		\label{eq_general_field_spin_sum} \\
			&\phantom{=} \times \left(\exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(A)}\right)\right)_{\bar{a}a} \left(\exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(B)}\right)\right)_{\bar{b}b} \nonumber \\
			&\phantom{=} \times \left(\exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(A')}\right)\right)_{\bar{a}'a'} \left(\exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(B')}\right)\right)_{\bar{b}'b'} \nonumber
	\end{align}

This looks horribly complicated, but it has been evaluated by the author in [Wei69]_. Without going into the actual calculations, we note the following two facts, which suffice our purposes. The first is that :math:`\pi_{ab,a'b'}(\pbf)` is a polynomial :math:`P` in :math:`p` on the mass shell as follows

.. math::
	:nowrap:

	\begin{equation}
		\pi_{ab,a'b'}(\pbf) = P_{ab,a'b'}\left( \sqrt{\pbf^2 + m^2}, \pbf \right)
		\label{eq_general_field_spin_sum_is_polynomial}
	\end{equation}

The second is that this polynomial is even or odd depending on the parity of :math:`2A + 2B'` as follows

.. math::
	:nowrap:

	\begin{equation}
		P_{ab,a'b'}(-p) = (-1)^{2A+2B'} P_{ab,a'b'}(p)
		\label{eq_general_field_spin_sum_polynomial_parity}
	\end{equation}

.. dropdown:: Evaluation of the spin sum
	:icon: unlock
	:animate: fade-in-slide-down

	We will evaluate :math:`\pi_{ab, a'b'}(\pbf)` in the simplest case where :math:`\pbf` is along the :math:`z`-axis, so that :math:`\hat{\pbf} \cdot \Jbf` is diagonal by :math:`\eqref{eq_j3_matrix}`. In this case we have the following

	.. math::
		:nowrap:

		\begin{equation}
			\pi_{ab, a'b'}(\pbf) = \sum_{\sigma} C^{AB}(\jfrak \sigma; ab) C^{A'B'}(\jfrak \sigma; a'b') \exp((-a+b-a'+b') \theta)
			\label{eq_general_field_spin_sum_along_z_axis}
		\end{equation}

	Since the Clebsch-Gordan coefficients vanish unless

	.. math::
		:nowrap:

		\begin{equation*}
			a + b = a' + b' = \sigma
		\end{equation*}

	we can eliminate :math:`b` and :math:`a'` from the exponential in :math:`\eqref{eq_general_field_spin_sum_along_z_axis}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			-a+b-a'+b' = -a+(\sigma - a)-(\sigma - b')+b' = 2b' - 2a
		\end{equation*}

	Moreover, recall from the definition of :math:`\theta` in :math:`\eqref{eq_general_field_defn_theta}` that

	.. math::
		:nowrap:

		\begin{equation*}
			\exp(\pm\theta) = \cosh\theta \pm \sinh\theta = \left( p_0 \pm p_3 \right) / m
		\end{equation*}

	where :math:`p_0 = \sqrt{\pbf^2 + m^2}`. Hence we can rewrite :math:`\eqref{eq_general_field_spin_sum_along_z_axis}` in a polynomial in :math:`p` in two cases as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\pi_{ab, a'b'}(\pbf) = \sum_{\sigma} C^{AB}(\jfrak \sigma; ab) C^{A'B'}(\jfrak \sigma; a'b') \times \!
				\begin{cases}
					\left( (p_0 + p_3) / m \right)^{2b'-2a} & \text{if } a \leq b' \\
					\left( (p_0 - p_3) / m \right)^{2a-2b'} & \text{if } a \geq b'
				\end{cases}
		\end{equation*}

	Finally, to verify :math:`\eqref{eq_general_field_spin_sum_polynomial_parity}`, it suffices to note that :math:`2a-2b'` differs from :math:`2A-2B'` by an even integer.

Assuming :math:`\eqref{eq_general_field_spin_sum_is_polynomial}`, we note that any :math:`P_{ab, a'b'}` can be written in such a way that it's (at most) linear in the first argument :math:`\sqrt{\pbf^2 + m^2}`. Changing the content of :math:`P_{ab, a'b'}` in :math:`\eqref{eq_general_field_spin_sum_is_polynomial}`, we may then rewrite it as follows

.. math::
	:nowrap:

	\begin{equation}
		\pi_{ab, a'b'}(\pbf) = P_{ab, a'b'}(\pbf) + 2 \sqrt{\pbf^2 + m^2} Q_{ab, a'b'}(\pbf)
		\label{eq_general_field_spin_sum_as_polynomial}
	\end{equation}

where :math:`P, Q` are polynomials in :math:`\pbf` that satisfy the following parity conditions

.. math::
	:nowrap:

	\begin{align*}
		P_{ab, a'b'}(-\pbf) &= (-1)^{2A+2B'} P_{ab, a'b'}(\pbf) \\
		Q_{ab, a'b'}(-\pbf) &= -(-1)^{2A+2B'} Q_{ab, a'b'}(\pbf)
	\end{align*}

Returning to the causality condition :math:`\eqref{eq_general_field_psi_commutator}`, let's consider space separated :math:`x` and :math:`y`. Up to a Lorentz transformation, we may assume that :math:`x-y = (0, \xbf-\ybf)`. Under this assumption, we can calculate as follows

.. math::
	:nowrap:

	\begin{align*}
		\left[ \psi_{ab}(x), \psi^{\prime~\dagger}_{a'b'}(y) \right]_{\pm} &= (2\pi)^{-3} \int \frac{d^3 p}{2p_0}~P_{ab, a'b'}(\pbf) \left( \kappa \kappa'^{~\ast} e^{\ifrak p \cdot (x-y)} \pm \lambda \lambda'^{~\ast} e^{-\ifrak p \cdot (x-y)} \right) \\
			&\phantom{=} + (2\pi)^{-3} \int d^3p~Q_{ab, a'b'}(\pbf) \left( \kappa \kappa'^{~\ast} e^{\ifrak p \cdot (x-y)} \pm \lambda \lambda'^{~\ast} e^{-\ifrak p \cdot (x-y)} \right) \\
			&= \kappa\kappa^{\prime~\ast} P_{ab, a'b'}(-\ifrak \nabla) \Delta_+(\xbf - \ybf) \pm \lambda\lambda^{\prime~\ast} P_{ab, a'b'}(\ifrak \nabla) \Delta_+(\ybf - \xbf) \\
			&\phantom{=} + \kappa\kappa^{\prime~\ast} Q_{ab, a'b'}(-\ifrak \nabla) \delta^3(\xbf - \ybf) \pm \lambda\lambda^{\prime~\ast} Q_{ab, a'b'}(\ifrak \nabla) \delta^3(\ybf - \xbf) \\
			&= \left( \kappa\kappa^{\prime~\ast} \pm (-1)^{2A+2B'} \lambda\lambda^{\prime~\ast} \right) P_{ab, a'b'}(-\ifrak \nabla) \Delta_+(\xbf - \ybf) \\
			&\phantom{=} + \left( \kappa\kappa^{\prime~\ast} \mp (-1)^{2A+2B'} \lambda\lambda^{\prime~\ast} \right) Q_{ab, a'b'}(-\ifrak \nabla) \delta^3(\xbf - \ybf)
	\end{align*}

where the derivative :math:`\nabla` is always taken with respect to :math:`x`. Here we've also used the fact that :math:`\Delta_+(x)` (for space-like :math:`x`) and the Dirac delta :math:`\delta^3(x)` are even functions. We see that for the (anti-)commutator to vanish for :math:`\xbf \neq \ybf`, i.e., when :math:`\delta^3(\xbf - \ybf) = 0`, we must have

.. math::
	:nowrap:

	\begin{equation}
		\kappa\kappa^{\prime~\ast} = \mp (-1)^{2A+2B'} \lambda\lambda^{\prime~\ast}
		\label{eq_general_field_causality_kappa_lambda_condition}
	\end{equation}

Now consider an important special case where :math:`\psi = \psi'`. It implies in particular that :math:`(A, B) = (A', B')` and :math:`(\kappa, \lambda) = (\kappa', \lambda')`. In this case we can rewrite :math:`\eqref{eq_general_field_causality_kappa_lambda_condition}` as follows

.. math::
	:nowrap:

	\begin{equation}
		|\kappa|^2 = \mp (-1)^{2A+2B} |\lambda|^2 = \mp (-1)^{2\jfrak} |\lambda|^2
		\label{eq_general_field_self_kappa_lambda_relation}
	\end{equation}

since :math:`\jfrak` differs from :math:`A+B` by an integer according to :math:`\eqref{eq_general_field_j_range}`. Hence in addition to the condition :math:`|\kappa| = |\lambda|`, the field (or rather the particle it describes) must be bosonic, i.e., the bottom sign is taken, if :math:`\jfrak` is an integer, and fermionic, i.e., the top sign is taken, if :math:`\jfrak` is a half-integer. This is consistent with the corresponding conclusions for scalar , vector, and Dirac fields found in previous sections, and is indeed a great clarification of the relationship between spin and statistics, e.g., `Pauli's exclusion principle <https://en.wikipedia.org/wiki/Pauli_exclusion_principle>`__.

Back to the general case. We know from :math:`\eqref{eq_general_field_self_kappa_lambda_relation}` that :math:`|\kappa'| = |\lambda'|` and :math:`(-1)^{2A+2B} = (-1)^{2\jfrak} = \mp`, which is the same sign as in :math:`\eqref{eq_general_field_causality_kappa_lambda_condition}`. Hence we can rewrite :math:`\eqref{eq_general_field_causality_kappa_lambda_condition}` by dividing both sides by :math:`|\kappa'|^2 = |\lambda'|^2` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\frac{\kappa}{\kappa'} = (-1)^{2B+2B'} \frac{\lambda}{\lambda'} \
		\implies (-1)^{2B} \frac{\kappa}{\lambda} = (-1)^{2B'} \frac{\kappa'}{\lambda'}
	\end{equation*}

Hence we conclude the following relationship between the coefficients :math:`\kappa` and :math:`\lambda`

.. math::
	:nowrap:

	\begin{equation*}
		\lambda = (-1)^{2B} c \kappa
	\end{equation*}

where :math:`c` is a constant that depends only on the field, or rather, the particle it describes, and not on the specific representation that gives rise to the field. Moreover we note that :math:`c` is a phase since :math:`|\kappa| = |\lambda|`. Hence by adjusting the phase of the creation operator (and correspondingly the annihilation operator), we can arrange so that :math:`c = 1`.

This marks the end of the discussion about the causality condition on general :math:`(A, B)`. As a result, we've obtained the following grand formula for a general (causal) field.

.. math::
	:nowrap:

	\begin{equation}
		\psi_{ab}(x) = (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \left( e^{\ifrak p \cdot x} u_{ab}(\pbf, \sigma) a(\pbf, \sigma) + \
			(-1)^{2B} e^{-\ifrak p \cdot x} v_{ab}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
		\label{eq_general_field_psi_field}
	\end{equation}

where the spinors :math:`u_{ab}` and :math:`v_{ab}` are given by :math:`\eqref{eq_general_field_u_at_finite_momentum}` and :math:`\eqref{eq_general_field_v_at_finite_momentum}`, respectively.

.. dropdown:: Same field given by different representations are physically indifferent
	:icon: unlock
	:animate: fade-in-slide-down

	As we've seen in :ref:`sec_scalar_field` and :ref:`sec_spin_zero_vector_field`, a spin-:math:`0` field can arise either as a scalar field, i.e., a :math:`(0, 0)` field, or from a vector field, i.e., a :math:`\left( \tfrac{1}{2}, \tfrac{1}{2} \right)` field. Moreover, the later turns out to be the (first) derivatives of the former, and hence doesn't produce anything really new. We'll see here that this is a rather general phenomenon.

	Indeed, according to :math:`\eqref{eq_general_field_j_range}`, a spin-:math:`0` field can only arise from an :math:`(A, A)` field. Given a scalar field :math:`\psi(x)`, one can construct a spin-:math:`0` :math:`(A, A)` field as follows

	.. math::
		:nowrap:

		\begin{equation}
			\left\{ \p_{\mu_1} \p_{\mu_2} \cdots \p_{\mu_{2A}} \right\} \psi(x)
			\label{eq_general_field_spin_zero_aa_field}
		\end{equation}

	where :math:`\left\{ \p_{\mu_1} \p_{\mu_2} \cdots \p_{\mu_{2A}} \right\}` is a traceless (symmetric) product of the partial derivatives. Here the trace of a symmetric tensor :math:`S_{\mu_1 \mu_2 \cdots \mu_n}` is defined to be

	.. math::
		:nowrap:

		\begin{equation}
			\op{tr}(S_{\mu_1 \mu_2 \cdots \mu_n}) \coloneqq \sum_{\mu_1, \mu_2 = 0}^3 \delta_{\mu_1 \mu_2} S_{\mu_1 \mu_2 \cdots \mu_n}
			\label{eq_trace_of_symmetric_tensor}
		\end{equation}

	We'll not verify that :math:`\eqref{eq_general_field_spin_zero_aa_field}` indeed transforms as a :math:`(A, A)` field, but we'll verify that it at least has the expected dimension :math:`(2A+1)^2`. To this end, we first note that a rank :math:`2A` symmetric tensor (in :math:`4` dimensions) has dimension

	.. math::
		:nowrap:

		\begin{equation*}
			\frac{(2A+1)(2A+2)(2A+3)}{3!}
		\end{equation*}

	which combinatorially is just the number of ways of having four nonnegative integers sum up to :math:`2A`. Next we must subtract from it the number of traces one can take on these tensors. According to :math:`\eqref{eq_trace_of_symmetric_tensor}`, each trace fixes two indexes, while the rest is still symmetric. Hence the number of traces is the same combinatorial number with :math:`2A` replaced by :math:`2A-2` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\frac{(2A-1)(2A)(2A+1)}{3!}
		\end{equation*}

	Now the difference of the two number above, which is the dimension of traceless symmetric tensors, is exactly :math:`(2A+1)^2` as expected.

	In general a spin-:math:`\jfrak` field can arise from an :math:`(A, B)` field as long as the triangle inequality

	.. math::
		:nowrap:

		\begin{equation}
			|A-B| \leq \jfrak \leq A+B
			\label{eq_spin_range_abj_triangle_inequality}
		\end{equation}

	is satisfied. We claim that the :math:`(A, B)` representation is the same as the tensor product of the :math:`(\jfrak, 0)` representation and the :math:`(B, B)` representation. Indeed, following :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`, we see that the later is a direct sum of :math:`(A, B)` representations as long as

	.. math::
		:nowrap:

		\begin{equation}
			|B-\jfrak| \leq A \leq B+\jfrak
		\end{equation}

	but this triangle inequality is exactly the same as :math:`\eqref{eq_spin_range_abj_triangle_inequality}`. Hence the claim is proved.

	By the same argument as in the spin-:math:`0` case, we conclude that any spin-:math:`\jfrak` field can be written as

	.. math::
		:nowrap:

		\begin{equation*}
			\left\{ \p_{\mu_1} \p_{\mu_2} \cdots \p_{\mu_{2B}} \right\} \psi_{\sigma}(x)
		\end{equation*}

	where :math:`\psi_{\sigma}(x)` is the :math:`(\jfrak, 0)` field. Swapping the role of :math:`A` and :math:`B`, one can also write it as rank :math:`2A` traceless derivatives of the :math:`(0, \jfrak)` field.


The CPT symmetries
++++++++++++++++++

The calculations of space, time, and charge conjugation transformations in the general case is essentially the same as for the Dirac field. In particular, instead of reverting the :math:`3`-momentum in Dirac spinors as in :math:`\eqref{eq_dirac_field_spatial_inversion_acts_on_psi_plus}`, we need to do it for general :math:`(A, B)` spinors :math:`\eqref{eq_general_field_u_at_finite_momentum}`, which involves the Clebsch-Gordan coefficients.

Without going to the details, we list the relevant symmetry properties of Clebsch-Gordan coefficients as follows

.. math::
	:nowrap:

	\begin{align}
		C^{AB}(\jfrak, \sigma; a, b) &= (-1)^{A+B-\jfrak} C^{BA}(\jfrak, \sigma; b, a)
		\label{eq_clebsch_gordan_symmetry_swap} \\
		C^{AB}(\jfrak, \sigma; a, b) &= (-1)^{A+B-\jfrak} C^{AB}(\jfrak, -\sigma; -a, -b)
		\label{eq_clebsch_gordan_symmetry_reverse}
	\end{align}

The first relation is proved in [Wei00]_ page 124, and the second relation can be deduced from the time reversal transformation law :math:`\eqref{eq_time_inversion_on_massive_general}`.

Consider first the spatial inversion. Combining :math:`\eqref{eq_clebsch_gordan_symmetry_swap}` with :math:`\eqref{eq_general_field_u_at_finite_momentum}`, one obtains the following relations on the spinors

.. math::
	:nowrap:

	\begin{align}
		u^{AB}_{ab}(-\pbf, \sigma) &= (-1)^{A+B-\jfrak} u^{BA}_{ba}(\pbf, \sigma)
		\label{eq_general_field_u_symmetry_swap} \\
		v^{AB}_{ab}(-\pbf, \sigma) &= (-1)^{A+B-\jfrak} v^{BA}_{ba}(\pbf, \sigma)
		\label{eq_general_field_v_symmetry_swap}
	\end{align}

We can calculate the spatial conjugation as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Pcal) \psi^{AB}_{ab}(x) U^{-1}(\Pcal) &\xlongequal{\eqref{eq_general_field_psi_field}} (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \left( \eta^{\ast} e^{\ifrak p \cdot x} u^{AB}_{ab}(\pbf, \sigma) a(-\pbf, \sigma) + (-1)^{2B} \eta^c e^{-\ifrak p \cdot x} v^{AB}_{ab}(\pbf, \sigma) a^{c \dagger}(-\pbf, \sigma) \right) \\
			&\xlongequal{\substack{\eqref{eq_general_field_u_symmetry_swap} \\ \eqref{eq_general_field_v_symmetry_swap}}} -(2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~(-1)^{A+B-\jfrak} \big( \blue{\eta^{\ast} e^{\ifrak p \cdot \Pcal x} u_{ba}(\pbf, \sigma) a(\pbf, \sigma)} \phantom{\big)} \\
   			&\phantom{\xlongequal{\eqref{eq_general_field_psi_field}}\big(}~\blue{+ (-1)^{2B} \eta^c e^{-\ifrak p \cdot \Pcal x} v_{ba}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma)} \big)
	\end{align*}

In order for the blue terms to be proportional to the corresponding terms in, in this case, a :math:`(B, A)` field, we must have

.. math::
	:nowrap:

	\begin{equation}
		(-1)^{2B} \eta^c = (-1)^{2A} \eta^{\ast} \iff \eta^c = (-1)^{2\jfrak} \eta^{\ast}
		\label{eq_general_field_space_inversion_parity_relation}
	\end{equation}

Under this assumption, we can complete the transformation law for spatial inversion as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \psi^{AB}_{ab}(x) U^{-1}(\Pcal) = -\eta^{\ast} (-1)^{A+B-\jfrak} \psi^{BA}_{ba}(\Pcal x)
		\label{eq_general_field_space_inversion_transformation}
	\end{equation}

which recovers the cases of scalar field :eq:`eq_scalar_field_spatial_inversion_transformation_law`, vector field :eq:`eq_vector_field_spatial_inversion_transformation_law`, and Dirac field :math:`\eqref{eq_dirac_field_spatial_inversion_transformation_law}` where :math:`\beta`, as defined by :eq:`eq_dirac_field_beta_matrix`, serves the function of swapping :math:`A` and :math:`B`.

Next consider the time inversion. As for the spatial inversion, we'll need the following identities

.. math::
	:nowrap:

	\begin{align}
		u^{AB \ast}_{ab}(-\pbf, -\sigma) &\xlongequal{\eqref{eq_general_field_u_at_finite_momentum}} \frac{1}{\sqrt{2p_0}} \sum_{a'b'} \left(\exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(A) \ast}\right)\right)_{-a', a} \left(\exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(B) \ast}\right)\right)_{-b', b} C^{AB}(\jfrak, -\sigma; -a', -b')
			\label{eq_general_field_u_symmetry_negation} \\
			&\xlongequal{\eqref{eq_clebsch_gordan_symmetry_reverse}} \frac{(-1)^{A+B-\jfrak}}{\sqrt{2p_0}} \sum_{a'b'} \left(\exp\left(\theta~\hat{\pbf} \cdot \Jbf^{(A) \ast}\right)\right)_{-a', a} \left(\exp\left(-\theta~\hat{\pbf} \cdot \Jbf^{(B) \ast}\right)\right)_{-b', b} C^{AB}(\jfrak, \sigma; a', b') \nonumber \\
			&\xlongequal{\substack{\eqref{eq_angular_momentum_representation_conjugate_formula} \\ \eqref{eq_angular_momentum_representation_conjugate_formula_as_matrix}}} \frac{(-1)^{A+B-\jfrak}}{\sqrt{2p_0}} \sum_{a'b'} (-1)^{a-a'+b-b'} \left( \exp(-\theta~\hat{\pbf} \cdot \Jbf^{(A)}) \right)_{a', -a} \left( \exp(\theta~\hat{\pbf} \cdot \Jbf^{(B)}) \right)_{b', -b} C^{AB}(\jfrak, \sigma; a', b') \nonumber \\
			&= (-1)^{A+B+a+b-\sigma - \jfrak} u^{AB}_{-a, -b}(\pbf, \sigma) \nonumber
	\end{align}

where it's convenient for the third equality to reformulate :math:`\eqref{eq_angular_momentum_representation_conjugate_formula}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\Jbf^{(\jfrak) \ast} = -C \Jbf^{(\jfrak)} C^{-1}, \quad \text{where}~~C_{\sigma' \sigma} = (-1)^{\jfrak - \sigma} \delta_{\sigma', -\sigma}
		\label{eq_angular_momentum_representation_conjugate_formula_as_matrix}
	\end{equation}

Since :math:`v^{AB}_{ab}(\pbf, \sigma)` is related to :math:`u^{AB}_{ab}(\pbf, \sigma)` by :math:`\eqref{eq_general_field_v_at_finite_momentum}`, we can derive the :math:`v`-counterpart of :math:`\eqref{eq_general_field_u_symmetry_negation}` as follows

.. math::
	:nowrap:

	\begin{align}
		v^{AB \ast}_{ab}(-\pbf, -\sigma) &\xlongequal{\eqref{eq_general_field_v_at_finite_momentum}} (-1)^{\jfrak - \sigma} u^{AB \ast}_{ab}(-\pbf, \sigma)
			\label{eq_general_field_v_symmetry_negation} \\
			&\xlongequal{\eqref{eq_general_field_u_symmetry_negation}} (-1)^{A+B+a+b} u^{AB}_{-a, -b}(\pbf, -\sigma) \nonumber \\
			&\xlongequal{\eqref{eq_general_field_v_at_finite_momentum}} (-1)^{A+B+a+b-\sigma - \jfrak} v^{AB}_{-a, -b}(\pbf, \sigma) \nonumber
	\end{align}

Remembering that :math:`U(\Tcal)` is anti-unitary, we can calculate the time inversion transformation as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Tcal) \psi^{AB}_{ab}(x) U^{-1}(\Tcal) &\xlongequal{\substack{\eqref{eq_general_field_psi_field} \\ \eqref{eq_creation_operator_time_inversion_conjugation_massive}}} (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~(-1)^{\jfrak - \sigma} \big( \zeta^{\ast} e^{-\ifrak p \cdot x} u^{AB \ast}_{ab}(\pbf, \sigma) a(-\pbf, -\sigma) \phantom{\big)}
		\label{eq_general_field_time_inversion_transformation} \\
	  		&\phantom{\xlongequal{\eqref{eq_general_field_psi_field}}\big(} + (-1)^{2B} \zeta^c e^{\ifrak p \cdot x} v^{AB \ast}_{ab}(\pbf, \sigma) a^{c \dagger}(-\pbf, -\sigma) \big) \nonumber \\
			&\xlongequal{\phantom{\eqref{eq_general_field_u_symmetry_negation}}} -(2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~(-1)^{\jfrak - \sigma} \big( \zeta^{\ast} e^{-\ifrak p \cdot \Pcal x} u^{AB \ast}_{ab}(-\pbf, -\sigma) a(\pbf, \sigma) \phantom{\big)} \nonumber \\
			&\phantom{\xlongequal{\eqref{eq_general_field_u_symmetry_negation}}\big(} + (-1)^{2B} \zeta^c e^{\ifrak p \cdot \Pcal x} v^{AB \ast}_{ab}(-\pbf, -\sigma) a^{c \dagger}(\pbf, \sigma) \big) \nonumber \\
			&\xlongequal{\substack{\eqref{eq_general_field_u_symmetry_negation} \\ \eqref{eq_general_field_v_symmetry_negation}}} (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~(-1)^{A+B+a+b-2\sigma} \big( \zeta^{\ast} e^{-\ifrak p \cdot \Pcal x} u^{AB}_{-a, -b}(\pbf, \sigma) \phantom{\big)} \nonumber \\
			&\phantom{\xlongequal{\eqref{eq_general_field_u_symmetry_negation}}\big(} + (-1)^{2B} \zeta^c e^{\ifrak p \cdot \Pcal x} v^{AB}_{-a, -b}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \big) \nonumber \\
			&\xlongequal{\phantom{\eqref{eq_general_field_u_symmetry_negation}}} (-1)^{A+B+a+b-2\jfrak} \zeta^{\ast} \psi^{AB}_{-a, -b}(-\Pcal x) \nonumber
	\end{align}

where the last equality assumes the following symmetry on the time-reversal parity

.. math::
	:nowrap:

	\begin{equation}
		\zeta^{\ast} = \zeta^c
		\label{eq_general_field_time_inversion_parity_relation}
	\end{equation}

At this point, we're pretty proficient at (and tired of) this kind of calculation. Hence we'll not spell out the (rather similar) details for the charge conjugation symmetry, but rather list the result as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Ccal) \psi^{AB}_{ab}(x) U^{-1}(\Ccal) = (-1)^{-2A-a-b-\jfrak} \xi^{\ast} \psi^{BA \dagger}_{-b, -a}(x)
		\label{eq_general_field_charge_inversion_transformation}
	\end{equation}

under the following assumption on the charge-reversal parity

.. math::
	:nowrap:

	\begin{equation}
		\xi^{\ast} = \xi^c
		\label{eq_general_field_charge_inversion_parity_relation}
	\end{equation}


The CPT Theorem
---------------

With all the hard work we've done in :ref:`sec_general_fields`, we can now reward ourselves a bit with the celebrated `CPT theorem <https://en.wikipedia.org/wiki/CPT_symmetry>`__ which is stated as follows

	For an appropriate choice of the inversion phases :math:`\eta` (space), :math:`\zeta` (time), and :math:`\xi` (charge), the product :math:`U(CPT)` is conserved.

We'll skip over the special case of scalar, vector, and Dirac fields, and jump directly into the general, and in fact simpler, case of :math:`(A, B)` fields.

.. math::
	:nowrap:

	\begin{align}
		U(CPT) \psi^{AB}_{ab}(x) U^{-1}(CPT) &\xlongequal{\eqref{eq_general_field_time_inversion_transformation}} (-1)^{A+B+a+b-2\jfrak} \zeta^{\ast} U(CP) \psi^{AB}_{-a, -b}(-\Pcal x) U^{-1}(CP)
		\label{eq_cpt_conjugation_general_field_calculation} \\
			&\xlongequal{\eqref{eq_general_field_space_inversion_transformation}} -(-1)^{a+b-\jfrak} \zeta^{\ast} \eta^{\ast} U(C) \psi^{BA}_{-b, -a}(-x) U^{-1}(C) \nonumber \\
			&\xlongequal{\eqref{eq_general_field_charge_inversion_transformation}} -(-1)^{-2B} \zeta^{\ast} \eta^{\ast} \xi^{\ast} \psi^{AB \dagger}_{ab}(-x) \nonumber
	\end{align}

Hence if we assume the following condition on the inversion parities

.. admonition:: Assumption on the inversion parities
	:class: Important

	.. math::
		:nowrap:

		\begin{equation}
			\zeta~\eta~\xi = 1
			\label{eq_cpt_parities_product_assumption}
		\end{equation}

then we can rewrite :math:`\eqref{eq_cpt_conjugation_general_field_calculation}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		U(CPT) \psi^{AB}_{ab}(x) U^{-1}(CPT) = -(-1)^{2B} \psi^{AB \dagger}_{ab}(-x)
	\end{equation*}

A few words are needed, however, to justify the seemingly strange assumption on the product of inversion parities. Indeed, it is physically meaningless to specify any inversion parity for a single species of particles because it's just a phase. The only conditions that we've seen on the parities come from pairs of particles and their antiparticles, notably :math:`\eqref{eq_general_field_space_inversion_parity_relation}, \eqref{eq_general_field_time_inversion_parity_relation}`, and :math:`\eqref{eq_general_field_charge_inversion_parity_relation}`. We saw that the time and charge inversion parities are the same between the particle and its antiparticle, respectively. However, a sign :math:`(-1)^{2\jfrak}` is involved in the space inversion parity. So if we impose :math:`\eqref{eq_cpt_parities_product_assumption}` on one particle species, then it will fail on its antiparticle species if the particle in question is a fermion! We're eventually saved by the fact that the interaction density must involve an even number of fermions (cf. discussions in :ref:`sec_causality_and_antiparticles`). In any case :math:`\eqref{eq_cpt_parities_product_assumption}` is a fairly sloppy assumption, which cannot hold in general, but it also doesn't make a difference to the CPT theorem.

Now suppose the interaction density :math:`\Hscr(x)` is defined by :math:`\eqref{eq_general_field_interaction_density}` as a linear combination of monomials like

.. math::
	:nowrap:

	\begin{equation*}
		\psi^{A_1 B_1}_{a_1 b_1}(x) \psi^{A_2 B_2}_{a_2 b_2}(x) \cdots \psi^{A_n B_n}_{a_n b_n}(x)
	\end{equation*}

Hence in light of :math:`\eqref{eq_general_field_g_coefficients_covariance}`, we know that both :math:`A_1 + A_2 + \cdots + A_n` and :math:`B_1 + B_2 + \cdots + B_n` must be integers, for otherwise they cannot be coupled to a spinless state. It follows then the following CPT transformation law on the interaction density

.. math::
	:nowrap:

	\begin{equation*}
		U(CPT) \Hscr(x) U^{-1}(CPT) = -\Hscr(-x)
	\end{equation*}

Recall from :math:`\eqref{eq_evolution_equation_of_u_operator}` and :math:`\eqref{eq_defn_v_by_density}` that the interaction term :math:`V = \int d^3 x~\Hscr(0, \xbf)` satisfies the following

.. math::
	:nowrap:

	\begin{equation*}
		U(CPT) V U^{-1}(CPT) = -\int d^3 x~\Hscr(-x) = V
	\end{equation*}

Since the CPT symmetry is clearly conserved for free particles, it is also conserved in interactions according to :math:`\eqref{eq_h_as_h0_plus_v}`.


Massless Fields
---------------

So far the story about quantum fields has been a 100% success. We've namely found the general formula :math:`\eqref{eq_general_field_psi_field}` for *any* field that represents a massive particle. However, such success will come to an end when we consider instead massless particles as we'll see in this section. This should not come as a surprise though since we've see in :eq:`eq_vector_field_defn_Pi` for example, that the spin sum blows up in the massless limit :math:`m \to 0`.

Let's nonetheless kickstart the routine of constructing fields as follows, and see where the problem should arise.

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\ell}(x) = (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \left( \kappa e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma) a(\pbf, \sigma) + \lambda e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
		\label{eq_massless_field_defn_psi_field}
	\end{equation}

This is reasonable because the translation symmetry is the same for massive and massless particles, and hence :eq:`eq_redefine_u_after_translation` and :eq:`eq_redefine_v_after_translation` apply.

Next, using the general transformation laws :math:`\eqref{eq_lorentz_transformation_formula_for_creation_operator}` and :math:`\eqref{eq_lorentz_transformation_formula_for_annihilation_operator}` for creation and annihilation operators, as well as the :math:`D` matrix :math:`\eqref{eq_little_group_d_matrix_massless}` for massless particles, we can infer the homogeneous Lorentz transformation laws as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda) a^{\dagger}(\pbf, \sigma) U^{-1}(\Lambda) &= \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(\ifrak \sigma \theta(\Lambda, p)) a^{\dagger}(\pbf_{\Lambda}, \sigma)
		\label{eq_massless_vector_field_creation_operator_conjugated_by_lorentz_transformation} \\
		U(\Lambda) a(\pbf, \sigma) U^{-1}(\Lambda) &= \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(-\ifrak \sigma \theta(\Lambda, p)) a(\pbf_{\Lambda}, \sigma)
		\label{eq_massless_vector_field_annihilation_operator_conjugated_by_lorentz_transformation}
	\end{align}

Just as in the massive case, we'd like :math:`\psi_{\ell}(x)` to satisfy the following transformation law

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \psi_{\ell}(x) U^{-1}(\Lambda) = \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}(\Lambda x)
		\label{eq_massless_field_psi_transform_by_d_matrix}
	\end{equation}

To see what conditions the spinors must satisfy, let's first expand the left-hand-side as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Lambda) \psi_{\ell}(x) U^{-1}(\Lambda) &= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \Big( \kappa e^{\ifrak p \cdot x} \blue{u_{\ell}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(-\ifrak \sigma \theta(\Lambda, p))} a(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + \lambda e^{-\ifrak p \cdot x} \blue{v_{\ell}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(\ifrak \sigma \theta(\Lambda, p))} a^{c \dagger}(\pbf_{\Lambda}, \sigma) \Big)
	\end{align*}

Then the right-hand-side as follows

.. math::
	:nowrap:

	\begin{align*}
		\sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}(\Lambda x) &= (2\pi)^{-3/2} \sum_{\sigma} \int \frac{d^3 p}{p_0}~p_0 \sum_{\ell'} \Big( \kappa e^{\ifrak p \cdot \Lambda x} D_{\ell \ell'}(\Lambda^{-1}) u_{\ell'}(\pbf, \sigma) a(\pbf, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + \lambda e^{-\ifrak p \cdot \Lambda x} D_{\ell \ell'}(\Lambda^{-1}) v_{\ell'}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \Big) \\
			&= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~\blue{\frac{(\Lambda p)_0}{p_0} \sum_{\ell'}} \Big( \kappa e^{\ifrak p \cdot x} \blue{D_{\ell \ell'}(\Lambda^{-1}) u_{\ell'}(\pbf_{\Lambda}, \sigma)} a(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + \lambda e^{-\ifrak p \cdot x} \blue{D_{\ell \ell'}(\Lambda^{-1}) v_{\ell'}(\pbf_{\Lambda}, \sigma)} a^{c \dagger}(\pbf_{\Lambda}, \sigma) \Big)
	\end{align*}

Equating the coefficients of :math:`a(\pbf_{\Lambda, \sigma})` and :math:`a^{c \dagger}(\pbf_{\Lambda}, \sigma)` (i.e., the blue terms), and inverting :math:`D_{\ell \ell'}(\Lambda^{-1})` as in :eq:`eq_annihilation_u_transformation`, we get the following conditions on the spinors

.. math::
	:nowrap:

	\begin{align}
		\exp(\ifrak \sigma \theta(\Lambda, p)) u_{\ell'}(\pbf_{\Lambda}, \sigma) &= \sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) u_{\ell}(\pbf, \sigma)
		\label{eq_massless_field_spinor_u_condition} \\
		\exp(-\ifrak \sigma \theta(\Lambda, p)) v_{\ell'}(\pbf_{\Lambda}, \sigma) &= \sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) v_{\ell}(\pbf, \sigma)
		\label{eq_massless_field_spinor_v_condition}
	\end{align}

The next step is to take the massless analogy to the boost operator in the massive case. Namely, if we let :math:`\Lambda = L(p)` be the chosen Lorentz transformation that takes the standard :math:`k = (1, 0, 0, 1)` to :math:`p`, then :math:`\theta(\Lambda, p) = 0`. Taking :math:`p = k` in :math:`\eqref{eq_massless_field_spinor_u_condition}` and :math:`\eqref{eq_massless_field_spinor_v_condition}`, we obtain the following

.. math::
	:nowrap:

	\begin{align}
		u_{\ell'}(\pbf, \sigma) &= \frac{1}{\sqrt{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) u_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_u_from_k} \\
		v_{\ell'}(\pbf, \sigma) &= \frac{1}{\sqrt{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) v_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_v_from_k}
	\end{align}

Next, in analogy to the rotation transformation, let's consider a little group element :math:`W` that fixes :math:`k`. In this case :math:`\eqref{eq_massless_field_spinor_u_condition}` and :math:`\eqref{eq_massless_field_spinor_v_condition}` take the following form

.. math::
	:nowrap:

	\begin{align}
		\exp(\ifrak \sigma \theta(W, k)) u_{\ell'}(\kbf, \sigma) &= \sum_{\ell} D_{\ell' \ell}(W) u_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_little_group_u_condition} \\
		\exp(-\ifrak \sigma \theta(W, k)) v_{\ell'}(\kbf, \sigma) &= \sum_{\ell} D_{\ell' \ell}(W) v_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_little_group_v_condition}
	\end{align}

Recall from :ref:`sec_massless_particle_states` that any :math:`W` can be written in the following form

.. math::
	:nowrap:

	\begin{equation}
		W(a, b, \theta) = S(a, b) R(\theta)
		\label{eq_massless_vector_field_w_equals_s_times_r}
	\end{equation}

where :math:`S(a, b)` is defined by :math:`\eqref{eq_massless_little_group_s_matrix}`, and :math:`R(\theta)` is defined by :math:`\eqref{eq_massless_little_group_r_matrix}`. Considering separately the two cases :math:`W(0, 0, \theta) = R(\theta)` and :math:`W(a, b, 0) = S(a, b)`, we get the following two consequences of :math:`\eqref{eq_massless_field_spinor_u_condition}`

.. math::
	:nowrap:

	\begin{align}
		e^{\ifrak \sigma \theta} u_{\ell'}(\kbf, \sigma) &= \sum_{\ell} D_{\ell' \ell}(R(\theta)) u_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_spinor_u_condition_from_r} \\
		u_{\ell'}(\kbf, \sigma) &= \sum_{\ell} D_{\ell' \ell}(S(a, b)) u_{\ell}(\kbf, \sigma)
		\label{eq_massless_field_spinor_u_condition_from_s}
	\end{align}

Similar constraints hold for :math:`v` as well, but we'll not bother to write them down.

It turns out, however, that these conditions can never be satisfied! To illustrate the difficulties, we'll first consider the case of vector fields, both as a warm-up and for later references when we'll analyze the electromagnetic theory. Then we'll show that the difficulties persist to the general case of arbitrary :math:`(A, B)` fields.

The failure for vector fields
+++++++++++++++++++++++++++++

For vector field :math:`D_{\mu}^{\nu}(\Lambda) = \Lambda_{\mu}^{\nu}` as in the massive case. As a convention, let's write

.. math::
	:nowrap:

	\begin{equation*}
		u_{\mu}(\pbf, \sigma) \eqqcolon \frac{1}{\sqrt{2p_0}} e_{\mu}(\pbf, \sigma)
	\end{equation*}

Since :math:`D_{\mu}^{\nu}(\Lambda)` is real, it follows from :math:`\eqref{eq_massless_field_spinor_u_condition}` and :math:`\eqref{eq_massless_field_spinor_v_condition}` that :math:`v` satisfies equations that are complex conjugate to those that :math:`u` satisfies. Hence :math:`v_{\mu}(\pbf, \sigma) = u^{\ast}_{\mu}(\pbf, \sigma)`.

Now we can translate :math:`\eqref{eq_massless_field_u_from_k}` from a boosting formula for :math:`u` to one for :math:`e` as follows

.. math::
	:nowrap:

	\begin{equation}
		e_{\mu}(\pbf, \sigma) = L(p)_{\mu}^{\nu} e_{\nu}(\kbf, \sigma)
		\label{eq_massless_vector_field_spinor_from_k_to_p}
	\end{equation}

Moreover :math:`\eqref{eq_massless_field_spinor_u_condition_from_r}` and :math:`\eqref{eq_massless_field_spinor_u_condition_from_s}` can be translated to conditions on :math:`e` as well as follows

.. math::
	:nowrap:

	\begin{align}
		e^{\ifrak \sigma \theta} e_{\mu}(\kbf, \sigma) &= R(\theta)_{\mu}^{\nu} e_{\nu}(\kbf, \sigma)
		\label{eq_massless_field_spinor_e_condition_from_r} \\
		e_{\mu}(\kbf, \sigma) &= S(a, b)_{\mu}^{\nu} e_{\nu}(\kbf, \sigma)
		\label{eq_massless_field_spinor_e_condition_from_s}
	\end{align}

Using the explicit formula :math:`\eqref{eq_massless_little_group_r_matrix}` for :math:`R(\theta)`, we can derive from :math:`\eqref{eq_massless_field_spinor_e_condition_from_r}` that the helicity :math:`\sigma = \pm 1`, and moreover,

.. math::
	:nowrap:

	\begin{equation}
		e_{\mu}(\kbf, \pm 1) = \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ \pm \ifrak \\ 0 \end{bmatrix*}
		\label{eq_massless_vector_field_e_at_k}
	\end{equation}

up to normalization. However, by the explicit formula :math:`\eqref{eq_massless_little_group_s_matrix}` for :math:`S(a, b)`, we see that :math:`\eqref{eq_massless_field_spinor_e_condition_from_s}` then requires

.. math::
	:nowrap:

	\begin{equation*}
		a \pm \ifrak b = 0
	\end{equation*}

which is impossible for any real :math:`a, b` that are not both zero.

For reasons that will be justified later, it's nonetheless legitimate to adopt the vectors :math:`e_{\mu}` as defined by :math:`\eqref{eq_massless_vector_field_e_at_k}` as the spinors, as well as the condition :math:`\kappa = \lambda = 1` as in the case of massive vector fields. With these assumptions, we can rename :math:`\psi` by :math:`a` (as it'll correspond to the `electromagnetic potential <https://en.wikipedia.org/wiki/Electromagnetic_four-potential>`__ which is conventionally named by :math:`a`), and rewrite :math:`\eqref{eq_massless_field_defn_psi_field}` as follows

.. math::
	:nowrap:

	\begin{equation}
		a_{\mu}(x) = (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \left( e^{\ifrak p \cdot x} e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) + e^{-\ifrak p \cdot x} e^{\ast}_{\mu}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
		\label{eq_massless_vector_field_a}
	\end{equation}

.. warning::

	The notations are getting slightly out of hands here. Namely, we've used :math:`a` for at least three different things in one place: the vector field :math:`a_{\mu}`, the parameter in :math:`S(a, b)`, and the creation operator :math:`a(\pbf, \sigma)`. There will actually be a fourth place where :math:`a` is used as the spin :math:`z`-component in an :math:`(A, B)` field. We can only hope that the context will make it clear what :math:`a` (or :math:`b`) really represents.

As for massive vector fields, we'll search for field equations that :math:`a_{\mu}(x)` must satisfy. First of all, it satisfies obviously the (massless) Klein-Gordon equation

.. math::
	:nowrap:

	\begin{equation}
		\square a_{\mu}(x) = 0
		\label{eq_massless_vector_field_klein_gordon}
	\end{equation}

which is nothing but an incarnation of the mass-shell condition :math:`p_0^2 = \pbf^2`. Then let's consider the massless analog to the gauge-fixing condition :eq:`eq_vector_field_gauge_fixing_condition`. To this end, we claim that :math:`e_0(\kbf, \pm 1) = 0` and :math:`\kbf \cdot \ebf(\kbf, \pm 1) = 0` imply the following

.. math::
	:nowrap:

	\begin{align}
		e_0(\pbf, \pm 1) &= 0
		\label{eq_massless_vector_field_spinor_zero_vanishes} \\
		\pbf \cdot \ebf(\pbf, \pm 1) &= 0
		\label{eq_massless_vector_field_spinor_orthogonal_to_momentum}
	\end{align}

in analogy to :eq:`eq_vector_field_spinor_orthogonal_to_momentum` by the following argument. First, note that :math:`e_{\mu}(\pbf, \pm 1)` can be obtained from :math:`e_{\mu}(\kbf, \pm 1)` by applying :math:`L(p)` as in :math:`\eqref{eq_massless_vector_field_spinor_from_k_to_p}`. Second, :math:`L(p)` can be decomposed into a boost along the :math:`z`-axis as in :math:`\eqref{eq_massless_boost}` followed by a :math:`3`-rotation. Finally, we conclude :math:`\eqref{eq_massless_vector_field_spinor_zero_vanishes}` and :math:`\eqref{eq_massless_vector_field_spinor_orthogonal_to_momentum}` by noting that :math:`e_{\mu}(\kbf, \pm 1)` is unaffected by any boost along the :math:`z`, and the dot product is preserved by any :math:`3`-rotation. In contrast to the massless case :eq:`eq_vector_field_spinor_orthogonal_to_momentum`, we have a stronger constraint :math:`\eqref{eq_massless_vector_field_spinor_zero_vanishes}` here because the helicity :math:`0` spinor is missing.

The corresponding constraints on :math:`a_{\mu}(x)` is the following

.. math::
	:nowrap:

	\begin{align}
		a_0(x) &= 0
		\label{eq_massless_vector_field_a0_vanishes} \\
		\nabla \cdot \abf(x) &= 0
		\label{eq_massless_vector_field_a_divergence_vanishes}
	\end{align}

which is clearly not Lorentz invariant.

But let's calculate :math:`U(\Lambda) a_{\mu}(x) U^{-1}(\Lambda)` anyway and see to some extent :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}` fails. To this end, we'll need to calculate the action of the :math:`D` matrix on the spinors, and we'll first do this for :math:`e_{\mu}(k, \pm 1)` as follows

.. math::
	:nowrap:

	\begin{align}
		{D_{\mu}}^{\nu}(W(a, b, \theta)) e_{\nu}(\kbf, \pm 1)
			&\xlongequal{\eqref{eq_massless_vector_field_w_equals_s_times_r}} {S(a, b)_{\mu}}^{\lambda} {R(\theta)_{\lambda}}^{\nu} e_{\nu}(\kbf, \pm 1)
			\label{eq_massless_vector_field_dw_acts_on_spinor} \\
			&\xlongequal{\eqref{eq_massless_field_spinor_e_condition_from_r}} e^{\pm \ifrak \theta} {S(a, b)_{\mu}}^{\lambda} e_{\lambda}(\kbf, \pm 1)
			\nonumber \\
			&\xlongequal{\eqref{eq_massless_little_group_s_matrix}} e^{\pm \ifrak \theta} \left( e_{\mu}(\kbf, \pm 1) + \frac{a \pm \ifrak b}{\sqrt{2}} k_{\mu} \right)
			\nonumber
	\end{align}

Next we recall the little group element :math:`W(\Lambda, p) = L^{-1}(\Lambda p) \Lambda L(p)` by definition. Plugging into :math:`\eqref{eq_massless_vector_field_dw_acts_on_spinor}`, we obtain the following

.. math::
	:nowrap:

	\begin{equation}
		{\Lambda_{\mu}}^{\nu} e_{\nu}(\pbf, \pm 1) = e^{\pm \ifrak \theta(\Lambda, p)} \left(e_{\mu}(\pbf_{\Lambda}, \pm 1) + (\Lambda p)_{\mu} \Omega_{\pm}(\Lambda, p)\right)
		\label{eq_massless_vector_field_spinor_e_condition}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{equation*}
		\Omega_{\pm}(\Lambda, p) \coloneqq \frac{a(\Lambda, p) \pm \ifrak b(\Lambda, p)}{\sqrt{2}}
	\end{equation*}

is the extra term that makes it different from :math:`\eqref{eq_massless_field_spinor_u_condition}` which would have been satisfied if :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}` holds.

To most conveniently utilize :math:`\eqref{eq_massless_vector_field_spinor_e_condition}`, let's calculate :math:`{D_{\mu}}^{\nu} \left( U(\Lambda) a_{\nu}(x) U^{-1}(\Lambda) \right)` using :math:`\eqref{eq_massless_vector_field_a}, \eqref{eq_massless_vector_field_creation_operator_conjugated_by_lorentz_transformation}`, and :math:`\eqref{eq_massless_vector_field_annihilation_operator_conjugated_by_lorentz_transformation}` as follows

.. math::
	:nowrap:

	\begin{align*}
		{D_{\mu}}^{\nu}(\Lambda) \left(U(\Lambda) a_{\nu}(x) U^{-1}(\Lambda)\right)
			&= (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot x} {\Lambda_{\mu}}^{\nu} e_{\nu}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} e^{-\ifrak \sigma \theta(\Lambda, p)} a(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} e^{-\ifrak p \cdot x} {\Lambda_{\mu}}^{\nu} e^{\ast}_{\nu}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} e^{\ifrak \sigma \theta} a^{c \dagger}(\pbf_{\Lambda}, \sigma) \Big) \\
			&\xlongequal{\eqref{eq_massless_vector_field_spinor_e_condition}} (2\pi)^{-3/2} \int \frac{d^3 p}{p_0} \sqrt{\frac{(\Lambda p)_0}{2}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot x} e_{\mu}(\pbf_{\Lambda}, \sigma) a(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{\xlongequal{\eqref{eq_massless_vector_field_spinor_e_condition}} \Big(} + e^{-\ifrak p \cdot x} e^{\ast}_{\mu}(\pbf_{\Lambda}, \sigma) a^{c \dagger}(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{\xlongequal{\eqref{eq_massless_vector_field_spinor_e_condition}} \Big(} + (\Lambda p)_{\mu} \left( e^{\ifrak p \cdot x} \Omega_{\sigma}(\Lambda, p) a(\pbf_{\Lambda}, \sigma) + e^{-\ifrak p \cdot x} \Omega^{\ast}_{\sigma}(\Lambda, p) a^{c \dagger}(\pbf_{\Lambda}, \sigma) \right) \Big) \\
			&= (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot \Lambda x} e_{\mu}(\pbf, \pm 1) a(\pbf, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + e^{-\ifrak p \cdot \Lambda x} e^{\ast}_{\mu}(\pbf, \pm 1) a^{c \dagger}(\pbf, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + p_{\mu} \left( e^{\ifrak p \cdot \Lambda x} \Omega_{\sigma}(\Lambda, \Lambda^{-1} p) + e^{-\ifrak p \cdot \Lambda x} \Omega^{\ast}_{\sigma}(\Lambda, \Lambda^{-1} p) a^{c \dagger}(\pbf, \sigma) \right) \Big) \\
			&= a_{\mu}(\Lambda x) + \frac{\p}{\p \left( (\Lambda x)_{\mu} \right)} \Omega(\Lambda, x)
	\end{align*}

where :math:`\Omega(\Lambda, x)` is a linear combination of creation and annihilation operators, whose precise form is not important here. Finally, moving :math:`D(\Lambda)` to the right-hand-side, we obtain the following variation of :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}`

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) a_{\mu}(x) U^{-1}(\Lambda) = {D(\Lambda^{-1})_{\mu}}^{\nu} a_{\nu}(\Lambda x) + \p_{\mu} \Omega(\Lambda, x)
		\label{eq_massless_vector_field_a_conjugation_by_lorentz_transformation}
	\end{equation}

which the massless vector field :math:`\eqref{eq_massless_vector_field_a}` actually satisfies.

It follows, using integration by parts, that one can construct interaction densities by coupling :math:`j^{\mu}(x) a_{\mu}(x)` as long as :math:`\p_{\mu} j^{\mu}(x) = 0`. This is completely parallel to :eq:`eq_vector_field_j_coupling` and :eq:`eq_vector_field_j_coupling_condition` for massive vector fields, and hence partially justifies the choice of the spinors in :math:`\eqref{eq_massless_vector_field_e_at_k}`, which satisfy :math:`\eqref{eq_massless_field_spinor_e_condition_from_r}` but not :math:`\eqref{eq_massless_field_spinor_e_condition_from_s}`.

Another byproduct of :math:`\eqref{eq_massless_vector_field_a_conjugation_by_lorentz_transformation}` is the observation that although :math:`a_{\mu}` fails to be a vector, one can easily construct a :math:`2`-tensor as follows

.. math::
	:nowrap:

	\begin{equation}
		f_{\mu \nu} = \p_{\mu} a_{\nu} - \p_{\nu} a_{\mu}
		\label{eq_massless_vector_field_curvature_tensor}
	\end{equation}

which obviously satisfies :math:`U(\Lambda) f_{\mu \nu} U^{-1}(\Lambda) = D_{\mu}^{\lambda}(\Lambda^{-1}) D_{\nu}^{\sigma}(\Lambda^{-1}) f_{\lambda \sigma}`. In fact, using :math:`\eqref{eq_massless_vector_field_klein_gordon}, \eqref{eq_massless_vector_field_a0_vanishes}` and :math:`\eqref{eq_massless_vector_field_a_divergence_vanishes}`, one can show that :math:`f_{\mu \nu}` satisfies the vacuum Maxwell equations

.. math::
	:nowrap:

	\begin{align*}
		\p_{\mu} f^{\mu \nu} &= 0 \\
		\epsilon^{\rho \tau \mu \nu} \p_{\tau} f_{\mu \nu} &= 0
	\end{align*}

where :math:`\epsilon^{\rho \tau \mu \nu}` is the totally anti-symmetric sign. Indeed :math:`f_{\mu \nu}` is the quantization of the electromagnetic field, while :math:`a_{\mu}` is the quantization of the electromagnetic potential.

.. dropdown:: Build a quantum field with :math:`f_{\mu \nu}`
	:icon: unlock
	:animate: fade-in-slide-down

	Since :math:`f_{\mu \nu}` is a tensor, one might try to build a (causal) quantum field with :math:`f_{\mu \nu}`, instead of :math:`a_{\mu}`. The causality condition :math:`\eqref{eq_h_commutativity_for_space_like_separations}` demands that the (anti-)commutator :math:`\left[ f_{\mu \nu}(x), f^{\dagger}_{\rho \tau}(y) \right]_{\pm}` must vanish for space-separated :math:`x` and :math:`y`. As in the other cases, some preliminary calculations are in order. First we calculate the spin sum as follows

	.. math::
		:nowrap:

		\begin{equation}
			\sum_{\sigma = \pm 1} \ebf_i(\kbf, \sigma) \ebf_j^{\ast}(\kbf, \sigma) = \delta_{ij} - \kbf_i \kbf_j ~~\xRightarrow{\eqref{eq_massless_vector_field_spinor_from_k_to_p}} \
				\sum_{\sigma = \pm 1} \ebf_i(\pbf, \sigma) \ebf_j^{\ast}(\pbf, \sigma) = \delta_{ij} - \frac{\pbf_i \pbf_j}{\pbf^2}
			\label{eq_massless_vector_field_spin_sum}
		\end{equation}

	where the first equality can be directly verified using :math:`\eqref{eq_massless_vector_field_e_at_k}`. In particular, we notice that the spin sum is real -- a neat fact that will be used in the next calculation. Moreover, the spin sum vanishes if any index is :math:`0` because of :math:`\eqref{eq_massless_vector_field_spinor_zero_vanishes}`.

	Next we calculate the (anti-)commutator between the derivatives of the components of an :math:`a`-field, with coefficients :math:`\kappa` and :math:`\lambda` as in :math:`\eqref{eq_massless_field_defn_psi_field}` restored, as follows

	.. math::
		:nowrap:

		\begin{align*}
			\left[ \p_{\mu} a_{\nu}(x), \p_{\rho} a^{\dagger}_{\tau}(y) \right]_{\pm} &\xlongequal{\eqref{eq_massless_vector_field_a}} (2\pi)^{-3} \int \frac{d^3 p~d^3 p'}{2\sqrt{p_0 p'_0}} \sum_{\sigma, \sigma' = \pm 1} \Big( |\kappa|^2 e^{\ifrak p \cdot x - \ifrak p' \cdot y} p_{\mu} p'_{\rho} e_{\nu}(\pbf, \sigma) e^{\ast}_{\tau}(\pbf', \sigma') \left[ a(\pbf, \sigma), a^{\dagger}(\pbf', \sigma') \right]_{\pm} \phantom{\Big)} \\
				&\phantom{\xlongequal{\eqref{eq_massless_vector_field_a}} \Big(} + |\lambda|^2 e^{-\ifrak p \cdot x + \ifrak p' \cdot y} p_{\mu} p'_{\rho} e^{\ast}_{\nu}(\pbf, \sigma) e_{\tau}(\pbf', \sigma') \left[ a^{c \dagger}(\pbf, \sigma), a^c(\pbf', \sigma') \right]_{\pm} \Big) \\
				&\xlongequal{\eqref{eq_creation_annihilation_commutator}} (2\pi)^{-3} \int \frac{d^3 p}{2p_0} \sum_{\sigma = \pm 1} p_{\mu} p_{\rho} \left( |\kappa|^2 e^{\ifrak p \cdot (x-y)} e_{\nu}(\pbf, \sigma) e^{\ast}_{\tau}(\pbf, \sigma) \pm |\lambda|^2 e^{-\ifrak p \cdot (x-y)} e^{\ast}_{\nu}(\pbf, \sigma) e_{\tau}(\pbf, \sigma) \right) \\
				&\xlongequal{\eqref{eq_massless_vector_field_spin_sum}} (2\pi)^{-3} \int \frac{d^3 p}{2p_0} \blue{p_{\mu} p_{\rho} (1 - \delta_{0\nu})(1 - \delta_{0\tau}) \left( \delta_{\nu \tau} - \frac{p_{\nu} p_{\tau}}{\pbf^2} \right)} \left( |\kappa|^2 e^{\ifrak p \cdot (x-y)} \pm |\lambda|^2 e^{-\ifrak p \cdot (x-y)} \right)
		\end{align*}

	Notice that the blue terms are the only terms that involve the indexes :math:`\mu, \nu, \rho`, and :math:`\tau`. For the convenience of notations, let's call it :math:`P_{\mu \nu \rho \tau}`. Now we can calculate the (anti-)commutator as follows

	.. math::
		:nowrap:

		\begin{align}
			\left[ f_{\mu \nu}(x), f^{\dagger}_{\rho \tau}(y) \right]_{\pm} &= \left[ \p_{\mu} a_{\nu}(x) - \p_{\nu} a_{\mu}(x), \p_{\rho} a^{\dagger}_{\tau}(y) - \p_{\tau} a^{\dagger}_{\rho}(y) \right]
			\label{eq_massless_curvature_commutator_first_calculation} \\
				&= \left[ \p_{\mu} a_{\nu}(x), \p_{\rho} a^{\dagger}_{\tau}(y) \right] \
	  				- \left[ \p_{\mu} a_{\nu}(x), \p_{\tau} a^{\dagger}_{\rho}(y) \right] \
	  				- \left[ \p_{\nu} a_{\mu}(x), \p_{\rho} a^{\dagger}_{\tau}(y) \right] \
	  				+ \left[ \p_{\nu} a_{\mu}(x), \p_{\tau} a^{\dagger}_{\rho}(y) \right] \nonumber \\
				&= (2\pi)^{-3} \int \frac{d^3 p}{2p_0} \left( P_{\mu \nu \rho \tau} - P_{\mu \nu \tau \rho} - P_{\nu \mu \rho \tau} + P_{\nu \mu \tau \rho} \right) \left( |\kappa|^2 e^{\ifrak p \cdot (x-y)} \pm |\lambda|^2 e^{-\ifrak p \cdot (x-y)} \right) \nonumber
		\end{align}

	It remains to calculate the following

	.. math::
		:nowrap:

		\begin{align*}
			P_{\mu \nu \rho \tau} - P_{\mu \nu \tau \rho} - P_{\nu \mu \rho \tau} + P_{\nu \mu \tau \rho} &= p_{\mu} p_{\rho} \left( \delta_{\nu \tau} - \delta_{0 \nu}\delta_{0 \tau} - (1 - \delta_{0 \nu})(1 - \delta_{0 \tau}) \frac{p_{\nu} p_{\tau}}{\pbf^2} \right) \\
				&\phantom{=} - p_{\mu} p_{\tau} \left( \delta_{\nu \rho} - \delta_{0 \nu}\delta_{0 \rho} - (1 - \delta_{0 \nu})(1 - \delta_{0 \rho}) \frac{p_{\nu} p_{\rho}}{\pbf^2} \right) \\
				&\phantom{=} - p_{\nu} p_{\rho} \left( \delta_{\mu \tau} - \delta_{0 \mu}\delta_{0 \tau} - (1 - \delta_{0 \mu})(1 - \delta_{0 \tau}) \frac{p_{\mu} p_{\tau}}{\pbf^2} \right) \\
				&\phantom{=} + p_{\nu} p_{\tau} \left( \delta_{\mu \rho} - \delta_{0 \mu}\delta_{0 \rho} - (1 - \delta_{0 \mu})(1 - \delta_{0 \rho}) \frac{p_{\mu} p_{\rho}}{\pbf^2} \right) \\
				&= p_{\mu} p_{\rho} (\delta_{\nu \tau} - 2\delta_{0 \nu}\delta_{0 \tau}) - p_{\mu} p_{\tau} (\delta_{\nu \rho} - 2\delta_{0 \nu}\delta_{0 \rho}) \\
				&\phantom{=} - p_{\nu} p_{\rho} (\delta_{\mu \tau} - 2\delta_{0 \mu}\delta_{0 \tau}) + p_{\nu} p_{\tau} (\delta_{\mu \rho} - 2\delta_{0 \mu}\delta_{0 \rho}) \\
				&= \eta_{\nu \tau} p_{\mu} p_{\rho} - \eta_{\nu \rho} p_{\mu} p_{\tau} - \eta_{\mu \tau} p_{\nu} p_{\rho} + \eta_{\mu \rho} p_{\nu} p_{\tau}
		\end{align*}

	Finally we can finish :math:`\eqref{eq_massless_curvature_commutator_first_calculation}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\left[ f_{\mu \nu}(x), f^{\dagger}_{\rho \tau}(y) \right]_{\pm} = -(2\pi)^{-3} \left( \eta_{\nu \tau} \p_{\mu} \p_{\rho} - \eta_{\nu \rho} \p_{\mu} \p_{\tau} - \eta_{\mu \tau} \p_{\nu} \p_{\rho} + \eta_{\mu \rho} \p_{\nu} \p_{\tau} \right) \int \frac{d^3 p}{2p_0} \left( |\kappa|^2 e^{\ifrak p \cdot (x-y)} \pm |\lambda|^2 e^{-\ifrak p \cdot (x-y)} \right)
		\end{equation*}

	For this expression to vanish, we must take the bottom sign, i.e., it's a commutator. Moreover, we must have :math:`|\kappa| = |\lambda|` since :math:`\Delta_+(x-y) = \Delta_+(y-x)` if :math:`x-y` is space-like. By further rescaling the creation operator if necessary, we may assume :math:`\kappa = \lambda = 1`, which partially justifies the same choice taken for the :math:`a`-field :math:`\eqref{eq_massless_vector_field_a}`. In particular, the massless helicity :math:`\pm 1` (vector) field must be bosonic, as expected.


The failure for general fields
++++++++++++++++++++++++++++++

The problem that massless fields cannot be made to satisfy the transformation law :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}` is not specific to vector fields. Indeed, we'll show in this section that the same problem persists to all :math:`(A, B)` representations.

Recall from :math:`\eqref{eq_general_field_a_from_jk}` -- :math:`\eqref{eq_general_field_b_from_jk}`, and :math:`\eqref{eq_general_field_a_repr}` -- :math:`\eqref{eq_general_field_b_repr}`, that we can explicitly write the :math:`\Jscr_{\mu \nu}` matrix as follows

.. math::
	:nowrap:

	\begin{align}
		\left( \Jscr_{ij} \right)_{a'b', ab} &= \epsilon_{ijk} \left( \delta_{b'b} \left( \Jbf^{(A)}_k \right)_{a'a} + \delta_{a'a} \left( \Jbf^{(B)}_k \right)_{b'b} \right)
		\label{eq_massless_general_field_jab} \\
		\left( \Jscr_{k0} \right)_{a'b', ab} &= -\ifrak \left( \delta_{b'b} \left( \Jbf^{(A)}_k \right)_{a'a} - \delta_{a'a} \left( \Jbf^{(B)}_k \right)_{b'b} \right)
		\label{eq_massless_general_field_kab}
	\end{align}

It follows from :eq:`eq_dirac_field_linearize_representation`, and the fact according to :math:`\eqref{eq_massless_little_group_r_matrix}` that :math:`R(\theta)` is the rotation in the :math:`xy`-plane, that

.. math::
	:nowrap:

	\begin{equation*}
		D(R(\theta)) = 1 + \ifrak \theta \Jscr_{12}
	\end{equation*}

for :math:`\theta` infinitesimal. The linearized version of :math:`\eqref{eq_massless_field_spinor_u_condition_from_r}`, and its counterpart for :math:`v`, in this case is given by the following

.. math::
	:nowrap:

	\begin{align*}
		\sigma u_{a'b'}(\kbf, \sigma) &\xlongequal{\eqref{eq_massless_field_spinor_u_condition_from_r}} \sum_{a'b'} \left( \Jscr_{12} \right)_{a'b', ab} u_{ab}(\kbf, \sigma) \
			\xlongequal{\eqref{eq_massless_general_field_jab}} \sum_{ab} (a' + b') \delta_{a'a} \delta_{b'b} u_{ab}(\kbf, \sigma) \
			= (a' + b') u_{a'b'}(\kbf, \sigma) \\
		-\sigma v_{a'b'}(\kbf, \sigma) &= (a' + b') v_{a'b'}(\kbf, \sigma)
	\end{align*}

It follows that

.. math::
	:nowrap:

	\begin{alignat}{2}
		a' + b' &= \sigma  \quad &&\text{if}~~ u_{a'b'}(\kbf, \sigma) \neq 0
		\label{eq_massless_general_field_ab_condition_1} \\
		a' + b' &= -\sigma \quad &&\text{if}~~ v_{a'b'}(\kbf, \sigma) \neq 0
		\label{eq_massless_general_field_ab_condition_for_v}
	\end{alignat}

Next let's linearize :math:`\eqref{eq_massless_field_spinor_u_condition_from_s}` using the explicit formula :math:`\eqref{eq_massless_little_group_s_matrix}` for :math:`S(a, b)` as follows

.. math::
	:nowrap:

	\begin{align*}
		0 &= \sum_{a, b} \left( \Jscr_{01} + \Jscr_{31} \right)_{a'b', ab} u_{ab}(\kbf, \sigma) \\
			&= \sum_{a, b} \left( \ifrak \delta_{b'b} \left( \Jbf^{(A)}_1 \right)_{a'a} - \ifrak \delta_{a'a} \left( \Jbf^{(B)}_1 \right)_{b'b} + \delta_{b'b} \left( \Jbf^{(A)}_2 \right)_{a'a} + \delta_{a'a} \left( \Jbf^{(B)}_2 \right)_{b'b} \right) u_{ab}(\kbf, \sigma) \\
			&= \sum_{a} \left( \ifrak \Jbf^{(A)}_1 + \Jbf^{(A)}_2 \right)_{a'a} u_{ab'}(\kbf, \sigma) + \sum_{b} \left( -\ifrak \Jbf^{(B)}_1 + \Jbf^{(B)}_2 \right)_{b'b} u_{a'b}(\kbf, \sigma) \\
		0 &= \sum_{a, b} \left( \Jscr_{02} + \Jscr_{32} \right)_{a'b', ab} u_{ab}(\kbf, \sigma) \\
			&= \sum_{a, b} \left( \ifrak \delta_{b'b} \left( \Jbf^{(A)}_2 \right)_{a'a} - \ifrak \delta_{a'a} \left( \Jbf^{(B)}_2 \right)_{b'b} - \delta_{b'b} \left( \Jbf^{(A)}_1 \right)_{a'a} - \delta_{a'a} \left( \Jbf^{(B)}_1 \right)_{b'b} \right) u_{ab}(\kbf, \sigma) \\
			&= \sum_{a} \left( -\Jbf^{(A)}_1 + \ifrak \Jbf^{(A)}_2 \right)_{a'a} u_{ab'}(\kbf, \sigma) - \sum_{b} \left( \Jbf^{(B)}_1 + \ifrak \Jbf^{(B)}_2 \right)_{b'b} u_{a'b}(\kbf, \sigma)
	\end{align*}

which correspond to taking :math:`a` infinitesimal and :math:`b = 0`, and :math:`b` infinitesimal and :math:`a = 0`, respectively. Since :math:`u_{ab}` are linearly independent, the above constraints are equivalent to the following

.. math::
	:nowrap:

	\begin{align*}
		\sum_{a} \left( \Jbf^{(A)}_1 - \ifrak \Jbf^{(A)}_2 \right)_{a'a} u_{ab'}(\kbf, \sigma) &= 0 \\
		\sum_{b} \left( \Jbf^{(B)}_1 + \ifrak \Jbf^{(B)}_2 \right)_{b'b} u_{a'b}(\kbf, \sigma) &= 0
	\end{align*}

In order for :math:`u_{ab'}` to be annihilated by the lowering operator, and for :math:`u_{a'b}` to be annihilated by the raising operator, we must have

.. math::
	:nowrap:

	\begin{equation}
		a' = -A, \quad b' = B
		\label{eq_massless_general_field_ab_condition_2}
	\end{equation}

if :math:`u_{a'b'}(\kbf, \sigma) \neq 0`. Combining :math:`\eqref{eq_massless_general_field_ab_condition_1}` and :math:`\eqref{eq_massless_general_field_ab_condition_2}`, we conclude that for the :math:`u`-spinor to not vanish, we must have

.. math::
	:nowrap:

	\begin{equation}
		\sigma = B - A
		\label{eq_massless_general_field_helicity_condition}
	\end{equation}

It follows that a general massless :math:`(A, B)` field, according to :math:`\eqref{eq_general_field_defn_psi_field}`, can only destroy particles of helicity :math:`B-A`. Similar argument can be applied to the :math:`v`-spinor, which, together with :math:`\eqref{eq_massless_general_field_ab_condition_for_v}`, implies that the field can only create antiparticles of helicity :math:`A-B`.

As a special case, we see once again that a massless helicity :math:`\pm 1` field cannot be constructed as a vector field, i.e., a :math:`\left( \tfrac{1}{2}, \tfrac{1}{2} \right)` field, because such vector field must be scalar by :math:`\eqref{eq_massless_general_field_helicity_condition}`. Indeed, the simplest massless helicity :math:`\pm 1` field must be a :math:`(1, 0) \oplus (0, 1)` field, which is nothing but the anti-symmetric :math:`2`-tensor :math:`f_{\mu \nu}` defined by :math:`\eqref{eq_massless_vector_field_curvature_tensor}`.

.. rubric:: Footnotes

.. [#two_ways_of_representation] This should be compared with :eq:`eq_d_repr_of_little_group`, which uses transpose instead of inverse to satisfy the group law. From the viewpoint of representation theory, it's more natural to use the inverse. But in the case of little group representations, since we're only interested in unitary representations, there is not much difference between the two choices.

.. [#wrong_integration_of_Delta_function] The evaluation of the integral eq. (5.2.8) on [Wei95]_ page 202 seems to be wrong, as the integrand oscillates as :math:`\sin(u)` for :math:`u` large enough, which will cause the integral to diverge.

.. [#charge_inversion_on_dirac_fields_sign] Our calculation differs from the calculation eq. (5.5.47) in [Wei95]_ by a sign.

.. [#clebsch_gordan_coefficients_orthonormality] Details of the argument can be found in [Wei15]_ page 121 -- 122.

.. [#clebsch_gordan_coefficient_zero_total_angular_momentum] An explicit evaluation of :math:`C^{AA}(0, 0;a, -a)`, with the only non-vanishing combination of spin :math:`z`-components, can be found in [Wei15]_ page 124 -- 125.