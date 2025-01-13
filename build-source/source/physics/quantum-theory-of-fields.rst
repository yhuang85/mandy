The Quantum Theory of Fields (S. Weinberg)
==========================================

.. warning::
	This note is work in progress.

This note covers the three volumes [Wei95]_, [Wei96]_, and [Wei00]_, written by S. Weinberg on quantum field theory, with occasional references to [Wei15]_ by the same author. What I like the most about these books is his attempt to make logical deductions from the most basic principles, in particular, the principle of symmetries, rather than to make analogies to experience, e.g., from classical physics (from a historical perspective). Such an endeavor may not always be possible because, after all, physics is about how we interpret Nature based on nothing but experience, and *not* about how Nature actually works. By the same token, the arguments that are considered logical here should really be interpreted as "reasonable-sounding", and have nothing to do with what mathematician would call "rigorous".

The order of the sections correspond roughly to the order of the chapters of the books.

What is a Quantum State?
------------------------

Quantum theory postulates that *any* physical state (of the world) can be represented by a *ray* in some complex Hilbert space. It's worth noting that it is the state, rather than the Hilbert space, that we actually care about. Let's write a state as

.. math::
	:nowrap:

	\begin{equation*}
		[\Psi] \coloneqq \{ c\Psi ~|~ c \in \Cbb \setminus 0 \}
	\end{equation*}

where :math:`\Psi` is a nonzero vector in the Hilbert space. It is, however, rather inconvenient to have to deal with :math:`[\Psi]` all the time. So instead, we will almost always pick a representative :math:`\Psi`, often out of a natural choice, and call it a *state vector*, and keep in mind that anything physically meaningful must not be sensitive to a scalar multiplication.

.. admonition:: Assumption
	:class: Important

	Throughout this post we always assume that state vectors are normalized so that :math:`||\Psi|| = 1`.

In fact, we don't really care about the states themselves either, because they are more of an abstraction rather than something one can physically measure. What we do care about are the (Hermitian) inner products between state vectors, denoted by :math:`(\Psi, \Phi)`. According to the so-called `Copenhagen interpretation <https://en.wikipedia.org/wiki/Copenhagen_interpretation>`_ of quantum mechanics, such inner product represents an *amplitude*, i.e., its squared norm gives the probability of finding a state :math:`[\Psi]` in :math:`[\Phi]` if we ever perform a measurement. We can write this statement as an equation as follows

.. math::
	:nowrap:

	\begin{equation*}
		P([\Psi] \to [\Phi]) = |(\Psi, \Phi)|^2
	\end{equation*}

In particular, the probability of finding any state in itself is one, due to the normalization above.

.. _sec_what_is_a_symmetry:

What is a Symmetry?
-------------------

We start with a *symmetry transformation*, by which we mean a transformation that preserves all quantities that one can ever measure about a system. Since it is the probabilities, rather than the states themselves, that are measurable, one is led to define a quantum symmetry transformation as a transformation of states :math:`T` such that

.. math::
	:nowrap:

	\begin{equation}
		P(T[\Psi] \to T[\Phi]) = P([\Psi] \to [\Phi])
		\label{eq_t_preserves_probability}
	\end{equation}

for any states :math:`[\Psi]` and :math:`[\Phi]`. Now a theorem of E. Wigner asserts that such :math:`T` can be realized either as a linear unitary or as an anti-linear anti-unitary transformation :math:`U = U(T)` of state vectors in the sense that :math:`[U\Psi] = T[\Psi]` for any :math:`\Psi`. In other words, :math:`U` satisfies either

.. math::
	:nowrap:

	\begin{equation*}
		U(c\Psi) = cU\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)
	\end{equation*}

or

.. math::
	:nowrap:

	\begin{equation*}
		U(c\Psi) = c^{\ast} U\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)^{\ast}
	\end{equation*}

where :math:`c` is any complex number.

.. dropdown:: Proof of Wigner's theorem
	:icon: unlock
	:animate: fade-in-slide-down

	The construction of a realization :math:`U` of :math:`T` takes the following steps.

	:Step 1: Fix an orthonormal basis :math:`\Psi_i, i \in \Nbb`, of the Hilbert space.

	:Step 2:  For each :math:`\Psi_i`, choose a unit vector :math:`\Psi'_i` such that :math:`P[\Psi_i] = [\Psi'_i]`. Then :math:`\Psi'_i, i \in \Nbb`, also form an orthonormal basis by :math:`\eqref{eq_t_preserves_probability}`. We'd like to define :math:`U` by asking

		.. math::
			:nowrap:

			\begin{equation*}
				U\Psi_i = \Psi'_i
			\end{equation*}

		for all :math:`i`, and extend by (anti-)linearity. But this is not going to realize :math:`T` in general because we haven't fixed the extra degrees of freedom -- the phases of :math:`\Psi'_i`.

	:Step 3: We fix the phases of :math:`\Psi'_k, k \geq 1`, relative to :math:`\Psi'_0`, by asking

		.. math::
			:nowrap:

			\begin{equation*}
				T[\Psi_0 + \Psi_k] = [\Psi'_0 + \Psi'_k].
			\end{equation*}

		To see why this is possible, note first that :math:`T[\Psi_0 + \Psi_k] = [\alpha \Psi'_0 + \beta \Psi'_k]`, where :math:`\alpha, \beta` are phase factors, due to :math:`\eqref{eq_t_preserves_probability}` and the basis being orthonormal. Now :math:`[\alpha \Psi'_0 + \beta \Psi'_k] = [\Psi'_0 + (\beta/\alpha) \Psi'_k]` and we can absorb the phase :math:`\beta/\alpha` into the definition of :math:`\Psi'_k`. This is indeed the best one can do, because the last one degree of freedom, which is to multiply all :math:`\Psi'_i` by a phase, cannot be fixed.

	:Step 4: We have so far specified the value of :math:`U` on all of :math:`\Psi_i, i \geq 0`, and :math:`\Psi_0 + \Psi_k, k \geq 1`. Notice that all the coefficients of :math:`\Psi` are real. It is therefore instructive to ask what :math:`\Psi_0 + \ifrak \Psi_1` should be. By the same argument as in the previous step, we can write

		.. math::
			:nowrap:

			\begin{equation*}
				T[\Psi_0 + \ifrak \Psi_1] = [\Psi'_0 + c \Psi'_1]
			\end{equation*}

		where :math:`c` is a phase. Let's apply :math:`\eqref{eq_t_preserves_probability}` once again as follows

		.. math::
			:nowrap:

			\begin{align*}
				\sqrt{2} &= \left| \left( [\Psi_0 + \ifrak \Psi_1], [\Psi_0 + \Psi_1] \right) \right| \\
					&= \left| \left( T[\Psi_0 + \ifrak \Psi_1], T[\Psi_0 + \Psi_1] \right) \right| \\
					&= \left| \left( [\Psi'_0 + c \Psi'_1], [\Psi'_0 + \Psi'_1] \right) \right| \\
					&= |1 + c|
			\end{align*}

		It follows that :math:`c = \pm \ifrak`, which correspond to :math:`U` being (complex) linear or anti-linear, respectively.

	At this point, we can extend :math:`U` to either a linear or anti-linear map of the Hilbert space. But we'll not be bothered about any further formal argument, including showing that (anti-)linearity must be coupled with (anti-)unitarity, respectively.

.. note::
	The *adjoint* of a linear operator :math:`A` is another linear operator :math:`A^{\dagger}` such that

	.. math::
		:nowrap:

		\begin{equation*}
			(\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)
		\end{equation*}

	for all any two state vectors :math:`\Psi` and :math:`\Phi`. On the other hand, the adjoint of an anti-linear :math:`A` is another anti-linear :math:`A^{\dagger}` such that

	.. math::
		:nowrap:

		\begin{equation*}
			(\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)^{\ast}
		\end{equation*}

	A (anti-)unitary operator :math:`U` thus satisfies :math:`U^{\dagger} = U^{-1}`.

In general we're not interested in just one symmetry transformation, but rather a group -- whether continuous or discrete -- of symmetry transformations, or just symmetry for short. In particular, if :math:`T_1, T_2` are two symmetry transformations, then we'd like :math:`T_2 T_1` to also be a symmetry transformation. In light of the :math:`U`-realization of symmetry transformations discussed above, we can rephrase this condition as

.. math::
	:nowrap:

	\begin{equation}
		U(T_2 T_1) \Psi = \exp(\ifrak \theta(T_1, T_2, \Psi)) U(T_2) U(T_1) \Psi
		\label{eq_u_depends_on_psi}
	\end{equation}

where :math:`\ifrak = \sqrt{-1}`, and :math:`\theta(T_1, T_2, \Psi)` is an angle, which depends a priori on :math:`T_1, T_2`, and :math:`\Psi`.

It turns out, however, the angle :math:`\theta(T_1, T_2, \Psi)` cannot depend on the state because if we apply :math:`\eqref{eq_u_depends_on_psi}` to the sum of two linearly independent state vectors :math:`\Psi_A + \Psi_B`, then we'll find

.. math::
	:nowrap:

	\begin{equation*}
		\exp(\pm \ifrak \theta(\Psi_A)) \Psi_A + \exp(\pm \ifrak \theta(\Psi_B)) \Psi_B = \exp(\pm \ifrak \theta(\Psi_A + \Psi_B)) (\Psi_A + \Psi_B)
	\end{equation*}

where we have suppressed the dependency of :math:`\theta` on :math:`T`, and the signs correspond to the cases of :math:`U` being linear or anti-linear, respectively. In any case, it follows that

.. math::
	:nowrap:

	\begin{equation*}
		\exp(\pm \ifrak \theta(\Psi_A)) = \exp(\pm \ifrak \theta(\Psi_B)) = \exp(\pm \ifrak \theta(\Psi_A + \Psi_B))
	\end{equation*}

which says nothing but the independence of :math:`\theta` on :math:`\Psi`.

.. todo::
	While the argument here appears to be purely mathematical, Weinberg pointed out in [Wei95]_ (page 53) the potential inabilities to create a state like :math:`\Psi_A + \Psi_B`. More precisely, he mentioned the general believe that it's impossible to prepare a superposition of two states, one with integer total angular momentum and the other with half-integer total angular momentum, in which case there will be a "super-selection rule" between different classes of states. After all, one Hilbert space may just not be enough to describe all states. It'd be nice to elaborate a bit more on the super-selection rules.

We can now simplify :math:`\eqref{eq_u_depends_on_psi}` to the following

.. math::
	:nowrap:

	\begin{equation}
		U(T_2 T_1) = \exp(\ifrak \theta(T_1, T_2)) U(T_2) U(T_1)
		\label{eq_u_not_depend_on_psi}
	\end{equation}

which, in mathematical terms, says that :math:`U` furnishes a *projective representation* of :math:`T`, or a representation up to a phase. It becomes a genuine representation if the phase is constantly one.

.. _assump_genuine_repr:

.. admonition:: Assumption
	:class: Important

	We will assume that :math:`U` furnishes a genuine representation of :math:`T` unless otherwise stated, because it's simpler and will be suffice for most scenarios of interest.

.. _sec_continuous_symmetry:

Continuous symmetry
^^^^^^^^^^^^^^^^^^^

Besides a handful of important discrete symmetries such as the time, charge, and parity conjugations, most of the interesting symmetries come in a continuous family, mathematically known as *Lie groups*. Note that continuous symmetries are necessarily unitary (and linear) because they can be continuously deformed into the identity, which is obviously unitary.

In fact, it will be of great importance to just look at the symmetry up to the first order at the identity transformation, mathematically known as the *Lie algebra*. Let :math:`\theta` be an element in the Lie algebra such that :math:`T(\theta) = 1 + \theta` up to the first order. We can expand :math:`U(T(\theta))` in a power series as follows

.. math::
	:nowrap:

	\begin{equation}
		U(T(\theta)) = 1 + \ifrak \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots
		\label{eq_u_expansion}
	\end{equation}

where :math:`\theta^a` are the (real) components of :math:`\theta`, and :math:`u_a` are operators independent of :math:`\theta`, and as a convention, repeated indexes are summed up. Here we put a :math:`\ifrak` in front of the linear term so that the unitarity of :math:`U` implies that :math:`u_a` are Hermitian.

.. note::
	We've implicitly used a summation convention in writing :math:`\eqref{eq_u_expansion}` that the same upper and lower indexes are automatically summed up. For example

	.. math::
		:nowrap:

		\begin{equation*}
			\theta^a \theta^b u_{ab} \equiv \sum_{a, b} \theta^a \theta^b u_{ab}
		\end{equation*}

	This convention will be used throughout this note, unless otherwise specified.

	Another noteworthy point is how one writes matrix or tensor elements using indexes. The point is that the indexes must come in certain order. This wouldn't really cause a problem if all indexes are lower or upper. However, care must be taken when both lower and upper indexes appear. For example, an element written as :math:`M^a_b` would be ambiguous as it's unclear whether it refers to :math:`M_{ab}` or :math:`M_{ba}` assuming that one can somehow raise/lower the indexes. To avoid such ambiguity, one writes either :math:`{M^a}_b` or :math:`{M_b}^a`.

	This is a particularly convenient convention when dealing with matrix or tensor multiplications. For example, one can multiply two matrices as follows

	.. math::
		:nowrap:

		\begin{equation*}
			{M^a}_b {N^b}_c = {(MN)^a}_c
		\end{equation*}

	while :math:`{M^a}_b {N_c}^b`, though still summed up over :math:`b`, wouldn't correspond to a matrix multiplication.


Now let :math:`\eta` be another element of the Lie algebra, and expand both sides of :math:`U(T(\eta)) U(T(\theta)) = U(T(\eta) T(\theta))` as follows

.. math::
	:nowrap:

	\begin{align*}
	U(T(\eta)) U(T(\theta)) &= \left( 1 + \ifrak \eta^a u_a + \tfrac{1}{2} \eta^a \eta^b u_{ab} + \cdots \right) \left( 1 + \ifrak \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots \right) \\
		&= 1 + \ifrak (\eta^a + \theta^a) u_a \blue{- \eta^a \theta^b u_a u_b} + \cdots \\
	\\
	U(T(\eta) T(\theta)) &= U \left( 1 + \eta + \theta + f_{ab} \eta^a \theta^b + \cdots \right) \\
		&= 1 + \blue{\ifrak} \left( \eta^c + \theta^c + \blue{{f^c}_{ab} \eta^a \theta^b} + \cdots \right) \blue{u_c} + \blue{\tfrac{1}{2}} \left( \blue{\eta^a + \theta^a} + \cdots \right) \left( \blue{\eta^b + \theta^b} + \cdots \right) \blue{u_{ab}} + \cdots
	\end{align*}

where :math:`{f^c}_{ab}` are the coefficients of the expansion of :math:`T(f(\eta, \theta)) = T(\eta) T(\theta)`. Equating the coefficients of :math:`\eta^a \theta^b`, i.e., the terms colored in blue, we get

.. math::
	:nowrap:

	\begin{equation*}
		-u_a u_b = \ifrak {f^c}_{ab} u_c + u_{ab} \implies u_{ab} = -u_a u_b - \ifrak {f^c}_{ab} u_c.
	\end{equation*}

It implies that one can calculate the higher-order operator :math:`u_{ab}` from the lower-order ones, assuming of course that we know the structure of the symmetry (Lie) group/algebra. In fact, this bootstrapping procedure can be continued to all orders, but we'll not be bothered about the details.

Next, note that :math:`u_{ab} = u_{ba}` since they are just partial derivatives. It follows that

.. math::
	:nowrap:

	\begin{equation*}
		[u_a, u_b] \coloneqq u_a u_b - u_b u_a = \ifrak ({f^c}_{ba} - {f^c}_{ab}) u_c \eqqcolon \ifrak {C^c}_{ab} u_c
	\end{equation*}

where the bracket is known as the *Lie bracket* and :math:`{C^c}_{ab}` are known as the *structure constants*.

We conclude the general discussion about continuous symmetry by considering a special, but important, case when :math:`T` is additive in the sense that :math:`T(\eta) T(\theta) = T(\eta + \theta)`. Notable examples of such symmetry include translations and rotations about a fixed axis. In this case :math:`f` vanishes, and it follows from :math:`\eqref{eq_u_expansion}` that

.. math::
	:nowrap:

	\begin{equation}
		U(T(\theta)) = \lim_{N \to \infty} (U(T(\theta / N)))^N = \lim_{N \to \infty} (1 + \ifrak \theta^a u_a / N)^N = \exp(\ifrak \theta^a u_a)
		\label{eq_additive_symmetry}
	\end{equation}

.. _sec_lorentz_symmetry:

Lorentz symmetry
^^^^^^^^^^^^^^^^

A particularly prominent continuous symmetry in our physical world is the Lorentz symmetry postulated by Einstein's special relativity, which supersedes the Galilean symmetry, which is respected by the Newtonian mechanics. We shall start from the classical theory of Lorentz symmetry, and then quantize it following the procedure discussed in the previous section.

Classical Lorentz symmetry
++++++++++++++++++++++++++

Classical Lorentz symmetry is a symmetry that acts on the (flat) spacetime and preserves the so-called *proper time*

.. math::
	:nowrap:

	\begin{equation}
		d\tau^2 \coloneqq dx_0^2 - dx_1^2 - dx_2^2 - dx_3^2 \eqqcolon -\eta^{\mu \nu} dx_{\mu} dx_{\nu}
		\label{eq_proper_time}
	\end{equation}

where

1. :math:`x_0` is also known as the time, and sometimes denoted by :math:`t`, and
2. the speed of light is set to :math:`1`, and
3. :math:`\eta = \op{diag}(-1, 1, 1, 1)` and the indexes :math:`\mu, \nu` run from :math:`0` to :math:`3`.

.. note::
	1. We will follow the common convention in physics that Greek letters such as :math:`\mu, \nu, \dots` run from :math:`0` to :math:`3`, while Roman letters such as :math:`i, j, \dots` run from :math:`1` to :math:`3`.
	2. We often write :math:`x` for a spacetime point :math:`(x_0, x_1, x_2, x_3)`, and :math:`\xbf` for a spatial point :math:`(x_1, x_2, x_3)`.
	3. A :math:`4`-index, i.e., those named by Greek letters, of a matrix or a tensor can be raised or lowered by :math:`\eta`. For example, one can raise an index of a matrix :math:`M_{\mu \nu}` by :math:`\eta^{\rho \mu} M_{\mu \nu} = {M^{\rho}}_{\nu}` or :math:`\eta^{\rho \nu} M_{\mu \nu} = {M_{\mu}}^{\rho}`, such that the order of (regardless of upper or lower) indexes are kept.

.. dropdown:: Einstein's special theory of relativity
	:icon: unlock
	:animate: fade-in-slide-down

	Using the notations introduced above, we can rewrite :math:`\eqref{eq_proper_time}` as :math:`d\tau^2 = dt^2 - d\xbf^2`, so that it's obvious that if a particle travels at the speed of light in one inertial frame, i.e., :math:`|d\xbf / dt| = 1`, and equivalently :math:`d\tau = 0`, then it travels at the speed of light in any other inertial frame, in direct contradiction with Newtonian mechanics.

	Instead of working with the spacetime coordinates, it can sometimes be convenient to work with the "dual" energy-momentum coordinates, also known as the *four momentum*. The transition can be done by imagining a particle of mass :math:`m`, and defining :math:`p = (E, \pbf) \coloneqq m dx / d\tau`. It follows from :math:`\eqref{eq_proper_time}` that

	.. math::
		:nowrap:

		\begin{equation}
			1 = (dt / d\tau)^2 - (d\xbf / d\tau)^2 \implies m^2 = (m dt / d\tau)^2 - (m d\xbf / d\tau)^2 = E^2 - \pbf^2
			\label{eq_four_momentum_mass_identity}
		\end{equation}

	which looks just like :math:`\eqref{eq_proper_time}`, and indeed, the mass (in our convention) is invariant in all inertial frames.

	One can also recover Newtonian mechanics at the low-speed limit (i.e., :math:`|\vbf| \ll 1`) using :math:`d\tau / dt = \sqrt{1 - \vbf^2}` as follows

	.. math::
		:nowrap:

		\begin{align}
			\pbf &= m d\xbf / d\tau = \frac{m \vbf}{\sqrt{1 - \vbf^2}} = m \vbf + O(|\vbf|^3)  \label{eq_p_from_v} \\
			E &= m dt / d\tau = m + \tfrac{1}{2} m \vbf^2 + O(|\vbf|^4)  \nonumber
		\end{align}

More precisely, by a Lorentz transformation we mean an inhomogeneous linear transformation

.. math::
	:nowrap:

	\begin{equation*}
		L(\Lambda, a)x \coloneqq \Lambda x + a
	\end{equation*}

which consists of a homogeneous part :math:`\Lambda` and a translation by :math:`a`. The proper time is obviously preserved by any translation, and also by :math:`\Lambda` if

.. math::
	:nowrap:

	\begin{equation}
		\eta^{\mu \nu} dx_{\mu} dx_{\nu} = \eta^{\mu \nu} {\Lambda_{\mu}}^{\rho} {\Lambda_{\nu}}^{\kappa} dx_{\rho} dx_{\kappa} \
			\implies \eta^{\mu \nu} = \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu}
		\label{eq_homogeneous_lorentz_transformation}
	\end{equation}

for any :math:`\mu` and :math:`\nu`. Moreover the group law is given by

.. math::
	:nowrap:

	\begin{equation*}
		L(\Lambda', a') L(\Lambda, a) x = L(\Lambda', a')(\Lambda x + a) = \Lambda' \Lambda x + \Lambda' a + a' = L(\Lambda' \Lambda, \Lambda' a + a') x
	\end{equation*}

For later use, let's also calculate the inverse matrix of :math:`\Lambda` as follows

.. math::
	:nowrap:

	\begin{equation}
		\delta_{\sigma}^{\nu} = \eta_{\sigma \mu} \eta^{\mu \nu} \xlongequal{\eqref{eq_homogeneous_lorentz_transformation}} \eta_{\sigma \mu} \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu} \
			\implies {(\Lambda^{-1})_{\sigma}}^{\kappa} = \eta_{\sigma\mu} \eta^{\rho\kappa} {\Lambda_{\rho}}^{\mu} = {\Lambda^{\kappa}}_{\sigma}
		\label{eq_lambda_inverse}
	\end{equation}

Now we'll take a look at the topology of the group of homogeneous Lorentz transformations. Taking determinant on both sides of :math:`\eqref{eq_homogeneous_lorentz_transformation}`, we see that :math:`\op{det}(\Lambda) = \pm 1`. Moreover, setting :math:`\mu = \nu = 0`, we have

.. math::
	:nowrap:

	\begin{equation*}
		1 = \left( {\Lambda_0}^0 \right)^2 - \sum_{i=1}^3 \left( {\Lambda_i}^0 \right) \implies \left| {\Lambda_0}^0 \right| \geq 1
	\end{equation*}

It follows that the homogeneous Lorentz group has four components. In particular, the one with :math:`\op{det}(\Lambda) = 1` and :math:`{\Lambda_0}^0 \geq 1` is the most common used and is given a name: *proper orthochronous* Lorentz group. Nonetheless, one can map one component to another by composing with either a time reversal transformation

.. math::
	:nowrap:

	\begin{equation}
		\Tcal: (t, \xbf) \mapsto (-t, \xbf)
		\label{eq_time_inversion}
	\end{equation}

or a space reversal transformation

.. math::
	:nowrap:

	\begin{equation}
		\Pcal: (t, \xbf) \mapsto (t, -\xbf)
		\label{eq_space_inversion}
	\end{equation}

or both.

So far everything have been rather abstract, but in fact, the (homogeneous) Lorentz group can be understood quite intuitively. There are basically two building blocks: one is a rotation in the :math:`3`-space, which says that the space is homogeneous in all (spatial) directions, and the other is a so-called *boost*, which says that, as G. Galileo originally noted, one cannot tell if a system is at rest or is moving in a constant velocity without making a reference to outside of the system. To spell out the details, let's consider a rest frame with :math:`d\xbf = 0` and a moving frame with :math:`d\xbf' / dt' = \vbf`. Then the transformation :math:`dx' = \Lambda dx` can be simplified as

.. math::
	:nowrap:

	\begin{equation*}
		dt' = {\Lambda_0}^0 dt, \quad dx'_i = {\Lambda_i}^0 dt \implies {\Lambda_i}^0 = v_i {\Lambda_0}^0
	\end{equation*}

Then using :math:`\eqref{eq_homogeneous_lorentz_transformation}`, we get

.. math::
	:nowrap:

	\begin{equation}
		1 = -\eta^{\mu \nu} {\Lambda_{\mu}}^0 {\Lambda_{\nu}}^0
			= \left( {\Lambda_0}^0 \right)^2 - \left( {\Lambda_i}^0 \right)^2
			= \left( 1 - \vbf^2 \right) \left( {\Lambda_0}^0 \right)^2
			\implies {\Lambda_0}^0 = \frac{1}{\sqrt{1 - \vbf^2}} \eqqcolon \gamma
			\label{eq_def_gamma}
	\end{equation}

assuming :math:`\Lambda` is proper orthochronous. It follows that

.. math::
	:nowrap:

	\begin{equation}
		{\Lambda_i}^0 = {\Lambda^0}_i = \gamma v_i
		\label{eq_lambda_boost}
	\end{equation}

The other components :math:`{\Lambda_i}^j, 1 \leq i, j \leq 3`, are not uniquely determined because a composition with a (spatial) rotation about the direction of :math:`\vbf` has no effect on :math:`\vbf`. To make it easier, one can apply a rotation so that :math:`\vbf` aligns with the :math:`3`-axis. Then an obvious choice of :math:`\Lambda` is given by

.. math::
	:nowrap:

	\begin{alignat}{2}
		t'   &= {\Lambda_0}^{\mu} x_{\mu} &&= \gamma (t + v_3 x_3)
		\label{eq_lambda_in_3_axis} \\
		x'_1 &= {\Lambda_1}^{\mu} x_{\mu} &&= x_1
		\nonumber \\
		x'_2 &= {\Lambda_2}^{\mu} x_{\mu} &&= x_2
		\nonumber \\
		x'_3 &= {\Lambda_3}^{\mu} x_{\mu} &&= \gamma (x_3 + v_3 t)
		\nonumber
	\end{alignat}

.. dropdown:: Time dilation and length contraction
	:icon: unlock
	:animate: fade-in-slide-down

	A few consequences can be drawn from the boost transformation, most notably the effects of `time dilation <https://en.wikipedia.org/wiki/Time_dilation>`__ and `length contraction <https://en.wikipedia.org/wiki/Length_contraction>`__. The time dilation, i.e., a clock ticks slower in a moving frame than in a rest frame, is quite obvious from :math:`\eqref{eq_lambda_boost}` and the fact that :math:`\gamma > 1`. But the length contraction requires some elaboration.

	To be more concrete, let's consider a rode of some fixed length. To measure the length, the measurement must be done *simultaneously* at the two ends of the rod. This constraint causes not much trouble in a rest frame, but must be taken care of in a moving frame since being simultaneous is not a Lorentz invariant property. Let :math:`x = (t, \xbf)` and :math:`y = (t', \ybf)` be the two endpoints of the rod in the rest frame, so that the length is :math:`|\xbf - \ybf|` regardless of whether :math:`t` and :math:`t'` are the same or not. Under the Lorentz transformation defined in :math:`\eqref{eq_lambda_in_3_axis}`, they become

	.. math::
		:nowrap:

		\begin{align*}
			\Lambda x &= (\gamma (t + v_3 x_3), x_1, x_2, \gamma(x_3 + v_3 t))  \\
			\Lambda y &= (\gamma (t' + v_3 y_3), y_1, y_2, \gamma(y_3 + v_3 t'))
		\end{align*}

	respectively. Setting the equal-time condition :math:`(\Lambda x)_0 = (\Lambda y)_0` gives :math:`t' = t + v_3 (x_3 - y_3)`. Substituting it into :math:`(\Lambda x)_3` and :math:`(\Lambda y)_3` then gives

	.. math::
		:nowrap:

		\begin{equation*}
			|(\Lambda x)_3 - (\Lambda y)_3| = \gamma \left| x_3 - y_3 - v_3^2 (x_3 - y_3) \right| = \frac{|x_3 - y_3|}{\gamma} < |x_3 - y_3|
		\end{equation*}

	This calculation says that the length of rod is contracted in the direction of movement.

	It should be emphasized that such contraction of length can only be observed in a frame where the rod is moving. Imagine for example a scenario where you're given a square box with equal sides while standing still, then after some unconscious period of time, e.g., sleeping, you wake up with the same box in hand, and you'd like to know if you're now moving. If you happen to have heard of such contraction of length, you might try to measure the sides of the box again. If one of the sides suddenly becomes shorter, then you know not only that you're moving, but also the direction of movement! This is of course absurd because the box is still at rest relative to you.

Finally, one can apply a rotation to :math:`\eqref{eq_lambda_in_3_axis}` to get the general formula

.. math::
	:nowrap:

	\begin{equation}
		\Lambda_{ij} = \delta_{ij} + \frac{v_i v_j}{\vbf^2} (\gamma - 1)
		\label{eq_general_lambda_in_spacetime}
	\end{equation}

for :math:`1 \leq i, j \leq 3`, which, together with :math:`\eqref{eq_lambda_boost}` and :math:`{\Lambda_0}^i = {\Lambda_i}^0,` gives the general formula for :math:`\Lambda`.

.. note::
	Any Lorentz transformation can be written as the composition of a boost followed by a rotation.

.. _sec_quantum_lorentz_symmetry:

Quantum Lorentz symmetry
++++++++++++++++++++++++

We will quantize the Lorentz symmetry :math:`L(\Lambda, a)` by looking for unitarity representations :math:`U(\Lambda, a)`. As discussed in :ref:`sec_continuous_symmetry`, we proceed by looking for infinitesimal symmetries. First of all, let's expand :math:`\Lambda` as

.. math::
	:nowrap:

	\begin{equation}
		{\Lambda_{\mu}}^{\nu} = {\delta_{\mu}}^{\nu} + {\omega_{\mu}}^{\nu} + \cdots
		\label{eq_expansion_of_Lambda}
	\end{equation}

where :math:`\delta` is the Kronecker delta, and *not* a tensor. It follows from :math:`\eta^{\mu \nu} = \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu}` that

.. math::
	:nowrap:

	\begin{align}
		\eta^{\mu \nu} &= \eta^{\rho \kappa} ({\delta_{\rho}}^{\mu} + {\omega_{\rho}}^{\mu} + \cdots) ({\delta_{\kappa}}^{\nu} + {\omega_{\kappa}}^{\nu} + \cdots)
			\label{eq_lorentz_lie_algebra_is_antisymmetric} \\
			&= \eta^{\mu \nu} + \eta^{\mu \kappa} {\omega_{\kappa}}^{\nu} + \eta^{\nu \rho} {\omega_{\rho}}^{\mu} + \cdots
			\nonumber \\
			&= \eta^{\mu \nu} + \omega^{\mu \nu} + \omega^{\nu \mu} + \cdots
			\nonumber
	\end{align}

Comparing the first order terms shows that :math:`\omega^{\mu \nu} = -\omega^{\nu \mu}` is anti-symmetric. It is therefore more convenient to use :math:`\omega^{\mu \nu}`, rather than :math:`\omega_{\mu}^{\nu}`, as the infinitesimal parameters in the expansion of :math:`\Lambda`.

.. note::
	A count of free parameters shows that the inhomogeneous Lorentz symmetry has :math:`10` degrees of freedom, :math:`4` of which come from the translation, and the rest :math:`6` come from the rank-:math:`2` anti-symmetric tensor :math:`\omega`.

We first postulate that :math:`U(1, 0) = 1` is the identity operator because the Lorentz transformation itself is the identity. Then we can write the power series expansion up to first order as follows

.. math::
	:nowrap:

	\begin{equation}
		U(1 + \omega, \epsilon) = 1 - \ifrak \epsilon^{\mu} P_{\mu} + \frac{\ifrak}{2} \omega^{\mu \nu} J_{\mu \nu} + \cdots
		\label{eq_u_lorentz_expansion}
	\end{equation}

Here we have inserted :math:`\ifrak` as usual so that the unitarity of :math:`U` implies that both :math:`P_{\mu}` and :math:`J_{\mu \nu}` are
Hermitian. Moreover, since :math:`\omega^{\mu \nu}` is anti-symmetric, we can assume the same holds for :math:`J_{\mu \nu}`.

.. note::
	Since we are expanding :math:`U(1 + \epsilon)` which is complex linear, the operators :math:`P` and :math:`J` are also complex linear. Hence we can freely move :math:`\ifrak` around these operators in calculations that follow. However, this will become an issue when we later consider other operators such as the space and time inversions, which can potentially be either complex linear or anti-linear. In the later case, a sign needs to be added when commuting with the multiplication by :math:`\ifrak`.

Let's evaluate how the expansion transformations under conjugation

.. math::
	:nowrap:

	\begin{align*}
		U(\Lambda, a) U(1 + \omega, \epsilon) U^{-1}(\Lambda, a) \
			&= U(\Lambda, a) U(1 + \omega, \epsilon) U(\Lambda^{-1}, -\Lambda^{-1} a) \\
			&= U(\Lambda, a) U((1 + \omega) \Lambda^{-1}, \epsilon - (1 + \omega) \Lambda^{-1} a) \\
			&= U(1 + \Lambda \omega \Lambda^{-1}, \Lambda \epsilon - \Lambda \omega \Lambda^{-1} a) \\
			&= 1 - \ifrak ({\Lambda^{\rho}}_{\mu} \epsilon^{\mu} - (\Lambda \omega \Lambda^{-1})^{\rho \kappa} a_{\kappa}) P_{\rho} \
				+ \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} J_{\rho \kappa} + \cdots \\
			&= 1 -\ifrak \epsilon^{\mu} {\Lambda^{\rho}}_{\mu} P_{\rho} + \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots \\
			&= 1 -\ifrak \epsilon^{\mu} {\Lambda^{\rho}}_{\mu} P_{\rho} + \tfrac{\ifrak}{2} {\Lambda^{\rho}}_{\mu} \omega^{\mu \nu} {\Lambda^{\kappa}}_{\nu} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots
	\end{align*}


where we have used :math:`\eqref{eq_lambda_inverse}` for :math:`\Lambda^{-1}`. Substituting :math:`U(1 + \omega, \epsilon)` with the expansion :math:`\eqref{eq_u_lorentz_expansion}` and equating the coefficients of :math:`\epsilon^{\mu}` and :math:`\omega_{\mu \nu}`, we have

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda, a) P_{\mu} U^{-1}(\Lambda, a) &= {\Lambda^{\rho}}_{\mu} P_{\rho}  \label{eq_p_conjugated_by_u} \\
		U(\Lambda, a) J_{\mu \nu} U^{-1}(\Lambda, a) &= {\Lambda^{\rho}}_{\mu} {\Lambda^{\kappa}}_{\nu} (J_{\rho \kappa} + a_{\kappa} P_{\rho} - a_{\rho} P_{\kappa})  \label{eq_j_conjugated_by_u}
	\end{align}

where in the second equation, we have also made the right-hand-side anti-symmetric with respect to :math:`\mu` and :math:`\nu`. It's now clear that :math:`P` transforms like a vector and is translation invariant, while :math:`J` transforms like a :math:`2`-tensor only for homogeneous Lorentz transformations and is not translation invariant in general. These are of course as expected since both :math:`P` and :math:`J` are quantization of rather familiar objects, which we now spell out.

We start with :math:`P` by writing :math:`H \coloneqq P_0` and :math:`\Pbf \coloneqq (P_1, P_2, P_3)`. Then :math:`H` is the energy operator, also know as the *Hamiltonian*, and :math:`\Pbf` is the momentum :math:`3`-vector. Similarly, let's write :math:`\Kbf \coloneqq (J_{01}, J_{02}, J_{03})` and :math:`\Jbf = (J_{23}, J_{31}, J_{12})`, as the *boost* :math:`3`-vector and the *angular momentum* :math:`3`-vector, respectively.

Now that we have named all the players (i.e., :math:`H, \Pbf, \Jbf, \Kbf`) in the game, it remains to find out their mutual commutation relations since they should form a Lie algebra of the (infinitesimal) Lorentz symmetry. This can be done by applying :math:`\eqref{eq_p_conjugated_by_u}` and :math:`\eqref{eq_j_conjugated_by_u}` to :math:`U(\Lambda, a)` that is itself infinitesimal. More precisely, keeping up to first order terms, we have :math:`{\Lambda^{\rho}}_{\mu} = {\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu}` and :math:`a_{\mu} = \epsilon_{\mu}`. It follows that :math:`\eqref{eq_p_conjugated_by_u}`, up to first order, can be written as follows

.. math::
	:nowrap:

	\begin{align*}
		\left( {\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu} \right) P_{\rho} &= \left( 1 - \ifrak \epsilon^{\nu} P_{\nu} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) P_{\mu} \left( 1 + \ifrak \epsilon^{\nu} P_{\nu} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
			&= P_{\mu} - \ifrak \epsilon^{\nu} [P_{\mu}, P_{\nu}] - \tfrac{\ifrak}{2} \omega^{\rho \kappa} [P_{\mu}, J_{\rho \kappa}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` and :math:`\omega` gives the following

.. math::
	:nowrap:

	\begin{align}
		[P_{\mu}, P_{\nu}] &= 0  \label{eq_bracket_p4_p4} \\
		[P_{\mu}, J_{\rho \kappa}] &= \ifrak (\eta_{\kappa \mu} P_{\rho} - \eta_{\rho \mu} P_{\kappa})  \label{eq_bracket_p4_j4}
	\end{align}

where for the second identity, we've also used the fact that :math:`J_{\rho \kappa} = -J_{\kappa \rho}`.

Similarly, expanding :math:`\eqref{eq_j_conjugated_by_u}` up to first order, we have

.. math::
	:nowrap:

	\begin{align*}
		J_{\mu \nu} + \epsilon_{\nu} P_{\mu} - \epsilon_{\mu} P_{\nu} + {\omega^{\rho}}_{\mu} J_{\rho \nu} + {\omega^{\kappa}}_{\nu} J_{\mu \kappa} \
			&= ({\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu}) ({\delta^{\kappa}}_{\nu} + {\omega^{\kappa}}_{\nu}) (J_{\rho \kappa} + \epsilon_{\kappa} P_{\rho} - \epsilon_{\rho} P_{\kappa}) \\
			&= \left( 1 - \ifrak \epsilon^{\rho} P_{\rho} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) J_{\mu \nu} \left( 1 + \ifrak \epsilon^{\rho} P_{\rho} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
			&= J_{\mu \nu} - \ifrak \epsilon^{\rho} [P_{\rho}, J_{\mu \nu}] + \tfrac{\ifrak}{2} \omega^{\rho \kappa} [J_{\rho \kappa}, J_{\mu \nu}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` reproduces :math:`\eqref{eq_bracket_p4_j4}`, but equating the coefficients of :math:`\omega` gives the following additional

.. math::
	:nowrap:

	\begin{equation}
		[J_{\rho \kappa}, J_{\mu \nu}] = \ifrak (\eta_{\rho \nu} J_{\mu \kappa} - \eta_{\rho \mu} J_{\nu \kappa} - \eta_{\kappa \mu} J_{\rho \nu} - \eta_{\kappa \nu} J_{\rho \mu})
		\label{eq_bracket_j4_j4}
	\end{equation}

Now that we have all the commutator relations, let's reorganize :math:`\eqref{eq_bracket_p4_p4}, \eqref{eq_bracket_p4_j4}, \eqref{eq_bracket_j4_j4}` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		\text{let } \mu = 0, \nu = i \text{ in \eqref{eq_bracket_p4_p4}} &\implies [H, P_i] &&= 0
		\label{eq_hp_commute} \\
		\text{let } \mu = 0, \rho = j, \kappa = k \text{ in \eqref{eq_bracket_p4_j4}} &\implies [H, J_i] &&= 0
		\label{eq_hj_commute} \\
		\text{let } \mu = 0, \rho = 0, \kappa = i \text{ in \eqref{eq_bracket_p4_j4}} &\implies [H, K_i] &&= \ifrak P_i
		\label{eq_hkp_commutation} \\
		\text{let } \mu = i, \nu = j \text{ in \eqref{eq_bracket_p4_p4}} &\implies [P_i, P_j] &&= 0
		\label{eq_pp_commute} \\
		\text{let } \mu = i, \rho = k, \kappa = i \text{ in \eqref{eq_bracket_p4_j4} and permutation (anti-)symmetry} &\implies [P_i, J_j] &&= \ifrak \epsilon_{ijk} P_k  \label{eq_pjp_commutation} \\
		\text{let } \mu = i, \rho = 0, \kappa = j \text{ in \eqref{eq_bracket_p4_j4}} &\implies [P_i, K_j] &&= \ifrak \delta_{ij} H  \label{eq_pkh_commutation} \\
		\text{let } \rho = j, \kappa = \mu = k, \nu = i \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} &\implies [J_i, J_j] &&= \ifrak \epsilon_{ijk} J_k  \label{eq_jjj_commutation} \\
		\text{let } \rho = \nu = j, \kappa = k, \mu = 0 \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} &\implies [J_i, K_j] &&= \ifrak \epsilon_{ijk} K_k  \label{eq_jkk_commutation} \\
		\text{let } \rho = \mu = 0, \kappa = i, \nu = j \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} &\implies [K_i, K_j] &&= -\ifrak \epsilon_{ijk} J_k
		\label{eq_kkj_commutation}
 	\end{alignat}

where :math:`\epsilon_{ijk}` is totally anti-symmetric with respect to permutations of indexes and satisfies :math:`\epsilon_{123} = 1`. [#tedious_calc_of_commutations]_

.. note::
	The Lie algebra generated by :math:`H, \Pbf, \Jbf, \Kbf` with commutation relations :math:`\eqref{eq_hp_commute}` -- :math:`\eqref{eq_kkj_commutation}` is known as the `Poincaré algebra <https://en.wikipedia.org/wiki/Poincar%C3%A9_group>`__.

Since the time evolution of a physical system is dictated by the Hamiltonian :math:`H`, quantities (i.e., observables) that commute with :math:`H` are conserved. In particular :math:`\eqref{eq_hp_commute}` and :math:`\eqref{eq_hj_commute}` imply that both momentum and angular momentum are conserved. Boosts, on the other hand, are *not* conserved, and therefore cannot be used to label (stable) physical states. Moreover :math:`\eqref{eq_pp_commute}` implies that translations commute with each other (as expected), which is *not* the case for the angular momenta according to :math:`\eqref{eq_jjj_commutation}`. Indeed, they furnish an infinitesimal representation of the :math:`3`-rotation group :math:`SO(3)`.

.. _sec_one_particle_states:

One-Particle States
-------------------

One neat application of our knowledge about Lorentz symmetry is to classify (free) one-particle states according to their transformation laws under (inhomogeneous) Lorentz transformations. Throughout this section, the Lorentz transformations will be assumed to be proper orthochronous, i.e., :math:`\op{det}(\Lambda) = 1` and :math:`{\Lambda_0}^0 \geq 1`.

In order to do so, we need some labels to identify states, which are typically conserved quantities. According to the commutation relations between :math:`H, \Pbf` and :math:`\Jbf` obtained in the previous section, we see that :math:`p = (H, \Pbf)` consists of mutually commutative conserved components, but not :math:`\Jbf`. Hence we can write our one-particle states as :math:`\Psi_{p, \sigma}` such that

.. math::
	:nowrap:

	\begin{equation*}
		P_{\mu} \Psi_{p, \sigma} = p_{\mu}
	\end{equation*}

where :math:`\sigma` are additional labels such as spin components that we will later specify.

Reduction to the little group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's first consider translations :math:`U(1, a)`. Since translations form an abelian group, it follows from :math:`\eqref{eq_additive_symmetry}` that

.. math::
	:nowrap:

	\begin{equation}
		U(1, a) \Psi_{p, \sigma} = \exp(-\ifrak a^{\mu} P_{\mu}) \Psi_{p, \sigma} = \exp(-\ifrak a^{\mu} p_{\mu}) \Psi_{p, \sigma}
		\label{eq_translation_formula_for_particle_state}
	\end{equation}

where the minus sign comes from our choice of expansion :math:`\eqref{eq_u_lorentz_expansion}`. Hence it remains to consider the action of homogeneous Lorentz transformations. For the convenience of notation, let's write :math:`U(\Lambda) \coloneqq U(\Lambda, 0)`. We would first like to know how :math:`U(\Lambda)` affects the :math:`4`-momentum. It follows from the following calculation

.. math::
	:nowrap:

	\begin{equation*}
		P_{\mu} U(\Lambda) \Psi_{p, \sigma}
			= U(\Lambda) (U^{-1} (\Lambda) P_{\mu} U(\Lambda)) \Psi_{p, \sigma}
			\xlongequal{\eqref{eq_p_conjugated_by_u}} U(\Lambda) {\Lambda_{\mu}}^{\nu} P_{\nu} \Psi_{p, \sigma}
			= \left( {\Lambda_{\mu}}^{\nu} p_{\nu} \right) U(\Lambda) \Psi_{p, \sigma}
	\end{equation*}

that :math:`U(\Lambda) \Psi_{p, \sigma}` has :math:`4`-momentum :math:`\Lambda p`. Therefore we can write

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \Psi_{p, \sigma} = C_{\sigma \sigma'} (\Lambda, p) \Psi_{\Lambda p, \sigma'}
		\label{eq_lorentz_acts_on_p_and_sigma}
	\end{equation}

where :math:`C_{\sigma \sigma'}` furnishes a representation of :math:`\Lambda` and :math:`p` under straightforward transformation rules, and an implicit summation over :math:`\sigma'` is assumed although it's not a :math:`4`-index.

Next we'd like to remove the dependency of :math:`C_{\sigma \sigma'}` on :math:`p` since, after all, it is :math:`\Lambda` that carries the symmetry. We can achieve this by noticing that :math:`U(\Lambda)` acts on the :math:`\Lambda`-orbits of :math:`p` transitively. The :math:`\Lambda`-orbits of :math:`p`, in turn, are uniquely determined by the value of :math:`p^2`, and in the case of :math:`p^2 \leq 0`, also by the sign of :math:`p_0`. In light of :math:`\eqref{eq_four_momentum_mass_identity}`, we can pick a convenient representative :math:`k` for each case as follows

+---------------------------------+-----------------------+----------+
| Case                            | Standard :math:`k`    | Physical |
+=================================+=======================+==========+
| :math:`p^2 = -M^2 < 0,~p_0 > 0` | :math:`(M, 0, 0, 0)`  | Yes      |
+---------------------------------+-----------------------+----------+
| :math:`p^2 = -M^2 < 0,~p_0 < 0` | :math:`(-M, 0, 0, 0)` | No       |
+---------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 > 0`        | :math:`(1, 0, 0, 1)`  | Yes      |
+---------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 = 0`        | :math:`(0, 0, 0, 0)`  | Yes      |
+---------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 < 0`        | :math:`(-1, 0, 0, 1)` | No       |
+---------------------------------+-----------------------+----------+
| :math:`p^2 = N^2 > 0`           | :math:`(0, N, 0, 0)`  | No       |
+---------------------------------+-----------------------+----------+

It turns out that only three of these cases are realized physically, and they correspond to the cases of a massive particle of mass :math:`M`, a massless particle and the vacuum, respectively. Since there is not much to say about the vacuum state, there are only two cases that we need to investigate.

With the choices of the standard :math:`k` in hand, we need to make one more set of choices. Namely, we will choose for each :math:`p` a standard Lorentz transformation :math:`L(p)` such that :math:`L(p) k = p`. Such :math:`L(p)` for a massive particle has been chosen in :math:`\eqref{eq_general_lambda_in_spacetime}`, albeit in spacetime coordinates, and we'll also handle the case of massless particles later. Once these choices have been made, we can *define*

.. math::
	:nowrap:

	\begin{equation}
		\Psi_{p, \sigma} \coloneqq N(p) U(L(p)) \Psi_{k, \sigma}
		\label{eq_def_of_one_particle_psi}
	\end{equation}

where :math:`N(p)` is a normalization factor to be determined later. In this way, we've also determined how :math:`\sigma` depends on :math:`p`. Applying :math:`\eqref{eq_lorentz_acts_on_p_and_sigma}` to :math:`\eqref{eq_def_of_one_particle_psi}` we can refactor the terms as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \Psi_{p, \sigma} = N(p) U(\Lambda) U(L(p)) \Psi_{k, \sigma} = \
		N(p) U(L(\Lambda p)) U(L(\Lambda p)^{-1} \Lambda L(p)) \Psi_{k, \sigma}
		\label{eq_def_of_one_particle_psi_refactored}
	\end{equation}

so that :math:`L(\Lambda p)^{-1} \Lambda L(p)` maps :math:`k` to itself, and hence :math:`U(L(\Lambda p)^{-1} \Lambda L(p))` acts solely on :math:`\sigma`.

At this point, we have reduced the problem to the classification of representations of the so-called *little group* defined as the subgroup of (proper orthochronous) Lorentz transformations :math:`W` that fixes :math:`k`, i.e., :math:`{W_{\mu}}^{\nu} k_{\nu} = k_{\mu}`. Element in the little group is known as `Wigner rotation <https://en.wikipedia.org/wiki/Wigner_rotation>`__ (and hence :math:`W`). More precisely, the task now is to find (unitary) representations :math:`D(W)` such that

.. math::
	:nowrap:

	\begin{equation*}
		\sum_{\sigma'} D_{\sigma \sigma'}(W_1) D_{\sigma' \sigma''}(W_2) \Psi_{k, \sigma''} = D_{\sigma \sigma''}(W_1 W_2) \Psi_{k, \sigma''}
	\end{equation*}

Once this is done, we can define

.. math::
	:nowrap:

	\begin{align}
		U(W) \Psi_{k, \sigma} &\coloneqq \sum_{\sigma'} D_{\sigma' \sigma}(W) \Psi_{k, \sigma'}  \label{eq_d_repr_of_little_group} \\
		W(\Lambda, p) &\coloneqq L(\Lambda p)^{-1} \Lambda L(p)  \label{eq_w_from_l}
	\end{align}

One can verify that :math:`\eqref{eq_d_repr_of_little_group}` indeed respects the group law as follows

.. math::
	:nowrap:

	\begin{equation*}
		U(W_2) U(W_1) \Psi_{k, \sigma} = U(W_2) \sum_{\sigma'} D_{\sigma' \sigma}(W_1) \Psi_{k, \sigma'} \
			= \sum_{\sigma' \sigma''} D_{\sigma' \sigma}(W_1) D_{\sigma'' \sigma'}(W_2) \Psi_{k, \sigma''} \
			= \sum_{\sigma''} D_{\sigma'' \sigma}(W_2 W_1) \Psi_{k, \sigma''}
	\end{equation*}

Now we can rewrite :math:`\eqref{eq_def_of_one_particle_psi_refactored}` as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda) \Psi_{p, \sigma} &\xlongequal{\eqref{eq_w_from_l}} N(p) U(L(\Lambda p)) U(W(\Lambda, p)) \Psi_{k, \sigma} \label{eq_little_group_acts_on_p_and_sigma} \\
			&\xlongequal{\eqref{eq_d_repr_of_little_group}} N(p) \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) U(L(\Lambda p)) \Psi_{k, \sigma'}  \nonumber \\
			&\xlongequal{\eqref{eq_def_of_one_particle_psi}} \frac{N(p)}{N(\Lambda p)} \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'} \nonumber
	\end{align}

which gives the sought-after coefficients :math:`C_{\sigma \sigma'}` in :math:`\eqref{eq_lorentz_acts_on_p_and_sigma}`.

It remains now, as far as the general discussion is concerned, to settle the normalization factor :math:`N(p)`. Indeed, it'd not have been needed at all if we'd like :math:`\Psi_{p, \sigma}` be to orthonormal in the sense that

.. math::
	:nowrap:

	\begin{equation}
		(\Psi_{p', \sigma'}, \Psi_{p, \sigma}) = \delta_{\sigma' \sigma} \delta(p' - p)
		\label{eq_psi_p4_sigma_orthonormal}
	\end{equation}

where the first delta is the Kronecker delta (for discrete indexes) and the second is the Dirac delta (for continuous indexes), since they are eigenvectors of the (Hermitian) operator :math:`P`. All we need is :math:`D_{\sigma \sigma'}` being unitary as is obvious from :math:`\eqref{eq_little_group_acts_on_p_and_sigma}`.

However, the Dirac delta in :math:`\eqref{eq_psi_p4_sigma_orthonormal}` is tricky to use since :math:`p` is constrained to the so-called *mass shell*, i.e., :math:`p_0 > 0` together with :math:`p^2 = -M^2` in the massive case and :math:`p^2 = 0` in the massless case, respectively. Hence the actual normalization we'd like to impose on the one-particle states is, instead of :math:`\eqref{eq_psi_p4_sigma_orthonormal}`, the following

.. math::
	:nowrap:

	\begin{equation}
		(\Psi_{p', \sigma'}, \Psi_{p, \sigma}) = \delta_{\sigma' \sigma} \delta(\pbf' - \pbf)
		\label{eq_psi_p3_sigma_orthonormal}
	\end{equation}

In fact, the problem eventually boils down to how to define the :math:`3`-momentum space Dirac delta in a Lorentz-invariant manner.

Since :math:`\Psi_{p, \sigma}` can be derived from :math:`\Psi_{k, \sigma}` by :math:`\eqref{eq_def_of_one_particle_psi}`, we can first ask :math:`\Psi_{k, \sigma}` to be orthonormal in the sense of :math:`\eqref{eq_psi_p3_sigma_orthonormal}`, where the Dirac delta plays no role, and then figure out how integration works on the mass shell (because Dirac delta is defined by integrals against test functions). As far as the mass shell integration is concerned, we can temporarily unify the massive and massless cases by allowing :math:`M \geq 0`. Consider a general mass shell integral of an arbitrary test function :math:`f(p)`

.. math::
	:nowrap:

	\begin{align*}
		\int d^4 p ~\delta(p^2 + M^2) \theta(p_0) f(p) &= \int d^3\pbf dp_0 ~\delta(p_0^2 - \pbf^2 - M^2) \theta(p_0) f(p_0, \pbf) \\
			&= \int d^3\pbf ~\frac{f\left( \sqrt{\pbf^2 + M^2}, \pbf \right)}{2 \sqrt{\pbf^2 + M^2}}
	\end{align*}

where :math:`\theta(p_0)` is the step function defined to be :math:`0` if :math:`p_0 \leq 0` and :math:`1` if :math:`p_0 > 1`. It follows that the Lorentz-invariant volume element in the :math:`3`-momentum space is

.. math::
	:nowrap:

	\begin{equation}
		\frac{d^3\pbf}{\sqrt{\pbf^2 + M^2}}  \label{eq_lorentz_invariant_3_momentum_volume_element}
	\end{equation}

We can use it to find the Lorentz-invariant Dirac delta (marked in blue) as follows

.. math::
	:nowrap:

	\begin{align*}
		f(\pbf') &\eqqcolon \int d^3\pbf ~\delta(\pbf' - \pbf) f(\pbf) \\
			&= \int \frac{d^3\pbf}{\sqrt{\pbf^2 + M^2}} \blue{p_0 \delta(\pbf' - \pbf)} f(\pbf)
	\end{align*}

It follows from Lorentz invariance that :math:`p_0 \delta(\pbf' - \pbf) = k_0 \delta(\kbf' - \kbf)`. Hence we can finally establish :math:`\eqref{eq_psi_p3_sigma_orthonormal}` as follows

.. math::
	:nowrap:

	\begin{align*}
		(\Psi_{p', \sigma'}, \Psi_{p, \sigma}) &= N(p) N(p')^{\ast} (U(L(p')) \Psi_{k', \sigma'}, U(L(p)) \Psi_{k, \sigma}) \\
			&= |N(p)|^2 \delta_{\sigma' \sigma} \delta(\kbf' - \kbf) \\
			&= \delta_{\sigma' \sigma} \delta(\pbf' - \pbf)
	\end{align*}

if we define :math:`N(p) = \sqrt{k_0 / p_0}`.

Putting everything together, we've obtained the following grand formula for the Lorentz transformation law

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \Psi_{p, \sigma} = \sqrt{\frac{(\Lambda p)_0}{p_0}} \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'}
		\label{eq_lorentz_transformation_formula_for_particle_state}
	\end{equation}

where :math:`D_{\sigma' \sigma}` is a unitary representation of the little group, and :math:`W(\Lambda, p)` is defined by :math:`\eqref{eq_w_from_l}`.

Massive particle states
^^^^^^^^^^^^^^^^^^^^^^^

Recall the standard :math:`4`-momentum :math:`k = (M, 0, 0, 0)` in this case. Obviously the little group here is nothing but the :math:`3`-rotation group :math:`SO(3)`. We can work out :math:`D_{\sigma \sigma'}(\Rcal)` by a rotation :math:`\Rcal \in SO(3)` up to first order as follows.

First write :math:`\Rcal^{ij} = \delta^{ij} + \Theta^{ij}` such that :math:`\Theta` is anti-symmetric. Then expand :math:`D_{\sigma \sigma'} (\Rcal)` similar to :math:`\eqref{eq_u_lorentz_expansion}` up to first order as follows

.. math::
	:nowrap:

	\begin{equation*}
		D_{\sigma \sigma'} (\Rcal) = \delta_{\sigma \sigma'} + \tfrac{\ifrak}{2} \Theta^{ij} (J_{ij})_{\sigma \sigma'}
	\end{equation*}

where :math:`J_{ij}` is a collection of Hermitian operators that satisfy :math:`J_{ij} = -J_{ji}` and the commutation relations :math:`\eqref{eq_jjj_commutation}`. It turns out that there exists an infinite number of such unitary representations indexed by nonnegative half-integers :math:`\jfrak = 0, \tfrac{1}{2}, 1, \tfrac{3}{2}, \cdots`, each of which has dimension :math:`2\jfrak + 1`. Choosing the :math:`3`-axis as the preferred axis of (definite) spin, we can summarize the result as follows

.. math::
	:nowrap:

	\begin{align}
		D^{(\jfrak)}_{\sigma \sigma'} (\Rcal) &= \delta_{\sigma \sigma'} + \tfrac{\ifrak}{2} \Theta^{ij} \left( J^{(\jfrak)}_{ij} \right)_{\sigma \sigma'}  \label{eq_representation_rotation_first_order} \\
		\left( J^{(\jfrak)}_{23} \pm \ifrak J^{(\jfrak)}_{31} \right)_{\sigma \sigma'} = \left( J^{(\jfrak)}_1 \pm \ifrak J^{(\jfrak)}_2 \right)_{\sigma \sigma'} &= \delta_{\sigma \pm 1, \sigma'} \sqrt{(\jfrak \mp \sigma)(\jfrak \pm \sigma + 1)}  \label{eq_j1_j2_matrix} \\
		\left( J^{(\jfrak)}_{12} \right)_{\sigma \sigma'} = \left( J^{(\jfrak)}_3 \right)_{\sigma \sigma'} &= \sigma \delta_{\sigma \sigma'}  \label{eq_j3_matrix}
	\end{align}

where :math:`\sigma, \sigma'` run through the values :math:`-\jfrak, -\jfrak + 1, \cdots, \jfrak - 1, \jfrak`.

.. _dropdown_repr_of_angular_momenta:

.. dropdown:: Representations of angular momenta
	:icon: unlock
	:animate: fade-in-slide-down

	Recall from :math:`\eqref{eq_jjj_commutation}` that the (quantum) angular momenta vector :math:`\Jbf` satisfy the commutation relations :math:`[J_i, J_j] = \ifrak \epsilon_{ijk} J_k`. Hence they cannot be simultaneously diagonalized. It's then a convention to use the angular momentum along the :math:`3`-axis to label the spin. The following two identities are straightforward but important

	.. math::
		:nowrap:

		\begin{align*}
			[\Jbf^2, J_i] &= 0, ~\forall i = 1, 2, 3 \\
			[J_3, J_1 \pm \ifrak J_2] &= \pm (J_1 \pm \ifrak J_2)
		\end{align*}

	where :math:`\Jbf^2 = J_1^2 + J_2^2 + J_3^2` as usual.

	Now if :math:`\Psi_{\sigma}` is an eigenstate of :math:`J_3` with eigenvalue :math:`\sigma`, then

	.. math::
		:nowrap:

		\begin{equation}
			J_3 (J_1 \pm \ifrak J_2) \Psi_{\sigma} = [J_3, J_1 \pm \ifrak J_2] \Psi_{\sigma} + (J_1 \pm \ifrak J_2) J_3 \Psi_{\sigma} = (\sigma \pm 1) \Psi_{\sigma}
			\label{eq_j1_j2_raises_or_lowers_state}
		\end{equation}

	In other words, applying :math:`J_1 \pm \ifrak J_2` to any eigenstate of :math:`J_3` raises or lowers the eigenvalue by one, and henceforth they are called *raising* and *lowering* operators, respectively. Moreover, since :math:`\Jbf^2` commutes with :math:`J_3`, we may assume that :math:`\Psi_{\sigma}` is also an eigenstate of :math:`\Jbf^2`, and since :math:`\Jbf^2` also commutes with :math:`J_1 \pm \ifrak J_2`, the whole series of :math:`J_3`-eigenstates obtained by applying the raising/lowering operators have the same :math:`\Jbf^2`-eigenvalue.

	We'll from now on focus on eigenstates with a fixed :math:`\Jbf^2`-eigenvalue. Moreover we'd like the eigenvalues of :math:`J_3` to be bounded, so both the raising and the lowering operations must stop after finite steps. Let :math:`\Psi_{\jfrak}` be the :math:`J_3`-eigenstate with the highest eigenvalue (if there are more than one, the representation is reducible). By repeatedly applying the lowering operator to :math:`\Psi_{\jfrak}`, we'll eventually reach the eigenstate :math:`\Psi_{\jfrak'}` with the lowest eigenvalue. Since the lowering operator decreases the eigenvalue by one, we know that :math:`\jfrak - \jfrak'` must be an integer.

	Consider the following two operators

	.. math::
		:nowrap:

		\begin{align}
			(J_1 - \ifrak J_2) (J_1 + \ifrak J_2) &= J_1^2 + J_2^2 + \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 - J_3  \label{eq_j1_j2_mixed_product_plus} \\
			(J_1 + \ifrak J_2) (J_1 - \ifrak J_2) &= J_1^2 + J_2^2 - \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 + J_3  \label{eq_j1_j2_mixed_product_minus}
		\end{align}

	Note that the first operator annihilates :math:`\Psi_{\jfrak}` and the second operator annihilates :math:`\Psi_{\jfrak'}` by assumption, which, together with the fact that :math:`\Jbf^2 \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'}`, implies

	.. math::
		:nowrap:

		\begin{equation}
			\Jbf^2 \Psi_{\jfrak} = (\jfrak^2 + \jfrak) \Psi_{\jfrak} = ((\jfrak')^2 - \jfrak') \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'} \implies \jfrak (\jfrak + 1) = \jfrak' (\jfrak' - 1)
			\label{eq_angular_momentum_squared_eigenvalue}
		\end{equation}

	The equation has two potential solutions: either :math:`\jfrak' = \jfrak + 1` or :math:`\jfrak = -\jfrak'`. The first option violates the maximality of :math:`\jfrak`, and so we must accept the second option. Since we also know :math:`\jfrak - \jfrak'` must be integral, we conclude that :math:`\jfrak` is itself a half-integer.

	As a piece of notation, we'll from now on write :math:`\Psi^{\jfrak}_{\sigma}` for the eigenstate of both :math:`\Jbf^2` and :math:`J_3` such that

	.. math::
		:nowrap:

		\begin{align}
			\Jbf^2 \Psi^{\jfrak}_{\sigma} &= \jfrak (\jfrak + 1) \Psi^{\jfrak}_{\sigma}  \label{eq_3j_square_eigenstate} \\
			J_3 \Psi^{\jfrak}_{\sigma} &= \sigma \Psi^{\jfrak}_{\sigma}  \nonumber
		\end{align}

	We'll also write :math:`J_i^{(\jfrak)} \coloneqq J_i` to explicitly indicate the dependency on :math:`\jfrak`.

	It remains to settle the constant term on the right-hand-side of :math:`\eqref{eq_j1_j2_matrix}`. By :math:`\eqref{eq_j1_j2_raises_or_lowers_state}` we can assume

	.. math::
		:nowrap:

		\begin{equation*}
			\left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} = \alpha_{\pm}(\jfrak, \sigma) \Psi^{\jfrak}_{\sigma \pm 1}
		\end{equation*}

	Applying :math:`\eqref{eq_j1_j2_mixed_product_plus}` and :math:`\eqref{eq_j1_j2_mixed_product_minus}` to :math:`\Psi^{\jfrak}_{\sigma}` then implies

	.. math::
		:nowrap:

		\begin{equation*}
			\alpha_{\mp} (\jfrak, \sigma \pm 1) \alpha_{\pm} (\jfrak, \sigma) = \jfrak^2 + \jfrak - \sigma^2 \mp \sigma
		\end{equation*}

	Now we use the fact that :math:`J^{(\jfrak)}_i, i = 1, 2, 3`, are Hermitian operators to calculate

	.. math::
		:nowrap:

		\begin{align*}
			|\alpha_{\pm} (\jfrak, \sigma)|^2 (\Psi^{\jfrak}_{\sigma}, \Psi^{\jfrak}_{\sigma}) &= \left( \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma}, \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} \right) \\
				&= \left( \Psi^{\jfrak}_{\sigma}, \left( J_1^{(\jfrak)} \mp \ifrak J_2^{(\jfrak)} \right) \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} \right) \\
				&= (\jfrak^2 + \jfrak - \sigma^2 \mp \sigma) (\Psi^{\jfrak}_{\sigma}, \Psi^{\jfrak}_{\sigma})
		\end{align*}

	It follows that, up to a choice of phase, :math:`\alpha_{\pm} (\jfrak, \sigma) = \sqrt{\jfrak^2 + \jfrak - \sigma^2 \mp \sigma} = \sqrt{(j \mp \sigma)(j \pm \sigma + 1)}`, which confirms :math:`\eqref{eq_j1_j2_matrix}`.

We end the discussion about massive particle states by working out the little group elements :math:`W(\Lambda, p)` defined by :math:`\eqref{eq_w_from_l}`. To this end, it suffices to work out the standard :math:`L(p)` such that :math:`L(p) k = p`, where :math:`k = (M, 0, 0, 0)`. We have already worked out such a transformation in :math:`\eqref{eq_lambda_boost}` and :math:`\eqref{eq_general_lambda_in_spacetime}` in spacetime coordinates, so we only need to translate it into :math:`4`-momentum coordinates.

Using :math:`\eqref{eq_p_from_v}`, we can rewrite :math:`\gamma` defined by :math:`\eqref{eq_def_gamma}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\pbf = \frac{M \vbf}{\sqrt{1 - \vbf^2}} \implies \gamma \coloneqq \frac{1}{\sqrt{1 - \vbf^2}} = \frac{\sqrt{M^2 + \pbf^2}}{M} \left( = \frac{p_0}{M} \right)
	\end{equation*}

It follows that

.. math::
	:nowrap:

	\begin{align}
		{L(p)_0}^0 &= \gamma \label{eq_L_transformation_for_massive_1} \\
		{L(p)_i}^0 = {L(p)^0}_i &= \frac{p_i}{M} \label{eq_L_transformation_for_massive_2} \\
		L(p)_{ij} &= \delta_{ij} + \frac{p_i p_j}{\pbf^2} (\gamma - 1) \label{eq_L_transformation_for_massive_3}
	\end{align}

Finally, we note an important fact that when :math:`\Lambda = \Rcal` is a :math:`3`-rotation, then

.. math::
	:nowrap:

	\begin{equation}
		W(\Rcal, p) = \Rcal  \label{eq_little_group_rotation}
	\end{equation}

for any :math:`p`. To see this, we'll work out how :math:`W(\Rcal, p)` acts on :math:`(1, \mathbf{0}), (0, \pbf)`, and :math:`(0, \qbf)`, respectively, where :math:`\qbf` is any :math:`3`-vector perpendicular to :math:`\pbf`, as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		W(\Rcal, p)(1, \mathbf{0}) &= L(\Rcal p)^{-1} \Rcal (\gamma, \pbf / M) &&= L(\Rcal p)^{-1} (\gamma, \Rcal p / M) &&= (1, \mathbf{0}) \\
		W(\Rcal, p)(0, \pbf) &= L(\Rcal p)^{-1} \Rcal (\pbf^2 / M, \gamma \pbf) &&= L(\Rcal p)^{-1} (\pbf^2 / M, \gamma \Rcal p) &&= (0, \Rcal \pbf) \\
		W(\Rcal, p)(0, \qbf) &= L(\Rcal p)^{-1} \Rcal (0, \qbf) &&= L(\Rcal p)^{-1} (0, \Rcal \qbf) &&= (0, \Rcal \qbf)
	\end{alignat*}

where we have used that fact that :math:`\gamma` is :math:`\Rcal`-invariant.

This observation is important since it implies that non-relativistic calculations about angular momenta, such as the `Clebsch-Gordan coefficients <https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients>`__, can be literally carried over to the relativistic setting.

.. _dropdown_clebsch_gordan_coefficients:

.. dropdown:: Clebsch-Gordan coefficients
	:icon: unlock
	:animate: fade-in-slide-down

	In a scenario where multiple particles present, or even just a single particle with both orbital angular momentum (i.e., the quantization of the classical angular momentum :math:`\xbf \times \Pbf`) and spin, it may happen that the full Hamiltonian doesn't commute with each individual :math:`3`-momentum :math:`\Jbf`, but commute with a "total" angular momentum. Therefore a formula, in terms of the so-called Clebsch-Gordan coefficients, that expresses the total angular momentum in terms of the individual ones is desirable. This section follows closely §4 from [Wei15]_. Note that the discussions that follow will be

	1. non-relativistic, which is justified by :math:`\eqref{eq_little_group_rotation}`, and
	2. applicable mostly (but not necessarily) to multi-particles states, rather than single-particle states.

	We shall focus on the composition of two angular momentum :math:`3`-vectors :math:`\Jbf'` and :math:`\Jbf''`, whether orbital or spin, that commute, i.e., they each satisfies :math:`\eqref{eq_jjj_commutation}` and in addition

	.. math::
		:nowrap:

		\begin{equation*}
			[J'_i, J''_j] = 0
		\end{equation*}

	for :math:`1 \leq i, j \leq 3`.

	Let's recollect the eigenstate representations :math:`\eqref{eq_j1_j2_matrix}, \eqref{eq_j3_matrix}` and :math:`\eqref{eq_3j_square_eigenstate}` as follows,

	.. math::
		:nowrap:

		\begin{align*}
			{\Jbf'}^2 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \jfrak' (\jfrak' + 1) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
			J'_3 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sigma' \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
			(J'_1 \pm \ifrak J'_2) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sqrt{{\jfrak'}^2 + \jfrak' - {\sigma'}^2 \mp \sigma'} ~\Psi^{\jfrak' ~\jfrak''}_{\sigma' \pm 1, \sigma''} \\
			{\Jbf''}^2 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \jfrak'' (\jfrak'' + 1) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
			J''_3 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sigma'' \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
			(J''_1 \pm \ifrak J''_2) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sqrt{{\jfrak''}^2 + \jfrak'' - {\sigma''}^2 \mp \sigma''} ~\Psi^{\jfrak' ~\jfrak''}_{\sigma', \sigma'' \pm 1}
		\end{align*}

	Without knowing exactly how the Hamiltonian :math:`H` looks like, we cannot really say what combinations of these angular momentum operators commute with :math:`H`, and therefore may be used to label states. However, one can imagine that a rotationally invariant Hamiltonian may contain terms like :math:`{\Jbf'}^2, {\Jbf''}^2` and interactions like :math:`\Jbf' \cdot \Jbf''`. In this case, we may choose to consider the following collection of (mutually commuting) operators

	.. math::
		:nowrap:

		\begin{equation*}
			{\Jbf'}^2, ~{\Jbf''}^2, ~\Jbf^2, \text{ and } J_3
		\end{equation*}

	where :math:`\Jbf \coloneqq \Jbf' + \Jbf''` is the total angular momentum, and :math:`J_i, i=1,2,3,` is its :math:`i`-th component.

	Now our goal is to express the eigenstates :math:`\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma}` which satisfy the following

	.. math::
		:nowrap:

		\begin{align}
			{\Jbf'}^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak' (\jfrak' + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \nonumber \\
			{\Jbf''}^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak'' (\jfrak'' + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \nonumber \\
			\Jbf^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak (\jfrak + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \nonumber \\
			J_3 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \sigma \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \nonumber \\
			(J_1 \pm \ifrak J_2) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \sqrt{\jfrak^2 + \jfrak - \sigma^2 \mp \sigma} ~\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma \pm 1} \label{eq_raising_lowering_am_pair}
		\end{align}

	in terms of :math:`\Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''}` as follows

	.. math::
		:nowrap:

		\begin{equation}
			\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} = \sum_{\sigma' \sigma''} C^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''}
			\label{eq_defn_clebsch_gordan_coefficients}
		\end{equation}

	where the coefficients are known as Clebsch-Gordan coefficients. For the clarity of exposition, let's divide the solution into a few steps.

	Step 1.
		First of all, note that since :math:`J_3 = J'_3 + J''_3`, we have the following constraint

		.. math::
			:nowrap:

			\begin{equation}
				C^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \neq 0 \implies \sigma = \sigma' + \sigma''
				\label{eq_sigma_additive}
			\end{equation}

		Moreover, we see that the maximum possible value of :math:`\sigma` is :math:`\jfrak' + \jfrak''`, and it's achieved exactly when :math:`\sigma' = \jfrak'` and :math:`\sigma'' = \jfrak''`. It follows, assuming the non-degeneracy of the representation at least, that

		.. math::
			:nowrap:

			\begin{equation}
				\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak''} =  \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''}
				\label{eq_highest_weight_am_pair}
			\end{equation}

		or equivalently

		.. math::
			:nowrap:

			\begin{equation}
				C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' ~\jfrak' + \jfrak''; \sigma' \sigma'') = \delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak''}
			\label{eq_clebsch_gordan_highest_weight}
			\end{equation}

	Step 2.
		Next consider a state with :math:`\sigma = \jfrak' + \jfrak'' - 1`. It follows from :math:`\eqref{eq_sigma_additive}` that it must be a superposition of :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''}` and :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1}` unless :math:`\jfrak'` and/or :math:`\jfrak''` vanishes, which leads to even simpler situations. Now we have two possible values of :math:`\jfrak`, namely :math:`\jfrak' + \jfrak''` and :math:`\jfrak' + \jfrak'' - 1`.

		In the former case, we can apply :math:`\eqref{eq_raising_lowering_am_pair}` by letting :math:`\sigma = \jfrak = \jfrak' + \jfrak''` as follows

		.. math::
			:nowrap:

			\begin{align*}
				\sqrt{2(\jfrak' + \jfrak'')} ~\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak'' - 1} \
					&\xlongequal{\eqref{eq_raising_lowering_am_pair}} (J_1 - \ifrak J_2) \Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak''} \\
					&\xlongequal{\eqref{eq_highest_weight_am_pair}} (J_1 - \ifrak J_2) \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''} \\
					&= (J'_1 - \ifrak J'_2 + J''_2 - \ifrak J''_2) \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''} \\
					&= \sqrt{2 \jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''} + \sqrt{2 \jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1}
			\end{align*}

		This gives us one of the :math:`\sigma = \jfrak' + \jfrak'' - 1` states

		.. math::
			:nowrap:

			\begin{equation}
				\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak'' - 1} = (\jfrak' + \jfrak'')^{-1/2} \left( \sqrt{\jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''} + \sqrt{\jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1} \right)
				\label{eq_second_highest_weight_am_pair_one}
			\end{equation}

		The other one, which has :math:`\jfrak = \jfrak' + \jfrak'' - 1`, must be orthogonal to :math:`\eqref{eq_second_highest_weight_am_pair_one}`. Therefore up to a normalization factor, we can write

		.. math::
			:nowrap:

			\begin{equation}
				\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak'' - 1}_{\jfrak' + \jfrak'' - 1} = (\jfrak' + \jfrak'')^{-1/2} \left( \sqrt{\jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak'-1 ~\jfrak''} - \sqrt{\jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1} \right)
			\label{eq_second_highest_weight_am_pair_two}
			\end{equation}

		We can translate :math:`\eqref{eq_second_highest_weight_am_pair_one}` and :math:`\eqref{eq_second_highest_weight_am_pair_two}` into Clebsch-Gordan coefficients as follows

		.. math::
			:nowrap:

			\begin{align}
				C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' ~\jfrak' + \jfrak'' - 1; \sigma' \sigma'') &= \
					\sqrt{\frac{\jfrak'}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak' - 1} \delta_{\sigma'' ~\jfrak''} + \
					\sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak'' - 1} \
					\label{eq_clebsch_gordan_second_highest_weight_one} \\

				C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' -1 ~\jfrak' + \jfrak'' - 1; \sigma' \sigma'') &= \
					\sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak' - 1} \delta_{\sigma'' ~\jfrak''} - \
					\sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} \delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak'' - 1} \
					\label{eq_clebsch_gordan_second_highest_weight_two}
			\end{align}

	Step 3.
		The pattern should now be clear. Namely for :math:`\sigma = \jfrak' + \jfrak'' - 2`, the :math:`J_3`-eigenspace must be :math:`3`-dimensional and spanned by :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 2 ~\jfrak''}, \Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak'' - 1}` and :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 2}`. Two of them come from :math:`\eqref{eq_second_highest_weight_am_pair_one}, \eqref{eq_second_highest_weight_am_pair_two}` by applying the lowering operator :math:`J_1 - \ifrak J_2`, and the third orthogonally complements the first two.

		This procedure can be continued but will eventually terminate because of the bounds :math:`|\sigma'| \leq \jfrak'` and :math:`|\sigma''| \leq \jfrak''`. It follows that :math:`\jfrak` can only take the following values

		.. math::
			:nowrap:

			\begin{equation}
				\jfrak = |\jfrak' - \jfrak''|, ~|\jfrak' - \jfrak''| + 1, \cdots, ~\jfrak' + \jfrak''
				\label{eq_composite_total_angular_momentum_range}
			\end{equation}

	To make the above rather abstract calculations more concrete, let's consider the example of a hydrogen atom, where the single electron is subject to the (radial) Coulomb force. Without actually spelling out the Hamiltonian, let's jump directly to the result. The result is that the energy levels can be labeled by positive integers :math:`n = 1, 2, \cdots`. For each :math:`n`, the orbital angular momentum :math:`\ell` can take any value from :math:`0, 1, \cdots, n - 1`, traditionally labeled by :math:`s, p, d, f, g, \cdots`. Finally the electron spin is :math:`\sfrak = 1/2`. We can then make the following table for the first few states with :math:`n \leq 2`.

	+------------------+-----------+------------------------+---------------------------+----------------+-----------------+-----------------+------------------+----------------------------------------------------------+
	| State            | :math:`n` | :math:`\ell(=\jfrak')` | :math:`\sfrak(=\jfrak'')` | :math:`\jfrak` | :math:`\sigma`  | :math:`\sigma'` | :math:`\sigma''` | :math:`C^{\ell~\sfrak}(\jfrak~\sigma;~\sigma' \sigma'')` |
	+==================+===========+========================+===========================+================+=================+=================+==================+==========================================================+
	| :math:`1s_{1/2}` | :math:`1` | :math:`0`              | :math:`1/2`               | :math:`1/2`    | :math:`\pm 1/2` | :math:`0`       | :math:`\pm 1/2`  | :math:`1`                                                |
	+------------------+-----------+------------------------+---------------------------+----------------+-----------------+-----------------+------------------+----------------------------------------------------------+
	| :math:`2s_{1/2}` | :math:`2` | :math:`0`              | :math:`1/2`               | :math:`1/2`    | :math:`\pm 1/2` | :math:`0`       | :math:`\pm 1/2`  | :math:`1`                                                |
	+------------------+-----------+------------------------+---------------------------+----------------+-----------------+-----------------+------------------+----------------------------------------------------------+
	| :math:`2p_{1/2}` | :math:`2` | :math:`1`              | :math:`1/2`               | :math:`1/2`    | :math:`\pm 1/2` | :math:`0`       | :math:`\pm 1/2`  | :math:`\pm \sqrt{1/3}`                                   |
	|                  |           |                        |                           |                |                 +-----------------+------------------+----------------------------------------------------------+
	|                  |           |                        |                           |                |                 | :math:`\pm 1`   | :math:`\mp 1/2`  | :math:`\mp \sqrt{2/3}`                                   |
	+------------------+-----------+------------------------+---------------------------+----------------+-----------------+-----------------+------------------+----------------------------------------------------------+
	| :math:`2p_{3/2}` | :math:`2` | :math:`1`              | :math:`1/2`               | :math:`3/2`    | :math:`\pm 3/2` | :math:`\pm 1`   | :math:`\pm 1/2`  | :math:`1`                                                |
	|                  |           |                        |                           |                +-----------------+-----------------+------------------+----------------------------------------------------------+
	|                  |           |                        |                           |                | :math:`\pm 1/2` | :math:`0`       | :math:`\pm 1/2`  | :math:`\sqrt{2/3}`                                       |
	|                  |           |                        |                           |                |                 +-----------------+------------------+----------------------------------------------------------+
	|                  |           |                        |                           |                |                 | :math:`\pm 1`   | :math:`\mp 1/2`  | :math:`\sqrt{1/3}`                                       |
	+------------------+-----------+------------------------+---------------------------+----------------+-----------------+-----------------+------------------+----------------------------------------------------------+

	where a state is structured as :math:`\{n\}\{\ell\}_{\{\jfrak\}}`. Interpreting the Clebsch-Gordan coefficients as probability amplitudes, we see for example that if one measures the :math:`3`-components of the orbital and spin of the electron in a hydrogen atom in state :math:`2p_{1/2}` with :math:`\sigma = 1/2`, then one gets either :math:`(0, 1/2)` or :math:`(1, -1/2)` with probabilities :math:`1/3` and :math:`2/3`, respectively.

	We end this (long) excursion with a few interesting facts that we'll come back to later.

	3. The energy difference between :math:`1s_{1/2}` and :math:`2s_{1/2}`, plotted on a spectrometer, is the famous `21-centimeter line <https://en.wikipedia.org/wiki/Hydrogen_line>`_.
	4. The energy difference between :math:`2p_{1/2}` and :math:`2p_{3/2}`, i.e., same orbital but different total angular momentum, is known as the `fine structure <https://en.wikipedia.org/wiki/Fine_structure>`_ of the hydrogen atom.
	5. The energy difference between :math:`2s_{1/2}` and :math:`2p_{1/2}`, i.e., same total but different orbital angular momentum, is known as the `Lamb shift <https://en.wikipedia.org/wiki/Lamb_shift>`_.
	6. The energy difference between states with the same orbital and total angular momentum, e.g., :math:`1s_{1/2}`, but different spin :math:`z`-component :math:`\sigma`, e.g., :math:`\pm 1/2`, due to the magnetic moment is known as the `hyperfine structure <https://en.wikipedia.org/wiki/Hyperfine_structure>`_.


.. _sec_massless_particle_states:

Massless particle states
^^^^^^^^^^^^^^^^^^^^^^^^

Recall the standard :math:`k = (1, 0, 0, 1)` for massless particles. Our first task is to work out the little group, i.e., Lorentz transformations :math:`W` such that :math:`Wk = k`. More precisely, we'll work out the column vectors of :math:`W` by thinking of them as the results of :math:`W` acting on the standard basis vectors. Let's start by :math:`v \coloneqq (1, 0, 0, 0)`, and perform the following calculations to :math:`Wv` using properties of Lorentz transformation

.. math::
	:nowrap:

	\begin{alignat}{2}
		(Wv)^{\mu} (Wv)_{\mu} &= v^{\mu} v_{\mu} &&= 1  \label{eq_vv_is_one} \\
		(Wv)^{\mu} k_\mu &= v^{\mu} k_{\mu} &&= 1  \label{eq_vk_is_one}
	\end{alignat}

It follows from :math:`\eqref{eq_vk_is_one}` that we can write :math:`Wv = (1 + c, a, b, c)`, and then from :math:`\eqref{eq_vv_is_one}` that :math:`c = (a^2 + b^2) / 2`. Playing similar games to the other basis vectors, we can engineer a particular Lorentz transformation as follows

.. math::
	:nowrap:

	\begin{equation}
		{S_{\mu}}^{\nu}(a, b) = \begin{bmatrix*}[r]
			1 + c & a & b & -c \\
			a & 1 & 0 & -a \\
			b & 0 & 1 & -b \\
			c & a & b & 1 - c
		\end{bmatrix*}
		\label{eq_massless_little_group_s_matrix}
	\end{equation}

which leaves :math:`k` invariant, and satisfies :math:`Sv = Wv`. It follows that :math:`S^{-1} W` must be a rotation about the :math:`3`-axis, which can be written as follows

.. math::
	:nowrap:

	\begin{equation}
		R(\theta) = \begin{bmatrix*}[r]
			1 & 0 & 0 & 0 \\
			0 & \cos\theta & \sin\theta & 0 \\
			0 & -\sin\theta & \cos\theta & 0 \\
			0 & 0 & 0 & 1
		\end{bmatrix*}
		\label{eq_massless_little_group_r_matrix}
	\end{equation}

Hence we can write any element in the little group as :math:`W(a, b, \theta) = S(a, b) R(\theta)`.

.. dropdown:: The little group is :math:`~E^+(2)`
	:icon: unlock
	:animate: fade-in-slide-down

	Although not necessary for our purposes here, we'd like to better understand the little group for :math:`k = (1, 0, 0, 1)` in terms of more familiar groups. It turns out that it's isomorphic to the :math:`2`-dimensional orientation-preserving `Euclidean group <https://en.wikipedia.org/wiki/Euclidean_group>`__ :math:`E^+(2)`, i.e., the group of rotations and translations on the plane.

	To see this, we go back to the defining property of :math:`W` that it fixes :math:`k`. It follows that it must also fix the orthogonal complement :math:`k^{\bot}` with respect to the bilinear form :math:`d\tau^2` defined in :math:`\eqref{eq_proper_time}`. Since :math:`k` is orthogonal to itself, we can uniquely determine :math:`W` by knowing its action on :math:`(1, 1, 0, 1)` and :math:`(1, 0, 1, 1)`. Letting :math:`S(a, b)` act on them, we see

	.. math::
		:nowrap:

		\begin{align*}
			S(a, b)(1, 1, 0, 1) &= (1 + a, 1, 0, 1 + a) \\
			S(a, b)(1, 0, 1, 1) &= (1 + b, 0, 1, 1 + b)
		\end{align*}

	Hence :math:`S` is isomorphic to a :math:`2`-dimensional translation group. Moreover, the direction of translation is determined by the rotation on the plane spanned by the second the the third coordinates, which is nothing but :math:`R`.

	Note that :math:`E^+ (2)` is not semisimple in the sense that it possesses an abelian normal subgroup. Indeed, it's obvious from the above discussion that a translation conjugated by a rotation is again a translation (in the rotated direction). The non-semisimplicity will have consequences on the representation as we will see below.

As in the massive case, we'll work out :math:`D_{\sigma \sigma'}` up to first order. To this end, note that up to first order

.. math::
	:nowrap:

	\begin{align*}
		W(a, b, \theta)_{\mu}^{\nu} &= \left(1 + \begin{bmatrix*}[r]
				0 & a & b & 0 \\
				a & 0 & 0 & -a \\
				b & 0 & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix*} \right) \left(1 + \begin{bmatrix*}[r]
				0 & 0 & 0 & 0 \\
				0 & 0 & \theta & 0 \\
				0 & -\theta & 0 & 0 \\
				0 & 0 & 0 & 0
			\end{bmatrix*} \right) + \cdots \\
			&= 1 + \begin{bmatrix*}[r]
				0 & a & b & 0 \\
				a & 0 & \theta & -a \\
				b & -\theta & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix*} + \cdots
	\end{align*}

where we've added the :math:`4`-indexes since we recall from discussions in :ref:`sec_quantum_lorentz_symmetry` that we must lift the :math:`\omega` index to make it anti-symmetric. We now rewrite

.. math::
	:nowrap:

	\begin{equation*}
		W(a, b , \theta)^{\mu \nu} = \eta^{\mu \sigma} W_{\sigma}^{\nu} = 1 + \begin{bmatrix*}[r]
				0 & -a & -b & 0 \\
				a & 0 & \theta & -a \\
				b & -\theta & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix*} + \cdots
	\end{equation*}

and spell out the expansion of :math:`D(a, b, \theta) \coloneqq D(W(a, b, \theta))` as follows

.. math::
	:nowrap:

	\begin{equation}
		D(a, b, \theta) = 1 + \ifrak aA + \ifrak bB + \ifrak \theta J_3
		\label{eq_massless_D_matrix_expansion}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{alignat*}{2}
		A &= -J_{01} - J_{13} &&= -K_1 + J_2 \\
		B &= -J_{02} - J_{23} &&= -K_2 - J_1
	\end{alignat*}

Next we use :math:`\eqref{eq_jjj_commutation}, \eqref{eq_jkk_commutation}` and :math:`\eqref{eq_kkj_commutation}` to calculate commutation relations between :math:`A, B` and :math:`J_3` as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		[J_3, A] &= -&&\ifrak K_2 &&- \ifrak J_1 &&= \ifrak B \\
		[J_3, B] &= &&\ifrak K_1 &&- \ifrak J_2 &&= -\ifrak A \\
		[A, B] &= -&&\ifrak J_3 &&+ \ifrak J_3 &&= 0
	\end{alignat*}

Since :math:`A, B` commute, we can use their eigenvalues to label states as follows

.. math::
	:nowrap:

	\begin{align*}
		A \Psi_{k, a, b} &= a \Psi_{k, a, b} \\
		B \Psi_{k, a, b} &= b \Psi_{k, a, b}
	\end{align*}

In fact, these states, corresponding to translation symmetries, come in continuous families as shown below

.. math::
	:nowrap:

	\begin{align*}
		AU^{-1}(R(\theta)) \Psi_{a, b, k} &= (a\cos\theta - b\sin\theta)U^{-1}(R(\theta)) \Psi_{a, b, k} \\
		BU^{-1}(R(\theta)) \Psi_{a, b, k} &= (a\sin\theta + b\cos\theta)U^{-1}(R(\theta)) \Psi_{a, b, k}
	\end{align*}

According to [Wei95]_ (page 72), massless particle states are not observed to come in such :math:`S^1`-families. Hence the only possibility is that :math:`a = b = 0` and the only symmetry left then is :math:`J_3`, which corresponds to a rotation about the :math:`3`-axis.

Unlike the :math:`SO(3)`-symmetry discussed in :ref:`Representations of angular momenta <dropdown_repr_of_angular_momenta>`, representations of :math:`J_3` alone cannot be characterized at the infinitesimal level, which would have resulted in a continuous spectrum. Instead, since a :math:`2\pi`-rotation about the :math:`3`-axis gives the identity transformation, one might expect an integer spectrum for :math:`J_3`. This is indeed the case if we :ref:`assume the representation is genuine <assump_genuine_repr>`. However, since the Lorentz group is not simplify connected (with fundamental group :math:`\Zbb/2`), one may encounter projective representations. Indeed, the :math:`2\pi`-rotation about the :math:`3`-axis represents a generator of the fundamental group, which has order :math:`2`, i.e., only the :math:`4\pi`-rotation about the :math:`3`-axis represents a contractible loop in the Lorentz group (see the `Plate trick <https://en.wikipedia.org/wiki/Plate_trick>`_). As a result, the :math:`J_3`-spectrum actually consists of half-integers, just like the spins. We can therefore write a general massless particle state as :math:`\Psi_{k, \sigma}` such that

.. math::
	:nowrap:

	\begin{equation*}
		J_3 \Psi_{k, \sigma} = \sigma \Psi_{k, \sigma}
	\end{equation*}

where :math:`\sigma` are half-integers, known as the *helicity*.

Combining the discussions so far, we can write down the :math:`D`-matrix defined by :math:`\eqref{eq_massless_D_matrix_expansion}` as follows

.. math::
	:nowrap:

	\begin{equation}
		D_{\sigma \sigma'}(W(a, b, \theta)) = \exp(\ifrak \theta \sigma) \delta_{\sigma \sigma'}
		\label{eq_little_group_d_matrix_massless}
	\end{equation}

where we recall :math:`W(a, b, \theta) = L(\Lambda p)^{-1} \Lambda L(p) = S(a, b)R(\theta)`. The Lorentz transformation formula :math:`\eqref{eq_lorentz_transformation_formula_for_particle_state}` for massless particles now becomes

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \Psi_{p, \sigma} = \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(\ifrak \theta(\Lambda, p) \sigma) \Psi_{\Lambda p, \sigma}
		\label{eq_lorentz_transformation_formula_for_massless}
	\end{equation}

In particular, we see that, unlike the spin :math:`z`-component of massive particles, helicity is Lorentz invariant (at least under genuine representations). It is reasonable, therefore, to think of massless particles of different helicity as different particle species. Examples include photons with :math:`\sigma = \pm 1` and gravitons with :math:`\sigma = \pm 2`, but *not* (anti-)neutrinos with hypothetical :math:`\sigma = \pm \tfrac{1}{2}` as otherwise stated in [Wei95]_ (page 73 -- 74), which are now known to have a nonzero mass. Here the :math:`\pm` signs are related to the space-inversion symmetry :math:`\eqref{eq_space_inversion}`, which will be discussed in detail later.

In order to use :math:`\eqref{eq_lorentz_transformation_formula_for_massless}` for a general :math:`(\Lambda, p)`, we first need to fix the choices of :math:`L(p)` that takes the standard :math:`k = (1, 0, 0, 1)` to :math:`p`. This can be done in two steps. First apply a (pure) boost along the :math:`3`-axis

.. math::
	:nowrap:

	\begin{equation}
		\begin{bmatrix}
			(p_0^2 + 1) / 2p_0 & 0 & 0 & (p_0^2 - 1) / 2p_0 \\
			0 & 1 & 0 & 0 \\
			0 & 0 & 1 & 0 \\
			(p_0^2 - 1) / 2p_0 & 0 & 0 & (p_0^2 + 1) / 2p_0
		\end{bmatrix}
		\begin{bmatrix}
			1 \\
			0 \\
			0 \\
			1
		\end{bmatrix} = \begin{bmatrix}
			p_0 \\
			0 \\
			0 \\
			p_0
		\end{bmatrix} = \begin{bmatrix}
			p_0 \\
			0 \\
			0 \\
			|\pbf|
		\end{bmatrix}
		\label{eq_massless_boost}
	\end{equation}

Then apply a (pure) rotation that takes :math:`(0, 0, |\pbf|)` to :math:`\pbf`. However, in contrast to the massive case :math:`\eqref{eq_L_transformation_for_massive_1}` -- :math:`\eqref{eq_L_transformation_for_massive_3}`, where :math:`L(p)` depends continuously on :math:`p`, there exists no continuous family of rotations that take :math:`(0, 0, |\pbf|)` to any other :math:`3`-vector (of the same length). Fortunately, any two choices of such rotations differ by (a pre-composition of) a rotation about the :math:`3`-axis, which, according to :math:`\eqref{eq_lorentz_transformation_formula_for_massless}`, only produces a physically immaterial phase factor.

.. dropdown:: Polarization of photons
	:icon: unlock
	:animate: fade-in-slide-down

	General photon states of definite momentum can be written as a superposition

	.. math::
		:nowrap:

		\begin{equation*}
			\Psi_{p, \alpha} \coloneqq \alpha_+ \Psi_{p, +1} + \alpha_- \Psi_{p, -1}
		\end{equation*}

	such that :math:`|\alpha_+|^2 + |\alpha_-|^2 = 1`. They are not in general Lorentz invariant due to the angle :math:`\theta` presented in the phase factor in :math:`\eqref{eq_lorentz_transformation_formula_for_massless}`. Here the base states :math:`\Psi_{p, \pm 1}` are known as (right and left-handedly) *circularly polarized*, and the other extreme cases where :math:`|\alpha_+| = |\alpha_-|` are known as *linearly polarized*. All other intermediate cases are then called *elliptically polarized*.

	It's not obvious at all why these states are named the way they are, if only viewed as abstract combinations of eigenstates of an abstract operator :math:`J_3`. They are named after analogies, either with `helicity <https://en.wikipedia.org/wiki/Helicity_(particle_physics)>`_ from classical mechanics or with classical electromagnetic fields (in vacuum) from Maxwell's theory.

.. _sec_space_and_time_inversions:

Space and time inversions
^^^^^^^^^^^^^^^^^^^^^^^^^

So far the discussions have been focused on orthochronous (and mostly homogeneous) Lorentz transformations, and in particular, the infinitesimal symmetries at the vicinity of the identity. Now it's time to take a look at the space and time inversions, defined in :math:`\eqref{eq_space_inversion}` and :math:`\eqref{eq_time_inversion}`, which takes us to the other components of the Lorentz group. The main goal is to understand their actions on the one-particle states, that have been worked out in the previous two sections.

Let's write

.. math::
	:nowrap:

	\begin{equation*}
		U(\Pcal) \coloneqq U(\Pcal, 0), \quad U(\Tcal) \coloneqq U(\Tcal, 0)
	\end{equation*}

for the corresponding quantum symmetry operators, which we haven't decided whether should be complex linear or anti-linear. The same calculations that led to :math:`\eqref{eq_p_conjugated_by_u}` and :math:`\eqref{eq_j_conjugated_by_u}` now give

.. math::
	:nowrap:

	\begin{align}
		U(\Pcal) \ifrak P_{\mu} U^{-1}(\Pcal) &= \ifrak \Pcal^{\rho}_{\mu} P_{\rho} \label{eq_p_conjugated_by_p} \\
		U(\Pcal) \ifrak J_{\mu \nu} U^{-1}(\Pcal) &= \ifrak \Pcal^{\rho}_{\mu} \Pcal^{\kappa}_{\nu} J_{\rho \kappa} \label{eq_j_conjugated_by_p} \\
		U(\Tcal) \ifrak P_{\mu} U^{-1}(\Tcal) &= \ifrak \Tcal^{\rho}_{\mu} P_{\rho} \label{eq_p_conjugated_by_t} \\
		U(\Tcal) \ifrak J_{\mu \nu} U^{-1}(\Tcal) &= \ifrak \Tcal^{\rho}_{\mu} \Tcal^{\kappa}_{\nu} J_{\rho \kappa} \label{eq_j_conjugated_by_t}
	\end{align}

The complex (anti-)linearity of :math:`U(\Pcal)` and :math:`U(\Tcal)` can then be decided by the postulation that physically meaningful energy must not be negative. More precisely, recall that :math:`P_0` is the energy operator. Then :math:`\eqref{eq_p_conjugated_by_p}` shows

.. math::
	:nowrap:

	\begin{equation*}
		U(\Pcal) \ifrak P_0 U^{-1}(\Pcal) = \ifrak P_0
	\end{equation*}

If :math:`U(\Pcal)` were anti-linear, then :math:`U(\Pcal) P_0 U^{-1}(\Pcal) = -P_0`. Then for any state :math:`\Psi` with positive energy, i.e., :math:`P_0 \Psi = p_0 \Psi`, we would have a state :math:`U^{-1}(\Pcal) \Psi` with negative energy :math:`-p_0`. Hence we conclude that :math:`U(\Pcal)` must be linear. The same argument shows also that :math:`U(\Tcal)` must be anti-linear (since :math:`\Tcal_0^0 = -1`).

As before, it'll be useful to rewrite :math:`\eqref{eq_p_conjugated_by_p}` -- :math:`\eqref{eq_j_conjugated_by_t}` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:nowrap:

	\begin{alignat}{3}
		U(\Pcal) &H U^{-1}(\Pcal) &&= &&H \nonumber \\
		U(\Pcal) &\Pbf U^{-1}(\Pcal) &&= -&&\Pbf \label{eq_p3_conjugated_by_p} \\
		U(\Pcal) &\Jbf U^{-1}(\Pcal) &&= &&\Jbf \label{eq_j3_conjugated_by_p} \\
		U(\Pcal) &\Kbf U^{-1}(\Pcal) &&= -&&\Kbf \label{eq_k3_conjugated_by_p} \\
		U(\Tcal) &H U^{-1}(\Tcal) &&= &&H \nonumber \\
		U(\Tcal) &\Pbf U^{-1}(\Tcal) &&= -&&\Pbf \label{eq_p3_conjugated_by_t} \\
		U(\Tcal) &\Jbf U^{-1}(\Tcal) &&= -&&\Jbf \label{eq_j3_conjugated_by_t} \\
		U(\Tcal) &\Kbf U^{-1}(\Tcal) &&= &&\Kbf \nonumber \\
	\end{alignat}

One can (and should) try to reconcile these implications with commonsense. For example, :math:`\eqref{eq_p3_conjugated_by_p}` and :math:`\eqref{eq_p3_conjugated_by_t}` say that the :math:`3`-momentum changes direction under either space or time inversion, which is of course as expected. Moreover :math:`\eqref{eq_j3_conjugated_by_p}` says that the spin (of for example a basketball) remains the same under space inversion because both the direction of the axis and the handedness of the rotation get reversed simultaneously, but it gets reversed under a time inversion according to :math:`\eqref{eq_j3_conjugated_by_t}` because the direction of rotation is reversed if time flows backwards.

In what follows we will work out the effects of space and time inversions on massive and massless particles, respectively.


.. _sec_space_inversion_for_massive_particles:

Space inversion for massive particles
+++++++++++++++++++++++++++++++++++++

We start by considering a state at rest :math:`\Psi_{k, \sigma}`, where :math:`k = (M, 0, 0, 0)` and :math:`\sigma` is an eigenvalue of :math:`J_3` under one of the spin representations discussed in :ref:`Representations of angular momenta <dropdown_repr_of_angular_momenta>`. Since the state is at rest and :math:`U(\Pcal)` commutes with :math:`J_3` according to :math:`\eqref{eq_j3_conjugated_by_p}`, we can write

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \Psi_{k, \sigma} = \eta \Psi_{k, \sigma}
		\label{eq_space_inversion_on_massive_standard}
	\end{equation}

where :math:`\eta` is a phase that depends a priori on :math:`\sigma`. It turns out, however, that :math:`\eta` is actually independent of :math:`\sigma`, and hence justifies the notation, since :math:`U(\Pcal)` commutes with the raising/lowering operators :math:`J_1 \pm \ifrak J_2` by :math:`\eqref{eq_j3_conjugated_by_p}`.

To move on to the general case, we recall that the general formula :math:`\eqref{eq_def_of_one_particle_psi}` takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		\Psi_{p, \sigma} = \sqrt{\frac{M}{p_0}} U(L(p)) \Psi_{k, \sigma}
	\end{equation*}

We can calculate as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \Psi_{p, \sigma} = \sqrt{\frac{M}{p_0}} U(\Pcal L(p) \Pcal^{-1}) U(\Pcal) \Psi_{k, \sigma} = \eta~\sqrt{\frac{M}{p_0}} U(L(\Pcal p)) \Psi_{k, \sigma} = \eta \Psi_{\Pcal p, \sigma}
		\label{eq_space_inversion_on_massive_general}
	\end{equation}

which generalizes :math:`\eqref{eq_space_inversion_on_massive_standard}`. Such :math:`\eta` is known as the *intrinsic parity*, which is intrinsic to a particle species.


.. _sec_time_inversion_for_massive_particles:

Time inversion for massive particles
++++++++++++++++++++++++++++++++++++

Consider the same :math:`\Psi_{k, \sigma}` as in the space inversion case. Now since :math:`U(\Tcal)` anti-commutes with :math:`J_3` according to :math:`\eqref{eq_j3_conjugated_by_t}`, we can write

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \Psi_{k, \sigma} = \zeta_{\sigma} \Psi_{k, -\sigma}
	\end{equation*}

where :math:`\zeta_{\sigma}` is a phase. Applying the raising/lowering operators and using :math:`\eqref{eq_j1_j2_matrix}`, we can calculate the left-hand-side, recalling that :math:`U(\Tcal)` is anti-linear, as follows

.. math::
	:nowrap:

	\begin{align*}
		(J_1 \pm \ifrak J_2) U(\Tcal) \Psi_{k, \sigma} &= -U(\Tcal) (J_1 \mp \ifrak J_2) \Psi_{k, \sigma} = -U(\Tcal) \sqrt{(\jfrak \pm \sigma)(\jfrak \mp \sigma + 1)} \Psi_{k, \sigma \mp 1} \\
		&= -\zeta_{\sigma \mp 1} \sqrt{(\jfrak \pm 1)(\jfrak \mp \sigma + 1)} \Psi_{k, -\sigma \pm 1}
	\end{align*}

where :math:`\jfrak` is the particle spin, and the right-hand-side as follows

.. math::
	:nowrap:

	\begin{equation*}
		(J_1 \pm \ifrak J_2) \zeta_{\sigma} \Psi_{k, -\sigma} = \zeta_{\sigma} \sqrt{(\jfrak \pm 1)(\jfrak \mp \sigma + 1)} \Psi_{k, -\sigma \pm 1}
	\end{equation*}

Equating the two sides, we see that :math:`\zeta_{\sigma} = -\zeta_{\sigma \pm 1}`. Up to an overall phase, we can set :math:`\zeta_{\sigma} = \zeta (-1)^{\jfrak - \sigma}` so that

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \Psi_{k, \sigma} = \zeta (-1)^{\jfrak - \sigma} \Psi_{k, -\sigma}.
	\end{equation*}

Here we have chosen to keep the option of a physically inconsequential phase :math:`\zeta` open. As in the case of space inversion, the formula generalizes to any :math:`4`-momentum :math:`p`

.. math::
	:nowrap:

	\begin{equation}
		U(\Tcal) \Psi_{p, \sigma} = \zeta (-1)^{\jfrak - \sigma} \Psi_{\Pcal p, -\sigma}
		\label{eq_time_inversion_on_massive_general}
	\end{equation}

since :math:`\Tcal L(p) \Tcal^{-1} = L(\Pcal p)`.

Space inversion for massless particles
++++++++++++++++++++++++++++++++++++++

Let's consider a state :math:`\Psi_{k, \sigma}` with :math:`k = (1, 0, 0, 1)` and :math:`\sigma` being the helicity, i.e., :math:`J_3 \Psi_{k, \sigma} = \sigma \Psi_{k, \sigma}`. Since :math:`U(\Pcal)` commutes with :math:`J_3`, the space inversion preserves :math:`\sigma`, just as in the massive case. However, since :math:`\Pcal` reverses the direction of motion, the helicity in the direction of motion actually reverses sign. It follows, in particular, that (massless) particles that respect the space inversion symmetry must come in companion with another particle of opposite helicity.

To spell out more details, note that since :math:`\Pcal` doesn't fix :math:`k`, it'll be convenient to introduce an additional rotation :math:`R_2`, which is defined to be a :math:`\pi`-rotation about the :math:`2`-axis, so that :math:`U(R_2) = \exp(\ifrak \pi J_2)` and :math:`R_2 \Pcal k = k`. Since :math:`U(R_2)` flips the sign of :math:`J_3`, as can be seen from the very definition of :math:`J_3` in :math:`\eqref{eq_u_lorentz_expansion}`, we have

.. math::
	:nowrap:

	\begin{equation*}
		U(R_2 \Pcal) \Psi_{k, \sigma} = \eta_{\sigma} \Psi_{k, -\sigma}
	\end{equation*}

where we see indeed that the helicity reverses sign (when :math:`k` is fixed).

To move on to the general case, recall that the :math:`L(p)` that takes :math:`k` to :math:`p` consists of a boost :math:`B` defined by :math:`\eqref{eq_massless_boost}` followed by a (chosen) pure rotation :math:`R(\pbf)` that takes :math:`(0, 0, |\pbf|)` to :math:`\pbf`. We calculate as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Pcal) \Psi_{p, \sigma} &= p_0^{-1/2} U(\Pcal R(\pbf)B) \Psi_{k, \sigma} \label{eq_space_inversion_on_massless_undetermined_phase} \\
			&= p_0^{-1/2} U(R(\pbf) B R_2^{-1}) U(R_2 \Pcal) \Psi_{k, \sigma} \nonumber \\
			&= p_0^{-1/2} \eta_{\sigma} U(R(\pbf) R_2^{-1} B) \Psi_{k, -\sigma} \nonumber \\
			&= \eta_{\sigma} \rho \Psi_{\Pcal p, -\sigma} \nonumber
	\end{align}

where :math:`\rho` is an extra phase due to the fact that although :math:`R(\pbf) R_2^{-1}` takes :math:`(0, 0, |\pbf|)` to :math:`-\pbf`, it may not be the chosen one.

To spell out :math:`\rho`, we need to be a bit more specific about :math:`R(\pbf)`. Following the usual convention of `spherical coordinates <https://en.wikipedia.org/wiki/Spherical_coordinate_system>`_, we can get from :math:`(0, 0, |\pbf|)` to :math:`\pbf` by first rotate (according to the right-handed rule) about the :math:`1`-axis at an angle :math:`0 \leq \phi \leq \pi`, known as the polar angle, and then rotate about the :math:`3`-axis at an angle :math:`0 \leq \theta < 2\pi`, known as the azimuthal angle. Now since we know that :math:`R(\pbf)R_2^{-1}` differs from :math:`R(-\pbf)` by a rotation about the :math:`3`-axis, we can figure the rotation out by examining their actions on some suitably generic :math:`\pbf`, for example, the unit vector along the :math:`2`-axis, which is fixed by :math:`R_2`. In this case :math:`R(-\pbf)` is a :math:`\pi/2`-rotation about the :math:`1`-axis, while :math:`R(\pbf)R_2^{-1}` is the same :math:`\pi/2`-rotation by the :math:`1`-axis, followed by a a :math:`\pi`-rotation about the :math:`3`-axis. Therefore we conclude that the difference is a :math:`\pi`-rotation about the :math:`3`-axis. In other words, we should have :math:`\rho = \exp(-\ifrak \pi \sigma)`. However, recalling that the helicity :math:`\sigma` may be a half-integer (thought not yet being found in Nature), there is a sign difference between :math:`\pm \pi`-rotations about the :math:`3`-axis. Without going into further details, we write down the general formula the space inversion as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \Psi_{p, \sigma} = \eta_{\sigma} \exp(\mp \ifrak \pi \sigma) \Psi_{\Pcal p, -\sigma}
		\label{eq_space_inversion_on_massless_general}
	\end{equation}

where the sign depends on the sign of :math:`p_2` (which can be seen by playing the same game as above with the negative unit vector along the :math:`2`-axis).

Time inversion for massless particles
+++++++++++++++++++++++++++++++++++++

Let :math:`k = (1, 0, 0, 1)` as usual and consider the state :math:`\Psi_{k, \sigma}`. Since :math:`U(\Tcal)` anti-commutes with both :math:`\Pbf` and :math:`J_3` by :math:`\eqref{eq_p3_conjugated_by_t}` and :math:`\eqref{eq_j3_conjugated_by_t}`, we have

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \Psi_{k, \sigma} = \Psi_{\Pcal k, -\sigma}.
	\end{equation*}

Composing with the rotation :math:`R_2` as in the previous section to fix :math:`k`, we have

.. math::
	:nowrap:

	\begin{equation*}
		U(R_2 \Tcal) \Psi_{k, \sigma} = \zeta_{\sigma} \Psi_{k, \sigma}
	\end{equation*}

where :math:`\zeta_{\sigma}` is (yet another) phase. We see that, unlike the space inversion, the time inversion doesn't produce a doublet of opposite helicity. Processing as in the space inversion case, one can derive the following general formula similar to :math:`\eqref{eq_space_inversion_on_massless_general}`

.. math::
	:nowrap:

	\begin{equation}
		U(\Tcal) \Psi_{p, \sigma} = \zeta_{\sigma} \exp(\pm \ifrak \pi \sigma) \Psi_{\Pcal p, \sigma}
		\label{eq_time_inversion_on_massless_general}
	\end{equation}

where the sign depends on the sign of :math:`p_2` as before.

Kramers' degeneracy
+++++++++++++++++++

We end our discussion about space and time inversions of one-particle states with an interesting observation on the squared time inversion :math:`U(\Tcal)^2`. It follows from both :math:`\eqref{eq_time_inversion_on_massive_general}` and :math:`\eqref{eq_time_inversion_on_massless_general}` that

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal)^2 \Psi_{p, \sigma} = (-1)^{2 s} \Psi_{p, \sigma}
	\end{equation*}

where :math:`s \in \tfrac{1}{2} \Zbb` equals the spin :math:`\jfrak` in the massive case and the absolute helicity :math:`|\sigma|` in the massless case.

Hence in a non-interacting system consisting of an odd number of half-integer spin/helicity particles and any number of integer spin/helicity particles, we have

.. math::
	:nowrap:

	\begin{equation}
		U(\Tcal)^2 \Psi = -\Psi. \label{eq_time_inversion_squared_reverse_sign}
	\end{equation}

Now for any eigenstate :math:`\Psi` of the Hamiltonian, there is an accompanying eigenstate :math:`U(\Tcal) \Psi` since :math:`U(\Tcal)` commutes with the Hamiltonian. The key observation then is that they are necessarily different states! To see this, let's suppose otherwise that :math:`U(\Tcal) \Psi = \zeta \Psi` represent the same state, where :math:`\zeta` is a phase. Then


.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal)^2 \Psi = U(\Tcal) \zeta \Psi = \zeta^{\ast} U(\Tcal) \Psi = |\zeta|^2 \Psi = \Psi
	\end{equation*}

contradicts :math:`\eqref{eq_time_inversion_squared_reverse_sign}`.

As a conclusion, we see that for such systems, any energy eigenvalue has at least a two-fold degeneracy. This is known as `Kramers' degeneracy <https://en.wikipedia.org/wiki/Kramers%27_theorem>`_.


.. _sec_scattering_theory:

Scattering Theory
-----------------

Physics would have been rather boring if nothing interacts, like the free particles that we have been studying so far. On the flip side, physics would have been impossible if we try to know exactly what happens in the interactions. The middle ground, where we assume that the particles are non-interacting long before and after the interaction, and something mysterious happened in between, is called scattering theory --  a place where theories meet experiments.

Non-interacting many-particles state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We shall, as always, start from the easiest part of the theory, which is clearly the non-interacting parts. Recall our grand formulae for the Lorentz transformation on one-particle state :math:`\eqref{eq_translation_formula_for_particle_state}` and :math:`\eqref{eq_lorentz_transformation_formula_for_particle_state}`. For a non-interacting many-particles system, it's conceivable to assume that the Lorentz transformation law is simply a direct product of the individual particles as follows (recall that :math:`U(\Lambda, a) = U(1, a) U(\Lambda, 0)`)

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda, a) \Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} =&~ \exp(-\ifrak a^{\mu} ((\Lambda p_1)_{\mu} + (\Lambda p_2)_{\mu} + \cdots)) \
		\label{eq_lorentz_transformation_formula_for_many_free_particles} \\
		&\times \sqrt{\frac{(\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p_1)_0 (p_2)_0 \cdots}} \nonumber \\
		&\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma'_1 \sigma_1}(W_1(\Lambda, p_1)) D_{\sigma'_2 \sigma_2}(W_2(\Lambda, p_2)) \cdots \nonumber \\
		&\times \Psi_{\Lambda p_1, \sigma'_1, n_1; ~\Lambda p_2, \sigma'_2, n_2; ~\cdots} \nonumber
	\end{align}

where the first component is the translation transformation :math:`\eqref{eq_translation_formula_for_particle_state}`, the second component is the normalization factor, and the third component is the little group representation, and the :math:`\sigma`'s are either the spin :math:`z`-component for massive particles or the helicity for massless particles, and the :math:`n`'s are additional (discrete) labels such as mass, charge, spin, etc.

Notice that by writing a many-particles state as :math:`\Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}`, we have given the particles an order, which is by no means unique. Hence the normalization of these states must take permutations into account as follows

.. math::
	:nowrap:

	\begin{equation}
		\left( \Psi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}, \Psi_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots} \right) = \delta^3(\pbf_1 - \pbf'_1) \delta_{\sigma_1 \sigma'_1} \delta_{n_1 n'_1} \delta^3(\pbf_2 - \pbf'_2) \delta_{\sigma_2 \sigma'_2} \delta_{n_2 n'_2} \pm \text{permutations}
		\label{eq_many_particles_state_normalization_rough}
	\end{equation}

The sign in front of the permutations has to do with the species of the particles, which will be discussed later. Note that although there are many terms in :math:`\eqref{eq_many_particles_state_normalization_rough}`, there is at most one nonzero term, which happens exactly when the two states differ by a permutation.

To suppress the annoyingly many sub-indexes in the states, we shall use letters such as :math:`\alpha, \beta, \cdots` to denote the compound index such as :math:`(p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots)`, so that, for example, :math:`\eqref{eq_many_particles_state_normalization_rough}` can be simplified as

.. math::
	:nowrap:

	\begin{equation*}
		\left( \Psi_{\alpha}, \Psi_{\alpha'} \right) = \delta(\alpha - \alpha')
	\end{equation*}

where the integral volume element reads

.. math::
	:nowrap:

	\begin{equation*}
		\int d\alpha \cdots = \sum_{\sigma_1, n_1; ~\sigma_2, n_2; ~\cdots} \int d^3 \pbf_1 d^3 \pbf_2 \cdots
	\end{equation*}

We have postulated that the transformation law :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}` works for non-interacting particles, but in fact, it's also only possible for non-interacting particles. One way to see this is through an energy calculation by letting :math:`\Lambda = 1` and :math:`a = (\tau, 0, 0, 0)` in :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}` to see that

.. math::
	:nowrap:

	\begin{equation*}
		\exp(\ifrak \tau E_{\alpha}) \Psi_{\alpha} = \exp(\ifrak \tau H) \Psi_{\alpha} = \exp(\ifrak \tau (E_1 + E_2 + \cdots)) \Psi_{\alpha} \implies E_\alpha = E_1 + E_2 + \cdots
	\end{equation*}

where :math:`E_i \coloneqq (p_i)_0` is the energy of the :math:`i`-th particle. There is obviously no energy left for any interaction.

In- and out-states
^^^^^^^^^^^^^^^^^^

As mentioned earlier, scattering theory is concerned with a scenario where interactions happen within a finite time period, long before and after which the system can be regarded as non-interacting. We can therefore define the in-state :math:`\Psi_{\alpha}^-` and the out-state :math:`\Psi_{\alpha}^+`, where :math:`\alpha` is the compound index as defined in the previous section, such that the states appear to be non-interacting with the prescribed particle states when *observed* at :math:`t \to \mp \infty`, respectively. [#in_out_state_sign_convention]_

Now it's time to bring forward an implicit assumption on the quantum states that we've been studying so far: they're defined in one chosen inertial frame. Indeed, the Lorentz transformation law :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}` tells us exactly how to transform the state to any other frame. States of this sort are called `Heisenberg picture <https://en.wikipedia.org/wiki/Heisenberg_picture>`_ states: they contain the entire history/future of the system and are not dynamical in time as opposed to the so-called `Schrödinger picture <https://en.wikipedia.org/wiki/Schr%C3%B6dinger_picture>`_ states.

Back to the scattering scenario, let's imagine a reference observer :math:`\Ocal`, who at :math:`t = 0` observes that the system is in a state :math:`\Psi`. Then imagine another observer :math:`\Ocal'` at rest with respect to :math:`\Ocal`, who sets his clock :math:`t' = 0` when :math:`t = \tau`, in other words :math:`t' = t - \tau`. Then from the viewpoint of :math:`\Ocal'`, the time-:math:`0` state should look like :math:`\exp(-\ifrak \tau H) \Psi`. It follows that the state :math:`\Psi`, viewed long before and long after the reference :math:`t = 0`, should look like :math:`\exp(-\ifrak \tau H) \Psi` for :math:`\tau \to \mp\infty`, respectively.

It follows that energy eigenstates such as :math:`\Psi_{\alpha}` will look the same at all time because :math:`\exp(-\ifrak \tau H) \Psi_{\alpha} = \exp(-\ifrak \tau E_{\alpha}) \Psi_{\alpha}` creates merely an inconsequential phase factor. This is one form of the uncertainty principle: if the energy is definitely known, then the time is completely unknown. Therefore we must consider a localized packet (or superposition) of states as follows

.. math::
	:nowrap:

	\begin{equation}
		\int d\alpha ~g(\alpha) \Psi_{\alpha}
		\label{eq_psi_packet}
	\end{equation}

where :math:`g(\alpha)` is a reasonably smooth function (e.g. without poles) which is non-vanishing within a finite range of energies. We can then demand that the time limits

.. math::
	:nowrap:

	\begin{equation*}
		\exp(-\ifrak \tau H) \int d\alpha ~g(\alpha) \Psi_{\alpha}^{\pm} = \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm}
	\end{equation*}

as :math:`\tau \to \pm\infty`, respectively, approach the corresponding superpositions of non-interacting particle states.

To be more precise, let's split the Hamiltonian into the free part and the interaction part as follows

.. math::
	:nowrap:

	\begin{equation}
		H = H_0 + V  \label{eq_h_as_h0_plus_v}
	\end{equation}

such that the energy eigenstates :math:`\Phi_{\alpha}` of :math:`H_0` (in the same frame as :math:`\Psi_{\alpha}^{\pm}`) transform according to :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`. Then the asymptotic freeness translates into the following conditions

.. math::
	:nowrap:

	\begin{equation}
		\lim_{\tau \to \pm\infty} \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm} = \lim_{\tau \to \pm\infty} \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha}
		\label{eq_in_out_states_asymptotic_by_energy}
	\end{equation}

or equivalently in terms of the Hamiltonians

.. math::
	:nowrap:

	\begin{equation}
		\lim_{\tau \to \pm\infty} \exp(-\ifrak \tau H) \int d\alpha ~g(\alpha) \Psi_{\alpha}^{\pm} = \lim_{\tau \to \pm\infty} \exp(-\ifrak \tau H_0) \int d\alpha ~g(\alpha) \Phi_{\alpha}
		\label{eq_in_out_states_asymptotic_by_hamiltonian}
	\end{equation}

This motivates the following definition

.. math::
	:nowrap:

	\begin{equation}
		\Omega(\tau) \coloneqq \exp(\ifrak \tau H) \exp(-\ifrak \tau H_0)
		\label{eq_defn_of_Omega}
	\end{equation}

so that :math:`\Psi_{\alpha}^{\pm} = \Omega(\pm\infty) \Phi_{\alpha}`, at least formally. Moreover, since :math:`\Omega` is unitary, the in- and out-states :math:`\Psi_{\alpha}^{\pm}` are normalized as long as :math:`\Phi_{\alpha}` are normalized.

In practice it will be assumed that the interaction term :math:`V` in :math:`\eqref{eq_h_as_h0_plus_v}` is relatively small so that a formal solution as power series in :math:`V` may be meaningful. As the first step, let's try to apply :math:`\eqref{eq_h_as_h0_plus_v}` to :math:`\Psi_{\alpha}^{\pm}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		E_{\alpha} \Psi_{\alpha}^{\pm} = H \Psi_{\alpha}^{\pm} = (H_0 + V) \Psi_{\alpha}^{\pm} \implies (E_{\alpha} - H_0) \Psi_{\alpha}^{\pm} = V \Psi_{\alpha}^{\pm}
	\end{equation*}

Note that :math:`\Phi_{\alpha}` is also annihilated by :math:`E_{\alpha} - H_0`. Considering the asymptotic :math:`\eqref{eq_in_out_states_asymptotic_by_energy}` or :math:`\eqref{eq_in_out_states_asymptotic_by_hamiltonian}`, it's reasonable to guess the following formal solution

.. math::
	:nowrap:

	\begin{equation}
		\Psi_{\alpha}^{\pm} = \Phi_{\alpha} + (E_{\alpha} - H_0 \mp \ifrak \epsilon)^{-1} V \Psi_{\alpha}^{\pm}
		\label{eq_lippmann_schwinger_mixed}
	\end{equation}

where the infinitesimal :math:`\mp \ifrak \epsilon` is a mathematical trick added to avoid division by zero, and the signs will be justified momentarily. One can obviously apply :math:`\eqref{eq_lippmann_schwinger_mixed}` recursively to get an expansion of :math:`\Psi_{\alpha}^{\pm}` as a power series in :math:`V`, and we shall come back to this point later. In order to express :math:`\Psi_{\alpha}^{\pm}` in terms of :math:`\Phi_{\alpha}`, let's expand the right-hand-side of :math:`\eqref{eq_lippmann_schwinger_mixed}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\Psi_{\alpha}^{\pm} = \Phi_{\alpha} + \int d\beta ~\frac{(\Phi_{\beta}, V \Psi_{\alpha}^{\pm}) \Phi_{\beta}}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon}
		\label{eq_lippmann_schwinger_pure}
	\end{equation}

Both :math:`\eqref{eq_lippmann_schwinger_mixed}` and :math:`\eqref{eq_lippmann_schwinger_pure}` are known as the `Lippmann-Schwinger equation <https://en.wikipedia.org/wiki/Lippmann%E2%80%93Schwinger_equation>`_.

Now let's justify the term :math:`\pm \ifrak \epsilon` by showing that :math:`\eqref{eq_lippmann_schwinger_pure}` indeed satisfies the asymptotic condition :math:`\eqref{eq_in_out_states_asymptotic_by_energy}` as follows

.. math::
	:nowrap:

	\begin{align}
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^{\pm} \
			&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha} \
	 			+ \int d\alpha d\beta ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm}) \Phi_{\beta}}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon} \label{eq_packet_expansion_by_lippmann_schwinger} \\
			&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha} \
   				+ \int d\beta ~\Phi_{\beta} \blue{\int d\alpha ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm})}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon}} \nonumber
 	\end{align}

Now the integral colored in blue can be integrated over :math:`E_{\alpha}` by a contour that runs from :math:`-\infty` to :math:`+\infty`, followed by a semicircle at infinity, in the upper-half-plane in the case of :math:`\Psi_{\alpha}^-` and the lower-half-plane in the case of :math:`\Psi_{\alpha}^+`, back to :math:`-\infty`. In either case, the sign in :math:`\mp \ifrak \epsilon` is chosen so that the integrant has no poles with infinitesimally small imaginary part, though both :math:`g(\alpha)` and :math:`(\Phi_{\beta}, V \Psi_{\alpha}^{\pm})`, viewed as complex functions, may have poles with finite imaginary parts. It follows then from the residual theorem and the damping factor :math:`\exp(-\ifrak \tau E_{\alpha})` as :math:`\tau \to \pm\infty` that the integral in blue vanishes, as desired.


.. _sec_s_matrix_and_its_symmetry:

S-matrix and its symmetry
^^^^^^^^^^^^^^^^^^^^^^^^^

The `S-matrix <https://en.wikipedia.org/wiki/S-matrix>`_ defined by

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} \coloneqq \left( \Psi_{\beta}^+, \Psi_{\alpha}^- \right)
		\label{eq_defn_s_matrix_by_in_and_out_states}
	\end{equation}

records the probability amplitude of finding the out-state :math:`\Psi_{\beta}^+` given the in-state :math:`\Psi_{\alpha}^-`. Note that since the in- and out-states both form an orthonormal basis of the same Hilbert space, the S-matrix is unitary. However, the way :math:`S` is defined in :math:`\eqref{eq_defn_s_matrix_by_in_and_out_states}` disqualifies it as an operator on the Hilbert space. Therefore it'll be convenient to convert both in- and out-states to the free states and define the *S-operator* by

.. math::
	:nowrap:

	\begin{equation}
		(\Phi_{\beta}, S \Phi_{\alpha}) \coloneqq S_{\beta \alpha} \label{eq_defn_s_operator}
	\end{equation}

Using :math:`\eqref{eq_defn_of_Omega}` we see that

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = (\Omega(\infty) \Phi_{\beta}, \Omega(-\infty) \Phi_{\alpha}) = (\Phi_{\beta}, \Omega^{\dagger}(\infty) \Omega(-\infty) \Phi_{\alpha}) \implies S = \Omega^{\dagger}(\infty) \Omega(-\infty) \eqqcolon U(\infty, -\infty)
		\label{eq_s_operator_by_u}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{equation}
		U(\tau_1, \tau_0) = \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)
		\label{eq_defn_u_operator}
	\end{equation}

The most straightforward way to calculate :math:`S_{\beta \alpha}` is probably to use :math:`\eqref{eq_lippmann_schwinger_pure}` directly. However. this turns out to be rather involved, and doesn't lead to a simple result. The issue is that we don't really want to convert both the in- and out-states to the non-interacting states, but rather to push, say, the in-states from the far past to the far future and compare with the out-states. To spell out the details, let's first calculate the asymptotic of the in-packet as :math:`\tau \to \infty` (but omitting the :math:`\lim_{\tau \to \infty}` symbol) using :math:`\eqref{eq_packet_expansion_by_lippmann_schwinger}`

.. math::
	:nowrap:

	\begin{align}
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- &= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta}
			\label{eq_positive_limit_of_in_state_by_lippmann_schwinger} \\
	  		&\phantom{=} + \int d\beta ~\Phi_{\beta} \int d\alpha \frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)}{E_{\alpha} - E_{\beta} + \ifrak \epsilon} \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} \nonumber \\
	 		&\phantom{=} - 2\pi\ifrak \int d\beta ~\Phi_{\beta} \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) \exp(-\ifrak \tau E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \left( g(\beta) - 2\pi\ifrak \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \right) \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \left( \blue{\delta(\alpha - \beta) - 2\pi\ifrak \delta(E_{\alpha} - E_{\beta}) (\Phi_{\beta}, V \Psi_{\alpha}^-)} \right) \nonumber
	\end{align}

where we've used the residue theorem again in the second equality. Next expand the left-hand-side of the equation in terms of the out-states and then let :math:`\tau \to \infty`

.. math::
	:nowrap:

	\begin{align}
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- \
			&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \int d\beta ~(\Psi_{\beta}^+, \Psi_{\alpha}^-) \Psi_{\beta}^+
			\label{eq_positive_limit_of_in_state_by_expanding_out_states} \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Psi_{\beta}^+ \int d\alpha ~g(\alpha) S_{\beta \alpha} \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \blue{S_{\beta \alpha}} \nonumber
	\end{align}

where we've used the fact that the S-matrix contains a :math:`\delta(E_{\alpha} - E_{\beta})` factor by energy conservation in the second equality, and the defining property :math:`\eqref{eq_in_out_states_asymptotic_by_energy}` of the out-state in the third equality.

Equating the blue terms from :math:`\eqref{eq_positive_limit_of_in_state_by_lippmann_schwinger}` and :math:`\eqref{eq_positive_limit_of_in_state_by_expanding_out_states}`, we've derived the following formula

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = \delta(\beta - \alpha) - 2\pi\ifrak \delta(E_{\beta} - E_{\alpha}) (\Phi_\beta, V \Psi_{\alpha}^-)
		\label{eq_s_matrix_pre_born_approx}
	\end{equation}

Up to the first order in :math:`V`, one can replace :math:`\Psi_{\alpha}^-` on the right-hand-side by :math:`\Phi_{\alpha}` and arrive at the so-called `Born approximation <https://en.wikipedia.org/wiki/Born_approximation>`_ of the S-matrix.

Lorentz symmetry
++++++++++++++++

Recall that in :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`, or really in :ref:`Lorentz symmetry of one-particle states <sec_lorentz_symmetry>`, we understood how Lorentz transformations act on particle states. Now we'd like to understand how they act on the S-matrix. Of course, since :math:`U(\Lambda, a)` is unitary, we always have

.. math::
	:nowrap:

	\begin{equation*}
		S_{\beta \alpha} = (\Psi_{\beta}^+, \Psi_{\alpha}^-) = \left( U(\Lambda, a) \Psi_{\beta}^+, U(\Lambda, a) \Psi_{\alpha}^- \right)
	\end{equation*}

but this is *not* what we mean by Lorentz symmetry. What we do want to know is, just like in :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`, how Lorentz transformation acts on the particle states, i.e., the (compound) indexes :math:`\alpha` and :math:`\beta`. Now although :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}` doesn't work for general (interacting) states, it does work for, say, :math:`\Psi_{\alpha}^-` in the :math:`\tau \to -\infty` limit because of the asymptotic freeness. By Lorentz we mean that :math:`U(\Lambda, a)` acts the same way on both in- and out-states. In other words, we'll be looking for some :math:`U(\Lambda, a)` such that the following general formula holds.

.. math::
	:nowrap:

	\begin{align}
		S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \
			&\blue{\exp\left( \ifrak a^{\mu} \Lambda_{\mu}^{\nu} \left( (p'_1)_{\nu} + (p'_2)_{\nu} + \cdots - (p_1)_{\nu} - (p_2)_{\nu} - \cdots \right) \right)}  \label{eq_lorentz_transformation_formula_for_s_matrix} \\
			&\times \sqrt{\frac{(\Lambda p'_1)_0 (\Lambda p'_2)_0 \cdots (\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p'_1)_0 (p'_2)_0 \cdots (p_1)_0 (p_2)_0 \cdots}}  \nonumber \\
			&\times \sum_{\bar{\sigma}'_1 \bar{\sigma}'_2 \cdots} D^{\ast}_{\bar{\sigma}'_1 \sigma'_1} (W(\Lambda, p'_1)) D^{\ast}_{\bar{\sigma}'_2 \sigma'_2} (W(\Lambda, p'_2)) \cdots  \nonumber \\
			&\times \sum_{\bar{\sigma}_1 \bar{\sigma}_2 \cdots} D_{\bar{\sigma}_1 \sigma_1} (W(\Lambda, p_1)) D_{\bar{\sigma}_2 \sigma_2} (W(\Lambda, p_2)) \cdots  \nonumber \\
			&\times S_{\Lambda p'_1, \bar{\sigma}'_1, n'_1; ~\Lambda p'_2, \bar{\sigma}'_2, n'_2; ~\cdots, ~~\Lambda p_1, \bar{\sigma}_1, n_1; ~\Lambda p_2, \bar{\sigma}_2, n_2, ~\cdots}  \nonumber
	\end{align}

where we've used primes to distinguish between labels from in- and out-states, and bars to distinguish between labels, specifically the spin-:math:`z` or helicity, before and after the Lorentz transformation.

Since the left-hand-side doesn't depend on the translation parameter :math:`a`, the blue term on the right-hand-side must be :math:`1`. In other words,

.. math::
	:nowrap:

	\begin{equation*}
		p_1 + p_2 + \cdots = p'_1 + p'_2 + \cdots
	\end{equation*}

which is nothing but the conservation of (total) momentum. Note that a special case, which is the energy conservation, has already been used in the derivation of :math:`\eqref{eq_positive_limit_of_in_state_by_expanding_out_states}` from the previous section.

As a consequence, we can now extract a delta function from the S-matrix as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} \eqqcolon \delta(\beta - \alpha) - 2\pi\ifrak M_{\beta \alpha} \delta^4 (p_{\beta} - p_{\alpha})
		\label{eq_s_matrix_with_m}
	\end{equation}

which should be compared with :math:`\eqref{eq_s_matrix_pre_born_approx}`.

Back to the core question of this section, how in the world can one engineer a magic :math:`U(\Lambda, a)` to satisfy the monstrous :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}`? One cannot. But remember that :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}` is readily satisfied for non-interacting particles. It follows that if we consider instead the S-operator defined by :math:`\eqref{eq_defn_s_operator}`, and let :math:`U_0(\Lambda, a)` be the Lorentz transformation on free particles defined by :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`, then :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}` would be satisfied if :math:`U_0(\Lambda, a)` commutes with :math:`S`. Indeed, using shorthand notations, we have

.. math::
	:nowrap:

	\begin{equation*}
		S_{\beta \alpha} = \left( \Phi_{\beta}, S \Phi_{\alpha} \right) \
			= \left( U_0 \Phi_{\beta}, U_0 S \Phi_{\alpha} \right) \
			= \left( U_0 \Phi_{\beta}, S U_0 \Phi_{\alpha} \right) \
			= \cdots \times S_{\underline{\beta \alpha}}
	\end{equation*}

where :math:`\cdots` denotes the coefficients on the right-hand-side of :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}` and the sub-indexes :math:`\underline{\beta \alpha}` is short for the compound index on the right-hand-side of :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}`.

In order for :math:`S` to commute with :math:`U_0(\Lambda, a)`, it suffices that it commutes with the infinitesimal generators of :math:`U_0(\Lambda, a)`, namely,

.. math::
	:nowrap:

	\begin{alignat}{2}
		&[H_0, S] &&= 0    \label{eq_h0_s_commute} \\
		&[\Pbf_0, S] &&= 0 \label{eq_p30_s_commute} \\
		&[\Jbf_0, S] &&= 0 \label{eq_j30_s_commute} \\
		&[\Kbf_0, S] &&= 0 \label{eq_k30_s_commute}
	\end{alignat}

where :math:`H_0, \Pbf_0, \Jbf_0, \Kbf_0` are discussed in :ref:`sec_quantum_lorentz_symmetry` and satisfy the commutation relations :math:`\eqref{eq_hp_commute}` -- :math:`\eqref{eq_kkj_commutation}`.

This shall be done in three steps, where :math:`\eqref{eq_p30_s_commute}, \eqref{eq_j30_s_commute}` will be handled first, followed by :math:`\eqref{eq_k30_s_commute}`, and finally :math:`\eqref{eq_h0_s_commute}`.

Step 1.
	Recall from :math:`\eqref{eq_s_operator_by_u}` and :math:`\eqref{eq_defn_u_operator}` that the S-operator can be understood as time translations dictated by :math:`H` and :math:`H_0`.  It's therefore necessary to understand how the free infinitesimal Lorentz transformations commute with :math:`H`. To this end, let's consider the in-states at :math:`\tau \to -\infty`, which is approximately free. We can similarly define infinitesimal operators :math:`\Pbf, \Jbf, \Kbf` that together with :math:`H` satisfy the same commutation relations :math:`\eqref{eq_hp_commute}` -- :math:`\eqref{eq_kkj_commutation}`.

	Now comes the crucial part, which is to make assumptions about :math:`H` so that :math:`\eqref{eq_h0_s_commute}` -- :math:`\eqref{eq_k30_s_commute}` are satisfied. Recall that :math:`H = H_0 + V` where :math:`V` describes the interactions. The first assumption we'll make is the following

	.. admonition:: Assumption on :math:`H` for Lorentz invariance of S-matrix #1
		:class: Important

		The interaction :math:`V` affects neither the momentum :math:`\Pbf` nor the angular momentum :math:`\Jbf`. In other words, we assume that

		.. math::
			:nowrap:

			\begin{equation}
				\Pbf = \Pbf_0, ~~\Jbf = \Jbf_0, ~~[V, \Pbf_0] = [V, \Jbf_0] = 0
				\label{eq_s_matrix_lorentz_invariance_assump_1}
			\end{equation}

	It follows from this assumption that :math:`\eqref{eq_p30_s_commute}` and :math:`\eqref{eq_j30_s_commute}` hold.

Step 2.
	Next we turn to :math:`\eqref{eq_k30_s_commute}`. This time we cannot "cheat" by assuming that :math:`\Kbf = \Kbf_0` because it would led to the undesirable consequence :math:`H = H_0` by :math:`\eqref{eq_pkh_commutation}`. So instead, let's write

	.. math::
		:nowrap:

		\begin{equation}
			\Kbf = \Kbf_0 + \Wbf
			\label{eq_k_as_k0_plus_w}
		\end{equation}

	where :math:`\Wbf` denotes the perturbation term. Let's calculate

	.. math::
		:nowrap:

		\begin{equation}
			[\Kbf_0, S] = \lim_{\substack{\tau_0 \to -\infty \\ \tau_1 \to \infty\phantom{-}}} [\Kbf_0, U(\tau_1, \tau_0)] \
				= \lim_{\substack{\tau_0 \to -\infty \\ \tau_1 \to \infty\phantom{-}}} [\Kbf_0, \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)]
			\label{eq_k0_s_commutator_by_u}
		\end{equation}

	as follow. First using :math:`\eqref{eq_hkp_commutation}` again

	.. math::
		:nowrap:

		\begin{alignat*}{2}
			[\Kbf_0, \exp(\ifrak \tau H_0)] &= [\Kbf_0, \ifrak \tau H_0] \exp(\ifrak \tau H_0) &&= \tau \Pbf_0 \exp(\ifrak \tau H_0) \\
			[\Kbf, \exp(\ifrak \tau H)] &= [\Kbf, \ifrak \tau H] \exp(\ifrak \tau H) &&= \tau \Pbf \exp(\ifrak \tau H) =\tau \Pbf_0 \exp(\ifrak \tau H)
		\end{alignat*}

	from which we can calculate

	.. math::
		:nowrap:

		\begin{align}
			[\Kbf_0, U(\tau_1, \tau_0)] &= [\Kbf_0, \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)] \label{eq_k30_u_commutation} \\
				&= [\Kbf_0, \exp(\ifrak \tau_1 H_0)] \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0) \nonumber \\
				&\phantom{=} + \exp(\ifrak \tau_1 H_0) [\Kbf - \Wbf, \exp(\ifrak (\tau_0 - \tau_1) H)] \exp(-\ifrak \tau_0 H_0) \nonumber \\
				&\phantom{=} + \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) [\Kbf_0, \exp(-\ifrak \tau_0 H_0)] \nonumber \\
				&= \blue{\tau_1 \Pbf_0 \exp(\ifrak \tau H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \nonumber \\
				&\phantom{=} \blue{+ (\tau_0 - \tau_1) \Pbf_0 \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \nonumber \\
				&\phantom{=} \blue{- \tau_0 \Pbf_0 \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)} \nonumber \\
				&\phantom{=} - \exp(\ifrak \tau_1 H_0) [\Wbf, \exp(\ifrak (\tau_0 - \tau_1) H)] \exp(-\ifrak \tau_0 H_0) \nonumber \\
				&= -\Wbf(\tau_1) U(\tau_1, \tau_0) + U(\tau_1, \tau_0) \Wbf(\tau_0) \nonumber
		\end{align}

	where :math:`\Wbf(\tau) \coloneqq \exp(\ifrak \tau H_0) \Wbf \exp(-\ifrak \tau H_0)`. Note that the three blue terms cancel out.

	We see that :math:`\eqref{eq_k0_s_commutator_by_u}`, and hence :math:`\eqref{eq_k30_s_commute}`, would follow if :math:`W(\tau) \to 0` as :math:`\tau \to \pm\infty`. The latter, in turn, would follow from the following assumption

	.. admonition:: Assumption on :math:`H` for Lorentz invariance of S-matrix #2
		:class: Important

		The matrix elements of :math:`W` with respect to the eigenstates :math:`\Phi_{\alpha}` of :math:`H_0` is smooth, so that :math:`W(\tau)` vanishes on any local packet of :math:`\Phi_{\alpha}` as in :math:`\eqref{eq_psi_packet}` as :math:`\tau \to \pm\infty`.

	This assumption should be compared with :math:`\eqref{eq_in_out_states_asymptotic_by_energy}` and :math:`\eqref{eq_in_out_states_asymptotic_by_hamiltonian}`, and can be justified by the asymptotic freeness of S-matrix theory.

Step 3.
	Finally let's handle :math:`\eqref{eq_h0_s_commute}`. Recall from :math:`\eqref{eq_s_operator_by_u}` that :math:`S = \Omega^{\dagger}(\infty) \Omega(-\infty)`. Hence the idea is to work out how :math:`H` and :math:`H_0` intertwine with :math:`\Omega(\pm\infty)`. To this end, let's use :math:`\eqref{eq_k30_u_commutation}` by setting :math:`\tau_1 = 0` and :math:`\tau_0 = \mp\infty` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			[\Kbf_0, \Omega(-\infty)] = -\Wbf \Omega(-\infty) \implies \Kbf \Omega(\mp\infty) = \Omega(\mp\infty) \Kbf_0
		\end{equation*}

	Moreover, by :math:`\eqref{eq_s_matrix_lorentz_invariance_assump_1}`, we have also :math:`\Pbf \Omega(\mp\infty) = \Omega(\mp\infty) \Pbf_0`. Now using the commutation relation :math:`\eqref{eq_pkh_commutation}` we conclude that

	.. math::
		:nowrap:

		\begin{equation*}
			H \Omega(\mp\infty) = \Omega(\mp\infty) H_0
		\end{equation*}

	which readily implies :math:`\eqref{eq_h0_s_commute}`.

.. note::
	Besides showing that :math:`\eqref{eq_h0_s_commute}` -- :math:`\eqref{eq_k30_s_commute}` hold, our calculations actually establish the following intertwining identities

	.. math::
		:nowrap:

		\begin{align*}
			H \Omega(\pm\infty) &= \Omega(\pm\infty) H_0 \\
			\Pbf \Omega(\pm\infty) &= \Omega(\pm\infty) \Pbf_0 \\
			\Jbf \Omega(\pm\infty) &= \Omega(\pm\infty) \Jbf_0
		\end{align*}

	which imply, in particular, that the standard commutation relations :math:`\eqref{eq_hp_commute}` -- :math:`\eqref{eq_kkj_commutation}` also hold in a frame where :math:`\tau \to \infty`, as expected.

Internal symmetry
+++++++++++++++++

An internal symmetry is a symmetry that leaves :math:`p` and :math:`\sigma` invariant and acts on the other labels such as charge, spin, and so on. We can write the general form of an internal symmetry on in- and out-states as follows

.. math::
	:nowrap:

	\begin{equation}
		U(T) \Psi^{\pm}_{p_1, \sigma_1, n_1;~p_2, \sigma_2, n_2;~\cdots} = \sum_{n'_1, n'_2, \cdots} \Dscr_{n'_1 n_1}(T) \Dscr_{n'_2 n_2}(T) \cdots \Psi^{\pm}_{p_1, \sigma_1, n'_1;~p_2, \sigma_2, n'_2;~\cdots}
		\label{eq_internal_symmetry_transformation_for_in_and_out_states}
	\end{equation}

where :math:`U(T)` is the unitary operator associated with the symmetry transformation :math:`T`, and the :math:`\Dscr`'s are analogs of the little group representations from :math:`\eqref{eq_d_repr_of_little_group}`.

Similar to :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}`, we can formulate the internal symmetry of S-matrix as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{n'_1, n'_2, \cdots, ~n_1, n_2, \cdots} = \
			\sum_{\bar{n}'_1, \bar{n}'_2, \cdots} \Dscr^{\ast}_{\bar{n}'_1 n'_1}(T) \Dscr^{\ast}_{\bar{n}'_2 n'_2}(T) \cdots \
			\sum_{\bar{n}_1, \bar{n}_2, \cdots} \Dscr_{\bar{n}_1 n_1}(T) \Dscr_{\bar{n}_2 n_2}(T) \cdots \
			S_{\bar{n}'_1, \bar{n}'_2, \cdots, ~\bar{n}_1, \bar{n}_2, \cdots}
		\label{eq_internal_symmetry_transformation_formula_for_s_matrix}
	\end{equation}

where we have suppressed the irrelevant :math:`p` and :math:`\sigma` labels.

For what kind of Hamiltonian :math:`H` does there exist an internal symmetry :math:`U(T)` that acts like :math:`\eqref{eq_internal_symmetry_transformation_for_in_and_out_states}`? The answer is similar to the case of Lorentz symmetry. Namely, if we can split :math:`H = H_0 + V` into the free and perturbation terms, such that the free symmetry transformation :math:`U_0(T)`, which satisfies :math:`\eqref{eq_internal_symmetry_transformation_for_in_and_out_states}` with :math:`\Phi` in place of :math:`\Psi^{\pm}`, commutes with both :math:`H_0` and :math:`V`.


Similar to the translations in Lorentz symmetry, let's consider a symmetry :math:`T(\theta)` parametrized by a real number. It follows from :math:`\eqref{eq_additive_symmetry}` that we can write

.. math::
	:nowrap:

	\begin{equation}
		U(T(\theta)) = \exp(\ifrak \theta Q)
		\label{eq_charge_internal_symmetry}
	\end{equation}

where :math:`Q` is a Hermitian operator called the charge. Probably the best known example of it is the electric charge. In this case, we can also write

.. math::
	:nowrap:

	\begin{equation}
		\Dscr_{n n'}(T(\theta)) = \delta_{n n'} \exp(\ifrak \theta q_n)
		\label{eq_infinitesimal_charge_d_matrix}
	\end{equation}

The general formula :math:`\eqref{eq_internal_symmetry_transformation_formula_for_s_matrix}` then translates into

.. math::
	:nowrap:

	\begin{equation*}
		q_1 + q_2 + \cdots = q'_1 + q'_2 + \cdots
	\end{equation*}

which is nothing about the conservation of charges. Besides the electric charge, there exist also other similar conserved, or approximately conserved, quantities, such as baryon number and lepton number.


.. _sec_parity_symmetry:

Parity symmetry
+++++++++++++++

Recall from :ref:`sec_space_inversion_for_massive_particles` that for non-interacting massive particles

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} \
			= \eta_{n_1} \eta_{n_2} \cdots \Psi^{\pm}_{U(\Pcal)p_1, \sigma_1, n_1; ~U(\Pcal)p_2, \sigma_2, n_2; ~\cdots}
		\label{eq_space_inversion_acts_on_in_and_out_states}
	\end{equation}

where :math:`\eta_n` denotes the intrinsic parity of particle :math:`n`. The S-matrix version of the parity symmetry is as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} \
			= \eta^{\ast}_{n'_1} \eta^{\ast}_{n'_2} \cdots \eta_{n_1} \eta_{n_2} \cdots \
			S_{U(\Pcal)p'_1, \sigma'_1, n'_1; ~U(\Pcal)p'_2, \sigma'_2, n'_2; ~\cdots, ~~U(\Pcal)p_1, \sigma_1, n_1; ~U(\Pcal)p_2, \sigma_2, n_2; ~\cdots}
		\label{eq_space_inversion_formula_for_s_matrix}
	\end{equation}

Although the space inversion operator :math:`\Pcal` is defined explicitly in :math:`\eqref{eq_space_inversion}`, the parity operator :math:`U(\Pcal)`, as far as the S-matrix is concerned, is completely determined by :math:`\eqref{eq_space_inversion_acts_on_in_and_out_states}` and :math:`\eqref{eq_space_inversion_formula_for_s_matrix}`. In particular, it's not uniquely determined if the particle species under question possesses internal symmetries as discussed in the previous section, because their composition with :math:`\Pcal` will also satisfy :math:`\eqref{eq_space_inversion_acts_on_in_and_out_states}` and :math:`\eqref{eq_space_inversion_formula_for_s_matrix}`, and therefore may equally well be called a parity operator.

Since :math:`\Pcal^2 = 1`, it's an obvious question to ask whether :math:`U(\Pcal)^2 = 1` necessarily. This would have been the case if :math:`U` furnishes a genuine representation, but it doesn't have to. In general, we have

.. math::
	:nowrap:

	\begin{equation*}
		U(\Pcal)^2 \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \eta_{n_1}^2 \eta_{n_2}^2 \cdots \
			\Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots}
	\end{equation*}

which looks just like an internal symmetry. Now if :math:`U(\Pcal)^2` belongs to a continuous family of internal symmetries, then it may be redefined, by suitably composing with internal symmetries, so that all :math:`\eta^2 = 1`. Examples of this kind include notably protons and neutrons. On the other hand, non-examples, i.e., those whose intrinsic parity cannot be reduced to :math:`\pm 1`, include only those hypothetical `Majorana fermions <https://en.wikipedia.org/wiki/Majorana_fermion>`_.

.. todo::
	Revise this part after I've learned more...

.. _dropdown_parities_of_elementary_particles:

.. dropdown:: Parities of elementary particles
	:icon: unlock
	:animate: fade-in-slide-down

	We shall in this section assume familiarity with angular momentum as discussed in :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`.

	Can parities be other than :math:`\pm 1`?
		Aside from electric charges, there exist other quantities that are (at least approximately) conserved by internal symmetries, for example, `baryon numbers <https://en.wikipedia.org/wiki/Baryon_number>`_ :math:`B` and `lepton numbers <https://en.wikipedia.org/wiki/Lepton_number>`_ :math:`L`. Examples of baryons include `protons <https://en.wikipedia.org/wiki/Proton>`_ and `neutrons <https://en.wikipedia.org/wiki/Neutron>`_. Examples of leptons include `electrons <https://en.wikipedia.org/wiki/Electron>`_, `muons <https://en.wikipedia.org/wiki/Muon>`_ and `neutrinos <https://en.wikipedia.org/wiki/Neutrino>`_. The internal symmetry operator generalizes :math:`\eqref{eq_charge_internal_symmetry}` in a straightforward way as follows

		.. math::
			:nowrap:

			\begin{equation}
				U(T(\alpha, \beta, \gamma)) = \exp(\ifrak(\alpha B + \beta L + \gamma Q))
				\label{eq_baryon_lepton_charge_internal_symmetry}
			\end{equation}

		so that :math:`T` is isomorphic to :math:`\Rbb^3` instead of :math:`\Rbb`. This will be the most general internal symmetry that will be considered here.

		By the conservation of angular momentum, the parity of the number of half-integer spin particles, which we denote by :math:`(-1)^F`, is conserved. Here :math:`F` stands for fermion. For all known (to Weinberg at least) particles, the following equality of parities holds

		.. math::
			:nowrap:

			\begin{equation}
				(-1)^F = (-1)^{B + L}
				\label{eq_fermion_count_eq_baryon_and_lepton_mod_2}
			\end{equation}

		In particular, the above mentioned protons, neutrons, electrons, neutrinos are all spin-:math:`1/2` particles.

		If, for whatever reason, the following holds

		.. math::
			:nowrap:

			\begin{equation}
				\orange{U(\Pcal)^2 = (-1)^F}
				\label{}
			\end{equation}

		and in addition :math:`\eqref{eq_fermion_count_eq_baryon_and_lepton_mod_2}` holds, then :math:`U(\Pcal)^2` is part of a continuous symmetry :math:`\eqref{eq_baryon_lepton_charge_internal_symmetry}` and hence can be set to one.  A hypothetical example that breaks :math:`\eqref{eq_fermion_count_eq_baryon_and_lepton_mod_2}` is the so-called Majorana fermions that are their own antiparticles, which implies :math:`B = L = 0`. For these particles, we have :math:`U(\Pcal)^4 = 1`, and hence the intrinsic parity may be :math:`\pm 1` or :math:`\pm \ifrak`.

	Can parities be :math:`-1`?
		The following reaction is observed experimentally

		.. math::
			:nowrap:

			\begin{equation}
				\pi^- + d \to n + n
				\label{eq_pion_deuteron_to_two_neutrons}
			\end{equation}

		where a negative pion is absorbed by a `deuteron <https://en.wikipedia.org/wiki/Deuterium>`_ to produce two neutrons. Moreover, the reaction assumes that the initial state, i.e., the left-hand-side of :math:`\eqref{eq_pion_deuteron_to_two_neutrons}` has orbital angular momentum :math:`\ell = 0` and total angular momentum :math:`j = 1`. Note that the spin of pion and deuteron is :math:`0` and :math:`1`, respectively.

		The conservation of angular momentum demands that the total angular momentum of the right-hand-side must also be :math:`1`, and this can be achieved, a priori, in a number of possibilities. Since neutrons have spin :math:`1/2`, the total spin :math:`\sfrak` of :math:`n + n` may be either :math:`0` or :math:`1` by :math:`\eqref{eq_composite_total_angular_momentum_range}`. But since neutrons are fermions and therefore the state :math:`n + n` must be anti-symmetric, we conclude that :math:`\sfrak = 0` by :math:`\eqref{eq_second_highest_weight_am_pair_two}`. [#pion_deuteron_reaction_final_state]_ Then it follows again from :math:`\eqref{eq_composite_total_angular_momentum_range}` that the orbital angular momentum of the right-hand-side of :math:`\eqref{eq_pion_deuteron_to_two_neutrons}` must be :math:`1`. We are left with exactly one choice.

		Now since the orbital angular momentum changes from :math:`0` in the initial state to :math:`1` in the final state, the S-matrix elements flip sign by the action of :math:`U(\Pcal)` (:red:`WHY? I guess I'm missing knowledge about how orbital angular momentum enters the S-matrix.`). It follows from :math:`\eqref{eq_space_inversion_formula_for_s_matrix}` that

		.. math::
			:nowrap:

			\begin{equation*}
				\eta_{\pi^-} \eta_d = -\eta_n^2
			\end{equation*}

		Deuteron is a nucleus consisting of a proton and a neutron. By the previous discussions about internal symmetries, one can arrange so that they have the same intrinsic parity and hence :math:`\eta_d = \eta_n^2`. It follows that :math:`\eta_{\pi^-} = -1` and the pion :math:`\pi^-` is what we set out to look for. Indeed, all its companions :math:`\pi^0` and :math:`\pi^+` also have parity :math:`-1` due to the isospin symmetry.

		The fact that :math:`\eta_{\pi} = -1` had led to a profound consequence because it was discovered through experiments that there are two spin-:math:`0` particles, now known as :math:`K`-mesons, one of which decays into two pions and the other into three pions. By rotational invariance one can exclude the effects of orbital angular momentum and conclude, assuming parity conservation, that they must have opposite intrinsic parities. However, as more experimental evidence pointing towards the fact that the two :math:`K`-mesons look alike, walk alike and quack alike, it was finally suggested by T. D. Lee and C. N. Yang that they're really the same particle and it's the parity conservation that fails to hold in these reactions, now known as the weak interactions. This suggestion was later verified more directly by an experiment of `C. S. Wu <https://en.wikipedia.org/wiki/Wu_experiment>`_.


Time inversion symmetry
+++++++++++++++++++++++

Recall from :ref:`sec_time_inversion_for_massive_particles` that for a single massive particle

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \Psi_{p, \sigma, n} = \zeta (-1)^{\jfrak - \sigma} \Psi_{\Pcal p, -\sigma, n}
	\end{equation*}

To generalize this to the in- and out-states, we need to remember that the time inversion also interchanges the very frame with respect to which the in- and out-states are defined. The result is as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Tcal) \Psi^{\pm}_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \
			\zeta_{n_1} (-1)^{\jfrak_1 - \sigma_1} \zeta_{n_2} (-1)^{\jfrak_2 - \sigma_2} \cdots \
			\Psi^{\mp}_{\Pcal p_1, -\sigma_1, n_1; ~\Pcal p_2, -\sigma_2, n_2; ~\cdots}
		\label{eq_time_inversion_acts_on_in_and_out_states}
	\end{equation}

The invariance of S-matrix can then be formulated as follows

.. math::
	:nowrap:

	\begin{align}
		S_{p'_1, \sigma'_1, n'_1; ~p'_2, \sigma'_2, n'_2; ~\cdots, ~p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} = \
			&\zeta_{n'_1} (-1)^{\jfrak'_1 - \sigma'_1} \zeta_{n'_2} (-1)^{\jfrak'_2 - \sigma'_2} \cdots \
				\zeta^{\ast}_{n_1} (-1)^{\jfrak_1 - \sigma_1} \zeta^{\ast}_{n_2} (-1)^{\jfrak_2 - \sigma_2}
			\label{eq_time_inversion_acts_on_s_matrix}  \\
			&\times S_{\Pcal p_1, -\sigma_1, n_1; ~\Pcal p_2, -\sigma_2, n_2; ~\cdots; ~\Pcal p'_1, -\sigma'_1, n'_1; ~\Pcal p'_2, -\sigma'_2, n'_2; ~\cdots} \nonumber
	\end{align}

Since we'll be mainly concerned with the rate of interactions in this section, the phase factors in front of :math:`\Psi` play little role. So let's simplify the notations in :math:`\eqref{eq_time_inversion_acts_on_in_and_out_states}` and :math:`\eqref{eq_time_inversion_acts_on_s_matrix}` using compound indexes as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Tcal) \Psi^{\pm}_{\alpha} &= \Psi^{\mp}_{\Tcal \alpha}  \nonumber \\
		S_{\beta, \alpha} &= S_{\Tcal\alpha, \Tcal\beta}  \label{eq_time_inversion_formula_for_s_matrix}
	\end{align}

where the phase factors have been "absorbed" in the right-hand-side.

Unlike the space inversions discussed in the previous section, time inversions don't directly lead to implications on reaction rates because, after all, we cannot turn time around in any experiment. However, under certain circumstances, one can use a trick to draw experimentally verifiable conclusions, which we now present.

The main assumption here is that one can expand the S-operator as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = S^{(0)}_{\beta \alpha} + S^{(1)}_{\beta \alpha} + \cdots
		\label{eq_s_operator_first_order_expansion}
	\end{equation}

such that :math:`S^{(1)} \ll S^{(0)}` can be regarded as the first-order perturbation. The unitarity of :math:`S` shows that

.. math::
	:nowrap:

	\begin{equation*}
		1 = S^{\dagger} S = {S^{(0)}}^{\dagger} S^{(0)} + {S^{(0)}}^{\dagger} S^{(1)} + {S^{(1)}}^{\dagger} S^{(0)} + \cdots
	\end{equation*}

which, in turn, implies

.. math::
	:nowrap:

	\begin{alignat*}{2}
		{S^{(0)}}^{\dagger} S^{(0)} &= 1 \quad &&\text{(Unitary)} \\
		S^{(1)} &= -S^{(0)} {S^{(1)}}^{\dagger} S^{(0)} \quad &&\text{((anti-)Hermitian)}
	\end{alignat*}

Using :math:`\eqref{eq_time_inversion_formula_for_s_matrix}`, the (anti-)Hermitian condition can be spelled out in matrix notations as follows

.. math::
	:nowrap:

	\begin{equation}
		S^{(1)}_{\beta \alpha} = -\int d\gamma' \int d\gamma ~S^{(0)}_{\beta \gamma'} ~S^{(1)~\ast}_{\Tcal \gamma' ~\Tcal \gamma} ~S^{(0)}_{\gamma \alpha}
		\label{eq_first_order_s_matrix_is_hermitian}
	\end{equation}

where we recall that the adjoint :math:`\dagger` equals the composition of the (complex) conjugation :math:`\ast` and transpose. Together with the unitarity of :math:`S^{(0)}`, we see that the rate of reaction :math:`\left| S^{(1)}_{\beta \alpha} \right|^2`, when summed up against a complete set of :math:`S^{(0)}` eigenstates, remains the same after applying :math:`\Tcal` to both initial and final states.

The simplest case where :math:`\eqref{eq_first_order_s_matrix_is_hermitian}` becomes applicable is obviously when both :math:`\alpha` and :math:`\beta` are eigenstates of :math:`S^{(0)}`, with eigenvalues, say, :math:`\exp(\ifrak \theta_{\alpha})` and :math:`\exp(\ifrak \theta_{\beta})`, respectively. In this case :math:`\eqref{eq_first_order_s_matrix_is_hermitian}` becomes

.. math::
	:nowrap:

	\begin{equation*}
		S^{(1)}_{\beta \alpha} = -\exp(\ifrak(\theta_{\alpha} + \theta_{\beta})) S^{(1)~\ast}_{\Tcal \beta ~\Tcal \alpha} \
			\implies \left| S^{(1)}_{\beta \alpha} \right|^2 = \left| S^{(1)}_{\Tcal \beta ~\Tcal \alpha} \right|^2
	\end{equation*}

This is to say that under the assumption that :math:`\eqref{eq_s_operator_first_order_expansion}` is valid, at least approximately, the rate of reaction :math:`S^{(1)}_{\beta \alpha}` should be invariant under a flip of the :math:`3`-momentum as well as the spin :math:`z`-component. This is *not* contradicted by Wu's experiment which disproved the parity conservation.


Rates and cross-sections
^^^^^^^^^^^^^^^^^^^^^^^^

As we already mentioned, the S-matrix entries :math:`S_{\beta \alpha}` can be interpreted as probability amplitudes of a reaction that turns an in-state :math:`\Psi^-_{\alpha}` into an out-state :math:`\Psi^+_{\beta}`. In other words, the probability :math:`P(\Psi^-_{\alpha} \to \Psi^+_{\beta}) = \left| S_{\beta \alpha} \right|^2`. It is, however, not completely straightforward to square S-matrix entries because, as we've seen in :math:`\eqref{eq_s_matrix_with_m}`, they contain Dirac delta functions.

Derivation in box model
+++++++++++++++++++++++

One trick that is often used in physics to deal with integration over an infinite space is to restrict the space to a (large) box, often with additional periodic boundary conditions, and hope that the final results will not depend on the size of the box, as long as it's large enough. This is exactly what we shall do.

Consider a cubic box whose sides have length :math:`L` and has volume :math:`V = L^3`. Imposing the periodic boundary condition on the cube, the :math:`3`-momentum is discretized as follows

.. math::
	:nowrap:

	\begin{equation}
		\pbf = \frac{2\pi}{L} (n_1, n_2, n_3)
		\label{eq_momentum_by_wave_number}
	\end{equation}

where :math:`n_1, n_2, n_3` are nonnegative integers. Of course, the higher the :math:`n`, the shorter the wave length if we interpret it as wave mechanics. By analogy with the continuous case, we can define the Dirac delta function as follows

.. math::
	:nowrap:

	\begin{equation}
		\delta^3_V (\pbf - \pbf') \coloneqq \frac{1}{(2\pi)^3} \int_V d^3 \xbf ~\exp(\ifrak (\pbf - \pbf') \cdot \xbf) = \frac{V}{(2\pi)^3} \delta_{\pbf \pbf'}
		\label{eq_3_momentum_delta_in_a_box}
	\end{equation}

where :math:`\delta_{\pbf \pbf'}` is the usual Kronecker delta. With this setup, the states inner product :math:`\eqref{eq_many_particles_state_normalization_rough}` will produce, from the Dirac deltas, an overall factor of :math:`\left( V/(2\pi)^3 \right)^N` where :math:`N` denotes the number of particles in the box. In order for the amplitudes to be independent of the size of the box, let's normalize the states as follows

.. math::
	:nowrap:

	\begin{equation*}
		\Psi_{\alpha}^{\square} \coloneqq \left( \frac{(2\pi)^3}{V} \right)^{N/2} \Psi_{\alpha}
	\end{equation*}

such that :math:`\left( \Psi^{\square}_{\beta}, \Psi^{\square}_{\alpha} \right) = \delta_{\beta \alpha}` is properly normalized. Correspondingly, we can express the S-matrix with respect to the box-normalized states as follows

.. math::
	:nowrap:

	\begin{equation*}
		S^{\square}_{\beta \alpha} = \left( \frac{(2\pi)^3}{V} \right)^{(N_{\alpha} + N_{\beta})/2} S_{\beta \alpha}
 	\end{equation*}

where :math:`N_{\alpha}, N_{\beta}` are the numbers of particles in the in- and out-states, respectively.

Now the transition probability in the box model takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		P(\alpha \to \beta) = \left| S^{\square}_{\beta \alpha} \right|^2 = \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha} + N_{\beta}} \left| S_{\beta \alpha} \right|^2
	\end{equation*}

which we can further turn to a differential form as follows

.. math::
	:nowrap:

	\begin{align}
		dP(\alpha \to \beta) \
			&= P(\alpha \to \beta) d\Nscr_{\beta}  \label{eq_differential_form_of_s_matrix_probability} \\
			&= P(\alpha \to \beta) \left( \frac{V}{(2\pi)^3} \right)^{N_{\beta}} d\beta  \nonumber \\
			&= \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha}} \left| S_{\beta \alpha} \right|^2 d\beta  \nonumber
	\end{align}

where :math:`d\beta` denotes an infinitesimal volume element around the state :math:`\beta`, or more precisely, a product of :math:`d^3 \pbf`, one for each particle. Then :math:`\Nscr_{\beta}` counts the number of states within the infinitesimal :math:`d\beta`, which can be readily calculated from :math:`\eqref{eq_momentum_by_wave_number}`.

Back to our core problem, which is to define :math:`\left| S_{\beta \alpha} \right|^2` as calculated by :math:`\eqref{eq_s_matrix_with_m}`. The first assumption we will make, at least for now, is a genericity condition

.. _assump_genericity_s_matrix:

.. admonition:: Genericity assumption on the S-matrix
	:class: Important

	No subset of particles in the state :math:`\beta` have exactly the same (total) :math:`4`-momentum as some subset in the state :math:`\alpha`.

Under this assumption, we can remove the term :math:`\delta(\beta - \alpha)` from :math:`\eqref{eq_s_matrix_with_m}` and write

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = -2 \pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha}
		\label{eq_generic_s_matrix_with_m}
	\end{equation}

and moreover, ensure that :math:`M_{\beta \alpha}` contains no more delta functions. Now the question becomes how to define :math:`\left| \delta^4(p_{\beta} - p_{\alpha}) \right|^2`. In fact, to align with the main theme of using in- and out-states to calculate the S-matrix, the interaction must be turned on for a finite period of time, say, :math:`T`. Hence the time-wise delta function becomes

.. math::
	:nowrap:

	\begin{equation}
		\delta_T(E_{\beta} - E_{\alpha}) \coloneqq \frac{1}{2 \pi} \int_{-T/2}^{T/2} dt ~\exp(\ifrak (E_{\beta} - E_{\alpha}) t)
		\label{eq_time_delta_in_a_period}
	\end{equation}

We can then modify :math:`\eqref{eq_generic_s_matrix_with_m}` in a "timed box" as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = -2\pi\ifrak \delta^3_V (\pbf_{\beta} - \pbf_{\alpha}) \delta_T(E_{\beta} - E_{\alpha}) M_{\beta \alpha}
		\label{eq_generic_s_matrix_in_time_box}
	\end{equation}

Now using :math:`\eqref{eq_3_momentum_delta_in_a_box}` and :math:`\eqref{eq_time_delta_in_a_period}`, we can calculate the squares as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		\left( \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \right)^2 &= \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \delta^3_V(0) &&= \delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) V/(2\pi)^3 \\
		\left( \delta_T(E_{\beta} - E_{\alpha}) \right)^2 &= \delta_T(E_{\beta} - E_{\alpha}) \delta_T(0) &&= \delta_T(E_{\beta} - E_{\alpha}) T/(2\pi)
	\end{alignat*}

All together, we can now rewrite :math:`\eqref{eq_differential_form_of_s_matrix_probability}` as follows

.. math::
	:nowrap:

	\begin{align*}
		dP(\alpha \to \beta) &= \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha}} \left| S_{\beta \alpha} \right|^2 d\beta \\
			&= (2\pi)^2 \left( \frac{(2\pi)^3}{V} \right)^{N_{\alpha} - 1} \frac{T}{2\pi} \
				\delta^3_V(\pbf_{\beta} - \pbf_{\alpha}) \delta_T(E_{\beta} - E_{\alpha}) \left| M_{\beta \alpha} \right|^2 d\beta \\
			&= (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} T \delta^4(p_{\beta} - p_{\alpha}) \left| M_{\beta \alpha} \right|^2 d\beta
	\end{align*}

where we have restored :math:`\delta^4(p_{\beta} - p_{\alpha})` by taking the large :math:`V` and :math:`T` limits.

If taking partial limits in :math:`V` and :math:`T` in the above derivation is not suspicious enough, then let's define the rate of transition by moving the :math:`T` factor from the right to the left as follows

.. math::
	:nowrap:

	\begin{equation}
		d\Gamma(\alpha \to \beta) \coloneqq dP(\alpha \to \beta) / T = (2\pi)^{3N_{\alpha}-2} V^{1-N_{\alpha}} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta
		\label{eq_rate_of_reaction_master_formula}
	\end{equation}

where :math:`M_{\beta \alpha}` is defined by :math:`\eqref{eq_generic_s_matrix_with_m}` instead of :math:`\eqref{eq_generic_s_matrix_in_time_box}` because we have restored :math:`\delta^4(p_{\beta} - p_{\alpha})`. This is totally wild because by taking rate one typically think of processes that happen within infinitesimal time periods, but we have at the same taken large time limit to recover the :math:`\delta^4(p_{\beta} - p_{\alpha})` factor. Despite the insane derivation, the end result seems reasonable and it will be the key formula that connects S-matrix to experimental measurement of probabilities.

Examples with few initial particles
+++++++++++++++++++++++++++++++++++

One special case of interest is when :math:`N_{\alpha} = 1`, or in other words, processes where one particle decays into multi-particles. In this case :math:`\eqref{eq_rate_of_reaction_master_formula}` becomes

.. math::
	:nowrap:

	\begin{equation}
		d\Gamma(\alpha \to \beta) = 2\pi \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta
		\label{eq_differential_reaction_rate_one_particle}
	\end{equation}

which becomes independent of the volume of the box. This is reasonable because the decay rate of one particle shouldn't care about the size of the containing box. However, the :math:`T \to \infty` limit in :math:`\delta^4(p_{\beta} - p_{\alpha})` is no longer valid. In fact, it cannot be longer than the (mean) lifetime :math:`\tau_{\alpha}` of the particle :math:`\alpha`, because the interaction wouldn't make sense if the particle itself already disintegrates. In this case, in order for :math:`\eqref{eq_time_delta_in_a_period}` to still approximate a delta function, we must assume that any characteristic energy of the interaction satisfies

.. math::
	:nowrap:

	\begin{equation*}
		|E_{\beta} - E_{\alpha}| \ll 1/\tau_{\alpha}
	\end{equation*}

where the right-hand-side is known as the total decay rate.

Another case of interest is when :math:`N_{\alpha} = 2`. In this case :math:`\eqref{eq_rate_of_reaction_master_formula}` takes the following form

.. math::
	:nowrap:

	\begin{equation}
		d\Gamma(\alpha \to \beta) = (2\pi)^4 V^{-1} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta
		\label{eq_differential_reaction_rate_two_particles}
	\end{equation}

It turns out that in the world of experimentalists, it's more common to use, instead of the transition rate, something called *cross-section*, or equivalently, rate per flux, where the flux is defined as [#abuse_of_phi_as_both_state_vector_and_flux]_

.. math::
	:nowrap:

	\begin{equation*}
		\Phi_{\alpha} \coloneqq u_{\alpha} / V
	\end{equation*}

and :math:`u_{\alpha}` is the (relativistic) relative velocity between the two particles, to be discussed in more detail in the next section by considering Lorentz symmetry. We can then rewrite :math:`\eqref{eq_differential_reaction_rate_two_particles}` in terms of the cross-section as follows

.. math::
	:nowrap:

	\begin{equation}
		d\sigma(\alpha \to \beta) \coloneqq d\Gamma(\alpha \to \beta) / \Phi_{\alpha} = (2\pi)^4 u_{\alpha}^{-1} \delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 d\beta
		\label{eq_cross_section_two_particles}
	\end{equation}

Note that :math:`d\sigma` has the dimension of an area.

Lorentz symmetry of rates and cross-sections
++++++++++++++++++++++++++++++++++++++++++++

We can investigate the Lorentz symmetry on the rates and cross-sections as follows. Squaring :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}`, and using the fact that the little group representations are unitary, we see that the following quantity

.. math::
	:nowrap:

	\begin{equation*}
		R_{\beta \alpha} \coloneqq \sum_{\text{spins}} |M_{\beta \alpha}|^2 \prod_{\beta} E \prod_{\alpha} E
	\end{equation*}

is Lorentz invariant, where :math:`E = p_0 = \sqrt{\pbf^2 + m^2}` for each particle in :math:`\alpha` and :math:`\beta`, respectively.

It follows that in the one-particle case, :math:`\eqref{eq_differential_reaction_rate_one_particle}` gives

.. math::
	:nowrap:

	\begin{equation*}
		\sum_{\text{spins}} d\Gamma(\alpha \to \beta) = 2\pi E_{\alpha}^{-1} R_{\beta \alpha} \delta^4(p_{\beta} - p_{\alpha}) \frac{d\beta}{\prod_{\beta} E}
	\end{equation*}

In particular, we recognize :math:`d\beta / \prod_{\beta} E` as a product of the Lorentz invariant :math:`3`-momentum volume elements constructed in :math:`\eqref{eq_lorentz_invariant_3_momentum_volume_element}`. Hence the only factor in the right-hand-side which is not Lorentz invariant is :math:`E_{\alpha}^{-1}`. It follows that the decay rate of a particle, summed up over all spins, is inverse proportional to its energy, or in other words, a faster moving particle decays slower, which is consistent with the special theory of relativity and experimentally observed slow decay rates of high energy particles coming from cosmic rays.

Next, let's turn to the two-particles case. In this case :math:`\eqref{eq_cross_section_two_particles}` gives

.. math::
	:nowrap:

	\begin{equation*}
		\sum_{\text{spins}} d\sigma(\alpha \to \beta) = (2\pi)^4 u_{\alpha}^{-1} E_1^{-1} E_2^{-1} R_{\beta \alpha} \delta^4(p_{\beta} - p_{\alpha}) \frac{d\beta}{\prod_{\beta} E}
	\end{equation*}

where :math:`E_1, E_2` are the energies of the two particles in state :math:`\alpha`. As in the one-particle case, in order for the cross-section to be Lorentz invariant, we must define the relative velocity :math:`u_{\alpha}` such that the product :math:`u_{\alpha} E_1 E_2` is Lorentz invariant. Indeed, such a quantity is uniquely determined by the requirement that when one of the particles stays still, then :math:`u_{\alpha}` should be the velocity of the other particle, and it takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		u_{\alpha} = \frac{\sqrt{(p_1 \cdot p_2)^2 - m_1^2 m_2^2}}{E_1 E_2}
	\end{equation*}

For later use, let's rewrite :math:`u_{\alpha}` in the center-of-mass frame as follows. In the center-of-mass frame, the total momentum vanishes, and therefore we can write :math:`p_1 = (E_1, \pbf)` and :math:`p_2 = (E_2, -\pbf)`. It follows that

.. math::
	:nowrap:

	\begin{align}
		u_{\alpha} &= \frac{\sqrt{(E_1 E_2 + \pbf^2)^2 - m_1^2 m_2^2}}{E_1 E_2} \label{eq_two_particles_relative_velocity_in_center_of_mass_frame} \\
			&= \frac{\sqrt{(E_1 E_2 + \pbf^2)^2 - (E_1^2 - \pbf^2)(E_2^2 - \pbf^2)}}{E_1 E_2} \nonumber \\
			&= \frac{|\pbf| (E_1 + E_2)}{E_1 E_2} \nonumber \\
			&= \left| \frac{\pbf_1}{E_1} - \frac{\pbf_2}{E_2} \right| \nonumber
	\end{align}

which indeed looks more like a relative velocity. Note, however, that this is *not* a physical velocity because its value may approach :math:`2` (i.e., faster than the speed of light) in relativistic limit.

The phase-space factor
++++++++++++++++++++++

By phase-space factor we mean the factor :math:`\delta^4(p_{\beta} - p_{\alpha}) d\beta` that appears in transition probabilities, rates and cross-sections discussed above. The goal of this section is to calculate it, particularly in the scenario where the final state consists of two particles. We'll use the center-of-mass frame with respect to the initial state so that :math:`\pbf_{\alpha} = 0`. Then the phase-space factors can be written as follows

.. math::
	:nowrap:

	\begin{equation*}
		\delta^4(p_{\beta} - p_{\alpha}) d\beta = \delta^3(\pbf'_1 + \pbf'_2 + \cdots) \delta(E'_1 + E'_2 + \cdots - E_{\alpha}) d^3 \pbf'_1 d^3 \pbf'_2 \cdots
	\end{equation*}

where we recall that the primes indicate that the quantities are taken from state :math:`\beta`, and :math:`E_{\alpha}` denotes the total energy of state :math:`\alpha`. In the case where the final state consists of exactly two particles, the phase-space factor can be further simplified as follows

.. math::
	:nowrap:

	\begin{align}
		\delta^4(p_{\beta} - p_{\alpha}) d\beta &= \delta(E'_1 + E'_2 - E_{\alpha}) d^3 \pbf'_1 \label{eq_simplified_two_final_particles_phase_space_factor} \\
			&= \delta \left( \sqrt{|\pbf'_1|^2 + {m'_1}^2} + \sqrt{|\pbf'_1|^2 + {m'_2}^2} - E_{\alpha} \right) |\pbf'_1|^2 d|\pbf'_1| d\Omega \nonumber
	\end{align}

where :math:`\Omega` is the solid angle in :math:`\pbf'_1`-space, if in the integration we replace any occurrence of :math:`\pbf'_2` with :math:`-\pbf'_1`.

To further simply the delta function in :math:`\eqref{eq_simplified_two_final_particles_phase_space_factor}`, we recall the following identity, which is an incarnation of integration by substitution,

.. math::
	:nowrap:

	\begin{equation*}
		\delta(f(x)) = \delta(x - x_0) / f'(x_0)
	\end{equation*}

where :math:`x_0` is a simple zero of :math:`f`. In the case of :math:`\eqref{eq_simplified_two_final_particles_phase_space_factor}`, we make the following choices

.. math::
	:nowrap:

	\begin{align}
		f(|\pbf'_1|) &= \sqrt{|\pbf'_1|^2 + {m'_1}^2} + \sqrt{|\pbf'_1|^2 + {m'_2}^2} - E_{\alpha} \nonumber \\
		k' &= \frac{1}{2E_{\alpha}} \sqrt{\left( E_{\alpha}^2 - {m'_1}^2 - {m'_2}^2 \right)^2 - 4 {m'_1}^2 {m'_2}^2} \label{eq_defn_root_k_prime}
	\end{align}

where :math:`k'` is the unique simple zero of :math:`f`. Differentiating :math:`f` at :math:`k'` we get

.. math::
	:nowrap:

	\begin{equation*}
		f'(k') = \frac{k'}{E'_1} + \frac{k'}{E'_2} = \frac{k' E_{\alpha}}{E_1 E_2}
	\end{equation*}

where

.. math::
	:nowrap:

	\begin{align}
		E'_1 &= \sqrt{{k'}^2 + {m'_1}^2} = \frac{E_{\alpha}^2 + {m'_1}^2 - {m'_2}^2}{2E_{\alpha}} \label{eq_defn_e1_prime} \\
		E'_2 &= \sqrt{{k'}^2 + {m'_2}^2} = \frac{E_{\alpha}^2 - {m'_1}^2 + {m'_2}^2}{2E_{\alpha}} \label{eq_defn_e2_prime}
	\end{align}

Putting all together, we can further simplify :math:`\eqref{eq_simplified_two_final_particles_phase_space_factor}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\delta^4(p_{\beta} - p_{\alpha}) d\beta = \frac{k' E'_1 E'_2}{E_{\alpha}} d\Omega
		\label{eq_two_particles_final_state_phase_factor_formula}
	\end{equation}

where :math:`k', E'_1` and :math:`E'_2` are defined by :math:`\eqref{eq_defn_root_k_prime}, \eqref{eq_defn_e1_prime}` and :math:`\eqref{eq_defn_e2_prime}`, respectively.

Substituting :math:`\eqref{eq_two_particles_final_state_phase_factor_formula}` into :math:`\eqref{eq_differential_reaction_rate_one_particle}`, we see that in the case of one particle decaying into two particles,

.. math::
	:nowrap:

	\begin{equation*}
		\frac{d\Gamma(\alpha \to \beta)}{d\Omega} = \frac{2\pi k' E'_1 E'_2}{E_{\alpha}} |M_{\beta \alpha}|^2
	\end{equation*}

The two-body scattering :math:`1~2 \to 1'~2'`, according to :math:`\eqref{eq_cross_section_two_particles}` and :math:`\eqref{eq_two_particles_relative_velocity_in_center_of_mass_frame}`, takes the following form

.. math::
	:nowrap:

	\begin{equation}
		\frac{d\sigma(\alpha \to \beta)}{d\Omega} = \frac{(2\pi)^4 k' E'_1 E'_2}{u_{\alpha} E_{\alpha}} |M_{\beta \alpha}|^2 \
			= \frac{(2\pi)^4 k' E'_1 E'_2 E_1 E_2}{k E_{\alpha}^2} |M_{\beta \alpha}|^2
		\label{eq_two_body_scattering_cross_section_per_solid_angle}
	\end{equation}

where :math:`k \coloneqq |\pbf_1| = |\pbf_2|`. These calculations will be used in the next section to get some insights into the scattering process.


Implications of the unitarity of S-matrix
+++++++++++++++++++++++++++++++++++++++++

In this section we'll no longer assume the :ref:`Genericity of the S-matrix <assump_genericity_s_matrix>`. This means that we'll get back to use :math:`\eqref{eq_s_matrix_with_m}`, instead of :math:`\eqref{eq_generic_s_matrix_with_m}`, which we recall as follows

.. math::
	:nowrap:

	\begin{equation*}
		S_{\beta \alpha} = \delta(\beta - \alpha) - 2\pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha}
	\end{equation*}

However, all the calculations from the previous sections can still be used here because we'll be caring about, for example, the *total* rates, which are integrations over all possible final states, and the degenerate ones will *not* contribute to such integrals.

First, let's spell out the consequence of the unitarity of the S-matrix, or more precisely :math:`S^{\dagger} S = 1`, as follows

.. math::
	:nowrap:

	\begin{align}
		\delta(\gamma - \alpha) &= \int d\beta ~S^{\ast}_{\beta \gamma} S_{\beta \alpha} \label{eq_s_matrix_unitarity_first_half} \\
			&= \int d\beta \left( \delta(\beta - \gamma) + 2\pi \ifrak \delta^4(p_{\beta} - p_{\gamma}) M^{\ast}_{\beta \gamma} \right) \
				\left( \delta(\beta - \alpha) - 2\pi \ifrak \delta^4(p_{\beta} - p_{\alpha}) M_{\beta \alpha} \right) \nonumber \\
			&= \delta(\gamma - \alpha) + 2\pi \ifrak \delta^4(p_{\alpha} - p_{\gamma}) M^{\ast}_{\alpha \gamma} - \
				2\pi \ifrak \delta^4(p_{\gamma} - p_{\alpha}) M_{\gamma \alpha} \nonumber \\
			&\phantom{=} + 4\pi^2 \delta^4(p_{\gamma} - p_{\alpha}) \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) M^{\ast}_{\beta \gamma} M_{\beta \alpha} \nonumber
	\end{align}

which implies

.. math::
	:nowrap:

	\begin{equation}
		\ifrak M^{\ast}_{\alpha \gamma} - \ifrak M_{\gamma \alpha} + 2\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) M^{\ast}_{\beta \gamma} M_{\beta \alpha} = 0
		\label{eq_s_matrix_unitarity_implication_on_m_general}
	\end{equation}

In the special case where :math:`\alpha = \gamma`, :math:`\eqref{eq_s_matrix_unitarity_implication_on_m_general}` gives the following key identity, known as the *generalized optical theorem*

.. math::
	:nowrap:

	\begin{equation}
		\op{Im} M_{\alpha \alpha} = -\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2
		\label{eq_s_matrix_unitarity_implication_on_m_special}
	\end{equation}

As an application we can calculate the total rate of all transitions produced by the initial state :math:`\alpha` using :math:`\eqref{eq_rate_of_reaction_master_formula}` as follows

.. math::
	:nowrap:

	\begin{align*}
		\Gamma_{\alpha} &\coloneqq \int d\beta ~\frac{d\Gamma(\alpha \to \beta)}{d\beta} \\
			&~= (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 \\
			&~= -\frac{1}{\pi} (2\pi)^{3N_{\alpha} - 2} V^{1 - N_{\alpha}} \op{Im} M_{\alpha \alpha}
	\end{align*}

.. dropdown:: An example of two-body scattering
	:icon: unlock
	:animate: fade-in-slide-down

	In the case where :math:`\alpha` is a two-particles state, we can use :math:`\eqref{eq_cross_section_two_particles}` to calculate the total cross-section as follows

	.. math::
		:nowrap:

		\begin{align}
			\sigma_{\alpha} &\coloneqq \int d\beta ~\frac{d\sigma(\alpha \to \beta)}{d\beta} \label{eq_two_body_total_cross_section_by_m} \\
				&~= (2\pi)^4 u_{\alpha}^{-1} \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\beta \alpha}|^2 \nonumber \\
				&~= -16\pi^3 u_{\alpha}^{-1} \op{Im} M_{\alpha \alpha} \nonumber
		\end{align}

	We then recall from :math:`\eqref{eq_two_body_scattering_cross_section_per_solid_angle}` that :math:`M_{\beta \alpha}` may also be expressed in terms of the differential cross-section by solid angle. Motivated by :math:`\eqref{eq_two_body_scattering_cross_section_per_solid_angle}`, let's define the *scattering amplitude* as follows

	.. math::
		:nowrap:

		\begin{equation*}
			f(\alpha \to \beta) \coloneqq -\frac{4\pi^2}{E_{\alpha}} \sqrt{\frac{k' E'_1 E'_2 E_1 E_2}{k}} ~M_{\beta \alpha}
		\end{equation*}

	so that :math:`d\sigma(\alpha \to \beta) / d\Omega = |f(\alpha \to \beta)|^2`, and the sign is just a convention. A particularly simple case is when the scattering is *elastic*, which means that :math:`E_i = E'_i` and :math:`|\pbf_i| = |\pbf'_i|` for :math:`i = 1, 2`, and consequentially :math:`k = k'`. In this case we have

	.. math::
		:nowrap:

		\begin{align*}
			f(\alpha \to \beta) &= -\frac{4\pi^2 E_1 E_2}{E_{\alpha}} M_{\beta \alpha} \\
			u_{\alpha} &= \frac{k E_{\alpha}}{E_1 E_2}
		\end{align*}

	These, together with :math:`\eqref{eq_two_body_total_cross_section_by_m}`, imply the following so-called `optical theorem <https://en.wikipedia.org/wiki/Optical_theorem>`_ for elastic two-body scattering

	.. math::
		:nowrap:

		\begin{equation}
			\op{Im} f(\alpha \to \alpha) = \frac{k}{4\pi} \sigma_{\alpha}
			\label{eq_optical_theorem}
		\end{equation}

	From the optical theorem we can derive an estimate on the forward diffraction angle as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\sigma_{\alpha} = \int d\Omega ~|f(\alpha \to \beta)|^2 > \frac{1}{2} |f(\alpha \to \alpha)|^2 \Delta\Omega \
				\geq \frac{1}{2} \left| \op{Im} f(\alpha \to \alpha) \right|^2 \Delta\Omega
		\end{equation*}

	where the second inequality follows from the assumption that :math:`f` is continuous, and the factor :math:`1/2` is a rather random choice and can be anything less than :math:`1`. Roughly speaking :math:`\Delta\Omega` measures the peak of the diffraction in the direction of :math:`\alpha`. Combining with :math:`\eqref{eq_optical_theorem}`, we conclude that

	.. math::
		:nowrap:

		\begin{equation*}
			\Delta\Omega < \frac{32 \pi^2}{k^2 \sigma_{\alpha}}
		\end{equation*}

	Assuming that total decay rate :math:`\sigma_{\alpha}` doesn't grow very fast at high energies, the peak in the direction of :math:`\alpha` shrink at the scale of :math:`1 / k^2`.

Another application of the unitary of the S-matrix is along the lines of statistical mechanics. Applying the same calculation in :math:`\eqref{eq_s_matrix_unitarity_first_half}` to :math:`S S^{\dagger} = 1`, we get the counterpart to :math:`\eqref{eq_s_matrix_unitarity_implication_on_m_special}`

.. math::
	:nowrap:

	\begin{equation*}
		\op{Im} M_{\alpha \alpha} = -\pi \int d\beta ~\delta^4(p_{\beta} - p_{\alpha}) |M_{\alpha \beta}|^2
	\end{equation*}

Combining with the master equation :math:`\eqref{eq_rate_of_reaction_master_formula}` we have

.. math::
	:nowrap:

	\begin{equation}
		\int d\beta ~c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} = \int d\beta ~c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha}
		\label{eq_s_matrix_unitarity_induced_symmetry_on_reaction_rates}
	\end{equation}

where :math:`c_{\alpha} \coloneqq \left( V / (2\pi)^3 \right)^{N_{\alpha}}`.

We shall carry out an equilibrium analysis for state :math:`\alpha`. To this end, let :math:`P_{\alpha} d\alpha` be the infinitesimal probability of finding the system in state :math:`\alpha`. Then we have

.. math::
	:nowrap:

	\begin{equation*}
		\frac{dP_{\alpha}}{dt} = \int d\beta ~P_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} - P_{\alpha} \int d\beta ~\frac{d\Gamma(\alpha \to \beta)}{d\beta}
	\end{equation*}

where the first term calculates the total rate that other states transit into :math:`\alpha`, and the second term calculates the total rate that the state :math:`\alpha` transits into other states. Recall that the *entropy* of the system is defined to be

.. math::
	:nowrap:

	\begin{equation*}
		-\int d\alpha ~P_{\alpha} \ln(P_{\alpha} / c_{\alpha})
	\end{equation*}

Its rate of change can be estimated as follows

.. math::
	:nowrap:

	\begin{align*}
		-\frac{d}{dt} \int d\alpha ~P_{\alpha} \ln(P_{\alpha} / c_{\alpha}) \
			&= -\int d\alpha ~(\ln(P_{\alpha} / c_{\alpha}) + 1) \frac{dP_{\alpha}}{dt} \\
			&= -\int d\alpha \int d\beta ~(\ln(P_{\alpha} / c_{\alpha}) + 1) \left( P_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \
	  			- P_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} \right) \\
			&= \int d\alpha \int d\beta ~P_{\beta} \ln\left( \frac{P_{\beta} c_{\alpha}}{P_{\alpha} c_{\beta}} \right) \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \\
			&\geq \int d\alpha \int d\beta ~\left( \frac{P_{\beta}}{c_{\beta}} - \frac{P_{\alpha}}{c_{\alpha}} \right) c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \\
			&= \int d\alpha \int d\beta ~\frac{P_{\beta}}{c_{\beta}} \left( c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} - c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} \right) \\
			&= \int d\alpha ~\frac{P_{\alpha}}{c_{\alpha}} \int d\beta \left( c_{\alpha} \frac{d\Gamma(\alpha \to \beta)}{d\beta} - c_{\beta} \frac{d\Gamma(\beta \to \alpha)}{d\alpha} \right) \xlongequal{\eqref{eq_s_matrix_unitarity_induced_symmetry_on_reaction_rates}} 0
	\end{align*}

where the fourth inequality follows from the general inequality :math:`\ln(x) \geq 1 - 1 / x` for any :math:`x > 0`. This is nothing but the famous slogan: entropy never decreases! And we see that as a consequence of the unitarity of the S-matrix.


Perturbation theory of S-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rather than being the epilogue of :ref:`sec_scattering_theory`, this section is more like a prelude to what comes next. In particular, we will work out a candidate Hamiltonian that satisfies the Lorentz invariance condition discussed in :ref:`sec_s_matrix_and_its_symmetry`.

One possible starting point of the perturbation theory is :math:`\eqref{eq_s_matrix_pre_born_approx}` together with the Lippmann-Schwinger formula :math:`\eqref{eq_lippmann_schwinger_pure}` which we recollect as follows

.. math::
	:nowrap:

	\begin{align}
		S_{\beta \alpha} &= \delta(\beta - \alpha) - 2\pi \ifrak \delta(E_{\beta} - E_{\alpha}) (\Phi_{\beta}, V\Psi_{\alpha}^-) \
		\label{eq_s_matrix_pre_born_approx_repeated} \\
		\Psi_{\alpha}^- &= \Phi_{\alpha} + \int d\beta ~\frac{(\Phi_{\beta}, V\Psi_{\alpha}^-) \Phi_{\beta}}{E_{\alpha} - E_{\beta} + \ifrak \epsilon} \
		\label{eq_lippmann_schwinger_repeated}
	\end{align}

Applying :math:`V` to :math:`\eqref{eq_lippmann_schwinger_repeated}` and taking scalar product with :math:`\Phi_{\beta}`, we get

.. math::
	:nowrap:

	\begin{equation}
		\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right) = V_{\beta \alpha} + \int d\beta ~\frac{\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right) V_{\beta \alpha}}{E_{\alpha} - E_{\beta} + \ifrak \epsilon}
		\label{eq_base_iter_old_fashioned_s_matrix_perturbation}
	\end{equation}

where :math:`V_{\beta \alpha} \coloneqq \left( \Phi_{\beta}, V\Phi_{\alpha} \right)`. One can apply :math:`\eqref{eq_base_iter_old_fashioned_s_matrix_perturbation}` iteratively to get the following

.. math::
	:nowrap:

	\begin{equation}
		\left( \Phi_{\beta}, V\Psi_{\alpha}^- \right) = V_{\beta \alpha} \
	  		+ \int d\gamma ~\frac{V_{\beta \gamma} V_{\gamma \alpha}}{E_{\alpha} - E_{\gamma} + \ifrak \epsilon} \
	  		+ \int d\gamma \int d\gamma' ~\frac{V_{\beta \gamma} V_{\gamma \gamma'} V_{\gamma' \alpha}}{(E_{\alpha} - E_{\gamma} + \ifrak \epsilon)(E_{\alpha} - E_{\gamma'} + \ifrak \epsilon)} \
	  		+ \cdots
		\label{eq_s_matrix_power_series_expansion_old_fashioned}
	\end{equation}

and therefore a power series expansion in :math:`V` of :math:`S_{\beta \alpha}` in view of :math:`\eqref{eq_s_matrix_pre_born_approx_repeated}`.

One obvious drawback of the expansion :math:`\eqref{eq_s_matrix_power_series_expansion_old_fashioned}` is that it obscures the Lorentz symmetry of the S-matrix because the denominators consist of only the energy terms. To overcome this, we shall use instead the other interpretation of the S-matrix in terms of the Hamiltonians given by :math:`\eqref{eq_s_operator_by_u}` and :math:`\eqref{eq_defn_u_operator}`, which we recall as follows

.. math::
	:nowrap:

	\begin{align}
		S &= U(\infty, -\infty) \nonumber \\
		U(\tau, \tau_0) &= \exp(\ifrak H_0 \tau) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \label{eq_defn_u_operator_repeated}
	\end{align}

Differentiating :math:`\eqref{eq_defn_u_operator_repeated}` in :math:`\tau` gives

.. math::
	:nowrap:

	\begin{align}
		\ifrak \frac{d}{d\tau} U(\tau, \tau_0) &= -H_0 \exp(\ifrak H_0 \tau) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \label{eq_evolution_equation_of_u_operator} \\
	  		&\phantom{=} + \exp(\ifrak H_0 \tau) H \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \nonumber \\
			&= \exp(\ifrak H_0 \tau) (H - H_0) \exp(-\ifrak H (\tau - \tau_0)) \exp(-\ifrak H_0 \tau_0) \nonumber \\
			&= \exp(\ifrak H_0 \tau) V \exp(-\ifrak H_0 \tau) U(\tau, \tau_0) \nonumber \\
			&\eqqcolon V(\tau) U(\tau, \tau_0) \nonumber
	\end{align}

Here

.. math::
	:nowrap:

	\begin{equation}
		V(\tau) = \exp(\ifrak H_0 \tau) V \exp(-\ifrak H_0 \tau)
		\label{eq_defn_interaction_perturbation_term}
	\end{equation}

is a time-dependent operator in the so-called *interaction picture*, to be distinguished from the Heisenberg picture operator where the true Hamiltonian :math:`H` should be used in place of :math:`H_0`. The differential equation :math:`\eqref{eq_evolution_equation_of_u_operator}` can be easily solved as follows

.. math::
	:nowrap:

	\begin{equation*}
		U(\tau, \tau_0) = 1 - \ifrak \int_{\tau_0}^{\tau} dt ~V(t) U(t, \tau_0)
	\end{equation*}

which can then be iterated to give the following

.. math::
	:nowrap:

	\begin{align*}
		U(\tau, \tau_0) &= 1 - \ifrak \int_{\tau_0}^{\tau} dt_1 ~V(t_1) \
				+ (-\ifrak)^2 \int_{\tau_0}^{\tau} dt_1 \int_{\tau_0}^{t_1} dt_2 ~V(t_1) V(t_2) \\
	  		&\phantom{=} + (-\ifrak)^3 \int_{\tau_0}^{\tau} dt_1 \int_{\tau_0}^{t_1} dt_2 \int_{\tau_0}^{t_2} dt_3 ~V(t_1) V(t_2) V(t_3) \
				+ \cdots
	\end{align*}

Letting :math:`\tau \to \infty` and :math:`\tau_0 \to -\infty` we get another power series expansion of :math:`S` in :math:`V` as follows

.. math::
	:nowrap:

	\begin{align}
		S &= 1 - \ifrak \int_{-\infty}^{\infty} dt_1 ~V(t_1) \
				+ (-\ifrak)^2 \int_{-\infty}^{\infty} dt_1 \int_{-\infty}^{t_1} dt_2 ~V(t_1) V(t_2) \label{eq_s_matrix_power_series_expansion_raw} \\
	  		&\phantom{=} + (-\ifrak)^3 \int_{-\infty}^{\infty} dt_1 \int_{-\infty}^{t_1} dt_2 \int_{-\infty}^{t_2} dt_3 ~V(t_1) V(t_2) V(t_3) \
				+ \cdots \nonumber
	\end{align}

It's somewhat inconvenient that the integral limits in :math:`\eqref{eq_s_matrix_power_series_expansion_raw}` ruins the permutation symmetry of the products of :math:`V`. But this can be fixed by introducing a *time-ordered product* as follows

.. math::
	:nowrap:

	\begin{align}
		T\{ V(t) \} &\coloneqq V(t)
		\label{eq_defn_time_ordered_product} \\
		T\{ V(t_1) V(t_2) \} &\coloneqq \theta(t_1 - t_2) V(t_1) V(t_2) + \theta(t_2 - t_1) V(t_2) V(t_1) \nonumber \\
		T\{ V(t_1) V(t_2) V(t_3) \} &\coloneqq \theta(t_1 - t_2) \theta(t_2 - t_3) V(t_1) V(t_2) V(t_3) + \cdots \nonumber \\
		&\cdots \nonumber
	\end{align}

where :math:`\theta(\tau)` is the step function which equals :math:`1` for :math:`\tau > 0` and :math:`0` for :math:`\tau < 0`, and it doesn't matter what the value at :math:`\tau = 0` is because it doesn't contribute to the integrals in :math:`\eqref{eq_s_matrix_power_series_expansion_raw}` anyway. With this definition, we can rewrite :math:`\eqref{eq_s_matrix_power_series_expansion_raw}` as follows

.. math::
	:nowrap:

	\begin{equation}
		S = 1 + \sum_{n=1}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 dt_2 \cdots dt_n ~T\{ V(t_1) V(t_2) \cdots V(t_n) \}
		\label{eq_s_matrix_power_series_expansion_time_ordered}
	\end{equation}

where the division by :math:`n!` is to account for the duplicated integrals introduced by the time-ordered product. Note that this power series looks much like the Taylor series of an exponential function. Indeed, in the unlikely event where :math:`V(t)` at different times all commute, one can remove the time-ordering and write :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered}` as an exponential function.

One great benefit of writing :math:`S` as in the form of :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered}` is that we can reformulate the condition of :math:`S` being Lorentz symmetric in terms of some condition on :math:`V`. Recall from :ref:`sec_s_matrix_and_its_symmetry` that a sufficient condition for a Lorentz invariant S-matrix is that the S-operator commutes with :math:`U_0(\Lambda, a)`, or equivalently in infinitesimal terms :math:`\eqref{eq_h0_s_commute}` -- :math:`\eqref{eq_k30_s_commute}` are satisfied. Now the main postulation is to express :math:`V` using a density function as follows

.. math::
	:nowrap:

	\begin{equation}
		V(t) = \int d^3 x ~\Hscr(t, \xbf)
		\label{eq_defn_v_by_density}
	\end{equation}

such that :math:`\Hscr(x)` is a scalar in the sense that

.. math::
	:nowrap:

	\begin{equation}
		U_0(\Lambda, a) \Hscr(x) U^{-1}_0(\Lambda, a) = \Hscr(\Lambda x + a)
		\label{eq_h_density_is_scalar}
	\end{equation}

Under these assumptions, we can further rewrite :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered}` in terms of :math:`\Hscr(x)` as follows

.. math::
	:nowrap:

	\begin{equation}
		S = 1 + \sum_{n=1}^{\infty} \frac{(-\ifrak)^n}{n!} \int d^4 x_1 \cdots d^4 x_n ~T\{ \Hscr(x_1) \cdots \Hscr(x_n) \}
		\label{eq_s_matrix_power_series_expansion_time_ordered_density}
	\end{equation}

This expression of :math:`S` is manifestly Lorentz invariant, except for the time-ordering part. In fact, the time-ordering between two spacetime points :math:`x_1, x_2` are Lorentz invariant if and only if :math:`x_1 - x_2` is time-like, namely, :math:`(x_1 - x_2)^2 \geq 0`. This is consistent with intuition because events with time-like (or light-like) separations may be observed by one observer, who definitely should know which event happened first. Therefore we obtain a sufficient condition for the Lorentz invariance of :math:`S` as follows

.. math::
	:nowrap:

	\begin{equation}
		[\Hscr(x_1), \Hscr(x_2)] = 0, \quad\forall ~(x_1 - x_2)^2 \leq 0
		\label{eq_h_commutativity_for_space_like_separations}
	\end{equation}

where we've also included the light-like case for technical reasons that will only become clear later. This condition is also referred to as the *causality* condition as it may be interpreted as saying that interactions happening at space-like separations should not be correlated.

.. dropdown:: A formal proof of the Lorentz invariance of the S-matrix
	:icon: unlock
	:animate: fade-in-slide-down

	We shall verify that the S-operator defined by :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered_density}` and satisfying :math:`\eqref{eq_h_commutativity_for_space_like_separations}` indeed satisfies :math:`\eqref{eq_h0_s_commute}` -- :math:`\eqref{eq_k30_s_commute}`, or rather, the most crucial one :math:`\eqref{eq_k30_s_commute}`. Using the definition of :math:`\Kbf_0` from :ref:`sec_quantum_lorentz_symmetry`, it follows from :math:`\eqref{eq_h_density_is_scalar}` that

	.. math::
		:nowrap:

		\begin{equation*}
			-\ifrak [\Kbf_0, \Hscr(t, \xbf)] = t \nabla \Hscr(t, \xbf) + \xbf \p_t \Hscr(t, \xbf)
		\end{equation*}

	Integrating over :math:`\xbf` and setting :math:`t = 0`, we get

	.. math::
		:nowrap:

		\begin{align*}
			[\Kbf_0, V] &= \left[ \Kbf_0, \int d^3 x ~\Hscr(0, \xbf) \right] \\
				&= \int d^3 x ~\xbf \left. \frac{\p}{\p t} \right\vert_{t=0} \Hscr(t, \xbf) \\
				&= \left[ H_0, \int d^3 x ~\xbf \Hscr(0, \xbf) \right] \\
				&\eqqcolon [H_0, \Wbf]
		\end{align*}

	where the third equality follows again from :math:`\eqref{eq_h_density_is_scalar}` by letting :math:`\Lambda` to be the infinitesimal time translation. Indeed the quantity

	.. math::
		\Wbf \coloneqq \int d^3 x ~\xbf \Hscr(0, \xbf)

	here is the same as the :math:`\Wbf` in :math:`\eqref{eq_k_as_k0_plus_w}`.

	Now for :math:`\Kbf = \Kbf_0 + \Wbf`, together with :math:`\Pbf` and :math:`H`, to satisfy the commutation relation :math:`\eqref{eq_hkp_commutation}`, which is readily satisfied by the free particle state Lorentz symmetry generators :math:`\Kbf_0, \Pbf_0` and :math:`H_0`, it suffices that

	.. math::
		:nowrap:

		\begin{equation}
			[\Wbf, V] = 0 \label{eq_w_v_commute}
		\end{equation}

	because of the following calculation

	.. math::
		:nowrap:

		\begin{align*}
			[H, \Kbf] &= [H_0 + V, \Kbf_0 + \Wbf] \\
				&= \ifrak \Pbf_0 + [V, \Kbf_0] + [H_0, W] + [V, \Wbf] \\
				&= \ifrak \Pbf_0 \\
				&= \ifrak \Pbf
		\end{align*}

	Finally, we note that :math:`\eqref{eq_w_v_commute}` is equivalent to the following

	.. math::
		:nowrap:

		\begin{equation*}
			\int d^3 x \int d^3 y ~\xbf [\Hscr(0, \xbf), \Hscr(0, \ybf)] = 0
		\end{equation*}

	which is clearly a weaker condition than :math:`\eqref{eq_h_commutativity_for_space_like_separations}`.

At last we've finally climbed the highest peak in scattering theory, namely :math:`\eqref{eq_h_commutativity_for_space_like_separations}`. It is specific to the relativistic theory because time-ordering is always preserved in Galilean symmetry. It is also this restriction that eventually leads us to a quantum field theory. In the words of the author

	"It is this condition that makes the combination of Lorentz invariance and quantum mechanics so restrictive." [#weinberg_quote_on_lorentz_invariance_and_quantum_mechanics]_

	-- S. Weinberg


.. _sec_the_cluster_decomposition_principle:

The Cluster Decomposition Principle
-----------------------------------

Nearly all modern texts on quantum field theory use the so-called *creation* and *annihilation* operators to describe Hamiltonians, but few explain why it is the way it is. Historically speaking, this formalism grew out of the so-called *canonical quantization* of electromagnetic fields. But of course we prefer logical reasons over historical ones, and it is the goal of this chapter to explain how this formalism leads to an S-matrix that satisfies the so-called cluster decomposition principle, or in plain words, that in effect distant experiments yield uncorrelated results. This is a quite fundamental assumption to keep because otherwise it'd be impossible to make an "isolated" experiment.

Bosons and fermions
^^^^^^^^^^^^^^^^^^^

We will now address an issue we left behind from :math:`\eqref{eq_many_particles_state_normalization_rough}`, namely, the permutations of particles in a many-particles state

.. math::
	:nowrap:

	\begin{equation*}
		\Phi_{\pbf_1, \sigma_1, n_1; ~\pbf_2, \sigma_2, n_2; ~\cdots}
	\end{equation*}

Note that we've used the :math:`3`-momenta instead of the :math:`4`-momenta to label the particles since we implicitly assume that the particles are all living on their mass shells. Moreover, we've decided to use the free particles states, which could equally well be in- or out-states. Since there is really no ordering of the particles, it's conceivable that any swap of two particles should just give the same state back. More precisely, there should exist a phase :math:`\alpha = \alpha(\pbf, \sigma, n, \pbf', \sigma', n')`, which depends a priori on the swapping particles, such that

.. math::
	:nowrap:

	\begin{equation}
		\Phi_{\cdots \pbf, \sigma, n; ~\cdots ~\pbf', \sigma', n'; ~\cdots} = \alpha \Phi_{\cdots \pbf', \sigma', n'; ~\cdots ~\pbf, \sigma, n; ~\cdots}
		\label{eq_factor_alpha_for_swapping_two_particles}
	\end{equation}

First of all, let's first argue why :math:`\alpha` should not depend on the other particles that appear in the states. This is in fact another essence of the cluster decomposition principle, namely, what happens between two particles should not depend on other unrelated particles, that in principle can be arbitrarily far away. Next, we argue that :math:`\alpha` should not depend on the spin :math:`z`-component (or the helicity for massless particles). This is because :math:`\alpha` would otherwise have to furnish a :math:`1`-dimensional representation of the rotation group, which, as we've seen in :ref:`Clebsch-Gordan coefficients <dropdown_clebsch_gordan_coefficients>`, doesn't exist. Finally, if :math:`\alpha` would depend on the momenta of the two swapping particles, then the Lorentz invariance would demand that the dependency takes the form of :math:`p^i p'_i` which is symmetric under the swap. Hence we can conclude, by applying :math:`\eqref{eq_factor_alpha_for_swapping_two_particles}` twice, that :math:`\alpha^2 = 1`.

.. warning::
	The argument above that led to the conclusion :math:`\alpha^2 = 1` neglected a possibility that :math:`\alpha` depends on the path that the particle are brought to the momenta :math:`\pbf_1, \pbf_2` and so on. We will come back to this point (much) later.

Now the question has become: should :math:`\alpha` be :math:`1` or :math:`-1`? At this point we shall just make up a story as follows. In this world there exist two types of particles, known as bosons and fermions, such that :math:`\alpha = -1` if the two swapping particles are both fermions and :math:`\alpha = 1` otherwise. This is really a convention rather than any sort of dark magic -- we could have from the beginning agreed upon a rule about how the particles should be ordered and always write states in that order. This convention, however, will turn out to be *mathematically* convenient when we have to deal with symmetries that involve multiple particles, such as the isospin symmetry.

We can now fix the signs in :math:`\eqref{eq_many_particles_state_normalization_rough}`. For the simplicity of notations, we shall write :math:`q \coloneqq (\pbf, \sigma, n)`, when details of the particle states are not important, so that a state can be shorthanded as :math:`\Phi_{q_1 q_2 \cdots q_N}`. In this notation :math:`\eqref{eq_many_particles_state_normalization_rough}` can be written as follows

.. math::
	:nowrap:

	\begin{equation*}
		\left( \Phi_{q'_1 q'_2 \cdots q'_M}, \Phi_{q_1 q_2 \cdots q_N} \right) = \delta_{MN} \sum_{\Pscr} \delta_{\Pscr} \prod_i \delta(q'_i - q_{\Pscr i})
	\end{equation*}

where :math:`\Pscr: \{1, 2, \cdots, N\} \to \{1, 2, \cdots, N\}` is a permutation and :math:`\delta_{\Pscr} = -1` if and only if :math:`\Pscr`, written as a product of swaps, contains an odd number of swaps of fermions.

Creation and annihilation operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a nutshell, creation and annihilation operators provide us a different way to write states like :math:`\Phi_{q_1 \cdots q_N}` and to write operators, e.g., the Hamiltonian, that act on the states. In this section, we shall

1. define the creation and annihilation operators by how they act on particle states,
2. calculate their (anti-)commutator,
3. show how to write *any* operator in terms of the creation and annihilation operators, and
4. work out the Lorentz and CPT transformation laws.

Definition of the operators
+++++++++++++++++++++++++++

Let's start with the creation operator that "creates" a particle with quantum numbers :math:`q` as follows

.. math::
	:nowrap:

	\begin{equation}
		a^{\dagger}(q) \Phi_{q_1 q_2 \cdots q_N} \coloneqq \Phi_{q q_1 q_2 \cdots q_N}
		\label{eq_defn_creation_operator}
	\end{equation}

By introducing a special state :math:`\Phi_{\VAC}`, called the *vacuum state*, which is a state with no particles, we can express any state as follows

.. math::
	:nowrap:

	\begin{equation}
		\Phi_{q_1 q_2 \cdots q_N} = a^{\dagger}(q_1) a^{\dagger}(q_2) \cdots a^{\dagger}(q_N) \Phi_{\VAC}
		\label{eq_particles_state_from_creation_operators}
	\end{equation}

The adjoint of :math:`a^{\dagger}(q)`, denoted by :math:`a(q)`, is then the annihilation operator, which "removes" a particle from the state. Unlike the the creation operator, which according to :math:`\eqref{eq_defn_creation_operator}` always add from the left to the list of existing particles, the annihilation operator necessarily needs to be able to remove the particle from anywhere in the state due to the permutation symmetry discussed in the previous section. To work out the formula for :math:`a(q)`, let's first write down the most general expression as follows

.. math::
	:nowrap:

	\begin{equation}
		a(q) \Phi_{q_1 \cdots q_N} = \sum_{i=1}^N \sigma(i) \delta(q - q_i) \Phi_{q_1 \cdots \hat{q}_i \cdots q_N}
		\label{eq_indeterminate_annihilation_formula}
	\end{equation}

where the hat means that the corresponding term is missing, and :math:`\sigma(i) = \pm 1` are the indeterminants that we need to solve for. Next we pair it with a test state :math:`\Phi_{q'_1 q'_2 \cdots q'_N}` and calculate the result in two ways. The first is a direct calculation using :math:`\eqref{eq_indeterminate_annihilation_formula}`

.. math::
	:nowrap:

	\begin{align*}
		\left( \Phi_{q'_1 \cdots q'_{N-1}}, a(q) \Phi_{q_1 \cdots q_N} \right) &= \sum_{i=1}^N \sigma(i) \delta(q - q_i) \left( \Phi_{q'_1 \cdots q'_{N-1}}, \Phi_{q_1 \cdots \hat{q}_i \cdots q_N} \right) \\
			&= \sum_{i=1}^N \sum _{\Pscr} \sigma(i) \delta(q - q_i) \delta_{\Pscr} \prod_{j=1}^{N-1} \delta(q'_j - q_{\Pscr j})
	\end{align*}

where :math:`\Pscr: \{1, 2, \cdots, N-1\} \to \{1, 2, \cdots, \hat{i}, \cdots, N\}` is a bijection. The second calculation uses the fact that :math:`a` and :math:`a^{\dagger}` are adjoint operators

.. math::
	:nowrap:

	\begin{align*}
		\left( \Phi_{q'_1 \cdots q'_{N-1}}, a(q) \Phi_{q_1 \cdots q_N} \right) &= \left( a^{\dagger}(q) \Phi_{q'_1 \cdots q'_{N-1}}, \Phi_{q_1 \cdots q_N} \right) \\
			&= \left( \Phi_{q q'_1 \cdots q'_{N-1}}, \Phi_{q_1 \cdots q_N} \right) \\
			&= \sum_{i=1}^N \sum_{\Pscr'} \delta_{\Pscr'} \delta(q - q_i) \prod_{j=1}^{N-1} \delta(q'_j - q_{\Pscr' j}) \\
			&= \sum_{i=1}^N \sum_{\Pscr} (\pm 1)^{c_i} \delta_{\Pscr} \delta(q - q_i) \prod_{j=1}^{N-1} \delta(q'_j - q_{\Pscr j})
	\end{align*}

A few notes are in order to explain the above calculation

1. If we think of :math:`q` in :math:`\Phi_{q q'_1 \cdots q'_{N-1}}` as having index :math:`0`, then :math:`\Pscr': \{0, 1, 2, \cdots, N-1\} \to \{1, 2, \cdots, N\}` is a bijection such that :math:`\Pscr' 0 = i` and the rest being the same as :math:`\Pscr`.
2. The sign in :math:`\pm 1` is positive if :math:`q` is a boson and negative if :math:`q` is a fermion.
3. The power :math:`c_i` counts the number of fermions among :math:`q_1, \cdots, q_{i-1}` because the map :math:`\Pscr' 0 = i` can be realized by a product of :math:`i` swaps :math:`(0 \leftrightarrow 1)(1 \leftrightarrow 2) \cdots (i-1 \leftrightarrow i)` and only those swaps with a fermion may contribute a :math:`-1`.

Comparing the results of the two ways of calculating the same quantity, we see that :math:`\sigma(i) = (\pm 1)^{c_i}` and therefore can rewrite :math:`\eqref{eq_indeterminate_annihilation_formula}` as follows

.. math::
	:nowrap:

	\begin{equation}
		a(q) \Phi_{q_1 \cdots q_N} = \sum_{i=1}^N (\pm 1)^{c_i} \delta(q - q_i) \Phi_{q_1 \cdots \hat{q}_i \cdots q_N}
		\label{eq_defn_annihilation_operator}
	\end{equation}

Note that :math:`a(q)` annihilates the vacuum state :math:`\Phi_{\VAC}`, whether :math:`q` is boson or fermion, since there is no state that contains :math:`-1` particles.

.. note::
	Although :math:`a^{\dagger}(q)` and :math:`a(q)` are called the creation and annihilation operators and they indeed appear to add and remove particles from a state, respectively, at least in our framework, it is really more of a mathematical convenience than anything physically realistic -- one should not imagine particles getting created or destroyed like magic.

The (anti-)commutation relation
+++++++++++++++++++++++++++++++

Just like all the operators we've encountered so far, it'll be important to calculate some sort of commutator :math:`\left[ a(q'), a^{\dagger}(q) \right]`. The calculation naturally splits into two halves: the first half is as follows

.. math::
	:nowrap:

	\begin{align*}
		a(q') a^{\dagger}(q) \Phi_{q_1 \cdots q_N} &\xlongequal{\eqref{eq_defn_creation_operator}} a(q') \Phi_{q q_1 \cdots q_N} \\
			&\xlongequal{\eqref{eq_defn_annihilation_operator}} \delta(q' - q) \Phi_{q_1 \cdots q_N} \
	  			+ \blue{\sum_{i=1}^N (\pm 1)^{c'_i} \delta(q' - q_i) \Phi_{q q_1 \cdots \hat{q}_i \cdots q_N}}
	\end{align*}

where the power :math:`c'_i = c_i + 1` if both :math:`q` and :math:`q'` are fermions, and :math:`c'_i = c_i` otherwise.

The second half is more straightforward and is done as follows

.. math::
	:nowrap:

	\begin{align*}
		a^{\dagger}(q) a(q') \Phi_{q_1 \cdots q_N} &\xlongequal{\eqref{eq_defn_annihilation_operator}} a^{\dagger}(q) \sum_{i=1}^N (-1)^{c_i} \delta(q' - q_i) \Phi_{q_1 \cdots \hat{q}_i \cdots q_N} \\
			&\xlongequal{\eqref{eq_defn_creation_operator}} \blue{\sum_{i=1}^N (\pm 1)^{c_i} \delta(q' - q_i) \Phi_{q q_1 \cdots \hat{q}_i \cdots q_N}}
	\end{align*}

Now we would like to combine the two halves to cancel the blue expressions. More precisely, we need to sum up the two if both :math:`q` and :math:`q'` are fermions, and subtract the two otherwise. The result can be formulated as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ a(q'), a^{\dagger}(q) \right]_{\pm} \coloneqq a(q') a^{\dagger}(q) \pm a^{\dagger}(q) a(q') = \delta(q' - q)
		\label{eq_creation_annihilation_commutator}
	\end{equation}

where the sign :math:`\pm` is positive if both :math:`q` and :math:`q'` are fermions, and negative otherwise.

Moreover, one can use the definitions :math:`\eqref{eq_defn_creation_operator}` and :math:`\eqref{eq_defn_annihilation_operator}` to show the following complementing identities

.. math::
	:nowrap:

	\begin{align*}
		\left[ a^{\dagger}(q'), a^{\dagger}(q) \right]_{\pm} &= 0 \\
		\left[ a(q'), a(q) \right]_{\pm} &= 0
	\end{align*}

with the same sign convention as in :math:`\eqref{eq_creation_annihilation_commutator}`.

A universal formula of operators
++++++++++++++++++++++++++++++++

We will show that any operator (on states) :math:`\Ocal` can be written as a sum of products of creation and annihilation operators as follows

.. math::
	:nowrap:

	\begin{align}
		\Ocal &= \sum_{N=0}^{\infty} \sum_{M=0}^{\infty} \int dq'_1 \cdots dq'_N dq_1 \cdots dq_M \label{eq_general_operator_expansion_in_creation_and_annihilation} \\
			&\phantom{=} \times a^{\dagger}(q'_1) \cdots a^{\dagger}(q'_N) a(q_M) \cdots a(q_1) \nonumber \\
			&\phantom{=} \times C_{NM}(q'_1, \cdots, q'_N, q_1, \cdots, q_M) \nonumber
	\end{align}

where :math:`C_{NM}(q'_1, \cdots, q'_N, q_1, \cdots, q_M)` are the coefficients to be determined. Indeed the coefficients :math:`C_{NM}` can be worked out by induction as follows. The base case if when :math:`N = M = 0`, where we simply define

.. math::
	:nowrap:

	\begin{equation}
		C_{00} \coloneqq \left( \Phi_{\VAC}, \Ocal \Phi_{\VAC} \right)
		\label{eq_defn_c00}
	\end{equation}

Now suppose inductively that :math:`C_{NM}` have been defined for all :math:`N < L, M \leq K` or :math:`N \leq L, M < K`. Then one calculates

.. math::
	:nowrap:

	\begin{align*}
		\left( \Phi_{q'_1 \cdots q'_L}, \Ocal \Phi_{q_1 \cdots q_K} \right) &= L!~K!~C_{LK}(q'_1 \cdots q'_L, q_1, \cdots q_K) \\
			&\phantom{=} + \text{ terms involving } C_{NM} \text{ with } N < L, M \leq K \text{ or } N \leq L, M < K
	\end{align*}

where the factorials :math:`L!` and :math:`K!` are included to account for the total permutations of :math:`q'_1, \cdots, q'_L` and :math:`q_1, \cdots q_K`, respectively. The coefficient :math:`C_{LK}` is therefore uniquely determined, and hence all :math:`C_{NM}` by induction.

.. note::
	In :math:`\eqref{eq_general_operator_expansion_in_creation_and_annihilation}` we've chosen a specific ordering of the creation and annihilation operators, known as the *normal* order. Namely, all the creation operators lie to the left of all the annihilation operators. It has at least one advantage of making :math:`\eqref{eq_defn_c00}` obvious. Finally we note that any ordering of a composition of creation and annihilation operators can be normally ordered by applying :math:`\eqref{eq_creation_annihilation_commutator}`.

.. _sec_the_lorentz_and_cpt_transformation_laws:

The Lorentz and CPT transformation laws
+++++++++++++++++++++++++++++++++++++++

Let's first work out how :math:`a^{\dagger}(\pbf, \sigma, n)` and :math:`a(\pbf, \sigma, n)` transform under Lorentz transformations. To this end, we recall the general Lorentz transformation law on free-particles state :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`, and use :math:`b` for the translational parameter to avoid conflict of notations, as follows

.. math::
	:nowrap:

	\begin{align*}
		U_0(\Lambda, b) \Phi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} =&~ \exp(-\ifrak b^{\mu} ((\Lambda p_1)_{\mu} + (\Lambda p_2)_{\mu} + \cdots)) \\
		&\times \sqrt{\frac{(\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p_1)_0 (p_2)_0 \cdots}} \\
		&\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma'_1 \sigma_1}(W_1(\Lambda, p_1)) D_{\sigma'_2 \sigma_2}(W_2(\Lambda, p_2)) \cdots \\
		&\times \Phi_{\Lambda p_1, \sigma'_1, n_1; ~\Lambda p_2, \sigma'_2, n_2; ~\cdots}
	\end{align*}

where we also recall that :math:`U_0` is the Lorentz transformation on free-particles state. Expanding :math:`\eqref{eq_particles_state_from_creation_operators}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\Phi_{\pbf_1, \sigma_1, n_1;~\pbf_2, \sigma_2, n_2;~\cdots} = a^{\dagger}(\pbf_1, \sigma_1, n_1) a^{\dagger}(\pbf_2, \sigma_2, n_2) \cdots \Phi_{\VAC}
	\end{equation*}

and imposing the assumption that the vacuum state is fixed by *any* Lorentz transformation

.. math::
	:nowrap:

	\begin{equation*}
		U_0(\Lambda, b) \Phi_{\VAC} = \Phi_{\VAC}
	\end{equation*}

we see that :math:`a^{\dagger}(\pbf, \sigma, n)` better transforms as follows

.. math::
	:nowrap:

	\begin{equation}
		U_0(\Lambda, b) a^{\dagger}(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) = \
			\exp(-\ifrak b^{\mu} (\Lambda p)_{\mu}) \sqrt{\frac{(\Lambda p)_0}{p_0}} D_{\sigma' \sigma}(W(\Lambda, p)) a^{\dagger}(\pbf_{\Lambda}, \sigma', n)
		\label{eq_lorentz_transformation_formula_for_creation_operator}
	\end{equation}

where :math:`\pbf_{\Lambda}` is the :math:`3`-momentum part of :math:`\Lambda p`. The corresponding transformation law for the annihilation operator :math:`a(\pbf, \sigma, n)` can be obtained by taking the adjoint of :math:`\eqref{eq_lorentz_transformation_formula_for_creation_operator}` and remembering that :math:`U_0` is unitary as follows

.. math::
	:nowrap:

	\begin{equation}
		U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) = \
			\exp(\ifrak b^{\mu} (\Lambda p)_{\mu}) \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{\ast}_{\sigma' \sigma}(W(\Lambda, p)) a(\pbf_{\Lambda}, \sigma, n)
		\label{eq_lorentz_transformation_formula_for_annihilation_operator}
	\end{equation}

The parity, time, and charge symmetry on the creation and annihilation operators can be worked out as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		U(\Pcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Pcal) &\xlongequal{\eqref{eq_space_inversion_on_massive_general}} \eta_n a^{\dagger}(-\pbf, \sigma, n) &&\quad\text{(massive particle)}
		\label{eq_creation_operator_space_inversion_conjugation_massive} \\
		U(\Pcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Pcal) &\xlongequal{\eqref{eq_space_inversion_on_massless_general}} \eta_{\sigma, n} \exp(\mp \ifrak \pi \sigma) a^{\dagger}(-\pbf, -\sigma, n) &&\quad\text{(massless particle)}
		\label{eq_creation_operator_space_inversion_conjugation_massless} \\
		U(\Tcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Tcal) &\xlongequal{\eqref{eq_time_inversion_on_massive_general}} \zeta_n (-1)^{\jfrak - \sigma} a^{\dagger}(-\pbf, -\sigma, n) &&\quad\text{(massive particle)}
		\label{eq_creation_operator_time_inversion_conjugation_massive} \\
		U(\Tcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Tcal) &\xlongequal{\eqref{eq_time_inversion_on_massless_general}} \zeta_{\sigma, n} \exp(\pm \ifrak \pi \sigma) a^{\dagger}(-\pbf, \sigma, n) &&\quad\text{(massless particle)}
		\label{eq_creation_operator_time_inversion_conjugation_massless} \\
		U(\Ccal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Ccal) &\xlongequal{\phantom{(00)}} \xi_n a^{\dagger}(\pbf, \sigma, n^c) &&\quad\text{(massive/massless particle)}
		\label{eq_creation_operator_charge_inversion_conjugation}
	\end{alignat}

where :math:`\Ccal` replaces a particle of species :math:`n` with its antiparticle -- a notion we haven't really explained yet -- :math:`n^c`. The corresponding transformation laws for :math:`a(\pbf, \sigma, n)` can be derived from the above identities by taking the adjoint, and are omitted here.

Cluster decomposition of S-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We are now ready to formalize the cluster decomposition principle, which states that experiments conducted at places far away from each other -- a notion which will be made more precise later in the section -- should yield uncorrelated results, in terms of the S-matrix.

Recall that :math:`S_{\beta \alpha} \coloneqq \left( \Psi^+_{\beta}, \Psi^-_{\alpha} \right)` is the amplitude of a reaction with in-state :math:`\alpha` and out-state :math:`\beta`, both of which are compound indexes of many-particles states. Now suppose the particles admit partitions :math:`\alpha = \alpha_1 \sqcup \cdots \sqcup \alpha_N` and :math:`\beta = \beta_1 \sqcup \cdots \sqcup \beta_N`, respectively, so that any particle from :math:`\alpha_i \cup \beta_i` is spatially far away from any particle from :math:`\alpha_j \cup \beta_j` for :math:`i \neq j`, or in plain words, the in- and out-state particles form :math:`N` clusters that are spatially far away from each other. Then the cluster decomposition principle demands a corresponding splitting of the S-matrix as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = \prod_{i=1}^N S_{\beta_i \alpha_i}
		\label{eq_s_matrix_splitting}
	\end{equation}

It is, however, not clear from :math:`\eqref{eq_s_matrix_splitting}` what conditions on :math:`S_{\beta \alpha}` would guarantee the cluster decomposition principle. To this end, we'll introduce a recursive expansion of :math:`S_{\beta \alpha}` in terms of the so-called "connected parts" :math:`S^C_{\beta \alpha}`, which can be best understood visually within the framework of Feynman diagrams to be introduced later. Roughly speaking, the idea is to decompose :math:`S_{\beta \alpha}` into components that involve only a (possibly full) subset of the particles.

Let's work out the recursive definitions of :math:`S^C_{\beta \alpha}` from bottom up. For simplicity, we'll in addition make the following assumption

.. _assump_on_particle_stability:

.. admonition:: Assumption on particle stability
	:class: Important

	The single-particle states are stable, i.e., it's not possible for a single-particle state to transit into any other state, including the vacuum state.

Under this assumption we can define the base case as follows

.. math::
	:nowrap:

	\begin{equation*}
		S^C_{q' q} \coloneqq S_{q' q} = \delta(q' - q)
	\end{equation*}

For a two-body interaction, we can define :math:`S^C_{q'_1 q'_2,~q_1 q_2}` by the following identity

.. math::
	:nowrap:

	\begin{equation}
		S_{q'_1 q'_2,~q_1 q_2} = S^C_{q'_1 q'_2,~q_1 q_2} + \delta(q'_1 - q_1) \delta(q'_2 - q_2) \pm \delta(q'_1 - q_2) \delta(q'_2 - q_1)
		\label{eq_defn_connected_part_two_body}
	\end{equation}

where the sign :math:`\pm` is negative if both :math:`q_1` and :math:`q_2` are fermions, and positive otherwise. This should be interpreted as saying that :math:`S^C_{q'_1 q'_2,~q_1 q_2}` is the amplitude that the two particles actually interact, and the rest is the amplitude that they do not interact at all. Indeed :math:`\eqref{eq_defn_connected_part_two_body}` can be seen as another incarnation of :math:`\eqref{eq_s_matrix_with_m}`.

This procedure can be continued to, say, a three-body problem as follows

.. math::
	:nowrap:

	\begin{align*}
		S_{q'_1 q'_2 q'_3,~q_1 q_2 q_3} &= S^C_{q'_1 q'_2 q'_3,~q_1 q_2 q_3} \\
			&\phantom{=} + \delta(q'_1 - q_1) S^C_{q'_2 q'_3,~q_2 q_3} \pm \text{ permutations} \\
			&\phantom{=} + \delta(q'_1 - q_1) \delta(q'_2 - q_2) \delta(q'_3 - q_3) \pm \text{ permutations}
	\end{align*}

and in general by recursion as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = S^C_{\beta \alpha} + \sum_{\substack{\alpha = \alpha_1 \sqcup \alpha_2 \sqcup \cdots \\ \beta = \beta_1 \sqcup \beta_2 \sqcup \cdots}} \
			(\pm) S^C_{\beta_1 \alpha_1} S^C_{\beta_2 \alpha_2} \cdots
		\label{eq_s_matrix_recursive_by_connected_parts}
	\end{equation}

where the sum is taken over all nontrivial partitions of :math:`\alpha` and :math:`\beta`.

The upshot of all the trouble of writing :math:`S_{\beta \alpha}` as a sum of connected parts is that the cluster decomposition principle :math:`\eqref{eq_s_matrix_splitting}` can now be rephrased as :math:`S^C_{\beta \alpha} = 0` if any two particles in :math:`\alpha \cup \beta` are (spatially) far apart. This can be illustrated by the following example of a four-body reaction :math:`1234 \to 1'2'3'4'` such that the particles :math:`\{1,2,1',2'\}` are far away from :math:`\{3,4,3',4'\}`. With obvious shorthand notations, we have

.. math::
	:nowrap:

	\begin{align*}
		S_{1'2'3'4',~1234} &= S^C_{1'2',~12} S^C_{3'4',~34} \\
			&\phantom{=} + (\delta_{1'1}\delta_{2'2} \pm \delta_{1'2}\delta_{2'1}) S^C_{3'4',~34} \\
			&\phantom{=} + (\delta_{3'3}\delta_{4'4} \pm \delta_{3'4}\delta_{4'3}) S^C_{1'2',~12} \\
			&\phantom{=} + (\delta_{1'1}\delta_{2'2} \pm \delta_{1'2}\delta_{2'1})(\delta_{3'3}\delta_{4'4} \pm \delta_{3'4}\delta_{4'3}) \\
			&= S_{1'2',~12} S_{3'4',~34}
	\end{align*}

as expected.

Finally, let's quantify the phrase "spatially far apart", which we've been using so far without much explanation. Since our states are defined in the momentum space, we first need to translate the amplitude into spatial coordinates using Fourier transform as follows

.. math::
	:nowrap:

	\begin{align*}
		S^C_{\xbf_1' \xbf'_2 \cdots,~\xbf_1 \xbf_2 \cdots} &\coloneqq \int d^3 \pbf'_1 d^3 \pbf'_2 \cdots d^3 \pbf_1 d^3 \pbf_2 \cdots S^C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots} \\
			&\phantom{\coloneqq} \times \exp(\ifrak \pbf'_1 \cdot \xbf'_1) \exp(\ifrak \pbf'_2 \cdot \xbf'_2) \cdots \exp(-\ifrak \pbf_1 \cdot \xbf_1) \exp(-\ifrak \pbf_2 \cdot \xbf_2) \cdots
	\end{align*}

where we've kept only the momentum of each particle because the other quantum numbers do not play a role in the current discussion.

Now if :math:`S^C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots}` were (reasonably) smooth, e.g., Lebesgue integrable, then the `Riemann-Lebesgue lemma <https://en.wikipedia.org/wiki/Riemann%E2%80%93Lebesgue_lemma>`_ asserts that the left-hand-side :math:`S^C_{\xbf'_1 \xbf'_2 \cdots,~\xbf_1 \xbf_2 \cdots}` vanishes as any linear combinations of :math:`\xbf'_1, \xbf'_2, \cdots, \xbf_1, \xbf_2, \cdots` go to infinity. This is a slightly too strong constraint because we know from the translational invariance that the left-hand-side shouldn't change by an overall translation, no matter how large it is. To remedy this defect, we introduce the factors of energy and momentum conservation delta functions in :math:`S^C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots}` as follows

.. math::
	:nowrap:

	\begin{align*}
		S^C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots} &= \delta^4(p' - p) C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots} \\
			&= \delta(E'_1 + E'_2 + \cdots - E_1 -E_2 - \cdots) \delta^3(\pbf'_1 + \pbf'_2 + \cdots - \pbf_1 - \pbf_2 - \cdots) C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots}
	\end{align*}

which guarantees that an overall translation will not change the integral. Moreover, we see that, in fact, the remaining :math:`C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots}` cannot contain any more delta functions of the momenta, such as e.g., :math:`\delta^3(\pbf'_1 - \pbf_1)`, because otherwise one could translate a subset of the particles, such as e.g., :math:`\{\xbf'_1, \xbf_1\}`, far away from the others, while keeping their relative position fixed, and henceforth leaving :math:`S^C_{\pbf'_1 \pbf'_2 \cdots,~\pbf_1 \pbf_2 \cdots}` unchanged. But this would violate the cluster decomposition principle. All in all, we've arrived at the following key conclusion

	The cluster decomposition principle is equivalent to the condition that every connected parts of the S-matrix contain exactly one momentum-conservation delta function.

at least under the :ref:`assumption on particle stability <assump_on_particle_stability>`.

.. dropdown:: Inevitability of quantum field theory
	:icon: unlock
	:animate: fade-in-slide-down

	One of the deepest consequences of the cluster decomposition principle is the necessity of a quantum field theory, or in other words, the inability to consistently define scattering amplitudes with a fixed number of particles. For example, let's consider a two-body problem with S-matrix amplitude :math:`S_{q'_1 q'_2,~q_1 q_2}`. Then it'd have been possible to adjust the Hamiltonian so that all higher-order amplitudes, and in particular, the three-body amplitude, vanish.

	.. math::
		:nowrap:

		\begin{equation*}
			S_{q'_1 q'_2 q'_3,~q_1 q_2 q_3} = S^C_{q'_1 q'_2 q'_3,~q_1 q_2 q_3} + S^C_{q'_1 q'_2,~q_1 q_2} \delta(q'_3 - q_3) \pm \text{ permutations} = 0
		\end{equation*}

	which implies that

	.. math::
		:nowrap:

		\begin{equation*}
			S^C_{q'_1 q'_2 q'_3,~q_1 q_2 q_3} = -S^C_{q'_1 q'_2,~q_1 q_2} \delta(q'_3 - q_3) \mp \text{ permutations}
		\end{equation*}

	But since the two-body amplitude :math:`S^C_{q'_1 q'_2,~q_1 q_2}` itself also contains a momentum-conservation delta function, the first summand of :math:`S^C_{q'_1 q'_2 q'_3,~q_1 q_2 q_3}` would contain a product of two delta functions, in violation against the cluster decomposition principle.

.. _sec_cluster_decomposable_hamiltonians:

Cluster decomposable Hamiltonians
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have seen from the previous section how the cluster decomposition principle is equivalent to a "single momentum-conservation delta function" condition on the connected parts of the S-matrix. The goal of this section is to further translate this condition to one on the Hamiltonian, as a prelude to the next chapter. Writing the Hamiltonian in its most general form as follows

.. math::
	:nowrap:

	\begin{equation}
		H = \sum_{N=0}^{\infty} \sum_{M=0}^{\infty} \int dq'_1 \cdots dq'_N dq_1 \cdots dq_M~a^{\dagger}(q'_1) \cdots a^{\dagger}(q'_N) a(q_M) \cdots a(q_1) \
			h_{NM}(q'_1, \cdots, q'_N, q_1 \cdots, q_M)
		\label{eq_general_expansion_of_hamiltonian}
	\end{equation}

we claim that the S-matrix corresponding to :math:`H` satisfies the clustering decomposition principle if each of the coefficients :math:`h_{NM}` contains exactly one momentum-conservation delta function.

To this end, we'll use the time-dependent perturbation of the S-matrix given by :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered}`, rephrased in terms of the matrix entries as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = \sum_{n=0}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 \cdots dt_n \
			\left( \Phi_{\beta}, T\{V(t_1) \cdots V(t_n)\} \Phi_{\alpha} \right)
		\label{eq_time_dependent_perturbative_s_matrix_element}
	\end{equation}

where we also recall :math:`V(t) \coloneqq \exp(\ifrak H_0 t) V \exp(-\ifrak H_0 t)` and :math:`H = H_0 + V`. Now we remember the following facts

1. Both :math:`\Phi_{\alpha}` and :math:`\Phi_{\beta}` can be written as a number of creation operators applied to :math:`\Phi_{\VAC}`.
2. The creation operators from :math:`\Phi_{\beta}` can be moved to the front of :math:`T\{V(t_1) \cdots V(t_n)\} \Phi_{\alpha}` by adjoint, and become annihilation operators.
3. Each :math:`V(t)` can be written as a polynomial in creation and annihilation operators as in :math:`\eqref{eq_general_operator_expansion_in_creation_and_annihilation}`.
4. Creation and annihilation operators satisfy the canonical commutation relations :math:`\eqref{eq_creation_annihilation_commutator}`, which can be reorganized as :math:`a(q') a^{\dagger}(q) = \pm a^{\dagger}(q) a(q') + \delta(q' - q)` to highlight the effect of moving a creation operator left over an annihilation operator -- it produces a delta function. Moreover, since we'll only care about the momentum in this section, we'll assume that all the Kronecker deltas of discrete quantum numbers have already been summed up , so that :math:`\delta(q' - q)` becomes really :math:`\delta(p' - p)`, or even :math:`\delta(\pbf' - \pbf)` assuming the particles are on their mass shells.

One particularly convenient way to calculate :math:`\eqref{eq_time_dependent_perturbative_s_matrix_element}` is to move all the creation operators to the left of the annihilation operators, while collecting all the delta functions along the way. In the end, the only nonzero term is a vacuum expectation value which is a polynomial of these delta functions. In other words, none of the creation and annihilation operators will survive because otherwise by construction, the rightmost operator is necessarily an annihilation operator, which then would be vanishing by acting on :math:`\Phi_{\VAC}`. This procedure can be rather messy, but luckily, there exists already a convenient bookkeeping device, known as the *Feynman diagrams*. We will encounter Feynman diagrams many times going forward as it is such a convenient tool, and here we'll just focus on the momentum conservations.

Let's write :math:`S^{(n)}_{\beta \alpha} \coloneqq \left( \Phi_{\beta}, T\{V(t_1) \cdots V(t_n)\} \Phi_{\alpha} \right)`. Then the general recipe for constructing a Feynman diagram to keep track of the delta functions in :math:`S^{(n)}_{\beta \alpha}` consists of the following steps.

.. warning::
	The Feynman diagram described here will *not* be very useful in quantum-theoretic calculations since, as we'll see in the next chapter, the interaction densities can only be constructed out of quantum fields, rather than the individual annihilation and creation operators. The only purpose of the diagrams to be described below is to illustrate the cluster decomposition principle in terms of the Hamiltonian. The generally useful Feynman diagrams will be introduced in :ref:`sec_the_feynman_rules`.

1. Orient the paper on which the diagram will be drawn so that the time direction goes upwards. (This is of course a random choice just to set up the scene.)
2. Place as many vertical strands as the number of particles in :math:`\Phi_{\alpha}` at the bottom of the page. Do the same to :math:`\Phi_{\beta}` but at the top of the page. In actual Feynman diagrams the strands will be oriented according to the flow of time (and a distinction between particles and antiparticles), but it's not important for our purposes here -- we will only use the orientations to guarantee that an annihilation operator is paired with a creation operator. For definiteness, let's orient these strands to point upwards.
3. Place one vertex (i.e., a fat dot) for each :math:`V(t)`, or more precisely, each monomial in annihilation and creation operators in :math:`V(t)`, with as many incoming strands as the number of annihilation operators, and outgoing strands as the number of creation operators. Note that since we'll only be doing a momentum conservation analysis, the details of :math:`V(t)`, e.g., the coefficients in front of the product of the annihilation and creation operators, which will be worked out in the next chapter, is not important here.
4. To each strand from the previous two steps, associate a :math:`3`-momentum :math:`\pbf` of the corresponding particle. To each vertex, associate a momentum-conservation delta function so that the total incoming momenta equals the total outgoing momenta. This is in fact a consequence of an assumption on :math:`V(t)`, which is related to our assumption on :math:`h_{NM}` in :math:`\eqref{eq_general_expansion_of_hamiltonian}`, and will be elaborated on later in this section.
5. Connect all the loose ends of the strands in a way that is compatible with the orientations on the strands.
6. To each edge from the previous step, which connects two strands with momenta :math:`\pbf` and :math:`\pbf'`, respectively, associate a delta function :math:`\delta(\pbf' - \pbf)`, which comes from the canonical commutation relation :math:`\eqref{eq_creation_annihilation_commutator}`.

Finally :math:`S^{(n)}_{\beta \alpha}` is simply a sum over products of delta functions, one for each such diagram.

As an example to illustrate the above steps, let's consider a five-body scattering with all particle identical to each other except for their momenta as follows

.. math::
	:nowrap:

	\begin{align*}
		\Phi_{\alpha} &= a^{\dagger}(\pbf_5) a^{\dagger}(\pbf_4) a^{\dagger_3}(\pbf) a^{\dagger}(\pbf_2) a^{\dagger}(\pbf_1) \Phi_{\VAC} \\
		\Phi_{\beta} &= a^{\dagger}(\pbf_{14}) a^{\dagger}(\pbf_{13}) a^{\dagger}(\pbf_{12}) a^{\dagger}(\pbf_{11}) a^{\dagger}(\pbf_{10}) \Phi_{\VAC} \\
		V(t) &= \delta^3(\pbf_8 + \pbf_9 - \pbf_6 - \pbf_7) a^{\dagger}(\pbf_9) a^{\dagger}(\pbf_8) a(\pbf_7) a(\pbf_6)
	\end{align*}

where the coefficient of :math:`V(t)`, except for the momentum conserving delta function, is suppressed. The figure below illustrates a few summands of the third order :math:`S^{(3)}_{\beta \alpha} = \left( \Phi_{\beta}, T\{V(t_1)V(t_2)V(t_3)\}\Phi_{\alpha} \right)`.

.. figure:: ./static/quantum-theory-of-fields/momenta-feynman-diagram.svg
	:align: center

	Figure: A few Feynman diagrams and the corresponding products of delta functions. All the connecting strands are directed upwards, i.e., in the positive time direction, and hence we've omitted the arrows.

Here we've used shorthand notations such as :math:`\delta_{1,6'} \coloneqq \delta^3(\pbf_1 - \pbf'_6)` and :math:`\delta_{8,9 - 6,7} \coloneqq \delta^3(\pbf_8 + \pbf_9 - \pbf_6 - \pbf_7)`.

Notice that some summands of :math:`S^{(n)}_{\beta \alpha}` are connected, in the sense that the graph, viewed as undirected, is connected, and some, for example the second one, are disconnected. For the disconnected ones, it's obvious that the product of the delta functions splits into a product of products, one for each connected components. Therefore we can rewrite :math:`S^{(n)}_{\beta \alpha}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		S^{(n)}_{\beta \alpha} = \sum_{\op{PART}'} (\pm) \prod_{j=1}^{\kappa} \left( \Phi_{\beta_j}, T\{V(t_{j1}) \cdots V(t_{jn_j})\} \Phi_{\alpha_j} \right)_C
	\end{equation*}

where the sum is taken over all partitions :math:`\op{PART}'` of :math:`\alpha = \sqcup_{j=1}^{\kappa} \alpha_j`, :math:`\beta = \sqcup_{j=1}^{\kappa} \beta_j` and :math:`\{1,2,\cdots,n\} = \sqcup_{j=1}^{\kappa} \{j1, \cdots, jn_j\}`, and the subscript :math:`C` indicates that only connected diagrams are allowed. Such refactorization is possible because in the eventual evaluation of a Feynman diagram, such as the ones worked out in the above example, different components are essentially independent of each other.

When evaluating :math:`\eqref{eq_time_dependent_perturbative_s_matrix_element}` using Feynman diagrams as explained above, we note that all the :math:`V(t)`'s in :math:`T\{V(t_1) \cdots V(t_n)\}` are interchangeable. It follows that a partition :math:`\op{PART}'` can be split up into a partition :math:`\op{PART}` of :math:`\alpha` and :math:`\beta` respectively into :math:`\kappa` clusters, and, modulo permutations, a partition :math:`n = n_1 + \cdots n_{\kappa}`. We can then evaluate :math:`\eqref{eq_time_dependent_perturbative_s_matrix_element}` as follows

.. math::
	:nowrap:

	\begin{align*}
		S_{\beta \alpha} &= \sum_{n=0}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 \cdots dt_n \left( \Phi_{\beta}, T\{V(t_1) \cdots V(t_n) \} \Phi_{\alpha} \right) \\
			&= \sum_{n=0}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 \cdots dt_n \sum_{\op{PART}'} (\pm) \prod_{j=1}^{\kappa} \left( \Phi_{\beta_j}, T\{V(t_{j1}) \cdots V(t_{jn_j})\} \Phi_{\alpha_j} \right)_C \\
			&= \sum_{n=0}^{\infty} \frac{(-\ifrak)^n}{n!} \sum_{\op{PART}} (\pm) \sum_{n_1 + \cdots n_{\kappa} = n} \frac{n!}{n_1! \cdots n_{\kappa}!} \prod_{j=1}^{\kappa} \int_{-\infty}^{\infty} dt_{j1} \cdots dt_{jn_j} \left( \Phi_{\beta_j}, T\{V(t_{j1}) \cdots V(t_{jn_j}) \Phi_{\alpha_j}\} \Phi_{\alpha_j} \right)_C \\
			&= \sum_{\op{PART}} (\pm) \prod_{j=1}^{\kappa} \sum_{n_j=0}^{\infty} \frac{(-\ifrak)^{n_j}}{n_j!} \int_{-\infty}^{\infty} dt_{j1} \cdots dt_{jn_j} \left( \Phi_{\beta_j}, T\{V(t_{j1}) \cdots V(t_{jn_j})\} \Phi_{\alpha_j} \right)_C
	\end{align*}

Comparing with :math:`\eqref{eq_s_matrix_recursive_by_connected_parts}`, we see that

.. math::
	:nowrap:

	\begin{equation*}
		S^C_{\beta \alpha} = \sum_{n=0}^{\infty} \frac{(-\ifrak)^n}{n!} \int_{-\infty}^{\infty} dt_1 \cdots dt_n \left( \Phi_{\beta}, T\{V(t_1) \cdots V(t_n)\} \Phi_{\alpha} \right)_C
	\end{equation*}

which also justifies calling :math:`S^C_{\beta \alpha}` a connected part of :math:`S_{\beta \alpha}`, because it corresponds to connected Feynman diagrams.

It remains to argue that :math:`S^C_{\beta \alpha}` contains exactly one momentum-conservation delta function, assuming that the same applies to the coefficients :math:`h_{NM}` in :math:`\eqref{eq_general_expansion_of_hamiltonian}`. Indeed, since :math:`H = H_0 + V` and the single momentum-conservation condition holds automatically true for :math:`H_0`, the same holds for :math:`V`. In other words, each vertex in the Feynman diagram produces one single momentum-conservation delta function, as promised earlier.

Finally, we note that the fact that each connected Feynman diagram gives rise to exactly one momentum-conservation delta function is the consequence of an elimination process. More precisely, we can first get rid of delta functions on the internal edges, i.e., edges between vertices, by equating the momenta on the two ends of the edge. So we're left with as many delta functions as there are vertices. Then for each internal edge, we can choose one of its endpoints, and solve the associated momentum in terms of the other momenta. Hence we can get rid of all but one delta function, as claimed.

An example elimination process is presented for the following Feynman diagram.

.. figure:: ./static/quantum-theory-of-fields/pseudo-feynman-diagram.svg
	:align: center

	Figure: An example Feynman diagram.

We skip the first step of equating momenta on the ends of edges, and continue with the conservation delta functions on the vertices as follows

.. math::
	:nowrap:

	\begin{align*}
		\delta^3(\pbf_3 - \pbf_1 - \pbf_4) \delta^3(\pbf_5 - \pbf_2 - \pbf_3) \delta^3(\pbf_6 + \pbf_4 - \pbf_5) \
			&= \delta^3(\pbf_5 - \pbf_1 - \pbf_2 - \pbf_4) \delta^3(\pbf_6 + \pbf_4 - \pbf_5) \\
			&= \delta^3(\pbf_6 - \pbf_1 - \pbf_2)
	\end{align*}

where we've first eliminated :math:`\pbf_3` using the first delta function, and then :math:`\pbf_5` using either one of the two remaining delta functions, to arrive at the final result.

.. _sec_quantum_fields_and_antiparticles:

Quantum Fields and Antiparticles
--------------------------------

In this chapter we will construct the Hamiltonians in the form of :math:`H = H_0 + V`, where :math:`H_0` is the Hamiltonian of free particles, and :math:`V = \int d^3 x~\Hscr(0, \xbf)` is a (small) interaction term in the form of :math:`\eqref{eq_defn_v_by_density}`, and the interaction density :math:`\Hscr(x)` is a Lorentz scalar in the sense of :math:`\eqref{eq_h_density_is_scalar}` and satisfies the cluster decomposition principle :math:`\eqref{eq_h_commutativity_for_space_like_separations}`. As a byproduct of the construction, we'll also demystify the so-called *antiparticles* which have been mentioned a number of times so far without definition.


Symmetries and quantum fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following the discussions in :ref:`sec_the_cluster_decomposition_principle`, we'll construct :math:`\Hscr(x)` out of creation and annihilation operators. However, as we've seen in :ref:`sec_the_lorentz_and_cpt_transformation_laws`, the Lorentz transformation of both :math:`a^{\dagger}(q)` and :math:`a(q)` involve in the coefficients a matrix element :math:`D_{\sigma \sigma'}(W(\Lambda, p))` that depend on the momenta of the particles, and hence are not scalars. The idea, then, is to construct :math:`\Hscr(x)` out of the so-called creation and annihilation fields defined as follows

.. math::
	:nowrap:

	\begin{align}
		\psi_{\ell}^+(x) &\coloneqq \sum_{\sigma, n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) a(\pbf, \sigma, n) \label{eq_defn_annihilation_field} \\
		\psi_{\ell}^-(x) &\coloneqq \sum_{\sigma, n} \int d^3 p~v_{\ell}(x;~\pbf, \sigma, n) a^{\dagger}(\pbf, \sigma, n) \label{eq_defn_creation_field}
	\end{align}

where :math:`\ell` is reserved for labeling particles later. We see, in particular, that the creation field :math:`\psi_{\ell}^-(x)` is a superposition of creation operators. Applying it to the vacuum state and let :math:`x` wander around the whole spacetime then creates a (quantum) field.

.. note::
	Just like particles, fields may be either bosonic or fermionic, but not mixed. For example :math:`\psi^-_{\ell}(x)` is a bosonic/fermionic if and only if all the particles created by :math:`a^{\dagger}` in :math:`\eqref{eq_defn_creation_field}` are bosonic/fermionic, respectively.

Lorentz symmetries
++++++++++++++++++

Now the hope is that the creation and annihilation fields transform, under (proper orthochronous) Lorentz symmetry, by a matrix that is independent of the spacetime coordinate :math:`x`. More precisely, we'd like the following to hold

.. math::
	:nowrap:

	\begin{align}
		U_0(\Lambda, b) \psi_{\ell}^+(x) U_0^{-1}(\Lambda, b) &= \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}^+ (\Lambda x + b) \label{eq_conjugate_annihilation_field} \\
		U_0(\Lambda, b) \psi_{\ell}^-(x) U_0^{-1}(\Lambda, b) &= \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi_{\ell'}^- (\Lambda x + b) \label{eq_conjugate_creation_field}
	\end{align}

Note that we've put :math:`\Lambda^{-1}` inside :math:`D_{\ell \ell'}` so that :math:`D` furnishes a representation of the homogeneous Lorentz transformations in the sense that :math:`D(\Lambda_1) D(\Lambda_2) = D(\Lambda_1 \Lambda_2)`. [#two_ways_of_representation]_ There is a priori no reason to use the same representation :math:`D` for both :math:`\psi^{\pm}_{\ell}`, but this turns out to be possible just by calculation. Moreover, the representation :math:`D` is not assumed to be irreducible. Indeed, as we'll see later, it generally decomposes into blocks fixed by further labels.

Then we can try to construct :math:`\Hscr(x)` by a formula similar to :math:`\eqref{eq_general_expansion_of_hamiltonian}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\Hscr(x) = \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \psi_{\ell'_1}^-(x) \cdots \psi_{\ell'_N}^-(x) \psi_{\ell_1}^+(x) \cdots \psi_{\ell_M}^+(x)
		\label{eq_construct_interaction_density_by_fields}
	\end{equation}

It follows from :math:`\eqref{eq_conjugate_annihilation_field}` and :math:`\eqref{eq_conjugate_creation_field}` that the interaction density defined by :math:`\eqref{eq_construct_interaction_density_by_fields}` is a scalar if the following holds

.. math::
	:nowrap:

	\begin{equation}
		g_{\bar{\ell}'_1, \cdots, \bar{\ell}'_N;~\bar{\ell}_1, \cdots, \bar{\ell}_M} = \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} D_{\ell'_1 \bar{\ell}'_1}(\Lambda^{-1}) \cdots D_{\ell'_N \bar{\ell}'_N}(\Lambda^{-1}) D_{\ell_1 \bar{\ell}_1}(\Lambda^{-1}) \cdots D_{\ell_M \bar{\ell}_M}(\Lambda^{-1}) g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M}
		\label{eq_coefficient_g_transformation_law}
	\end{equation}

The solution to the last problem relies on a classification of the representations of the Lorentz group, which has been discussed (in terms of the little group representations) in :ref:`sec_one_particle_states`, and shall be dealt with at a later point. The main goal of this section is to pin down the conditions :math:`u_{\ell}` and :math:`v_{\ell}` must satisfy so that they can be explicitly solved in the following sections.

For this section, we'll focus on the massive case. Recall the Lorentz transformation laws for the creation and annihilation of massive particles :math:`\eqref{eq_lorentz_transformation_formula_for_creation_operator}, \eqref{eq_lorentz_transformation_formula_for_annihilation_operator}` as follows

.. math::
	:nowrap:

	\begin{align}
		U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) &= e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \label{eq_lorentz_transformation_formula_for_annihilation_operator_revisited} \\
		U_0(\Lambda, b) a^{\dagger}(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) &= e^{-\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n) \ast}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a^{\dagger}(\pbf_{\Lambda}, \sigma', n) \label{eq_lorentz_transformation_formula_for_creation_operator_revisited}
	\end{align}

where we've used the fact that :math:`D` is unitary in the sense that :math:`D^{\dagger} = D^{-1}` to invert :math:`W(\Lambda, p)` (and flip the indexes) for later convenience -- mostly because of the use of :math:`\Lambda^{-1}` in :math:`\eqref{eq_conjugate_annihilation_field}` and :math:`\eqref{eq_conjugate_creation_field}`.

Applying :math:`\eqref{eq_lorentz_transformation_formula_for_annihilation_operator_revisited}` to :math:`\eqref{eq_defn_annihilation_field}`, respectively, we get

.. math::
	:nowrap:

	\begin{align*}
		U_0(\Lambda, b) \psi_{\ell}^+(x) U_0^{-1}(\Lambda, b) &= \sum_{\sigma, n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) \\
			&= \sum_{\sigma, \sigma', n} \int d^3 p~u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{(j_n)}_{\sigma \sigma'} (W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \\
			&= \sum_{\sigma, \sigma', n} \int d^3 (\Lambda p)~u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{p_0}{(\Lambda p)_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p)) a(\pbf_{\Lambda}, \sigma', n) \\
			&= \sum_{\sigma', n} \int d^3 (\Lambda p) \blue{\sum_{\sigma} u_{\ell}(x;~\pbf, \sigma, n) e^{\ifrak b \cdot \Lambda p} \sqrt{\frac{p_0}{(\Lambda p)_0}} D^{(j_n)}_{\sigma \sigma'}(W^{-1}(\Lambda, p))} a(\pbf_{\Lambda}, \sigma', n)
	\end{align*}

where in the last equality we've used that the fact :math:`\eqref{eq_lorentz_invariant_3_momentum_volume_element}` that :math:`d^3 p / p_0` is Lorentz invariant. Comparing this with the right-hand-side of :math:`\eqref{eq_conjugate_annihilation_field}`

.. math::
	:nowrap:

	\begin{align*}
		\sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \psi^+_{\ell'}(\Lambda x + b) &= \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \sum_{\sigma, n} \int d^3 p~u_{\ell'}(\Lambda x + b;~\pbf, \sigma, n) a(\pbf, \sigma, n) \\
			&= \sum_{\sigma', n} \sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) \int d^3 (\Lambda p)~u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) a(\pbf_{\Lambda}, \sigma', n) \\
			&= \sum_{\sigma', n} \int d^3 (\Lambda p)~\blue{\sum_{\ell'} D_{\ell \ell'}(\Lambda^{-1}) u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n)} a(\pbf_{\Lambda}, \sigma', n)
	\end{align*}

Equating the blue parts the two calculations, and inverting :math:`D^{(j_n)}_{\sigma \sigma'} (W^{-1}(\Lambda, p))` and :math:`D_{\ell \ell'}(\Lambda^{-1})`, we get

.. math::
	:nowrap:

	\begin{equation}
		\sqrt{\frac{p_0}{(\Lambda p)_0}} e^{\ifrak b \cdot \Lambda p} \sum_{\ell} D_{\ell' \ell}(\Lambda) u_{\ell}(x;~\pbf, \sigma, n) \
			= \sum_{\sigma'} u_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) D_{\sigma' \sigma}^{(j_n)} (W(\Lambda, p))
		\label{eq_annihilation_u_transformation}
	\end{equation}

A parallel calculation for the :math:`v_{\ell}` in :math:`\eqref{eq_conjugate_creation_field}`, which we'll omit, gives the following

.. math::
	:nowrap:

	\begin{equation}
		\sqrt{\frac{p_0}{(\Lambda p)_0}} e^{-\ifrak b \cdot \Lambda p} \sum_{\ell} D_{\ell' \ell}(\Lambda) v_{\ell}(x;~\pbf, \sigma, n) \
			= \sum_{\sigma'} v_{\ell'}(\Lambda x + b;~\pbf_{\Lambda}, \sigma', n) D^{(j_n) \ast}_{\sigma' \sigma}(W(\Lambda, p))
		\label{eq_creation_v_transformation}
	\end{equation}

The identities :math:`\eqref{eq_annihilation_u_transformation}` and :math:`\eqref{eq_creation_v_transformation}` pose the fundamental conditions on :math:`u_{\ell}` and :math:`v_{\ell}`, respectively, which we'll utilize to eventually solve for their solutions. Currently both :math:`u_{\ell}` and :math:`v_{\ell}` depend on :math:`x, \pbf, \sigma` and :math:`n`, and the goal is to use the Lorentz symmetry to reduce the dependencies. This will be carried out in three steps, corresponding to the three types of Lorentz symmetries: translations, boosts, and rotations, as follows.


Translations
	Taking :math:`\Lambda = 1` in :math:`\eqref{eq_creation_v_transformation}` gives :math:`\exp(\ifrak b \cdot p) u_{\ell}(x;~\pbf, \sigma, n) = u_{\ell}(x + b;~\pbf, \sigma, n)`, which then implies

	.. math::
		:nowrap:

		\begin{equation}
			u_{\ell}(x;~\pbf, \sigma, n) = (2\pi)^{-3/2} e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n)
			\label{eq_redefine_u_after_translation}
		\end{equation}

	and similarly

	.. math::
		:nowrap:

		\begin{equation}
			v_{\ell}(x;~\pbf, \sigma, n) = (2\pi)^{-3/2} e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n)
			\label{eq_redefine_v_after_translation}
		\end{equation}

	where we've slightly abused notations by keeping the names of :math:`u_{\ell}` and :math:`v_{\ell}`, while changing their arguments. Here the seemingly redundant :math:`(2\pi)^{-3/2}` is inserted so that the fields

	.. math::
		:nowrap:

		\begin{align}
			\psi^+_{\ell}(x) &= \sum_{\sigma, n} (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n) a(\pbf, \sigma, n) \label{eq_annihilation_field_simplified_by_translation} \\
			\psi^-_{\ell}(x) &= \sum_{\sigma, n} (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n) a^{\dagger}(\pbf, \sigma, n) \label{eq_creation_field_simplified_by_translation}
		\end{align}

	look like the usual Fourier transforms.

	Plugging :math:`\eqref{eq_redefine_u_after_translation}` and :math:`\eqref{eq_redefine_v_after_translation}` into :math:`\eqref{eq_annihilation_u_transformation}` and :math:`\eqref{eq_creation_v_transformation}`, respectively, they can be simplified as follows

	.. math::
		:nowrap:

		\begin{align}
			\sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) u_{\ell}(\pbf, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(\pbf_{\Lambda}, \sigma', n) D^{(j_n)}_{\sigma' \sigma}(W(\Lambda, p)) \label{eq_annihilation_u_transformation_simplified_by_translation} \\
			\sqrt{\frac{p_0}{(\Lambda p)_0}} \sum_{\ell} D_{\ell' \ell}(\Lambda) v_{\ell}(\pbf, \sigma, n) &= \sum_{\sigma'} v_{\ell'}(\pbf_{\Lambda}, \sigma', n) D^{(j_n) \ast}_{\sigma' \sigma}(W(\Lambda, p)) \label{eq_creation_v_transformation_simplified_by_translation}
		\end{align}

	for any homogeneous Lorentz transformation :math:`\Lambda`.

Boosts
	Taking :math:`\pbf = 0` and :math:`\Lambda = L(q)` which takes a particle at rest to one with (arbitrary) momentum :math:`q`, we see that

	.. math::
		:nowrap:

		\begin{equation*}
			W(\Lambda, p) \xlongequal{\eqref{eq_w_from_l}} L(\Lambda p)^{-1} \Lambda L(p) = L(q)^{-1} L(q) = 1
		\end{equation*}

	In this case :math:`\eqref{eq_annihilation_u_transformation_simplified_by_translation}` and :math:`\eqref{eq_creation_v_transformation_simplified_by_translation}` take the following form (with :math:`\qbf` substituted by :math:`\pbf`)

	.. math::
		:nowrap:

		\begin{align}
			\sqrt{\frac{m}{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) u_{\ell}(0, \sigma, n) &= u_{\ell'}(\pbf, \sigma, n)
			\label{eq_annihilation_u_transformation_simplified_by_boost} \\
			\sqrt{\frac{m}{p_0}} \sum_{\ell} D_{\ell' \ell}(L(p)) v_{\ell}(0, \sigma, n) &= v_{\ell'}(\pbf, \sigma, n)
			\label{eq_creation_v_transformation_simplified_by_boost}
		\end{align}

	It follows that one can calculate :math:`u_{\ell}(\pbf, \sigma, n)` for any :math:`\pbf` from the special case of :math:`\pbf = 0` given a representation :math:`D`.

Rotations
	Taking :math:`\pbf = 0` and :math:`\Lambda = \Rcal` a :math:`3`-rotation, and recalling from :math:`\eqref{eq_little_group_rotation}` that :math:`W(\Lambda, p) = \Rcal`, we get special cases of :math:`\eqref{eq_annihilation_u_transformation_simplified_by_translation}` and :math:`\eqref{eq_creation_v_transformation_simplified_by_translation}` as follows

	.. math::
		:nowrap:

		\begin{align*}
			\sum_{\ell} D_{\ell' \ell}(\Rcal) u_{\ell}(0, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(0, \sigma', n) D_{\sigma' \sigma}^{(j_n)}(\Rcal) \\
			\sum_{\ell} D_{\ell' \ell}(\Rcal) v_{\ell}(0, \sigma, n) &= \sum_{\sigma'} v_{\ell'}(0, \sigma', n) D_{\sigma' \sigma}^{(j_n) \ast}(\Rcal)
		\end{align*}

	Using :math:`\eqref{eq_representation_rotation_first_order}` we can further reduce it to the first order as follows

	.. math::
		:nowrap:

		\begin{align}
			\sum_{\ell} \hat{\Jbf}_{\ell' \ell} u_{\ell}(0, \sigma, n) &= \sum_{\sigma'} u_{\ell'}(0, \sigma', n) \Jbf^{(j_n)}_{\sigma' \sigma} \label{eq_j_intertwines_u} \\
			\sum_{\ell} \hat{\Jbf}_{\ell' \ell} v_{\ell}(0, \sigma, n) &= -\sum_{\sigma'} v_{\ell'}(0, \sigma', n) \Jbf^{(j_n) \ast}_{\sigma' \sigma} \label{eq_j_intertwines_v}
		\end{align}

	where :math:`\hat{\Jbf}` denotes the angular momentum vector for the representation :math:`D_{\ell' \ell}(\Rcal)`, in analogy with the usual angular momentum :math:`\Jbf^{(\jfrak)}` for :math:`D^{(j)}(\Rcal)`.


The cluster decomposition principle
+++++++++++++++++++++++++++++++++++

Let's verify that the fields defined by :math:`\eqref{eq_annihilation_field_simplified_by_translation}` and :math:`\eqref{eq_creation_field_simplified_by_translation}`, when plugged into :math:`\eqref{eq_construct_interaction_density_by_fields}`, indeed satisfy the cluster decomposition principle as discussed in :ref:`sec_the_cluster_decomposition_principle`. It's really just a straightforward but tedious calculation which we spell out as follows

.. math::
	:nowrap:

	\begin{align*}
		V(t) &= \int d^3 x~\Hscr(x) \\
			&= \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \int d^3 x~\psi^-_{\ell'_1}(x) \cdots \psi^-_{\ell'_N}(x) \psi^+_{\ell_1}(x) \cdots \psi^+_{\ell_M}(x) \\
			&= \sum_{N,M=0}^{\infty} \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} \sum_{\sigma'_1, \cdots, \sigma'_N} \sum_{n'_1, \cdots, n'_N} \sum_{\sigma_1, \cdots, \sigma_M} \sum_{n_1, \cdots, n_M} (2\pi)^{-3N/2-3M/2} \\
			&\phantom{=} \times \int d^3 p'_1 \cdots d^3 p'_N d^3 p_1 \cdots d^3 p_M~\exp(\ifrak (E_1 + \cdots E_M - E'_1 - \cdots - E'_N)t) \\
			&\phantom{=} \times \blue{\int d^3 x~\exp(\ifrak (\pbf_1 + \cdots \pbf_M - \pbf'_1 - \cdots - \pbf'_N) \cdot \xbf)} \\
			&\phantom{=} \times v_{\ell'_1}(\pbf'_1, \sigma'_1, n'_1) \cdots v_{\ell'_N}(\pbf'_N, \sigma'_N, n'_N) u_{\ell_1}(\pbf_1, \sigma_1, n_1) \cdots u_{\ell_M}(\pbf_M, \sigma_M, n_M) \\
			&\phantom{=} \times a^{\dagger}(\pbf'_1, \sigma'_1, n'_1) \cdots a^{\dagger}(\pbf'_N, \sigma'_N, n'_N) a(\pbf_1, \sigma_1, n_1) \cdots a(\pbf_M, \sigma_M, n_M) \\
			&= \sum_{N,M=0}^{\infty} \sum_{\sigma'_1, \cdots, \sigma'_N} \sum_{n'_1, \cdots, n'_N} \sum_{\sigma_1, \cdots, \sigma_M} \sum_{n_1, \cdots, n_M} \int d^3 p'_1 \cdots d^3 p'_N d^3 p_1 \cdots d^3 p_M \\
			&\phantom{=} \times (2\pi)^{3 - 3N/2 - 3M/2} \exp(\ifrak (E_1 + \cdots E_M - E'_1 - \cdots - E'_N)t) \\
			&\phantom{=} \times a^{\dagger}(\pbf'_1, \sigma'_1, n'_1) \cdots a^{\dagger}(\pbf'_N, \sigma'_N, n'_N) a(\pbf_1, \sigma_1, n_1) \cdots a(\pbf_M, \sigma_M, n_M) \blue{\delta^3(\pbf_1 + \cdots + \pbf_M - \pbf'_1 - \cdots - \pbf'_N)} \\
			&\phantom{=} \times \left( \sum_{\ell'_1, \cdots, \ell'_N} \sum_{\ell_1, \cdots, \ell_M} g_{\ell'_1, \cdots, \ell'_N;~\ell_1, \cdots, \ell_M} v_{\ell'_1}(\pbf'_1, \sigma'_1, n'_1) \cdots v_{\ell'_N}(\pbf'_N, \sigma'_N, n'_N) u_{\ell_1}(\pbf_1, \sigma_1, n_1) \cdots u_{\ell_M}(\pbf_M, \sigma_M, n_M) \right)
	\end{align*}

Besides re-ordering the terms, the only actual calculation is highlighted in the two blue terms, where the second one is the integral of the first. One can compare this calculation with :math:`\eqref{eq_general_expansion_of_hamiltonian}` and see that the cluster decomposition principle is indeed satisfied because there is a unique momentum conservation delta function in each coefficient, as long as :math:`g, u, v` are reasonably smooth, i.e., it's fine to have poles and/or branching singularities but no delta functions.


.. _sec_causality_and_antiparticles:

Causality and antiparticles
+++++++++++++++++++++++++++

We now turn to the other crucial condition on the Hamiltonian, namely, the causality condition :math:`\eqref{eq_h_commutativity_for_space_like_separations}`. Given the general formula :math:`\eqref{eq_construct_interaction_density_by_fields}` of the interaction density, we are forced to require that :math:`[\psi^+_{\ell}(x), \psi^-_{\ell'}(y)] = 0` whenever :math:`x - y` is space-like. However, according to :math:`\eqref{eq_annihilation_field_simplified_by_translation}` and :math:`\eqref{eq_creation_field_simplified_by_translation}`, we have

.. math::
	:nowrap:

	\begin{align*}
		[\psi^+_{\ell}(x), \psi^-_{\ell'}(y)]_{\pm} &= \sum_{\sigma, n} (2\pi)^{-3} \int d^3 p~d^3 p'~e^{\ifrak (p \cdot x - p' \cdot y)} u_{\ell}(\pbf, \sigma, n) v_{\ell'}(\pbf', \sigma, n) [a(\pbf, \sigma, n), a^{\dagger}(\pbf', \sigma, n)]_{\pm} \\
			&= \sum_{\sigma, n} (2\pi)^{-3} \int d^3 p~e^{\ifrak p \cdot (x - y)} u_{\ell}(\pbf, \sigma, n) v_{\ell'}(\pbf, \sigma, n)
	\end{align*}

where the sign :math:`\pm` is positive if the field is fermionic, and negative otherwise. This quantity is not necessarily vanishing even if :math:`x - y` is space-like.

Therefore in order to construct :math:`\Hscr` in the form of :math:`\eqref{eq_construct_interaction_density_by_fields}` that satisfies :math:`\eqref{eq_h_commutativity_for_space_like_separations}`, we must not just use :math:`\psi^{\pm}(x)` as the building blocks. It turns out that one may consider a linear combination of the two as follows

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\ell}(x) \coloneqq \kappa_{\ell} \psi^+_{\ell}(x) + \lambda_{\ell} \psi^-_{\ell}(x)
		\label{eq_defn_psi_field}
	\end{equation}

as well as its adjoint :math:`\psi^{\dagger}_{\ell}(x)`, and hope that they satisfy

.. math::
	:nowrap:

	\begin{equation}
		[\psi_{\ell}(x), \psi_{\ell}(y)]_{\pm} = [\psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y)]_{\pm} = 0
		\label{eq_space_like_commutativity_for_combined_field}
	\end{equation}

whenever :math:`x-y` is space-like, and replace :math:`\psi^{\pm}_{\ell}(x)` with :math:`\psi_{\ell}(x), \psi^{\dagger}_{\ell}(x)` in :math:`\eqref{eq_construct_interaction_density_by_fields}`. Under these assumptions, we can then construct the interaction density :math:`\Hscr` as a polynomial in :math:`\psi_{\ell}(x), \psi_{\ell}^{\dagger}(x)` with an even number of fermionic fields (so that the sign in :math:`\eqref{eq_space_like_commutativity_for_combined_field}` is negative).

There remains, however, one issue with field like :math:`\eqref{eq_space_like_commutativity_for_combined_field}` that mixes creation and annihilation fields. Namely, special conditions must hold in order for such fields to play well with conserved quantum numbers. To be more specific, let :math:`Q` be a conserved quantum number, e.g., the electric charge. Then the following hold

.. math::
	:nowrap:

	\begin{align}
		[Q, a(\pbf, \sigma, n)] &= -q(n) a(\pbf, \sigma, n) \label{eq_charge_of_creation_operator} \\
		[Q, a^{\dagger}(\pbf, \sigma, n)] &= q(n) a^{\dagger}(\pbf, \sigma, n) \label{eq_charge_of_annihilation_operator}
	\end{align}

where :math:`q(n)` denotes the quantum number of the particle species :math:`n`. These identities can be verified by applying both sides to :math:`\Psi_{\pbf, \sigma, n}` and :math:`\Psi_{\VAC}`, respectively.

Now in order for :math:`Q` to commute with :math:`\Hscr`, which is constructed as a polynomial of :math:`\psi_{\ell}(x)` and :math:`\psi_{\ell}^{\dagger}(x)`, we better have

.. math::
	:nowrap:

	\begin{equation}
		[Q, \psi_{\ell}(x)] = -q_{\ell} \psi_{\ell}(x)
		\label{eq_charge_of_psi_field}
	\end{equation}

so that each monomial (with coefficient neglected) :math:`\psi^{\dagger}_{\ell'_1}(x) \cdots \psi^{\dagger}_{\ell'_M}(x) \psi_{\ell_1}(x) \cdots \psi_{\ell_N}(x)` in :math:`\Hscr` will commute with :math:`Q` if

.. math::
	:nowrap:

	\begin{equation*}
		q_{\ell_1} + \cdots + q_{\ell_N} = q_{\ell'_1} + \cdots + q_{\ell'_M}
	\end{equation*}

Note that the negative sign in :math:`\eqref{eq_charge_of_psi_field}` is a formal analogy to :math:`\eqref{eq_charge_of_creation_operator}`, where we think of :math:`\psi_{\ell}(x)` as an annihilation field even though it's really not. Since :math:`\psi_{\ell}(x)` is a linear combination of :math:`\psi^+_{\ell}(x)` and :math:`\psi^-_{\ell}(x)`, which in turn are superpositions of annihilation and creation operators, respectively, it follows from :math:`\eqref{eq_charge_of_creation_operator}` and :math:`\eqref{eq_charge_of_annihilation_operator}` that in order for :math:`\eqref{eq_charge_of_psi_field}` to hold, the following conditions must be satisfied

1. all particles annihilated by :math:`\psi^+_{\ell}(x)` must have the same charge :math:`q(n) = q_{\ell}`,
2. all particles created by :math:`\psi^-_{\ell}(x)` must have the same charge :math:`q(n) = -q_{\ell}`, and
3. for any particle of species :math:`n`, which is annihilated by :math:`\psi^+_{\ell}(x)`, there exists a particle of species :math:`\bar{n}`, which is created by :math:`\psi^-_{\ell}(x)`, such that :math:`q(n) = -q(\bar{n})`.

The particles of species :math:`n` and :math:`\bar{n}` are called *antiparticles* of each other -- they are exactly the same except for the charges which are opposite. It is the last condition that demands the existence of particle-antiparticle pairs so that one can formulate a consistent (relativistic) quantum field theory.

.. dropdown:: The Klein-Gordon equation
	:icon: unlock
	:animate: fade-in-slide-down

	It follows from the definition :math:`\eqref{eq_defn_psi_field}`, together with :math:`\eqref{eq_annihilation_field_simplified_by_translation}` and :math:`\eqref{eq_creation_field_simplified_by_translation}`, that the field :math:`\psi_{\ell}(x)` satisfies the following so-called `Klein-Gordon equation <https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation>`_

	.. math::
		:nowrap:

		\begin{equation}
			\left( \square - m^2 \right) \psi_{\ell}(x) = 0
			\label{eq_klein_gordon}
		\end{equation}

	where :math:`\square \coloneqq \eta^{\mu \nu} \p_{\mu} \p_{\nu}` is the `d'Alembert operator <https://en.wikipedia.org/wiki/D%27Alembert_operator>`_ and :math:`m` is the (definite) mass of the field. This equation is traditionally one of the starting points of quantum field theory, from which creation/annihilation operators can be derived through the so-called canonical quantization formalism. However, we've derived the equation here from the other way around, namely, the creation/annihilation operators, which in turn come from the first principles of quantum mechanics and Lorentz symmetry.


.. _sec_scalar_field:

Scalar fields
^^^^^^^^^^^^^

We'll start, as always, with the simplest case of scalar fields, namely, when :math:`\psi^+(x) = \psi^+_{\ell}(x)` and :math:`\psi^-(x) = \psi^-_{\ell}(x)` are scalar functions. We argue first that such fields can only create/annihilate spinless particles. Indeed, since :math:`\hat{\Jbf}` necessarily vanishes, it follows from :math:`\eqref{eq_j_intertwines_u}` and :math:`\eqref{eq_j_intertwines_v}` that :math:`u` and :math:`v` may be nonzero if and only if :math:`j_n = 0`. If we, for the moment, are concerned with just one particle species, then we can write :math:`u(\pbf, \sigma, n) = u(\pbf)` and :math:`v(\pbf, \sigma, n) = v(\pbf)`. Lastly, we note that since :math:`D = 1` in this case, :math:`\eqref{eq_annihilation_u_transformation_simplified_by_translation}` and :math:`\eqref{eq_creation_v_transformation_simplified_by_translation}` become

.. math::
	:nowrap:

	\begin{align*}
		\sqrt{p_0}~u(\pbf) &= \sqrt{(\Lambda p)_0}~u(\pbf_{\Lambda}) \\
		\sqrt{p_0}~v(\pbf) &= \sqrt{(\Lambda p)_0}~v(\pbf_{\Lambda})
	\end{align*}

It follows that

.. math::
	:nowrap:

	\begin{equation}
		u(\pbf)	= v(\pbf) = (2p_0)^{-1/2}
		\label{eq_scalar_u_and_v}
	\end{equation}

where the factor :math:`2` is purely conventional. In particular :math:`u(0) = v(0) = (2m)^{-1/2}`.

Plugging :math:`\eqref{eq_scalar_u_and_v}` into :math:`\eqref{eq_annihilation_field_simplified_by_translation}, \eqref{eq_creation_field_simplified_by_translation}`, we get

.. math::
	:nowrap:

	\begin{align}
		\psi^+(x) &= \int d^3 p~(2\pi)^{-3/2} e^{\ifrak p \cdot x} (2p_0)^{-1/2} a(\pbf)
		\label{eq_scalar_field_psi_plus} \\
		\psi^-(x) &= \int d^3 p~(2\pi)^{-3/2} e^{-\ifrak p \cdot x} (2p_0)^{-1/2} a^{\dagger}(\pbf) = \psi^{+ \dagger}(x)
		\label{eq_scalar_field_psi_plus_and_minus_are_adjoints}
	\end{align}

To construct :math:`\Hscr`, we first note that :math:`\eqref{eq_coefficient_g_transformation_law}` trivially holds in this case because :math:`D = 1` and :math:`g` is simply a scalar. Next let's consider the causality condition which demands that :math:`\left[ \psi^+(x), \psi^-(y) \right] = 0` whenever :math:`x - y` is space-like. Using the canonical commutation relation :math:`\eqref{eq_creation_annihilation_commutator}` we calculate

.. math::
	:nowrap:

	\begin{align}
		\left[ \psi^+(x), \psi^-(y) \right]_{\pm} &= \int d^3 p~d^3 q~(2\pi)^{-3} e^{\ifrak (p \cdot x - q \cdot y)} (4 p_0 q_0)^{-1/2} \left[ a(\pbf), a^{\dagger}(\qbf) \right]_{\pm}
		\label{eq_scalar_field_commutator_as_Delta} \\
			&= \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0}~e^{\ifrak p \cdot (x - y)} \eqqcolon \Delta_+(x - y) \nonumber
	\end{align}

where

.. math::
	:nowrap:

	\begin{equation}
		\Delta_+(x) \coloneqq \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0}~e^{\ifrak p \cdot x}
		\label{eq_defn_Delta_plus}
	\end{equation}

We notice that :math:`\Delta_+(x)` is manifestly (proper orthochronous) Lorentz invariant -- the invariance of the volume element comes from :math:`\eqref{eq_lorentz_invariant_3_momentum_volume_element}`. It is, however, not in general invariant under transformations like :math:`x \to -x`. But, as we'll see, such invariance holds assuming :math:`x` is space-like.

.. note::
	The plus subscript in :math:`\Delta_+(x)` is there to distinguish it from an anti-symmetrized version :math:`\Delta(x)` to be introduced later.

Now we'll restrict ourselves to the special case of a space-like :math:`x` which, up to a Lorentz transformation, can be assumed to take the form :math:`x = (0, \xbf)` with :math:`|\xbf| > 0`. In this case, we can then calculate :math:`\Delta_+(x)` as follows [#wrong_integration_of_Delta_function]_

.. math::
	:nowrap:

	\begin{align*}
		\Delta_+(x) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2\sqrt{\pbf^2 + m^2}}~\exp(\ifrak \pbf \cdot \xbf) \\
			&= \frac{4\pi}{(2\pi)^3} \int_0^{\infty} \frac{\pbf^2 d|\pbf|}{2\sqrt{\pbf^2 + m^2}} \int_{S^2} d^2 \hat{\pbf}~\exp(\ifrak |\pbf| |\xbf| \hat{\pbf} \cdot \hat{\xbf}) \\
			&= \frac{1}{2\pi} \int_0^{\infty} \frac{\pbf^2 d|\pbf|}{2\sqrt{\pbf^2 + m^2}} \int_0^{\pi} d\theta~\exp(\ifrak |\pbf| |\xbf| \cos\theta)
	\end{align*}

The last integral cannot be easily evaluated, at least without some knowledge about special functions. Nonetheless, we observe that :math:`\Delta_+(x) \neq 0`, which means that :math:`\Hscr` cannot be just any polynomial in :math:`\psi^{\pm}(x)`. Moreover, we note that :math:`\Delta_+(x) = \Delta_+(-x)` as promised earlier. As already mentioned in :math:`\eqref{eq_defn_psi_field}`, let's try

.. math::
	:nowrap:

	\begin{equation}
		\psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^-(x)
		\label{eq_scalar_field_first_defn_of_psi}
	\end{equation}

We can then try to make :math:`\eqref{eq_space_like_commutativity_for_combined_field}` hold by the following calculations

.. math::
	:nowrap:

	\begin{align*}
		\left[ \psi(x), \psi(y) \right]_{\pm} &\xlongequal{\phantom{(000)}} \kappa\lambda \left(\left[ \psi^+(x), \psi^-(y) \right]_{\pm} + \left[ \psi^-(x), \psi^+(y) \right]_{\pm} \right) \\
			&\xlongequal{\eqref{eq_scalar_field_commutator_as_Delta}} \kappa\lambda (1 \pm 1) \Delta(x - y) \\
		\left[ \psi(x), \psi^{\dagger}(y) \right]_{\pm} &\xlongequal{\eqref{eq_scalar_field_psi_plus_and_minus_are_adjoints}} \left[ \kappa \psi^+(x) + \lambda \psi^-(x), \kappa^{\ast} \psi^-(y) + \lambda^{\ast} \psi^+(y) \right]_{\pm} \\
			&\xlongequal{\phantom{(000)}} |\kappa|^2 \left[ \psi^+(x), \psi^-(y) \right]_{\pm} + |\lambda|^2 \left[ \psi^-(x), \psi^+(y) \right]_{\pm} \\
			&\xlongequal{\eqref{eq_scalar_field_commutator_as_Delta}} \left( |\kappa|^2 \pm |\lambda|^2 \right) \Delta(x - y)
	\end{align*}

We see that :math:`\eqref{eq_space_like_commutativity_for_combined_field}` holds for scalar fields if the fields are bosonic, i.e., the bottom sign in :math:`\pm` applies, and :math:`|\kappa| = |\lambda|`. By adjust the phase of :math:`a(\pbf)`, we can actually arrange so that :math:`\kappa = \lambda`, in which case we have

.. math::
	:nowrap:

	\begin{equation}
		\psi(x) = \psi^+(x) + \psi^-(x) = \psi^+(x) + \psi^{+ \dagger}(x) = \psi^{\dagger}(x)
		\label{eq_scalar_field_psi_fixed_phase}
	\end{equation}

.. note::
	Although the arrangement of phase so that :math:`\kappa = \lambda` is a mere convention, it's a convention that needs to be applied to *all* scalar fields appearing in :math:`\Hscr`. Namely, one cannot have both :math:`\psi(x)` as in :math:`\eqref{eq_scalar_field_psi_fixed_phase}` and another

	.. math::
		:nowrap:

		\begin{equation*}
			\psi'(x) = e^{\ifrak \theta} \psi^+(x) + e^{-\ifrak \theta} \psi^{+ \dagger}(x)
		\end{equation*}

	for some :math:`\theta`, because :math:`\psi(x)` won't commute with :math:`\psi'(y)` even if :math:`x - y` is space-like.

Now if the particle created and annihilated by :math:`\psi(x)` carries a (non-vanishing) conserved quantum number :math:`Q`, then by the discussions on the charge conservation from the previous section, a density :math:`\Hscr` made up of :math:`\psi(x)` as defined by :math:`\eqref{eq_scalar_field_psi_fixed_phase}` will not commute with :math:`Q`. Instead, one must assume the existence of a field :math:`\psi^{+ c}(x)` that creates and annihilates the corresponding antiparticle, in the sense that

.. math::
	:nowrap:

	\begin{align*}
		\left[ Q, \psi^+(x) \right] &= -q \psi^+(x) \\
		\left[ Q, \psi^{+ c}(x) \right] &= q \psi^{+ c}(x)
	\end{align*}

Here the supscript :math:`c` stands for charge (conjugation). Now instead of :math:`\eqref{eq_scalar_field_first_defn_of_psi}`, let's try

.. math::
	:nowrap:

	\begin{equation*}
		\psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x)
	\end{equation*}

so that :math:`[Q, \psi(x)] = -q \psi(x)`. We calculate the commutators, assuming the antiparticle is different from the particle, just as before as follows

.. math::
	:nowrap:

	\begin{align*}
		\left[\psi(x), \psi(y) \right]_{\pm} &= \left[ \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x), \kappa \psi^+(y) + \lambda \psi^{+ c \dagger}(y) \right]_{\pm} = 0 \\
		\left[\psi(x), \psi^{\dagger}(y) \right]_{\pm} &= \left[ \kappa \psi^+(x) + \lambda \psi^{+ c \dagger}(x), \kappa^{\ast} \psi^{+ \dagger}(y) + \lambda^{\ast} \psi^{+ c}(y) \right]_{\pm} \\
			&= |\kappa|^2 \left[ \psi^+(x), \psi^{+ \dagger}(y) \right]_{\pm} + |\lambda|^2 \left[ \psi^{+ c \dagger}(x), \psi^{+ c}(y) \right]_{\pm} \\
			&= (|\kappa|^2 \pm |\lambda|^2) \Delta(x - y)
	\end{align*}

where we've assumed, in particular that the particle and its particle share the same mass so that :math:`\eqref{eq_scalar_field_commutator_as_Delta}` equally applies.

By the same argument as in the case where no quantum number is involved, we see that a scalar field can satisfy the causality condition if it describes a boson. Moreover, by adjusting the phase of :math:`a(\pbf)`, one can arrange so that :math:`\kappa = \lambda` so that

.. math::
	:nowrap:

	\begin{equation}
		\psi(x) = \psi^+(x) + \psi^{+ c \dagger}(x)
		\label{eq_scalar_field_psi_fixed_phase_with_antiparticle}
	\end{equation}

Note that this is compatible with :math:`\eqref{eq_scalar_field_psi_fixed_phase}` in the case where the particle is its own antiparticle.

Using :math:`\eqref{eq_scalar_field_psi_plus}` and :math:`\eqref{eq_scalar_field_psi_plus_and_minus_are_adjoints}`, we can write :math:`\psi(x)` in terms of the creation and annihilation operators as follows

.. math::
	:nowrap:

	\begin{equation}
		\psi(x) = \int \frac{d^3 p}{(2\pi)^{3/2} (2p_0)^{1/2}}~\left[ e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{c \dagger}(\pbf) \right]
		\label{eq_scalar_field_psi_by_creation_and_annihilation_operators}
	\end{equation}

with the possibility of :math:`a^{c \dagger}(\pbf) = a^{\dagger}(\pbf)` in the case where the created particle is its own antiparticle.

For later use (e.g., the evaluation of Feynman diagrams), we note the following identity which holds for any, and not just space-like, :math:`x` and :math:`y`.

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi(x), \psi^{\dagger}(y) \right] = \Delta(x - y)
		\label{eq_scalar_field_commutator}
	\end{equation}

where :math:`\Delta(x)` is defined as follows

.. math::
	:nowrap:

	\begin{equation}
		\Delta(x) \coloneqq \Delta_+(x) - \Delta_+(-x) = \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0} \left( e^{\ifrak p \cdot x} - e^{-\ifrak p \cdot x} \right)
		\label{eq_defn_Delta}
	\end{equation}


The CPT symmetries
++++++++++++++++++

Let's investigate how a scalar field transforms under spatial inversion :math:`\Pcal`, time inversion :math:`\Tcal`, and charge conjugation :math:`\Ccal`. This follows essentially from :math:`\eqref{eq_scalar_field_psi_by_creation_and_annihilation_operators}` together with our knowledge about how creation/annihilation operators transform under CPT transformations in :ref:`sec_the_lorentz_and_cpt_transformation_laws`. Recall that we consider the case of massive particles here, leaving the massless case to a later section.

We start with the spatial inversion :math:`\Pcal` by recalling the following transformation rules

.. math::
	:nowrap:

	\begin{align*}
		U(\Pcal) a(\pbf) U^{-1}(\Pcal) &= \eta^{\ast} a(-\pbf) \\
		U(\Pcal) a^{c \dagger}(\pbf) U^{-1}(\Pcal) &= \eta^c a^{c \dagger}(-\pbf)
	\end{align*}

where :math:`\eta` and :math:`\eta^c` are the intrinsic parities of the particle and antiparticle, respectively. In order for the scalar field :math:`\eqref{eq_scalar_field_psi_by_creation_and_annihilation_operators}` to transform nicely with :math:`\Pcal`, one must have :math:`\eta^{\ast} = \eta^c` (or :math:`\eta^{\ast} = \eta` in the case where the particle is its own antiparticle). As a result, we have

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \psi(x) U^{-1}(\Pcal) = \eta^{\ast} \psi(\Pcal x)
		\label{eq_scalar_field_spatial_inversion_transformation_law}
	\end{equation}

Next let's consider the time inversion :math:`\Tcal`. We recall the transformation rules as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Tcal) a(\pbf) U^{-1}(\Tcal) &= \zeta^{\ast} a(-\pbf) \\
		U(\Tcal) a^{c \dagger}(\pbf) U^{-1}(\Tcal) &= \zeta^c a^{c \dagger}(-\pbf)
	\end{align*}

Similar to the case of spatial inversions, in order for :math:`\psi(x)` to transform nicely with :math:`U(\Tcal)`, one must have :math:`\zeta^{\ast} = \zeta^c`. Moreover, since :math:`U(\Tcal)` is anti-unitary, we have

.. math::
	:nowrap:

	\begin{equation*}
		U(\Tcal) \psi(x) U^{-1}(\Tcal) = \zeta^{\ast} \psi(-\Tcal x)
	\end{equation*}

Finally let's consider the charge conjugation :math:`\Ccal` with the following transformation laws

.. math::
	:nowrap:

	\begin{align*}
		U(\Ccal) a(\pbf) U^{-1}(\Ccal) &= \xi^{\ast} a^c(\pbf) \\
		U(\Ccal) a^{c \dagger}(\pbf) U^{-1}(\Ccal) &= \xi^c a^{\dagger}(\pbf)
	\end{align*}

As before, we must have :math:`\xi^{\ast} = \xi^c` and therefore

.. math::
	:nowrap:

	\begin{equation*}
		U(\Ccal) \psi(x) U^{-1}(\Ccal) = \xi^{\ast} \psi^{\dagger}(x)
	\end{equation*}


Vector fields
^^^^^^^^^^^^^

The next simplest scenario after scalar field is vector field, where the representation :math:`D(\Lambda) = \Lambda`. Once again, let's consider particles of one species so that we can drop the :math:`n` label from, for example, :math:`a(\pbf, \sigma, n)`. In this case, we can rewrite :math:`\eqref{eq_annihilation_field_simplified_by_translation}` and :math:`\eqref{eq_creation_field_simplified_by_translation}` as follows

.. math::
	:nowrap:

	\begin{align}
		\psi^+_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} u_{\mu}(\pbf, \sigma) a(\pbf, \sigma)
		\label{eq_vector_field_psi_plus} \\
		\psi^-_{\nu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} v_{\nu}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma)
		\label{eq_vector_field_psi_minus}
	\end{align}

where :math:`\mu, \nu` are the :math:`4`-indexes. Moreover, the boost transformation formulae :math:`\eqref{eq_annihilation_u_transformation_simplified_by_boost}` and :math:`\eqref{eq_creation_v_transformation_simplified_by_boost}` take the following form

.. math::
	:nowrap:

	\begin{align}
		u_{\mu}(\pbf, \sigma) &= (m / p_0)^{1/2} L(p)_{\mu}^{\nu} u_{\nu}(0, \sigma)
		\label{eq_vector_field_boost_u} \\
		v_{\mu}(\pbf, \sigma) &= (m / p_0)^{1/2} L(p)_{\mu}^{\nu} v_{\nu}(0, \sigma)
		\label{eq_vector_field_boost_v}
	\end{align}

Finally the (linearized) rotation transformation formulae :math:`\eqref{eq_j_intertwines_u}` and :math:`\eqref{eq_j_intertwines_v}` take the following form

.. math::
	:nowrap:

	\begin{align}
		\sum_{\sigma'} u_{\mu}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_{\nu} \hat{\Jbf}_{\mu \nu} u_{\nu}(0, \sigma)
		\label{eq_vector_field_angular_momentum_intertwines_u} \\
		-\sum_{\sigma'} v_{\mu}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma} &= \sum_{\nu} \hat{\Jbf}_{\mu \nu} v_{\nu}(0, \sigma)
		\label{eq_vector_field_angular_momentum_intertwines_v}
	\end{align}

where :math:`\hat{\Jbf}` is the angular momentum vector associated with the (tautological) representation :math:`\Lambda`. It follows from :math:`\eqref{eq_expansion_of_Lambda}` and :math:`\eqref{eq_u_lorentz_expansion}` that

.. math::
	:nowrap:

	\begin{align}
		\left( \hat{\Jbf}_k \right)_{00} = \left( \hat{\Jbf}_k \right)_{0i} = \left( \hat{\Jbf}_k \right)_{i0} &= 0
		\label{eq_vector_field_j_intertwines_u} \\
		\left( \hat{\Jbf}_k \right)_{ij} &= -\ifrak \epsilon_{ijk}
		\label{eq_vector_field_j_intertwines_v}
	\end{align}

where :math:`\{i,j,k\} = \{1,2,3\}`. From this one can then calculate the following squares

.. math::
	:nowrap:

	\begin{align*}
		\left( \hat{\Jbf}^2 \right)_{00} &= \left( \hat{\Jbf}^2 \right)_{0i} = \left( \hat{\Jbf}^2 \right)_{i0} = 0 \\
		\left( \hat{\Jbf}^2 \right)_{ij} &= \sum_{k,m=1}^3 \left( \hat{\Jbf}_k \right)_{im} \left( \hat{\Jbf}_k \right)_{mj} = \sum_{k,m=1}^3 -\epsilon_{imk} \epsilon_{mjk} = 2\delta_{ij}
	\end{align*}

It follows then from :math:`\eqref{eq_vector_field_j_intertwines_u}` and :math:`\eqref{eq_vector_field_j_intertwines_v}` that

.. math::
	:nowrap:

	\begin{align}
		\sum_{\sigma'} u_0(0, \sigma') \left(\Jbf^{(\jfrak)}\right)^2_{\sigma' \sigma} \
			&= \sum_{\sigma' \sigma''} u_0(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma''} \Jbf^{(\jfrak)}_{\sigma'' \sigma} \label{eq_vector_field_u0_intertwines_j_sq} \\
			&= \sum_{\nu \sigma'} \Jbf^{(\jfrak)}_{\sigma' \sigma''} \hat{\Jbf}_{0 \nu} u_{\nu}(0, \sigma'') = 0 \nonumber \\
		\sum_{\sigma'} u_i(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma} \
			&= \sum_{\sigma' \sigma''} u_i(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma''} \Jbf^{(\jfrak)}_{\sigma'' \sigma} \label{eq_vector_field_ui_intertwines_j_sq} \\
			&= \sum_{\nu \sigma''} \hat{\Jbf}_{i \nu} u_{\nu}(0, \sigma'') \Jbf^{(\jfrak)}_{\sigma'' \sigma} \nonumber \\
			&= \sum_{\nu \rho} \hat{\Jbf}_{i \nu} \hat{\Jbf}_{\nu \rho} u_{\rho}(0, \sigma) \nonumber \\
			&= 2 \delta_{i \rho} u_{\rho}(0, \sigma) = 2 u_i(0, \sigma) \nonumber \\
		\sum_{\sigma'} v_0(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma} &= 0
		\label{eq_vector_field_v0_intertwines_j_sq} \\
		\sum_{\sigma'} v_i(0, \sigma') \left( \Jbf^{(\jfrak)} \right)^2_{\sigma' \sigma} &=2v_i(0, \sigma)
		\label{eq_vector_field_vi_intertwines_j_sq}
	\end{align}

where we've worked out the details of the calculations for :math:`u`, but not :math:`v` because they are essentially the same.

The reason for squaring the angular momentum :math:`3`-vector as above is because it has a particularly simple matrix form

.. math::
	:nowrap:

	\begin{equation*}
		\left( \Jbf^{(\jfrak)} \right)^2_{\sigma \sigma'} = \jfrak (\jfrak + 1) \delta_{\sigma \sigma'}
	\end{equation*}

by :math:`\eqref{eq_angular_momentum_squared_eigenvalue}`. It follows that in order for :math:`\eqref{eq_vector_field_u0_intertwines_j_sq}` -- :math:`\eqref{eq_vector_field_vi_intertwines_j_sq}` to have nonzero solutions, one must have either :math:`\jfrak = 0`, in which case only the time-components :math:`u_0(0)` and :math:`v_0(0)` may be nonzero, where we've also suppressed :math:`\sigma` because spin vanishes, or :math:`\jfrak = 1`, in which case only the space-components :math:`u_i(0, \sigma)` and :math:`v_i(0, \sigma)` may be nonzero. These two cases are discussed in more details as follows.

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

It follows from :math:`\eqref{eq_vector_field_boost_u}` and :math:`\eqref{eq_vector_field_boost_v}` (see also :math:`\eqref{eq_L_transformation_for_massive_1}, \eqref{eq_L_transformation_for_massive_2}`) that

.. math::
	:nowrap:

	\begin{align*}
		u_{\mu}(\pbf) &= (m / p_0)^{1/2} L(p)^0_{\mu} u_0(0) \\
			&= (m / p_0)^{1/2} (p_{\mu} / m) \ifrak (m / 2)^{1/2} \\
			&= \ifrak p_{\mu} (2p_0)^{-1/2} \\
		v_{\mu}(\pbf) &= -\ifrak p_{\mu} (2p_0)^{-1/2}
	\end{align*}

where we once again have omitted the details of the calculation of :math:`v` because it's similar to that of :math:`u`. Plugging into :math:`\eqref{eq_vector_field_psi_plus}` and :math:`\eqref{eq_vector_field_psi_minus}`, we see that the field components take the following form

.. math::
	:nowrap:

	\begin{align*}
		\psi^+_{\mu}(x) &= (2\pi)^{-3/2} \int d^3 p~e^{\ifrak p \cdot x} \ifrak p_{\mu} (2p_0)^{-1/2} a(\pbf) \\
		\psi^-_{\mu}(x) &= (2\pi)^{-3/2} \int d^3 p~e^{-\ifrak p \cdot x} (-\ifrak p_{\mu}) (2p_0)^{-1/2} a^{\dagger}(\pbf)
	\end{align*}

Comparing these with :math:`\eqref{eq_scalar_field_psi_plus}` and :math:`\eqref{eq_scalar_field_psi_plus_and_minus_are_adjoints}`, and thanks to the choices of :math:`u_0(0)` and :math:`v_0(0)` above, we see that

.. math::
	:nowrap:

	\begin{equation*}
		\psi^+_{\mu}(x) = \p_{\mu} \psi^+(x), \quad \psi^-_{\mu}(x) = \p_{\mu} \psi^-(x)
	\end{equation*}

It follows that in fact a spinless vector field defined by :math:`\psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^-_{\mu}(x)` as usual is nothing but the gradient vector field of a (spinless) scalar field. Hence we get nothing new from spinless vector fields.

.. _sec_spin_1_vector_fields:

Spin-:math:`1` vector fields
++++++++++++++++++++++++++++

In this case :math:`\jfrak = 1`. We start with the state whose spin :math:`z`-component vanishes, i.e., :math:`u_{\mu}(0,0)` and :math:`v_{\mu}(0,0)`. First we claim that they are both in the :math:`z`-direction, i.e., :math:`u_{\mu}(0,0) = v_{\mu}(0,0) = 0` unless :math:`\mu=3`. Indeed, taking the :math:`z`-components of both sides of :math:`\eqref{eq_vector_field_angular_momentum_intertwines_u}` and recalling that :math:`\left( J_3^{(1)} \right)_{0 \sigma} = 0`, we have for :math:`\mu = 1`

.. math::
	:nowrap:

	\begin{equation*}
		0 = \sum_{\nu} \left( \hat{\Jbf}_3 \right)_{1 \nu} u_{\nu}(0, 0) = -\ifrak u_2(0, 0) \implies u_2(0, 0) = 0
	\end{equation*}

and for :math:`\mu = 2`

.. math::
	:nowrap:

	\begin{equation*}
		0 = \sum_{\nu} \left( \hat{\Jbf}_3 \right)_{2 \nu} u_{\nu}(0, 0) = \ifrak u_1(0, 0) \implies u_1(0, 0) = 0
	\end{equation*}

These, together with :math:`\eqref{eq_vector_field_u0_intertwines_j_sq}`, imply that only :math:`u_3(0, 0)` is nonzero. The same conclusion can also be drawn for :math:`v_3(0, 0)`. Therefore up to a normalization factor, we can write

.. math::
	:nowrap:

	\begin{equation}
		u_{\mu}(0, 0) = v_{\mu}(0, 0) = (2m)^{-1/2} \begin{bmatrix*}[r] 0 \\ 0 \\ 0 \\ 1 \end{bmatrix*}
		\label{eq_vector_field_uv_spin_z_0}
	\end{equation}

Now to calculate :math:`u` and :math:`v` for the other spin :math:`z`-components, we'll try to use :math:`\eqref{eq_j1_j2_matrix}` as follows. According to :math:`\eqref{eq_vector_field_angular_momentum_intertwines_u}` we have

.. math::
	:nowrap:

	\begin{align*}
		\sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{\mu \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{\mu \nu} \right) u_{\nu}(0, \sigma) \
			&= \sum_{\sigma'} u_{\mu}(0, \sigma') \left( J^{(1)}_1 + \ifrak J^{(1)}_2 \right)_{\sigma \sigma'} \\
			&= \sum_{\sigma'} u_{\mu}(0, \sigma') \delta_{\sigma+1, \sigma'} \sqrt{(1 - \sigma)(2 + \sigma)}
	\end{align*}

Letting :math:`\sigma=0` and :math:`\mu=1`, we have

.. math::
	:nowrap:

	\begin{equation*}
		\sqrt{2}~u_1(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{1 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{1 \nu} \right) u_{\nu}(0, 0) = \ifrak \left( \hat{\Jbf}_2 \right)_{13} u_3(0, 0) = -(2m)^{-1/2}
	\end{equation*}

Changing to :math:`\mu=2`, we have

.. math::
	:nowrap:

	\begin{equation*}
		\sqrt{2}~u_2(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{2 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{2 \nu} \right) u_{\nu}(0, 0) = (\hat{\Jbf}_1)_{23} u_3(0, 0) = -\ifrak (2m)^{-1/2}
	\end{equation*}

Finally take :math:`\mu=3`, we have

.. math::
	:nowrap:

	\begin{equation*}
		\sqrt{2}~u_3(0, 1) = \sum_{\nu} \left(\left( \hat{\Jbf}_1 \right)_{3 \nu} + \ifrak \left( \hat{\Jbf}_2 \right)_{3 \nu} \right) u_{\nu}(0, 0) \
			= \left( \hat{\Jbf}_1 \right)_{32} u_2(0, 0) + \ifrak \left( \hat{\Jbf}_2 \right)_{31} u_1(0, 0) = 0
	\end{equation*}

Putting these all together, we have calculated :math:`u_{\mu}(0, 1)` as follows

.. math::
	:nowrap:

	\begin{equation*}
		u_{\mu}(0, 1) = -\frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*}
	\end{equation*}

Calculations for :math:`\sigma = -1` as well as for :math:`v` are similar and hence omitted. The results are listed for future reference as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		u_{\mu}(0, 1) &= -v_{\mu}(0, -1) &&= -\frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*}
		\label{eq_vector_field_uv_spin_z_1} \\
		u_{\mu}(0, -1) &= -v_{\mu}(0, 1) &&= \frac{1}{2\sqrt{m}} \begin{bmatrix*}[r] 0 \\ 1 \\ -\ifrak \\ 0 \end{bmatrix*}
		\label{eq_vector_field_uv_spin_z_2}
	\end{alignat}

Applying the boosting formulae :math:`\eqref{eq_vector_field_boost_u}, \eqref{eq_vector_field_boost_v}` to :math:`\eqref{eq_vector_field_uv_spin_z_0}, \eqref{eq_vector_field_uv_spin_z_1}, \eqref{eq_vector_field_uv_spin_z_2}`, we obtain the formulae for :math:`u, v` at arbitrary momentum as follows

.. math::
	:nowrap:

	\begin{equation}
		u_{\mu}(\pbf, \sigma) = v_{\mu}^{\ast}(\pbf, \sigma) = (2p_0)^{-1/2} L(p)_{\mu}^{\nu} e_{\nu}(0, \sigma) \eqqcolon (2p_0)^{-1/2} e_{\mu}(\pbf, \sigma)
		\label{eq_vector_field_defn_e_vector_at_p}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{equation}
		e_{\mu}(0, 0) = \begin{bmatrix*}[r] 0 \\ 0 \\ 0 \\ 1 \end{bmatrix*}, \quad \
		e_{\mu}(0, 1) = -\frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ \ifrak \\ 0 \end{bmatrix*}, \quad \
		e_{\mu}(0, -1) = \frac{1}{\sqrt{2}} \begin{bmatrix*}[r] 0 \\ 1 \\ -\ifrak \\ 0 \end{bmatrix*}
		\label{eq_vector_field_defn_e_vector_at_rest}
	\end{equation}

Now we can rewrite the fields :math:`\eqref{eq_vector_field_psi_plus}` and :math:`\eqref{eq_vector_field_psi_minus}` as follows

.. math::
	:nowrap:

	\begin{align}
		\psi^+_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}}~\exp(\ifrak p \cdot x) e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) \nonumber \\
		\psi^-_{\mu}(x) &= \sum_{\sigma} (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}}~\exp(-\ifrak p \cdot x) e_{\mu}^{\ast}(\pbf, \sigma) a^{\dagger}(\pbf, \sigma) = \psi^{+ \dagger}_{\mu}(x)
		\label{eq_vector_field_psi_minus_adjoint_to_plus}
	\end{align}

Similar to the calculation :math:`\eqref{eq_scalar_field_commutator_as_Delta}` for scalar field, the (anti-)commutator can be calculated as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi^+_{\mu}(x), \psi^-_{\nu}(y) \right]_{\pm} = \int \frac{d^3 p}{(2\pi)^3 2p_0}~\exp(\ifrak p \cdot (x - y)) \Pi_{\mu \nu}(\pbf)
		\label{eq_vector_field_commutator_by_Pi}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{equation}
		\Pi_{\mu \nu}(\pbf) \coloneqq \sum_{\sigma} e_{\mu}(\pbf, \sigma) e^{\ast}_{\nu}(\pbf, \sigma)
		\label{eq_vector_field_Pi_matrix}
	\end{equation}

To better understand the quantity :math:`\Pi_{\mu \nu}(\pbf)`, let's first evaluate it at :math:`\pbf = 0` as follows

.. math::
	:nowrap:

	\begin{equation*}
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
	\end{equation*}

which is nothing but the projection to the spatial :math:`3`-space, or in other words, the orthogonal complement of the time direction :math:`\pbf = 0`. Considering the definition :math:`e_{\mu}(\pbf, \sigma) \coloneqq L(p)_{\mu}^{\nu} e_{\nu}(0, \sigma)` as in :math:`\eqref{eq_vector_field_defn_e_vector_at_p}`, we see that the general :math:`\Pi_{\mu \nu}(\pbf)` is really just a projection to the orthogonal complement of :math:`p`, and therefore can be written as

.. math::
	:nowrap:

	\begin{equation}
		\Pi_{\mu \nu}(\pbf) = \eta_{\mu \nu} + \frac{p_{\mu} p_{\nu}}{m^2}
		\label{eq_vector_field_defn_Pi}
	\end{equation}

because of the mass-shell condition :math:`p^2 + m^2 = 0`.

In light of :math:`\eqref{eq_scalar_field_commutator_as_Delta}`, we can rewrite :math:`\eqref{eq_vector_field_commutator_by_Pi}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi^+_{\mu}(x), \psi^-_{\nu}(y) \right]_{\pm} = \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x - y)
		\label{eq_vector_field_commutator_by_Delta}
	\end{equation}

where :math:`\Delta_+(x - y)` is defined by :math:`\eqref{eq_defn_Delta_plus}`. As in the case of scalar fields, this (anti-)commutator doesn't vanish even for space-like :math:`x - y`. Nonetheless, it's still an even function for space-like separations. The trick, as usual, is to consider a linear combination of :math:`\psi^+_{\mu}(x)` and :math:`\psi^-_{\mu}(x)` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\psi_{\mu}(x) \coloneqq \kappa \psi^+_{\mu}(x) + \lambda \psi^-_{\mu}(x)
	\end{equation*}

Now for space-separated :math:`x` and :math:`y`, we can calculate using :math:`\eqref{eq_vector_field_commutator_by_Delta}` and :math:`\eqref{eq_vector_field_psi_minus_adjoint_to_plus}` as follows

.. math::
	:nowrap:

	\begin{align*}
		\left[ \psi_{\mu}(x), \psi_{\nu}(y) \right]_{\pm} &= \kappa\lambda(1 \pm 1) \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x-y) \\
		\left[ \psi_{\mu}(x), \psi^{\dagger}_{\nu}(y) \right]_{\pm} &= (|\kappa|^2 \pm |\lambda|^2) \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta_+(x-y)
	\end{align*}

For them to vanishes, we see that first of all, we must adopt the top sign, i.e., take the commutator, or in other words, the vector field of spin :math:`1` must be bosonic. In addition, we must have :math:`|\kappa| = |\lambda|`. In fact, by adjusting the phase of the creation/annihilation operators, we can arrange so that :math:`\kappa = \lambda = 1`. To summarize, we can write a general vector field in the following form

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^-_{\mu}(x) = \psi^+_{\mu}(x) + \psi^{+ \dagger}_{\mu}(x)
		\label{eq_vector_field_psi_fixed_phase}
	\end{equation}

just like :math:`\eqref{eq_scalar_field_psi_fixed_phase}`. It's also obvious that :math:`\psi_{\mu}(x)` is Hermitian.

Now if the vector field carries a nonzero (conserved) quantum charge, then one must adjust :math:`\eqref{eq_vector_field_psi_fixed_phase}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\psi_{\mu}(x) \coloneqq \psi^+_{\mu}(x) + \psi^{+ c \dagger}_{\mu}(x)
	\end{equation*}

in analogy with :math:`\eqref{eq_scalar_field_psi_fixed_phase_with_antiparticle}` for scalar fields. Finally, we can express the vector field in terms of creation and annihilation operators as follows

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\mu}(x) = \sum_{\sigma} \int \frac{d^3 p}{(2\pi)^{3/2} (2p_0)^{1/2}}~\left[ e^{\ifrak p \cdot x} e_{\mu}(\pbf, \sigma) a(\pbf, \sigma) \
			+ e^{-\ifrak p \cdot x} e_{\mu}^{\ast}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right]
		\label{eq_vector_field_psi_by_creation_and_annihilation_operators}
	\end{equation}

in analogy with :math:`\eqref{eq_scalar_field_psi_by_creation_and_annihilation_operators}` for scalar fields. Finally, let's calculate the commutator (for general :math:`x` and :math:`y`) for later use as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi_{\mu}(x), \psi^{\dagger}_{\nu}(y) \right] = \left( \eta_{\mu \nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta(x-y)
		\label{eq_vector_field_commutator}
	\end{equation}

where :math:`\Delta(x-y)` as defined by :math:`\eqref{eq_defn_Delta}`.

So far, besides the introduction of the vectors :math:`e_{\mu}(\pbf, \sigma)` in :math:`\eqref{eq_vector_field_defn_e_vector_at_p}` and :math:`\eqref{eq_vector_field_defn_e_vector_at_rest}`, the discussion on vector fields looks very much like scalar fields. A key difference, however, stems from the following observation

.. math::
	:nowrap:

	\begin{equation}
		e^{\mu}(\pbf, \sigma) p_{\mu} = 0
		\label{eq_vector_field_spinor_orthogonal_to_momentum}
	\end{equation}

which, in turn, implies that

.. math::
	:nowrap:

	\begin{equation}
		\p_{\mu} \psi^{\mu}(x) = 0
		\label{eq_vector_field_gauge_fixing_condition}
	\end{equation}

This condition turns out to be coincide with a so-called "gauge fixing" condition for spin-:math:`1` photons in quantum electrodynamics. However, it's known that photons are massless particles. Therefore we may wonder if a vanishing mass limit :math:`m \to 0` may be applied. Now the simplest way to construct a (scalar) interaction density :math:`\Hscr(x)` using :math:`\psi_{\mu}(x)` is

.. math::
	:nowrap:

	\begin{equation}
		\Hscr(x) = J^{\mu}(x) \psi_{\mu}(x)
		\label{eq_vector_field_j_coupling}
	\end{equation}

where :math:`J^{\mu}(x)` is a :math:`4`-vector current. Suppose we fix the in- and out-states in the interaction. Then according to :math:`\eqref{eq_vector_field_psi_by_creation_and_annihilation_operators}`, the rate of (anti-)particle emission is proportional to

.. math::
	:nowrap:

	\begin{equation*}
		\sum_{\sigma} \left| \langle J^{\mu} \rangle e_{\mu}^{\ast}(\pbf, \sigma) \right|^2 \
			= \langle J^{\mu} \rangle \langle J^{\mu} \rangle^{\ast} \Pi_{\mu \nu}(\pbf) \
			\xlongequal{\eqref{eq_vector_field_defn_Pi}} \langle J^{\mu} \rangle \langle J^{\mu} \rangle^{\ast} \left( \eta_{\mu \nu} - p_{\mu} p_{\nu} / m^2 \right)
	\end{equation*}

where :math:`\langle J^{\mu} \rangle` denotes the matrix element of the current between the fixed in- and out-states. Now this rate blows up at :math:`m \to 0` limit unless :math:`p_{\mu} \langle J^{\mu} \rangle = 0`. This last condition can be translated to spacetime coordinates as follows

.. math::
	:nowrap:

	\begin{equation}
		\p_{\mu} J^{\mu}(x) = 0
		\label{eq_vector_field_j_coupling_condition}
	\end{equation}

or in other words :math:`J^{\mu}(x)` is a conserved current.

The CPT symmetries
++++++++++++++++++

Let's start with the spatial inversion. First recall from :ref:`sec_the_lorentz_and_cpt_transformation_laws`

.. math::
	:nowrap:

	\begin{align*}
		U(\Pcal) a(\pbf, \sigma) U^{-1}(\Pcal) &= \eta^{\ast} a(-\pbf, \sigma) \\
		U(\Pcal) a^{c \dagger}(\pbf, \sigma) U^{-1}(\Pcal) &= \eta^c a^{c \dagger}(-\pbf, \sigma)
	\end{align*}

It follows that we need to express :math:`e_{\mu}(-\pbf, \sigma)` in terms of :math:`e_{\mu}(\pbf, \sigma)`. To this end, let's calculate

.. math::
	:nowrap:

	\begin{equation}
		e_{\mu}(-\pbf, \sigma) = L(-\pbf)_{\mu}^{\nu} e_{\nu}(0, \sigma) \
			= \Pcal_{\mu}^{\rho} L(p)_{\rho}^{\tau} \Pcal_{\tau}^{\nu} e_{\nu}(0, \sigma) \
			= -\Pcal_{\mu}^{\rho} L(p)_{\rho}^{\tau} e_{\tau}(0, \sigma) \
			= -\Pcal_{\mu}^{\rho} e_{\rho}(p, \sigma)
		\label{eq_vector_field_revert_momentum_transformation}
	\end{equation}

It follows that the spatial inversion transformation law is given as follows

.. math::
	:nowrap:

	\begin{equation}
		U(\Pcal) \psi_{\mu}(x) U^{-1}(\Pcal) = -\eta^{\ast} \Pcal_{\mu}^{\nu} \psi_{\nu}(\Pcal x)
		\label{eq_vector_field_spatial_inversion_transformation_law}
	\end{equation}

under the following assumption :math:`\eta^c = \eta^{\ast}`. Omitting further details, the transformation laws for time inversion and charge conjugation is given as follows

.. math::
	:nowrap:

	\begin{align*}
		U(\Tcal) \psi_{\mu}(x) U^{-1}(\Tcal) &= \zeta^{\ast} \Pcal_{\mu}^{\nu} \psi_{\nu}(x) \\
		U(\Ccal) \psi_{\mu}(x) U^{-1}(\Ccal) &= \xi^{\ast} \psi_{\mu}^{\dagger}(x)
	\end{align*}

under the assumptions :math:`\zeta^c = \zeta^{\ast}` and :math:`\xi^c = \xi^{\ast}`.


Dirac fields
^^^^^^^^^^^^

Here we'll encounter the first nontrivial representation of the (homogeneous orthochronous) Lorentz group, first discovered by P. Dirac in a completely different (and more physical) context. Our treatment here will be purely mathematical, and will serve as a warm-up for the general representation theory.

.. _sec_dirac_representation_and_gamma_matrices:

Dirac representation and gamma matrices
+++++++++++++++++++++++++++++++++++++++

Let :math:`D` be a representation of the Lorentz group in the sense that :math:`D(\Lambda_1) D(\Lambda_2) = D(\Lambda_1 \Lambda_2)`. By the discussion in :ref:`sec_quantum_lorentz_symmetry` and ignoring the translation part, we can write :math:`D(\Lambda)` up to first order as follows

.. math::
	:nowrap:

	\begin{equation}
		D(\Lambda) = 1 + \frac{\ifrak}{2} \omega^{\mu \nu} \Jscr_{\mu \nu}
		\label{eq_dirac_field_linearize_representation}
	\end{equation}

where :math:`\Jscr_{\mu \nu} = -\Jscr_{\nu \mu}` are (Hermitian) matrices that, according to :math:`\eqref{eq_bracket_j4_j4}`, satisfy in addition the following Lie-algebraic condition

.. math::
	:nowrap:

	\begin{equation}
		\left[ \Jscr_{\mu \nu}, \Jscr_{\rho \kappa} \right] = \ifrak\left( \eta_{\mu \rho} \Jscr_{\nu \kappa} - \eta_{\nu \rho} \Jscr_{\mu \kappa} + \eta_{\kappa \mu} \Jscr_{\rho \nu} - \eta_{\kappa \nu} \Jscr_{\rho \mu} \right)
		\label{eq_bracket_repr_j}
	\end{equation}

The idea is to turn :math:`\eqref{eq_bracket_repr_j}` into a Jacobi identity. To this end, let's assume the existence of matrices :math:`\gamma_{\mu}` such that

.. math::
	:nowrap:

	\begin{equation}
		\left\{ \gamma_{\mu}, \gamma_{\nu} \right\} = 2\eta_{\mu \nu}
		\label{eq_dirac_field_clifford_algebra}
	\end{equation}

where the curly bracket denotes the anti-commutator and is equivalent to the notation :math:`[\cdots]_+` used in the previous chapters. Here the right-hand-side, being a bear number, should be interpreted as a multiple of the identity matrix of the same rank as :math:`\gamma_{\mu}`. Such matrices :math:`\gamma_{\mu}` form a so-called `Clifford algebra <https://en.wikipedia.org/wiki/Clifford_algebra>`_ of the symmetric bilinear form :math:`\eta_{\mu \nu}`. Then we can define

.. math::
	:nowrap:

	\begin{equation}
		\Jscr_{\mu \nu} \coloneqq -\frac{\ifrak}{4} \left[ \gamma_{\mu}, \gamma_{\nu} \right]
		\label{eq_dirac_field_defn_j}
	\end{equation}

Now we can calculate

.. math::
	:nowrap:

	\begin{align}
		\left[ \Jscr_{\mu \nu}, \gamma_{\rho} \right] &= -\frac{\ifrak}{4} \left[ \left[ \gamma_{\mu}, \gamma_{\nu} \right], \gamma_{\rho} \right] \label{eq_dirac_field_j_gamma_commutator} \\
			&= -\frac{\ifrak}{4} \left( \gamma_{\mu}\gamma_{\nu}\gamma_{\rho} - \gamma_{\nu}\gamma_{\mu}\gamma_{\rho} - \gamma_{\rho}\gamma_{\mu}\gamma_{\nu} + \gamma_{\rho}\gamma_{\nu}\gamma_{\mu} \right) \nonumber \\
			&= -\frac{\ifrak}{4} \left( (\gamma_{\mu}\gamma_{\nu}\gamma_{\rho} + \gamma_{\mu}\gamma_{\rho}\gamma_{\nu}) - (\gamma_{\mu}\gamma_{\rho}\gamma_{\nu} + \gamma_{\rho}\gamma_{\mu}\gamma_{\nu}) - (\gamma_{\nu}\gamma_{\mu}\gamma_{\rho} + \gamma_{\nu}\gamma_{\rho}\gamma_{\mu}) + (\gamma_{\nu}\gamma_{\rho}\gamma_{\mu} + \gamma_{\rho}\gamma_{\nu}\gamma_{\mu}) \right) \nonumber \\
			&= -\frac{\ifrak}{4} \left( 2\gamma_{\mu}\eta_{\nu \rho} - 2\eta_{\mu \rho}\gamma_{\nu} - 2\gamma_{\nu}\eta_{\mu \rho} + 2\eta_{\nu \rho}\gamma_{\mu} \right) \nonumber \\
			&= -\ifrak \eta_{\nu \rho}\gamma_{\mu} + \ifrak \eta_{\mu \rho}\gamma_{\nu} \nonumber
	\end{align}

Then we can verify :math:`\eqref{eq_bracket_repr_j}`, starting from the left-hand-side, as follows

.. math::
	:nowrap:

	\begin{align*}
		\left[ \Jscr_{\mu \nu}, \Jscr_{\rho \kappa} \right] &= -\frac{\ifrak}{4} \left[ \Jscr_{\mu \nu}, \left[ \gamma_{\rho}, \gamma_{\kappa} \right] \right] \\
			&= \frac{\ifrak}{4} \left[ \gamma_{\rho}, \left[ \gamma_{\kappa}, \Jscr_{\mu \nu} \right] \right] + \frac{\ifrak}{4} \left[ \gamma_{\kappa}, \left[ \Jscr_{\mu \nu}, \gamma_{\rho} \right] \right] \\
			&= \frac{\ifrak}{4} \left[ \gamma_{\rho}, \left( \ifrak \eta_{\nu \kappa} \gamma_{\mu} - \ifrak \eta_{\mu \kappa} \gamma_{\nu} \right) \right] + \frac{\ifrak}{4} \left[ \gamma_{\kappa}, \left( -\ifrak \eta_{\nu \rho} \gamma_{\mu} + \ifrak \eta_{\mu \rho} \gamma_{\nu} \right) \right] \\
			&= -\ifrak \eta_{\nu \kappa} \Jscr_{\rho \mu} + \ifrak \eta_{\mu \kappa} \Jscr_{\rho \nu} + \ifrak \eta_{\nu \rho} \Jscr_{\kappa \mu} - \ifrak \eta_{\mu \rho} \Jscr_{\kappa \nu}
	\end{align*}

The last expression is easily seen to be equal to the right-hand-side of :math:`\eqref{eq_bracket_repr_j}` using the anti-symmetry of :math:`\Jscr_{\mu \nu}`. Hence :math:`\eqref{eq_bracket_repr_j}` is indeed satisfied for :math:`\Jscr_{\mu \nu}` defined by :math:`\eqref{eq_dirac_field_defn_j}`.

In fact, the calculation :math:`\eqref{eq_dirac_field_j_gamma_commutator}` may be rephrased more conveniently as follows

.. math::
	:nowrap:

	\begin{equation}
		D(\Lambda) \gamma_{\mu} D^{-1}(\Lambda) = \left( \Lambda^{-1} \right)_{\mu}^{\nu} \gamma_{\nu}
		\label{eq_dirac_field_gamma_is_vector}
	\end{equation}

or in plain words, :math:`\gamma_{\mu}` is a vector. Indeed, one can calculate the left-hand-side of :math:`\eqref{eq_dirac_field_gamma_is_vector}` up to first order as follows [#dirac_gamma_is_vector]_

.. math::
	:nowrap:

	\begin{align*}
		D(\Lambda) \gamma_{\mu} D^{-1}(\Lambda) &= \left( 1 + \frac{\ifrak}{2} \omega^{\nu \rho} \Jscr_{\nu \rho} \right) \gamma_{\mu} \left( 1 - \frac{\ifrak}{2} \omega^{\nu \rho} \Jscr_{\nu \rho} \right) \\
			&= \gamma_{\mu} + \frac{\ifrak}{2} \omega^{\nu \rho} \left[ \Jscr_{\nu \rho}, \gamma_{\mu} \right] \\
			&= \gamma_{\mu} + \frac{\ifrak}{2} \omega^{\nu \rho} \left( -\ifrak \eta_{\rho \mu} \gamma_{\nu} + \ifrak \eta_{\nu \mu} \gamma_{\rho} \right) \\
			&= \left( \delta_{\mu}^{\nu} - \omega_{\mu}^{\nu} \right) \gamma_{\nu} = \left( \Lambda^{-1} \right)_{\mu}^{\nu} \gamma_{\nu}
	\end{align*}

Using the very definition :math:`\eqref{eq_dirac_field_defn_j}`, one then sees that :math:`\Jscr_{\mu \nu}` is an anti-symmetric tensor in the sense that

.. math::
	:nowrap:

	\begin{equation*}
		D(\Lambda) \Jscr_{\mu \nu} D^{-1}(\Lambda) = \left( \Lambda^{-1} \right)_{\mu}^{\rho} \left( \Lambda^{-1} \right)_{\nu}^{\kappa} \Jscr_{\rho \kappa}
	\end{equation*}

Indeed, one can construct even more anti-symmetric tensors as follows

.. math::
	:nowrap:

	\begin{align*}
		\Ascr_{\mu \nu \rho} &\coloneqq \gamma_{\mu}\gamma_{\nu}\gamma_{\rho} \pm \text{ signed permutations} \\
		\Pscr_{\mu \nu \rho \kappa} &\coloneqq \gamma_{\mu}\gamma_{\nu}\gamma_{\rho}\gamma_{\kappa} \pm \text{ signed permutations}
	\end{align*}

but no more because we're constrained by the :math:`4`-dimensional spacetime. Here we remark also that the consideration of only anti-symmetric tensors loses no generality because of :math:`\eqref{eq_dirac_field_clifford_algebra}`. Now we claim that the matrices :math:`1, \gamma_{\mu}, \Jscr_{\mu \nu}, \Ascr_{\mu \nu \rho}` and :math:`\Pscr_{\mu \nu \rho \kappa}` are all linearly independent. This can be seen either by observing that they transform differently under conjugation by :math:`D(\Lambda)`, or, more directly, by observing that they are orthogonal to each other under the inner product defined by the trace of the product (of two matrices). For example, one easily sees that :math:`\Jscr_{\mu \nu}` is traceless because the trace of a commutator vanishes, and :math:`\op{tr}(\gamma_{\mu} \gamma_{\nu}) = 0` for :math:`\mu \neq \nu` using :math:`\eqref{eq_dirac_field_clifford_algebra}`. It's slightly tricker to see that :math:`\gamma_{\mu}` itself is also traceless, but this is again a consequence of the Clifford algebra relations :math:`\eqref{eq_dirac_field_clifford_algebra}`, which we demonstrate as follows

.. math::
	:nowrap:

	\begin{equation*}
		\op{tr}(\gamma_{\mu}) = \op{tr}(\gamma_{\nu} \gamma_{\mu} \gamma^{-1}_{\nu}) = -\op{tr}(\gamma_{\mu} \gamma_{\nu} \gamma^{-1}_{\nu}) = -\op{tr}(\gamma_{\mu}) \
			\implies \op{tr}(\gamma_{\mu}) = 0
	\end{equation*}

Counting these linearly independent matrices, we see that there are :math:`1 + \binom{4}{1} + \binom{4}{2} + \binom{4}{3} + \binom{4}{4} = 16` of them. It means that the size of the :math:`\gamma_{\mu}` matrices is at least :math:`4 \times 4`.

It turns out that there exists indeed a solution of :math:`\eqref{eq_dirac_field_clifford_algebra}` in terms of :math:`4 \times 4` matrices, conveniently known as the `gamma matrices <https://en.wikipedia.org/wiki/Gamma_matrices>`_, which we define as follows

.. math::
	:nowrap:

	\begin{equation}
		\gamma_0 \coloneqq -\ifrak \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad \
			\bm{\gamma} \coloneqq -\ifrak \begin{bmatrix*}[r] 0 & \bm{\sigma} \\ -\bm{\sigma} & 0 \end{bmatrix*}
		\label{eq_dirac_field_defn_gamma_matrices}
	\end{equation}

where :math:`\bm{\sigma} = (\sigma_1, \sigma_2, \sigma_3)` is made up of the so-called `Pauli matrices <https://en.wikipedia.org/wiki/Pauli_matrices>`_ defined as follows

.. math::
	:nowrap:

	\begin{equation}
		\sigma_1 \coloneqq \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad \
			\sigma_2 \coloneqq \begin{bmatrix*}[r] 0 & -\ifrak \\ \ifrak & 0 \end{bmatrix*}, \quad \
			\sigma_3 \coloneqq \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}
		\label{eq_pauli_matrices}
	\end{equation}

Indeed the Pauli matrices make up a solution to not only a :math:`3`-dimensional Clifford algebra with respect to the Euclidean inner product, but also an angular momentum representation if multiplied by :math:`1/2`, or more precisely, a spin-:math:`1/2` representation. Note that this representation, in terms of Hermitian matrices, is apparently different from the one given by :math:`\eqref{eq_representation_rotation_first_order} - \eqref{eq_j3_matrix}` with :math:`\jfrak = 1/2`, in terms of real matrices. One can verify that they differ by just a change of basis. As long as the terminology is concerned, one often uses Hermitian and real interchangeably.

Now using the definition :math:`\eqref{eq_dirac_field_defn_j}`, we can calculate, using the Clifford relations, as follows

.. math::
	:nowrap:

	\begin{align}
		\Jscr_{ij} &= -\frac{\ifrak}{4} \left[ \gamma_i, \gamma_j \right] \
			= -\frac{\ifrak}{2} \epsilon_{ij} \gamma_i \gamma_j \
			= -\frac{\ifrak}{2} \epsilon_{ij} \begin{bmatrix} \sigma_i \sigma_j & 0 \\ 0 & \sigma_i \sigma_j \end{bmatrix} \
			\xlongequal{\eqref{eq_jjj_commutation}} \frac{1}{2} \epsilon_{ijk} \begin{bmatrix} \sigma_k & 0 \\ 0 & \sigma_k \end{bmatrix}  \label{eq_dirac_field_jscr_matrix} \\
		\Jscr_{i0} &= -\frac{\ifrak}{4} \left[ \gamma_i, \gamma_0 \right] = -\frac{\ifrak}{2} \gamma_i \gamma_0 = \frac{\ifrak}{2} \begin{bmatrix} \sigma_i & 0 \\ 0 & \sigma_i \end{bmatrix}  \nonumber
	\end{align}

.. _paragraph_dirac_field_representation_not_unitary:

where :math:`i, j \in \{1,2,3\}` and :math:`\epsilon` is the totally anti-symmetric sign. We see that the representation :math:`\Jscr_{\mu \nu}` is in fact reducible. Moreover, we see that that the corresponding representation :math:`D` of the Lorentz group given by :math:`\eqref{eq_dirac_field_linearize_representation}` is not unitary, since while :math:`\Jscr_{ij}` are Hermitian, :math:`\Jscr_{i0}` are anti-Hermitian. The fact that :math:`D` is not unitary will have consequences when we try to construct the interaction density as in :math:`\eqref{eq_construct_interaction_density_by_fields}`, because products like :math:`\psi^{\dagger} \psi` will not be a scalar in light of :math:`\eqref{eq_conjugate_annihilation_field}` or :math:`\eqref{eq_conjugate_creation_field}`.

Next let's consider the parity transformation in the formalism of gamma matrices. In light of the transformation laws :math:`\eqref{eq_j3_conjugated_by_p}` and :math:`\eqref{eq_k3_conjugated_by_p}`, we can define

.. math::
	:nowrap:

	\begin{equation}
		\beta \coloneqq \ifrak \gamma_0 = \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}
		\label{eq_dirac_field_beta_matrix}
	\end{equation}

as the parity transformation. It's clear that :math:`\beta^2 = 1`. Moreover, it follows from the Clifford relations :math:`\eqref{eq_dirac_field_clifford_algebra}` that

.. math::
	:nowrap:

	\begin{equation*}
		\beta \gamma_i \beta^{-1} = -\gamma_i, \quad \beta \gamma_0 \beta^{-1} = \gamma_0
	\end{equation*}

It follows then from the definition :math:`\eqref{eq_dirac_field_defn_j}` that

.. math::
	:nowrap:

	\begin{equation*}
		\beta \Jscr_{ij} \beta^{-1} = \Jscr_{ij}, \quad \beta \Jscr_{i0} \beta^{-1} = - \Jscr_{i0}
	\end{equation*}

which are in complete analogy with :math:`\eqref{eq_j3_conjugated_by_p}` and :math:`\eqref{eq_k3_conjugated_by_p}` if we think of :math:`\Jscr_{ij}` as the "angular momenta" and :math:`\Jscr_{i0}` as the "boosts".

In connection to the non-unitarity of :math:`D(\Lambda)`, let's note that since :math:`\beta \gamma_{\mu}^{\dagger} \beta^{-1} = -\gamma_{\mu}`, we have :math:`\beta \Jscr_{\mu \nu}^{\dagger} \beta^{-1} = \Jscr_{\mu \nu}`, and therefore

.. math::
	:nowrap:

	\begin{equation}
		\beta D^{\dagger}(\Lambda) \beta^{-1} = D^{-1}(\Lambda)
		\label{eq_dirac_field_pseudo_unitarity_of_d_matrix}
	\end{equation}

in light of :math:`\eqref{eq_dirac_field_linearize_representation}`. This identity will be useful when we later construct the interaction density.

At last we'll introduce yet another special element to the family of gamma matrices, namely :math:`\gamma_5`, defined as follows

.. math::
	:nowrap:

	\begin{equation}
		\gamma_5 \coloneqq -\ifrak \gamma_0 \gamma_1 \gamma_2 \gamma_3 = \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}
		\label{eq_dirac_field_defn_gamma_5}
	\end{equation}

One nice thing about :math:`\gamma_5` is that it anti-commutes with all other :math:`\gamma` matrices, and in particular

.. math::
	:nowrap:

	\begin{equation}
		\beta \gamma_5 = -\gamma_5 \beta
		\label{eq_dirac_field_gamma_5_anti_commutes_beta}
	\end{equation}

In fact, the collection :math:`\gamma_0, \gamma_1, \gamma_2, \gamma_3, \gamma_5` makes exactly a :math:`5`-dimensional spacetime Clifford algebra.


Construction of Dirac fields
++++++++++++++++++++++++++++

As in the case of scalar and vector fields, let's write the Dirac fields as follows

.. math::
	:nowrap:

	\begin{align}
		\psi^+_{\ell}(x) &= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma) a(\pbf, \sigma)
		\label{eq_dirac_field_psi_plus} \\
		\psi^{-c}_{\ell}(x) &= (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p~e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma)
		\label{eq_dirac_field_psi_minus}
	\end{align}

Moreover, using the particularly convenient form of :math:`\Jscr_{ij}` in :math:`\eqref{eq_dirac_field_jscr_matrix}`, we can write the zero-momentum conditions :math:`\eqref{eq_j_intertwines_u}` and :math:`\eqref{eq_j_intertwines_v}` as follows

.. math::
	:nowrap:

	\begin{align}
		\frac{1}{2} \sum_m \bm{\sigma}_{m' m} u_{m \pm}(0, \sigma) &= \sum_{\sigma'} u_{m' \pm}(0, \sigma') \Jbf^{(\jfrak)}_{\sigma' \sigma}
		\label{eq_dirac_field_sigma_intertwines_j_by_u} \\
		-\frac{1}{2} \sum_m \bm{\sigma}_{m' m} v_{m \pm}(0, \sigma) &= \sum_{\sigma'} v_{m' \pm}(0, \sigma') \Jbf^{(\jfrak) \ast}_{\sigma' \sigma}
		\label{eq_dirac_field_sigma_intertwines_j_by_v}
	\end{align}

Here the signs :math:`\pm` correspond to the two (identical) irreducible representations of :math:`\Jscr_{ij}` in the form of :math:`\eqref{eq_dirac_field_jscr_matrix}`, and :math:`m` and :math:`m'` index the Pauli matrices.

Now if we think of :math:`u_{m \pm}(0, \sigma)` as matrix elements of a matrix :math:`U_{\pm}` and similarly for :math:`v`, then :math:`\eqref{eq_dirac_field_sigma_intertwines_j_by_u}` and :math:`\eqref{eq_dirac_field_sigma_intertwines_j_by_v}` can be rewritten as

.. math::
	:nowrap:

	\begin{align}
		\frac{1}{2} \bm{\sigma} U_{\pm} &= U_{\pm} \Jbf^{(\jfrak)}
		\label{eq_dirac_field_sigma_intertwines_j_by_u_matrix_form} \\
		-\frac{1}{2} \bm{\sigma} V_{\pm} &= V_{\pm} \Jbf^{(\jfrak) \ast}
		\label{eq_dirac_field_sigma_intertwines_j_by_v_matrix_form}
	\end{align}

We recall that both :math:`\tfrac{1}{2} \bm{\sigma}` and :math:`\Jbf^{(\jfrak)}` (as well as :math:`-\Jbf^{(\jfrak) \ast}`) are irreducible representations of the rotation group (or rather, its Lie algebra). We first claim that :math:`U_{\pm}` must be isomorphism. Indeed, the kernel of :math:`U_{\pm}` is easily seen to be an invariant subspace under the action of :math:`\Jbf^{(\jfrak)}`, and hence must be null if :math:`U_{\pm} \neq 0`. On the other hand, the image of :math:`U_{\pm}` is an invariant subspace under the action of :math:`\bm{\sigma}`, and hence must be the whole space if :math:`U_{\pm} \neq 0`. It follows then that the rank of :math:`\Jbf^{(\jfrak)}` and :math:`U_{\pm}` must be the same as :math:`\bm{\sigma}`, which is :math:`2`. The same argument applies also to :math:`V_{\pm}`. In particular we must have :math:`\jfrak = 1/2`, or in other words, the Dirac fields describe :math:`1/2`-spin particles.

.. note::
	The mathematical argument above is commonly known as `Schur's lemma <https://en.wikipedia.org/wiki/Schur%27s_lemma>`_.

The matrix form of :math:`\Jbf^{(1/2)}` according to :math:`\eqref{eq_j1_j2_matrix}` and :math:`\eqref{eq_j3_matrix}` is given By

.. math::
	:nowrap:

	\begin{equation*}
		J_1^{(1/2)} = \frac{1}{2} \begin{bmatrix*}[r] 0 & 1 \\ 1 & 0 \end{bmatrix*}, \quad \
			J_2^{(1/2)} = \frac{\ifrak}{2} \begin{bmatrix*}[r] 0 & -1 \\ 1 & 0 \end{bmatrix*}, \quad \
			J_3^{(1/2)} = \frac{1}{2} \begin{bmatrix*}[r] 1 & 0 \\ 0 & -1 \end{bmatrix*}
	\end{equation*}

Comparing with the Pauli matrices :math:`\eqref{eq_pauli_matrices}`, we see that

.. math::
	:nowrap:

	\begin{equation*}
		\Jbf^{(1/2)} = \frac{1}{2} \bm{\sigma}, \quad -\Jbf^{(1/2) \ast} = \frac{1}{2} \sigma_2 \bm{\sigma} \sigma_2
	\end{equation*}

Hence :math:`\eqref{eq_dirac_field_sigma_intertwines_j_by_u_matrix_form}` may be rewritten as :math:`\bm{\sigma} U_{\pm} = U_{\pm} \bm{\sigma}`. One can apply Schur's lemma once again to conclude that :math:`U_{\pm}` must be a scalar. To spell out more details, note that since :math:`\bm{\sigma}` commutes with any scalar it commutes with :math:`U_{\pm} - \lambda` for any :math:`\lambda \in \Cbb`. It follows that :math:`U_{\pm} - \lambda` is either an isomorphism or zero. The later must be the case if :math:`\lambda` is an eigenvalue of :math:`U_{\pm}`. Hence :math:`U_{\pm}` must be a scalar as claimed. A similar argument can be applied to :math:`V_{\pm}` by rewriting :math:`\eqref{eq_dirac_field_sigma_intertwines_j_by_v_matrix_form}` as :math:`\bm{\sigma} (V_{\pm} \sigma_2) = (V_{\pm} \sigma_2) \bm{\sigma}`. Hence :math:`V_{\pm}` must be proportional to :math:`\sigma_2`.

Going back to :math:`u_{m \pm}(0, \sigma)` (resp. :math:`v_{m \pm}(0, \sigma)`) from :math:`U_{\pm}` (resp. :math:`V_{\pm}`), we have concluded that

.. math::
	:nowrap:

	\begin{align*}
		u_{m \pm}(0, \sigma) &= c_{\pm} \delta_{m \sigma} \\
		v_{m \pm}(0, \sigma) &= -\ifrak d_{\pm} (\sigma_2)_{m \sigma}
	\end{align*}

where we've inserted the extra factor :math:`-\ifrak` in front of :math:`d_{\pm}` to make the final results look more uniform. We can further unwrap this result in matrix notations as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		u(0, 1/2) &= \begin{bmatrix*}[l] c_+ \\ 0 \\ c_- \\ 0 \end{bmatrix*}, \quad u(0, -1/2) &&= \begin{bmatrix*}[l] 0 \\ c_+ \\ 0 \\ c_- \end{bmatrix*}
		\label{eq_dirac_field_u_matrix_by_c} \\
		v(0, 1/2) &= \begin{bmatrix*}[l] 0 \\ d_+ \\ 0 \\ d_- \end{bmatrix*}, \quad v(0, -1/2) &&= -\begin{bmatrix*}[l] d_+ \\ 0 \\ d_- \\ 0 \end{bmatrix*}
		\label{eq_dirac_field_v_matrix_by_d}
	\end{alignat}

These are known as the (Dirac) spinors. Spinors at finite momentum can be determined by :math:`\eqref{eq_annihilation_u_transformation_simplified_by_boost}` and :math:`\eqref{eq_creation_v_transformation_simplified_by_boost}` as follows

.. math::
	:nowrap:

	\begin{align}
		u(\pbf, \sigma) &= \sqrt{m/p_0} D(L(p)) u(0, \sigma)
		\label{eq_dirac_field_spinor_u_at_finite_momentum} \\
		v(\pbf, \sigma) &= \sqrt{m/p_0} D(L(p)) v(0, \sigma)
		\label{eq_dirac_field_spinor_v_at_finite_momentum}
	\end{align}

In general the constants :math:`c_{\pm}` and :math:`d_{\pm}` may be arbitrary. However, if we assume the conservation of parity, or in other words, that the fields transform covariantly under the spatial inversion, then we can determine the constants, and henceforth the spinors uniquely as follows. Applying the spatial inversion transformation laws of creation and annihilation operators from :ref:`sec_the_lorentz_and_cpt_transformation_laws` to :math:`\eqref{eq_dirac_field_psi_plus}` and :math:`\eqref{eq_dirac_field_psi_minus}`, we can calculate as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Pcal) \psi^+(x) U^{-1}(\Pcal) &= (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot x} u(\pbf, \sigma) a(-\pbf, \sigma)
			\label{eq_dirac_field_spatial_inversion_acts_on_psi_plus} \\
			&= (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot \Pcal x} u(-\pbf, \sigma) a(\pbf, \sigma) \nonumber \\
			&= (2\pi)^{-3/2} \eta^{\ast} \sum_{\sigma} \int d^3 p~e^{\ifrak p \cdot \Pcal x} \sqrt{m/p_0}~\beta D(L(p)) \beta u(0, \sigma) a(\pbf, \sigma) \nonumber \\
		U(\Pcal) \psi^{- c}(x) U^{-1}(\Pcal) &= (2\pi)^{-3/2} \eta^c \sum_{\sigma} \int d^3 p~e^{-\ifrak p \cdot \Pcal x} \sqrt{m/p_0}~\beta D(L(p))\beta v(0, \sigma) a^{c \dagger}(\pbf, \sigma)
		\label{eq_dirac_field_spatial_inversion_acts_on_psi_minus}
	\end{align}

Here we've evaluated :math:`u(-\pbf, \sigma)` (and :math:`v(-\pbf, \sigma)`) in a way similar to :math:`\eqref{eq_vector_field_revert_momentum_transformation}`, where, instead of applying :math:`\Pcal` directly on vectors, we recall for Dirac fields that the spatial inversion acts as :math:`\beta`, defined by :math:`\eqref{eq_dirac_field_beta_matrix},` on spinors.

In order for :math:`U(\Pcal) \psi^+(x) U^{-1}(\Pcal)` to be proportional to :math:`\psi^+(x)`, we observe that it would be sufficient (and most likely also necessary) if :math:`u(0, \sigma)` is an eigenvector of :math:`\beta`. Similar argument applies to :math:`\psi^-(x)` as well, so we can write

.. math::
	:nowrap:

	\begin{align}
		\beta u(0, \sigma) &= b_+ u(0, \sigma)
		\label{eq_dirac_field_beta_eigenvalue_of_u} \\
		\beta v(0, \sigma) &= b_- v(0, \sigma)
		\label{eq_dirac_field_beta_eigenvalue_of_v}
	\end{align}

where :math:`b^2_{\pm} = 1` since :math:`\beta^2 = 1`. Assuming this, we can rewrite :math:`\eqref{eq_dirac_field_spatial_inversion_acts_on_psi_plus}` and :math:`\eqref{eq_dirac_field_spatial_inversion_acts_on_psi_minus}` as follows

.. math::
	:nowrap:

	\begin{align}
		U(\Pcal) \psi^+(x) U^{-1}(\Pcal) &= \eta^{\ast} b_+ \beta \psi^+(\Pcal x)
		\label{eq_dirac_field_psi_plus_conjugated_by_u} \\
		U(\Pcal) \psi^{- c}(x) U^{-1}(\Pcal) &= \eta^c b_- \beta \psi^{- c}(\Pcal x)
		\label{eq_dirac_field_psi_minus_conjugated_by_u}
	\end{align}

so that the parity is indeed conserved.

Now using :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_u}` and :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_v}` (and an appropriate rescaling), we can rewrite :math:`\eqref{eq_dirac_field_u_matrix_by_c}` and :math:`\eqref{eq_dirac_field_v_matrix_by_d}` as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		u(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 1 \\ 0 \\ b_+ \\ 0 \end{bmatrix*}, \quad u(0, -1/2) &&= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 0 \\ 1 \\ 0 \\ b_+ \end{bmatrix*} \\
		v(0, 1/2) &= \frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 0 \\ 1 \\ 0 \\ b_- \end{bmatrix*}, \quad v(0, -1/2) &&= -\frac{1}{\sqrt{2}} \begin{bmatrix*}[l] 1 \\ 0 \\ b_- \\ 0 \end{bmatrix*}
	\end{alignat*}

So far we haven't really achieved much by assuming the parity conservation. What eventually pins down the spinors is, once again, the causality condition. To see this, let's try to construct the Dirac field following :math:`\eqref{eq_defn_psi_field}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\psi(x) \coloneqq \kappa \psi^+(x) + \lambda \psi^{- c}(x)
		\label{eq_dirac_field_psi_field_raw}
	\end{equation}

The (anti-)commutator can be calculated using :math:`\eqref{eq_dirac_field_psi_plus}` and :math:`\eqref{eq_dirac_field_psi_minus}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_{\pm} = (2\pi)^{-3} \int d^3 p~\left( |\kappa|^2 N_{\ell \ell'}(\pbf) \exp(\ifrak p \cdot (x-y)) \pm |\lambda|^2 M_{\ell \ell'}(\pbf) \exp(-\ifrak p \cdot (x-y)) \right)
		\label{eq_dirac_field_commutator_raw}
	\end{equation}

where :math:`N_{\ell \ell'}` and :math:`M_{\ell \ell'}` are the spin sums defined as follows

.. math::
	:nowrap:

	\begin{align}
		N_{\ell \ell'}(\pbf) &= \sum_{\sigma} u_{\ell}(\pbf, \sigma) u^{\ast}_{\ell'}(\pbf, \sigma)
		\label{eq_dirac_field_n_matrix_as_spinor_sum} \\
		M_{\ell \ell'}(\pbf) &= \sum_{\sigma} v_{\ell}(\pbf, \sigma) v^{\ast}_{\ell'}(\pbf, \sigma)
		\label{eq_dirac_field_m_matrix_as_spinor_sum}
	\end{align}

To evaluate the spin sums, we first turn back to the zero-momentum case and use :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_u}` and :math:`\eqref{eq_dirac_field_beta_eigenvalue_of_v}` to express the values in terms of :math:`\beta` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		N(0) &= \sum_{\sigma} u(0, \sigma) u^{\dagger}(0, \sigma) &&= \frac{1 + b_+ \beta}{2}
		\label{eq_dirac_field_n_matrix_at_zero} \\
		M(0) &= \sum_{\sigma} v(0, \sigma) v^{\dagger}(0, \sigma) &&= \frac{1 + b_- \beta}{2}
		\label{eq_dirac_field_m_matrix_at_zero}
	\end{alignat}

where :math:`\dagger` here means transpose conjugation. The easiest way to see this is probably to momentarily forget about the spin :math:`z`-component :math:`\sigma` so that :math:`\beta` as in :math:`\eqref{eq_dirac_field_beta_matrix}` behaves like a :math:`2 \times 2` matrix with obvious eigenvectors :math:`[1, 1]^T` for :math:`b_{\pm} = 1` and :math:`[1, -1]^T` for :math:`b_{\pm} = -1`. Then :math:`\eqref{eq_dirac_field_n_matrix_at_zero}` and :math:`\eqref{eq_dirac_field_m_matrix_at_zero}` can be verified by a direct calculation. Here the superscript :math:`T` means taking transpose so we have column vectors.

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

To go further, we need to invoke the the transformation laws of the gamma matrices and their relatives, e.g., :math:`\beta`, under :math:`D(\Lambda)` from :ref:`sec_dirac_representation_and_gamma_matrices`. More precisely, using the pseudo-unitarity :math:`\eqref{eq_dirac_field_pseudo_unitarity_of_d_matrix}` we have the following

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

Now that we've finished evaluating the spin sums, we can plug them into :math:`\eqref{eq_dirac_field_commutator_raw}` to get the following evaluation of the (anti-)commutator

.. math::
	:nowrap:

	\begin{align}
		\left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_{\pm} &= (2\pi)^{-3} \int \frac{d^3 p}{2p_0}~\big( |\kappa|^2 (-\ifrak p_{\mu} \gamma^{\mu} + b_+ m) e^{\ifrak p \cdot (x-y)} \beta \phantom{\big)}
		\label{eq_dirac_field_commutator_first_evaluation} \\
			&\phantom{= \phantom{\big(}} \pm |\lambda|^2 (-\ifrak p_{\mu} \gamma^{\mu} + b_- m) e^{-\ifrak p \cdot (x-y)} \beta \big)_{\ell \ell'} \nonumber \\
			&= \left( |\kappa|^2 \left( -\ifrak \gamma^{\mu} \p_{x_{\mu}} + b_+ m \right) \Delta_+(x-y) \beta \pm |\lambda|^2 \left( -\ifrak \gamma^{\mu} \p_{y_{\mu}} + b_- m \right) \Delta_+(y-x) \beta \right)_{\ell \ell'} \nonumber
	\end{align}

where :math:`\Delta_+` is defined by :math:`\eqref{eq_defn_Delta_plus}`. Recall that :math:`\Delta_+(x) = \Delta_+(-x)` for space-like :math:`x`. Hence for space-like :math:`x-y`, the following holds

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

By the usual phase adjustments on the creation and annihilation operators and rescaling, we can arrange so that :math:`\kappa = \lambda = 1`. Recalling :math:`\eqref{eq_dirac_field_gamma_5_anti_commutes_beta}` and replacing :math:`\psi` with :math:`\gamma_5 \psi` if necessary, we can arrange so that :math:`b_{\pm} = \pm 1`. Putting these all together, we have evaluated the Dirac field :math:`\eqref{eq_dirac_field_psi_field_raw}` as follows

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

and the anti-commutator, calculated by plugging the spin sums into :math:`\eqref{eq_dirac_field_commutator_raw}`, is

.. math::
	:nowrap:

	\begin{equation}
		\left[ \psi_{\ell}(x), \psi_{\ell'}^{\dagger}(y) \right]_+ = \left( \left( -\gamma^{\mu}\p_{\mu} + m \right) \beta \right)_{\ell \ell'} \Delta(x-y)
		\label{eq_dirac_field_commutator}
	\end{equation}

where :math:`\p_{\mu} = \p_{x_{\mu}}` and :math:`\Delta(x-y)` is defined by :math:`\eqref{eq_defn_Delta}`.

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

	However, unlike the original derivation of Dirac, we've derived it here as a consequence of parity conservation. In fact, we note that the Dirac equation is not something completely different from the Klein-Gordon equation, because using the Clifford algebra relations :math:`\eqref{eq_dirac_field_clifford_algebra}`, we see that :math:`\gamma^{\mu} \p_{\mu}` is actually a square root of the D'Alembert operator :math:`\square`. This is allegedly one of the motivations of Dirac to find these gamma matrices in the first place.

The CPT symmetries
++++++++++++++++++

The transformation law of Dirac fields under spatial inversion has already been worked out in :math:`\eqref{eq_dirac_field_spatial_inversion_transformation_law}`, so we're left to work out the transformation laws under time inversion and charge conjugation.

Recall from :ref:`sec_space_and_time_inversions` that the time-inversion operator :math:`U(\Tcal)` is complex anti-linear. Hence to work out the transformation law under time inversion, we'll need to work out the complex-conjugated spinors :math:`u^{\ast}(\pbf, \sigma)` and :math:`v^{\ast}(\pbf, \sigma)`. Now in light of :math:`\eqref{eq_dirac_field_spinor_u_at_finite_momentum}` and :math:`\eqref{eq_dirac_field_spinor_v_at_finite_momentum}`, and the fact that the spinors are real at zero-momentum, we just need to work out the complex conjugate :math:`D^{\ast}(L(p))` in terms of :math:`D(L(p))` and the gamma matrices. Now according to :math:`\eqref{eq_dirac_field_linearize_representation}`, it suffices to work out :math:`\Jscr^{\ast}_{\mu \nu}`. Finally according to :math:`\eqref{eq_dirac_field_defn_j}`, it suffices to work out :math:`\gamma^{\ast}_{\mu}`.

Inspecting the explicit forms of the gamma matrices given by :math:`\eqref{eq_dirac_field_defn_gamma_matrices}` and :math:`\eqref{eq_pauli_matrices}` we see that :math:`\gamma_0, \gamma_1, \gamma_3` are anti-Hermitian while :math:`\gamma_2` is Hermitian, or more explicitly

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

where we also note that :math:`(\Cscr^{-1} \beta)^{-1} = \beta \Cscr`. It follows from :math:`\eqref{eq_dirac_field_defn_j}` that

.. math::
	:nowrap:

	\begin{equation*}
		\Jscr^{\ast}_{\mu \nu} = -\beta\Cscr \Jscr_{\mu \nu} \Cscr^{-1}\beta
	\end{equation*}

and hence from :math:`\eqref{eq_dirac_field_linearize_representation}` that

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

These relations turns out to be useful for the charge conjugation transformation, but not for the time inversion because the spinors :math:`u` and :math:`v` are swapped. To remedy this, we notice from :math:`\eqref{eq_dirac_field_u_spinor_zero_momentum}` and :math:`\eqref{eq_dirac_field_v_spinor_zero_momentum}` that the :math:`u` and :math:`v` spinors at zero-momentum are related by :math:`\gamma_5` defined by :math:`\eqref{eq_dirac_field_defn_gamma_5}`. Moreover, in order to cancel the :math:`\beta` matrices at the two ends of right-hand-side of :math:`\eqref{eq_dirac_field_d_matrix_conjugation_1}`, we can replace :math:`\pbf` with :math:`-\pbf`, which is also desirable as far as the time inversion is concerned in light of :math:`\eqref{eq_creation_operator_time_inversion_conjugation_massive}`. Putting all these considerations together, let's try the following

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

Construction of the interaction density
+++++++++++++++++++++++++++++++++++++++

As mentioned in :ref:`Dirac representation and gamma matrices <paragraph_dirac_field_representation_not_unitary>`, the fact that the Dirac representation is not unitary means that we cannot construct the interaction density using :math:`\psi^{\dagger} \psi` because it won't be a scalar. Indeed, let's work out how :math:`\psi^{\dagger}` transforms under a (homogeneous orthochronous) Lorentz transformation using :math:`\eqref{eq_dirac_field_pseudo_unitarity_of_d_matrix}` as follows

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

General fields
^^^^^^^^^^^^^^

We've now seen how scalar, vector, and Dirac fields can be constructed out of specific representations of the (homogeneous orthochronous) Lorentz group. These constructions can be generalized and unified by understanding the general representation theory of the Lorentz group.

General representation theory of the Lorentz group
++++++++++++++++++++++++++++++++++++++++++++++++++

The starting point, as in the case of Dirac fields, is the general commutation relation :math:`\eqref{eq_bracket_repr_j}` that the :math:`\Jscr_{\mu \nu}` matrices must satisfy. As explained in :ref:`sec_quantum_lorentz_symmetry`, we can rename the :math:`\Jscr_{\mu \nu}` matrices as follows

.. math::
	:nowrap:

	\begin{alignat*}{3}
		\Jscr_1 &\coloneqq \Jscr_{23}, \quad \Jscr_2 &&\coloneqq \Jscr_{31}, \quad \Jscr_3 &&\coloneqq \Jscr_{12} \\
		\Kscr_1 &\coloneqq \Jscr_{10}, \quad \Kscr_2 &&\coloneqq \Jscr_{20}, \quad \Kscr_3 &&\coloneqq \Jscr_{30}
	\end{alignat*}

and rewrite :math:`\eqref{eq_bracket_repr_j}` as a set of equations as follows

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
3. Dirac field: :math:`A = \tfrac{1}{2}, B = 0` or :math:`A = 0, B = \tfrac{1}{2}`. In either case :math:`\jfrak` must be :math:`\tfrac{1}{2}`. Indeed, they correspond to the two irreducible components of (the angular momentum part of) the Dirac representation :math:`\eqref{eq_dirac_field_jscr_matrix}`. Therefore the Dirac field may be written in short hand as :math:`\left( \tfrac{1}{2}, 0 \right) \oplus \left( 0, \tfrac{1}{2} \right)`.

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

As usual, let's translate :math:`\eqref{eq_j_intertwines_u}` and :math:`\eqref{eq_j_intertwines_v}` in the context of :math:`(A, B)` representations as follows

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

In light of :math:`\eqref{eq_dirac_field_linearize_representation}`, we can then write

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

	According to :math:`\eqref{eq_construct_interaction_density_by_fields}` and :math:`\eqref{eq_coefficient_g_transformation_law}`, a general interaction density can be constructed as follows

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

which recovers the cases of scalar field :math:`\eqref{eq_scalar_field_spatial_inversion_transformation_law}`, vector field :math:`\eqref{eq_vector_field_spatial_inversion_transformation_law}`, and Dirac field :math:`\eqref{eq_dirac_field_spatial_inversion_transformation_law}` where :math:`\beta`, as defined by :math:`\eqref{eq_dirac_field_beta_matrix}`, serves the function of swapping :math:`A` and :math:`B`.

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


The CPT theorem
^^^^^^^^^^^^^^^

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

Recall from :math:`\eqref{eq_evolution_equation_of_u_operator}` and :math:`\eqref{eq_defn_v_by_density}` that the interaction term :math:`V = \int d^3 x~\Hscr(0, \xbf)` satisfies the following [#volume_element_inversion_sign_cancels_for_cpt_symmetry]_

.. math::
	:nowrap:

	\begin{equation*}
		U(CPT) V U^{-1}(CPT) = -\int d^3 x~\Hscr(-x) = V
	\end{equation*}

Since the CPT symmetry is clearly conserved for free particles, it is also conserved in interactions according to :math:`\eqref{eq_h_as_h0_plus_v}`.


Massless fields
^^^^^^^^^^^^^^^

So far the story about quantum fields has been a 100% success. We've namely found the general formula :math:`\eqref{eq_general_field_psi_field}` for *any* field that represents a massive particle. However, such success will come to an end when we consider instead massless particles as we'll see in this section. This should not come as a surprise though since we've see in :math:`\eqref{eq_vector_field_defn_Pi}` for example, that the spin sum blows up in the massless limit :math:`m \to 0`.

Let's nonetheless kickstart the routine of constructing fields as follows, and see where the problem should arise.

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\ell}(x) = (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \left( \kappa e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma) a(\pbf, \sigma) + \lambda e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma) a^{c \dagger}(\pbf, \sigma) \right)
		\label{eq_massless_field_defn_psi_field}
	\end{equation}

This is reasonable because the translation symmetry is the same for massive and massless particles, and hence :math:`\eqref{eq_redefine_u_after_translation}` and :math:`\eqref{eq_redefine_v_after_translation}` apply.

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

Equating the coefficients of :math:`a(\pbf_{\Lambda, \sigma})` and :math:`a^{c \dagger}(\pbf_{\Lambda}, \sigma)` (i.e., the blue terms), and inverting :math:`D_{\ell \ell'}(\Lambda^{-1})` as in :math:`\eqref{eq_annihilation_u_transformation}`, we get the following conditions on the spinors

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

which is nothing but an incarnation of the mass-shell condition :math:`p_0^2 = \pbf^2`. Then let's consider the massless analog to the gauge-fixing condition :math:`\eqref{eq_vector_field_gauge_fixing_condition}`. To this end, we claim that :math:`e_0(\kbf, \pm 1) = 0` and :math:`\kbf \cdot \ebf(\kbf, \pm 1) = 0` imply the following

.. math::
	:nowrap:

	\begin{align}
		e_0(\pbf, \pm 1) &= 0
		\label{eq_massless_vector_field_spinor_zero_vanishes} \\
		\pbf \cdot \ebf(\pbf, \pm 1) &= 0
		\label{eq_massless_vector_field_spinor_orthogonal_to_momentum}
	\end{align}

in analogy to :math:`\eqref{eq_vector_field_spinor_orthogonal_to_momentum}` by the following argument. First, note that :math:`e_{\mu}(\pbf, \pm 1)` can be obtained from :math:`e_{\mu}(\kbf, \pm 1)` by applying :math:`L(p)` as in :math:`\eqref{eq_massless_vector_field_spinor_from_k_to_p}`. Second, :math:`L(p)` can be decomposed into a boost along the :math:`z`-axis as in :math:`\eqref{eq_massless_boost}` followed by a :math:`3`-rotation. Finally, we conclude :math:`\eqref{eq_massless_vector_field_spinor_zero_vanishes}` and :math:`\eqref{eq_massless_vector_field_spinor_orthogonal_to_momentum}` by noting that :math:`e_{\mu}(\kbf, \pm 1)` is unaffected by any boost along the :math:`z`, and the dot product is preserved by any :math:`3`-rotation. In contrast to the massless case :math:`\eqref{eq_vector_field_spinor_orthogonal_to_momentum}`, we have a stronger constraint :math:`\eqref{eq_massless_vector_field_spinor_zero_vanishes}` here because the helicity :math:`0` spinor is missing.

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
		D_{\mu}^{\nu}(W(a, b, \theta)) e_{\nu}(\kbf, \pm 1) &\xlongequal{\eqref{eq_massless_vector_field_w_equals_s_times_r}} S(a, b)_{\mu}^{\lambda} R(\theta)_{\lambda}^{\nu} e_{\nu}(\kbf, \pm 1)
		\label{eq_massless_vector_field_dw_acts_on_spinor} \\
			&\xlongequal{\eqref{eq_massless_field_spinor_e_condition_from_r}} e^{\pm \ifrak \theta} S(a, b)_{\mu}^{\lambda} e_{\lambda}(\kbf, \pm 1)
			\nonumber \\
			&\xlongequal{\eqref{eq_massless_little_group_s_matrix}} e^{\pm \ifrak \theta} \left( e_{\mu}(\kbf, \pm 1) + \frac{a \pm \ifrak b}{\sqrt{2}} k_{\mu} \right)
			\nonumber
	\end{align}

Next we recall the little group element :math:`W(\Lambda, p) = L^{-1}(\Lambda p) \Lambda L(p)` by definition. Plugging into :math:`\eqref{eq_massless_vector_field_dw_acts_on_spinor}`, we obtain the following

.. math::
	:nowrap:

	\begin{equation}
		\Lambda_{\mu}^{\nu} e_{\nu}(\pbf, \pm 1) = e^{\pm \ifrak \theta(\Lambda, p)} \left(e_{\mu}(\pbf_{\Lambda}, \pm 1) + (\Lambda p)_{\mu} \Omega_{\pm}(\Lambda, p)\right)
		\label{eq_massless_vector_field_spinor_e_condition}
	\end{equation}

where

.. math::
	:nowrap:

	\begin{equation*}
		\Omega_{\pm}(\Lambda, p) \coloneqq \frac{a(\Lambda, p) \pm \ifrak b(\Lambda, p)}{\sqrt{2}}
	\end{equation*}

is the extra term that makes it different from :math:`\eqref{eq_massless_field_spinor_u_condition}` which would have been satisfied if :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}` holds.

To most conveniently utilize :math:`\eqref{eq_massless_vector_field_spinor_e_condition}`, let's calculate :math:`D_{\mu}^{\nu} \left( U(\Lambda) a_{\nu}(x) U^{-1}(\Lambda) \right)` using :math:`\eqref{eq_massless_vector_field_a}, \eqref{eq_massless_vector_field_creation_operator_conjugated_by_lorentz_transformation}`, and :math:`\eqref{eq_massless_vector_field_annihilation_operator_conjugated_by_lorentz_transformation}` as follows

.. math::
	:nowrap:

	\begin{align*}
		D_{\mu}^{\nu}(\Lambda) \left(U(\Lambda) a_{\nu}(x) U^{-1}(\Lambda)\right) &= (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot x} \Lambda_{\mu}^{\nu} e_{\nu}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} e^{-\ifrak \sigma \theta(\Lambda, p)} a(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} e^{-\ifrak p \cdot x} \Lambda_{\mu}^{\nu} e^{\ast}_{\nu}(\pbf, \sigma) \sqrt{\frac{(\Lambda p)_0}{p_0}} e^{\ifrak \sigma \theta} a^{c \dagger}(\pbf_{\Lambda}, \sigma) \Big) \\
			&\xlongequal{\eqref{eq_massless_vector_field_spinor_e_condition}} (2\pi)^{-3/2} \int \frac{d^3 p}{p_0} \sqrt{\frac{(\Lambda p)_0}{2}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot x} e_{\mu}(\pbf_{\Lambda}, \sigma) a(\pbf_{\Lambda}, \sigma) + e^{-\ifrak p \cdot x} e^{\ast}_{\mu}(\pbf_{\Lambda}, \sigma) a^{c \dagger}(\pbf_{\Lambda}, \sigma) \phantom{\Big)} \\
			&\phantom{\xlongequal{\eqref{eq_massless_vector_field_spinor_e_condition}} \Big(} + (\Lambda p)_{\mu} \left( e^{\ifrak p \cdot x} \Omega_{\sigma}(\Lambda, p) a(\pbf_{\Lambda}, \sigma) + e^{-\ifrak p \cdot x} \Omega^{\ast}_{\sigma}(\Lambda, p) a^{c \dagger}(\pbf_{\Lambda}, \sigma) \right) \Big) \\
			&= (2\pi)^{-3/2} \int \frac{d^3 p}{\sqrt{2p_0}} \sum_{\sigma = \pm 1} \Big( e^{\ifrak p \cdot \Lambda x} e_{\mu}(\pbf, \pm 1) a(\pbf, \sigma) + e^{-\ifrak p \cdot \Lambda x} e^{\ast}_{\mu}(\pbf, \pm 1) a^{c \dagger}(\pbf, \sigma) \phantom{\Big)} \\
			&\phantom{= \Big(} + p_{\mu} \left( e^{\ifrak p \cdot \Lambda x} \Omega_{\sigma}(\Lambda, \Lambda^{-1} p) + e^{-\ifrak p \cdot \Lambda x} \Omega^{\ast}_{\sigma}(\Lambda, \Lambda^{-1} p) a^{c \dagger}(\pbf, \sigma) \right) \Big) \\
			&= a_{\mu}(\Lambda x) + \frac{\p}{\p (\Lambda_{\mu}^{\nu} x_{\nu})} \Omega(\Lambda, x)
	\end{align*}

where :math:`\Omega(\Lambda, x)` is a linear combination of creation and annihilation operators, whose precise form is not important here. Finally, moving :math:`D(\Lambda)` to the right-hand-side, we obtain the following variation of :math:`\eqref{eq_massless_field_psi_transform_by_d_matrix}`

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) a_{\mu}(x) U^{-1}(\Lambda) = D(\Lambda^{-1})_{\mu}^{\nu} a_{\nu}(\Lambda x) + \p_{\mu} \Omega(\Lambda, x)
		\label{eq_massless_vector_field_a_conjugation_by_lorentz_transformation}
	\end{equation}

which the massless vector field :math:`\eqref{eq_massless_vector_field_a}` actually satisfies.

It follows, using integration by parts, that one can construct interaction densities by coupling :math:`j^{\mu}(x) a_{\mu}(x)` as long as :math:`\p_{\mu} j^{\mu}(x) = 0`. This is completely parallel to :math:`\eqref{eq_vector_field_j_coupling}` and :math:`\eqref{eq_vector_field_j_coupling_condition}` for massive vector fields, and hence partially justifies the choice of the spinors in :math:`\eqref{eq_massless_vector_field_e_at_k}`, which satisfy :math:`\eqref{eq_massless_field_spinor_e_condition_from_r}` but not :math:`\eqref{eq_massless_field_spinor_e_condition_from_s}`.

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

It follows from :math:`\eqref{eq_dirac_field_linearize_representation}`, and the fact according to :math:`\eqref{eq_massless_little_group_r_matrix}` that :math:`R(\theta)` is the rotation in the :math:`xy`-plane, that

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


.. _sec_the_feynman_rules:

The Feynman Rules
-----------------

In :ref:`sec_cluster_decomposable_hamiltonians`, we've discussed the condition that the Hamiltonian must satisfy in order for the cluster decomposition principle to hold. Namely, the Hamiltonian can be written as a polynomial :math:`\eqref{eq_general_expansion_of_hamiltonian}` of creation and annihilation operators in normal order, such that the coefficients contains exactly one momentum conversing delta function. In order to derive this condition, we've encountered the idea of Feynman diagrams, which is a bookkeeping device for the evaluation of S-matrix. We couldn't say more about the coefficients besides the delta function because we had not introduced the building blocks of the Hamiltonian, namely, the quantum fields. Now that we've seen how to construct even the most general fields in the previous chapter, we're ready to spell out the full details of Feynman diagrams, first in spacetime coordinates, and then in momentum space coordinates.

.. warning::

	We'll consider only fields of massive particles in this chapter.


.. _sec_spacetime_feynman_rules:

Spacetime Feynman rules
^^^^^^^^^^^^^^^^^^^^^^^

In this section, we'll derive the Feynman rules in spacetime coordinates, which will establish a diagrammatic framework for calculating S-matrices. In the next section, we'll translate these rules to the momentum space.

First, recall the S-matrix formulated in terms of the S-operator :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered_density}` as follows

.. math::
	:nowrap:

	\begin{align}
		S_{\pbf'_1, \sigma'_1, n'_1;~\pbf'_2, \sigma'_2, n'_2;~\cdots,~\pbf_1, \sigma_1, n_1;~\pbf_2, \sigma_2, n_2;~\cdots} &= \sum_{N=0}^{\infty} \frac{(-\ifrak)^N}{N!} \int d^4 x_1 d^4 x_2 \cdots d^4 x_n \Big( \Phi_{\VAC} \cdots a(\pbf'_2, \sigma'_2, n'_2) a(\pbf'_1, \sigma'_1, n'_2) \phantom{\Big)}
		\label{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density} \\
			&\phantom{= \Big(} \times T\left\{ \Hscr(x_1) \cdots \Hscr(x_N) \right\} a^{\dagger}(\pbf_1, \sigma_1, n_1) a^{\dagger}(\pbf_2, \sigma_2, n_2) \cdots \Phi_{\VAC} \Big) \nonumber
	\end{align}

where :math:`T\{ \cdots \}` is the time-ordered product defined by :math:`\eqref{eq_defn_time_ordered_product}`. As we've seen from the previous chapter, the interaction density :math:`\Hscr(x)` may be written as a polynomial in fields and their adjoint as follows

.. math::
	:nowrap:

	\begin{equation}
		\Hscr(x) = \sum_{i} g_i \Hscr_i(x)
		\label{eq_interaction_density_as_sum_of_monomials}
	\end{equation}

where :math:`\Hscr_i(x)` is monomial of certain fields and their adjoint. Finally, we recall the general formula :math:`\eqref{eq_general_field_psi_field}` for a quantum field as follows

.. math::
	:nowrap:

	\begin{equation}
		\psi_{\ell}(x) = (2\pi)^{-3/2} \sum_{\sigma} \int d^3 p \left( e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n) a(\pbf, \sigma, n) + e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n) a^{\dagger}(\pbf, \sigma, n^c) \right)
		\label{eq_generic_field_expression}
	\end{equation}

where we've restored the particle species index :math:`n`, and absorbed the sign :math:`(-1)^{2B}` in :math:`\eqref{eq_general_field_psi_field}` into the :math:`v`-spinor. Unlike the notation used in :ref:`sec_quantum_lorentz_symmetry`, the sub-index :math:`\ell` here includes not only the running indexes of the Lorentz representation, but also the representation itself, as well as the particle species :math:`n`.

.. warning::

	Any specific field of interest, whether it's scalar, vector, or Dirac, must be first put into the form :math:`\eqref{eq_generic_field_expression}` in order to use the Feynman rules, that we'll introduce below.

According to :math:`\eqref{eq_defn_annihilation_field}` and :math:`\eqref{eq_defn_creation_field}`, we can also write, with obvious modifications, :math:`\psi_{\ell}(x) = \psi^+_{\ell}(x) + \psi^-_{\ell}(x)` such as :math:`\psi^+_{\ell}(x)` is a linear combination of annihilation operators, and :math:`\psi^-_{\ell}(x)` is a linear combination of creation operators.

Now the idea of the Feynman rules to calculate the S-matrix is same as what has been discussed in :ref:`sec_cluster_decomposable_hamiltonians`. Namely, we'd like to move any annihilation operator to the right of a creation operator using the standard commutation rule :math:`\eqref{eq_creation_annihilation_commutator}`. To be more specific, we'll list all the possible scenarios as follows

.. _listing_feynman_rules:

#. Paring a final particle (in out-state) :math:`(\pbf, \sigma, n)` with a field adjoint :math:`\psi^{\dagger}_{\ell}(x)` gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{a(\pbf, \sigma, n) \psi^{\dagger}_{\ell}(x)}
				\coloneqq \left[ a(\pbf, \sigma, n), \psi^{\dagger}_{\ell}(x) \right]_{\pm}
				= (2\pi)^{-3/2} e^{-\ifrak p \cdot x} u^{\ast}_{\ell}(\pbf, \sigma, n)
			\label{eq_feynman_rule_a_psi_dagger}
		\end{equation}

#. Paring a final antiparticle :math:`(\pbf, \sigma, n^c)` with a field :math:`\psi_{\ell}(x)` gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{a(\pbf, \sigma, n^c) \psi_{\ell}(x)}
				\coloneqq \left[ a(\pbf, \sigma, n^c), \psi_{\ell}(x) \right]_{\pm}
				= (2\pi)^{-3/2} e^{-\ifrak p \cdot x} v_{\ell}(\pbf, \sigma, n)
			\label{eq_feynman_rule_a_psi}
		\end{equation}

#. Paring a field :math:`\psi_{\ell}(x)` with an initial particle (in in-state) :math:`(\pbf, \sigma, n)` gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{\psi_{\ell}(x) a^{\dagger}(\pbf, \sigma, n)}
			\coloneqq \left[ \psi_{\ell}(x), a^{\dagger}(\pbf, \sigma, n) \right]_{\pm}
			= (2\pi)^{-3/2} e^{\ifrak p \cdot x} u_{\ell}(\pbf, \sigma, n)
			\label{eq_feynman_rule_psi_a_dagger}
		\end{equation}

#. Paring a field adjoint :math:`\psi^{\dagger}_{\ell}(x)` with an initial antiparticle :math:`(\pbf, \sigma, n^c)` gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{\psi^{\dagger}(x) a^{\dagger}(\pbf, \sigma, n^c)}
			\coloneqq \left[ \psi^{\dagger}(x), a^{\dagger}(\pbf, \sigma, n^c) \right]_{\pm}
			= (2\pi)^{-3/2} e^{\ifrak p \cdot x} v^{\ast}_{\ell}(\pbf, \sigma, n)
			\label{eq_feynman_rule_psi_dagger_a_dagger}
		\end{equation}

#. Paring a final particle :math:`(\pbf, \sigma, n)` (or antiparticle) with an initial particle :math:`(\pbf', \sigma', n')` (or antiparticle) gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{a(\pbf', \sigma', n') a^{\dagger}(\pbf, \sigma, n)}
			\coloneqq \left[ a(\pbf', \sigma', n'), a^{\dagger}(\pbf, \sigma, n) \right]_{\pm}
			= \delta^3(\pbf' - \pbf) \delta_{\sigma' \sigma} \delta_{n' n}
			\label{eq_feynman_rule_a_a_dagger}
		\end{equation}

#. Paring a field :math:`\psi_{\ell}(x)` in :math:`\Hscr_i(x)` with a field adjoint :math:`\psi_m^{\dagger}(y)` in :math:`\Hscr_j(y)` gives
	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{\psi_{\ell}(x) \psi^{\dagger}_m(y)}
			\coloneqq \theta(x_0 - y_0) \left[ \psi^+_{\ell}(x), \psi^{+ \dagger}_m(y) \right]_{\pm} \mp \theta(y_0 - x_0) \left[ \psi^{- \dagger}_m(y), \psi^-_{\ell}(x) \right]_{\pm} \eqqcolon -\ifrak \Delta_{\ell m}(x, y)
			\label{eq_feynman_rule_propagator}
		\end{equation}

   where :math:`\theta(\tau)` is the step function which equals :math:`1` for :math:`\tau > 0` and vanishes for :math:`\tau < 0`. Here we remind ourselves once again that the Feynman rule is all about moving annihilation operators, e.g. :math:`\psi^+_{\ell}(x)` and :math:`\psi^{- \dagger}_m(y)`, to the right of creation operators, e.g. :math:`\psi^-_{\ell}(x)` and :math:`\psi^{- \dagger}_m(y)`. The sign :math:`\mp` in the middle is due to the fact that when the top sign should to be used, the particles are fermions, and hence the interchange of the fields due to time ordering requires an extra minus sign.

   This quantity is known as a *propagator*, which will be evaluated in the next section.

.. note::
	1. In the above listing, we've assumed that the item on the left in the (anti-)commutators also lies to the left of the item on the right in the (anti-)commutator in the vacuum expectation value in :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}` if we ignore the time ordering operator. The same applies to case (6) where :math:`\psi_{\ell}(x)` is assumed to lie to the left of :math:`\psi_m^{\dagger}(y)`. In particular, if we assume that the interaction density :math:`\Hscr(x)` is normally ordered in the sense that all the field adjoints lie to the left of the fields, then :math:`\Hscr_i(x)` necessarily lies to the left of :math:`\Hscr_j(y)`, and we don't have to to define :math:`\Delta_{\ell m}(x, x)` which would require certain regulation to not blow up integrals.
	2. The parings listed above are commonly known as `Wick contractions <https://en.wikipedia.org/wiki/Wick%27s_theorem>`__.

A great invention of Feynman is the following diagrammatic representation of the above rules, known as the Feynman diagrams.

.. _fig_spacetime_feynman_diagrams:

.. figure:: ./static/quantum-theory-of-fields/space-propagators.svg
	:align: center

	Figure: All possible edges in Feynman diagrams.

A few comments are in order to clarify the meaning of these diagrams

* The arrow points towards the (positive) direction of time, which is upwards for particles and downwards for antiparticles. In other words, we interpret an antiparticle as a particle that moves backwards in time. An exceptional case is (6), where the edge is placed horizontally. The reason is that a field or its adjoint doesn't just create or destroy (anti-)particles -- they create/destroy a particle and at the same time destroy/create the corresponding antiparticle, respectively. In light of :math:`\eqref{eq_feynman_rule_propagator}`, there is no reason to prefer either an upward or a downward arrow.
* The arrow in (6) points from :math:`(m, y)` to :math:`(\ell, x)` since :math:`\psi_{\ell}(x)` is a field and :math:`\psi^{\dagger}_m(y)` is a field adjoint. Two processes happen in this scenario, namely, a particle created by :math:`\psi^{+ \dagger}_m(y)` is absorbed by :math:`\psi^+_{\ell}(x)`, and an antiparticle created by :math:`\psi^-_{\ell}(x)` is absorbed by :math:`\psi^{- \dagger}_m(y)`. The arrow is compatible with both processes.
* In the case where the particle is its own antiparticle, the arrows in (1) -- (6) will be omitted because one cannot tell apart a field and a field adjoint according to :math:`\eqref{eq_general_field_charge_inversion_transformation}`.
* We didn't draw the other scenario in (5) where an antiparticle is created and then destroyed without any interaction. In this case we need to flip the direction of the arrow.
* Every nodes in the diagram, marked by a fat dot, correspond to a monomial :math:`\Hscr_i(x)` in :math:`\eqref{eq_interaction_density_as_sum_of_monomials}`. Moreover, for each node, there are as many incoming edges as there are fields, and as many outgoing edges as there are field adjoints.

With these basic building blocks at hand, we're ready to evaluate :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}` using Feynman diagrams. The first thing to notice is that the S-matrix, as defined by :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}`, should be viewed as a power series, where a term of order :math:`N` corresponds to a monomial given as a product of powers, each of which is :math:`N_i`-th power of an interaction type :math:`g_i \Hscr_i(x)` (see :math:`\eqref{eq_interaction_density_as_sum_of_monomials}`), such that :math:`N = \sum_i N_i`.

To a given term of order :math:`N`, one can draw the associated Feynman diagrams in two steps as follows. The first step is to draw (on a piece of paper) one upward-pointing and downward-pointing strand for each in-state particle and antiparticle, respectively, at the bottom; do the same to the out-state particles and antiparticles at the top; and draw in the middle :math:`N_i` vertices for each interaction types :math:`i`, which corresponds to a monomial :math:`\Hscr_i(x)`, with incoming and outgoing edges corresponding to its fields and field adjoints, respectively. The second step is to connect any pair of open strands by an oriented edge if the particle in question is different from its antiparticle, and by an unoriented edge otherwise. Moreover, the vertices and the edges are labelled as illustrated in the figure above.

Now knowing how to draw a Feynman diagram of any given order for an interaction, we can spell out the recipe for a diagrammatic calculation of the S-matrix given by :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}` in the following steps.

1. Draw all (distinct) Feynman diagrams (in reality, to a finite order) following the rules described above. We'll come back to what we mean by "distinct" right after the recipe.
2. For each diagram, we assign a factor :math:`-\ifrak` to each vertex, corresponding to one factor of the power :math:`(-\ifrak)^n` in :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}`; and a factor :math:`g_i` to each vertex corresponding to the coefficient of :math:`\Hscr_i(x)` in :math:`\eqref{eq_interaction_density_as_sum_of_monomials}`; and a factor in the :ref:`listing of Feynman rules <listing_feynman_rules>` to each edge. Multiplying all the factors together and integrating over all the coordinates :math:`x_1, x_2, \cdots`, one for each vertex, we obtain a (numeric) valuation of the Feynman diagram.
3. The S-matrix is the "sum" over all the evaluations of the Feynman diagrams. Here the sum is put in quotation marks because we might do subtraction instead of addition when there are fermionic fields involved in the interaction. More precisely, for each Feynman diagram, one can move the fields and field adjoints over each other that the two ends of every edge are next to each other (in the right order). Then we add a minus sign in front of the evaluation if such rearrangement involves an odd number of swaps between fermionic fields (or field adjoints).

These are really all one needs to evaluate S-matrices using Feynman diagrams, but a few further notes may be necessary to make it completely clear. Firstly, note that we've ignored the factor :math:`1 / N!` in :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}` from our recipe above. The reason lies in the word "distinct" from the first step. More precisely, we consider two Feynman diagrams, which differ by a re-labeling of the vertices, to be the same. Since there are :math:`N!` ways of labeling :math:`N` vertices, we have already taken the fraction :math:`1 / N!` into account by including only distinct diagrams.

Secondly, by the discussion in :ref:`sec_cluster_decomposable_hamiltonians`, only connected Feynman diagrams will be included so that the resulting S-matrix satisfies the cluster decomposition principle.

.. _paragraph_interaction_density_symmetry_factor:

The last note is more of a convention (for convenience), which aims at further remove duplications among Feynman diagrams. Here by duplication we mean diagrams that are not exactly the same but whose evaluations are the same. A basic example is when a monomial in the interaction density :math:`\Hscr(x)` contains a power of the same field (or field adjoint). From the viewpoint of Feynman diagrams, it means that a vertex may have more than one identical attached (incoming or outgoing) strands. Now when other strands want to connect to these identical ones, they may choose which one to connect first, and second, and so on, but the result will be the same regardless of the choices. Hence it's a convention to write the coefficient of :math:`\Hscr_i(x)` as :math:`g_i / k!` if it contains :math:`k` identical fields (or field adjoints), so that in a diagrammatic calculation, one only need to include one of such diagrams. Other numerical factors might be inserted in more complex situations, such as when two vertices with identical attached strands try to connect to each other, or when there is a loop of identical vertices. These cases are discussed in [Wei95]_ page 265 -- 267, and we'll come back to them when they become relevant in calculations.

To make things concrete and to prepare for the calculations in the next sections, we conclude the discussion of Feynman rules with two prototypical examples.

Example 1: :math:`\psi^{\dagger} \psi \phi`-interaction
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Consider the following interaction density

.. math::
	:nowrap:

	\begin{equation}
		\Hscr(x) = \sum_{\ell m k} g_{\ell m k} \psi_{\ell}^{\dagger}(x) \psi_m(x) \phi_k(x)
		\label{eq_psi_dagger_psi_phi_interaction_density}
	\end{equation}

where :math:`\psi(x)` is a (complex) fermionic field, and :math:`\phi(x)` is a real, i.e., :math:`\phi(x) = \phi^{\dagger}(x)`, bosonic field. This is an important type of interaction since it shows up not only in quantum electrodynamics, but in fact in the whole Standard Model of electromagnetic, weak, and strong interactions.

This type of interaction allows three kinds of scattering processes, namely, fermion-fermion, fermion-boson, and boson-boson scattering, which we'll discuss one by one.

.. _listing_fermion_fermion_scattering:

Fermion-fermion scattering
	The scattering is represented as :math:`12 \to 1'2'`, where all in- and out-state particles :math:`1, 2, 1', 2'` are fermions. Up to the second order, there are two (connected) Feynman diagrams

	.. figure:: ./static/quantum-theory-of-fields/psi-dagger-psi-phi-fermion-fermion-scattering.svg
		:align: center

		Figure: Two second-order fermion-fermion scattering diagrams :math:`\psi^{\dagger} \psi \phi`-theory.

	where the solid (directed) edges represent the fermions, and the dashed (undirected) edges represent the (neutral) boson. More explicitly, the two diagrams correspond to the following two contractions

	.. math::
		:nowrap:

		\begin{align}
			\underbracket{ a(1') \psi_{\ell}^{\dagger}(x) }
			\underbracket{ a(2') \psi_{\ell'}^{\dagger}(y) }
			\underbracket{ \phi_k(x) \phi_{k'}(y) }
			\underbracket{ \psi_m(x) a^{\dagger}(1) }
			\underbracket{ \psi_{m'}(y) a^{\dagger}(2) }
			\label{eq_fermion_fermion_scattering_1212} \\

			\underbracket{ a(2') \psi_{\ell}^{\dagger}(x) }
			\underbracket{ a(1') \psi_{\ell'}^{\dagger}(y) }
			\underbracket{ \phi_k(x) \phi_{k'}(y) }
			\underbracket{ \psi_m(x) a^{\dagger}(1) }
			\underbracket{ \psi_{m'}(y) a^{\dagger}(2) }
			\label{eq_fermion_fermion_scattering_1221}
		\end{align}

	respectively. Moreover, comparing with the original order

	.. math::
		:nowrap:

		\begin{equation}
			a(2') a(1') \psi^{\dagger}(x) \psi(x) \phi(x) \psi^{\dagger}(y) \psi(y) \phi(y) a^{\dagger}(1) a^{\dagger}(2)
			\label{eq_two_particles_scattering_original_order}
		\end{equation}

	we see that the contractions :math:`\eqref{eq_fermion_fermion_scattering_1212}` and :math:`\eqref{eq_fermion_fermion_scattering_1221}` require an even and odd number of swaps between fermionic operators, respectively. We note that whether a given diagram requires an even or odd fermionic swaps is rather arbitrary, and depends on many conventions. However, the fact that the two diagrams corresponding to :math:`\eqref{eq_fermion_fermion_scattering_1212}` and :math:`\eqref{eq_fermion_fermion_scattering_1221}` carry opposite signs is independent of the conventions and hence meaningful. Indeed, it's another incarnation of the Fermi statistics in the sense that the S-matrix switches sign if either the in-state fermions :math:`1 \leftrightarrow 2` or the out-state fermions :math:`1' \leftrightarrow 2'` are swapped.

	Now let's use the Feynman rules :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_propagator}` to evaluate the fermion-fermion scattering S-matrix up to second order as follows [#connected_s_matrix_by_feynman_diagrams]_

	.. math::
		:nowrap:

		\begin{align*}
			S^C_{\pbf'_1 \sigma'_1 n'_1,~\pbf'_2 \sigma'_2 n'_2;~~\pbf_1 \sigma_1 n_1,~\pbf_2 \sigma_2 n_2} &=
					\sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y
					\big( \eqref{eq_fermion_fermion_scattering_1212} - \eqref{eq_fermion_fermion_scattering_1221} \big) \\
				&= (2\pi)^{-6} \sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y~(-\ifrak) \Delta_{k k'}(x, y) \times \\
				&\phantom{=} \times e^{\ifrak p_1 \cdot x + \ifrak p_2 \cdot y} u_m(\pbf_1, \sigma_1, n_1) u_{m'}(\pbf_2, \sigma_2, n_2) \\
				&\phantom{=} \times \left( e^{-\ifrak p'_1 \cdot x - \ifrak p'_2 \cdot y} u^{\ast}_{\ell}(\pbf'_1, \sigma'_1, n'_1) u^{\ast}_{\ell'}(\pbf'_2, \sigma'_2, n'_2)
	  				- e^{-\ifrak p'_2 \cdot x - \ifrak p'_1 \cdot y} u^{\ast}_{\ell}(\pbf'_2, \sigma'_2, n'_2) u^{\ast}_{\ell'}(\pbf'_1, \sigma'_1, n'_1) \right)
		\end{align*}

.. _listing_fermion_boson_scattering:

Fermion-boson scattering
	Consider the scattering :math:`12 \to 1'2'`, where particles :math:`1, 1'` are fermions and :math:`2, 2'` are bosons, under interaction density :math:`\eqref{eq_psi_dagger_psi_phi_interaction_density}`. Up to second order, there are again two Feynman diagrams as follows

	.. figure:: ./static/quantum-theory-of-fields/psi-dagger-psi-phi-fermion-boson-scattering.svg
		:align: center

		Figure: Two second-order fermion-boson scattering diagrams :math:`\psi^{\dagger} \psi \phi`-theory.

	They correspond to the following two contractions

	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{ a(2') \phi_k(x) }
			\underbracket{ a(1') \psi^{\dagger}_{\ell}(x) }
			\underbracket{ \psi_m(x) \psi_{\ell'}^{\dagger}(y) }
			\underbracket{ \psi_{m'}(y) a^{\dagger}(1) }
			\underbracket{ \phi_{k'}(y) a^{\dagger}(2) }
			\label{eq_fermion_boson_scattering_1212}
		\end{equation}

	and

	.. math::
		:nowrap:

		\begin{equation}
			\underbracket{ a(2') \phi_{k'}(y) }
			\underbracket{ a(1') \psi_{\ell}^{\dagger}(x) }
			\underbracket{ \psi_m(x) \psi_{\ell'}^{\dagger}(y) }
			\underbracket{ \psi_{m'}(y) a^{\dagger}(1) }
			\underbracket{ \phi_k(x) a^{\dagger}(2) }
			\label{eq_fermion_boson_scattering_1221}
		\end{equation}

	respectively. Comparing with the ordering of operators :math:`\eqref{eq_two_particles_scattering_original_order}`, we see that neither :math:`\eqref{eq_fermion_boson_scattering_1212}` nor :math:`\eqref{eq_fermion_boson_scattering_1221}` require any fermionic swap, and hence no extra signs are need in this case, in contrast to the previous case of fermion-fermion scattering.

	Next let's use :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_propagator}` to evaluate the second order S-matrix as follows

	.. math::
		:nowrap:

		\begin{align*}
			S^C_{\pbf'_1 \sigma'_1 n'_1,~\pbf'_2 \sigma'_2 n'_2;~\pbf_1 \sigma_1 n_1,~\pbf_2 \sigma_2 n_2} &=
					\sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y \
					\big( \eqref{eq_fermion_boson_scattering_1212} + \eqref{eq_fermion_boson_scattering_1221} \big) \\
				&= (2\pi)^{-6} \sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y~(-\ifrak) \Delta_{m \ell'}(x, y) \times \\
				&\phantom{=} \times e^{-\ifrak p'_1 \cdot x + \ifrak p_1 \cdot y} u^{\ast}_{\ell}(\pbf'_1, \sigma'_1, n'_1) u_{m'}(\pbf_1, \sigma_1, n_1) \\
				&\phantom{=} \times \left( e^{-\ifrak p'_2 \cdot x + \ifrak p_2 \cdot y} u^{\ast}_k(\pbf'_2, \sigma'_2, n'_2) u_{k'}(\pbf_2, \sigma_2, n_2) \
					+ e^{-\ifrak p'_2 \cdot y + \ifrak p_2 \cdot x} u^{\ast}_{k'}(\pbf'_2, \sigma'_2, n'_2) u_k(\pbf_2, \sigma_2, n_2) \right)
		\end{align*}

.. _listing_boson_boson_scattering:

Boson-boson scattering
	It turns out that the lowest order boson-boson scattering under the interaction density :math:`\eqref{eq_psi_dagger_psi_phi_interaction_density}` is four. An example of such scattering is given by the following Feynman diagram

	.. figure:: ./static/quantum-theory-of-fields/psi-dagger-psi-phi-boson-boson-scattering.svg
		:align: center

		Figure: The fourth order boson-boson scattering diagram in :math:`\psi^{\dagger} \psi \phi`-theory.

	which involves a fermionic loop. We'll not evaluate the corresponding S-matrix here, but note that the corresponding (fermionic) contraction

	.. math::
		:nowrap:

		\begin{equation*}
			\underbracket{ \psi(x_1) \psi^{\dagger}(x_2) }
			\underbracket{ \psi(x_2) \psi^{\dagger}(x_3) }
			\underbracket{ \psi(x_3) \psi^{\dagger}(x_4) }
			\underbracket{ \psi(x_4) \psi^{\dagger}(x_1) }
		\end{equation*}

	where we've ignored the terms involving the bosonic operators, requires an odd number of fermionic swaps, which, in turn, requires an extra minus sign. This is a general phenomenon for any diagram that involves a fermionic loop.


Example 2: :math:`\phi^3`-interaction
+++++++++++++++++++++++++++++++++++++

Now let's consider an interaction density that involves a power of the same field as follows

.. math::
	:nowrap:

	\begin{equation}
		\Hscr(x) = \frac{1}{3!} \sum_{\ell m k} g_{\ell m k} \phi_{\ell}(x) \phi_m(x) \phi_k(x)
		\label{eq_phi3_interaction_density}
	\end{equation}

where :math:`g_{\ell m k}` is totally symmetric, and the extra factor :math:`1/3!` is to :ref:`account for this symmetry <paragraph_interaction_density_symmetry_factor>`.

We'll evaluate the scattering amplitude of two (identical) bosons :math:`12 \to 1'2'`. In this case, there are three connected Feynman diagrams of second order as follows

.. figure:: ./static/quantum-theory-of-fields/phi3-boson-boson-scattering.svg
	:align: center

	Figure: Three second-order boson-boson scattering diagrams in :math:`\phi^3`-theory.

They correspond to the following three contractions

.. math::
	:nowrap:

	\begin{align}
		&\underbracket{ a(2') \phi_{\ell}(x) }
		\underbracket{ a(1') \phi_{\ell'}(y) }
		\underbracket{ \phi_m(x) \phi_{m'}(y) }
		\underbracket{ \phi_{k'}(y) a^{\dagger}(1) }
		\underbracket{ \phi_k(x) a^{\dagger}(2) }
		\label{eq_phi3_h1212} \\

		&\underbracket{ a(2') \phi_{\ell'}(y) }
		\underbracket{ a(1') \phi_{\ell}(x) }
		\underbracket{ \phi_m(x) \phi_{m'}(y) }
		\underbracket{ \phi_{k'}(y) a^{\dagger}(1) }
		\underbracket{ \phi_k(x) a^{\dagger}(2) }
		\label{eq_phi3_h1221} \\

		&\underbracket{ a(2') \phi_{\ell}(x) }
		\underbracket{ a(1') \phi_k (x)}
		\underbracket{ \phi_m(x) \phi_{m'}(y) }
		\underbracket{ \phi_{k'}(y) a^{\dagger}(1) }
		\underbracket{ \phi_{\ell'}(y) a^{\dagger}(2) }
		\label{eq_phi3_v1212}
	\end{align}

respectively. We note once again that the factor :math:`1/3!` in :math:`\eqref{eq_phi3_interaction_density}` allows us to include just the above three contractions, and ignore the other ones obtained by, say, permuting :math:`\ell, m`, and :math:`k` (or equivalently :math:`\ell', m'`, and :math:`k'`). With this in mind, we can now calculate the (connected) S-matrix using :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_propagator}` as follows

.. math::
	:nowrap:

	\begin{align*}
		S^C_{\pbf'_1 \sigma'_1 n'_1,~\pbf'_2 \sigma'_2 n'_2;~\pbf_1 \sigma_1 n_1,~\pbf_2 \sigma_2 n_2} &=
				\sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y \
				\big( \eqref{eq_phi3_h1212} + \eqref{eq_phi3_h1221} + \eqref{eq_phi3_v1212} \big) \\
			&= (2\pi)^{-6} \sum_{\ell m k, \ell' m' k'} (-\ifrak)^2 g_{\ell m k} g_{\ell' m' k'} \int d^4 x d^4 y \
				(-\ifrak) \Delta_{m m'}(x, y) \times \\
			&\phantom{=} \times \Big( e^{\ifrak \left( p_2 - p'_2 \right) \cdot x + \ifrak \left( p_1 - p'_1 \right) \cdot y} u^{\ast}_{\ell}(\pbf'_2, \sigma'_2, n'_2) u^{\ast}_{\ell'}(\pbf'_1, \sigma'_1, n'_1) u_{k'}(\pbf_1, \sigma_1, n_1) u_k(\pbf_2, \sigma_2, n_2) \\
			&\phantom{= \times} + e^{\ifrak (p_2 - p'_1) \cdot x + \ifrak (p_1 - p'_2) \cdot y} u^{\ast}_{\ell'}(\pbf'_2, \sigma'_2, n'_2) u^{\ast}_{\ell}(\pbf'_1, \sigma'_1, n'_1) u_{k'}(\pbf_1, \sigma_1, n_1) u_k(\pbf_2, \sigma_2, n_2) \\
			&\phantom{= \times} + e^{-\ifrak (p'_1 + p'_2) \cdot x + \ifrak (p_1 + p_2) \cdot y} u^{\ast}_{\ell}(\pbf'_2, \sigma'_2, n'_2) u^{\ast}_k(\pbf'_1, \sigma'_1, n'_1) u_{k'}(\pbf_1, \sigma_1, n_1) u_{\ell'}(\pbf_2, \sigma_2, n_2) \Big)
	\end{align*}

A special case is when :math:`\phi` is a scalar field so that :math:`\eqref{eq_phi3_interaction_density}` takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		\Hscr(x) = \frac{g}{3!} \phi^3(x)
	\end{equation*}

In this case, the S-matrix can be simplified as follows

.. math::
	:nowrap:

	\begin{align*}
		S^C_{\pbf'_1,~\pbf'_2;~\pbf_1,~\pbf_2} &= \frac{\ifrak g^2}{(2\pi)^6 \sqrt{16 E'_1 E'_2 E_1 E_2}} \int d^4 x d^4 y~\Delta_F(x, y) \times \\
			&\phantom{=} \times \left( e^{\ifrak \left( p_2 - p'_2 \right) \cdot x + \ifrak \left( p_1 - p'_1 \right) \cdot y}
			+ e^{\ifrak (p_2 - p'_1) \cdot x + \ifrak (p_1 - p'_2) \cdot y}
			+ e^{-\ifrak (p'_1 + p'_2) \cdot x + \ifrak (p_1 + p_2) \cdot y} \right)
	\end{align*}

where :math:`\Delta_F` is the so-called `Feynman propagator <https://en.wikipedia.org/wiki/Propagator#Feynman_propagator>`__ and will be discussed in detail in the next section.


Momentum space Feynman rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There turns out to be many advantages in working in the :math:`4`-momentum space, instead of the spacetime. Hence the goal of this section is to translate the Feynman rules :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_propagator}` to the momentum space. The biggest challenge, however, lies in the fact that so far we've been working exclusively under the assumption that the :math:`4`-momentum lies on the mass shell, whether in the spinors appearing in :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_psi_a_dagger}` or in the propagator :math:`\eqref{eq_feynman_rule_propagator}`. We'll tackle this challenge in three steps. First, we'll rewrite the propagator as an integral on the momentum space, as opposed to the integral :math:`\eqref{eq_defn_Delta_plus}` defined on the mass shell. This will then allow us to translate the Feynman rules to the momentum space, except for the external edges which still live on the mass shell because the in- and out-state particles do so. Finally, we'll discuss how to generalize the Feynman rules so that the "external lines" do not necessarily have to live on the mass shell.

.. _sec_propagator_in_momentum_space:

Propagator in momentum space
++++++++++++++++++++++++++++

Plugging :math:`\eqref{eq_generic_field_expression}` into :math:`\eqref{eq_feynman_rule_propagator}`, we can rewrite the propagator in terms of the spin sums as follows

.. math::
	:nowrap:

	\begin{align}
		-\ifrak \Delta_{\ell m}(x, y) &= \theta(x_0 - y_0) (2\pi)^{-3} \int d^3 p~e^{\ifrak p \cdot (x-y)} \sum_{\sigma} u_{\ell}(\pbf, \sigma, n) u^{\ast}_m(\pbf, \sigma, n)
		\label{eq_propagator_as_spin_sums} \\
			&\phantom{=} \pm \theta(y_0 - x_0) (2\pi)^{-3} \int d^3 p~e^{\ifrak p \cdot (y-x)} \sum_{\sigma} v^{\ast}_m(\pbf, \sigma, n) v_{\ell}(\pbf, \sigma, n)
			\nonumber
	\end{align}

where the top and bottom signs apply to bosonic and fermionic fields, respectively.

Then, in light of :math:`\eqref{eq_general_field_spin_sums_as_pi}` and :math:`\eqref{eq_general_field_spin_sum_is_polynomial}`, we can write the spin sums as polynomials in (on-mass-shell) :math:`4`-momenta as follows

.. math::
	:nowrap:

	\begin{align}
		\sum_{\sigma} u_{\ell}(\pbf, \sigma, n) u^{\ast}_m(\pbf, \sigma, n) &= \left( 2\sqrt{\pbf^2 + M_n^2} \right)^{-1} P_{\ell m} \left( \sqrt{\pbf^2 + M_n^2}, \pbf \right)
		\label{eq_spin_sum_u_as_polynomial} \\
		\sum_{\sigma} v_{\ell}(\pbf, \sigma, n) v^{\ast}_m(\pbf, \sigma, n) &= \pm \left( 2\sqrt{\pbf^2 + M_n^2} \right)^{-1} P_{\ell m} \left( -\sqrt{\pbf^2 + M_n^2}, -\pbf \right)
		\label{eq_spin_sum_v_as_polynomial}
	\end{align}

where the top and bottom signs correspond to bosonic and fermionic fields, respectively. Note that although according to :math:`\eqref{eq_general_field_spin_sums_as_pi}`, the spin sums for :math:`u` and :math:`v` can be made the same (with respect to some basis of representation), we've introduced here an additional sign for :math:`v` using the symmetry property :math:`\eqref{eq_general_field_spin_sum_polynomial_parity}`. In particular, as we'll see right below, this extra sign is also needed to be consistent with our previous calculations for Dirac fields.

Before we continue evaluating the propagator, let's take a moment to figure out how :math:`P_{\ell m}(p)` looks like in the cases of scalar, vector, Dirac, and general :math:`(A, B)` fields. The case of scalar fields is the simplest, and it follows from :math:`\eqref{eq_scalar_u_and_v}` that

.. math::
	:nowrap:

	\begin{equation}
		P^{\text{scalar}}(p) = 1
		\label{eq_p_polynomial_scalar}
	\end{equation}

Next, for spin-:math:`1` vector fields, it follows from :math:`\eqref{eq_vector_field_Pi_matrix}` and :math:`\eqref{eq_vector_field_defn_Pi}` that

.. math::
	:nowrap:

	\begin{equation}
		P^{\text{vector}}_{\mu \nu}(p) = \eta_{\mu \nu} + p_{\mu} p_{\nu} / M^2
		\label{eq_p_polynomial_vector}
	\end{equation}

which is even in :math:`p`, and hence consistent with :math:`\eqref{eq_spin_sum_v_as_polynomial}`.

Then for spin-:math:`1/2` Dirac fields, it follows from :math:`\eqref{eq_dirac_field_n_matrix_as_spinor_sum}` -- :math:`\eqref{eq_dirac_field_m_matrix_as_spinor_sum}` and :math:`\eqref{eq_dirac_field_spin_sum_u}` -- :math:`\eqref{eq_dirac_field_spin_sum_v}` that

.. math::
	:nowrap:

	\begin{equation}
		P^{\text{Dirac}}_{\ell m}(p) = \big( (-\ifrak p^{\mu} \gamma_{\mu} + M) \beta \big)_{\ell m}
		\label{eq_p_polynomial_dirac}
	\end{equation}

We see that :math:`\eqref{eq_dirac_field_spin_sum_v}` is consistent with :math:`\eqref{eq_spin_sum_v_as_polynomial}`, thanks to the sign convention.

Finally, for any :math:`(A, B)` field, the polynomial :math:`P_{ab, a'b'}(p) = \pi_{ab, a'b'}(p)` (on the mass shell) is given by :math:`\eqref{eq_general_field_spin_sum}`. We note that :math:`\eqref{eq_spin_sum_v_as_polynomial}` is consistent with :math:`\eqref{eq_general_field_spin_sums_as_pi}` due to the parity symmetry :math:`\eqref{eq_general_field_spin_sum_polynomial_parity}`, which we only proved in a special case.

Back to the evaluation of the propagator :math:`\eqref{eq_propagator_as_spin_sums}`. By plugging :math:`\eqref{eq_spin_sum_u_as_polynomial}` and :math:`\eqref{eq_spin_sum_v_as_polynomial}` into :math:`\eqref{eq_propagator_as_spin_sums}`, we can further rewrite the propagator as follows

.. math::
	:nowrap:

	\begin{equation}
		-\ifrak \Delta_{\ell m}(x, y) = \theta(x_0 - y_0) P_{\ell m}(-\ifrak \p_x) \Delta_+(x-y) + \theta(y_0 - x_0) P_{\ell m}(-\ifrak \p_x) \Delta_+(y-x)
		\label{eq_propagator_as_delta_plus}
	\end{equation}

where :math:`\Delta_+(x)` is defined by :math:`\eqref{eq_defn_Delta_plus}`, and :math:`P_{\ell m}(p)` takes the form of :math:`\eqref{eq_general_field_spin_sum_as_polynomial}` and is linear in :math:`p_0 = \sqrt{\pbf^2 + M^2}`. In particular, we see that :math:`\Delta_{\ell m}(x, y)` is really a function of :math:`x - y`, and henceforth can be written as

.. math::
	:nowrap:

	\begin{equation*}
		\Delta_{\ell m}(x - y) = \Delta_{\ell m}(x, y)
	\end{equation*}

Now we must remember that although :math:`P_{\ell m}(p)` for scalar, vector, and Dirac fields, look like a polynomial defined generally on the momentum space, they're really only defined on the mass shell, as shown in :math:`\eqref{eq_general_field_spin_sum_as_polynomial}` for general fields. We'll first try the poor man's extension of :math:`P_{\ell m}(p)` to a genuine momentum space polynomial by the following definition

.. math::
	:nowrap:

	\begin{equation*}
		P_{\ell m}^{(L)}(p) \coloneqq P_{\ell m}^{(0)}(\pbf) + p_0 P_{\ell m}^{(1)}(\pbf)
	\end{equation*}

where :math:`P_{\ell m}^{(0)}, P_{\ell m}^{(1)}` correspond to :math:`P_{ab, a'b'}, 2Q_{ab, a'b'}` in :math:`\eqref{eq_general_field_spin_sum_as_polynomial}`, respectively. Clearly :math:`P^{(L)}(p) = P(p)` when :math:`p` is on the mass shell, or equivalently :math:`p_0 = \sqrt{\pbf^2 + M^2}`. Here the supscript :math:`L` stands for linear, since :math:`P^{(L)}(p)` is linear in :math:`p_0`. This linearity, as we'll now demonstrate, turns out to be the key feature of this rather naive extension.

Since the step function :math:`\theta(\tau)`, which equals :math:`1` for :math:`\tau > 0` and vanishes for :math:`\tau < 0`, has the following property

.. math::
	:nowrap:

	\begin{equation*}
		\p_{x_0} \theta(x_0 - y_0) = -\p_{x_0} \theta(y_0 - x_0) = \delta(x_0 - y_0)
	\end{equation*}

we can rewrite :math:`\eqref{eq_propagator_as_delta_plus}` as follows

.. math::
	:nowrap:

	\begin{align*}
		-\ifrak \Delta_{\ell m}(x - y) &= P^{(L)}_{\ell m}(-\ifrak \p_x) \big( \theta(x_0-y_0) \Delta_+(x-y) + \theta(y_0-x_0) \Delta_+(y-x) \big) \\
			&\phantom{=} - \big( \ifrak \p_{x_0} \theta(x_0-y_0) \big) P^{(1)}_{\ell m}(-\ifrak \p_x) \Delta_+(x-y) - \big( \ifrak \p_{x_0} \theta(y_0-x_0) \big) P^{(1)}_{\ell m}(-\ifrak \p_x) \Delta_+(y-x) \\
			&= P^{(L)}_{\ell m}(-\ifrak \p_x) \big( \theta(x_0-y_0) \Delta_+(x-y) + \theta(y_0-x_0) \Delta_+(y-x) \big) \\
			&\phantom{=} - \blue{\ifrak  \delta(x_0-y_0) P^{(1)}_{\ell m}(-\ifrak \nabla) \left( \Delta_+(x-y) - \Delta_+(y-x) \right)} \\
			&= P^{(L)}_{\ell m}(-\ifrak \p_x) \big( \theta(x_0-y_0) \Delta_+(x-y) + \theta(y_0-x_0) \Delta_+(y-x) \big)
	\end{align*}

where the blue terms vanish since :math:`\Delta_+(x)` is even if :math:`x_0 = 0`.

Define the *Feynman propagator* :math:`\Delta_F(x)` by

.. math::
	:nowrap:

	\begin{equation}
		-\ifrak \Delta_F(x) \coloneqq \theta(x_0) \Delta_+(x) + \theta(-x_0) \Delta_+(-x)
		\label{eq_defn_feynman_propagator}
	\end{equation}

which, by the way, is also the propagator for scalar fields. Then the general propagator can be derived from it as follows

.. math::
	:nowrap:

	\begin{equation}
		\Delta_{\ell m}(x-y) = P^{(L)}_{\ell m}(-\ifrak \p_x) \Delta_F(x-y)
		\label{eq_general_propagator_from_feynman_propagator}
	\end{equation}

Now we'll evaluate :math:`\eqref{eq_defn_feynman_propagator}` as an integral over the momentum space. The key trick is to express the step function :math:`\theta(t)` as follows

.. math::
	:nowrap:

	\begin{equation}
		\theta(t) = \frac{-1}{2\pi\ifrak} \int_{-\infty}^{\infty} ds~\frac{e^{-\ifrak st}}{s + \ifrak \epsilon}
		\label{eq_theta_step_function_as_contour_integral}
	\end{equation}

where :math:`\epsilon > 0` is infinitesimally small. The fact that :math:`\theta(t)` vanishes for :math:`t < 0` and equals :math:`1` for :math:`t > 0` can be easily verified using the `residue theorem <https://en.wikipedia.org/wiki/Residue_theorem>`__. Plugging :math:`\eqref{eq_theta_step_function_as_contour_integral}` and :math:`\eqref{eq_defn_Delta_plus}` into :math:`\eqref{eq_defn_feynman_propagator}`, we can calculate as follows

.. math::
	:nowrap:

	\begin{align*}
		-\ifrak \Delta_F(x) &= \frac{-1}{\ifrak (2\pi)^4} \int d^3 q ds~\frac
				{\exp(-\ifrak (s + \sqrt{\qbf^2 + M^2}) x_0 + \ifrak \qbf \cdot \xbf) + \exp(\ifrak (s + \sqrt{\qbf^2 + M^2}) x_0 - \ifrak \qbf \cdot \xbf)}
				{2(s + \ifrak \epsilon) \sqrt{\qbf^2 + M^2}} \\
			&\xlongequal{\substack{p_0 = s + \sqrt{\qbf^2 + M^2} \\ \pbf = \qbf}} \frac{-1}{\ifrak (2\pi)^4} \int d^4 p~\frac
				{\exp(\ifrak p \cdot x) + \exp(-\ifrak p \cdot x)}
				{2\left( p_0 - \sqrt{\pbf^2 + M^2} + \ifrak \epsilon \right) \sqrt{\pbf^2 + M^2}} \\
			&= \frac{-1}{\ifrak (2\pi)^4} \int d^4 p~\frac{\exp(\ifrak p \cdot x)}{2\sqrt{\pbf^2 + M^2}} \left( \frac{1}
				{p_0 - \sqrt{\pbf^2 + M^2} + \ifrak \epsilon} - \frac{1}{p_0 + \sqrt{\pbf^2 + M^2} - \ifrak \epsilon} \right) \\
			&= \frac{1}{\ifrak (2\pi)^4} \int d^4 p~\frac{\exp(\ifrak p \cdot x)}{p^2 + M^2 - \ifrak \epsilon}
	\end{align*}

where in the last equality we've replaced :math:`2\sqrt{\pbf^2 + M^2} \epsilon` with :math:`\epsilon`, and ignored the :math:`\epsilon^2` term. It follows that

.. math::
	:nowrap:

	\begin{equation}
		\Delta_F(x) = (2\pi)^{-4} \int d^4 p~\frac{\exp(\ifrak p \cdot x)}{p^2 + M^2 - \ifrak \epsilon}
		\label{eq_feynman_propagator_as_momentum_space_integral}
	\end{equation}

.. dropdown:: Green's function for the Klein-Gordon equation
	:animate: fade-in-slide-down

	The Feynman propagator :math:`\Delta_F(x)` defined by :math:`\eqref{eq_feynman_propagator_as_momentum_space_integral}` turns out to be also a Green's function for the Klein-Gordon equation is the sense that

	.. math::
		:nowrap:

		\begin{equation}
			\left( \square - M^2 \right) \Delta_F(x) = -\delta^4(x)
			\label{eq_klein_gordon_greens_function}
		\end{equation}

	Indeed, it follows readily from the following formal expansion of the (one-dimensional) Dirac delta function

	.. math::
		:nowrap:

		\begin{equation*}
			\delta(x) = \frac{1}{2\pi} \int dy~e^{\pm \ifrak xy}
		\end{equation*}

	Although the factor :math:`-\ifrak \epsilon` in the denominator in :math:`\eqref{eq_feynman_propagator_as_momentum_space_integral}` plays no role in the defining equation :math:`\eqref{eq_klein_gordon_greens_function}` of the Green's function, it does determine the boundary condition of the Green's function in the follow sense

	.. math::
		:nowrap:

		\begin{align*}
			\Delta_F(x) &\to \frac{\ifrak}{(2\pi)^3} \int \frac{d^3 p}{2\sqrt{\pbf^2 + M^2}}~\exp\left( -\ifrak x_0 \sqrt{\pbf^2 + M^2} + \ifrak \pbf \cdot \xbf \right) \quad \text{as}~~x_0 \to +\infty \\
			\Delta_F(x) &\to \frac{\ifrak}{(2\pi)^3} \int \frac{d^3 p}{2\sqrt{\pbf^2 + M^2}}~\exp\left(+\ifrak x_0 \sqrt{\pbf^2 + M^2} - \ifrak \pbf \cdot \xbf \right) \quad \text{as}~~x_0 \to -\infty
		\end{align*}

	This can be seen either directly from :math:`\eqref{eq_defn_feynman_propagator}` or by factoring the denominator :math:`p^2 + M^2 = \left( \sqrt{\pbf^2+M^2} - p_0 \right)\left( \sqrt{\pbf^2+M^2} + p_0 \right)` and using the residue theorem.

Plugging :math:`\eqref{eq_feynman_propagator_as_momentum_space_integral}` into :math:`\eqref{eq_general_propagator_from_feynman_propagator}`, we can express the propagator as an integral over the entire momentum space (without the mass-shell constraint) as follows

.. math::
	:nowrap:

	\begin{equation}
		\Delta_{\ell m}(x, y) = \Delta_{\ell m}(x-y) = (2\pi)^{-4} \int d^4 p~\frac{P^{(L)}_{\ell m}(p) e^{\ifrak p \cdot (x-y)}}{p^2 + M^2 - \ifrak \epsilon}
		\label{eq_propagator_as_momentum_space_integral_linear}
	\end{equation}

This expression is, however, not Lorentz covariant in general since :math:`P^{(L)}_{\ell m}(p)`, being linear in :math:`p_0`, is not. More precisely, we'd like :math:`P_{\ell m}(p)` to satisfy the following

.. math::
	:nowrap:

	\begin{equation*}
		P_{\ell m}(\Lambda p) = D_{\ell \ell'}(\Lambda) D^{\ast}_{m m'}(\Lambda) P_{\ell' m'}(p)
	\end{equation*}

where we recall the :math:`m` index correspond to a field adjoint according to :math:`\eqref{eq_feynman_rule_propagator}`. Nonetheless, we see from :math:`\eqref{eq_p_polynomial_scalar}` -- :math:`\eqref{eq_p_polynomial_dirac}` that the polynomial :math:`P` is Lorentz covariant in the case of scalar, vector, and Dirac fields. Among these three cases, the only case where :math:`P` is not linear in :math:`p_0` is vector field. Indeed, we have the following from :math:`\eqref{eq_p_polynomial_vector}`

.. math::
	:nowrap:

	\begin{align*}
		P^{(L)}_{\mu \nu}(p) &= \eta_{\mu \nu} + M^{-2} \left( p_{\mu}p_{\nu} + \delta_{\mu}^0 \delta_{\nu}^0 \left( -p_0^2 + \pbf^2 + M^2 \right) \right) \\
			&= P_{\mu \nu}(p) + M^{-2} (p^2 + M^2) \delta_{\mu}^0 \delta_{\nu}^0
	\end{align*}

Plugging into :math:`\eqref{eq_propagator_as_momentum_space_integral_linear}`, we have

.. math::
	:nowrap:

	\begin{equation}
		\Delta^{\text{vector}}_{\mu \nu}(x, y) = (2\pi)^{-4} \int d^4 p~\frac{P_{\mu \nu}(p) e^{\ifrak p \cdot (x-y)}}{p^2 + M^2 - \ifrak \epsilon} + \blue{M^{-2}\delta^4(x-y) \delta_{\mu}^0 \delta_{\nu}^0}
		\label{eq_vector_field_propagator_needs_local_term}
	\end{equation}

We see that the price to pay for making :math:`P_{\mu \nu}(p)` Lorentz covariant is the blue term, which is local in the sense that it's non-vanishing only when :math:`x = y`.

.. todo::
	It's claimed in [Wei95]_ page 278 -- 279 that such an extra non-covariant term may be countered by a modification to the interaction density, but I have had a hard time seeing why that's the case. In order to not get stuck at this point, we'll ignore the difference between :math:`P^{(L)}` and :math:`P` at momentum so that :math:`\eqref{eq_propagator_as_momentum_space_integral_linear}` and write

	.. math::
		:nowrap:

		\begin{equation}
			\Delta_{\ell m}(x, y) = (2\pi)^{-4} \int d^4 p~\frac{P_{\ell m}(p) e^{\ifrak p \cdot (x-y)}}{p^2 + M^2 - \ifrak \epsilon}
			\label{eq_propagator_as_momentum_space_integral}
		\end{equation}

	and remember to be extra careful when, in a concrete case, :math:`P` is not linear in :math:`p_0`.

.. _sec_feynman_rules_in_momentum_space:

Feynman rules in momentum space
+++++++++++++++++++++++++++++++

To turn the Feynman rules derived in :ref:`sec_spacetime_feynman_rules`, specifically :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_propagator}`, from an integral over spacetime coordinates as in :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}` to an integral over the momentum space coordinates, we need to integrate out the :math:`x` coordinates. Indeed, as can be seen from :math:`\eqref{eq_feynman_rule_a_psi_dagger}` -- :math:`\eqref{eq_feynman_rule_psi_dagger_a_dagger}`, as well as :math:`\eqref{eq_propagator_as_momentum_space_integral}`, the :math:`x` variables appear only in the exponential terms. More precisely, at each vertex, the only term that involves :math:`x` is the following

.. math::
	:nowrap:

	\begin{equation*}
		\exp\left( \sum \ifrak x \cdot p_{\text{in}} - \sum \ifrak x \cdot p_{\text{out}} + \sum \ifrak x \cdot p_{\text{entering}} - \sum \ifrak x \cdot p_{\text{leaving}} \right)
	\end{equation*}

where :math:`p_{\text{in}}` and :math:`p_{\text{out}}` denote the momenta of the in- and out-state (anti-)particles that connect the vertex, respectively, if there is any, and :math:`p_{\text{entering}}` and :math:`p_{\text{leaving}}` denote the momenta on the internal edges that entering and leaving the vertex, respectively. Integrating :math:`x` out, we get the following momentum-conservation factor

.. math::
	:nowrap:

	\begin{equation*}
		(2\pi)^4 \delta^4 \left( \sum p_{\text{in}} - \sum p_{\text{out}} + \sum p_{\text{entering}} - \sum p_{\text{leaving}} \right)
	\end{equation*}

assigned to each vertex.

Now with all the :math:`x` coordinates integrated out, we can reformulate the :ref:`diagrammatic Feynman rules <fig_spacetime_feynman_diagrams>` as follows

.. _fig_momentum_space_feynman_diagrams:

.. figure:: ./static/quantum-theory-of-fields/momenta-propagators.svg
	:align: center

	Figure. Feynman rules in momentum space.

To evaluate the (connected) S-matrix of the form :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}`, we follow the steps explained in :ref:`sec_spacetime_feynman_rules` with one modification:

* Associate to each vertex a factor :math:`-\ifrak (2\pi)^4 \delta^4 \left( \sum p_{\text{in}} - \sum p_{\text{out}} + \sum p_{\text{entering}} - \sum p_{\text{leaving}} \right)`, and to each edge a factor as indicated in the figure above.

Indeed, the momentum-conservation delta functions at the vertices correspond precisely to the delta functions appeared in :ref:`sec_cluster_decomposable_hamiltonians`. Hence the same argument implies, once again, that the connected S-matrix contains exactly one momentum-conservation delta function.

Before moving onto the discussion about the external edges, let's revisit the example calculations of S-matrices in :math:`\psi^{\dagger} \psi \phi`-interaction in momentum space as follows.

:ref:`Fermion-boson scattering <listing_fermion_boson_scattering>` in momentum space
	The goal is to re-evaluate the (connected) S-matrix :math:`S^C_{1'2',12}`, where :math:`1` stands for the fermion, and :math:`2` stands for the boson, in momentum space. Using the momentum space Feynman rules, together with obvious abbreviations, the calculation goes as follows.

	.. math::
		:nowrap:

		\begin{align*}
			S^C_{1'2',12} &=
					\sum_{\ell m k, \ell' m' k'} (-\ifrak)^{2+1} (2\pi)^{8-6-4} g_{\ell m k} g_{\ell' m' k'} \int d^4 p~\frac{P_{m \ell'}(p)}{p^2 + M^2 - \ifrak \epsilon} u_{m'}(1) u^{\ast}_{\ell}(1') \\
				&\phantom{=~} \times \left( u_{k'}(2) u^{\ast}_k(2') \delta^4(p-p_1-p_2) \delta^4(p'_1+p'_2-p) + u_k(2) u^{\ast}_{k'}(2') \delta^4(p+p'_2-p_1) \delta^4(p'_1-p_2-p) \right) \\
				&= \ifrak (2\pi)^{-2} \delta^4(p'_1+p'_2-p_1-p_2) \sum_{\ell m k, \ell' m' k'} g_{\ell m k} g_{\ell' m' k'} u_{m'}(1) u^{\ast}_{\ell}(1') \\
				&\phantom{=~} \times \left( \frac{P_{m \ell'}(p_1+p_2)}{(p_1+p_2)^2+M^2-\ifrak \epsilon} u_{k'}(2) u^{\ast}_k(2') + \frac{P_{m \ell'}(p_1-p'_2)}{(p_1-p'_2)^2+M^2-\ifrak \epsilon} u_k(2) u^{\ast}_{k'}(2') \right) \\
				&= \ifrak (2\pi)^{-2} \delta^4(p'_1+p'_2-p_1-p_2) \sum_{k, k'} \bigg( \phantom{)} \\
				&\phantom{=~} \left( \blue{u^{\dagger}(1') \Gamma_k \frac{P(p_1+p_2)}{(p_1+p_2)^2+M^2-\ifrak \epsilon} \Gamma_{k'} u(1)} \right) u_{k'}(2) u^{\ast}_k(2') \\
				&\phantom{=~(} + \left( \blue{u^{\dagger}(1') \Gamma_k \frac{P(p_1-p'_2)}{(p_1-p'_2)^2+M^2-\ifrak \epsilon} \Gamma_{k'} u(1)} \right) u_k(2) u^{\ast}_{k'}(2') \bigg)
		\end{align*}

	where the blue terms are supposed to be understood as matrix multiplications with

	.. math::
		:nowrap:

		\begin{equation}
			(\Gamma_k)_{\ell m} \coloneqq g_{\ell m k}
			\label{eq_three_way_coupling_constant_as_matrix}
		\end{equation}

	and :math:`u^{\dagger}` as a transpose-conjugated row vector, and :math:`M` is the mass of the intermediate fermion.

	One can see from the calculation above that a general S-matrix can readily be read off from the Feynman diagrams, and is made up of the following pieces

	1. An appropriate power of :math:`\ifrak` and :math:`\pi` determined by the number of edges and vertices.

	2. A momentum-conservation delta function equating the total in- and out-state particles.

	3. A summation of (matrix) products of field coefficients, propagator integrands, and coupling constants, one for each Feynman diagram. In particular, the momentum carried by a propagator is determined by the momentum-conservation law at each vertex.

:ref:`Fermion-fermion scattering <listing_fermion_fermion_scattering>` in momentum space
	Using the observation we made from the previous calculation, we can immediately write down the S-matrix :math:`S^C_{1'2',12}`, where both :math:`1` and :math:`2` are fermions, as follows

	.. math::
		:nowrap:

		\begin{align*}
			S^C_{1'2',12} &= \ifrak (2\pi)^{-2} \delta^4(p'_1+p'_2-p_1-p_2) \sum_{k, k'} \bigg( \phantom{)} \\
				&\phantom{=~} \frac{P_{k k'}(p'_1 - p_1)}{(p'_1-p_1)^2 + M^2 - \ifrak\epsilon} \left( u^{\dagger}(1) \Gamma_k u(1') \right) \left( u^{\dagger}(2) \Gamma_{k'} u(2') \right) \\
				&\phantom{=(} - \frac{P_{k k'}(p'_2-p_1)}{(p'_2-p_1)^2 + M^2 - \ifrak\epsilon} \left( u^{\dagger}(1) \Gamma_k u(2') \right) \left( u^{\dagger}(2) \Gamma_{k'} u(1') \right) \bigg)
		\end{align*}

:ref:`Boson-boson scattering <listing_boson_boson_scattering>` in momentum space
	We didn't actually calculate the (4th order) boson-boson scattering S-matrix :math:`S^C_{1'2',12}` in spacetime coordinates due to its complexity. But it becomes much simpler in momentum space coordinates, and can be calculated as follows. First, let's figure out the powers of :math:`\ifrak` and :math:`\pi`, respectively. Since there are :math:`4` vertices and :math:`4` internal edges, each of which contribute one :math:`-\ifrak`, we get a contribution of :math:`(-\ifrak)^8 = 1`. Moreover, since there are equal numbers of vertices and internal edges, each of which contribute :math:`(2\pi)^4` and :math:`(2\pi)^{-4}`, respectively, and :math:`4` external edges, each of which contribute :math:`(2\pi)^{-3/2}`, we get a total contribution of :math:`(2\pi)^{-6}`. Remembering an additional minus sign coming from the fermionic loop, we have the following

	.. math::
		:nowrap:

		\begin{align*}
			S^C_{1'2',12} &= -(2\pi)^{-6} \delta^4(p'_1+p'_2-p_1-p_2) \sum_{k_1, k_2, k'_1, k'_2} u^{\ast}_{k'_1}(1') u^{\ast}_{k'_2}(2') u_{k_1}(1) u_{k_2}(2) \\
				&\phantom{=~} \times \int d^4 p~\op{Tr} \bigg( \frac{P(p)}{p^2+M^2-\ifrak\epsilon} \Gamma_{k'_1} \frac{P(p-p'_1)}{(p-p'_1)^2+M^2-\ifrak\epsilon}
					\Gamma_{k'_2} \frac{P(p-p'_1-p'_2)}{(p-p'_1-p'_2)^2+M^2-\ifrak\epsilon} \Gamma_{k_2} \\
				&\phantom{=~} \times \frac{P(p-p'_1-p'_2+p_2)}{(p-p'_1-p'_2+p_2)^2+M^2-\ifrak\epsilon} \Gamma_{k_1} \bigg) + \cdots
		\end{align*}

	where :math:`\cdots` denotes valuations of other (4th order) Feynman diagrams.


External edges off the mass shell
+++++++++++++++++++++++++++++++++

By transforming propagators into integrals over the :math:`4`-momentum space as in :math:`\eqref{eq_propagator_as_momentum_space_integral}`, and integrating out the spacetime coordinates, we've been able to express the S-matrix as an integral over a product of :math:`4`-momentum space coordinates, one for each internal edge that doesn't split the Feynman diagram into disconnected components. In particular, as we've seen in several examples from the previous section, the S-matrix involves no momentum space integrals at all when the Feynman diagram is a tree.

However, the external edges are still confined to the mass shell since the in- and out-state particles are. We'd like to relax such confinement to remove any mass-shell constraints. It will allow us to calculate contributions from a large Feynman diagram in terms of its smaller local pieces, and facilitate the derivation of the path integral formalism later. The idea of a "local" Feynman diagram, i.e., a diagram without external edges, is to add additional "vertices at infinity" at the places of in- and out-state particles, which turn the external edges into internal ones. Moreover, unlike the (internal) vertices, whose spacetime coordinates will be integrated as in :math:`\eqref{eq_s_matrix_fully_expanded_by_timed_ordered_interaction_density}`, the spacetime coordinates of the added vertices will be kept intact. In particular, no momentum conservation law will be imposed on these vertices.

It turns out that this procedure can be generalized to a field theory with external interaction, in which case the interaction-picture perturbation term takes the following form

.. math::
	:nowrap:

	\begin{equation}
		V_{\epsilon}(t) = V(t) + \sum_a \int d^3 x~\epsilon_a(t, \xbf) o_a(t, \xbf)
		\label{eq_defn_interaction_perturbation_term_with_external_fields}
	\end{equation}

Here :math:`\epsilon_a(t, \xbf)` are the infinitesimal parameters, and :math:`o_a(t, \xbf)` are the external currents in the interaction picture in the following sense

.. math::
	:nowrap:

	\begin{equation*}
		o_a(t) = \exp(\ifrak H_0 t) o_a(0) \exp(-\ifrak H_0 t)
	\end{equation*}

in analogy to :math:`\eqref{eq_defn_interaction_perturbation_term}`.

Now under the perturbation :math:`V_{\epsilon}(t)`, the S-matrix, whose entry :math:`S_{\beta\alpha}` is a complex number, becomes a functional :math:`S_{\beta\alpha}[\epsilon]` in complex functions :math:`\epsilon_a`. The functional derivative

.. math::
	:nowrap:

	\begin{equation}
		\left. \frac{\delta^r S_{\beta\alpha}[\epsilon]}{\delta\epsilon_a(x) \delta\epsilon_b(y) \cdots} \right|_{\epsilon=0}
		\label{eq_functional_derivative_of_s_matrix}
	\end{equation}

can be evaluated using the usual Feynman rules. Indeed, besides the internal vertices coming from the interaction density :math:`\Hscr(t, \xbf)`, we must also include external vertices coming from :math:`o(t, \xbf)`. For example, if :math:`o(t, \xbf)` are monomials of fields and field adjoints, then there will be :math:`r` external vertices, each of which has as many incoming and outgoing edges as there are fields and field adjoints in the corresponding :math:`o(t, \xbf)`, respectively. In particular, if :math:`o(t, \xbf)` are all degree one, then we recover the off-mass-shell Feynman diagrams, as promised.

.. dropdown:: Gell-Mann and Low's Theorem
	:icon: unlock
	:animate: fade-in-slide-down

	It's turns out that :math:`\eqref{eq_functional_derivative_of_s_matrix}` may be evaluated to an amplitude as follows

	.. math::
		:nowrap:

		\begin{equation}
			\left. \frac{\delta^r S_{\beta\alpha}[\epsilon]}{\delta\epsilon_{a_1}(x_1) \cdots \delta\epsilon_{a_r}(x_r)} \right|_{\epsilon=0}
				= (-\ifrak)^r \left( \Psi_{\beta}^+, T\left\{ O_{a_1}(x_1) \cdots O_{a_r}(x_r) \right\} \Psi_{\alpha}^- \right)
			\label{eq_gell_mann_low_theorem}
		\end{equation}

	where

	.. math::
		:nowrap:

		\begin{equation}
			O_a(t, \xbf) = \exp(\ifrak H t) o_a(0, \xbf) \exp(-\ifrak H t) = \Omega(t) o_a(t, \xbf) \Omega^{-1}(t)
			\label{eq_heisenberg_picture_external_field}
		\end{equation}

	using the definition :math:`\eqref{eq_defn_of_Omega}` represent the Heisenberg picture operators, and :math:`\Psi_{\alpha}^-, \Psi_{\beta}^+` are the in- and out-states, respectively.

	The formula :math:`\eqref{eq_gell_mann_low_theorem}` is know as `Gell-Mann and Low's theorem <https://en.wikipedia.org/wiki/Gell-Mann_and_Low_theorem>`_. To see this, let's write out the functional derivative, following :math:`\eqref{eq_s_matrix_power_series_expansion_time_ordered}`, as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\left. \frac{\delta^r S_{\beta\alpha[\epsilon]}}{\delta\epsilon_{a_1}(x_1) \cdots \delta\epsilon_{a_r}(x_r)} \right|_{\epsilon=0}
				= \sum_{N=0}^{\infty} \frac{(-\ifrak)^{N+r}}{N!} \int_{-\infty}^{\infty} d\tau_1 \cdots d\tau_N
					\left( \Phi_{\beta}, T\left\{ V(\tau_1) \cdots V(\tau_N) o_{a_1}(x_1) \cdots o_{a_r}(x_r) \right\} \Phi_{\alpha} \right)
		\end{equation*}

	where we note that each external vertex corresponding to some :math:`o_a(x)` comes with a factor :math:`-\ifrak`, but the spacetime coordinates :math:`x` are not integrated. Assume without loss of generality that :math:`x_1, \cdots, x_r` are already time-ordered in the sense that

	.. math::
		:nowrap:

		\begin{equation}
			(x_1)_0 \geq (x_2)_0 \geq \cdots \geq (x_r)_0
			\label{eq_gell_mann_low_external_fields_are_time_ordered}
		\end{equation}

	Then we can expand the time-ordered product above as follows

	.. math::
		:nowrap:

		\begin{align*}
			\left. \frac{\delta^r S_{\beta\alpha[\epsilon]}}{\delta\epsilon_{a_1}(x_1) \cdots \delta\epsilon_{a_r}(x_r)} \right|_{\epsilon=0}
				&= \sum_{N=0}^{\infty} \frac{(-\ifrak)^{N+r}}{N!} \sum_{N_0 + \cdots N_r = N} \frac{N!}{N_0! \cdots N_r!} \\
				&\phantom{=~} \times \int_{(x_1)_0}^{\infty} d\tau_{01} \cdots d\tau_{0N_0} \int_{(x_2)_0}^{(x_1)_0} d\tau_{11} \cdots d\tau_{1N1} \cdots \int_{-\infty}^{(x_r)_0} d\tau_{r1} \cdots d\tau_{rN_r} \\
				&\phantom{=~} \times \Big( \Phi_{\beta}, T\left\{ V(\tau_{01}) \cdots V(\tau_{0N_0}) \right\} o_{a_1}(x_1) T\left\{ V(\tau_{11}) \cdots V(\tau_{1N1}) \right\} o_{a_2}(x_2) \cdots \\
				&\phantom{=~\times} \times \cdots o_{a_r}(x_r) T\left\{ V(\tau_{r1}) \cdots V(\tau_{rN_r}) \right\} \Phi_{\alpha} \Big)
		\end{align*}

	where :math:`\tau_{01}, \cdots, \tau_{0N0}` are greater than :math:`(x_1)_0`, and :math:`\tau(11), \cdots, \tau(1N1)` are between :math:`(x_2)_0` and :math:`(x_1)_0`, and so on. In other words, we're summing over a partition of :math:`\tau_1, \cdots, \tau_N` into :math:`r+1` unordered groups, separated by :math:`(x_1)_0, \cdots, (x_r)_0`. The combinatorial coefficient :math:`N! / (N_0! \cdots N_r!)` counts exactly the number of such partitions.

	It follows that

	.. math::
		:nowrap:

		\begin{equation*}
			\left. \frac{\delta^r S_{\beta\alpha[\epsilon]}}{\delta\epsilon_{a_1}(x_1) \cdots \delta\epsilon_{a_r}(x_r)} \right|_{\epsilon=0}
				= (-\ifrak)^r \left( \Phi_{\beta}, U(\infty, (x_1)_0) o_{a_1}(x_1) U((x_1)_0, (x_2)_0) o_{a_2}(x_2) \cdots o_{a_r}(x_r) U((x_r)_0, -\infty) \Phi_{\alpha} \right)
		\end{equation*}

	where

	.. math::
		:nowrap:

		\begin{equation*}
			U(\tau', \tau) = \sum_{N=0}^{\infty} \frac{(-\ifrak)^N}{N!} \int_{\tau}^{\tau'} d\tau_1 \cdots d\tau_N T\left\{ V(\tau_1) \cdots V(\tau_N) \right\}
				\xlongequal{\eqref{eq_s_matrix_power_series_expansion_time_ordered},~\eqref{eq_s_operator_by_u}}
				\Omega^{-1}(\tau') \Omega(\tau)
		\end{equation*}

	Combining with :math:`\eqref{eq_heisenberg_picture_external_field}` and :math:`\eqref{eq_defn_of_Omega}`, we conclude that

	.. math::
		:nowrap:

		\begin{align*}
			\left. \frac{\delta^r S_{\beta\alpha[\epsilon]}}{\delta\epsilon_{a_1}(x_1) \cdots \delta\epsilon_{a_r}(x_r)} \right|_{\epsilon=0}
				&= (-\ifrak)^r \left( \Omega(\infty)\Phi_{\beta}, O_{a_1}(x_1) \cdots O_{a_r}(x_r) \Omega(-\infty)\Phi_{\alpha} \right) \\
				&= (-\ifrak)^r \left( \Psi_{\beta}^+, O_{a_1}(x_1) \cdots O_{a_r}(x_r) \Psi_{\alpha}^- \right)
		\end{align*}

	under the assumption :math:`\eqref{eq_gell_mann_low_external_fields_are_time_ordered}`. Finally the Gell-Mann and Low theorem :math:`\eqref{eq_gell_mann_low_theorem}` follows from the symmetry property that both sides are symmetric under swaps of bosonic (external) fields and anti-symmetric under swaps of fermionic (external) fields.


The Canonical Formalism
-----------------------

The quantum theory we've been developing so far has been based almost solely on the symmetry principles, especially Lorentz symmetries. This is a very satisfying approach since it's logically clean and relies only on the most fundamental principles, however, this is not the way quantum theory historically had been developed. Not surprisingly, the original development of quantum theory is much messier and requires substantial experience in "classical" physics. It's largely based on the so-called *Lagrangian formalism*, which is a readily well-established principle in classical physics and can be "quantized". The main goal of this chapter is to go through this formalism, not for historical sake, but because it offers a particularly convenient way to construct Hamiltonians that generate Lorentz-invariant S-matrices, which has been difficult for us as can be seen in :ref:`sec_feynman_rules_in_momentum_space`.

Canonical variables
^^^^^^^^^^^^^^^^^^^

We've seen in :ref:`sec_quantum_fields_and_antiparticles` a few ways of constructing (Lorentz-invariant) interaction densities. However, we don't have a systematic way to do so. The so-called Lagrangian formalism will not provide a systematic solution either, but it'll allow us to construct more interesting interaction densities (from classical physics theories), to the extent that all known quantum field theories arise in this way! In addition, it'll shed light on the mysterious local terms as for example in :math:`\eqref{eq_vector_field_propagator_needs_local_term}`, that are needed to compensate for a Lorentz-invariant momentum space propagator.

The offer from the Lagrangian formalism regarding constructing a quantum field theory is the following. Instead of using the fields :math:`\eqref{eq_defn_annihilation_field}` -- :math:`\eqref{eq_defn_creation_field}` to construct the Hamiltonians, we'll use the so-called *canonical variables*, which have particularly simple (equal time) commutation relations. More precisely, it consists of a collection of quantum operators :math:`q_n(t, \xbf)` and its canonical conjugates :math:`p_n(t, \xbf)`, which satisfy the following (anti-)commutation relations

.. math::
	:nowrap:

	\begin{align}
		\left[ q_n(t, \xbf), p_{n'}(t, \ybf) \right]_{\pm} &= \ifrak \delta^3(\xbf - \ybf) \delta_{n n'}
		\label{eq_canonical_commutation_relation_1} \\
		\left[ q_n(t, \xbf), q_{n'}(t, \ybf) \right]_{\pm} &= 0
		\label{eq_canonical_commutation_relation_2} \\
		\left[ p_n(t, \xbf), p_{n'}(t, \ybf) \right]_{\pm} &= 0
		\label{eq_canonical_commutation_relation_3}
	\end{align}

where :math:`\pm` correspond to when the particle under question is fermionic or bosonic, respectively.

To see how canonical variables may be constructed from fields considered in :ref:`sec_quantum_fields_and_antiparticles`, let's consider a few examples.

Scalar fields
	Let's start by considering scalar fields of particles that are their own antiparticles. Using notations from :ref:`sec_scalar_field`, it means that :math:`\psi(x) = \psi^{
	\dagger}(x)`, i.e., the field is Hermitian. It follows then from :math:`\eqref{eq_scalar_field_commutator}` and :math:`\eqref{eq_defn_Delta}` that

	.. math::
		:nowrap:

		\begin{equation}
			\left[ \psi(x), \psi(y) \right] = \Delta(x-y) = \frac{1}{(2\pi)^3} \int \frac{d^3 p}{2p_0} \left( e^{\ifrak p \cdot (x-y)} - e^{-\ifrak p \cdot (x-y)} \right)
			\label{eq_scalar_field_commutation_relation_reproduced}
		\end{equation}

	where :math:`p_0 = \sqrt{\pbf^2 + m^2}`.

	We claim that the canonical commutation relations :math:`\eqref{eq_canonical_commutation_relation_1}` -- :math:`\eqref{eq_canonical_commutation_relation_3}` are satisfied by

	.. math::
		:nowrap:

		\begin{align}
			q(t, \xbf) &\coloneqq \psi(t, \xbf)
			\label{eq_defn_q_scalar_field_self_dual} \\
			p(t, \xbf) &\coloneqq \dot{\psi}(t, \xbf)
			\label{eq_defn_p_scalar_field_self_dual}
		\end{align}

	Indeed, it follows from the following calculations

	.. math::
		:nowrap:

		\begin{alignat}{2}
			\left[ q(t, \xbf), p(t, \xbf) \right] &= \left[ \psi(t, \xbf), \dot{\psi}(t, \ybf) \right] &&= -\dot{\Delta}(0, \xbf-\ybf) = \ifrak \delta^3(\xbf-\ybf)
			\nonumber \\
			\left[ q(t, \xbf), q(t, \xbf) \right] &= \left[ \psi(t, \xbf), \psi(t, \ybf) \right] &&= \Delta(0, \xbf-\ybf) = 0
			\label{eq_canonical_commutator_scalar_field_self_dual_qq} \\
			\left[ p(t, \xbf), p(t, \xbf) \right] &= \left[ \dot{\psi}(t, \xbf), \dot{\psi}(t, \ybf) \right] &&= -\ddot{\Delta}(0, \xbf-\ybf) = 0
			\nonumber
		\end{alignat}

	Now for particles that are different from their antiparticles, we must modify :math:`\eqref{eq_defn_q_scalar_field_self_dual}` -- :math:`\eqref{eq_defn_p_scalar_field_self_dual}` as follows

	.. math::
		:nowrap:

		\begin{align*}
			q(t, \xbf) &= \psi(t, \xbf) \\
			p(t, \xbf) &= \dot{\psi}^{\dagger}(t, \xbf)
		\end{align*}

	and note that in this case :math:`\left[ \psi(t, \xbf), \psi(t', \ybf) \right] = 0`, in contrast to :math:`\eqref{eq_canonical_commutator_scalar_field_self_dual_qq}`.

Spin-:math:`1` vector fields
	Consider once again particles that are self-charge-dual. Using notations from :ref:`sec_spin_1_vector_fields`, we recall the commutation relation :math:`\eqref{eq_vector_field_commutator}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\left[ \psi_{\mu}(x), \psi_{\nu}(y) \right] = \left( \eta_{\mu\nu} - \frac{\p_{\mu} \p_{\nu}}{m^2} \right) \Delta(x-y)
		\end{equation*}

	The canonical variables in this case can be defined as follows

	.. math::
		:nowrap:

		\begin{align}
			q_i(t, \xbf) &= \psi_i(t, \xbf)
			\label{eq_defn_q_vector_field_self_dual} \\
			p_i(t, \xbf) &= \dot{\psi}_i(t, \xbf) - \frac{\p \psi_0(t, \xbf)}{\p x_i}
			\label{eq_defn_p_vector_field_self_dual}
		\end{align}

	where :math:`i=1,2,3`. Indeed, let's calculate the equal-time commutators as follows

	.. math::
		:nowrap:

		\begin{align*}
			\left[ q_i(t, \xbf), p_j(t, \ybf) \right] &= \left[ \psi_i(t, \xbf), \dot{\psi}_j(t, \ybf) \right] - \left[ \psi_i(t, \xbf), \frac{\p \psi_0(t, \ybf)}{\p y_j} \right] \\
				&= -\left( \eta_{ij} -\frac{\p_i \p_j}{m^2} \right) \dot{\Delta}(0, \xbf-\ybf) - \left. \frac{\p_i \p_0}{m^2} \right|_{t=0} \left( \p_j \Delta(t, \xbf-\ybf) \right) \\
				&= \ifrak \delta^3(\xbf-\ybf) \delta_{ij} \\
			\left[ q_i(t, \xbf), q_j(t, \ybf) \right] &= \left( \eta_{ij} - \frac{\p_i \p_j}{m^2} \right) \Delta(0, \xbf-\ybf) = 0 \\
			\left[ p_i(t, \xbf), p_j(t, \ybf) \right] &= \left[ \dot{\psi}_i(t, \xbf), \dot{\psi}_j(t, \ybf)\right] + \p_{x_i} \p_{y_j} \left[ \psi_0(t, \xbf), \psi_0(t, \ybf) \right] \\
			&\phantom{=} - \p_{x_i} \left[ \psi_0(t, \xbf), \dot{\psi}_j(t, \ybf) \right] - \p_{y_j} \left[ \dot{\psi}_i(t, \xbf), \psi_0(t, \ybf) \right] = 0
		\end{align*}

	We've omitted some details about the vanishing of the last quantities -- it turns out that the the first and second terms cancel out, and the third and the fourth terms also cancel out.

	In any case, we've constructed three pairs of canonical variables, one for each spatial index. But what about the time index? It turns out that :math:`\psi_0` is *not* an independent variable. Indeed, we can derive from :math:`\eqref{eq_defn_p_vector_field_self_dual}` an expression of :math:`\psi_0` as follows

	.. math::
		:nowrap:

		\begin{align*}
			& & p_i & = \p_0 \psi_i - \p_i \psi_0 \\
			& \xRightarrow{\phantom{\eqref{eq_klein_gordon}}} & \p_i p_i & = \p_0 \p_i \psi_i - \p^2_i \psi_0 \\
			& \xRightarrow{\phantom{\eqref{eq_klein_gordon}}} & \nabla \cdot \pbf & = \p_0 \sum_{i=1}^3 \p_i \psi_i - \sum_{i=1}^3 \p^2_i \psi_0 \\
			& \xRightarrow{\eqref{eq_vector_field_gauge_fixing_condition}} & \nabla \cdot \pbf & = \p_0^2 \psi_0 - \sum_{i=1}^3 \p_i^2 \psi_0 = -\square \psi_0 \\
			& \xRightarrow{\eqref{eq_klein_gordon}} & \psi_0 & = -m^{-2} \nabla \cdot \pbf
		\end{align*}

Spin-:math:`1/2` Dirac fields
	Recall the anti-commutator of Dirac fields :math:`\eqref{eq_dirac_field_commutator}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\left[ \psi_{\ell}(x), \psi^{\dagger}_{\ell'}(y) \right]_+ = \left( (-\gamma^{\mu} \p_{\mu} + m) \beta \right)_{\ell \ell'} \Delta(x-y)
		\end{equation*}

	where :math:`\ell, \ell'` are indexes corresponding to the two spin :math:`z`-component :math:`\pm 1/2`. Assuming that particle under question has distinct antiparticle, i.e., it's not a Majorana fermion, the following holds trivially

	.. math::
		:nowrap:

		\begin{equation*}
			\left[ \psi_{\ell}(x), \psi_{\ell'}(y) \right]_+ = 0
		\end{equation*}

	It follows that the canonical variables can be defined by

	.. math::
		:nowrap:

		\begin{align*}
			q_{\ell}(x) &= \psi_{\ell}(x) \\
			p_{\ell}(x) &= \ifrak \psi^{\dagger}_{\ell}(x)
		\end{align*}

	Indeed, the only nontrivial (and non-vanishing) anti-commutator can be calculated as follows

	.. math::
		:nowrap:

		\begin{align*}
			\left[ q_{\ell}(t, \xbf), p_{\ell'}(t, \ybf) \right]_+ &= \ifrak \left[ \psi_{\ell}(t, \xbf), \psi_{\ell'}^{\dagger}(t, \ybf) \right]_+ \\
				&= -\ifrak \left( \gamma^0 \beta \right)_{\ell \ell'} \dot{\Delta}(0, \xbf-\ybf) \\
				&= \ifrak \delta^3(\xbf-\ybf) \delta_{\ell \ell'}
		\end{align*}

Through these examples, we see that there is no particular pattern in how one may define canonical variables. In fact, one doesn't really define canonical variables in this way either -- they are simply given for granted in the Lagrangian formalism as we will see.

We begin by a general discussion on functionals :math:`F[q(t), p(t)]` of canonical variables, since both Hamiltonians and Lagrangians will be such functionals. A few notes are in order. First we've used a shorthand notation :math:`q(t)` and :math:`p(t)` to denote a collection of canonical variables. Moreover, in writing :math:`q(t)` (and similarly for :math:`p(t)`) we implicitly think of them as fields at a given time. Indeed, as we'll see, the time variable plays an exceptional role in the Lagrangian formalism, in contrast to our mindset so far that space and time are all mixed up in a Lorentz invariant theory. Finally, we've used square bracket to differentiate it from regular functions of spacetime or momentum variables.

At the heart of the Lagrangian formalism lies a variational principle. Hence it's crucial to be able to take infinitesimal variations on :math:`F[q(t), p(t)]`, which we write as follows

.. math::
	:nowrap:

	\begin{equation}
		\delta F[q(t), p(t)] = \int d^3 x \sum_n \left( \delta q_n(t, \xbf) \frac{\delta F[q(t), p(t)]}{\delta q_n(t, \xbf)} + \frac{\delta F[q(t), p(t)]}{\delta p_n(t, \xbf)} \delta p_n(t, \xbf) \right)
		\label{eq_infinitesimal_variation_of_functional_of_canonical_variables}
	\end{equation}

Here the infinitesimal fields :math:`\delta q_n` and :math:`\delta p_n` are assumed to (anti-)commute with all other fields. Now assuming :math:`F[q(t), p(t)]` is written so that all the :math:`q` fields lie to the left of all the :math:`p` fields, then :math:`\eqref{eq_infinitesimal_variation_of_functional_of_canonical_variables}` can be realized by the following definition of variational derivatives

.. math::
	:nowrap:

	\begin{align*}
		\frac{\delta F[q(t), p(t)]}{\delta q_n(t, \xbf)} \coloneqq \ifrak \big[ p_n(t, \xbf), F[q(t), p(t)] \big] \\
		\frac{\delta F[q(t), p(t)]}{\delta p_n(t, \xbf)} \coloneqq \ifrak \big[ F[q(t), p(t)], q_n(t, \xbf) \big]
	\end{align*}

Hamiltonian and Lagrangian for free fields
++++++++++++++++++++++++++++++++++++++++++

For free fields we have

.. math::
	:nowrap:

	\begin{align}
		q_n(t, \xbf) &= e^{\ifrak H_0 t} q_n(0, \xbf) e^{-\ifrak H_0 t}
		\label{eq_free_field_q_time_evolution} \\
		p_n(t, \xbf) &= e^{\ifrak H_0 t} p_n(0, \xbf) e^{-\ifrak H_0 t}
		\label{eq_free_field_p_time_evolution}
	\end{align}

where :math:`H_0` is the free field Hamiltonian, also known as the symmetry generator for the time translation, or the energy operator. However, rather than thinking of it as an abstract operator as we've done so far, we'll (momentarily) make it a functional of canonical variables. With this in mind, we can take the time derivative of :math:`\eqref{eq_free_field_q_time_evolution}` and :math:`\eqref{eq_free_field_p_time_evolution}` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		\dot{q}_n(t, \xbf) &= \ifrak \left[ H_0, q_n(t, \xbf) \right] &&= \frac{\delta H_0}{\delta p_n(t, \xbf)}
		\label{eq_free_field_hamilton_equation_q_dot} \\
		\dot{p}_n(t, \xbf) &= \ifrak \left[ H_0, p_n(t, \xbf) \right] &&= -\frac{\delta H_0}{\delta q_n(t, \xbf)}
		\label{eq_free_field_hamilton_equation_p_dot}
	\end{alignat}

We recognize these as the quantum analog of `Hamilton's equation of motion <https://en.wikipedia.org/wiki/Hamiltonian_mechanics>`__.

To turn :math:`H_0` into a functional of canonical variables, we first make it a functional of creation and annihilation operators. Remembering that :math:`H_0` is the energy operator, and :math:`p_0 = \sqrt{\pbf^2 + m^2}` is the energy in the :math:`4`-momentum, we can write :math:`H_0` as a diagonal operator as follows

.. math::
	:nowrap:

	\begin{equation}
		H_0 = \sum_{n, \sigma} \int d^3 p~a^{\dagger}(\pbf, \sigma, n) a(\pbf, \sigma, n) \sqrt{\pbf^2 + m^2}
		\label{eq_free_field_hamiltonian_diagonal}
	\end{equation}

For simplicity, let's consider the case of a real scalar field :math:`\psi(x)` given by :math:`\eqref{eq_scalar_field_psi_by_creation_and_annihilation_operators}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		q(t, \xbf) = \psi(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} \left( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)
	\end{equation*}

The canonical conjugate variable is

.. math::
	:nowrap:

	\begin{equation*}
		p(t, \xbf) = \dot{\psi}(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} (-\ifrak p_0) \left( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)
	\end{equation*}

These look a bit far from :math:`\eqref{eq_free_field_hamiltonian_diagonal}`. But since :math:`H_0` involves products like :math:`a^{\dagger}(\pbf, \sigma, n) a(\pbf, \sigma, n)`, let's try to square the canonical variables as follows

.. math::
	:nowrap:

	\begin{align*}
		\int d^3 x~q^2(t, \xbf) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}}
				\Big( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big)
				\Big( e^{\ifrak p' \cdot x} a(\pbf') + e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
			&= \int \frac{d^3 p}{2p_0} \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } + \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)
	\end{align*}

and

.. math::
	:nowrap:

	\begin{align*}
		\int d^3 x~p^2(t, \xbf) &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}} (-p_0 p'_0)
				\Big( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big)
				\Big( e^{\ifrak p' \cdot x} a(\pbf') - e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
			&= \int \frac{d^3 p}{2p_0} \left( -p_0^2 \right) \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } - \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)
	\end{align*}

and finally, inspired by the calculations above

.. math::
	:nowrap:

	\begin{align*}
		\int d^3 x~\left( \nabla q(t, \xbf) \right)^2 &= \frac{1}{(2\pi)^3} \int \frac{d^3 p~d^3 p'~d^3 x}{2\sqrt{p_0 p'_0}} \left( -\pbf \cdot \pbf' \right)
				\Big( e^{\ifrak p \cdot x} a(\pbf) - e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \Big)
				\Big( e^{\ifrak p' \cdot x} a(\pbf') - e^{-\ifrak p' \cdot x} a^{\dagger}(\pbf') \Big) \\
			&= \int \frac{d^3 p}{2p_0} \pbf^2 \left( \blue{ e^{-2\ifrak p_0 t} a(\pbf) a(-\pbf) + e^{2\ifrak p_0 t} a^{\dagger}(\pbf) a^{\dagger}(-\pbf) } + \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+ \right)
	\end{align*}

Putting these calculations together in a specific way, and using the identity :math:`p_0^2 - \pbf^2 = m^2`, we can eliminate the blue terms as follows

.. math::
	:nowrap:

	\begin{align}
		\frac{1}{2} \int d^3 x \left( p^2 + \left( \nabla q \right)^2 + m^2 q^2 \right)
			&= \frac{1}{2} \int d^3 p~p_0 \left[ a(\pbf), a^{\dagger}(\pbf) \right]_+
			\label{eq_calculate_free_real_scalar_field_hamiltonian} \\
			&= \int d^3 p~p_0 \left( a^{\dagger}(\pbf) a(\pbf) + \frac{1}{2} \delta^3(\pbf-\pbf) \right)
			\nonumber \\
			&= H_0 + \blue{ \frac{1}{2} \int d^3 p~p_0 \delta^3(0) }
			\nonumber
	\end{align}

Here we've encountered for the first time an infinite term (which we've marked in blue). As long as the Hamiltonian dynamics :math:`\eqref{eq_free_field_hamilton_equation_q_dot}` -- :math:`\eqref{eq_free_field_hamilton_equation_p_dot}` is concerned, it makes no difference adding a constant to the Hamiltonian. Hence we can write the free Hamiltonian for real scalar fields as follows

.. math::
	:nowrap:

	\begin{equation}
		H_0^{\text{RSF}} = \frac{1}{2} \int d^3 x \left( p^2 + \left( \nabla q \right)^2 + m^2 q^2 \right)
		\label{eq_free_scalar_field_hamiltonian}
	\end{equation}

.. warning::
	Throwing away the infinite term in :math:`\eqref{eq_calculate_free_real_scalar_field_hamiltonian}` is an instance of a well-known criticism in quantum field theory: "just because something is infinite doesn't mean it's zero". Indeed, Weinberg mentioned in page 297 [Wei95]_ that such "infinities" shouldn't be thrown away when, for example, the fields are constrained within a finite space, or there is an involvement of gravity.

Now it's time to introduce the rather mysterious Lagrangian, which can be derived from the Hamiltonian via the so-called `Legendre transformation <https://en.wikipedia.org/wiki/Legendre_transformation>`__ as follows

.. math::
	:nowrap:

	\begin{equation}
		L_0\left[ q(t), \dot{q}(t) \right] \coloneqq \sum_n \int d^3 x~p_n(t, \xbf) \dot{q}_n(t, \xbf) - H_0
		\label{eq_legendre_transformation_lagrangian_from_hamiltonian}
	\end{equation}

where each occurrence of :math:`p_n(t)` is replaced by its expression in :math:`q_n(t)` and :math:`\dot{q}_n(t)`.

As a concrete example, let's consider again the real scalar field, where :math:`p = \dot{q}`. It follows that

.. math::
	:nowrap:

	\begin{align}
		L_0^{\text{RSF}} &= \int d^3 x \left( p\dot{q} - \frac{1}{2} p^2 - \frac{1}{2} \left( \nabla q \right)^2 - \frac{1}{2} m^2 q^2 \right)
			\label{eq_free_real_scalar_field_lagrangian} \\
			&= \frac{1}{2} \int d^3 x \left( \dot{q}^2 - \left( \nabla q \right)^2 - m^2 q^2 \right)
			\nonumber \\
			&= -\frac{1}{2} \int d^3 x \left( \p_{\mu} \psi \p^{\mu} \psi + m^2 \psi^2 \right)
			\nonumber
	\end{align}

It should be noted that expressing :math:`p` in terms of :math:`q` and :math:`\dot{q}` isn't always easy. Indeed, it's far from obvious how the :math:`p_i` defined by :math:`\eqref{eq_defn_p_vector_field_self_dual}` could be expressed in the corresponding :math:`q_i` and :math:`\dot{q}_i`. (Un)Fortunately, we'd never really need to do so -- writing down a Lagrangian turns out to be mostly a guess work.

.. _sec_hamiltonian_and_lagrangian_for_interacting_fields:

Hamiltonian and Lagrangian for interacting fields
+++++++++++++++++++++++++++++++++++++++++++++++++

Let :math:`H` be the full Hamiltonian. Then the Heisenberg picture canonical variables can be defined as follows

.. math::
	:nowrap:

	\begin{align}
		Q_n(t, \xbf) &\coloneqq e^{\ifrak Ht} q_n(0, \xbf) e^{-\ifrak Ht}
		\label{eq_defn_heisenberg_canonical_q} \\
		P_n(t, \xbf) &\coloneqq e^{\ifrak Ht} p_n(0, \xbf) e^{-\ifrak Ht}
		\nonumber
	\end{align}

Then obviously these canonical variables also satisfy the canonical (anti-)commutation relations

.. math::
	:nowrap:

	\begin{align*}
		\left[ Q_n(t, \xbf), P_{n'}(t, \ybf) \right]_{\pm} &= \ifrak \delta^3(\xbf-\ybf) \delta_{n n'} \\
		\left[ Q_n(t, \xbf), Q_{n'}(t, \ybf) \right]_{\pm} &= 0 \\
		\left[ P_n(t, \xbf), P_{n'}(t, \ybf) \right]_{\pm} &= 0
	\end{align*}

Moreover, the analog of :math:`\eqref{eq_free_field_hamilton_equation_q_dot}` and :math:`\eqref{eq_free_field_hamilton_equation_p_dot}` holds as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		\dot{Q}_n(t, \xbf) &= \ifrak \left[ H, Q_n(t, \xbf) \right] &&= \frac{\delta H}{\delta P_n(t, \xbf)} \\
		\dot{P}_n(t, \xbf) &= \ifrak \left[ H, P_n(t, \xbf) \right] &&= -\frac{\delta H}{ \delta Q_n(t, \xbf)}
	\end{alignat*}

As an example, we note that, in light of :math:`\eqref{eq_free_scalar_field_hamiltonian}`, the full Hamiltonian for real scalar fields may be written as

.. math::
	:nowrap:

	\begin{equation*}
		H^{RSF} = \int d^3 x \left( \tfrac{1}{2} P^2 + \tfrac{1}{2} \left( \nabla Q \right)^2 + \tfrac{1}{2} m^2 Q^2 + \Hscr(Q) \right)
	\end{equation*}

where :math:`\Hscr(Q)` is the perturbation term giving rise to the interaction.

The Lagrangian formalism
^^^^^^^^^^^^^^^^^^^^^^^^

We'll leave aside the discussion of canonical variables for a bit to introduce the Lagrangian formalism in its most general form. After that we'll play the game backwards. Namely, instead of constructing canonical variables out of the free fields that we've been exclusively considering since :ref:`sec_quantum_fields_and_antiparticles`, we'll get canonically conjugate fields out of the (magically appearing) Lagrangians, and then *impose* the canonical commutation relations :math:`\eqref{eq_canonical_commutation_relation_1}` -- :math:`\eqref{eq_canonical_commutation_relation_3}` on them -- a procedure generally known as "quantization".

In the classical physical theory of fields, a Lagrangian is a functional :math:`L[\Psi(t), \dot{\Psi}(t)]`, where :math:`\Psi(t)` is any field and :math:`\dot{\Psi}(t)` is its time derivative. Here we've capitalized the field variables to distinguish them from the free fields considered in the previous section. Define the conjugate fields as follows

.. math::
	:nowrap:

	\begin{equation}
		\Pi_n(t, \xbf) \coloneqq \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)}
		\label{eq_general_lagrangian_conjugate_pi}
	\end{equation}

so that the field equations are given by

.. math::
	:nowrap:

	\begin{equation}
		\dot{\Pi}_n(t, \xbf) = \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi_n(t, \xbf)}
		\label{eq_equation_of_motion_for_fields}
	\end{equation}

.. warning::
	Unlike the functional derivatives considered in :math:`\eqref{eq_infinitesimal_variation_of_functional_of_canonical_variables}` for canonical variables, the functional derivative :math:`\eqref{eq_general_lagrangian_conjugate_pi}`, interpreted quantum mechanically, is not really well-defined since :math:`\Psi(t)` and :math:`\dot{\Psi}(t)` don't in general satisfy a simple (same time) commutation relation. According to Weinberg (see footnote on page 299 in [Wei95]_), "no important issues hinge on the details here". So we'll pretend that it behaves just like usual derivatives.

Indeed, recall that in the classical Lagrangian formalism, the field equations are given by a variational principle applied to the so-called *action*, defined as follows

.. math::
	:nowrap:

	\begin{equation}
		I[\Psi] \coloneqq \int_{-\infty}^{\infty} dt~L[\Psi(t), \dot{\Psi}(t)]
		\label{eq_defn_action_of_fields}
	\end{equation}

The infinitesimal variation of :math:`I[\Psi]` is given by

.. math::
	:nowrap:

	\begin{align*}
		\delta I[\Psi] &= \sum_n \int_{-\infty}^{\infty} dt \int d^3 x \left(
				\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi_n(t, \xbf)} \delta \Psi_n(t, \xbf) +
				\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)} \delta \dot{\Psi}_n(t, \xbf) \right) \\
			&= \sum_n \int_{-\infty}^{\infty} dt \int d^3 x \left(
				\frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \Psi(t, \xbf)} - \frac{d}{dt} \frac{\delta L[\Psi(t), \dot{\Psi}(t)]}{\delta \dot{\Psi}_n(t, \xbf)} \right) \delta \Psi_n(t, \xbf)
	\end{align*}

where for the last equality, integration by parts is used under the assumption that the infinitesimal variation :math:`\delta \Psi_n(t, \xbf)` vanishes at :math:`t \to \pm\infty`. Obviously :math:`\delta I[\Psi]` vanishes for any :math:`\delta \Psi_n(t, \xbf)` if and only if :math:`\eqref{eq_equation_of_motion_for_fields}` is satisfied.

Now we're interested in constructing Lorentz invariant theories, but an action defined by :math:`\eqref{eq_defn_action_of_fields}` apparently distinguishes the time from space variables. This motivates the hypothesis that the Lagrangian itself is given by a spatial integral of a so-called *Lagrangian density* as follows

.. math::
	:nowrap:

	\begin{equation}
		L[\Psi(t), \dot{\Psi}(t)] = \int d^3 x~\Lscr(\Psi(t, \xbf), \nabla\Psi(t, \xbf), \dot{\Psi}(t, \xbf))
		\label{eq_defn_lagrangian_density}
	\end{equation}

In terms of the Lagrangian density, we can rewrite the action :math:`\eqref{eq_defn_action_of_fields}` as a :math:`4`-integral as follows

.. math::
	:nowrap:

	\begin{equation*}
		I[\Psi] = \int d^4 x~\Lscr(\Psi(x), \p_{\mu} \Psi(x))
	\end{equation*}

We'd also like to reexpress the field equations :math:`\eqref{eq_equation_of_motion_for_fields}` in terms of the Lagrangian density. To this end, let's first calculate the variation of :math:`\eqref{eq_defn_lagrangian_density}` by an amount :math:`\delta \Psi_n(t, \xbf)` as follows

.. math::
	:nowrap:

	\begin{align*}
		\delta L &= \sum_n \int d^3 x \left( \frac{\delta\Lscr}{\delta\Psi_n} \delta\Psi_n + \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \cdot \nabla \delta\Psi_n + \frac{\delta\Lscr}{\delta\dot{\Psi}_n} \delta\dot{\Psi}_n \right) \\
			&= \sum_n \int d^3 x \left( \left( \frac{\delta\Lscr}{\delta\Psi_n} - \nabla \cdot \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \right) \delta\Psi_n + \frac{\delta\Lscr}{\delta\dot{\Psi}_n} \delta\dot{\Psi}_n \right)
	\end{align*}

It follows that

.. math::
	:nowrap:

	\begin{align*}
		\frac{\delta L}{\delta\Psi_n} &= \frac{\delta\Lscr}{\delta\Psi_n} - \nabla \cdot \frac{\delta\Lscr}{\delta(\nabla \Psi_n)} \\
		\frac{\delta L}{\delta\dot{\Psi}_n} &= \frac{\delta \Lscr}{\delta \dot{\Psi}_n}
	\end{align*}

Combining these with :math:`\eqref{eq_equation_of_motion_for_fields}` and :math:`\eqref{eq_equation_of_motion_for_fields}`, we've derived the so-called Euler-Lagrange equations for the Lagrangian density

.. math::
	:nowrap:

	\begin{equation}
		\frac{\delta \Lscr}{\delta \Psi_n} = \p_{\mu} \frac{\delta \Lscr}{\delta(\p_{\mu} \Psi_n)}
		\label{eq_euler_lagrange}
	\end{equation}

Note that the summing :math:`4`-index :math:`\mu` here represents :math:`x_{\mu}`. Most importantly, the field equations given by :math:`\eqref{eq_euler_lagrange}` will be Lorentz invariant if :math:`\Lscr` is. Indeed, guessing such :math:`\Lscr` will be more or less the only way to construct Lorentz invariant (quantum) field theories.

.. note::
	The Lagrangian density :math:`\Lscr` is assumed to be real for two reasons. First, if :math:`\Lscr` were complex, then splitting it into the real and imaginary parts, :math:`\eqref{eq_euler_lagrange}` would contain twice as many equations as there are fields, regardless whether real or complex. This is undesirable because generically there will be no solutions. The second reason has to wait until the next section, where symmetries will be discussed. It turns out that the reality of :math:`\Lscr` will guarantee that the symmetry generators are Hermitian.

Now recall from the previous section that the anchor of our knowledge is the Hamiltonian -- we know how it must look like, at least for free fields. To go from the Lagrangian to the Hamiltonian, we use again the Legendre transformation (cf. :math:`\eqref{eq_legendre_transformation_lagrangian_from_hamiltonian}`) to *define* the Hamiltonian as follows

.. math::
	:nowrap:

	\begin{equation}
		H[\Psi, \Pi] \coloneqq \sum_n \int d^3 x~\Pi_n(t, \xbf) \dot{\Psi}_n(t, \xbf) - L[\Psi(t), \dot{\Psi}(t)]
		\label{eq_legendre_transformation_hamiltonian_from_lagrangian}
	\end{equation}

.. warning::
	In order to realize :math:`H` as a functional of :math:`\Psi` and :math:`\Pi`, one must in principle be able to solve for :math:`\dot{\Psi}_n` in terms of :math:`\Psi_n` and :math:`\Pi_n` from :math:`\eqref{eq_general_lagrangian_conjugate_pi}`. This isn't always easy, if at all possible, but it rarely pose serious difficulties in applications either.

As a double check, let's verify that the Hamiltonian defined by :math:`\eqref{eq_legendre_transformation_hamiltonian_from_lagrangian}` also satisfies Hamilton's equations (cf. :math:`\eqref{eq_free_field_hamilton_equation_q_dot}` -- :math:`\eqref{eq_free_field_hamilton_equation_p_dot}`). Indeed, the variational derivatives are calculated as follows

.. math::
	:nowrap:

	\begin{align*}
		\frac{\delta H}{\delta \Pi_n(t, \xbf)} &= \sum_m \int d^3 y \left( \frac{\delta \Pi_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \dot{\Psi}_m(t, \ybf) + \Pi_m(t, \ybf) \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \right) - \sum_m \int d^3 y \frac{\delta L}{\delta \dot{\Psi}_m(t, \ybf)} \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Pi_n(t, \xbf)} \\
			&= \sum_m \int d^3 y~\delta_{m,n} \delta^3(\ybf-\xbf) \dot{\Psi}_m(t, \ybf) \\
			&= \dot{\Psi}_n(t, \xbf) \\
		\frac{\delta H}{\delta \Psi_n(t, \xbf)} &= \sum_m \int d^3 y~\Pi_m(t, \ybf) \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} - \sum_m \int d^3 y \left( \frac{\delta L}{\delta \Psi_m(t, \ybf)} \frac{\delta \Psi_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} + \frac{\delta L}{\delta \dot{\Psi}_m(t, \ybf)} \frac{\delta \dot{\Psi}_m(t, \ybf)}{\delta \Psi_n(t, \xbf)} \right) \\
			&\xlongequal{\eqref{eq_general_lagrangian_conjugate_pi},~\eqref{eq_equation_of_motion_for_fields}} -\sum_m \int d^3 y~\delta_{m, n} \delta^3(\ybf-\xbf) \dot{\Pi}_m(t, \ybf) \\
			&= -\dot{\Pi}_n(t, \xbf)
	\end{align*}

It's therefore attempting to demand, in the Lagrangian formalism, that :math:`\Psi_n` and :math:`\Pi_n`, defined by :math:`\eqref{eq_general_lagrangian_conjugate_pi}`, satisfy the canonical commutation relations. In other words, they are (Heisenberg picture) canonically conjugate fields. But this is not true in general, as it turns out.

The issue is that the Lagrangian :math:`L[\Psi(t), \dot{\Psi}(t)]` may contain certain field, but not its time derivative. One example is spin-:math:`1` vector fields, where we see from :math:`\eqref{eq_defn_q_vector_field_self_dual}` that the spatial fields :math:`\psi_i` are part of the canonical variables, but not :math:`\psi_0`, which nonetheless should present in the Lagrangian by Lorentz invariance. It turns out that what's missing from the Lagrangian is :math:`\dot{\psi}_0`, which causes its conjugate variable defined by :math:`\eqref{eq_general_lagrangian_conjugate_pi}` to vanish.

But instead of dealing with vector fields further, we'll turn back to the general ground to establish the fundamental principles. Inspired by above discussion, we can rewrite the Lagrangian as

.. math::
	:nowrap:

	\begin{equation}
		L[Q(t), \dot{Q}(t), C(t)]
		\label{eq_general_quantum_lagrangian}
	\end{equation}

where each :math:`Q_n(t)` has a corresponding :math:`\dot{Q}_n(t)`, but not for :math:`C(t)`. It follows that one can define the canonical conjugates by

.. math::
	:nowrap:

	\begin{equation*}
		P_n(t, \xbf) \coloneqq \frac{\delta L[Q(t), \dot{Q}(t), C(t)]}{\delta \dot{Q}_n(t, \xbf)}
	\end{equation*}

and hence the Hamiltonian takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		H[Q, P] = \sum_n \int d^3 x~P_n \dot{Q}_n - L[Q(t), \dot{Q}(t), C(t)]
	\end{equation*}

.. dropdown:: Quantization of free scalar fields
	:animate: fade-in-slide-down
	:icon: unlock

	We'll illustrate how "quantization" works in the simplest case of free scalar fields :math:`\phi(t, \xbf)`. Namely, we'll reverse our earlier approach to the quantum theory by starting from a Lagrangian, and then work out the field equations, solve them for the fields, impose canonical commutation relations, and finally arrive at the familiar commutation relations between creation and annihilation operators introduced in :ref:`sec_the_cluster_decomposition_principle`.

	Following :math:`\eqref{eq_free_real_scalar_field_lagrangian}`, let's consider the following Lagrangian

	.. math::
		:nowrap:

		\begin{equation*}
			L_0[\phi, \dot{\phi}] = -\frac{1}{2} \int d^3 x \left( \p_{\mu} \phi \p^{\mu} \phi + m^2 \phi^2 \right)
		\end{equation*}

	The canonical conjugate of :math:`\phi(t, \xbf)` is then

	.. math::
		:nowrap:

		\begin{equation}
			\pi(t, \xbf) \coloneqq \frac{\delta L_0}{\delta \dot{\phi}(t, \xbf)} = \dot{\phi}(t, \xbf)
			\label{eq_free_scalar_field_pi_equals_dot_phi}
		\end{equation}

	Hence the Hamiltonian takes the following form

	.. math::
		:nowrap:

		\begin{align*}
			H_0[\phi, \pi] &= \int d^3 x~\pi(t, \xbf) \dot{\phi}(t, \xbf) - L_0 \\
				&= \frac{1}{2} \int d^3 x \left( \pi^2(t, \xbf) + \big( \nabla \phi(t, \xbf) \big)^2 + m^2 \phi^2 \right)
		\end{align*}

	The field equations :math:`\eqref{eq_free_field_hamilton_equation_q_dot}` -- :math:`\eqref{eq_free_field_hamilton_equation_p_dot}` are then given as follows

	.. math::
		:nowrap:

		\begin{alignat*}{2}
			\dot{\phi}(t, \xbf) &= \frac{\delta H_0}{\delta \pi(t, \xbf)} &&= \pi(t, \xbf) \\
			\dot{\pi}(t, \xbf) &= -\frac{\delta H_0}{\delta \phi(t, \xbf)} &&= \nabla^2 \phi(t, \xbf) - m^2 \phi(t, \xbf)
		\end{alignat*}

	Together, it implies that the field equation is precisely the Klein-Gordon equation

	.. math::
		:nowrap:

		\begin{equation*}
			\left( \square - m^2 \right) \phi(x) = 0
		\end{equation*}

	Using Fourier transform, the general Hermitian solution, up to a scalar, can be written as follows

	.. math::
		:nowrap:

		\begin{equation*}
			\phi(x) = \frac{1}{(2\pi)^{3/2}} \int \frac{d^3 p}{\sqrt{2p_0}} \left( e^{\ifrak p \cdot x} a(\pbf) + e^{-\ifrak p \cdot x} a^{\dagger}(\pbf) \right)
		\end{equation*}

	where :math:`p_0 = \sqrt{\pbf^2+m^2}` and :math:`a(\pbf)` is, at the moment, just any operator function of :math:`\pbf`. Using :math:`\eqref{eq_free_scalar_field_pi_equals_dot_phi}` we have

	.. math::
		:nowrap:

		\begin{equation*}
			\pi(x) = \frac{-\ifrak}{(2\pi)^{3/2}} \int d^3 p \sqrt{\frac{p_0}{2}} \left( e^{\ifrak p \cdot x}a(\pbf) - e^{-\ifrak p \cdot x}a^{\dagger}(\pbf) \right)
		\end{equation*}

	One can then verify that if we impose the canonical commutation relations on the conjugate fields :math:`\phi(t, \xbf)` and :math:`\pi(t, \xbf)` as follows

	.. math::
		:nowrap:

		\begin{align*}
			\left[ \phi(t, \xbf), \pi(y, \ybf) \right] &= \ifrak \delta^3(\xbf-\ybf) \\
			\left[ \phi(t, \xbf), \phi(t, \ybf) \right] &= 0 \\
			\left[ \pi(t, \xbf), \pi(t, \ybf) \right] &= 0
		\end{align*}

	then the familiar commutation relations

	.. math::
		:nowrap:

		\begin{align*}
			\left[ a(\pbf), a^{\dagger}(\pbf') \right] &= \delta^3(\pbf-\pbf') \\
			\left[ a(\pbf), a(\pbf') \right] &= 0
		\end{align*}

	must hold. In this way we've completely reversed the process of deriving a Lagrangian from free fields made up of creation and annihilation operators.


Global symmetries
^^^^^^^^^^^^^^^^^

Of course, the reason for introducing the Lagrangian formalism is not to reproduce the Hamiltonians and the fields that we already knew. The main motivation is that, as we'll see, the Lagrangian formalism provides a framework for studying symmetries. Recall from :ref:`sec_what_is_a_symmetry` that a symmetry was defined to be a(n anti-)unitary transformation on the Hilbert space of states, i.e., a transformation that preserves amplitudes. Now in the Lagrangian formalism, field equations come out of the stationary action condition. Therefore in this context, we'll redefine a symmetry as an infinitesimal variation of the fields that leaves the action invariant. As it turns out, symmetries in this sense lead to conserved currents, which are nothing but the symmetry operators considered earlier. Hence besides a slight abuse of terminology, the notion of symmetries will be consistent.

.. note::
	Throughout this section, repeated indexes like :math:`n`, which are used to index various fields, in an equation are not automatically summed up. On the other hand, repeated :math:`4`-indexes like :math:`\mu` do follow the Einstein summation convention.

Consider an infinitesimal variation

.. math::
	:nowrap:

	\begin{equation}
		\Psi_n(x) \to \Psi_n(x) + \ifrak \epsilon \Fscr_n(x)
		\label{eq_infinitesimal_variation_of_field}
	\end{equation}

which leaves the action :math:`I[\Psi]`  invariant

.. math::
	:nowrap:

	\begin{equation}
		0 = \delta I = \ifrak \epsilon \sum_n \int d^4 x~\frac{\delta I[\Psi]}{\delta \Psi_n(x)} \Fscr_n(x)
		\label{eq_vanishing_of_action_under_infinitesimal_variation}
	\end{equation}

A few remarks are in order. First of all, if we think of :math:`\eqref{eq_infinitesimal_variation_of_field}` as an infinitesimal (unitary) symmetry transformation, then the coefficient :math:`\ifrak` can be justified by then intention of making :math:`\Fscr_n(x)` Hermitian. Next, although :math:`\eqref{eq_vanishing_of_action_under_infinitesimal_variation}` *always* holds when :math:`\Psi_n(x)` is stationary, the infinitesimal :math:`\Fscr_n(x)` being a symmetry demands that :math:`\eqref{eq_vanishing_of_action_under_infinitesimal_variation}` holds true for *any* :math:`\Psi_n(x)`. Finally, we emphasize the fact that :math:`\epsilon` is an infinitesimal *constant*, rather than a function of :math:`x`, is the defining property for the symmetry to be called "global". Indeed, we'll be dealing with symmetries that are not global in the next chapter, namely, the gauge symmetries.

.. _sec_from_symmetries_to_conservation_laws:

From symmetries to conservation laws
++++++++++++++++++++++++++++++++++++

The general principle that "symmetries imply conservation laws" is mathematically known as `Noether's theorem <https://en.wikipedia.org/wiki/Noether%27s_theorem>`__, but we'll not bother with any mathematical formality here. To see how to derive conserved quantitites from an assumed symmetry, let's change :math:`\eqref{eq_infinitesimal_variation_of_field}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\Psi_n(x) \to \Psi_n(x) + \ifrak \epsilon(x) \Fscr_n(x)
		\label{eq_functional_infinitesimal_variation_of_field}
	\end{equation}

where :math:`\epsilon(x)` now is an infinitesimal function of :math:`x`. Under this variation, the corresponding :math:`\delta I` may not vanish. But it must take the following form

.. math::
	:nowrap:

	\begin{equation}
		\delta I = -\int d^4 x J^{\mu}(x) \p_{\mu} \epsilon(x)
		\label{eq_variation_of_action_by_functional_deformation}
	\end{equation}

because it must vanishe when :math:`\epsilon(x)` is constant. Here :math:`J^{\mu}(x)` is a function(al) to be determined in individual cases, and is usually known as *current*. Now if :math:`\Psi_n(x)` satisfies the field equations, i.e., it's a stationary point of the action, then :math:`\eqref{eq_variation_of_action_by_functional_deformation}` must vanishes for any :math:`\epsilon(x)`. Applying integration by parts (and assuming :math:`\Fscr_n(x)` vanishes at infinity), we must have

.. math::
	:nowrap:

	\begin{equation}
		\p_{\mu} J^{\mu}(x) = 0
		\label{eq_general_conservation_of_current}
	\end{equation}

which is the conservation law for :math:`J`, which then can be called a conserved current. One gets also a conserved quantity, i.e., a quantity that doesn't change by time, by integrating :math:`\eqref{eq_general_conservation_of_current}` over the :math:`3`-space as follows

.. math::
	:nowrap:

	\begin{equation*}
		\dot{J}^0(x) = -\nabla \cdot \Jbf(x)
			\implies \int d^3 x~\dot{J}^0(x) = -\int d^3 x~\nabla \cdot \Jbf(x) = 0
			\implies F \coloneqq \int d^3 x~J^0(x) \text{ is conserved.}
	\end{equation*}

Unfortunately, not much more can be said about the conserved current :math:`J` at this level of generality. This is, however, not the case if one imposes stronger assumptions on the symmetry, as we now explain.

Lagrangian-preserving symmetry
	This is the first strengthening of the symmetry assumption. Namely, instead of assuming that the variation :math:`\eqref{eq_infinitesimal_variation_of_field}` fixes the action, we assume that it fixes the Lagrangian itself. Namely,

	.. math::
		:nowrap:

		\begin{equation}
			\delta L = \ifrak \epsilon \sum_n \int d^3 x \left( \frac{\delta L}{\delta \Psi_n(t, \xbf)} \Fscr_n(t, \xbf) + \frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \dot{\Fscr}_n(t, \xbf) \right) = 0
			\label{eq_stationary_lagrangian}
		\end{equation}

	Now let :math:`\epsilon(t)` be a time-dependent infinitesimal in :math:`\eqref{eq_functional_infinitesimal_variation_of_field}`. Then we can calculate :math:`\delta I` under such variation as follows

	.. math::
		:nowrap:

		\begin{align*}
			\delta I &= \ifrak \sum_n \int dt \int d^3 x \left( \frac{\delta L}{\delta \Psi_n(t, \xbf)} \epsilon(t) \Fscr_n(t, \xbf) + \frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \frac{d}{dt} \big( \epsilon(t) \Fscr_n(t, \xbf) \big) \right) \\
				&= \ifrak \sum_n \int dt \int d^3 x~\frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \dot{\epsilon}(t) \Fscr_n(t, \xbf)
		\end{align*}

	Comparing with :math:`\eqref{eq_variation_of_action_by_functional_deformation}`, we can derive an explicit formula for the conserved quantity as follows

	.. math::
		:nowrap:

		\begin{equation}
			F = -\ifrak \sum_n \int d^3 x~\frac{\delta L}{\delta \dot{\Psi}_n(t, \xbf)} \Fscr_n(t, \xbf)
			\label{eq_lagrangian_preserving_symmetry_conserved_quantity}
		\end{equation}

	Indeed, one can verify directly that :math:`\dot{F}(t) = 0` using :math:`\eqref{eq_stationary_lagrangian}` together with the field equations :math:`\eqref{eq_general_lagrangian_conjugate_pi}` and :math:`\eqref{eq_equation_of_motion_for_fields}`.

Lagrangian-density-preserving symmetry
	Taking the previous assumption further, let's impose the even stronger condition that the Lagrangian density is invariant under :math:`\eqref{eq_infinitesimal_variation_of_field}`. It means that

	.. math::
		:nowrap:

		\begin{equation}
			\delta \Lscr = \ifrak \epsilon \sum_n \left( \frac{\delta \Lscr}{\delta \Psi_n(x)} \Fscr_n(x) + \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \p_{\mu} \Fscr_n(x) \right) = 0
			\label{eq_stationary_lagrangian_density}
		\end{equation}

	Now under :math:`\eqref{eq_functional_infinitesimal_variation_of_field}`, we can calculate the variation of the action as follows

	.. math::
		:nowrap:

		\begin{align*}
			\delta I &= \ifrak \sum_n \int d^4 x~\left( \frac{\delta \Lscr}{\delta \Psi_n(x)} \epsilon(x) \Fscr_n(x) + \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \p_{\mu} \big( \epsilon(x) \Fscr_n(x) \big) \right) \\
				&= \ifrak \sum_n \int d^4 x~\frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \Fscr_n(x) \p_{\mu}\epsilon(x)
		\end{align*}

	Comparing with :math:`\eqref{eq_variation_of_action_by_functional_deformation}` as before, we can derive an explicit formula for the conserved current as follows

	.. math::
		:nowrap:

		\begin{equation}
			J^{\mu}(x) = -\ifrak \sum_n \frac{\delta \Lscr}{\delta (\p_{\mu} \Psi_n(x))} \Fscr_n(x)
			\label{eq_lagrangian_density_preserving_symmetry_conserved_density}
		\end{equation}

	Once again, one can directly verify that :math:`\p_{\mu} J^{\mu}(x) = 0` using :math:`\eqref{eq_stationary_lagrangian_density}` together with the Euler-Lagrange equation :math:`\eqref{eq_euler_lagrange}`.

So far everything has been completely classical. To make it a quantum theory, we'll involve the canonical fields introduced in :ref:`sec_hamiltonian_and_lagrangian_for_interacting_fields`. More precisely, instead of any :math:`\Fscr_n(t, \xbf)`, we'll suppose that it takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		\Fscr_n(Q(t), \xbf)
	\end{equation*}

where :math:`Q(t)` is defined by :math:`\eqref{eq_defn_heisenberg_canonical_q}`. Next, recall from :math:`\eqref{eq_general_quantum_lagrangian}` that the field :math:`\Psi_n` is either a :math:`Q_n`, in which case :math:`\delta L / \delta \dot{Q}_n = P_n`, or a :math:`C_n`, in which case the functional derivative vaninshes.

Now in the case of a Lagrangian-preserving symmetry, the conserved quantity :math:`\eqref{eq_lagrangian_preserving_symmetry_conserved_quantity}` takes the following form

.. math::
	:nowrap:

	\begin{equation}
		F = -\ifrak \sum_n \int d^3 x~P_n(t, \xbf) \Fscr_n(Q(t), \xbf)
		\label{eq_lagrangian_preserving_symmetry_generator_formula}
	\end{equation}

which of course is time-independent. Moreover, one can show that :math:`F` in fact generates the quantum symmetry in the following sense

.. math::
	:nowrap:

	\begin{equation}
		\left[ F, Q_n(t, \xbf) \right] = -\ifrak \sum_m \int d^3 y~\left[ P_m(t, \ybf), Q_n(t, \xbf) \right] \Fscr_m(Q(t), \ybf) = -\Fscr_n(Q(t), \xbf)
		\label{eq_lagrangian_preserving_symmetry_generator}
	\end{equation}

where we've taken advantage of the time-independency of :math:`F` to arrange the same-time commutator.

.. _sec_spacetime_translations:

Spacetime translations
++++++++++++++++++++++

So far the symmetries have been rather abstract, to make it more explicit, and also to get warmed up for the general case, let's assume the Lagrangian is invariant under the (spacetime) translation transformation given as follows

.. math::
	:nowrap:

	\begin{equation*}
		\Psi_n(x) \to \Psi_n(x + \epsilon) = \Psi_n(x) + \epsilon^{\mu} \p_{\mu} \Psi_n(x)
	\end{equation*}

Comparing with :math:`\eqref{eq_infinitesimal_variation_of_field}` we see that

.. math::
	:nowrap:

	\begin{equation*}
		\Fscr_{\mu} = -\ifrak \p_{\mu} \Psi_n
	\end{equation*}

It follows from :math:`\eqref{eq_variation_of_action_by_functional_deformation}` and :math:`\eqref{eq_general_conservation_of_current}` that there exists a conserved :math:`4`-current :math:`T^{\mu}_{\nu}`, which is known as the `energy-momentum tensor <https://en.wikipedia.org/wiki/Stress%E2%80%93energy_tensor>`__, such that

.. math::
	:nowrap:

	\begin{equation*}
		\p_{\mu} T^{\mu}_{\nu} = 0
	\end{equation*}

The corresponding conserved currents then take the form

.. math::
	:nowrap:

	\begin{equation}
		P_{\nu} \coloneqq \int d^3 x~T^0_{\nu}
		\label{eq_spacetime_translation_conserved_quantity_is_momentum}
	\end{equation}

such that :math:`\dot{P}_{\nu} = 0`. Here it's important to not confuse :math:`P_{\nu}` with a canonical variable -- it's just a conserved quantity which turns out to be the :math:`4`-momentum.

Now recall from :math:`\eqref{eq_defn_lagrangian_density}` that the Lagrangian is usually the spatial integral of a density functional. Hence it's not unreasonable to suppose that the Lagrangian is indeed invariant under spatial translations. Under this assumption, we can rewrite :math:`\eqref{eq_lagrangian_preserving_symmetry_generator_formula}` as follows

.. math::
	:nowrap:

	\begin{equation}
		\Pbf \coloneqq -\sum_n \int d^3 x~P_n(t, \xbf) \nabla Q_n(t, \xbf)
		\label{eq_spatial_translation_conserved_quantity}
	\end{equation}

with the understanding that :math:`\Psi_n = Q_n`.

To verify that :math:`\Pbf` indeed generates spatial translations, let's calculate using the fact that :math:`\Pbf` is time-independent as follows

.. math::
	:nowrap:

	\begin{align*}
		\left[ \Pbf, Q_n(t, \xbf) \right] &= \ifrak \nabla Q_n(t, \xbf) \\
		\left[ \Pbf, P_n(t, \xbf) \right] &= \ifrak \nabla P_n(t, \xbf)
	\end{align*}

It follows that

.. math::
	:nowrap:

	\begin{equation}
		\left[ \Pbf, \Gscr \right] = \ifrak \nabla \Gscr
		\label{eq_momenta_act_as_spatial_derivative}
	\end{equation}

for any functional :math:`\Gscr` that doesn't explicitly involve :math:`\xbf`. This verifies that :math:`\Pbf` indeed generates the spatial translation.

In contrast, one cannot hope that the Lagrangian to be invariant under time translation, if there should be any interaction. But we already know the operator that generates time translation, namely, the Hamiltonian. In other words, we define :math:`P_0 \coloneqq -H` such that

.. math::
	:nowrap:

	\begin{equation}
		\left[ H, \Gscr \right] = -\ifrak \dot{\Gscr}
		\label{eq_hamiltonian_acts_as_time_derivative}
	\end{equation}

for any functional :math:`\Gscr` that doesn't explicitly involve :math:`t`.

In general, the Lagrangian density is not invariant under spacetime translations. However, it turns out that the conserved current, which in this case is :math:`T^{\mu}_{\nu}`, can nonetheless be calculated. To spell out the details, let's consider the following variation

.. math::
	:nowrap:

	\begin{equation*}
		\Psi_n(x) \to \Psi_n(x + \epsilon(x)) = \Psi_n(x) + \epsilon^{\mu}(x) \p_{\mu} \Psi_n(x)
	\end{equation*}

The corresponding variation of the action is given as follows

.. math::
	:nowrap:

	\begin{align*}
		\delta I[\Psi] &= \sum_n \int d^4 x \left( \frac{\delta \Lscr}{\delta \Psi_n} \epsilon^{\mu} \p_{\mu} \Psi_n + \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\nu}(\epsilon^{\mu} \p_{\mu} \Psi_n) \right) \\
			&= \int d^4 x \left( \epsilon^{\mu} \p_{\mu} \Lscr + \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu}\Psi_n \p_{\nu} \epsilon^{\mu} \right) \\
			&= -\int d^4 x \left( \delta^{\nu}_{\mu} \Lscr - \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu} \Psi_n \right) \p_{\nu} \epsilon^{\mu}
	\end{align*}

where we've used the chain rule for derivatives in the second equality, and integration by parts in the third. Comparing with :math:`\eqref{eq_variation_of_action_by_functional_deformation}`, we see that

.. math::
	:nowrap:

	\begin{equation*}
		T^{\nu}_{\mu} = \delta^{\nu}_{\mu} \Lscr - \sum_n \frac{\delta \Lscr}{\delta (\p_{\nu} \Psi_n)} \p_{\mu} \Psi_n
	\end{equation*}

.. _note_energy_momentum_tensor_not_symmetric:

.. note::
	The energy-momentum tensor :math:`T_{\mu}^{\nu}` is not yet suitable for general relativity since it's not symmetric. As we'll see in :ref:`sec_lorentz_symmetry`, when taking homogeneous Lorentz transformation symmetry into account, one can supplement :math:`T_{\mu}^{\nu}` with some extra terms to make it both conserved and symmetric.

Indeed, this calculation recovers :math:`\eqref{eq_spatial_translation_conserved_quantity}` by letting :math:`\nu = 0` and :math:`\mu \neq 0`. Moreover, it recovers the Hamiltonian by letting :math:`\mu = \nu = 0` as follows

.. math::
	:nowrap:

	\begin{equation*}
		H = -P_0 = \int d^3 x \left( \sum_n P_n \dot{Q}_n - \Lscr \right)
	\end{equation*}


Linear transformations
++++++++++++++++++++++

As another example, let's consider linear variations as follows

.. math::
	:nowrap:

	\begin{align*}
		Q_n(x) &\to Q_n(x) + \ifrak \epsilon^a (t_a)_n^m Q_m(x) \\
		C_r(x) &\to C_r(x) + \ifrak \epsilon^a (\tau_a)_r^s C_s(x)
	\end{align*}

where we've adopted the Einstein summation convention for repeated upper and lower indexes because it'd otherwise be too tedious to write out the summations. Here :math:`(t_{\square})^{\square}_{\square}` should furnish a representation of the Lie algebra of the symmetry group.

As before, the invariance of action under such variations implies the existstence of conserved currents :math:`J_a^{\mu}` such that

.. math::
	:nowrap:

	\begin{equation*}
		\p_{\mu} J_a^{\mu} = 0
	\end{equation*}

as well as the conserved quantity

.. math::
	:nowrap:

	\begin{equation*}
		T_a \coloneqq \int d^3 x~J^0_a
	\end{equation*}

If, in addition, the Lagrangian is invariant under such variations, then :math:`T_a` takes the following form by :math:`\eqref{eq_lagrangian_preserving_symmetry_generator_formula}`

.. math::
	:nowrap:

	\begin{equation*}
		T_a = -\ifrak \int d^3 x~P_n(t, \xbf) (t_a)^n_m Q^m(t, \xbf)
	\end{equation*}

It follows that

.. math::
	:nowrap:

	\begin{align*}
		\left[ T_a, Q^n(x) \right] &= -(t_a)^n_m Q^m(x) \\
		\left[ T_a, P_n(x) \right] &= (t_a)^m_n P_m(x)
	\end{align*}

In particular, when :math:`t_a` is diagonal (e.g., in electrodynamics), the operators :math:`Q^n` and :math:`P_n` may be regarded as raising/lowering operators. In fact, we claim that :math:`T_a` form a Lie algebra by the following calculation

.. math::
	:nowrap:

	\begin{align*}
		\left[ T_a, T_b \right] &= -\left[ \int d^3 x~P_n(t, \xbf) (t_a)^n_m Q^m(t, \xbf), \int d^3 y~P_r(t, \ybf) (t_b)^r_s Q^s(t, \ybf) \right] \\
			&= -\int d^3 x~d^3 y~(t_a)^n_m (t_b)^r_s \left[ P_n(t, \xbf) Q^m(t, \xbf), P_r(t, \ybf) Q^s(t, \ybf) \right] \\
			&= -\int d^3 x~d^3 y~(t_a)^n_m (t_b)^r_s \Big( P_n(t, \xbf) \left[ Q^m(t, \xbf), P_r(t, \ybf) \right] Q^s(t, \ybf) - P_r(t, \ybf) \left[ Q^s(t, \ybf), P_n(t, \xbf) \right] Q^m(t, \xbf) \Big) \\
			&= -\ifrak \int d^3 x \Big( (t_a)^n_m (t_b)^m_s P_n(t, \xbf) Q^s(t, \xbf) - (t_a)^n_m (t_b)^r_n P_r(t, \xbf) Q^m(t, \xbf) \Big) \\
			&= -\ifrak \int d^3 x~\left[ t_a, t_b \right]^n_m P_n(t, \xbf) Q^m(t, \xbf)
	\end{align*}

Now if :math:`t_a` form a Lie algebra with structure constants :math:`f^c_{ab}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\left[ t_a, t_b \right] = \ifrak f^c_{ab} t_c
	\end{equation*}

then

.. math::
	:nowrap:

	\begin{equation*}
		\left[ T_a, T_b \right] = \ifrak f^c_{ab} T_c
	\end{equation*}

In other words, the conserved quantities also form the same Lie algebra.

Now if, in addition, the Lagrangian density is also invariant, then :math:`\eqref{eq_lagrangian_density_preserving_symmetry_conserved_density}` takes the following form

.. math::
	:nowrap:

	\begin{equation}
		J^{\mu}_a = -\ifrak \frac{\delta \Lscr}{\delta (\p_{\mu} Q_n)} (t_a)^m_n Q_m - \ifrak \frac{\delta \Lscr}{\delta (\p_{\mu} C_r)} (\tau_a)^s_r C_s
		\label{eq_lagrangian_density_invariant_linear_transformation_conserved_current}
	\end{equation}

.. dropdown:: Interacting equal-mass real scalar fields
	:animate: fade-in-slide-down
	:icon: unlock

	As a specific example, let's consider the following Lagrangian density for two interacting equal-mass real scalar fields

	.. math::
		:nowrap:

		\begin{equation}
			\Lscr = -\frac{1}{2} \p_{\mu} \Phi_1 \p^{\mu} \Phi_1 - \frac{1}{2} m^2 \Phi_1^2 - \frac{1}{2} \p_{\mu} \Phi_2 \p^{\mu} \Phi_2 - \frac{1}{2} m^2 \Phi_2^2 - \Hscr(\Phi_1^2 + \Phi_2^2)
		\end{equation}

	This density is invariant under the following linear transformation

	.. math::
		:nowrap:

		\begin{align*}
			\Phi_1 &\to \Phi_1 - \epsilon \Phi_2 \\
			\Phi_2 &\to \Phi_2 + \epsilon \Phi_1
		\end{align*}

	In this case, we can evaluate :math:`\eqref{eq_lagrangian_density_invariant_linear_transformation_conserved_current}` as follows

	.. math::
		:nowrap:

		\begin{equation*}
			J^{\mu} = -\Psi_2 \p^{\mu} \Phi_1 + \Phi_1 \p^{\mu} \Phi_2
		\end{equation*}

Note that since :math:`\Lscr` doesn't have :math:`\dot{C}_r` dependencies, we have the following by letting :math:`\mu = 0` in :math:`\eqref{eq_lagrangian_density_invariant_linear_transformation_conserved_current}`

.. math::
	:nowrap:

	\begin{equation*}
		J^0_a = -\ifrak P^n (t_a)^m_n Q_m
	\end{equation*}

whose equal-time commutation relations with canonical variables :math:`P` and :math:`Q` can be easily calculated.


.. _sec_lorentz_invariance:

Lorentz invariance
^^^^^^^^^^^^^^^^^^

The goal of this section is to show that the Lorentz invariance of the Lagrangian density implies the Lorentz invariance of the S-matrix, which justifies our interest in the Lagrangian formalism in the first place.

Recall from :math:`\eqref{eq_expansion_of_Lambda}` and :math:`\eqref{eq_lorentz_lie_algebra_is_antisymmetric}` that

.. math::
	:nowrap:

	\begin{align}
		\Lambda^{\mu}_{\nu} &= \delta^{\mu}_{\nu} + \omega^{\mu}_{\nu}
		\nonumber \\
		\omega_{\mu \nu} &= -\omega_{\nu \mu}
		\label{eq_lorentz_omega_is_antisymmetric}
	\end{align}

is a :math:`(\mu, \nu)`-parametrized anti-symmetric variation. It follows then from :math:`\eqref{eq_variation_of_action_by_functional_deformation}` that there exist :math:`(\mu, \nu)`-parametrized anti-symmetric conserved currents as follows

.. math::
	:nowrap:

	\begin{align}
		\p_{\rho} \Mscr^{\rho \mu \nu} &= 0
		\label{eq_lorentz_invariance_m_conservation} \\
		\Mscr^{\rho \mu \nu} &= -\Mscr^{\rho \nu \mu}
		\label{eq_lorentz_invariance_conserved_m_antisymmetric}
	\end{align}

which, in turn, make conversed quantities

.. math::
	:nowrap:

	\begin{equation}
		J^{\mu \nu} \coloneqq \int d^3 x~\Mscr^{0 \mu \nu}
		\label{eq_lorentz_invariance_conserved_j}
	\end{equation}

such that :math:`\dot{J}^{\mu \nu} = 0`. These, as we'll see, turn out to be rather familiar objects that we've encountered as early as in :math:`\eqref{eq_u_lorentz_expansion}`.

In light of :math:`\eqref{eq_lagrangian_density_preserving_symmetry_conserved_density}`, one can work out an explicit formula for :math:`\Mscr^{\rho \mu \nu}` if the Lagrangian density is invariant under the symmetry transformation. Now since the Lagrangian density is expressed in terms of quantum fields, one'd like to know how they transform under Lorentz transformations. Since the translation symmetry has already been dealt with in :ref:`sec_spacetime_translations`, we'll consider here homogeneous Lorentz transformations. Luckily this has been worked out already in :ref:`sec_quantum_fields_and_antiparticles`. More precisely, recall from :math:`\eqref{eq_dirac_field_linearize_representation}` that the variation term can be written as follows

.. math::
	:nowrap:

	\begin{equation*}
		\delta \Psi_n = \frac{\ifrak}{2} \omega^{\mu \nu} (\Jscr_{\mu \nu})_n^m \Psi_m
	\end{equation*}

where :math:`\Jscr` are matrices satisfying :math:`\eqref{eq_bracket_repr_j}`. The corresponding derivatives then have the following variation term

.. math::
	:nowrap:

	\begin{equation*}
		\delta (\p_{\kappa} \Psi_n) = \frac{\ifrak}{2} \omega^{\mu \nu} (\Jscr_{\mu \nu})_n^m \p_{\kappa} \Psi_m + \omega_{\kappa}^{\lambda} \p_{\lambda} \Psi_n
	\end{equation*}

where the second summand on the right-hand-side corresponds to the fact the the Lorentz transformation also acts on the spacetime coordinates.

Now the invariance of the Lagrangian density under such variation can be written as follows

.. math::
	:nowrap:

	\begin{equation*}
		\frac{\delta \Lscr}{\delta \Psi_n} \frac{\ifrak}{2} \omega^{\mu \nu} (\Jscr_{\mu \nu})_n^m \Psi_m
      		+ \frac{\delta \Lscr}{\delta (\p_{\kappa} \Psi_n)} \left( \frac{\ifrak}{2} \omega^{\mu \nu} (\Jscr_{\mu \nu})_n^m \p_{\kappa} \Psi_m + \omega_{\kappa}^{\lambda} \p_{\lambda} \Psi_n \right) = 0
	\end{equation*}

Since :math:`\omega^{\mu \nu}` is not in general zero, its coefficient must be zero, which, taking :math:`\eqref{eq_lorentz_omega_is_antisymmetric}` into account, implies the following

.. math::
	:nowrap:

	\begin{equation}
		\frac{\ifrak}{2} \frac{\delta\Lscr}{\delta\Psi_n} (\Jscr_{\mu\nu})^m_n \Psi_m + \frac{\ifrak}{2} \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} (\Jscr_{\mu\nu})^m_n \p_{\kappa}\Psi_m + \frac{1}{2} \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} \left( \eta_{\kappa \mu} \p_{\nu} - \eta_{\kappa \nu} \p_{\mu} \right) \Psi_n = 0
		\label{eq_lorentz_invariance_current_raw_identity}
	\end{equation}

Using :math:`\eqref{eq_euler_lagrange}`, we can get rid of the :math:`\delta\Lscr / \delta\Psi_n` term in :math:`\eqref{eq_lorentz_invariance_current_raw_identity}` to arrive at the following

.. math::
	:nowrap:

	\begin{equation}
		\ifrak \p_{\kappa} \left( \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} (\Jscr_{\mu\nu})_n^m \Psi_m \right) - \frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} (T_{\mu\nu} - T_{\nu\mu}) = 0
		\label{eq_lorentz_invariance_current_identity}
	\end{equation}

Now we can address the issue of :ref:`energy-momentum tensor not being symmetric <note_energy_momentum_tensor_not_symmetric>` by introducing the following so-called `Belinfante tensor <https://en.wikipedia.org/wiki/Belinfante%E2%80%93Rosenfeld_stress%E2%80%93energy_tensor>`__

.. math::
	:nowrap:

	\begin{equation}
		\Theta_{\mu\nu} \coloneqq T_{\mu\nu} - \frac{\ifrak}{2} \p_{\kappa} \left(
			\frac{\delta\Lscr}{\delta(\p_{\kappa} \Psi_n)} (\Jscr_{\mu\nu})_n^m \Psi_m -
			\frac{\delta\Lscr}{\delta(\p_{\mu} \Psi_n)} (\Jscr_{\kappa\nu})_n^m \Psi_m -
			\frac{\delta\Lscr}{\delta(\p_{\nu} \Psi_n)} (\Jscr_{\kappa\mu})_n^m \Psi_m \right)
		\label{eq_defn_belinfante_tensor}
	\end{equation}

which is both conserved in the sense that

.. math::
	:nowrap:

	\begin{equation}
		\p^{\mu} \Theta_{\mu\nu} = 0
		\label{eq_belinfante_tensor_is_conserved}
	\end{equation}

and symmetric in the sense that

.. math::
	:nowrap:

	\begin{equation}
		\Theta_{\mu\nu} = \Theta_{\nu\mu}
		\label{eq_belinfante_tensor_is_symmetric}
	\end{equation}

Indeed :math:`\eqref{eq_belinfante_tensor_is_conserved}` follows from the observation that the term inside the parenthesis of :math:`\eqref{eq_defn_belinfante_tensor}` is anti-symmetric in :math:`\mu` and :math:`\kappa`, and :math:`\eqref{eq_belinfante_tensor_is_symmetric}` is a direct consequence of :math:`\eqref{eq_lorentz_invariance_current_identity}`.

The conserved quantities corresponding to :math:`\Theta_{\mu\nu}` are

.. math::
	:nowrap:

	\begin{equation}
		\int d^3 x~\Theta_{0 \nu} = \int d^3 x~T_{0 \nu} \xlongequal{\eqref{eq_spacetime_translation_conserved_quantity_is_momentum}} P_{\nu}
		\label{eq_p_as_integral_of_belinfante_tensor}
	\end{equation}

where the first first equality holds because, again, the item in the parenthesis of :math:`\eqref{eq_defn_belinfante_tensor}` is anti-symmetric is :math:`\mu` and :math:`\nu`, and therefore :math:`\nu \neq 0` given :math:`\mu = 0`. Hence it's at least equally legitimate to call :math:`\Theta_{\mu \nu}` the energy-momentum tensor. Indeed, the fact that :math:`\Theta_{\mu \nu}` is the symmetric makes it the right choice in general relatively.

Unlike the other conserved currents, which are derived under the general principles explained in :ref:`sec_from_symmetries_to_conservation_laws`, we'll construct the anti-symmetric :math:`\Mscr^{\rho \mu \nu}` declared in :math:`\eqref{eq_lorentz_invariance_m_conservation}` and :math:`\eqref{eq_lorentz_invariance_conserved_m_antisymmetric}` by hand as follows

.. math::
	:nowrap:

	\begin{equation*}
		\Mscr^{\rho\mu\nu} \coloneqq x^{\mu} \Theta^{\rho\nu} - x^{\nu} \Theta^{\rho\mu}
	\end{equation*}

While :math:`\eqref{eq_lorentz_omega_is_antisymmetric}` is automatically satisfied by definition, we can verify :math:`\eqref{eq_lorentz_invariance_m_conservation}` as follows

.. math::
	:nowrap:

	\begin{equation*}
		\p_{\rho} \Mscr^{\rho\mu\nu} \xlongequal{\eqref{eq_belinfante_tensor_is_conserved}} \Theta^{\mu\nu} - \Theta^{\nu\mu} \xlongequal{\eqref{eq_belinfante_tensor_is_symmetric}} 0
	\end{equation*}

Moreover :math:`\eqref{eq_lorentz_invariance_conserved_j}` takes the following form

.. math::
	:nowrap:

	\begin{equation*}
		J^{\mu\nu} = \int d^3 x \left( x^{\mu} \Theta^{0\nu} - x^{\nu} \Theta^{0\mu} \right)
	\end{equation*}

Now if we consider the rotation generators :math:`J_i \coloneqq \tfrac{1}{2} \epsilon_{ijk} J^{jk}`, then it follows from :math:`\eqref{eq_hamiltonian_acts_as_time_derivative}` that

.. math::
	:nowrap:

	\begin{equation*}
		[H, \Jbf] = -\ifrak \dot{\Jbf} = 0
	\end{equation*}

since :math:`\Jbf` doesn't implicitly involve :math:`t`. This recovers one of the commutation relations :math:`\eqref{eq_hj_commute}` for the Poincaré algebra. Next, let's verify :math:`\eqref{eq_pjp_commutation}` as follows

.. math::
	:nowrap:

	\begin{align*}
		[P_i, J_j] &\xlongequal{\phantom{\eqref{eq_momenta_act_as_spatial_derivative}}} \frac{1}{2} \epsilon_{jk\ell} \left[ P_i, J^{k\ell} \right] \\
			&\xlongequal{\eqref{eq_momenta_act_as_spatial_derivative}} \frac{\ifrak}{2} \epsilon_{jk\ell} \int d^3x \left( x^k \p_i \Theta^{0\ell} - x^{\ell} \p_i \Theta^{0k} \right) \\
			&\xlongequal{\phantom{\eqref{eq_momenta_act_as_spatial_derivative}}} \frac{\ifrak}{2} \epsilon_{jk\ell} \int d^3x \left( -\delta^k_i \Theta^{0\ell} + \delta^{\ell}_i \Theta^{0k} \right) \\
			&\xlongequal{\phantom{\eqref{eq_momenta_act_as_spatial_derivative}}} \ifrak \epsilon_{ijk} \int d^3x~\Theta^{0k} \\
			&\xlongequal{\eqref{eq_p_as_integral_of_belinfante_tensor}} \ifrak \epsilon_{ijk} P^k
	\end{align*}

What come next are the boost operators defined as follows

.. math::
	:nowrap:

	\begin{equation*}
		K_i = K^i \coloneqq J^{0i} = \int d^3x \left( x^0 \Theta^{0i} - x^i \Theta^{00} \right)
	\end{equation*}



.. rubric:: Footnotes

.. [#tedious_calc_of_commutations] These are some rather tedious and error-prone calculations, but in the end, we manage to arrive at the same results as stated in [Wei95]_ page 61.

.. [#in_out_state_sign_convention] I really don't like the sign convention of the in- and out-states in [Wei95]_, even though Weinberg said that it has become a tradition. So our sign convention, as far as the definitions of in- and out-states are concerned, is the opposite to the one used in [Wei95]_.

.. [#pion_deuteron_reaction_final_state] The final state claimed in [Wei95]_ page 126 has total spin :math:`1` rather than :math:`0`. I suspect that the claim in [Wei95]_ is wrong, but otherwise it doesn't affect any subsequent arguments anyway.

.. [#abuse_of_phi_as_both_state_vector_and_flux] It's really unfortunate that one constantly runs out symbols to represent physical quantities, and it's uncommon to use anything other than just one letter. So we have to live with the fact that the meaning of the symbol will depend on the context.

.. [#weinberg_quote_on_lorentz_invariance_and_quantum_mechanics] See [Wei95]_ page 145.

.. [#two_ways_of_representation] This should be compared with :math:`\eqref{eq_d_repr_of_little_group}`, which uses transpose instead of inverse to satisfy the group law. From the viewpoint of representation theory, it's more natural to use the inverse. But in the case of little group representations, since we're only interested in unitary representations, there is not much difference between the two choices.

.. [#wrong_integration_of_Delta_function] The evaluation of the integral eq. (5.2.8) on [Wei95]_ page 202 seems to be wrong, as the integrand oscillates as :math:`\sin(u)` for :math:`u` large enough, which will cause the integral to diverge.

.. [#dirac_gamma_is_vector] Our calculation is consistent with eq. (5.4.8) in [Wei95]_ because of the convention :math:`{(\Lambda^{-1})^{\nu}}_{\mu} = {\Lambda_{\mu}}^{\nu}` introduced in formula (2.3.10).

.. [#charge_inversion_on_dirac_fields_sign] Our calculation differs from the calculation eq. (5.5.47) in [Wei95]_ by a sign.

.. [#clebsch_gordan_coefficients_orthonormality] Details of the argument can be found in [Wei15]_ page 121 -- 122.

.. [#clebsch_gordan_coefficient_zero_total_angular_momentum] An explicit evaluation of :math:`C^{AA}(0, 0;a, -a)`, with the only non-vanishing combination of spin :math:`z`-components, can be found in [Wei15]_ page 124 -- 125.

.. [#volume_element_inversion_sign_cancels_for_cpt_symmetry] Comparing with eq. (5.8.9) in [Wei95]_ page 246, we get probably the most important result right, even though we've insisted in reversing the sign of the volume elements such as :math:`d^3 p` and :math:`d^3 x`, when :math:`p` and :math:`x` are reversed, respectively.

.. [#connected_s_matrix_by_feynman_diagrams] It's stated in [Wei95]_ page 271 and afterwards that it's the full S-matrices that are calculated. However, it seems that only the connected part of the S-matrices are calculated.