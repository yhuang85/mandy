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

Now let :math:`\eta` be another element of the Lie algebra, and expand both sides of :math:`U(T(\eta)) U(T(\theta)) = U(T(\eta) T(\theta))` as follows

.. math::
	:nowrap:

	\begin{align*}
	U(T(\eta)) U(T(\theta)) &= \left( 1 + \ifrak \eta^a u_a + \tfrac{1}{2} \eta^a \eta^b u_{ab} + \cdots \right) \left( 1 + \ifrak \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots \right) \\
		&= 1 + \ifrak (\eta^a + \theta^a) u_a \blue{- \eta^a \theta^b u_a u_b} + \cdots \\
	\\
	U(T(\eta) T(\theta)) &= U \left( 1 + \eta + \theta + f_{ab} \eta^a \theta^b + \cdots \right) \\
		&= 1 + \blue{\ifrak} \left( \eta^c + \theta^c + \blue{f^c_{ab} \eta^a \theta^b} + \cdots \right) \blue{u_c} + \blue{\tfrac{1}{2}} \left( \blue{\eta^a + \theta^a} + \cdots \right) \left( \blue{\eta^b + \theta^b} + \cdots \right) \blue{u_{ab}} + \cdots
	\end{align*}

where :math:`f^c_{ab}` are the coefficients of the expansion of :math:`T(f(\eta, \theta)) = T(\eta) T(\theta)`. Equating the coefficients of :math:`\eta^a \theta^b`, i.e., the terms colored in blue, we get

.. math::
	:nowrap:

	\begin{equation*}
		-u_a u_b = \ifrak f^c_{ab} u_c + u_{ab} \quad \Longrightarrow \quad u_{ab} = -u_a u_b - \ifrak f^c_{ab} u_c.
	\end{equation*}

It implies that one can calculate the higher-order operator :math:`u_{ab}` from the lower-order ones, assuming of course that we know the structure of the symmetry (Lie) group/algebra. In fact, this bootstrapping procedure can be continued to all orders, but we'll not be bothered about the details.

Next, note that :math:`u_{ab} = u_{ba}` since they are just partial derivatives. It follows that

.. math::
	:nowrap:

	\begin{equation*}
		[u_a, u_b] \coloneqq u_a u_b - u_b u_a = \ifrak (f^c_{ba} - f^c_{ab}) u_c \eqqcolon \ifrak \Gamma^c_{ab} u_c
	\end{equation*}

where the bracket is known as the *Lie bracket* and :math:`\Gamma^c_{ab}` are known as the *structure constants*.

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
	1. We will follow the common convention in physics that greek letters such as :math:`\mu, \nu, \dots` run from :math:`0` to :math:`3`, while roman letters such as :math:`i, j, \dots` run from :math:`1` to :math:`3`.
	2. We often write :math:`x` for a spacetime point :math:`(x_0, x_1, x_2, x_3)`, and :math:`\xbf` for a spatial point :math:`(x_1, x_2, x_3)`.
	3. Matrix/Tensor indexes can will be raised or lowered using :math:`\eta`. For example :math:`dx^0 = \eta^{0 \mu} dx_{\mu} = -dx_0`.
	4. Einstein's summation convention is implicitly applied, i.e., upper and lower indexes of the same label are summed up. *The same label that appears both as upper or lower indexes are considered illegal expression.* This is particularly important because otherwise we'd run into troublesome expressions such as :math:`x_{\mu} x_{\mu} = x_0^2 + x_1^2 + x_2^2 + x_3^2`.

.. dropdown:: Einstein's special theory of relativity
	:animate: fade-in-slide-down

	Using the notations introduced above, we can rewrite :math:`\eqref{eq_proper_time}` as :math:`d\tau^2 = dt^2 - d\xbf^2`, so that it's obvious that if a particle travels at the speed of light in one inertial frame, i.e., :math:`|d\xbf / dt| = 1`, and equivalently :math:`d\tau = 0`, then it travels at the speed of light in any other inertial frame, in direct contradiction with Newtonian mechanics.

	Instead of working with the spacetime coordinates, it can sometimes be convenient to work with the "dual" energy-momentum coordinates, also known as the *four momentum*. The transition can be done by imagining a particle of mass :math:`m`, and defining :math:`p = (E, \pbf) \coloneqq m dx / d\tau`. It follows from :math:`\eqref{eq_proper_time}` that

	.. math::
		:nowrap:

		\begin{equation*}
			1 = (dt / d\tau)^2 - (d\xbf / d\tau)^2 ~\Longrightarrow~ m^2 = (m dt / d\tau)^2 - (m d\xbf / d\tau)^2 = E^2 - \pbf^2
		\end{equation*}

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
		\eta^{\mu \nu} dx_{\mu} dx_{\nu} = \eta^{\mu \nu} \Lambda^{\rho}_{\mu} \Lambda^{\kappa}_{\nu} dx_{\rho} dx_{\kappa} \
			~\Longrightarrow~ \eta^{\mu \nu} = \eta^{\rho \kappa} \Lambda^{\mu}_{\rho} \Lambda^{\nu}_{\kappa}
		\label{eq_homogeneous_lorentz_transformation}
	\end{equation}

for any :math:`\mu, \nu`. Moreover the group law is given by

.. math::
	:nowrap:

	\begin{equation*}
		L(\Lambda', a') L(\Lambda, a) x = L(\Lambda', a')(\Lambda x + a) = \Lambda' \Lambda x + \Lambda' a + a' = L(\Lambda' \Lambda, \Lambda' a + a') x
	\end{equation*}

For later use, let's also calculate the inverse transformation using :math:`\eqref{eq_homogeneous_lorentz_transformation}`

.. math::
	:nowrap:

	\begin{equation}
		\delta_{\sigma}^{\nu} = \eta_{\sigma \mu} \eta^{\mu \nu} = \left(\eta_{\sigma \mu} \eta^{\rho \kappa} \Lambda_{\rho}^{\mu}\right) \Lambda^{\nu}_{\kappa} \
			~\Longrightarrow~ (\Lambda^{-1})_{\mu}^{\nu} = \eta_{\mu \sigma} \eta^{\nu \rho} \Lambda_{\rho}^{\sigma} \
			~\Longleftrightarrow~ (\Lambda^{-1})^{\mu \nu} = \Lambda^{\nu \mu}
		\label{eq_lambda_inverse}
	\end{equation}

Taking determinant on :math:`\eqref{eq_homogeneous_lorentz_transformation}` implies that :math:`\op{det}(\Lambda) = \pm 1`, and setting :math:`\mu = \nu = 0` implies that

.. math::
	:nowrap:

	\begin{equation*}
		1 = (\Lambda^0_0)^2 - \sum_{i=1}^3 \Lambda^0_i \Lambda^0_i ~\Longrightarrow~ \left| \Lambda^0_0 \right| \geq 1
	\end{equation*}

It follows that the homogeneous Lorentz group has four components. In particular, the one with :math:`\op{det}(\Lambda) = 1` and :math:`\Lambda^0_0 \geq 1` is the most common used and is given a name: *proper orthochronous* Lorentz group. Nonetheless, one can map one component to another by composing with either a time reversal transformation

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
		dt' = \Lambda^0_0 dt, \quad dx'_i = \Lambda^0_i dt ~\Longrightarrow~ \Lambda^0_i = v_i \Lambda^0_0
	\end{equation*}

Then using :math:`\eqref{eq_homogeneous_lorentz_transformation}`, we get

.. math::
	:nowrap:

	\begin{equation}
		1 = -\eta^{\mu \nu} \Lambda^0_{\mu} \Lambda^0_{\nu} = (\Lambda^0_0)^2 - \Lambda^0_i \Lambda^0_i = (1 - \vbf^2) (\Lambda^0_0)^2 \
			~\Longrightarrow~ \Lambda^0_0 = \frac{1}{\sqrt{1 - \vbf^2}} \eqqcolon \gamma
			\label{eq_def_gamma}
	\end{equation}

assuming :math:`\Lambda` is proper orthochronous. It follows that

.. math::
	:nowrap:

	\begin{equation}
		\Lambda^0_0 = \gamma, \quad \Lambda^0_i = \gamma v_i
		\label{eq_lambda_boost}
	\end{equation}

The other components :math:`\Lambda^i_j, 1 \leq i, j \leq 3`, are not uniquely determined because a composition with a (spatial) rotation about the direction of :math:`\vbf` has no effect on :math:`\vbf`. To make it easier, one can apply a rotation so that :math:`\vbf` aligns with the :math:`3`-axis. Then an obvious choice of :math:`\Lambda` is given by

.. math::
	:nowrap:

	\begin{alignat}{2}
		t' &=  \Lambda_0^{\mu} x_{\mu} &&= \gamma (t + v_3 x_3)  \label{eq_lambda_in_3_axis} \\
		x'_1 &= \Lambda_1^{\mu} x_{\mu} &&= x_1  \nonumber \\
		x'_2 &= \Lambda_2^{\mu} x_{\mu} &&= x_2  \nonumber \\
		x'_3 &= \Lambda_3^{\mu} x_{\mu} &&= \gamma (x_3 + v_3 t) \nonumber
	\end{alignat}

.. dropdown:: Time dilation and length contraction
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
			|(\Lambda x)_3 - (\Lambda y)_3| = \gamma |x_3 - y_3 - v_3^2 (x_3 - y_3)| = \frac{|x_3 - y_3|}{\gamma} < |x_3 - y_3|
		\end{equation*}

	This calculation says that the length of rod is contracted in the direction of movement.

	It should be emphasized that such contraction of length can only be observed in a frame where the rod is moving. Imagine for example a scenario where you're given a square box with equal sides while standing still, then after some unconscious period of time, e.g., sleeping, you wake up with the same box in hand, and you'd like to know if you're now moving. If you happen to have heard of such contraction of length, you might try to measure the sides of the box again. If one of the sides suddenly becomes shorter, then you know not only that you're moving, but also the direction of movement! This is of course absurd because the box is still at rest relative to you.

Finally, one can apply a rotation to :math:`\eqref{eq_lambda_in_3_axis}` to get the general formulae

.. math::
	:nowrap:

	\begin{equation}
		\Lambda_{ij} = \delta_{ij} + \frac{v_i v_j}{\vbf^2} (\gamma - 1)
		\label{eq_general_lambda_in_spacetime}
	\end{equation}

for :math:`1 \leq i, j \leq 3`, which, together with :math:`\eqref{eq_lambda_boost}` and :math:`\Lambda_0^i = \Lambda_i^0,` gives the general formula for :math:`\Lambda`.

.. note::
	Any Lorentz transformation can be written as the composition of a boost followed by a rotation.

.. _sec_quantum_lorentz_symmetry:

Quantum Lorentz symmetry
++++++++++++++++++++++++

We will quantize the Lorentz symmetry :math:`L(\Lambda, a)` by looking for unitarity representations :math:`U(\Lambda, a)`. As discussed in :ref:`sec_continuous_symmetry`, we proceed by looking for infinitesimal symmetries. First of all, let's expand :math:`\Lambda` as

.. math::
	:nowrap:

	\begin{equation*}
		\Lambda_{\mu}^{\nu} = \delta_{\mu}^{\nu} + \omega_{\mu}^{\nu} + \cdots
	\end{equation*}

where :math:`\delta` is the Kronecker delta, and *not* a tensor. It follows from :math:`\eta^{\mu \nu} = \eta^{\rho \kappa} \Lambda^{\mu}_{\rho} \Lambda^{\nu}_{\kappa}` that

.. math::
	:nowrap:

	\begin{align*}
		\eta^{\mu \nu} &= \eta^{\rho \kappa} (\delta_{\rho}^{\mu} + \omega_{\rho}^{\mu} + \cdots) (\delta_{\kappa}^{\nu} + \omega_{\kappa}^{\nu} + \cdots) \\
			&= \eta^{\mu \nu} + \eta^{\mu \kappa} \omega^{\nu}_{\kappa} + \eta^{\rho \nu} \omega^{\mu}_{\rho} + \cdots \\
			&= \eta^{\mu \nu} + \omega^{\mu \nu} + \omega^{\nu \mu} + \cdots
	\end{align*}

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
			&= 1 - \ifrak (\Lambda^{\rho}_{\mu} \epsilon^{\mu} - (\Lambda \omega \Lambda^{-1})^{\rho \kappa} a_{\kappa}) P_{\rho} \
				+ \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} J_{\rho \kappa} + \cdots \\
			&= 1 -\ifrak \epsilon^{\mu} \Lambda_{\mu}^{\rho} P_{\rho} + \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots \\
			&= 1 -\ifrak \epsilon^{\mu} \Lambda_{\mu}^{\rho} P_{\rho} + \tfrac{\ifrak}{2} \Lambda^{\rho \mu} \omega_{\mu \nu} \Lambda^{\kappa \nu} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots \\
			&= 1 -\ifrak \epsilon^{\mu} \Lambda_{\mu}^{\rho} P_{\rho} + \tfrac{\ifrak}{2} \omega^{\mu \nu} (\Lambda^{-1})^{\rho}_{\mu} (\Lambda^{-1})^{\kappa}_{\nu} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots
	\end{align*}


where we have used :math:`\eqref{eq_lambda_inverse}` for :math:`\Lambda^{-1}`. Substituting :math:`U(1 + \omega, \epsilon)` with the expansion :math:`\eqref{eq_u_lorentz_expansion}` and equating the coefficients of :math:`\epsilon^{\mu}` and :math:`\omega_{\mu \nu}`, we have

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda, a) P_{\mu} U^{-1}(\Lambda, a) &= \Lambda^{\rho}_{\mu} P_{\rho}  \label{eq_p_conjugated_by_u} \\
		U(\Lambda, a) J_{\mu \nu} U^{-1}(\Lambda, a) &= (\Lambda^{-1})^{\rho}_{\mu} (\Lambda^{-1})^{\kappa}_{\nu} (J_{\rho \kappa} + a_{\kappa} P_{\rho} - a_{\rho} P_{\kappa})  \label{eq_j_conjugated_by_u}
	\end{align}

where in the second equation, we have also made the right-hand-side anti-symmetric with respect to :math:`\mu` and :math:`\nu`. It's now clear that :math:`P` transforms like a vector and is translation invariant, while :math:`J` transforms like a :math:`2`-tensor only for homogeneous Lorentz transformations and is not translation invariant in general. These are of course as expected since both :math:`P` and :math:`J` are quantization of rather familiar objects, which we now spell out.

We start with :math:`P` by writing :math:`H \coloneqq P_0` and :math:`\Pbf \coloneqq (P_1, P_2, P_3)`. Then :math:`H` is the energy operator, also know as the *Hamiltonian*, and :math:`\Pbf` is the momentum :math:`3`-vector. Similarly, let's write :math:`\Kbf \coloneqq (J_{01}, J_{02}, J_{03})` and :math:`\Jbf = (J_{23}, J_{31}, J_{12})`, as the *boost* :math:`3`-vector and the *angular momentum* :math:`3`-vector, respectively.

Now that we have named all the players (i.e., :math:`H, \Pbf, \Jbf, \Kbf`) in the game, it remains to find out their mutual commutation relations since they should form a Lie algebra of the (infinitesimal) Lorentz symmetry. This can be done by applying :math:`\eqref{eq_p_conjugated_by_u}` and :math:`\eqref{eq_j_conjugated_by_u}` to :math:`U(\Lambda, a)` that is itself infinitesimal. More precisely, keeping up to first order terms, write :math:`\Lambda_{\mu}^{\nu} = \delta_{\mu}^{\nu} + \omega_{\mu}^{\nu}` and :math:`a^{\mu} = \epsilon^{\mu}` so that :math:`\eqref{eq_p_conjugated_by_u}` becomes

.. math::
	:nowrap:

	\begin{align*}
		\left( \delta_{\mu}^{\rho} + \omega_{\mu}^{\rho} \right) P_{\rho} &= \left( 1 - \ifrak \epsilon^{\nu} P_{\nu} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) P_{\mu} \left( 1 + \ifrak \epsilon^{\nu} P_{\nu} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
			&= P_{\mu} - \ifrak \epsilon^{\nu} [P_{\mu}, P_{\nu}] - \tfrac{\ifrak}{2} \omega^{\rho \kappa} [P_{\mu}, J_{\rho \kappa}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` and :math:`\omega` gives the following

.. math::
	:nowrap:

	\begin{align}
		[P_{\mu}, P_{\nu}] &= 0  \label{eq_bracket_p4_p4} \\
		[P_{\mu}, J_{\rho \kappa}] &= \ifrak (\eta_{\mu \rho} P_{\kappa} - \eta_{\mu \kappa} P_{\rho})  \label{eq_bracket_p4_j4}
	\end{align}

where we've used the identity :math:`\omega_{\mu}^{\rho} P_{\rho} = \eta_{\mu \kappa} \omega^{\kappa \rho} P_{\rho} = \tfrac{1}{2} \omega^{\rho \kappa} (\eta_{\mu \rho} P_{\kappa} - \eta_{\mu \kappa} P_{\rho})`. Next :math:`\eqref{eq_j_conjugated_by_u}` (up to first order) becomes

.. math::
	:nowrap:

	\begin{align*}
		J_{\mu \nu} + \epsilon_{\nu} P_{\mu} - \epsilon_{\mu} P_{\nu} - \omega_{\mu}^{\rho} J_{\rho \nu} - \omega_{\nu}^{\kappa} J_{\mu \kappa} \
			&= (\delta_{\mu}^{\rho} - \omega_{\mu}^{\rho}) (\delta_{\nu}^{\kappa} - \omega_{\nu}^{\kappa}) (J_{\rho \kappa} + \epsilon_{\kappa} P_{\rho} - \epsilon_{\rho} P_{\kappa}) \\
			&= \left( 1 - \ifrak \epsilon^{\rho} P_{\rho} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) J_{\mu \nu} \left( 1 + \ifrak \epsilon^{\rho} P_{\rho} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
			&= J_{\mu \nu} - \ifrak \epsilon^{\rho} [P_{\rho}, J_{\mu \nu}] + \tfrac{\ifrak}{2} \omega^{\rho \kappa} [J_{\rho \kappa}, J_{\mu \nu}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` reproduces :math:`\eqref{eq_bracket_p4_j4}`, but equating the coefficients of :math:`\omega` gives the following additional

.. math::
	:nowrap:

	\begin{equation}
		[J_{\rho \kappa}, J_{\mu \nu}] = \ifrak (\eta_{\rho \mu} J_{\kappa \nu} - \eta_{\kappa \mu} J_{\rho \nu} + \eta_{\nu \rho} J_{\mu \kappa} - \eta_{\nu \kappa} J_{\mu \rho})
		\label{eq_bracket_j4_j4}
	\end{equation}

Now that we have all the commutator relations, let's reorganize :math:`\eqref{eq_bracket_p4_p4}, \eqref{eq_bracket_p4_j4}, \eqref{eq_bracket_j4_j4}` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		\text{let } \mu = 0, \nu = i \text{ in \eqref{eq_bracket_p4_p4}} ~&\Longrightarrow~ [H, P_i] &&= 0  \label{eq_hp_commute} \\
		\text{let } \mu = 0, \rho = j, \kappa = k \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [H, J_i] &&= 0  \label{eq_hj_commute} \\
		\text{let } \mu = 0, \rho = 0, \kappa = i \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [H, K_i] &&= \ifrak P_i  \label{eq_hkp_commutation} \\
		\text{let } \mu = i, \nu = j \text{ in \eqref{eq_bracket_p4_p4}} ~&\Longrightarrow~ [P_i, P_j] &&= 0  \label{eq_pp_commute} \\
		\text{let } \mu = i, \rho = k, \kappa = i \text{ in \eqref{eq_bracket_p4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [P_i, J_j] &&= -\ifrak \epsilon_{ijk} P_k  \nonumber \\
		\text{let } \mu = i, \rho = 0 \text{ and enumerate } \kappa \in \{1, 2, 3\} \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [P_i, K_j] &&= \ifrak \delta_{ij} H  \label{eq_pkh_commutation} \\
		\text{let } \rho = j, \kappa = \mu = k, \nu = i \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [J_i, J_j] &&= \ifrak \epsilon_{ijk} J_k  \label{eq_jjj_commutation} \\
		\text{let } \rho = \nu = j, \kappa = k, \mu = 0 \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [J_i, K_j] &&= \ifrak \epsilon_{ijk} K_k  \label{eq_jkk_commutation} \\
		\text{let } \rho = \mu = 0, \kappa = i, \nu = j \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [K_i, K_j] &&= -\ifrak \epsilon_{ijk} J_k  \label{eq_kkj_commutation}
 	\end{alignat}

where :math:`\epsilon_{ijk}` is totally anti-symmetric with respect to permutations of indexes and satisfies :math:`\epsilon_{123} = 1`. [#tedious_calc_of_commutations]_

Since the time evolution of a physical system is dictated by the Hamiltonian :math:`H`, quantities (i.e., observables) that commute with :math:`H` are conserved. In particular :math:`\eqref{eq_hp_commute}` and :math:`\eqref{eq_hj_commute}` imply that both momentum and angular momentum are conserved. Boosts, on the other hand, are *not* conserved, and therefore cannot be used to label (stable) physical states. Moreover :math:`\eqref{eq_pp_commute}` implies that translations commute with each other (as expected), which is *not* the case for the angular momenta according to :math:`\eqref{eq_jjj_commutation}`. Indeed, they furnish an infinitesimal representation of the :math:`3`-rotation group :math:`SO(3)`.

One-Particle States
-------------------

One neat application of our knowledge about Lorentz symmetry is to classify (free) one-particle states according to their transformation laws under (inhomogeneous) Lorentz transformations. Throughout this section, the Lorentz transformations will be assumed to be proper orthochronous, i.e., :math:`\op{det}(\Lambda) = 1` and :math:`\Lambda_0^0 \geq 1`.

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
		P_{\mu} U(\Lambda) \Psi_{p, \sigma} = U(\Lambda) (U^{-1} (\Lambda) P_{\mu} U(\Lambda)) \Psi_{p, \sigma} \
		\xlongequal{\eqref{eq_p_conjugated_by_u}} U(\Lambda) \Lambda_{\mu}^{\nu} P_{\nu} \Psi_{p, \sigma} = (\Lambda_{\mu}^{\nu} p_{\nu}) U(\Lambda) \Psi_{p, \sigma}
	\end{equation*}

that :math:`U(\Lambda) \Psi_{p, \sigma}` has :math:`4`-momentum :math:`\Lambda p`. Therefore we can write

.. math::
	:nowrap:

	\begin{equation}
		U(\Lambda) \Psi_{p, \sigma} = C_{\sigma \sigma'} (\Lambda, p) \Psi_{\Lambda p, \sigma'}
		\label{eq_lorentz_acts_on_p_and_sigma}
	\end{equation}

where :math:`C_{\sigma \sigma'}` furnishes a representation of :math:`\Lambda` and :math:`p` under straightforward transformation rules, and an implicit summation over :math:`\sigma'` is assumed although it's not a :math:`4`-index.

Next we'd like to remove the dependency of :math:`C_{\sigma \sigma'}` on :math:`p` since, after all, it is :math:`\Lambda` that carries the symmetry. We can achieve this by noticing that :math:`U(\Lambda)` acts on the :math:`\Lambda`-orbits of :math:`p` transitively. The :math:`\Lambda`-orbits of :math:`p`, in turn, are uniquely determined by the value of :math:`p^2`, and in the case of :math:`p^2 \geq 0`, also by the sign of :math:`p_0`. We can therefore pick a convenient representative :math:`k` for each case as follows

+--------------------------------+-----------------------+----------+
| Case                           | Standard :math:`k`    | Physical |
+================================+=======================+==========+
| :math:`p^2 = M^2 > 0,~p_0 > 0` | :math:`(M, 0, 0, 0)`  | Yes      |
+--------------------------------+-----------------------+----------+
| :math:`p^2 = M^2 > 0,~p_0 < 0` | :math:`(-M, 0, 0, 0)` | No       |
+--------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 > 0`       | :math:`(1, 0, 0, 1)`  | Yes      |
+--------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 = 0`       | :math:`(0, 0, 0, 0)`  | Yes      |
+--------------------------------+-----------------------+----------+
| :math:`p^2 = 0,~p_0 < 0`       | :math:`(-1, 0, 0, 1)` | No       |
+--------------------------------+-----------------------+----------+
| :math:`p^2 = -N^2 < 0`         | :math:`(0, N, 0, 0)`  | No       |
+--------------------------------+-----------------------+----------+

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

At this point, we have reduced the problem to the classification of representations of the so-called *little group* defined as the subgroup of (proper orthochronous) Lorentz transformations :math:`W` that fixes :math:`k`, i.e., :math:`W_{\mu}^{\nu} k_{\nu} = k_{\mu}`. Element in the little group is known as `Wigner rotation <https://en.wikipedia.org/wiki/Wigner_rotation>`__ (and hence :math:`W`). More precisely, the task now is to find (unitary) representations :math:`D(W)` such that

.. math::
	:nowrap:

	\begin{equation*}
		D_{\sigma \sigma'}(W_1) D_{\sigma' \sigma''}(W_2) \Psi_{k, \sigma''} = D_{\sigma \sigma''}(W_1 W_2) \Psi_{k, \sigma''}
	\end{equation*}

Once this is done, we can define

.. math::
	:nowrap:

	\begin{align}
		U(W) \Psi_{k, \sigma} &\coloneqq D_{\sigma \sigma'}(W) \Psi_{k, \sigma'}  \label{eq_d_repr_of_little_group} \\
		W(\Lambda, p) &\coloneqq L(\Lambda p)^{-1} \Lambda L(p)  \label{eq_w_from_l}
	\end{align}


and substitute them into :math:`\eqref{eq_def_of_one_particle_psi_refactored}` to get

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda) \Psi_{p, \sigma} &\xlongequal{\eqref{eq_w_from_l}} N(p) U(L(\Lambda p)) U(W(\Lambda, p)) \Psi_{k, \sigma} \label{eq_little_group_acts_on_p_and_sigma} \\
			&\xlongequal{\eqref{eq_d_repr_of_little_group}} N(p) D_{\sigma \sigma'}(W(\Lambda, p)) U(L(\Lambda p)) \Psi_{k, \sigma'}  \nonumber \\
			&\xlongequal{\eqref{eq_def_of_one_particle_psi}} \frac{N(p)}{N(\Lambda p)} D_{\sigma \sigma'}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'} \nonumber
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

However, the Dirac delta in :math:`\eqref{eq_psi_p4_sigma_orthonormal}` is tricky to use since :math:`p` is constrained to the so-called *mass shell*, i.e., :math:`p_0 > 0` together with :math:`p^2 = M^2` in the massive case and :math:`p^2 = 0` in the massless case, respectively. Hence the actual normalization we'd like to impose on the one-particle states is, instead of :math:`\eqref{eq_psi_p4_sigma_orthonormal}`, the following

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
		\int d^4 p ~\delta(p^2 - M^2) \theta(p_0) f(p) &= \int d^3\pbf dp_0 ~\delta(p_0^2 - \pbf^2 - M^2) \theta(p_0) f(p_0, \pbf) \\
			&= \int d^3\pbf ~\frac{f\left( \sqrt{\pbf^2 + M^2}, \pbf \right)}{2 \sqrt{\pbf^2 + M^2}}
	\end{align*}

where :math:`\theta(p_0)` is the step function defined to be :math:`0` if :math:`p_0 \leq 0` and :math:`1` if :math:`p_0 > 1`. It follows that the Lorentz-invariant volume element in the :math:`3`-momentum space is

.. math::
	:nowrap:

	\begin{equation}
		d^3\pbf / \sqrt{\pbf^2 + M^2}  \label{eq_lorentz_invariant_3_momentum_volume_element}
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
		U(\Lambda) \Psi_{p, \sigma} = \sqrt{\frac{(\Lambda p)_0}{p_0}} D_{\sigma \sigma'}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'}
		\label{eq_lorentz_transformation_formula_for_particle_state}
	\end{equation}

where :math:`D_{\sigma \sigma'}` is a unitary representation of the little group, and :math:`W(\Lambda, p)` is defined by :math:`\eqref{eq_w_from_l}`.

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
		D^{(\jfrak)}_{\sigma \sigma'} (\Rcal) &= \delta_{\sigma \sigma'} + \tfrac{\ifrak}{2} \Theta^{ij} \left( J^{(\jfrak)}_{ij} \right)_{\sigma \sigma'}  \nonumber \\
		\left( J^{(\jfrak)}_{23} \pm \ifrak J^{(\jfrak)}_{31} \right)_{\sigma \sigma'} = \left( J^{(\jfrak)}_1 \pm \ifrak J^{(\jfrak)}_2 \right)_{\sigma \sigma'} &= \delta_{\sigma \pm 1, \sigma'} \sqrt{(\jfrak \mp \sigma)(\jfrak \pm \sigma + 1)}  \label{eq_j1_j2_matrix} \\
		\left( J^{(\jfrak)}_{12} \right)_{\sigma \sigma'} = \left( J^{(\jfrak)}_3 \right)_{\sigma \sigma'} &= \sigma \delta_{\sigma \sigma'}  \label{eq_j3_matrix}
	\end{align}

where :math:`\sigma, \sigma'` run through the values :math:`-\jfrak, -\jfrak + 1, \cdots, \jfrak - 1, \jfrak`.

.. _dropdown_repr_of_angular_momenta:

.. dropdown:: Representations of angular momenta
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
			\Jbf^2 \Psi_{\jfrak} = (\jfrak^2 + \jfrak) \Psi_{\jfrak} = ((\jfrak')^2 - \jfrak') \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'} ~\Longrightarrow~ \jfrak (\jfrak + 1) = \jfrak' (\jfrak' - 1)
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
		\pbf = \frac{M \vbf}{\sqrt{1 - \vbf^2}} ~\Longrightarrow~ \gamma = \frac{1}{\sqrt{1 - \vbf^2}} = \frac{\sqrt{M^2 + \pbf^2}}{M}
	\end{equation*}

It follows that [#boost_in_p_formula]_

.. math::
	:nowrap:

	\begin{align}
		L(p)_0^0 &= \gamma \label{eq_L_transformation_for_massive_1} \\
		L(p)_i^0 = L(p)_0^i &= \frac{p_i}{M} \label{eq_L_transformation_for_massive_2} \\
		L(p)_i^j &= \delta_i^j + \frac{p_i p_j}{\pbf^2} (\gamma - 1) \label{eq_L_transformation_for_massive_3}
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

		\begin{equation*}
			\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} = \sum_{\sigma' \sigma''} C^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''}
		\end{equation*}

	where the coefficients are known as Clebsch-Gordan coefficients. For the clarity of exposition, let's divide the solution into a few steps.

	Step 1.
		First of all, note that since :math:`J_3 = J'_3 + J''_3`, we have the following constraint

		.. math::
			:nowrap:

			\begin{equation}
				\Psi^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \neq 0 ~\Longrightarrow~ \sigma = \sigma' + \sigma''
				\label{eq_sigma_additive}
			\end{equation}

		Moreover, we see that the maximum possible value of :math:`\sigma` is :math:`\jfrak' + \jfrak''`, and it is achieved exactly when :math:`\sigma' = \jfrak'` and :math:`\sigma'' = \jfrak''`. It follows, assuming the non-degeneracy of the representation at least, that

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
		Next consider a state with :math:`\sigma = \jfrak' + \jfrak'' - 1`. It follows from :math:`\eqref{eq_sigma_additive}` that it must be a superposition of :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''}` and :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1}`, unless :math:`\jfrak'` and/or :math:`\jfrak''` vanishes, which leads to even simpler situations. Now we have two possible values of :math:`\jfrak`, namely :math:`\jfrak' + \jfrak''` and :math:`\jfrak' + \jfrak'' - 1`.

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

	1. The energy difference between :math:`1s_{1/2}` and :math:`2s_{1/2}`, plotted on a spectrometer, is the famous `21-centimeter line <https://en.wikipedia.org/wiki/Hydrogen_line>`_.
	2. The energy difference between :math:`2p_{1/2}` and :math:`2p_{3/2}`, i.e., same orbital but different total angular momentum, is known as the `fine structure <https://en.wikipedia.org/wiki/Fine_structure>`_ of the hydrogen atom.
	3. The energy difference between :math:`2s_{1/2}` and :math:`2p_{1/2}`, i.e., same total but different orbital angular momentum, is known as the `Lamb shift <https://en.wikipedia.org/wiki/Lamb_shift>`_.
	4. The energy difference between states with the same orbital and total angular momentum, e.g., :math:`1s_{1/2}`, but different spin :math:`z`-component :math:`\sigma`, e.g., :math:`\pm 1/2`, due to the magnetic moment is known as the `hyperfine structure <https://en.wikipedia.org/wiki/Hyperfine_structure>`_.

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

	\begin{equation*}
		S(a, b) = \begin{bmatrix}
			1 + c & a & b & -c \\
			a & 1 & 0 & -a \\
			b & 0 & 1 & -b \\
			c & a & b & 1 - c
		\end{bmatrix}
	\end{equation*}

which leaves :math:`k` invariant, and satisfies :math:`Sv = Wv`. It follows that :math:`S^{-1} W` must be a rotation about the :math:`3`-axis, which can be written as follows

.. math::
	:nowrap:

	\begin{equation*}
		R(\theta) = \begin{bmatrix}
			1 & 0 & 0 & 0 \\
			0 & \cos\theta & \sin\theta & 0 \\
			0 & -\sin\theta & \cos\theta & 0 \\
			0 & 0 & 0 & 1
		\end{bmatrix}
	\end{equation*}

Hence we can write any element in the little group as :math:`W(a, b, \theta) = S(a, b) R(\theta)`.

.. dropdown:: The little group is :math:`~E^+(2)`
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
		W(a, b, \theta)_{\mu}^{\nu} &= \left(1 + \begin{bmatrix}
				0 & a & b & 0 \\
				a & 0 & 0 & -a \\
				b & 0 & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix} \right) \left(1 + \begin{bmatrix}
				0 & 0 & 0 & 0 \\
				0 & 0 & \theta & 0 \\
				0 & -\theta & 0 & 0 \\
				0 & 0 & 0 & 0
			\end{bmatrix} \right) + \cdots \\
			&= 1 + \begin{bmatrix}
				0 & a & b & 0 \\
				a & 0 & \theta & -a \\
				b & -\theta & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix} + \cdots
	\end{align*}

where we've added the :math:`4`-indexes since we recall from discussions in :ref:`sec_quantum_lorentz_symmetry` that we must lift the :math:`\omega` index to make it anti-symmetric. We now rewrite

.. math::
	:nowrap:

	\begin{equation*}
		W(a, b , \theta)^{\mu \nu} = \eta^{\mu \sigma} W_{\sigma}^{\nu} = 1 + \begin{bmatrix}
				0 & -a & -b & 0 \\
				a & 0 & \theta & -a \\
				b & -\theta & 0 & -b \\
				0 & a & b & 0
			\end{bmatrix} + \cdots
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

	\begin{equation*}
		D_{\sigma \sigma'}(W(a, b, \theta)) = \exp(\ifrak \theta \sigma) \delta_{\sigma \sigma'}
	\end{equation*}

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
		U(\Pcal) &\Kbf U^{-1}(\Pcal) &&= -&&\Kbf \nonumber \\
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
		&\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma_1 \sigma'_1}(W_1(\Lambda, p_1)) D_{\sigma_2 \sigma'_2}(W_2(\Lambda, p_2)) \cdots \nonumber \\
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
		\exp(\ifrak \tau E_{\alpha}) \Psi_{\alpha} = \exp(\ifrak \tau H) \Psi_{\alpha} = \exp(\ifrak \tau (E_1 + E_2 + \cdots)) \Psi_{\alpha} ~\Longrightarrow~ E_\alpha = E_1 + E_2 + \cdots
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
		E_{\alpha} \Psi_{\alpha}^{\pm} = H \Psi_{\alpha}^{\pm} = (H_0 + V) \Psi_{\alpha}^{\pm} ~\Longrightarrow~ (E_{\alpha} - H_0) \Psi_{\alpha}^{\pm} = V \Psi_{\alpha}^{\pm}
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
		S_{\beta \alpha} = (\Omega(\infty) \Phi_{\beta}, \Omega(-\infty) \Phi_{\alpha}) = (\Phi_{\beta}, \Omega^{\dagger}(\infty) \Omega(-\infty) \Phi_{\alpha}) ~\Longrightarrow~ S = \Omega^{\dagger}(\infty) \Omega(-\infty) \eqqcolon U(\infty, -\infty)
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
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- \
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} + \int d\beta ~\Phi_{\beta} \int d\alpha \
				\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)}{E_{\alpha} - E_{\beta} + \ifrak \epsilon}  \label{eq_positive_limit_of_in_state_by_lippmann_schwinger} \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} \
     			- 2\pi\ifrak \int d\beta ~\Phi_{\beta} \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) \exp(-\ifrak \tau E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \left( g(\beta) - 2\pi\ifrak \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \right)  \nonumber \\
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
			&\times \sum_{\underline{\sigma}'_1 \underline{\sigma}'_2 \cdots} D^{\ast}_{\sigma'_1 \underline{\sigma}'_1} (W(\Lambda, p'_1)) D^{\ast}_{\sigma'_2 \underline{\sigma}'_2} (W(\Lambda, p'_2)) \cdots  \nonumber \\
			&\times \sum_{\underline{\sigma}_1 \underline{\sigma}_2 \cdots} D_{\sigma_1 \underline{\sigma}_1} (W(\Lambda, p_1)) D_{\sigma_2 \underline{\sigma}_2} (W(\Lambda, p_2)) \cdots  \nonumber \\
			&\times S_{\Lambda p'_1, \underline{\sigma}'_1, n'_1; ~\Lambda p'_2, \underline{\sigma}'_2, n'_2; ~\cdots, ~~\Lambda p_1, \underline{\sigma}_1, n_1; ~\Lambda p_2, \underline{\sigma}_2, n_2, ~\cdots}  \nonumber
	\end{align}

where we've used primes to distinguish between labels from in- and out-states, and underlines to distinguish between labels, specifically the spin-:math:`z` or helicity, before and after the Lorentz transformation.

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
			[\Kbf_0, \Omega(-\infty)] = -\Wbf \Omega(-\infty) ~\Longrightarrow~ \Kbf \Omega(\mp\infty) = \Omega(\mp\infty) \Kbf_0
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
		U(T) \Psi^{\pm}_{p_1, \sigma_1, n_1;~p_2, \sigma_2, n_2;~\cdots} = \sum_{n'_1, n'_2, \cdots} \Dscr_{n_1 n'_1} \Dscr_{n_2 n'_2} \cdots \Psi^{\pm}_{p_1, \sigma_1, n'_1;~p_2, \sigma_2, n'_2;~\cdots}
		\label{eq_internal_symmetry_transformation_for_in_and_out_states}
	\end{equation}

where :math:`U(T)` is the unitary operator associated with the symmetry transformation :math:`T`, and the :math:`\Dscr`'s are analogs of the little group representations from :math:`\eqref{eq_d_repr_of_little_group}`.

Similar to :math:`\eqref{eq_lorentz_transformation_formula_for_s_matrix}`, we can formulate the internal symmetry of S-matrix as follows

.. math::
	:nowrap:

	\begin{equation}
		S_{n'_1, n'_2, \cdots, ~n_1, n_2, \cdots} = \
			\sum_{\underline{n}'_1, \underline{n}'_2, \cdots} \Dscr^{\ast}_{n'_1 \underline{n}'_1}(T) \Dscr^{\ast}_{n'_2 \underline{n}'_2}(T) \cdots \
			\sum_{\underline{n}_1, \underline{n}_2, \cdots} \Dscr_{n_1 \underline{n}_1}(T) \Dscr_{n_2 \underline{n}_2}(T) \cdots \
			S_{\underline{n}'_1, \underline{n}'_2, \cdots, ~\underline{n}_1, \underline{n}_2, \cdots}
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

.. dropdown:: Parities of elementary particles
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

		and in addition :math:`\eqref{eq_fermion_count_eq_baryon_and_lepton_mod_2}` holds, then :math:`U(\Pcal)^2` is part of a continuous symmetry :math:`\eqref{eq_baryon_lepton_charge_internal_symmetry}` and hence can be set to one.  A hypothetical example that breaks :math:`\eqref{eq_fermion_count_eq_baryon_and_lepton_mod_2}` is the so-called Majorana fermions that are their own anti-particles, which implies :math:`B = L = 0`. For these particles, we have :math:`U(\Pcal)^4 = 1`, and hence the intrinsic parity may be :math:`\pm 1` or :math:`\pm \ifrak`.

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
			~\Longrightarrow~ \left| S^{(1)}_{\beta \alpha} \right|^2 = \left| S^{(1)}_{\Tcal \beta ~\Tcal \alpha} \right|^2
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
		u_{\alpha} &= \frac{\sqrt{(E_1 E_2 + |\pbf|^2)^2 - m_1^2 m_2^2}}{E_1 E_2} \label{eq_two_particles_relative_velocity_in_center_of_mass_frame} \\
			&= \frac{\sqrt{(E_1 E_2 + |\pbf|^2)^2 - (E_1^2 - |\pbf|^2)(E_2^2 - |\pbf|^2)}}{E_1 E_2} \nonumber \\
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
		k' &= \frac{\sqrt{\left( E_{\alpha}^2 - {m'_1}^2 - {m'_2}^2 \right)^2 - 4 {m'_1}^2 {m'_2}^2}}{2E_{\alpha}} \label{eq_defn_root_k_prime}
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

Here :math:`V(\tau) = \exp(\ifrak H_0 \tau) V \exp(-\ifrak H_0 \tau)` is a time-dependent operator in the so-called *interaction picture*, to be distinguished from the Heisenberg picture operator where the true Hamiltonian :math:`H` should be used in place of :math:`H_0`. The differential equation :math:`\eqref{eq_evolution_equation_of_u_operator}` can be easily solved as follows

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

	\begin{align*}
		T\{ V(t) \} &\coloneqq V(t) \\
		T\{ V(t_1) V(t_2) \} &\coloneqq \theta(t_1 - t_2) V(t_1) V(t_2) + \theta(t_2 - t_1) V(t_2) V(t_1) \\
		T\{ V(t_1) V(t_2) V(t_3) \} &\coloneqq \theta(t_1 - t_2) \theta(t_2 - t_3) V(t_1) V(t_2) V(t_3) + \cdots \\
		&\cdots
	\end{align*}

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

	\begin{equation*}
		V(t) = \int d^3 x ~\Hscr(t, \xbf)
	\end{equation*}

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
		S = 1 + \sum_{n=1}^{\infty} \frac{(-1)^n}{n!} \int d^4 x_1 \cdots d^4 x_n ~T\{ \Hscr(x_1) \cdots \Hscr(x_n) \}
		\label{eq_s_matrix_power_series_expansion_time_ordered_density}
	\end{equation}

This expression of :math:`S` is manifestly Lorentz invariant, except for the time-ordering part. In fact, the time-ordering between two spacetime points :math:`x_1, x_2` are Lorentz invariant if and only if :math:`x_1 - x_2` is time-like, namely, :math:`(x_1 - x_2)^2 \geq 0`. This is consistent with intuition because events with time-like (or light-like) separations may be observed by one observer, who definitely should know which event happened first. Therefore we obtain a sufficient condition for the Lorentz invariance of :math:`S` as follows

.. math::
	:nowrap:

	\begin{equation}
		[\Hscr(x_1), \Hscr(x_2)] = 0, \quad\forall ~(x_1 - x_2)^2 \leq 0
		\label{eq_h_commutativity_for_space_like_separations}
	\end{equation}

where we've also included the light-like case for technical reasons that will only become clear later.

.. dropdown:: A formal proof of the Lorentz invariance of the S-matrix
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

The Lorentz and CPT transformation laws
+++++++++++++++++++++++++++++++++++++++

Let's first work out how :math:`a^{\dagger}(\pbf, \sigma, n)` and :math:`a(\pbf, \sigma, n)` transform under Lorentz transformations. To this end, we recall the general Lorentz transformation law on free-particles state :math:`\eqref{eq_lorentz_transformation_formula_for_many_free_particles}`, and use :math:`b` for the translational parameter to avoid conflict of notations, as follows

.. math::
	:nowrap:

	\begin{align*}
        U_0(\Lambda, b) \Phi_{p_1, \sigma_1, n_1; ~p_2, \sigma_2, n_2; ~\cdots} =&~ \exp(-\ifrak b^{\mu} ((\Lambda p_1)_{\mu} + (\Lambda p_2)_{\mu} + \cdots)) \\
        &\times \sqrt{\frac{(\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p_1)_0 (p_2)_0 \cdots}} \\
        &\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma_1 \sigma'_1}(W_1(\Lambda, p_1)) D_{\sigma_2 \sigma'_2}(W_2(\Lambda, p_2)) \cdots \\
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
			\exp(-\ifrak b^{\mu} (\Lambda p)_{\mu}) \sqrt{\frac{(\Lambda p)_0}{p_0}} D_{\sigma \sigma'}(W(\Lambda, p)) a^{\dagger}(\pbf_{\Lambda}, \sigma', n)
		\label{eq_lorentz_transformation_formula_for_creation_operator}
	\end{equation}

where :math:`\pbf_{\Lambda}` is the :math:`3`-momentum part of :math:`\Lambda p`. The corresponding transformation law for the annihilation operator :math:`a(\pbf, \sigma, n)` can be obtained by taking the adjoint of :math:`\eqref{eq_lorentz_transformation_formula_for_creation_operator}` and remembering that :math:`U_0` is unitary as follows

.. math::
	:nowrap:

	\begin{equation}
		U_0(\Lambda, b) a(\pbf, \sigma, n) U_0^{-1}(\Lambda, b) = \
			\exp(\ifrak b^{\mu} (\Lambda p)_{\mu}) \sqrt{\frac{(\Lambda p)_0}{p_0}} D^{\ast}_{\sigma \sigma'}(W(\Lambda, p)) a(\pbf_{\Lambda}, \sigma, n)
		\label{eq_lorentz_transformation_formula_for_annihilation_operator}
	\end{equation}

The parity, time, and charge symmetry on the creation and annihilation operators can be worked out as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		U(\Pcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Pcal) &\xlongequal{\eqref{eq_space_inversion_on_massive_general}} \eta_n a^{\dagger}(-\pbf, \sigma, n) &&\quad\text{(massive particle)} \\
		U(\Pcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Pcal) &\xlongequal{\eqref{eq_space_inversion_on_massless_general}} \eta_{\sigma, n} \exp(\mp \ifrak \pi \sigma) a^{\dagger}(-\pbf, -\sigma, n) &&\quad\text{(massless particle)} \\
		U(\Tcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Tcal) &\xlongequal{\eqref{eq_time_inversion_on_massive_general}} \zeta_n (-1)^{\jfrak - \sigma} a^{\dagger}(-\pbf, -\sigma, n) &&\quad\text{(massive particle)} \\
		U(\Tcal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Tcal) &\xlongequal{\eqref{eq_time_inversion_on_massless_general}} \zeta_{\sigma, n} \exp(\pm \ifrak \pi \sigma) a^{\dagger}(-\pbf, \sigma, n) &&\quad\text{(massless particle)} \\
		U(\Ccal) a^{\dagger}(\pbf, \sigma, n) U^{-1}(\Ccal) &\xlongequal{\phantom{(00)}} \zeta_n a^{\dagger}(\pbf, \sigma, n^c) &&\quad\text{(massive/massless particle)}
	\end{alignat*}

where :math:`\Ccal` replaces a particle of species :math:`n` with its anti-particle -- a notion we haven't really explained yet -- :math:`n^c`. The corresponding transformation laws for :math:`a(\pbf, \sigma, n)` can be derived from the above identities by taking the adjoint, and are omitted here.

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

Let's write :math:`S^{(n)}_{\beta \alpha} \coloneqq \left( \Phi_{\beta}, T\{V(t_1) \cdots V(t_n)\} \Phi_{\alpha} \right)`. Then the general recipe for constructing a Feynman diagram to keep track of the delta functions in :math:`S^{(n)}_{\beta \alpha}` consists of the following steps

1. Orient the paper on which the diagram will be drawn so that the time direction goes upwards. (This is of course a random choice just to set up the scene.)
2. Place as many vertical strands as the number of particles in :math:`\Phi_{\alpha}` at the bottom of the page. Do the same to :math:`\Phi_{\beta}` but at the top of the page. Moreover, all such strands are oriented to flow upwards.
3. Place one vertex (or a fat dot) for each :math:`V(t)` so that their are sorted vertically by time. Moreover, each vertex has the same number of incoming strands from below as the total number of annihilation operators in the corresponding :math:`V(t)`, written in the form :math:`\eqref{eq_general_operator_expansion_in_creation_and_annihilation}`, and the same number of outgoing strands pointing upwards as the total number of creation operators. Note that since we'll be doing a rather crude analysis of just momentum conservations, the actual structure of :math:`V(t)` plays no role here.
4. To each strand from the previous step, associate a :math:`3`-momentum :math:`\pbf` corresponding to the created/annihilated particle. To each vertex, associate a momentum-conservation delta function so that the total incoming momenta equals the total outgoing momenta. This is in fact a consequence of an assumption on :math:`V(t)`, which is related to our assumption on :math:`h_{NM}` in :math:`\eqref{eq_general_expansion_of_hamiltonian}`, and will be elaborated on later in this section.
5. Connect all the loose ends of the strands so that the flow always goes in the positive time direction, i.e., upwards. All loose ends must be connected to qualify as a valid (directed) graph. As we'll see, invalid graphs don't contribute to :math:`S^{(n)}_{\beta \alpha}`.
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

where we've also suppressed the time factor from :math:`V(t)` as it won't really play a role in the diagram. The figure below illustrates a few summands of the third order :math:`S^{(3)}_{\beta \alpha} = \left( \Phi_{\beta}, T\{V(t_1)V(t_2)V(t_3)\}\Phi_{\alpha} \right)`.

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

To verify that each connected Feynman diagram gives rise to exactly one momentum-conservation delta function, we note a crucial fact that the diagram contains no *directed* loops, which allows us to start from the bottom, i.e., the total momentum of the in-state particles, and work our way up to the top, i.e., the total momentum of the out-state particles, so that there is exactly one momentum-conservation delta function at each step when a(n intermediate) vertex is added, which guarantees that the incoming total momentum equals the outgoing total momentum.

As a side remark, we point out that the argument in [Wei95]_ (page 186 --187) in terms of a Euler characteristic calculation of the diagram appears to be wrong. A counterexample (to Weinberg's argument) can be illustrated in the following pseudo Feynman diagram which contains a directed loop

.. figure:: ./static/quantum-theory-of-fields/pseudo-feynman-diagram.svg
	:align: center

	Figure: A pseudo Feynman diagram where time orientation is not respected, and the edge-wise delta functions have been evaluated so that the momentum of an intermediate particle can be put on an edge. In particular, there is a directed closed loop.

Evaluating the momentum-conservation delta function of the diagram then gives

.. math::
	:nowrap:

	\begin{equation*}
		\delta^3(\pbf_3 - \pbf_1 - \pbf_4) \delta^3(\pbf_5 - \pbf_2 - \pbf_3) \delta^3(\pbf_6 + \pbf_4 - \pbf_5) \
			= \delta^3(\pbf_6 - \pbf_1 - \pbf_2) \delta^3(\pbf_4 + \pbf_6 - \pbf_2 - \pbf_3)
	\end{equation*}

where the first factor demands the usual conservation of momentum between in- and out-state particles, and the second demands a conservation of momentum on the edge labeled by :math:`\pbf_5`.




.. rubric:: Footnotes

.. [#tedious_calc_of_commutations] These are some rather tedious and error-prone calculations, but in the end, we manage to arrive at the same results as stated in [Wei95]_ page 61.

.. [#boost_in_p_formula] The formula in [Wei95]_ page 68, eq. (2.5.24) is wrong.

.. [#in_out_state_sign_convention] I really don't like the sign convention of the in- and out-states in [Wei95]_, even though Weinberg said that it has become a tradition. So our sign convention, as far as the definitions of in- and out-states are concerned, is the opposite to the one used in [Wei95]_.

.. [#pion_deuteron_reaction_final_state] The final state claimed in [Wei95]_ page 126 has total spin :math:`1` rather than :math:`0`. I suspect that the claim in [Wei95]_ is wrong, but otherwise it doesn't affect any subsequent arguments anyway.

.. [#abuse_of_phi_as_both_state_vector_and_flux] It's really unfortunate that one constantly runs out symbols to represent physical quantities, and it's uncommon to use anything other than just one letter. So we have to live with the fact that the meaning of the symbol will depend on the context.

.. [#weinberg_quote_on_lorentz_invariance_and_quantum_mechanics] See [Wei95]_ page 145.