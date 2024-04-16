The Quantum Theory of Fields (S. Weinberg)
==========================================

.. warning::
	This note is work in progress.

This note covers the three volumes [Wei95]_, [Wei96]_, and [Wei00]_, written by S. Weinberg on quantum field theory. What I like the most about these books is his attempt to make logical deductions from the most basic principles, in particular, the principle of symmetries, rather than to make analogies to experience, e.g., from classical physics (from a historical perspective). Such an endeavor may not always be possible because, after all, physics is about how we interpret Nature based on nothing but experience, and *not* about how Nature actually works. By the same token, the arguments that are considered logical here should really be interpreted as "reasonable-sounding", and have nothing to do with what mathematician would call "rigorous".

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

		It follows that :math:`c = \pm\i`, which correspond to :math:`U` being (complex) linear or anti-linear, respectively.

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
	While the argument here appears to be purely mathematical, Weinberg pointed out in the book the potential inabilities to create a state like :math:`\Psi_A + \Psi_B`. More precisely, he mentioned the general believe that it's impossible to prepare a superposition of two states, one with integer total angular momentum and the other with half-integer total angular momentum, in which case there will be a "super-selection rule" between different classes of states. After all, one Hilbert space may just not be enough to describe all states. It'd be nice to elaborate a bit more on the super-selection rules.

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

where :math:`\theta^a` are the (real) components of :math:`\theta`, and :math:`u_a` are operators independent of :math:`\theta`, and as a convention, repeated indexes are summed up. Here we put a :math:`\i` in front of the linear term so that the unitarity of :math:`U` implies that :math:`u_a` are Hermitian.

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
		t' &=  \Lambda_0^{\mu} x_{\mu} &&= \gamma (t + v_3 x_3)  \nonumber \\
		x'_1 &= \Lambda_1^{\mu} x_{\mu} &&= x_1  \nonumber \\
		x'_2 &= \Lambda_2^{\mu} x_{\mu} &&= x_2  \nonumber \\
		x'_3 &= \Lambda_3^{\mu} x_{\mu} &&= \gamma (x_3 + v_3 t)
		\label{eq_lambda_in_3_axis}
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
		\text{let } \mu = 0, \rho = 0, \kappa = i \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [H, K_i] &&= \ifrak P_i  \nonumber \\
		\text{let } \mu = i, \nu = j \text{ in \eqref{eq_bracket_p4_p4}} ~&\Longrightarrow~ [P_i, P_j] &&= 0  \label{eq_pp_commute} \\
		\text{let } \mu = i, \rho = k, \kappa = i \text{ in \eqref{eq_bracket_p4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [P_i, J_j] &&= -\ifrak \epsilon_{ijk} P_k  \nonumber \\
		\text{let } \mu = i, \rho = 0 \text{ and enumerate } \kappa \in \{1, 2, 3\} \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [P_i, K_j] &&= \ifrak \delta_{ij} H  \nonumber \\
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
		U(\Lambda) \Psi_{p, \sigma} &\xlongequal{\eqref{eq_w_from_l}} N(p) U(L(\Lambda p)) U(W(\Lambda, p)) \Psi_{k, \sigma}  \nonumber \\
			&\xlongequal{\eqref{eq_d_repr_of_little_group}} N(p) D_{\sigma \sigma'}(W(\Lambda, p)) U(L(\Lambda p)) \Psi_{k, \sigma'}  \nonumber \\
			&\xlongequal{\eqref{eq_def_of_one_particle_psi}} \frac{N(p)}{N(\Lambda p)} D_{\sigma \sigma'}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'}
			\label{eq_little_group_acts_on_p_and_sigma}
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

where :math:`\theta(p_0)` is the step function defined to be :math:`0` if :math:`p_0 \leq 0` and :math:`1` if :math:`p_0 > 1`. It follows that the Lorentz-invariant volume element in the :math:`3`-momentum space is :math:`d^3\pbf / \sqrt{\pbf^2 + M^2}`. We can use it to find the Lorentz-invariant Dirac delta (marked in blue) as follows

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

	In other words, applying :math:`J_1 \pm \ifrak J_2` to any eigenstate of :math:`J_3` raises or lowers the eigenvalue by one, and henceforth they are called *raising* and *lowering* operators, respectively. Moreover, since :math:`\Jbf^2` commutes with :math:`J_3`, we may assume that :math:`\Psi_m` is also an eigenstate of :math:`\Jbf^2`, and since :math:`\Jbf^2` also commutes with :math:`J_1 \pm \ifrak J_2`, the whole series of :math:`J_3`-eigenstates obtained by applying raising and/or lowering operators has the same :math:`\Jbf^2`-eigenvalue.

	We are only interested in eigenstates with finite :math:`\Jbf^2`-eigenvalue. So both the raising and the lowering operations must stop after finite steps. Let :math:`\Psi_{\jfrak}` be the :math:`J_3`-eigenstate with the highest eigenvalue (if there are more than one, the representation is reducible). By repeatedly applying the lowering operator to :math:`\Psi_{\jfrak}`, we'll eventually reach the eigenstate :math:`\Psi_{\jfrak'}` with the lowest eigenvalue. Since the lowering operator decreases the eigenvalue by one, we know that :math:`\jfrak - \jfrak'` must be an integer.

	Let's consider the following two operators

	.. math::
		:nowrap:

		\begin{align}
			(J_1 - \ifrak J_2) (J_1 + \ifrak J_2) &= J_1^2 + J_2^2 + \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 - J_3  \label{eq_j1_j2_mixed_product_plus} \\
			(J_1 + \ifrak J_2) (J_1 - \ifrak J_2) &= J_1^2 + J_2^2 - \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 + J_3  \label{eq_j1_j2_mixed_product_minus}
		\end{align}

	Note that the first operator annihilates :math:`\Psi_{\jfrak}` and the second operator annihilates :math:`\Psi_{\jfrak'}` by assumption, which, together with the fact that :math:`\Jbf^2 \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'}`, implies

	.. math::
		:nowrap:

		\begin{equation*}
			\Jbf^2 \Psi_{\jfrak} = (\jfrak^2 + \jfrak) \Psi_{\jfrak} = ((\jfrak')^2 - \jfrak') \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'} ~\Longrightarrow~ \jfrak (\jfrak + 1) = \jfrak' (\jfrak' - 1)
		\end{equation*}

	The equation has two potential solutions: either :math:`\jfrak' = \jfrak + 1` or :math:`\jfrak = -\jfrak'`. The first option violates the maximality of :math:`\jfrak`, and so we must accept the second option. Since we also know :math:`\jfrak - \jfrak'` must be integral, we conclude that :math:`\jfrak` is itself a half-integer.

	It remains to settle the constant term on the right-hand-side of :math:`\eqref{eq_j1_j2_matrix}`. By :math:`\eqref{eq_j1_j2_raises_or_lowers_state}` we can assume

	.. math::
		:nowrap:

		\begin{equation*}
			\left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi_{\sigma} = \alpha_{\pm}(\jfrak, \sigma) \Psi_{\sigma \pm 1}
		\end{equation*}

	Applying :math:`\eqref{eq_j1_j2_mixed_product_plus}` and :math:`\eqref{eq_j1_j2_mixed_product_minus}` to :math:`\Psi_{\sigma}` then implies

	.. math::
		:nowrap:

		\begin{equation*}
			\alpha_{\mp} (\jfrak, \sigma \pm 1) \alpha_{\pm} (\jfrak, \sigma) = \jfrak^2 + \jfrak - \sigma^2 \mp \sigma
		\end{equation*}

	Now we use the fact that :math:`J_i, i = 1, 2, 3`, are Hermitian operators to calculate

	.. math::
		:nowrap:

		\begin{align*}
			|\alpha_{\pm} (\jfrak, \sigma)|^2 (\Psi_{\sigma}, \Psi_{\sigma}) &= \left( \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi_{\sigma}, \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi_{\sigma} \right) \\
				&= \left( \Psi_{\sigma}, \left( J_1^{(\jfrak)} \mp \ifrak J_2^{(\jfrak)} \right) \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi_{\sigma} \right) \\
				&= (\jfrak^2 + \jfrak - \sigma^2 \mp \sigma) (\Psi_{\sigma}, \Psi_{\sigma})
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

Finally, we note an important fact that when :math:`\Lambda = \Rcal` is a :math:`3`-rotation, then :math:`W(\Rcal, p) = \Rcal` for any :math:`p`. To see this, we'll work out how :math:`W(\Rcal, p)` acts on :math:`(1, \mathbf{0}), (0, \pbf)`, and :math:`(0, \qbf)`, respectively, where :math:`\qbf` is any :math:`3`-vector perpendicular to :math:`\pbf`, as follows

.. math::
	:nowrap:

	\begin{alignat*}{2}
		W(\Rcal, p)(1, \mathbf{0}) &= L(\Rcal p)^{-1} \Rcal (\gamma, \pbf / M) &&= L(\Rcal p)^{-1} (\gamma, \Rcal p / M) &&= (1, \mathbf{0}) \\
		W(\Rcal, p)(0, \pbf) &= L(\Rcal p)^{-1} \Rcal (\pbf^2 / M, \gamma \pbf) &&= L(\Rcal p)^{-1} (\pbf^2 / M, \gamma \Rcal p) &&= (0, \Rcal \pbf) \\
		W(\Rcal, p)(0, \qbf) &= L(\Rcal p)^{-1} \Rcal (0, \qbf) &&= L(\Rcal p)^{-1} (0, \Rcal \qbf) &&= (0, \Rcal \qbf)
	\end{alignat*}

where we have used that fact that :math:`\gamma` is :math:`\Rcal`-invariant.

This observation is important since it implies that non-relativistic calculations about spins, such as the `Clebsch-Gordan coefficients <https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients>`__, can be literally carried over to the quantum setting.

.. todo::
	Review Clebsch-Gordan coefficients.


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

According to the book, massless particle states are not observed to come in such :math:`S^1`-families. Hence the only possibility is that :math:`a = b = 0` and the only symmetry left then is :math:`J_3`, which corresponds to a rotation about the :math:`3`-axis.

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

In particular, we see that, unlike the spin of massive particles, helicity is Lorentz invariant (at least under genuine representations). It is reasonable, therefore, to think of massless particles of different helicity as different particle species. Examples include photons with :math:`\sigma = \pm 1` and gravitons with :math:`\sigma = \pm 2`, but *not* (anti-)neutrinos with hypothetical :math:`\sigma = \pm \tfrac{1}{2}` as otherwise stated in the book, which are now known to have a nonzero mass. Here the :math:`\pm` signs are related to the space-inversion symmetry :math:`\eqref{eq_space_inversion}`, which will be discussed in detail later.

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
		U(\Tcal) \Psi_{p, \sigma} = \zeta (-1)^{\jfrak - \sigma} \Psi_{p, -\sigma}
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
		U(\Pcal) \Psi_{p, \sigma} &= p_0^{-1/2} U(\Pcal R(\pbf)B) \Psi_{k, \sigma} \nonumber \\
			&= p_0^{-1/2} U(R(\pbf) B R_2^{-1}) U(R_2 \Pcal) \Psi_{k, \sigma} \nonumber \\
			&= p_0^{-1/2} \eta_{\sigma} U(R(\pbf) R_2^{-1} B) \Psi_{k, -\sigma} \nonumber \\
			&= \eta_{\sigma} \rho \Psi_{\Pcal p, -\sigma}
			\label{eq_space_inversion_on_massless_undetermined_phase}
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

Scattering Theory
-----------------

Physics would have been rather boring if nothing interacts, like the free particles that we have been studying so far. On the flip side, physics would have been impossible if we try to know exactly what happens in the interactions. The middle ground, where we assume that the particles are non-interacting long before and after the interaction, and something mysterious happened in between, is called scattering theory.

Non-interacting many-particles state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We shall, as always, start from the easiest part of the theory, which is clearly the non-interacting parts. Recall our grand formulae for the Lorentz transformation on one-particle state :math:`\eqref{eq_translation_formula_for_particle_state}` and :math:`\eqref{eq_lorentz_transformation_formula_for_particle_state}`. For a non-interacting many-particles system, it's conceivable to assume that the Lorentz transformation law is simply a direct product of the individual particles as follows (recall that :math:`U(\Lambda, a) = U(1, a) U(\Lambda, 0)`)

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda, a) \Psi_{p_1, \sigma_1, n_1; p_2, \sigma_2, n_2; \cdots} =&~ \exp(-\ifrak a^{\mu} ((\Lambda p_1)_{\mu} + (\Lambda p_2)_{\mu} + \cdots)) \nonumber \\
		&\times \sqrt{\frac{(\Lambda p_1)_0 (\Lambda p_2)_0 \cdots}{(p_1)_0 (p_2)_0 \cdots}} \nonumber \\
		&\times \sum_{\sigma'_1 \sigma'_2 \cdots} D_{\sigma_1 \sigma'_1}(W_1(\Lambda, p_1)) D_{\sigma_2 \sigma'_2}(W_2(\Lambda, p_2)) \cdots \Psi_{\Lambda p_1, \sigma'_1, n_1; \Lambda p_2, \sigma'_2, n_2; \cdots}
		\label{eq_lorentz_transformation_formula_for_many_free_particles}
	\end{align}

where the first component is the translation transformation :math:`\eqref{eq_translation_formula_for_particle_state}`, the second component is the normalization factor, and the third component is the little group representation, and the :math:`\sigma`'s are either the spin-:math:`z` component for massive particles or the helicity for massless particles, and the :math:`n`'s are additional (discrete) labels such as mass, charge, spin, etc.

Notice that by writing a many-particles state as :math:`\Psi_{p_1, \sigma_1, n_1; p_2, \sigma_2, n_2; \cdots}`, we have given the particles an order, which is by no means unique. Hence the normalization of these states must take permutations into account as follows

.. math::
	:nowrap:

	\begin{equation}
		\left( \Psi_{p_1, \sigma_1, n_1; p_2, \sigma_2, n_2; \cdots}, \Psi_{p'_1, \sigma'_1, n'_1; p'_2, \sigma'_2, n'_2; \cdots} \right) = \delta^3(\pbf_1 - \pbf'_1) \delta_{\sigma_1 \sigma'_1} \delta_{n_1 n'_1} \delta^3(\pbf_2 - \pbf'_2) \delta_{\sigma_2 \sigma'_2} \delta_{n_2 n'_2} \pm \text{permutations}
		\label{eq_many_particles_state_normalization_rough}
	\end{equation}

The sign in front of the permutations has to do with the species of the particles, which will be discussed later. Note that although there are many terms in :math:`\eqref{eq_many_particles_state_normalization_rough}`, there is at most one nonzero term, which happens exactly when the two states differ by a permutation.

To suppress the annoyingly many subscripts in the states, we shall use letters such as :math:`\alpha, \beta, \cdots` to denote the compound index such as :math:`(p_1, \sigma_1, n_1; p_2, \sigma_2, n_2; \cdots)`, so that, for example, :math:`\eqref{eq_many_particles_state_normalization_rough}` can be simplified as

.. math::
	:nowrap:

	\begin{equation*}
		\left( \Psi_{\alpha}, \Psi_{\alpha'} \right) = \delta(\alpha - \alpha')
	\end{equation*}

where the integral volume element reads

.. math::
	:nowrap:

	\begin{equation*}
		\int d\alpha \cdots = \sum_{\sigma_1, n_1; \sigma_2, n_2; \cdots} \int d^3 \pbf_1 d^3 \pbf_2 \cdots
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

	\begin{equation*}
		\int d\alpha ~g(\alpha) \Psi_{\alpha}
	\end{equation*}

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
     			+ \int d\alpha d\beta ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm}) \Phi_{\beta}}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon} \nonumber \\
		    &= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Phi_{\alpha} \
   			    + \int d\beta ~\Phi_{\beta} \blue{\int d\alpha ~\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^{\pm})}{E_{\alpha} - E_{\beta} \mp \ifrak \epsilon}}  \label{eq_packet_expansion_by_lippmann_schwinger}
 	\end{align}

Now the integral colored in blue can be integrated over :math:`E_{\alpha}` by a contour that runs from :math:`-\infty` to :math:`+\infty`, followed by a semicircle at infinity, in the upper-half-plane in the case of :math:`\Psi_{\alpha}^-` and the lower-half-plane in the case of :math:`\Psi_{\alpha}^+`, back to :math:`-\infty`. In either case, the sign in :math:`\mp \ifrak \epsilon` is chosen so that the integrant has no poles with infinitesimally small imaginary part, though both :math:`g(\alpha)` and :math:`(\Phi_{\beta}, V \Psi_{\alpha}^{\pm})`, viewed as complex functions, may have poles with finite imaginary parts. It follows then from the residual theorem and the damping factor :math:`\exp(-\ifrak \tau E_{\alpha})` as :math:`\tau \to \pm\infty` that the integral in blue vanishes, as desired.

S-matrix and its symmetry
^^^^^^^^^^^^^^^^^^^^^^^^^

The `S-matrix <https://en.wikipedia.org/wiki/S-matrix>`_ defined by

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} \coloneqq \left( \Psi_{\beta}^+, \Psi_{\alpha}^- \right)
		\label{eq_defn_s_matrix_by_in_and_out_states}
	\end{equation}

records the probability amplitude of finding the out-state :math:`\Psi_{\beta}^+` given the in-state :math:`\Psi_{\alpha}^-`. Note that since the in- and out-states both form an orthonormal basis of the same Hilbert space, the S-matrix is unitary. However, the way :math:`S` is defined in :math:`\eqref{eq_defn_s_matrix_by_in_and_out_states}` disqualifies it as an operator on the Hilbert space. Therefore it's convenient to convert both in- and out-states to the free states and define the *S-operator* by

.. math::
	:nowrap:

	\begin{equation}
		(\Phi_{\beta}, S \Phi_{\alpha}) \coloneqq S_{\beta \alpha}
	\end{equation}

Using :math:`\eqref{eq_defn_of_Omega}` we see that

.. math::
	:nowrap:

	\begin{equation*}
		S_{\beta \alpha} = (\Omega(\infty) \Phi_{\beta}, \Omega(-\infty) \Phi_{\alpha}) = (\Phi_{\beta}, \Omega^{\dagger}(\infty) \Omega(-\infty) \Phi_{\alpha}) ~\Longrightarrow~ S = \Omega^{\dagger}(\infty) \Omega(-\infty) \eqqcolon U(\infty, -\infty)
	\end{equation*}

where

.. math::
	:nowrap:

	\begin{equation}
		U(\tau_1, \tau_0) = \exp(\ifrak \tau_1 H_0) \exp(\ifrak (\tau_0 - \tau_1) H) \exp(-\ifrak \tau_0 H_0)
		\label{eq_defn_u_operator}
	\end{equation}

The most straightforward way to calculate :math:`S_{\beta \alpha}` is probably to use :math:`\eqref{eq_lippmann_schwinger_pure}` directly. However. this turns out to be rather complicated, and doesn't lead to a simple result. Instead, we shall follow a rather smart trick from the book as follows. First let's calculate the asymptotic of the in-packet as :math:`\tau \to \infty` (but omitting the :math:`\lim_{\tau \to \infty}` symbol) using :math:`\eqref{eq_packet_expansion_by_lippmann_schwinger}`

.. math::
	:nowrap:

	\begin{align}
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- \
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} + \int d\beta ~\Phi_{\beta} \int d\alpha \
				\frac{\exp(-\ifrak \tau E_{\alpha}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)}{E_{\alpha} - E_{\beta} + \ifrak \epsilon}  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) g(\beta) \Phi_{\beta} \
     			- 2\pi\ifrak \int d\beta ~\Phi_{\beta} \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) \exp(-\ifrak \tau E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-)  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \left( g(\beta) - 2\pi\ifrak \int d\alpha ~\delta(E_{\alpha} - E_{\beta}) g(\alpha) (\Phi_{\beta}, V \Psi_{\alpha}^-) \right)  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \left( \blue{\delta(\alpha - \beta) - 2\pi\ifrak (\Phi_{\beta}, V \Psi_{\alpha}^-)} \right)  \label{eq_positive_limit_of_in_state_by_lippmann_schwinger}
	\end{align}

where we've used the residue theorem again in the second equality. Next expand the left-hand-side of the equation, i.e., the :math:`\tau \to \infty` limit of the in-packet in terms of the out-states

.. math::
	:nowrap:

	\begin{align}
		\int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \Psi_{\alpha}^- \
			&= \int d\alpha ~\exp(-\ifrak \tau E_{\alpha}) g(\alpha) \int d\beta ~(\Psi_{\beta}^+, \Psi_{\alpha}^-) \Psi_{\beta}^+  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Psi_{\beta}^+ \int d\alpha ~g(\alpha) S_{\beta \alpha}  \nonumber \\
			&= \int d\beta ~\exp(-\ifrak \tau E_{\beta}) \Phi_{\beta} \int d\alpha ~g(\alpha) \blue{S_{\beta \alpha}}
			\label{eq_positive_limit_of_in_state_by_expanding_out_states}
	\end{align}

where we've used the fact that the S-matrix contains a :math:`\delta(E_{\alpha} - E_{\beta})` factor by energy conservation in the second equality, and the defining property :math:`\eqref{eq_in_out_states_asymptotic_by_energy}` of the out-state in the third equality.

Equating the blue terms from :math:`\eqref{eq_positive_limit_of_in_state_by_lippmann_schwinger}` and :math:`\eqref{eq_positive_limit_of_in_state_by_expanding_out_states}`, we've derived the following formula

.. math::
	:nowrap:

	\begin{equation}
		S_{\beta \alpha} = \delta(\alpha - \beta) - 2\pi\ifrak (\Phi_\beta, V \Psi_{\alpha}^-)
	\end{equation}

Up to the first order of :math:`V`, one can replace the :math:`\Psi_{\alpha}^-` on the right-hand-side by :math:`\Phi_{\alpha}`, and arrive at the `Born approximation <https://en.wikipedia.org/wiki/Born_approximation>`_ of the S-matrix.

.. rubric:: Footnotes

.. [#tedious_calc_of_commutations] These are some rather tedious and error-prone calculations, but in the end, we manage to arrive at the same results as stated in the book.

.. [#boost_in_p_formula] The formula in [Wei95]_ page 68, eq. (2.5.24) is wrong.

.. [#in_out_state_sign_convention] I really don't like the sign convention of the in- and out-states in the book, even though Weinberg said that it's become a tradition. So our convention is the opposite to the one used in the book.