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

	Throughout this post we always assume that state vectors are normalized so that :math:`||\Psi|| = 1`.

In fact, we don't really care about the states themselves either, because they are more of an abstraction rather than something one can physically measure. What we do care about are the (Hermitian) inner products between state vectors, denoted by :math:`(\Psi, \Phi)`. According to the so-called `Copenhagen interpretation <https://en.wikipedia.org/wiki/Copenhagen_interpretation>`_ of quantum mechanics, such inner product represents an *amplitude*, i.e., its squared norm gives the probability of finding a state :math:`[\Psi]` in :math:`[\Phi]` if we ever perform a measurement. We can write this statement as an equation as follows:

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

	:Step 4: We have so far specified the value of :math:`U` on all of :math:`\Psi_i, i \geq 0`, and :math:`\Psi_0 + \Psi_k, k \geq 1`. Notice that all the coefficients of :math:`\Psi` are real. It is therefore instructive to ask what :math:`\Psi_0 + \i \Psi_1` should be. By the same argument as in the previous step, we can write

		.. math::
			:nowrap:

			\begin{equation*}
				T[\Psi_0 + \i \Psi_1] = [\Psi'_0 + c \Psi'_1]
			\end{equation*}

		where :math:`c` is a phase. Let's apply :math:`\eqref{eq_t_preserves_probability}` once again as follows

		.. math::
			:nowrap:

			\begin{align*}
				\sqrt{2} &= \left| \left( [\Psi_0 + \i \Psi_1], [\Psi_0 + \Psi_1] \right) \right| \\
					&= \left| \left( T[\Psi_0 + \i \Psi_1], T[\Psi_0 + \Psi_1] \right) \right| \\
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
		U(T_2 T_1) \Psi = e^{\i \theta(T_1, T_2, \Psi)} U(T_2) U(T_1) \Psi
		\label{eq_u_depends_on_psi}
	\end{equation}

where :math:`\theta(T_1, T_2, \Psi)` is an angle, which depends a priori on :math:`T_1, T_2`, and :math:`\Psi`.

It turns out, however, the angle :math:`\theta(T_1, T_2, \Psi)` cannot depend on the state because if we apply :math:`\eqref{eq_u_depends_on_psi}` to the sum of two linearly independent state vectors :math:`\Psi_A + \Psi_B`, then we'll find

.. math::
	:nowrap:

	\begin{equation*}
		e^{\pm \i \theta(\Psi_A)} \Psi_A + e^{\pm \i \theta(\Psi_B)} \Psi_B = e^{\pm \i \theta(\Psi_A + \Psi_B)} (\Psi_A + \Psi_B)
	\end{equation*}

where we have suppressed the dependency of :math:`\theta` on :math:`T`, and the signs correspond to the cases of :math:`U` being linear or anti-linear, respectively. In any case, it follows that

.. math::
	:nowrap:

	\begin{equation*}
		e^{\pm \i \theta(\Psi_A)} = e^{\pm \i \theta(\Psi_B)} = e^{\pm \i \theta(\Psi_A + \Psi_B)}
	\end{equation*}

which says nothing but the independence of :math:`\theta` on :math:`\Psi`.

.. todo::
	While the argument here appears to be purely mathematical, Weinberg pointed out in the book the potential inabilities to create a state like :math:`\Psi_A + \Psi_B`. More precisely, he mentioned the general believe that it's impossible to prepare a superposition of two states, one with integer total angular momentum and the other with half-integer total angular momentum, in which case there will be a "super-selection rule" between different classes of states. After all, one Hilbert space may just not be enough to describe all states. It'd be nice to elaborate a bit more on the super-selection rules.

We can now simplify :math:`\eqref{eq_u_depends_on_psi}` to the following

.. math::
	:nowrap:

	\begin{equation}
		U(T_2 T_1) = e^{\i \theta(T_1, T_2)} U(T_2) U(T_1)
		\label{eq_u_not_depend_on_psi}
	\end{equation}

which, in mathematical terms, says that :math:`U` furnishes a *projective representation* of :math:`T`, or a representation up to a phase. It becomes a genuine representation if the phase is constantly one.

.. note::
	We always assume that :math:`U` furnishes a genuine representation of :math:`T` since, as it turns out, the phase factor in :math:`\eqref{eq_u_not_depend_on_psi}` is more of a mathematical artifact than something that bears any physical significance.

.. _sec_continuous_symmetry:

Continuous symmetry
^^^^^^^^^^^^^^^^^^^

Besides a handful of important discrete symmetries such as the time, charge, and parity conjugations, most of the interesting symmetries come in a continuous family, mathematically known as *Lie groups*. Note that continuous symmetries are necessarily unitary (and linear) because they can be continuously deformed into the identity, which is obviously unitary.

In fact, it will be of great importance to just look at the symmetry up to the first order at the identity transformation, mathematically known as the *Lie algebra*. Let :math:`\theta` be an element in the Lie algebra such that :math:`T(\theta) = 1 + \theta` up to the first order. We can expand :math:`U(T(\theta))` in a power series as follows

.. math::
	:nowrap:

	\begin{equation}
		U(T(\theta)) = 1 + \i \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots
		\label{eq_u_expansion}
	\end{equation}

where :math:`\theta^a` are the (real) components of :math:`\theta`, and :math:`u_a` are operators independent of :math:`\theta`, and as a convention, repeated indexes are summed up. Here we put a :math:`\i` in front of the linear term so that the unitarity of :math:`U` implies that :math:`u_a` are Hermitian.

Now let :math:`\eta` be another element of the Lie algebra, and expand both sides of :math:`U(T(\eta)) U(T(\theta)) = U(T(\eta) T(\theta))` as follows

.. math::
	:nowrap:

	\begin{align*}
	U(T(\eta)) U(T(\theta)) &= \left( 1 + \i \eta^a u_a + \tfrac{1}{2} \eta^a \eta^b u_{ab} + \cdots \right) \left( 1 + \i \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots \right) \\
		&= 1 + \i (\eta^a + \theta^a) u_a \blue{- \eta^a \theta^b u_a u_b} + \cdots \\
	\\
	U(T(\eta) T(\theta)) &= U \left( 1 + \eta + \theta + f_{ab} \eta^a \theta^b + \cdots \right) \\
		&= 1 + \blue{\i} \left( \eta^c + \theta^c + \blue{f^c_{ab} \eta^a \theta^b} + \cdots \right) \blue{u_c} + \blue{\tfrac{1}{2}} \left( \blue{\eta^a + \theta^a} + \cdots \right) \left( \blue{\eta^b + \theta^b} + \cdots \right) \blue{u_{ab}} + \cdots
	\end{align*}

where :math:`f^c_{ab}` are the coefficients of the expansion of :math:`T(f(\eta, \theta)) = T(\eta) T(\theta)`. Equating the coefficients of :math:`\eta^a \theta^b`, i.e., the terms colored in blue, we get

.. math::
	:nowrap:

	\begin{equation*}
		-u_a u_b = \i f^c_{ab} u_c + u_{ab} \quad \Longrightarrow \quad u_{ab} = -u_a u_b - \i f^c_{ab} u_c.
	\end{equation*}

It implies that one can calculate the higher-order operator :math:`u_{ab}` from the lower-order ones, assuming of course that we know the structure of the symmetry (Lie) group/algebra. In fact, this bootstrapping procedure can be continued to all orders, but we'll not be bothered about the details.

Next, note that :math:`u_{ab} = u_{ba}` since they are just partial derivatives. It follows that

.. math::
	:nowrap:

	\begin{equation*}
		[u_a, u_b] \coloneqq u_a u_b - u_b u_a = \i (f^c_{ba} - f^c_{ab}) u_c \eqqcolon \i \Gamma^c_{ab} u_c
	\end{equation*}

where the bracket is known as the *Lie bracket* and :math:`\Gamma^c_{ab}` are known as the *structure constants*.

We conclude the general discussion about continuous symmetry by considering a special, but important, case when :math:`T` is additive in the sense that :math:`T(\eta) T(\theta) = T(\eta + \theta)`. Notable examples of such symmetry include translations and rotations about a fixed axis. In this case :math:`f` vanishes, and it follows from :math:`\eqref{eq_u_expansion}` that

.. math::
	:nowrap:

	\begin{equation}
		U(T(\theta)) = \lim_{N \to \infty} (U(T(\theta / N)))^N = \lim_{N \to \infty} (1 + \i \theta^a u_a / N)^N = \op{exp}(\i \theta^a u_a)
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

		\begin{align*}
			\pbf &= m d\xbf / d\tau = \frac{m \vbf}{\sqrt{1 - \vbf^2}} = m \vbf + O(|\vbf|^3) \\
			E &= m dt / d\tau = m + \tfrac{1}{2} m \vbf^2 + O(|\vbf|^4)
		\end{align*}

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

Taking determinant on :math:`\eqref{eq_homogeneous_lorentz_transformation}` implies that :math:`\op{det}(\Lambda) = \pm 1`, and setting :math:`\mu = \nu = 0` implies that

.. math::
	:nowrap:

	\begin{equation*}
		1 = (\Lambda^0_0)^2 - \Lambda^0_i \Lambda^0_i ~\Longrightarrow~ \left| \Lambda^0_0 \right| \geq 1
	\end{equation*}

It follows that the homogeneous Lorentz group has four components. In particular, the one with :math:`\op{det}(\Lambda) = 1` and :math:`\Lambda^0_0 \geq 1` is the most common used and is given a name: *proper orthochronous* Lorentz group. Nonetheless, one can map one component to another by composing with either a time reversal transformation

.. math::
	:nowrap:

	\begin{equation*}
		\Tcal: (t, \xbf) \mapsto (-t, \xbf)
	\end{equation*}

or a space reversal transformation

.. math::
	:nowrap:

	\begin{equation*}
		\Pcal: (t, \xbf) \mapsto (t, -\xbf)
	\end{equation*}

or both.

So far everything have been rather abstract, but in fact, the (homogeneous) Lorentz group can be understood quite intuitively. There are basically two building blocks: one is a rotation in the :math:`3`-space, which says that the space is homogeneous in all (spatial) directions, and the other is a so-called *boost*, which says that, as G. Galileo originally noted, one cannot tell if a system is stationary or is moving in a constant velocity without making a reference to outside of the system. To spell out the details, let's consider a stationary frame with :math:`d\xbf = 0` and a moving frame with :math:`d\xbf' / dt' = \vbf`. Then the transformation :math:`dx' = \Lambda dx` can be simplified as

.. math::
	:nowrap:

	\begin{equation*}
		dt' = \Lambda^0_0 dt, \quad dx'_i = \Lambda^0_i dt ~\Longrightarrow~ \Lambda^0_i = v_i \Lambda^0_0
	\end{equation*}

Then using :math:`\eqref{eq_homogeneous_lorentz_transformation}`, we get

.. math::
	:nowrap:

	\begin{equation*}
		1 = -\eta^{\mu \nu} \Lambda^0_{\mu} \Lambda^0_{\nu} = (\Lambda^0_0)^2 - \Lambda^0_i \Lambda^0_i = (1 - \vbf^2) (\Lambda^0_0)^2 \
			~\Longrightarrow~ \Lambda^0_0 = \frac{1}{\sqrt{1 - \vbf^2}}
	\end{equation*}

assuming :math:`\Lambda` is proper orthochronous. It's customary to write :math:`\gamma \coloneqq 1 / \sqrt{1 - \vbf^2}` so that

.. math::
	:nowrap:

	\begin{equation}
		\Lambda^0_0 = \gamma, \quad \Lambda^0_i = \gamma v_i
		\label{eq_lambda_boost}
	\end{equation}

The other components :math:`\Lambda^i_j, 1 \leq i, j \leq 3`, are not uniquely determined because a composition with a (spatial) rotation about the direction of :math:`\vbf` has no effect on :math:`\vbf`. To make it easier, one can apply a rotation so that :math:`\vbf` aligns with the :math:`3`-axis. Then an obvious choice of :math:`\Lambda` is given by

.. math::
	:nowrap:

	\begin{equation}
		\Lambda = \begin{bmatrix}
				\gamma & 0 & 0 & -\gamma v_3 \\
				0 & 1 & 0 & 0 \\
				0 & 0 & 1 & 0 \\
				\gamma v_3 & 0 & 0 & \gamma
			\end{bmatrix}
		\label{eq_lambda_in_3_axis}
	\end{equation}

which can be spelled out to the following possibly more intuitive expressions

.. math::
	:nowrap:

	\begin{align*}
		t' &= \gamma (t - v_3 x_3), \\
		x'_1 &= x_1, \\
		x'_2 &= x_2, \\
		x'_3 &= \gamma (x_3 + v_3 t).
	\end{align*}

It's now obvious that

1. a stationary point in the original frame is moving with constant velocity :math:`\vbf` (which is :math:`v_3` in this special case) in the prime-frame as planned, and
2. time appears to flow slower in the moving frame because :math:`\gamma < 1`, and
3. length appears to be shorter in the moving frame, again because :math:`\gamma < 1`.

Finally, one can apply a rotation to :math:`\eqref{eq_lambda_in_3_axis}` to get the general formulae

.. math::
	:nowrap:

	\begin{equation*}
		\Lambda_{ij} = \delta_{ij} + \frac{v_i v_j}{\vbf^2} (\gamma - 1)
	\end{equation*}

for :math:`1 \leq i, j \leq 3`, which together with :math:`\eqref{eq_lambda_boost}` gives the general formula for :math:`\Lambda`.

.. note::
	Any Lorentz transformation can be written as the composition of a boost followed by a rotation.

Quantum Lorentz symmetry
++++++++++++++++++++++++

We will quantize the Lorentz symmetry :math:`L(\Lambda, a)` by looking for unitarity representations :math:`U(\Lambda, a)`. As discussed in the section of :ref:`Continuous symmetry <sec_continuous_symmetry>`, we proceed by looking for infinitesimal symmetries. First of all, let's expand :math:`\Lambda` as

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

Comparing the first order terms shows that :math:`\omega^{\mu \nu} = -\omega^{\nu \mu}` is anti-symmetric. It is therefore more convenient to use :math:`\omega^{\mu \nu}`, rather than :math:`\omega_{\mu}^{\nu}`, as the infinitesimal parameters of the expansion of :math:`\Lambda`.

.. note::
	A count of free parameters shows that the inhomogeneous Lorentz symmetry has :math:`10` degrees of freedom, :math:`4` of which come from the translation, and the rest :math:`6` come from the rank-:math:`2` anti-symmetric tensor :math:`\omega`.

We first postulate that :math:`U(1, 0) = 1` is the identity operator because the Lorentz transformation itself is the identity. Then we can write the power series expansion up to first order as follows

.. math::
	:nowrap:

	\begin{equation}
		U(1 + \omega, \epsilon) = 1 - \i \epsilon^{\mu} P_{\mu} + \tfrac{\i}{2} \omega^{\mu \nu} J_{\mu \nu} + \cdots
		\label{eq_u_lorentz_expansion}
	\end{equation}

Here we have inserted :math:`\i` as usual so that the unitarity of :math:`U` implies that both :math:`P_{\mu}` and :math:`J_{\mu \nu}` are Hermitian. Moreover, since :math:`\omega^{\mu \nu}` is anti-symmetric, we can assume the same holds for :math:`J_{\mu \nu}`.

Let's evaluate how the expansion transformations under conjugation

.. math::
	:nowrap:

	\begin{align*}
		U(\Lambda, a) U(1 + \omega, \epsilon) U^{-1}(\Lambda, a) \
			&= U(\Lambda, a) U(1 + \omega, \epsilon) U(\Lambda^{-1}, -\Lambda^{-1} a) \\
			&= U(\Lambda, a) U((1 + \omega) \Lambda^{-1}, \epsilon - (1 + \omega) \Lambda^{-1} a) \\
			&= U(1 + \Lambda \omega \Lambda^{-1}, \Lambda \epsilon - \Lambda \omega \Lambda^{-1} a) \\
			&= 1 - \i (\Lambda^{\rho}_{\mu} \epsilon^{\mu} - \Lambda^{\rho}_{\mu} \omega^{\mu \nu} \Lambda_{\nu}^{\kappa} a_{\kappa}) P_{\rho} \
				+ \tfrac{\i}{2} \Lambda^{\rho}_{\mu} \omega^{\mu \nu} \Lambda_{\nu}^{\kappa} J_{\rho \kappa} + \cdots
	\end{align*}


where we have used the identity :math:`(\Lambda^{-1})_{\mu \nu} = \Lambda_{\mu \nu}`. Substituting :math:`U(1 + \omega, \epsilon)` with the expansion :math:`\eqref{eq_u_lorentz_expansion}` and equating the coefficients of :math:`\epsilon^{\mu}` and :math:`\omega^{\mu \nu}`, we have

.. math::
	:nowrap:

	\begin{align}
		U(\Lambda, a) P_{\mu} U^{-1}(\Lambda, a) &= \Lambda^{\rho}_{\mu} P_{\rho}  \label{eq_p_conjugated_by_u} \\
		U(\Lambda, a) J_{\mu \nu} U^{-1}(\Lambda, a) &= \Lambda^{\rho}_{\mu} \Lambda^{\kappa}_{\nu} (J_{\rho \kappa} + a_{\kappa} P_{\rho} - a_{\rho} P_{\kappa})  \label{eq_j_conjugated_by_u}
	\end{align}

where in the second equation, we have also made the right-hand-side anti-symmetric with respect to :math:`\mu` and :math:`\nu`. It's now clear that :math:`P` transforms like a vector and is translation invariant, while :math:`J` transforms like a :math:`2`-tensor only for homogeneous Lorentz transformations and is not translation invariant in general. These are of course as expected since both :math:`P` and :math:`J` are quantization of rather familiar objects, which we now spell out.

We start with :math:`P` by writing :math:`H \coloneqq P_0` and :math:`\Pbf \coloneqq (P_1, P_2, P_3)`. Then :math:`H` is the energy operator, also know as the *Hamiltonian*, and :math:`\Pbf` is the momentum :math:`3`-vector. Similarly, let's write :math:`\Kbf \coloneqq (J_{01}, J_{02}, J_{03})` and :math:`\Jbf = (J_{23}, J_{31}, J_{12})`, as the *boost* :math:`3`-vector and the *angular momentum* :math:`3`-vector, respectively.

Now that we have named all the players (i.e., :math:`H, \Pbf, \Jbf, \Kbf`) in the game, it remains to find out their mutual commutation relations since they should form a Lie algebra of the (infinitesimal) Lorentz symmetry. This can be done by applying :math:`\eqref{eq_p_conjugated_by_u}` and :math:`\eqref{eq_j_conjugated_by_u}` to :math:`U(\Lambda, a)` that is itself infinitesimal. More precisely, keeping up to first order terms, write :math:`\Lambda_{\mu}^{\nu} = \delta_{\mu}^{\nu} + \omega_{\mu}^{\nu}` and :math:`a^{\mu} = \epsilon^{\mu}` so that :math:`\eqref{eq_p_conjugated_by_u}` becomes

.. math::
	:nowrap:

	\begin{align*}
		\left( \delta_{\mu}^{\rho} + \omega_{\mu}^{\rho} \right) P_{\rho} &= \left( 1 - \i \epsilon^{\nu} P_{\nu} + \tfrac{\i}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) P_{\mu} \left( 1 + \i \epsilon^{\nu} P_{\nu} - \tfrac{\i}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
			&= P_{\mu} - \i \epsilon^{\nu} [P_{\mu}, P_{\nu}] - \tfrac{\i}{2} \omega^{\rho \kappa} [P_{\mu}, J_{\rho \kappa}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` and :math:`\omega` gives the following

.. math::
	:nowrap:

	\begin{align}
		[P_{\mu}, P_{\nu}] &= 0  \label{eq_bracket_p4_p4} \\
		[P_{\mu}, J_{\rho \kappa}] &= -\i (\eta_{\mu \rho} P_{\kappa} - \eta_{\mu \kappa} P_{\rho})  \label{eq_bracket_p4_j4}
	\end{align}

where we've used the identity :math:`\omega_{\mu}^{\rho} P_{\rho} = \eta_{\mu \kappa} \omega^{\rho \kappa} P_{\rho} = \tfrac{1}{2} \omega^{\rho \kappa} (\eta_{\mu \kappa} P_{\rho} - \eta_{\mu \rho} P_{\kappa})`. Now :math:`\eqref{eq_j_conjugated_by_u}` (up to first order) becomes

.. math::
	:nowrap:

	\begin{align*}
		J_{\mu \nu} + \epsilon_{\nu} P_{\mu} - \epsilon_{\mu} P_{\nu} + \omega_{\mu}^{\rho} J_{\rho \nu} + \omega_{\nu}^{\kappa} J_{\mu \kappa} &= (\delta_{\mu}^{\rho} + \omega_{\mu}^{\rho}) (\delta_{\nu}^{\kappa} + \omega_{\nu}^{\kappa}) (J_{\rho \kappa} + \epsilon_{\kappa} P_{\rho} - \epsilon_{\rho} P_{\kappa}) \\
		&= \left( 1 - \i \epsilon^{\rho} P_{\rho} + \tfrac{\i}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) J_{\mu \nu} \left( 1 + \i \epsilon^{\rho} P_{\rho} - \tfrac{\i}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
		&= J_{\mu \nu} - \i \epsilon^{\rho} [P_{\rho}, J_{\mu \nu}] + \tfrac{\i}{2} \omega^{\rho \kappa} [J_{\rho \kappa}, J_{\mu \nu}]
	\end{align*}

Equating the coefficients of :math:`\epsilon` reproduces :math:`\eqref{eq_bracket_p4_j4}`, but equating the coefficients of :math:`\omega` gives the following additional

.. math::
	:nowrap:

	\begin{equation}
		[J_{\rho \kappa}, J_{\mu \nu}] = -\i (\eta_{\kappa \mu} J_{\rho \nu} - \eta_{\rho \mu} J_{\kappa \nu} + \eta_{\nu \rho} J_{\mu \kappa} - \eta_{\nu \kappa} J_{\mu \rho})
		\label{eq_bracket_j4_j4}
	\end{equation}

Now that we have all the commutator relations, let's reorganize :math:`\eqref{eq_bracket_p4_p4}, \eqref{eq_bracket_p4_j4}, \eqref{eq_bracket_j4_j4}` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:nowrap:

	\begin{alignat}{2}
		\text{let } \mu = 0, \nu = i \text{ in \eqref{eq_bracket_p4_p4}} ~&\Longrightarrow~ [H, P_i] &&= 0  \label{eq_hp_commute} \\
		\text{let } \mu = 0, \rho = j, \kappa = k \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [H, J_i] &&= 0  \label{eq_hj_commute} \\
		\text{let } \mu = 0, \rho = 0, \kappa = i \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [H, K_i] &&= \i P_i  \nonumber \\
		\text{let } \mu = i, \nu = j \text{ in \eqref{eq_bracket_p4_p4}} ~&\Longrightarrow~ [P_i, P_j] &&= 0  \label{eq_pp_commute} \\
		\text{let } \mu = i, \rho = k, \kappa = i \text{ in \eqref{eq_bracket_p4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [P_i, J_j] &&= \i \epsilon_{ijk} P_k  \nonumber \\
		\text{let } \mu = i, \rho = 0 \text{ and enumerate } \kappa \in \{1, 2, 3\} \text{ in \eqref{eq_bracket_p4_j4}} ~&\Longrightarrow~ [P_i, K_j] &&= \i \delta_{ij} H  \nonumber \\
		\text{let } \rho = j, \kappa = \mu = k, \nu = i \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [J_i, J_j] &&= \i \epsilon_{ijk} J_k  \label{eq_jjj_commutation} \\
		\text{let } \rho = \nu = j, \kappa = k, \mu = 0 \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [J_i, K_j] &&= -\i \epsilon_{ijk} K_k  \nonumber \\
		\text{let } \rho = \mu = 0, \kappa = i, \nu = j \text{ in \eqref{eq_bracket_j4_j4} and permutation (anti-)symmetry} ~&\Longrightarrow~ [K_i, K_j] &&= -\i \epsilon_{ijk} J_k  \nonumber
 	\end{alignat}

where :math:`\epsilon_{ijk}` is totally anti-symmetric with respect to permutations of indexes and satisfies :math:`\epsilon_{123} = 1`.

These are some rather tedious and error-prone calculations. But in the end, we seem to at least get the important pieces right, namely :math:`\eqref{eq_hp_commute}, \eqref{eq_hj_commute}, \eqref{eq_pp_commute}` and :math:`\eqref{eq_jjj_commutation}`. Since the time evolution of a physical system is dictated by the Hamiltonian :math:`H`, quantities (i.e., observables) that commute with :math:`H` are conserved. In particular :math:`\eqref{eq_hp_commute}` and :math:`\eqref{eq_hj_commute}` imply that both momentum and angular momentum are conserved. Boosts, on the other hand, are *not* conserved, and therefore cannot be used to label (stable) physical states. Moreover :math:`\eqref{eq_pp_commute}` implies that translations commute with each other (as expected), which is *not* the case for the angular momenta according to :math:`\eqref{eq_jjj_commutation}`. Indeed, they furnish an infinitesimal representation of the :math:`3`-rotation group :math:`SO(3)`.

One-Particle States
-------------------

One neat application of our knowledge about Lorentz symmetry is to classify (free) one-particle states according to their transformation laws under (inhomogeneous) Lorentz transformations.

In order to do so, we need some labels to identify states, which are typically conserved quantities. According to the commutation relations between :math:`H, \Pbf` and :math:`\Jbf` obtained in the previous section, we see that :math:`p = (H, \Pbf)` consists of mutually commutative conserved components, but not :math:`\Jbf`. Hence we can write our one-particle states as :math:`\Psi_{p, \sigma}` such that

.. math::
	:nowrap:

	\begin{equation*}
		P_{\mu} \Psi_{p, \sigma} = p_{\mu}
	\end{equation*}

where :math:`\sigma` are additional labels such as spin components that we will later specify.

Let's first consider translations :math:`U(1, a)`. Since translations form an abelian group, it follows from :math:`\eqref{eq_additive_symmetry}` that

.. math::
	:nowrap:

	\begin{equation*}
		U(1, a) \Psi_{p, \sigma} = \op{exp}^{\i a^{\mu} P_{\mu}} \Psi_{p, \sigma} = \op{exp}^{\i a^{\mu} p_{\mu}} \Psi_{p, \sigma}
	\end{equation*}

Hence it remains to consider the action of homogeneous Lorentz transformations. For the convenience of notation, let's write :math:`U(\Lambda) \coloneqq U(\Lambda, 0)`. We would first like to know how :math:`U(\Lambda)` affects the :math:`4`-momentum. It follows from the following calculation

.. math::
	:nowrap:

	\begin{equation*}
		P_{\mu} U(\Lambda) \Psi_{p, \sigma} = U(\Lambda) (U^{-1} (\Lambda) P_{\mu} U(\Lambda)) \Psi_{p, \sigma} \
		\xlongequal{\eqref{eq_p_conjugated_by_u}} U(\Lambda) \Lambda_{\mu}^{\nu} P_{\nu} \Psi_{p, \sigma} = (\Lambda_{\mu}^{\nu} p_{\nu}) U(\Lambda) \Psi_{p, \sigma}
	\end{equation*}

that :math:`U(\Lambda) \Psi_{p, \sigma}` has :math:`4`-momentum :math:`\Lambda p`. Therefore we can write

.. math::
	:nowrap:

	\begin{equation*}
		U(\Lambda) \Psi_{p, \sigma} = C_{\sigma}^{\sigma'} (\Lambda, p) \Psi_{\Lambda p, \sigma'}
	\end{equation*}

where :math:`C(\Lambda, p)` furnishes, for each fixed :math:`p`, a representation of :math:`\Lambda`.

Next we'd like to remove the dependency of :math:`C` on :math:`p` since, after all, it is :math:`\Lambda` that carries the symmetry. We can achieve this by noticing that :math:`U(\Lambda)` acts on the :math:`p`-orbits of :math:`\Lambda` transitively. The :math:`p`-orbits of :math:`\Lambda`, in turn, is uniquely determined by the value of :math:`p^2`.