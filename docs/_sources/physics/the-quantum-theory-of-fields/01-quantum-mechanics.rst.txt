Quantum Mechanics
=================

The goal of this chapter it to cover the basics of quantum mechanics, i.e., the quantum theory of particles.

What is a Quantum State?
------------------------

Quantum theory postulates that *any* physical state (of the world) can be represented by a *ray* in some complex Hilbert space. It's worth noting that it is the state, rather than the Hilbert space, that we actually care about. Let's write a state as

.. math:: [\Psi] \coloneqq \{ c\Psi ~|~ c \in \Cbb \setminus 0 \}

where :math:`\Psi` is a nonzero vector in the Hilbert space. It is, however, rather inconvenient to have to deal with :math:`[\Psi]` all the time. So instead, we will almost always pick a representative :math:`\Psi`, often out of a natural choice, and call it a *state vector*, and keep in mind that anything physically meaningful must not be sensitive to a scalar multiplication.

.. admonition:: Assumption
	:class: Important

	Throughout this post we always assume that state vectors are normalized so that :math:`||\Psi|| = 1`.

In fact, we don't really care about the states themselves either, because they are more of an abstraction rather than something one can physically measure. What we do care about are the (Hermitian) inner products between state vectors, denoted by :math:`(\Psi, \Phi)`. According to the so-called `Copenhagen interpretation <https://en.wikipedia.org/wiki/Copenhagen_interpretation>`_ of quantum mechanics, such inner product represents an *amplitude*, i.e., its squared norm gives the probability of finding a state :math:`[\Psi]` in :math:`[\Phi]` if we ever perform a measurement. We can write this statement as an equation as follows

.. math:: P([\Psi] \to [\Phi]) = |(\Psi, \Phi)|^2

In particular, the probability of finding any state in itself is one, due to the normalization above.

.. _sec_what_is_a_symmetry:

What is a Symmetry?
-------------------

We start with a *symmetry transformation*, by which we mean a transformation that preserves all quantities that one can ever measure about a system. Since it is the probabilities, rather than the states themselves, that are measurable, one is led to define a quantum symmetry transformation as a transformation of states :math:`T` such that

.. math::
	:label: eq_t_preserves_probability

	P(T[\Psi] \to T[\Phi]) = P([\Psi] \to [\Phi])


for any states :math:`[\Psi]` and :math:`[\Phi]`. Now a theorem of E. Wigner asserts that such :math:`T` can be realized either as a linear unitary or as an anti-linear anti-unitary transformation :math:`U = U(T)` of state vectors in the sense that :math:`[U\Psi] = T[\Psi]` for any :math:`\Psi`. In other words, :math:`U` satisfies either

.. math:: U(c\Psi) = cU\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)

or

.. math:: U(c\Psi) = c^{\ast} U\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)^{\ast}

where :math:`c` is any complex number.

.. dropdown:: Proof of Wigner's theorem
	:icon: unlock
	:animate: fade-in-slide-down

	The construction of a realization :math:`U` of :math:`T` takes the following steps.

	:Step 1: Fix an orthonormal basis :math:`\Psi_i, i \in \Nbb`, of the Hilbert space.

	:Step 2:  For each :math:`\Psi_i`, choose a unit vector :math:`\Psi'_i` such that :math:`P[\Psi_i] = [\Psi'_i]`. Then :math:`\Psi'_i, i \in \Nbb`, also form an orthonormal basis by :eq:`eq_t_preserves_probability`. We'd like to define :math:`U` by asking

		.. math:: U\Psi_i = \Psi'_i

		for all :math:`i`, and extend by (anti-)linearity. But this is not going to realize :math:`T` in general because we haven't fixed the extra degrees of freedom -- the phases of :math:`\Psi'_i`.

	:Step 3: We fix the phases of :math:`\Psi'_k, k \geq 1`, relative to :math:`\Psi'_0`, by asking

		.. math:: T[\Psi_0 + \Psi_k] = [\Psi'_0 + \Psi'_k]

		To see why this is possible, note first that :math:`T[\Psi_0 + \Psi_k] = [\alpha \Psi'_0 + \beta \Psi'_k]`, where :math:`\alpha, \beta` are phase factors, due to :eq:`eq_t_preserves_probability` and the basis being orthonormal. Now :math:`[\alpha \Psi'_0 + \beta \Psi'_k] = [\Psi'_0 + (\beta/\alpha) \Psi'_k]` and we can absorb the phase :math:`\beta/\alpha` into the definition of :math:`\Psi'_k`. This is indeed the best one can do, because the last one degree of freedom, which is to multiply all :math:`\Psi'_i` by a phase, cannot be fixed.

	:Step 4: We have so far specified the value of :math:`U` on all of :math:`\Psi_i, i \geq 0`, and :math:`\Psi_0 + \Psi_k, k \geq 1`. Notice that all the coefficients of :math:`\Psi` are real. It is therefore instructive to ask what :math:`\Psi_0 + \ifrak \Psi_1` should be. By the same argument as in the previous step, we can write

		.. math:: T[\Psi_0 + \ifrak \Psi_1] = [\Psi'_0 + c \Psi'_1]

		where :math:`c` is a phase. Let's apply :eq:`eq_t_preserves_probability` once again as follows

		.. math::

			\sqrt{2} &= \left| \left( [\Psi_0 + \ifrak \Psi_1], [\Psi_0 + \Psi_1] \right) \right| \\
				&= \left| \left( T[\Psi_0 + \ifrak \Psi_1], T[\Psi_0 + \Psi_1] \right) \right| \\
				&= \left| \left( [\Psi'_0 + c \Psi'_1], [\Psi'_0 + \Psi'_1] \right) \right| \\
				&= |1 + c|

		It follows that :math:`c = \pm \ifrak`, which correspond to :math:`U` being (complex) linear or anti-linear, respectively.

	At this point, we can extend :math:`U` to either a linear or anti-linear map of the Hilbert space. But we'll not be bothered about any further formal argument, including showing that (anti-)linearity must be coupled with (anti-)unitarity, respectively.

.. note::
	The *adjoint* of a linear operator :math:`A` is another linear operator :math:`A^{\dagger}` such that

	.. math:: (\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)

	for all any two state vectors :math:`\Psi` and :math:`\Phi`. On the other hand, the adjoint of an anti-linear :math:`A` is another anti-linear :math:`A^{\dagger}` such that

	.. math:: (\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)^{\ast}

	A (anti-)unitary operator :math:`U` thus satisfies :math:`U^{\dagger} = U^{-1}`.

In general we're not interested in just one symmetry transformation, but rather a group -- whether continuous or discrete -- of symmetry transformations, or just symmetry for short. In particular, if :math:`T_1, T_2` are two symmetry transformations, then we'd like :math:`T_2 T_1` to also be a symmetry transformation. In light of the :math:`U`-realization of symmetry transformations discussed above, we can rephrase this condition as

.. math::
	:label: eq_u_depends_on_psi

	U(T_2 T_1) \Psi = \exp(\ifrak \theta(T_1, T_2, \Psi)) U(T_2) U(T_1) \Psi


where :math:`\ifrak = \sqrt{-1}`, and :math:`\theta(T_1, T_2, \Psi)` is an angle, which depends a priori on :math:`T_1, T_2`, and :math:`\Psi`.

It turns out, however, the angle :math:`\theta(T_1, T_2, \Psi)` cannot depend on the state because if we apply :eq:`eq_u_depends_on_psi` to the sum of two linearly independent state vectors :math:`\Psi_A + \Psi_B`, then we'll find

.. math::

	\exp(\pm \ifrak \theta(\Psi_A)) \Psi_A + \exp(\pm \ifrak \theta(\Psi_B)) \Psi_B
		= \exp(\pm \ifrak \theta(\Psi_A + \Psi_B)) (\Psi_A + \Psi_B)

where we have suppressed the dependency of :math:`\theta` on :math:`T`, and the signs correspond to the cases of :math:`U` being linear or anti-linear, respectively. In any case, it follows that

.. math::

	\exp(\pm \ifrak \theta(\Psi_A)) = \exp(\pm \ifrak \theta(\Psi_B))
		= \exp(\pm \ifrak \theta(\Psi_A + \Psi_B))

which says nothing but the independence of :math:`\theta` on :math:`\Psi`.

.. todo::
	While the argument here appears to be purely mathematical, Weinberg pointed out in [Wei95]_ (page 53) the potential inabilities to create a state like :math:`\Psi_A + \Psi_B`. More precisely, he mentioned the general believe that it's impossible to prepare a superposition of two states, one with integer total angular momentum and the other with half-integer total angular momentum, in which case there will be a "super-selection rule" between different classes of states. After all, one Hilbert space may just not be enough to describe all states. It'd be nice to elaborate a bit more on the super-selection rules.

We can now simplify :eq:`eq_u_depends_on_psi` to the following

.. math:: U(T_2 T_1) = \exp(\ifrak \theta(T_1, T_2)) U(T_2) U(T_1)

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
	:label: eq_u_expansion

	U(T(\theta)) = 1 + \ifrak \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots

where :math:`\theta^a` are the (real) components of :math:`\theta`, and :math:`u_a` are operators independent of :math:`\theta`, and as a convention, repeated indexes are summed up. Here we put a :math:`\ifrak` in front of the linear term so that the unitarity of :math:`U` implies that :math:`u_a` are Hermitian.

.. note::
	We've implicitly used a summation convention in writing :eq:`eq_u_expansion` that the same upper and lower indexes are automatically summed up. For example

	.. math:: \theta^a \theta^b u_{ab} \equiv \sum_{a, b} \theta^a \theta^b u_{ab}

	This convention will be used throughout this note, unless otherwise specified.

	Another noteworthy point is how one writes matrix or tensor elements using indexes. The point is that the indexes must come in certain order. This wouldn't really cause a problem if all indexes are lower or upper. However, care must be taken when both lower and upper indexes appear. For example, an element written as :math:`M^a_b` would be ambiguous as it's unclear whether it refers to :math:`M_{ab}` or :math:`M_{ba}` assuming that one can somehow raise/lower the indexes. To avoid such ambiguity, one writes either :math:`{M^a}_b` or :math:`{M_b}^a`.

	This is a particularly convenient convention when dealing with matrix or tensor multiplications. For example, one can multiply two matrices as follows

	.. math:: {M^a}_b {N^b}_c = {(MN)^a}_c

	while :math:`{M^a}_b {N_c}^b`, though still summed up over :math:`b`, wouldn't correspond to a matrix multiplication.


Now let :math:`\eta` be another element of the Lie algebra, and expand both sides of the equality :math:`U(T(\eta)) U(T(\theta)) = U(T(\eta) T(\theta))` as follows

.. math::

	U(T(\eta)) U(T(\theta))
		&= \left( 1 + \ifrak \eta^a u_a + \tfrac{1}{2} \eta^a \eta^b u_{ab} + \cdots \right) \left( 1 + \ifrak \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots \right) \\
		&= 1 + \ifrak (\eta^a + \theta^a) u_a \blue{- \eta^a \theta^b u_a u_b} + \cdots \\
	U(T(\eta) T(\theta))
		&= U \left( 1 + \eta + \theta + f_{ab} \eta^a \theta^b + \cdots \right) \\
		&= 1 + \blue{\ifrak} \left( \eta^c + \theta^c + \blue{{f^c}_{ab} \eta^a \theta^b} + \cdots \right) \blue{u_c} \\
		&\phantom{=} + \blue{\tfrac{1}{2}} \left( \blue{\eta^a + \theta^a} + \cdots \right) \left( \blue{\eta^b + \theta^b} + \cdots \right) \blue{u_{ab}} + \cdots

where :math:`{f^c}_{ab}` are the coefficients of the expansion of :math:`T(f(\eta, \theta)) = T(\eta) T(\theta)`. Equating the coefficients of :math:`\eta^a \theta^b`, i.e., the terms colored in blue, we get

.. math::

	-u_a u_b = \ifrak {f^c}_{ab} u_c + u_{ab} \implies u_{ab} = -u_a u_b - \ifrak {f^c}_{ab} u_c.

It implies that one can calculate the higher-order operator :math:`u_{ab}` from the lower-order ones, assuming of course that we know the structure of the symmetry (Lie) group/algebra. In fact, this bootstrapping procedure can be continued to all orders, but we'll not be bothered about the details.

Next, note that :math:`u_{ab} = u_{ba}` since they are just partial derivatives. It follows that

.. math::

	[u_a, u_b] \coloneqq u_a u_b - u_b u_a = \ifrak ({f^c}_{ba} - {f^c}_{ab}) u_c \eqqcolon \ifrak {C^c}_{ab} u_c

where the bracket is known as the *Lie bracket* and :math:`{C^c}_{ab}` are known as the *structure constants*.

We conclude the general discussion about continuous symmetry by considering a special, but important, case when :math:`T` is additive in the sense that :math:`T(\eta) T(\theta) = T(\eta + \theta)`. Notable examples of such symmetry include translations and rotations about a fixed axis. In this case :math:`f` vanishes, and it follows from :eq:`eq_u_expansion` that

.. math::
	:label: eq_additive_symmetry

	U(T(\theta)) = \lim_{N \to \infty} (U(T(\theta / N)))^N = \lim_{N \to \infty} (1 + \ifrak \theta^a u_a / N)^N = \exp(\ifrak \theta^a u_a)

.. _sec_lorentz_symmetry:

Lorentz symmetry
^^^^^^^^^^^^^^^^

A particularly prominent continuous symmetry in our physical world is the Lorentz symmetry postulated by Einstein's special relativity, which supersedes the Galilean symmetry, which is respected by the Newtonian mechanics. We shall start from the classical theory of Lorentz symmetry, and then quantize it following the procedure discussed in the previous section.

Classical Lorentz symmetry
++++++++++++++++++++++++++

Classical Lorentz symmetry is a symmetry that acts on the (flat) spacetime and preserves the so-called *proper time*

.. math::
	:label: eq_proper_time

	d\tau^2 \coloneqq dx_0^2 - dx_1^2 - dx_2^2 - dx_3^2 \eqqcolon -\eta^{\mu \nu} dx_{\mu} dx_{\nu}

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

	Using the notations introduced above, we can rewrite :eq:`eq_proper_time` as :math:`d\tau^2 = dt^2 - d\xbf^2`, so that it's obvious that if a particle travels at the speed of light in one inertial frame, i.e., :math:`|d\xbf / dt| = 1`, and equivalently :math:`d\tau = 0`, then it travels at the speed of light in any other inertial frame, in direct contradiction with Newtonian mechanics.

	Instead of working with the spacetime coordinates, it can sometimes be convenient to work with the "dual" energy-momentum coordinates, also known as the *four momentum*. The transition can be done by imagining a particle of mass :math:`m`, and defining :math:`p = (E, \pbf) \coloneqq m dx / d\tau`. It follows from :eq:`eq_proper_time` that

	.. math::
		:label: eq_four_momentum_mass_identity

		1 = (dt / d\tau)^2 - (d\xbf / d\tau)^2 \implies m^2 = (m dt / d\tau)^2 - (m d\xbf / d\tau)^2 = E^2 - \pbf^2

	which looks just like :eq:`eq_proper_time`, and indeed, the mass (in our convention) is invariant in all inertial frames.

	One can also recover Newtonian mechanics at the low-speed limit (i.e., :math:`|\vbf| \ll 1`) using :math:`d\tau / dt = \sqrt{1 - \vbf^2}` as follows

	.. math::
		:label: eq_p_from_v

		\begin{alignat*}{2}
			\pbf &= m d\xbf / d\tau &&= \frac{m \vbf}{\sqrt{1 - \vbf^2}} = m \vbf + O(|\vbf|^3) \\
			E &= m dt / d\tau &&= m + \tfrac{1}{2} m \vbf^2 + O(|\vbf|^4)
		\end{alignat*}

More precisely, by a Lorentz transformation we mean an inhomogeneous linear transformation

.. math:: L(\Lambda, a)x \coloneqq \Lambda x + a

which consists of a homogeneous part :math:`\Lambda` and a translation by :math:`a`. The proper time is obviously preserved by any translation, and also by :math:`\Lambda` if

.. math::
	:label: eq_homogeneous_lorentz_transformation

	\eta^{\mu \nu} dx_{\mu} dx_{\nu} = \eta^{\mu \nu} {\Lambda_{\mu}}^{\rho} {\Lambda_{\nu}}^{\kappa} dx_{\rho} dx_{\kappa} \
		\implies \eta^{\mu \nu} = \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu}

for any :math:`\mu` and :math:`\nu`. Moreover the group law is given by

.. math::

	L(\Lambda', a') L(\Lambda, a) x = L(\Lambda', a')(\Lambda x + a) = \Lambda' \Lambda x + \Lambda' a + a' = L(\Lambda' \Lambda, \Lambda' a + a') x

For later use, let's also calculate the inverse matrix of :math:`\Lambda` using :eq:`eq_homogeneous_lorentz_transformation` as follows

.. math::
	:label: eq_lambda_inverse

	\delta_{\sigma}^{\nu}
		= \eta_{\sigma \mu} \eta^{\mu \nu} = \eta_{\sigma \mu} \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu}
		\implies {(\Lambda^{-1})_{\sigma}}^{\kappa} = \eta_{\sigma\mu} \eta^{\rho\kappa} {\Lambda_{\rho}}^{\mu} = {\Lambda^{\kappa}}_{\sigma}

Now we'll take a look at the topology of the group of homogeneous Lorentz transformations. Taking determinant on both sides of :eq:`eq_homogeneous_lorentz_transformation`, we see that :math:`\op{det}(\Lambda) = \pm 1`. Moreover, setting :math:`\mu = \nu = 0`, we have

.. math::

	1 = \left( {\Lambda_0}^0 \right)^2 - \sum_{i=1}^3 \left( {\Lambda_i}^0 \right) \implies \left| {\Lambda_0}^0 \right| \geq 1

It follows that the homogeneous Lorentz group has four components. In particular, the one with :math:`\op{det}(\Lambda) = 1` and :math:`{\Lambda_0}^0 \geq 1` is the most common used and is given a name: *proper orthochronous* Lorentz group. Nonetheless, one can map one component to another by composing with either a time reversal transformation

.. math::
	:label: eq_time_inversion

	\Tcal: (t, \xbf) \mapsto (-t, \xbf)

or a space reversal transformation

.. math::
	:label: eq_space_inversion

	\Pcal: (t, \xbf) \mapsto (t, -\xbf)

or both.

So far everything have been rather abstract, but in fact, the (homogeneous) Lorentz group can be understood quite intuitively. There are basically two building blocks: one is a rotation in the :math:`3`-space, which says that the space is homogeneous in all (spatial) directions, and the other is a so-called *boost*, which says that, as G. Galileo originally noted, one cannot tell if a system is at rest or is moving in a constant velocity without making a reference to outside of the system. To spell out the details, let's consider a rest frame with :math:`d\xbf = 0` and a moving frame with :math:`d\xbf' / dt' = \vbf`. Then the transformation :math:`dx' = \Lambda dx` can be simplified as

.. math:: dt' = {\Lambda_0}^0 dt, \quad dx'_i = {\Lambda_i}^0 dt \implies {\Lambda_i}^0 = v_i {\Lambda_0}^0

Then using :eq:`eq_homogeneous_lorentz_transformation`, we get

.. math::
	:label: eq_def_gamma

	1 &= -\eta^{\mu \nu} {\Lambda_{\mu}}^0 {\Lambda_{\nu}}^0 \\
		&= \left( {\Lambda_0}^0 \right)^2 - \left( {\Lambda_i}^0 \right)^2 \\
		&= \left( 1 - \vbf^2 \right) \left( {\Lambda_0}^0 \right)^2
			\implies {\Lambda_0}^0 = \frac{1}{\sqrt{1 - \vbf^2}} \eqqcolon \gamma

assuming :math:`\Lambda` is proper orthochronous. It follows that

.. math::
	:label: eq_lambda_boost

	{\Lambda_i}^0 = -{\Lambda^0}_i = \gamma v_i

The other components :math:`{\Lambda_i}^j, 1 \leq i, j \leq 3`, are not uniquely determined because a composition with a (spatial) rotation about the direction of :math:`\vbf` has no effect on :math:`\vbf`. To make it easier, one can apply a rotation so that :math:`\vbf` aligns with the :math:`3`-axis. Then an obvious choice of :math:`\Lambda` is given by

.. math::
	:label: eq_lambda_in_3_axis

	\begin{alignat*}{2}
		t'   &= {\Lambda_0}^{\mu} x_{\mu} &&= \gamma (t + v_3 x_3) \\
		x'_1 &= {\Lambda_1}^{\mu} x_{\mu} &&= x_1 \\
		x'_2 &= {\Lambda_2}^{\mu} x_{\mu} &&= x_2 \\
		x'_3 &= {\Lambda_3}^{\mu} x_{\mu} &&= \gamma (x_3 + v_3 t)
	\end{alignat*}

.. dropdown:: Time dilation and length contraction
	:icon: unlock
	:animate: fade-in-slide-down

	A few consequences can be drawn from the boost transformation, most notably the effects of `time dilation <https://en.wikipedia.org/wiki/Time_dilation>`__ and `length contraction <https://en.wikipedia.org/wiki/Length_contraction>`__. The time dilation, i.e., a clock ticks slower in a moving frame than in a rest frame, is quite obvious from :eq:`eq_lambda_boost` and the fact that :math:`\gamma > 1`. But the length contraction requires some elaboration.

	To be more concrete, let's consider a rode of some fixed length. To measure the length, the measurement must be done *simultaneously* at the two ends of the rod. This constraint causes not much trouble in a rest frame, but must be taken care of in a moving frame since being simultaneous is not a Lorentz invariant property. Let :math:`x = (t, \xbf)` and :math:`y = (t', \ybf)` be the two endpoints of the rod in the rest frame, so that the length is :math:`|\xbf - \ybf|` regardless of whether :math:`t` and :math:`t'` are the same or not. Under the Lorentz transformation defined by :eq:`eq_lambda_in_3_axis`, they become

	.. math::

		\Lambda x &= (\gamma (t + v_3 x_3), x_1, x_2, \gamma(x_3 + v_3 t))  \\
		\Lambda y &= (\gamma (t' + v_3 y_3), y_1, y_2, \gamma(y_3 + v_3 t'))

	respectively. Setting the equal-time condition :math:`(\Lambda x)_0 = (\Lambda y)_0` gives :math:`t' = t + v_3 (x_3 - y_3)`. Substituting it into :math:`(\Lambda x)_3` and :math:`(\Lambda y)_3` then gives

	.. math::

		|(\Lambda x)_3 - (\Lambda y)_3| = \gamma \left| x_3 - y_3 - v_3^2 (x_3 - y_3) \right| = \frac{|x_3 - y_3|}{\gamma} < |x_3 - y_3|

	This calculation says that the length of rod is contracted in the direction of movement.

	It should be emphasized that such contraction of length can only be observed in a frame where the rod is moving. Imagine for example a scenario where you're given a square box with equal sides while standing still, then after some unconscious period of time, e.g., sleeping, you wake up with the same box in hand, and you'd like to know if you're now moving. If you happen to have heard of such contraction of length, you might try to measure the sides of the box again. If one of the sides suddenly becomes shorter, then you know not only that you're moving, but also the direction of movement! This is of course absurd because the box is still at rest relative to you.

Finally, one can apply a rotation to :eq:`eq_lambda_in_3_axis` to get the general formula

.. math::
	:label: eq_general_lambda_in_spacetime

	\Lambda_{ij} = \delta_{ij} + \frac{v_i v_j}{\vbf^2} (\gamma - 1)

for :math:`1 \leq i, j \leq 3`, which, together with :eq:`eq_lambda_boost` and :math:`{\Lambda_0}^i = {\Lambda_i}^0,` gives the general formula for :math:`\Lambda`.

.. note::
	Any Lorentz transformation can be written as the composition of a boost followed by a rotation.

.. _sec_quantum_lorentz_symmetry:

Quantum Lorentz symmetry
++++++++++++++++++++++++

We will quantize the Lorentz symmetry :math:`L(\Lambda, a)` by looking for unitarity representations :math:`U(\Lambda, a)`. As discussed in :ref:`sec_continuous_symmetry`, we proceed by looking for infinitesimal symmetries. First of all, let's expand :math:`\Lambda` as

.. math::
	:label: eq_expansion_of_Lambda

	{\Lambda_{\mu}}^{\nu} = {\delta_{\mu}}^{\nu} + {\omega_{\mu}}^{\nu} + \cdots

where :math:`\delta` is the Kronecker delta, and *not* a tensor. It follows from :math:`\eta^{\mu \nu} = \eta^{\rho \kappa} {\Lambda_{\rho}}^{\mu} {\Lambda_{\kappa}}^{\nu}` that

.. math::
	:label: eq_lorentz_lie_algebra_is_antisymmetric

	\eta^{\mu \nu} &= \eta^{\rho \kappa} ({\delta_{\rho}}^{\mu} + {\omega_{\rho}}^{\mu} + \cdots) ({\delta_{\kappa}}^{\nu} + {\omega_{\kappa}}^{\nu} + \cdots) \\
		&= \eta^{\mu \nu} + \eta^{\mu \kappa} {\omega_{\kappa}}^{\nu} + \eta^{\nu \rho} {\omega_{\rho}}^{\mu} + \cdots \\
		&= \eta^{\mu \nu} + \omega^{\mu \nu} + \omega^{\nu \mu} + \cdots

Comparing the first order terms shows that :math:`\omega^{\mu \nu} = -\omega^{\nu \mu}` is anti-symmetric. It is therefore more convenient to use :math:`\omega^{\mu \nu}`, rather than :math:`\omega_{\mu}^{\nu}`, as the infinitesimal parameters in the expansion of :math:`\Lambda`.

.. note::
	A count of free parameters shows that the inhomogeneous Lorentz symmetry has :math:`10` degrees of freedom, :math:`4` of which come from the translation, and the rest :math:`6` come from the rank-:math:`2` anti-symmetric tensor :math:`\omega`.

We first postulate that :math:`U(1, 0) = 1` is the identity operator because the Lorentz transformation itself is the identity. Then we can write the power series expansion up to first order as follows

.. math::
	:label: eq_u_lorentz_expansion

	U(1 + \omega, \epsilon) = 1 - \ifrak \epsilon^{\mu} P_{\mu} + \frac{\ifrak}{2} \omega^{\mu \nu} J_{\mu \nu} + \cdots

Here we have inserted :math:`\ifrak` as usual so that the unitarity of :math:`U` implies that both :math:`P_{\mu}` and :math:`J_{\mu \nu}` are
Hermitian. Moreover, since :math:`\omega^{\mu \nu}` is anti-symmetric, we can assume the same holds for :math:`J_{\mu \nu}`.

.. note::
	Since we are expanding :math:`U(1 + \epsilon)` which is complex linear, the operators :math:`P` and :math:`J` are also complex linear. Hence we can freely move :math:`\ifrak` around these operators in calculations that follow. However, this will become an issue when we later consider other operators such as the space and time inversions, which can potentially be either complex linear or anti-linear. In the later case, a sign needs to be added when commuting with the multiplication by :math:`\ifrak`.

Let's evaluate how the expansion transformations under conjugation

.. math::

	U(\Lambda, a) U(1 + \omega, \epsilon) U^{-1}(\Lambda, a)
		&= U(\Lambda, a) U(1 + \omega, \epsilon) U(\Lambda^{-1}, -\Lambda^{-1} a) \\
		&= U(\Lambda, a) U((1 + \omega) \Lambda^{-1}, \epsilon - (1 + \omega) \Lambda^{-1} a) \\
		&= U(1 + \Lambda \omega \Lambda^{-1}, \Lambda \epsilon - \Lambda \omega \Lambda^{-1} a) \\
		&= 1 - \ifrak ({\Lambda^{\rho}}_{\mu} \epsilon^{\mu} - (\Lambda \omega \Lambda^{-1})^{\rho \kappa} a_{\kappa}) P_{\rho} \
			+ \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} J_{\rho \kappa} + \cdots \\
		&= 1 -\ifrak \epsilon^{\mu} {\Lambda^{\rho}}_{\mu} P_{\rho} + \tfrac{\ifrak}{2} (\Lambda \omega \Lambda^{-1})^{\rho \kappa} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots \\
		&= 1 -\ifrak \epsilon^{\mu} {\Lambda^{\rho}}_{\mu} P_{\rho} + \tfrac{\ifrak}{2} {\Lambda^{\rho}}_{\mu} \omega^{\mu \nu} {\Lambda^{\kappa}}_{\nu} (J_{\rho \kappa} + 2a_{\kappa} P_{\rho}) + \cdots

where we have used :eq:`eq_lambda_inverse` for :math:`\Lambda^{-1}`. Substituting :math:`U(1 + \omega, \epsilon)` with the expansion :eq:`eq_u_lorentz_expansion` and equating the coefficients of :math:`\epsilon^{\mu}` and :math:`\omega_{\mu \nu}`, we have

.. math::
	:label: eq_p_and_j_conjugated_by_u

	U(\Lambda, a) P_{\mu} U^{-1}(\Lambda, a) &= {\Lambda^{\rho}}_{\mu} P_{\rho} \\
	U(\Lambda, a) J_{\mu \nu} U^{-1}(\Lambda, a) &= {\Lambda^{\rho}}_{\mu} {\Lambda^{\kappa}}_{\nu} (J_{\rho \kappa} + a_{\kappa} P_{\rho} - a_{\rho} P_{\kappa})

where in the second equation, we have also made the right-hand-side anti-symmetric with respect to :math:`\mu` and :math:`\nu`. It's now clear that :math:`P` transforms like a vector and is translation invariant, while :math:`J` transforms like a :math:`2`-tensor only for homogeneous Lorentz transformations and is not translation invariant in general. These are of course as expected since both :math:`P` and :math:`J` are quantization of rather familiar objects, which we now spell out.

We start with :math:`P` by writing :math:`H \coloneqq P_0` and :math:`\Pbf \coloneqq (P_1, P_2, P_3)`. Then :math:`H` is the energy operator, also know as the *Hamiltonian*, and :math:`\Pbf` is the momentum :math:`3`-vector. Similarly, let's write :math:`\Kbf \coloneqq (J_{01}, J_{02}, J_{03})` and :math:`\Jbf = (J_{23}, J_{31}, J_{12})`, as the *boost* :math:`3`-vector and the *angular momentum* :math:`3`-vector, respectively.

Now that we have named all the players (i.e., :math:`H, \Pbf, \Jbf, \Kbf`) in the game, it remains to find out their mutual commutation relations since they should form a Lie algebra of the (infinitesimal) Lorentz symmetry. This can be done by applying :eq:`eq_p_and_j_conjugated_by_u` to :math:`U(\Lambda, a)` that is itself infinitesimal. More precisely, keeping up to first order terms, we have :math:`{\Lambda^{\rho}}_{\mu} = {\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu}` and :math:`a_{\mu} = \epsilon_{\mu}`. It follows that :eq:`eq_p_and_j_conjugated_by_u`, up to first order, can be written as follows

.. math::

	\left( {\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu} \right) P_{\rho}
		&= \left( 1 - \ifrak \epsilon^{\nu} P_{\nu} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) P_{\mu} \left( 1 + \ifrak \epsilon^{\nu} P_{\nu} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
		&= P_{\mu} - \ifrak \epsilon^{\nu} [P_{\mu}, P_{\nu}] - \tfrac{\ifrak}{2} \omega^{\rho \kappa} [P_{\mu}, J_{\rho \kappa}]

Equating the coefficients of :math:`\epsilon` and :math:`\omega` gives the following

.. math::
	:label: eq_bracket_pp_and_pj

	[P_{\mu}, P_{\nu}] &= 0  \label{eq_bracket_p4_p4} \\
	[P_{\mu}, J_{\rho \kappa}] &= \ifrak (\eta_{\kappa \mu} P_{\rho} - \eta_{\rho \mu} P_{\kappa})

where for the second identity, we've also used the fact that :math:`J_{\rho \kappa} = -J_{\kappa \rho}`.

Similarly, expanding :eq:`eq_p_and_j_conjugated_by_u` up to first order, we have

.. math::

	J_{\mu \nu} + \epsilon_{\nu} P_{\mu} - \epsilon_{\mu} P_{\nu} + {\omega^{\rho}}_{\mu} J_{\rho \nu} + {\omega^{\kappa}}_{\nu} J_{\mu \kappa}
		&= ({\delta^{\rho}}_{\mu} + {\omega^{\rho}}_{\mu}) ({\delta^{\kappa}}_{\nu} + {\omega^{\kappa}}_{\nu}) (J_{\rho \kappa} + \epsilon_{\kappa} P_{\rho} - \epsilon_{\rho} P_{\kappa}) \\
		&= \left( 1 - \ifrak \epsilon^{\rho} P_{\rho} + \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) J_{\mu \nu} \left( 1 + \ifrak \epsilon^{\rho} P_{\rho} - \tfrac{\ifrak}{2} \omega^{\rho \kappa} J_{\rho \kappa} \right) \\
		&= J_{\mu \nu} - \ifrak \epsilon^{\rho} [P_{\rho}, J_{\mu \nu}] + \tfrac{\ifrak}{2} \omega^{\rho \kappa} [J_{\rho \kappa}, J_{\mu \nu}]

Equating the coefficients of :math:`\epsilon` reproduces :eq:`eq_bracket_pp_and_pj`, but equating the coefficients of :math:`\omega` gives the following additional

.. math::
	:label: eq_bracket_j4_j4

	[J_{\rho \kappa}, J_{\mu \nu}] = \ifrak (\eta_{\rho \nu} J_{\mu \kappa} - \eta_{\rho \mu} J_{\nu \kappa} - \eta_{\kappa \mu} J_{\rho \nu} - \eta_{\kappa \nu} J_{\rho \mu})

Now that we have all the commutator relations, let's reorganize :eq:`eq_bracket_pp_and_pj` and :eq:`eq_bracket_j4_j4` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:label: eq_poincare_algebra

	[H, P_i] &= 0 \\
	[H, J_i] &= 0 \\
	[H, K_i] &= \ifrak P_i \\
	[P_i, P_j] &= 0 \\
	[P_i, J_j] &= \ifrak \epsilon_{ijk} P_k \\
	[P_i, K_j] &= \ifrak \delta_{ij} H \\
	[J_i, J_j] &= \ifrak \epsilon_{ijk} J_k \\
	[J_i, K_j] &= \ifrak \epsilon_{ijk} K_k \\
	[K_i, K_j] &= -\ifrak \epsilon_{ijk} J_k

where :math:`\epsilon_{ijk}` is totally anti-symmetric with respect to permutations of indexes and satisfies :math:`\epsilon_{123} = 1`. [#tedious_calc_of_commutations]_

.. note::
	The Lie algebra generated by :math:`H, \Pbf, \Jbf, \Kbf` with commutation relations :eq:`eq_poincare_algebra` is known as the `Poincaré algebra <https://en.wikipedia.org/wiki/Poincar%C3%A9_group>`__.

Since the time evolution of a physical system is dictated by the Hamiltonian :math:`H`, quantities that commute with :math:`H` are conserved. In particular we see from :eq:`eq_poincare_algebra` that both momentum :math:`\Pbf` and angular momentum :math:`\Jbf` are conserved. Boosts :math:`\Kbf`, on the other hand, are *not* conserved, and therefore cannot be used to label (stable) physical states. Moreover, momenta (which generate translations) commute with each other, while angular momenta (which generate rotation) do not, indeed, they furnish an infinitesimal representation of the :math:`3`-rotation group :math:`SO(3)`. This should be all consistent with our intuition.

.. _sec_one_particle_states:

One-Particle States
-------------------

One neat application of our knowledge about Lorentz symmetry is to classify (free) one-particle states according to their transformation laws under (inhomogeneous) Lorentz transformations. Throughout this section, the Lorentz transformations will be assumed to be proper orthochronous, i.e., :math:`\op{det}(\Lambda) = 1` and :math:`{\Lambda_0}^0 \geq 1`.

In order to do so, we need some labels to identify states, which are typically conserved quantities. According to the commutation relations between :math:`H, \Pbf` and :math:`\Jbf` obtained in the previous section, we see that :math:`p = (H, \Pbf)` consists of mutually commutative conserved components, but not :math:`\Jbf`. Hence we can write our one-particle states as :math:`\Psi_{p, \sigma}` such that

.. math:: P_{\mu} \Psi_{p, \sigma} = p_{\mu}

where :math:`\sigma` are additional labels such as spin components that we will later specify.

Reduction to the little group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's first consider translations :math:`U(1, a)`. Since translations form an abelian group, it follows from :eq:`eq_additive_symmetry` that

.. math::
	:label: eq_translation_formula_for_particle_state

	U(1, a) \Psi_{p, \sigma} = \exp(-\ifrak a^{\mu} P_{\mu}) \Psi_{p, \sigma} = \exp(-\ifrak a^{\mu} p_{\mu}) \Psi_{p, \sigma}

where the minus sign comes from our choice of expansion :eq:`eq_u_lorentz_expansion`. Hence it remains to consider the action of homogeneous Lorentz transformations. For the convenience of notation, let's write :math:`U(\Lambda) \coloneqq U(\Lambda, 0)`. We would first like to know how :math:`U(\Lambda)` affects the :math:`4`-momentum. It follows from the following calculation (using :eq:`eq_p_and_j_conjugated_by_u`)

.. math::

	P_{\mu} U(\Lambda) \Psi_{p, \sigma}
		= U(\Lambda) (U^{-1} (\Lambda) P_{\mu} U(\Lambda)) \Psi_{p, \sigma}
		= U(\Lambda) {\Lambda_{\mu}}^{\nu} P_{\nu} \Psi_{p, \sigma}
		= \left( {\Lambda_{\mu}}^{\nu} p_{\nu} \right) U(\Lambda) \Psi_{p, \sigma}

that :math:`U(\Lambda) \Psi_{p, \sigma}` has :math:`4`-momentum :math:`\Lambda p`. Therefore we can write

.. math::
	:label: eq_lorentz_acts_on_p_and_sigma

	U(\Lambda) \Psi_{p, \sigma} = C_{\sigma \sigma'} (\Lambda, p) \Psi_{\Lambda p, \sigma'}

where :math:`C_{\sigma \sigma'}` furnishes a representation of :math:`\Lambda` and :math:`p` under straightforward transformation rules, and an implicit summation over :math:`\sigma'` is assumed although it's not a :math:`4`-index.

Next we'd like to remove the dependency of :math:`C_{\sigma \sigma'}` on :math:`p` since, after all, it is :math:`\Lambda` that carries the symmetry. We can achieve this by noticing that :math:`U(\Lambda)` acts on the :math:`\Lambda`-orbits of :math:`p` transitively. The :math:`\Lambda`-orbits of :math:`p`, in turn, are uniquely determined by the value of :math:`p^2`, and in the case of :math:`p^2 \leq 0`, also by the sign of :math:`p_0`. In light of :eq:`eq_four_momentum_mass_identity`, we can pick a convenient representative :math:`k` for each case as follows

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

With the choices of the standard :math:`k` in hand, we need to make one more set of choices. Namely, we will choose for each :math:`p` a standard Lorentz transformation :math:`L(p)` such that :math:`L(p) k = p`. Such :math:`L(p)` for a massive particle has been chosen in :eq:`eq_general_lambda_in_spacetime`, albeit in spacetime coordinates, and we'll also handle the case of massless particles later. Once these choices have been made, we can *define*

.. math::
	:label: eq_def_of_one_particle_psi

	\Psi_{p, \sigma} \coloneqq N(p) U(L(p)) \Psi_{k, \sigma}

where :math:`N(p)` is a normalization factor to be determined later. In this way, we've also determined how :math:`\sigma` depends on :math:`p`. Applying :eq:`eq_lorentz_acts_on_p_and_sigma` to :eq:`eq_def_of_one_particle_psi` we can refactor the terms as follows

.. math::
	:label: eq_def_of_one_particle_psi_refactored

	U(\Lambda) \Psi_{p, \sigma}
		&= N(p) U(\Lambda) U(L(p)) \Psi_{k, \sigma} \\
		&= N(p) U(L(\Lambda p)) U(L(\Lambda p)^{-1} \Lambda L(p)) \Psi_{k, \sigma}

so that :math:`L(\Lambda p)^{-1} \Lambda L(p)` maps :math:`k` to itself, and hence :math:`U(L(\Lambda p)^{-1} \Lambda L(p))` acts solely on :math:`\sigma`.

At this point, we have reduced the problem to the classification of representations of the so-called *little group* defined as the subgroup of (proper orthochronous) Lorentz transformations :math:`W` that fixes :math:`k`, i.e., :math:`{W_{\mu}}^{\nu} k_{\nu} = k_{\mu}`. Element in the little group is known as `Wigner rotation <https://en.wikipedia.org/wiki/Wigner_rotation>`__ (and hence :math:`W`). More precisely, the task now is to find (unitary) representations :math:`D(W)` such that

.. math::

	\sum_{\sigma'} D_{\sigma \sigma'}(W_1) D_{\sigma' \sigma''}(W_2) \Psi_{k, \sigma''} = D_{\sigma \sigma''}(W_1 W_2) \Psi_{k, \sigma''}

Once this is done, we can define

.. math::
	:label: eq_d_repr_of_little_group

	U(W) \Psi_{k, \sigma} \coloneqq \sum_{\sigma'} D_{\sigma' \sigma}(W) \Psi_{k, \sigma'},  \quad\text{where}~~
		W(\Lambda, p) \coloneqq L(\Lambda p)^{-1} \Lambda L(p)


.. dropdown:: Validation of :eq:`eq_d_repr_of_little_group`
	:animate: fade-in-slide-down
	:icon: unlock

	One can verify that :eq:`eq_d_repr_of_little_group` indeed respects the group law as follows

	.. math::

		U(W_2) U(W_1) \Psi_{k, \sigma}
			&= U(W_2) \sum_{\sigma'} D_{\sigma' \sigma}(W_1) \Psi_{k, \sigma'} \\
			&= \sum_{\sigma' \sigma''} D_{\sigma' \sigma}(W_1) D_{\sigma'' \sigma'}(W_2) \Psi_{k, \sigma''} \\
			&= \sum_{\sigma''} D_{\sigma'' \sigma}(W_2 W_1) \Psi_{k, \sigma''}

Now we can rewrite :eq:`eq_def_of_one_particle_psi_refactored` (using :eq:`eq_def_of_one_particle_psi` and :eq:`eq_d_repr_of_little_group`) as follows

.. math::
	:label: eq_little_group_acts_on_p_and_sigma

	U(\Lambda) \Psi_{p, \sigma}
		&= N(p) U(L(\Lambda p)) U(W(\Lambda, p)) \Psi_{k, \sigma} \\
		&= N(p) \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) U(L(\Lambda p)) \Psi_{k, \sigma'} \\
		&= \frac{N(p)}{N(\Lambda p)} \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'}

which gives the sought-after coefficients :math:`C_{\sigma \sigma'}` in :eq:`eq_lorentz_acts_on_p_and_sigma`.

It remains now, as far as the general discussion is concerned, to settle the normalization factor :math:`N(p)`. Indeed, it'd not have been needed at all if we'd like :math:`\Psi_{p, \sigma}` be to orthonormal in the sense that

.. math::
	:label: eq_psi_p4_sigma_orthonormal

	(\Psi_{p', \sigma'}, \Psi_{p, \sigma}) = \delta_{\sigma' \sigma} \delta(p' - p)

where the first delta is the Kronecker delta (for discrete indexes) and the second is the Dirac delta (for continuous indexes), since they are eigenvectors of the (Hermitian) operator :math:`P`. All we need is :math:`D_{\sigma \sigma'}` being unitary as is obvious from :eq:`eq_little_group_acts_on_p_and_sigma`.

However, the Dirac delta in :eq:`eq_psi_p4_sigma_orthonormal` is tricky to use since :math:`p` is constrained to the so-called *mass shell*, i.e., :math:`p_0 > 0` together with :math:`p^2 = -M^2` in the massive case and :math:`p^2 = 0` in the massless case, respectively. Hence the actual normalization we'd like to impose on the one-particle states is, instead of :eq:`eq_psi_p4_sigma_orthonormal`, the following

.. math::
	:label: eq_psi_p3_sigma_orthonormal

	(\Psi_{p', \sigma'}, \Psi_{p, \sigma}) = \delta_{\sigma' \sigma} \delta(\pbf' - \pbf)

In fact, the problem eventually boils down to how to define the :math:`3`-momentum space Dirac delta in a Lorentz-invariant manner.

Since :math:`\Psi_{p, \sigma}` can be derived from :math:`\Psi_{k, \sigma}` by :eq:`eq_def_of_one_particle_psi`, we can first ask :math:`\Psi_{k, \sigma}` to be orthonormal in the sense of :eq:`eq_psi_p3_sigma_orthonormal`, where the Dirac delta plays no role, and then figure out how integration works on the mass shell (because Dirac delta is defined by integrals against test functions). As far as the mass shell integration is concerned, we can temporarily unify the massive and massless cases by allowing :math:`M \geq 0`. Consider a general mass shell integral of an arbitrary test function :math:`f(p)`

.. math::

	\int d^4 p ~\delta(p^2 + M^2) \theta(p_0) f(p)
		&= \int d^3\pbf dp_0 ~\delta(p_0^2 - \pbf^2 - M^2) \theta(p_0) f(p_0, \pbf) \\
		&= \int d^3\pbf ~\frac{f\left( \sqrt{\pbf^2 + M^2}, \pbf \right)}{2 \sqrt{\pbf^2 + M^2}}

where :math:`\theta(p_0)` is the step function defined to be :math:`0` if :math:`p_0 \leq 0` and :math:`1` if :math:`p_0 > 1`. It follows that the Lorentz-invariant volume element in the :math:`3`-momentum space is

.. math::
	:label: eq_lorentz_invariant_3_momentum_volume_element

	\frac{d^3\pbf}{\sqrt{\pbf^2 + M^2}}

We can use it to find the Lorentz-invariant Dirac delta (marked in blue) as follows

.. math::

	f(\pbf') &\eqqcolon \int d^3\pbf ~\delta(\pbf' - \pbf) f(\pbf) \\
		&= \int \frac{d^3\pbf}{\sqrt{\pbf^2 + M^2}} \blue{p_0 \delta(\pbf' - \pbf)} f(\pbf)

It follows from Lorentz invariance that :math:`p_0 \delta(\pbf' - \pbf) = k_0 \delta(\kbf' - \kbf)`. Hence we can finally establish :eq:`eq_psi_p3_sigma_orthonormal` as follows

.. math::

	(\Psi_{p', \sigma'}, \Psi_{p, \sigma})
		&= N(p) N(p')^{\ast} (U(L(p')) \Psi_{k', \sigma'}, U(L(p)) \Psi_{k, \sigma}) \\
		&= |N(p)|^2 \delta_{\sigma' \sigma} \delta(\kbf' - \kbf) \\
		&= \delta_{\sigma' \sigma} \delta(\pbf' - \pbf)

if we define :math:`N(p) = \sqrt{k_0 / p_0}`.

Putting everything together, we've obtained the following grand formula for the Lorentz transformation law

.. math::
	:label: eq_lorentz_transformation_formula_for_particle_state

	U(\Lambda) \Psi_{p, \sigma} = \sqrt{\frac{(\Lambda p)_0}{p_0}} \sum_{\sigma'} D_{\sigma' \sigma}(W(\Lambda, p)) \Psi_{\Lambda p, \sigma'}

where :math:`D_{\sigma' \sigma}` is a unitary representation of the little group, and :math:`W(\Lambda, p)` is defined by :eq:`eq_d_repr_of_little_group`.

Massive particle states
^^^^^^^^^^^^^^^^^^^^^^^

Recall the standard :math:`4`-momentum :math:`k = (M, 0, 0, 0)` in this case. Obviously the little group here is nothing but the :math:`3`-rotation group :math:`SO(3)`. We can work out :math:`D_{\sigma \sigma'}(\Rcal)` by a rotation :math:`\Rcal \in SO(3)` up to first order as follows.

First write :math:`\Rcal^{ij} = \delta^{ij} + \Theta^{ij}` such that :math:`\Theta` is anti-symmetric. Then expand :math:`D_{\sigma \sigma'} (\Rcal)` similar to :eq:`eq_u_lorentz_expansion` up to first order as follows

.. math::

	D_{\sigma \sigma'} (\Rcal) = \delta_{\sigma \sigma'} + \tfrac{\ifrak}{2} \Theta^{ij} (J_{ij})_{\sigma \sigma'}

where :math:`J_{ij}` is a collection of Hermitian operators that satisfy :math:`J_{ij} = -J_{ji}` and the commutation relations :eq:`eq_poincare_algebra`. It turns out that there exists an infinite number of such unitary representations indexed by nonnegative half-integers :math:`\jfrak = 0, \tfrac{1}{2}, 1, \tfrac{3}{2}, \cdots`, each of which has dimension :math:`2\jfrak + 1`. Choosing the :math:`3`-axis as the preferred axis of (definite) spin, we can summarize the result as follows

.. math::
	:label: eq_representation_rotation_first_order

	D^{(\jfrak)}_{\sigma \sigma'} (\Rcal) = \delta_{\sigma \sigma'} + \tfrac{\ifrak}{2} \Theta^{ij} \left( J^{(\jfrak)}_{ij} \right)_{\sigma \sigma'}

where :math:`J^{\jfrak}_{ij}` satisfy the following commutation relations

.. math::
	:label: eq_rotation_j_matrix

	\left( J^{(\jfrak)}_{23} \pm \ifrak J^{(\jfrak)}_{31} \right)_{\sigma \sigma'}
		\equiv \left( J^{(\jfrak)}_1 \pm \ifrak J^{(\jfrak)}_2 \right)_{\sigma \sigma'}
		&= \delta_{\sigma \pm 1, \sigma'} \sqrt{(\jfrak \mp \sigma)(\jfrak \pm \sigma + 1)} \\
	\left( J^{(\jfrak)}_{12} \right)_{\sigma \sigma'}
		\equiv \left( J^{(\jfrak)}_3 \right)_{\sigma \sigma'}
		&= \sigma \delta_{\sigma \sigma'}  \label{eq_j3_matrix}

where :math:`\sigma, \sigma'` run through the values :math:`-\jfrak, -\jfrak + 1, \cdots, \jfrak - 1, \jfrak`.

.. _dropdown_repr_of_angular_momenta:

.. dropdown:: Representations of angular momenta
	:icon: unlock
	:animate: fade-in-slide-down

	Recall from :eq:`eq_poincare_algebra` that the (quantum) angular momenta vector :math:`\Jbf` satisfy the commutation relations :math:`[J_i, J_j] = \ifrak \epsilon_{ijk} J_k`. Hence they cannot be simultaneously diagonalized. It's then a convention to use the angular momentum along the :math:`3`-axis to label the spin. The following two identities are straightforward but important

	.. math::

		[\Jbf^2, J_i] &= 0, ~\forall i = 1, 2, 3 \\
		[J_3, J_1 \pm \ifrak J_2] &= \pm (J_1 \pm \ifrak J_2)

	where :math:`\Jbf^2 = J_1^2 + J_2^2 + J_3^2` as usual.

	Now if :math:`\Psi_{\sigma}` is an eigenstate of :math:`J_3` with eigenvalue :math:`\sigma`, then

	.. math::
		:label: eq_j1_j2_raises_or_lowers_state

		J_3 (J_1 \pm \ifrak J_2) \Psi_{\sigma} = [J_3, J_1 \pm \ifrak J_2] \Psi_{\sigma} + (J_1 \pm \ifrak J_2) J_3 \Psi_{\sigma} = (\sigma \pm 1) \Psi_{\sigma}

	In other words, applying :math:`J_1 \pm \ifrak J_2` to any eigenstate of :math:`J_3` raises or lowers the eigenvalue by one, and henceforth they are called *raising* and *lowering* operators, respectively. Moreover, since :math:`\Jbf^2` commutes with :math:`J_3`, we may assume that :math:`\Psi_{\sigma}` is also an eigenstate of :math:`\Jbf^2`, and since :math:`\Jbf^2` also commutes with :math:`J_1 \pm \ifrak J_2`, the whole series of :math:`J_3`-eigenstates obtained by applying the raising/lowering operators have the same :math:`\Jbf^2`-eigenvalue.

	We'll from now on focus on eigenstates with a fixed :math:`\Jbf^2`-eigenvalue. Moreover we'd like the eigenvalues of :math:`J_3` to be bounded, so both the raising and the lowering operations must stop after finite steps. Let :math:`\Psi_{\jfrak}` be the :math:`J_3`-eigenstate with the highest eigenvalue (if there are more than one, the representation is reducible). By repeatedly applying the lowering operator to :math:`\Psi_{\jfrak}`, we'll eventually reach the eigenstate :math:`\Psi_{\jfrak'}` with the lowest eigenvalue. Since the lowering operator decreases the eigenvalue by one, we know that :math:`\jfrak - \jfrak'` must be an integer.

	Consider the following two operators

	.. math::
		:label: eq_j1_j2_mixed_product

		(J_1 - \ifrak J_2) (J_1 + \ifrak J_2) &= J_1^2 + J_2^2 + \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 - J_3 \\
		(J_1 + \ifrak J_2) (J_1 - \ifrak J_2) &= J_1^2 + J_2^2 - \ifrak [J_1, J_2] = \Jbf^2 - J_3^2 + J_3

	Note that the first operator annihilates :math:`\Psi_{\jfrak}` and the second operator annihilates :math:`\Psi_{\jfrak'}` by assumption, which, together with the fact that :math:`\Jbf^2 \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'}`, implies

	.. math::
		:label: eq_angular_momentum_squared_eigenvalue

		\Jbf^2 \Psi_{\jfrak} = (\jfrak^2 + \jfrak) \Psi_{\jfrak} = ((\jfrak')^2 - \jfrak') \Psi_{\jfrak} = \Jbf^2 \Psi_{\jfrak'} \implies \jfrak (\jfrak + 1) = \jfrak' (\jfrak' - 1)

	The equation has two potential solutions: either :math:`\jfrak' = \jfrak + 1` or :math:`\jfrak = -\jfrak'`. The first option violates the maximality of :math:`\jfrak`, and so we must accept the second option. Since we also know :math:`\jfrak - \jfrak'` must be integral, we conclude that :math:`\jfrak` is itself a half-integer.

	As a piece of notation, we'll from now on write :math:`\Psi^{\jfrak}_{\sigma}` for the eigenstate of both :math:`\Jbf^2` and :math:`J_3` such that

	.. math::
		:label: eq_3j_square_eigenstate

		\Jbf^2 \Psi^{\jfrak}_{\sigma} &= \jfrak (\jfrak + 1) \Psi^{\jfrak}_{\sigma} \\
		J_3 \Psi^{\jfrak}_{\sigma} &= \sigma \Psi^{\jfrak}_{\sigma}

	We'll also write :math:`J_i^{(\jfrak)} \coloneqq J_i` to explicitly indicate the dependency on :math:`\jfrak`.

	It remains to settle the constant term on the right-hand-side of :eq:`eq_rotation_j_matrix`. By :eq:`eq_j1_j2_raises_or_lowers_state` we can assume

	.. math:: \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} = \alpha_{\pm}(\jfrak, \sigma) \Psi^{\jfrak}_{\sigma \pm 1}

	Applying :eq:`eq_j1_j2_mixed_product` to :math:`\Psi^{\jfrak}_{\sigma}` then implies

	.. math:: \alpha_{\mp} (\jfrak, \sigma \pm 1) \alpha_{\pm} (\jfrak, \sigma) = \jfrak^2 + \jfrak - \sigma^2 \mp \sigma

	Now we use the fact that :math:`J^{(\jfrak)}_i, i = 1, 2, 3`, are Hermitian operators to calculate

	.. math::

		|\alpha_{\pm} (\jfrak, \sigma)|^2 (\Psi^{\jfrak}_{\sigma}, \Psi^{\jfrak}_{\sigma})
			&= \left( \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma}, \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} \right) \\
			&= \left( \Psi^{\jfrak}_{\sigma}, \left( J_1^{(\jfrak)} \mp \ifrak J_2^{(\jfrak)} \right) \left( J_1^{(\jfrak)} \pm \ifrak J_2^{(\jfrak)} \right) \Psi^{\jfrak}_{\sigma} \right) \\
			&= (\jfrak^2 + \jfrak - \sigma^2 \mp \sigma) (\Psi^{\jfrak}_{\sigma}, \Psi^{\jfrak}_{\sigma})

	It follows that, up to a choice of phase, :math:`\alpha_{\pm} (\jfrak, \sigma) = \sqrt{\jfrak^2 + \jfrak - \sigma^2 \mp \sigma} = \sqrt{(j \mp \sigma)(j \pm \sigma + 1)}`, which confirms :eq:`eq_rotation_j_matrix`.

We end the discussion about massive particle states by working out the little group elements :math:`W(\Lambda, p)` defined by :eq:`eq_d_repr_of_little_group`. To this end, it suffices to work out the standard :math:`L(p)` such that :math:`L(p) k = p`, where :math:`k = (M, 0, 0, 0)`. We have already worked out such a transformation in :eq:`eq_lambda_boost` and :eq:`eq_general_lambda_in_spacetime` in spacetime coordinates, so we only need to translate it into :math:`4`-momentum coordinates.

Using :eq:`eq_p_from_v`, we can rewrite :math:`\gamma` defined by :eq:`eq_def_gamma` as follows

.. math::

	\pbf = \frac{M \vbf}{\sqrt{1 - \vbf^2}} \implies \gamma \coloneqq \frac{1}{\sqrt{1 - \vbf^2}} = \frac{\sqrt{M^2 + \pbf^2}}{M} \left( = \frac{p_0}{M} \right)

It follows that

.. math::
	:label: eq_L_transformation_for_massive

	{L(p)_0}^0 = {L(p)^0}_0 &= \gamma \\
	{L(p)_i}^0 = -{L(p)^0}_i &= \frac{p_i}{M} \\
	L(p)_{ij} &= \delta_{ij} + \frac{p_i p_j}{\pbf^2} (\gamma - 1)

Finally, we note an important fact that when :math:`\Lambda = \Rcal` is a :math:`3`-rotation, then

.. math::
	:label: eq_little_group_rotation

	W(\Rcal, p) = \Rcal

for any :math:`p`. To see this, we'll work out how :math:`W(\Rcal, p)` acts on :math:`(1, \mathbf{0}), (0, \pbf)`, and :math:`(0, \qbf)`, respectively, where :math:`\qbf` is any :math:`3`-vector perpendicular to :math:`\pbf`, as follows

.. math::

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

	1. non-relativistic, which is justified by :eq:`eq_little_group_rotation`, and
	2. applicable mostly (but not necessarily) to multi-particles states, rather than single-particle states.

	We shall focus on the composition of two angular momentum :math:`3`-vectors :math:`\Jbf'` and :math:`\Jbf''`, whether orbital or spin, that commute, i.e., they each satisfies :eq:`eq_poincare_algebra` and in addition

	.. math:: [J'_i, J''_j] = 0

	for :math:`1 \leq i, j \leq 3`.

	Let's recollect the eigenstate representations :eq:`eq_rotation_j_matrix` and :eq:`eq_3j_square_eigenstate` as follows,

	.. math::

		{\Jbf'}^2 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \jfrak' (\jfrak' + 1) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
		J'_3 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sigma' \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
		(J'_1 \pm \ifrak J'_2) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sqrt{{\jfrak'}^2 + \jfrak' - {\sigma'}^2 \mp \sigma'} ~\Psi^{\jfrak' ~\jfrak''}_{\sigma' \pm 1, \sigma''} \\
		{\Jbf''}^2 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \jfrak'' (\jfrak'' + 1) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
		J''_3 \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sigma'' \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} \\
		(J''_1 \pm \ifrak J''_2) \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''} &= \sqrt{{\jfrak''}^2 + \jfrak'' - {\sigma''}^2 \mp \sigma''} ~\Psi^{\jfrak' ~\jfrak''}_{\sigma', \sigma'' \pm 1}

	Without knowing exactly how the Hamiltonian :math:`H` looks like, we cannot really say what combinations of these angular momentum operators commute with :math:`H`, and therefore may be used to label states. However, one can imagine that a rotationally invariant Hamiltonian may contain terms like :math:`{\Jbf'}^2, {\Jbf''}^2` and interactions like :math:`\Jbf' \cdot \Jbf''`. In this case, we may choose to consider the following collection of (mutually commuting) operators

	.. math:: {\Jbf'}^2, ~{\Jbf''}^2, ~\Jbf^2, \text{ and } J_3

	where :math:`\Jbf \coloneqq \Jbf' + \Jbf''` is the total angular momentum, and :math:`J_i, i=1,2,3,` is its :math:`i`-th component.

	Now our goal is to express the eigenstates :math:`\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma}` which satisfy the following

	.. math::
		:label: eq_raising_lowering_am_pair

		{\Jbf'}^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak' (\jfrak' + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \\
		{\Jbf''}^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak'' (\jfrak'' + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \\
		\Jbf^2 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \jfrak (\jfrak + 1) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \\
		J_3 \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \sigma \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} \\
		(J_1 \pm \ifrak J_2) \Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} &= \sqrt{\jfrak^2 + \jfrak - \sigma^2 \mp \sigma} ~\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma \pm 1}

	in terms of :math:`\Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''}` as follows

	.. math::
		:label: eq_defn_clebsch_gordan_coefficients

		\Psi^{\jfrak' ~\jfrak'' ~\jfrak}_{\sigma} = \sum_{\sigma' \sigma''} C^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \Psi^{\jfrak' ~\jfrak''}_{\sigma' \sigma''}

	where the coefficients are known as Clebsch-Gordan coefficients. For the clarity of exposition, let's divide the solution into a few steps.

	Step 1.
		First of all, note that since :math:`J_3 = J'_3 + J''_3`, we have the following constraint

		.. math::
			:label: eq_sigma_additive

			C^{\jfrak' ~\jfrak''}(\jfrak ~\sigma; \sigma' \sigma'') \neq 0 \implies \sigma = \sigma' + \sigma''

		Moreover, we see that the maximum possible value of :math:`\sigma` is :math:`\jfrak' + \jfrak''`, and it's achieved exactly when :math:`\sigma' = \jfrak'` and :math:`\sigma'' = \jfrak''`. It follows, assuming the non-degeneracy of the representation at least, that

		.. math::
			:label: eq_highest_weight_am_pair

			\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak''} =  \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''}


		or equivalently

		.. math::

			C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' ~\jfrak' + \jfrak''; \sigma' \sigma'') = \delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak''}

	Step 2.
		Next consider a state with :math:`\sigma = \jfrak' + \jfrak'' - 1`. It follows from :eq:`eq_sigma_additive` that it must be a superposition of :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''}` and :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1}` unless :math:`\jfrak'` and/or :math:`\jfrak''` vanishes, which leads to even simpler situations. Now we have two possible values of :math:`\jfrak`, namely :math:`\jfrak' + \jfrak''` and :math:`\jfrak' + \jfrak'' - 1`.

		In the former case, we can use :eq:`eq_raising_lowering_am_pair` with :math:`\sigma = \jfrak = \jfrak' + \jfrak''` and :eq:`eq_highest_weight_am_pair` to calculate as follows

		.. math::

			\sqrt{2(\jfrak' + \jfrak'')} ~\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak'' - 1} \
				&= (J_1 - \ifrak J_2) \Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak''} \\
				&= (J_1 - \ifrak J_2) \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''} \\
				&= (J'_1 - \ifrak J'_2 + J''_2 - \ifrak J''_2) \Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak''} \\
				&= \sqrt{2 \jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''} + \sqrt{2 \jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1}

		This gives us one of the :math:`\sigma = \jfrak' + \jfrak'' - 1` states

		.. math::
			:label: eq_second_highest_weight_am_pair_one

			\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak''}_{\jfrak' + \jfrak'' - 1}
				= (\jfrak' + \jfrak'')^{-1/2} \left( \sqrt{\jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak''} + \sqrt{\jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1} \right)

		The other one, which has :math:`\jfrak = \jfrak' + \jfrak'' - 1`, must be orthogonal to :eq:`eq_second_highest_weight_am_pair_one`. Therefore up to a normalization factor, we can write

		.. math::
			:label: eq_second_highest_weight_am_pair_two

			\Psi^{\jfrak' ~\jfrak'' ~\jfrak' + \jfrak'' - 1}_{\jfrak' + \jfrak'' - 1}
				= (\jfrak' + \jfrak'')^{-1/2} \left( \sqrt{\jfrak''} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak'-1 ~\jfrak''} - \sqrt{\jfrak'} ~\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 1} \right)

		We can translate :eq:`eq_second_highest_weight_am_pair_one` and :eq:`eq_second_highest_weight_am_pair_two` into Clebsch-Gordan coefficients as follows

		.. math::
			:label: eq_clebsch_gordan_second_highest_weight

			& C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' ~\jfrak' + \jfrak'' - 1; \sigma' \sigma'') \\
			& \quad = \sqrt{\frac{\jfrak'}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak' - 1} \delta_{\sigma'' ~\jfrak''} +
					\sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak'' - 1} \\
			& C^{\jfrak' ~\jfrak''}(\jfrak' + \jfrak'' -1 ~\jfrak' + \jfrak'' - 1; \sigma' \sigma'') \\
			& \quad = \sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} ~\delta_{\sigma' ~\jfrak' - 1} \delta_{\sigma'' ~\jfrak''} -
				\sqrt{\frac{\jfrak''}{\jfrak' + \jfrak''}} \delta_{\sigma' ~\jfrak'} \delta_{\sigma'' ~\jfrak'' - 1}

	Step 3.
		The pattern should now be clear. Namely for :math:`\sigma = \jfrak' + \jfrak'' - 2`, the :math:`J_3`-eigenspace must be :math:`3`-dimensional and spanned by :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 2 ~\jfrak''}, \Psi^{\jfrak' ~\jfrak''}_{\jfrak' - 1 ~\jfrak'' - 1}` and :math:`\Psi^{\jfrak' ~\jfrak''}_{\jfrak' ~\jfrak'' - 2}`. Two of them come from :eq:`eq_second_highest_weight_am_pair_one` and :eq:`eq_second_highest_weight_am_pair_two` by applying the lowering operator :math:`J_1 - \ifrak J_2`, and the third orthogonally complements the first two.

		This procedure can be continued but will eventually terminate because of the bounds :math:`|\sigma'| \leq \jfrak'` and :math:`|\sigma''| \leq \jfrak''`. It follows that :math:`\jfrak` can only take the following values

		.. math::
			:label: eq_composite_total_angular_momentum_range

			\jfrak = |\jfrak' - \jfrak''|, ~|\jfrak' - \jfrak''| + 1, \cdots, ~\jfrak' + \jfrak''

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

	\begin{alignat*}{2}
		(Wv)^{\mu} (Wv)_{\mu} &= v^{\mu} v_{\mu} &&= 1 \\
		(Wv)^{\mu} k_\mu &= v^{\mu} k_{\mu} &&= 1
	\end{alignat*}

It follows that we can write :math:`Wv = (1 + c, a, b, c)` with :math:`a^2 + b^2 = 2c`. Playing similar games to the other basis vectors, we can engineer a particular Lorentz transformation as follows

.. math::
	:label: eq_massless_little_group_s_matrix

	{S_{\mu}}^{\nu}(a, b) = \begin{bmatrix*}[r]
		1 + c & a & b & -c \\
		a & 1 & 0 & -a \\
		b & 0 & 1 & -b \\
		c & a & b & 1 - c
	\end{bmatrix*}

which leaves :math:`k` invariant, and satisfies :math:`Sv = Wv`. It follows that :math:`S^{-1} W` must be a rotation about the :math:`3`-axis, which can be written as follows

.. math::
	:label: eq_massless_little_group_r_matrix

	R(\theta) = \begin{bmatrix*}[r]
		1 & 0 & 0 & 0 \\
		0 & \cos\theta & \sin\theta & 0 \\
		0 & -\sin\theta & \cos\theta & 0 \\
		0 & 0 & 0 & 1
	\end{bmatrix*}

Hence we can write any element in the little group as :math:`W(a, b, \theta) = S(a, b) R(\theta)`.

.. dropdown:: The little group is :math:`~E^+(2)`
	:icon: unlock
	:animate: fade-in-slide-down

	Although not necessary for our purposes here, we'd like to better understand the little group for :math:`k = (1, 0, 0, 1)` in terms of more familiar groups. It turns out that it's isomorphic to the :math:`2`-dimensional orientation-preserving `Euclidean group <https://en.wikipedia.org/wiki/Euclidean_group>`__ :math:`E^+(2)`, i.e., the group of rotations and translations on the plane.

	To see this, we go back to the defining property of :math:`W` that it fixes :math:`k`. It follows that it must also fix the orthogonal complement :math:`k^{\bot}` with respect to the bilinear form :math:`d\tau^2` defined in :eq:`eq_proper_time`. Since :math:`k` is orthogonal to itself, we can uniquely determine :math:`W` by knowing its action on :math:`(1, 1, 0, 1)` and :math:`(1, 0, 1, 1)`. Letting :math:`S(a, b)` act on them, we see

	.. math::

		S(a, b)(1, 1, 0, 1) &= (1 + a, 1, 0, 1 + a) \\
		S(a, b)(1, 0, 1, 1) &= (1 + b, 0, 1, 1 + b)

	Hence :math:`S` is isomorphic to a :math:`2`-dimensional translation group. Moreover, the direction of translation is determined by the rotation on the plane spanned by the second the the third coordinates, which is nothing but :math:`R`.

	Note that :math:`E^+ (2)` is not semisimple in the sense that it possesses an abelian normal subgroup. Indeed, it's obvious from the above discussion that a translation conjugated by a rotation is again a translation (in the rotated direction). The non-semisimplicity will have consequences on the representation as we will see below.

As in the massive case, we'll work out :math:`D_{\sigma \sigma'}` up to first order. To this end, note that up to first order

.. math::

	{W(a, b, \theta)_{\mu}}^{\nu} &= \left(1 + \begin{bmatrix*}[r]
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

where we've added the :math:`4`-indexes since we recall from discussions in :ref:`sec_quantum_lorentz_symmetry` that we must lift the :math:`\omega` index to make it anti-symmetric. We now rewrite

.. math::

	W(a, b , \theta)^{\mu \nu} = \eta^{\mu \sigma} {W_{\sigma}}^{\nu} = 1 + \begin{bmatrix*}[r]
			0 & -a & -b & 0 \\
			a & 0 & \theta & -a \\
			b & -\theta & 0 & -b \\
			0 & a & b & 0
		\end{bmatrix*} + \cdots

and spell out the expansion of :math:`D(a, b, \theta) \coloneqq D(W(a, b, \theta))` as follows

.. math::
	:label: eq_massless_D_matrix_expansion

	D(a, b, \theta) = 1 + \ifrak aA + \ifrak bB + \ifrak \theta J_3

where

.. math::

	\begin{alignat*}{2}
		A &= -J_{01} - J_{13} &&= -K_1 + J_2 \\
		B &= -J_{02} - J_{23} &&= -K_2 - J_1
	\end{alignat*}

Next we use :eq:`eq_poincare_algebra` to calculate commutation relations between :math:`A, B` and :math:`J_3` as follows

.. math::

	\begin{alignat*}{2}
		[J_3, A] &= -&&\ifrak K_2 &&- \ifrak J_1 &&= \ifrak B \\
		[J_3, B] &= &&\ifrak K_1 &&- \ifrak J_2 &&= -\ifrak A \\
		[A, B] &= -&&\ifrak J_3 &&+ \ifrak J_3 &&= 0
	\end{alignat*}

Since :math:`A, B` commute, we can use their eigenvalues to label states as follows

.. math::

	A \Psi_{k, a, b} &= a \Psi_{k, a, b} \\
	B \Psi_{k, a, b} &= b \Psi_{k, a, b}

In fact, these states, corresponding to translation symmetries, come in continuous families as shown below

.. math::

	AU^{-1}(R(\theta)) \Psi_{a, b, k} &= (a\cos\theta - b\sin\theta)U^{-1}(R(\theta)) \Psi_{a, b, k} \\
	BU^{-1}(R(\theta)) \Psi_{a, b, k} &= (a\sin\theta + b\cos\theta)U^{-1}(R(\theta)) \Psi_{a, b, k}

According to [Wei95]_ (page 72), massless particle states are not observed to come in such :math:`S^1`-families. Hence the only possibility is that :math:`a = b = 0` and the only symmetry left then is :math:`J_3`, which corresponds to a rotation about the :math:`3`-axis.

Unlike the :math:`SO(3)`-symmetry discussed in :ref:`Representations of angular momenta <dropdown_repr_of_angular_momenta>`, representations of :math:`J_3` alone cannot be characterized at the infinitesimal level, which would have resulted in a continuous spectrum. Instead, since a :math:`2\pi`-rotation about the :math:`3`-axis gives the identity transformation, one might expect an integer spectrum for :math:`J_3`. This is indeed the case if we :ref:`assume the representation is genuine <assump_genuine_repr>`. However, since the Lorentz group is not simplify connected (with fundamental group :math:`\Zbb/2`), one may encounter projective representations. Indeed, the :math:`2\pi`-rotation about the :math:`3`-axis represents a generator of the fundamental group, which has order :math:`2`, i.e., only the :math:`4\pi`-rotation about the :math:`3`-axis represents a contractible loop in the Lorentz group (see the `Plate trick <https://en.wikipedia.org/wiki/Plate_trick>`_). As a result, the :math:`J_3`-spectrum actually consists of half-integers, just like the spins. We can therefore write a general massless particle state as :math:`\Psi_{k, \sigma}` such that

.. math::

	J_3 \Psi_{k, \sigma} = \sigma \Psi_{k, \sigma}

where :math:`\sigma` are half-integers, known as the *helicity*.

Combining the discussions so far, we can write down the :math:`D`-matrix defined by :eq:`eq_massless_D_matrix_expansion` as follows

.. math::
	:label: eq_little_group_d_matrix_massless

	D_{\sigma \sigma'}(W(a, b, \theta)) = \exp(\ifrak \theta \sigma) \delta_{\sigma \sigma'}

where we recall :math:`W(a, b, \theta) = L(\Lambda p)^{-1} \Lambda L(p) = S(a, b)R(\theta)`. The Lorentz transformation formula :eq:`eq_lorentz_transformation_formula_for_particle_state` for massless particles now becomes

.. math::
	:label: eq_lorentz_transformation_formula_for_massless

	U(\Lambda) \Psi_{p, \sigma} = \sqrt{\frac{(\Lambda p)_0}{p_0}} \exp(\ifrak \theta(\Lambda, p) \sigma) \Psi_{\Lambda p, \sigma}

In particular, we see that, unlike the spin :math:`z`-component of massive particles, helicity is Lorentz invariant (at least under genuine representations). It is reasonable, therefore, to think of massless particles of different helicity as different particle species. Examples include photons with :math:`\sigma = \pm 1` and gravitons with :math:`\sigma = \pm 2`, but *not* (anti-)neutrinos with hypothetical :math:`\sigma = \pm \tfrac{1}{2}` as otherwise stated in [Wei95]_ (page 73 -- 74), which are now known to have a nonzero mass. Here the :math:`\pm` signs are related to the space-inversion symmetry :eq:`eq_space_inversion`, which will be discussed in detail later.

In order to use :eq:`eq_little_group_d_matrix_massless` for a general :math:`(\Lambda, p)`, we first need to fix the choices of :math:`L(p)` that takes the standard :math:`k = (1, 0, 0, 1)` to :math:`p`. This can be done in two steps. First apply a (pure) boost along the :math:`3`-axis

.. math::
	:label: eq_massless_boost

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

Then apply a (pure) rotation that takes :math:`(0, 0, |\pbf|)` to :math:`\pbf`. However, in contrast to the massive case :eq:`eq_L_transformation_for_massive`, where :math:`L(p)` depends continuously on :math:`p`, there exists no continuous family of rotations that take :math:`(0, 0, |\pbf|)` to any other :math:`3`-vector (of the same length). Fortunately, any two choices of such rotations differ by (a pre-composition of) a rotation about the :math:`3`-axis, which, according to :eq:`eq_little_group_d_matrix_massless`, only produces a physically immaterial phase factor.

.. dropdown:: Polarization of photons
	:icon: unlock
	:animate: fade-in-slide-down

	General photon states of definite momentum can be written as a superposition

	.. math:: \Psi_{p, \alpha} \coloneqq \alpha_+ \Psi_{p, +1} + \alpha_- \Psi_{p, -1}

	such that :math:`|\alpha_+|^2 + |\alpha_-|^2 = 1`. They are not in general Lorentz invariant due to the angle :math:`\theta` presented in the phase factor in :eq:`eq_little_group_d_matrix_massless`. Here the base states :math:`\Psi_{p, \pm 1}` are known as (right and left-handedly) *circularly polarized*, and the other extreme cases where :math:`|\alpha_+| = |\alpha_-|` are known as *linearly polarized*. All other intermediate cases are then called *elliptically polarized*.

	It's not obvious at all why these states are named the way they are, if only viewed as abstract combinations of eigenstates of an abstract operator :math:`J_3`. They are named after analogies, either with `helicity <https://en.wikipedia.org/wiki/Helicity_(particle_physics)>`_ from classical mechanics or with classical electromagnetic fields (in vacuum) from Maxwell's theory.

.. _sec_space_and_time_inversions:

Space and time inversions
^^^^^^^^^^^^^^^^^^^^^^^^^

So far the discussions have been focused on orthochronous (and mostly homogeneous) Lorentz transformations, and in particular, the infinitesimal symmetries at the vicinity of the identity. Now it's time to take a look at the space and time inversions, defined in :eq:`eq_space_inversion` and :eq:`eq_time_inversion`, respectively, which takes us to the other components of the Lorentz group. The main goal is to understand their actions on the one-particle states, that have been worked out in the previous two sections.

Let's write

.. math:: U(\Pcal) \coloneqq U(\Pcal, 0), \quad U(\Tcal) \coloneqq U(\Tcal, 0)

for the corresponding quantum symmetry operators, which we haven't decided whether should be complex linear or anti-linear. The same calculations that led to :eq:`eq_p_and_j_conjugated_by_u` now give

.. math::
	:label: eq_p_and_j_conjugated_by_space_and_time_inversions

	U(\Pcal) \ifrak P_{\mu} U^{-1}(\Pcal) &= \ifrak {\Pcal_{\mu}}^{\rho} P_{\rho} \\
	U(\Pcal) \ifrak J_{\mu \nu} U^{-1}(\Pcal) &= \ifrak {\Pcal_{\mu}}^{\rho} {\Pcal_{\nu}}^{\kappa} J_{\rho \kappa} \\
	U(\Tcal) \ifrak P_{\mu} U^{-1}(\Tcal) &= \ifrak {\Tcal_{\mu}}^{\rho} P_{\rho} \\
	U(\Tcal) \ifrak J_{\mu \nu} U^{-1}(\Tcal) &= \ifrak {\Tcal_{\mu}}^{\rho} {\Tcal_{\nu}}^{\kappa} J_{\rho \kappa}

The complex (anti-)linearity of :math:`U(\Pcal)` and :math:`U(\Tcal)` can then be decided by the postulation that physically meaningful energy must not be negative. More precisely, recall that :math:`P_0` is the energy operator. Then :eq:`eq_p_and_j_conjugated_by_space_and_time_inversions` shows

.. math:: U(\Pcal) \ifrak P_0 U^{-1}(\Pcal) = \ifrak P_0

If :math:`U(\Pcal)` were anti-linear, then :math:`U(\Pcal) P_0 U^{-1}(\Pcal) = -P_0`. Then for any state :math:`\Psi` with positive energy, i.e., :math:`P_0 \Psi = p_0 \Psi`, we would have a state :math:`U^{-1}(\Pcal) \Psi` with negative energy :math:`-p_0`. Hence we conclude that :math:`U(\Pcal)` must be linear. The same argument shows also that :math:`U(\Tcal)` must be anti-linear (since :math:`{\Tcal_0}^0 = -1`).

As before, it'll be useful to rewrite :eq:`eq_p_and_j_conjugated_by_space_and_time_inversions` in terms of :math:`H, \Pbf, \Jbf, \Kbf` as follows

.. math::
	:label: eq_hpjk_conjugated_by_space_and_time_inversions

	\begin{alignat*}{3}
		U(\Pcal) &H U^{-1}(\Pcal) &&= &&H \\
		U(\Pcal) &\Pbf U^{-1}(\Pcal) &&= -&&\Pbf \\
		U(\Pcal) &\Jbf U^{-1}(\Pcal) &&= &&\Jbf \\
		U(\Pcal) &\Kbf U^{-1}(\Pcal) &&= -&&\Kbf \\
		U(\Tcal) &H U^{-1}(\Tcal) &&= &&H \\
		U(\Tcal) &\Pbf U^{-1}(\Tcal) &&= -&&\Pbf \\
		U(\Tcal) &\Jbf U^{-1}(\Tcal) &&= -&&\Jbf \\
		U(\Tcal) &\Kbf U^{-1}(\Tcal) &&= &&\Kbf \\
	\end{alignat*}

One can (and should) try to reconcile these implications with commonsense. For example, the :math:`3`-momentum :math:`\Pbf` changes direction under either space or time inversion as expected. Moreover, the spin (of for example a basketball) remains the same under space inversion because both the direction of the axis and the handedness of the rotation get reversed simultaneously, but it gets reversed under time inversion because the direction of rotation is reversed if time flows backwards.

In what follows we will work out the effects of space and time inversions on massive and massless particles, respectively.


.. _sec_space_inversion_for_massive_particles:

Space inversion for massive particles
+++++++++++++++++++++++++++++++++++++

We start by considering a state at rest :math:`\Psi_{k, \sigma}`, where :math:`k = (M, 0, 0, 0)` and :math:`\sigma` is an eigenvalue of :math:`J_3` under one of the spin representations discussed in :ref:`Representations of angular momenta <dropdown_repr_of_angular_momenta>`. Since the state is at rest and :math:`U(\Pcal)` commutes with :math:`J_3` according to :eq:`eq_p_and_j_conjugated_by_space_and_time_inversions`, we can write

.. math::
	:label: eq_space_inversion_on_massive_standard

	U(\Pcal) \Psi_{k, \sigma} = \eta \Psi_{k, \sigma}

where :math:`\eta` is a phase that depends a priori on :math:`\sigma`. It turns out, however, that :math:`\eta` is actually independent of :math:`\sigma`, and hence justifies the notation, since :math:`U(\Pcal)` commutes with the raising/lowering operators :math:`J_1 \pm \ifrak J_2` by :eq:`eq_hpjk_conjugated_by_space_and_time_inversions`.

To move on to the general case, we recall that the general formula :eq:`eq_def_of_one_particle_psi` takes the following form

.. math:: \Psi_{p, \sigma} = \sqrt{\frac{M}{p_0}} U(L(p)) \Psi_{k, \sigma}

We can calculate as follows

.. math::
	:label: eq_space_inversion_on_massive_general

	U(\Pcal) \Psi_{p, \sigma} = \sqrt{\frac{M}{p_0}} U(\Pcal L(p) \Pcal^{-1}) U(\Pcal) \Psi_{k, \sigma} = \eta~\sqrt{\frac{M}{p_0}} U(L(\Pcal p)) \Psi_{k, \sigma} = \eta \Psi_{\Pcal p, \sigma}

which generalizes :eq:`eq_space_inversion_on_massive_standard`. Such :math:`\eta` is known as the *intrinsic parity*, which is intrinsic to a particle species.


.. _sec_time_inversion_for_massive_particles:

Time inversion for massive particles
++++++++++++++++++++++++++++++++++++

Consider the same :math:`\Psi_{k, \sigma}` as in the space inversion case. Now since :math:`U(\Tcal)` anti-commutes with :math:`J_3` according to :eq:`eq_p_and_j_conjugated_by_space_and_time_inversions`, we can write

.. math:: U(\Tcal) \Psi_{k, \sigma} = \zeta_{\sigma} \Psi_{k, -\sigma}

where :math:`\zeta_{\sigma}` is a phase. Applying the raising/lowering operators and using :eq:`eq_rotation_j_matrix`, we can calculate the left-hand-side, recalling that :math:`U(\Tcal)` is anti-linear, as follows

.. math::

	(J_1 \pm \ifrak J_2) U(\Tcal) \Psi_{k, \sigma}
		&= -U(\Tcal) (J_1 \mp \ifrak J_2) \Psi_{k, \sigma} \\
		&= -U(\Tcal) \sqrt{(\jfrak \pm \sigma)(\jfrak \mp \sigma + 1)} \Psi_{k, \sigma \mp 1} \\
		&= -\zeta_{\sigma \mp 1} \sqrt{(\jfrak \pm 1)(\jfrak \mp \sigma + 1)} \Psi_{k, -\sigma \pm 1}

where :math:`\jfrak` is the particle spin, and the right-hand-side as follows

.. math::

	(J_1 \pm \ifrak J_2) \zeta_{\sigma} \Psi_{k, -\sigma} = \zeta_{\sigma} \sqrt{(\jfrak \pm 1)(\jfrak \mp \sigma + 1)} \Psi_{k, -\sigma \pm 1}

Equating the two sides, we see that :math:`\zeta_{\sigma} = -\zeta_{\sigma \pm 1}`. Up to an overall phase, we can set :math:`\zeta_{\sigma} = \zeta (-1)^{\jfrak - \sigma}` so that

.. math:: U(\Tcal) \Psi_{k, \sigma} = \zeta (-1)^{\jfrak - \sigma} \Psi_{k, -\sigma}

Here we have chosen to keep the option of a physically inconsequential phase :math:`\zeta` open. As in the case of space inversion, the formula generalizes to any :math:`4`-momentum :math:`p`

.. math::
	:label: eq_time_inversion_on_massive_general

	U(\Tcal) \Psi_{p, \sigma} = \zeta (-1)^{\jfrak - \sigma} \Psi_{\Pcal p, -\sigma}

since :math:`\Tcal L(p) \Tcal^{-1} = L(\Pcal p)`.

Space inversion for massless particles
++++++++++++++++++++++++++++++++++++++

Let's consider a state :math:`\Psi_{k, \sigma}` with :math:`k = (1, 0, 0, 1)` and :math:`\sigma` being the helicity, i.e., :math:`J_3 \Psi_{k, \sigma} = \sigma \Psi_{k, \sigma}`. Since :math:`U(\Pcal)` commutes with :math:`J_3`, the space inversion preserves :math:`\sigma`, just as in the massive case. However, since :math:`\Pcal` reverses the direction of motion, the helicity in the direction of motion actually reverses sign. It follows, in particular, that (massless) particles that respect the space inversion symmetry must come in companion with another particle of opposite helicity.

To spell out more details, note that since :math:`\Pcal` doesn't fix :math:`k`, it'll be convenient to introduce an additional rotation :math:`R_2`, which is defined to be a :math:`\pi`-rotation about the :math:`2`-axis, so that :math:`U(R_2) = \exp(\ifrak \pi J_2)` and :math:`R_2 \Pcal k = k`. Since :math:`U(R_2)` flips the sign of :math:`J_3`, as can be seen from the very definition of :math:`J_3` in :eq:`eq_u_lorentz_expansion`, we have

.. math:: U(R_2 \Pcal) \Psi_{k, \sigma} = \eta_{\sigma} \Psi_{k, -\sigma}

where we see indeed that the helicity reverses sign (when :math:`k` is fixed).

To move on to the general case, recall that the :math:`L(p)` that takes :math:`k` to :math:`p` consists of a boost :math:`B` defined by :eq:`eq_massless_boost` followed by a (chosen) pure rotation :math:`R(\pbf)` that takes :math:`(0, 0, |\pbf|)` to :math:`\pbf`. We calculate as follows

.. math::
	:label: eq_space_inversion_on_massless_undetermined_phase

	U(\Pcal) \Psi_{p, \sigma} &= p_0^{-1/2} U(\Pcal R(\pbf)B) \Psi_{k, \sigma} \\
		&= p_0^{-1/2} U(R(\pbf) B R_2^{-1}) U(R_2 \Pcal) \Psi_{k, \sigma} \\
		&= p_0^{-1/2} \eta_{\sigma} U(R(\pbf) R_2^{-1} B) \Psi_{k, -\sigma} \\
		&= \eta_{\sigma} \rho \Psi_{\Pcal p, -\sigma}

where :math:`\rho` is an extra phase due to the fact that although :math:`R(\pbf) R_2^{-1}` takes :math:`(0, 0, |\pbf|)` to :math:`-\pbf`, it may not be the chosen one.

To spell out :math:`\rho`, we need to be a bit more specific about :math:`R(\pbf)`. Following the usual convention of `spherical coordinates <https://en.wikipedia.org/wiki/Spherical_coordinate_system>`_, we can get from :math:`(0, 0, |\pbf|)` to :math:`\pbf` by first rotate (according to the right-handed rule) about the :math:`1`-axis at an angle :math:`0 \leq \phi \leq \pi`, known as the polar angle, and then rotate about the :math:`3`-axis at an angle :math:`0 \leq \theta < 2\pi`, known as the azimuthal angle. Now since we know that :math:`R(\pbf)R_2^{-1}` differs from :math:`R(-\pbf)` by a rotation about the :math:`3`-axis, we can figure the rotation out by examining their actions on some suitably generic :math:`\pbf`, for example, the unit vector along the :math:`2`-axis, which is fixed by :math:`R_2`. In this case :math:`R(-\pbf)` is a :math:`\pi/2`-rotation about the :math:`1`-axis, while :math:`R(\pbf)R_2^{-1}` is the same :math:`\pi/2`-rotation by the :math:`1`-axis, followed by a a :math:`\pi`-rotation about the :math:`3`-axis. Therefore we conclude that the difference is a :math:`\pi`-rotation about the :math:`3`-axis. In other words, we should have :math:`\rho = \exp(-\ifrak \pi \sigma)`. However, recalling that the helicity :math:`\sigma` may be a half-integer (thought not yet being found in Nature), there is a sign difference between :math:`\pm \pi`-rotations about the :math:`3`-axis. Without going into further details, we write down the general formula the space inversion as follows

.. math::
	:label: eq_space_inversion_on_massless_general

	U(\Pcal) \Psi_{p, \sigma} = \eta_{\sigma} \exp(\mp \ifrak \pi \sigma) \Psi_{\Pcal p, -\sigma}

where the sign depends on the sign of :math:`p_2` (which can be seen by playing the same game as above with the negative unit vector along the :math:`2`-axis).

Time inversion for massless particles
+++++++++++++++++++++++++++++++++++++

Let :math:`k = (1, 0, 0, 1)` as usual and consider the state :math:`\Psi_{k, \sigma}`. Since :math:`U(\Tcal)` anti-commutes with both :math:`\Pbf` and :math:`J_3` by :eq:`eq_hpjk_conjugated_by_space_and_time_inversions`, we have

.. math:: U(\Tcal) \Psi_{k, \sigma} = \Psi_{\Pcal k, -\sigma}.

Composing with the rotation :math:`R_2` as in the previous section to fix :math:`k`, we have

.. math:: U(R_2 \Tcal) \Psi_{k, \sigma} = \zeta_{\sigma} \Psi_{k, \sigma}

where :math:`\zeta_{\sigma}` is (yet another) phase. We see that, unlike the space inversion, the time inversion doesn't produce a doublet of opposite helicity. Processing as in the space inversion case, one can derive the following general formula similar to :eq:`eq_space_inversion_on_massless_general`

.. math::
	:label: eq_time_inversion_on_massless_general

	U(\Tcal) \Psi_{p, \sigma} = \zeta_{\sigma} \exp(\pm \ifrak \pi \sigma) \Psi_{\Pcal p, \sigma}

where the sign depends on the sign of :math:`p_2` as before.

Kramers' degeneracy
+++++++++++++++++++

We end our discussion about space and time inversions of one-particle states with an interesting observation on the squared time inversion :math:`U(\Tcal)^2`. It follows from both :eq:`eq_time_inversion_on_massive_general` and :eq:`eq_time_inversion_on_massless_general` that

.. math:: U(\Tcal)^2 \Psi_{p, \sigma} = (-1)^{2 s} \Psi_{p, \sigma}

where :math:`s \in \tfrac{1}{2} \Zbb` equals the spin :math:`\jfrak` in the massive case and the absolute helicity :math:`|\sigma|` in the massless case.

Hence in a non-interacting system consisting of an odd number of half-integer spin/helicity particles and any number of integer spin/helicity particles, we have

.. math::
	:label: eq_time_inversion_squared_reverse_sign

	U(\Tcal)^2 \Psi = -\Psi

Now for any eigenstate :math:`\Psi` of the Hamiltonian, there is an accompanying eigenstate :math:`U(\Tcal) \Psi` since :math:`U(\Tcal)` commutes with the Hamiltonian. The key observation then is that they are necessarily different states! To see this, let's suppose otherwise that :math:`U(\Tcal) \Psi = \zeta \Psi` represent the same state, where :math:`\zeta` is a phase. Then


.. math::

	U(\Tcal)^2 \Psi = U(\Tcal) \zeta \Psi = \zeta^{\ast} U(\Tcal) \Psi = |\zeta|^2 \Psi = \Psi

contradicts :eq:`eq_time_inversion_squared_reverse_sign`.

As a conclusion, we see that for such systems, any energy eigenvalue has at least a two-fold degeneracy. This is known as `Kramers' degeneracy <https://en.wikipedia.org/wiki/Kramers%27_theorem>`_.


.. rubric:: Footnotes

.. [#tedious_calc_of_commutations] These are some rather tedious and error-prone calculations, but in the end, we manage to arrive at the same results as stated in [Wei95]_ page 61.