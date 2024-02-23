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
    [\Psi] \coloneqq \{ c\Psi ~|~ c \in \Cbb \setminus 0 \}

where :math:`\Psi` is a nonzero vector in the Hilbert space. It is, however, rather inconvenient to have to deal with :math:`[\Psi]` all the time. So instead, we will almost always pick a representative :math:`\Psi`, often out of a natural choice, and call it a *state vector*, and keep in mind that anything physically meaningful must not be sensitive to a scalar multiplication.

    Throughout this post we always assume that state vectors are normalized so that :math:`||\Psi|| = 1`.

In fact, we don't really care about the states themselves either, because they are more of an abstraction rather than something one can physically measure. What we do care about are the (Hermitian) inner products between state vectors, denoted by :math:`(\Psi, \Phi)`. According to the so-called `Copenhagen interpretation <https://en.wikipedia.org/wiki/Copenhagen_interpretation>`_ of quantum mechanics, such inner product represents an *amplitude*, i.e., its squared norm gives the probability of finding a state :math:`[\Psi]` in :math:`[\Phi]` if we ever perform a measurement. We can write this statement as an equation as follows:

.. math::
    P([\Psi] \to [\Phi]) = |(\Psi, \Phi)|^2

In particular, the probability of finding any state in itself is one, due to the normalization above.

What is a Symmetry?
-------------------

We start with a *symmetry transformation*, by which we mean a transformation that preserves all quantities that one can ever measure about a system. Since it is the probabilities, rather than the states themselves, that are measurable, one is led to define a quantum symmetry transformation as a transformation of states :math:`T` such that

.. math::
    :label: t_preserves_probability

    P(T[\Psi] \to T[\Phi]) = P([\Psi] \to [\Phi])

for any states :math:`[\Psi]` and :math:`[\Phi]`. Now a theorem of E. Wigner asserts that such :math:`T` can be realized either as a linear unitary or as an anti-linear anti-unitary transformation :math:`U = U(T)` of state vectors in the sense that :math:`[U\Psi] = T[\Psi]` for any :math:`\Psi`. In other words, :math:`U` satisfies either

.. math::
    U(c\Psi) = cU\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)

or

.. math::
    U(c\Psi) = c^{\ast} U\Psi, \quad (U\Psi, U\Phi) = (\Psi, \Phi)^{\ast}

where :math:`c` is any complex number.

.. dropdown:: Proof of Wigner's theorem
    :animate: fade-in-slide-down

    The construction of a realization :math:`U` of :math:`T` takes the following steps.

    :Step 1: Fix an orthonormal basis :math:`\Psi_i, i \in \Nbb`, of the Hilbert space.

    :Step 2:  For each :math:`\Psi_i`, choose a unit vector :math:`\Psi'_i` such that :math:`P[\Psi_i] = [\Psi'_i]`. Then :math:`\Psi'_i, i \in \Nbb`, also form an orthonormal basis by :eq:`t_preserves_probability`. We'd like to define :math:`U` by asking

        .. math::
            :label: u_definition

            U\Psi_i = \Psi'_i

        for all :math:`i`, and extend by (anti-)linearity. But this is not going to realize :math:`T` in general because we haven't fixed the extra degrees of freedom -- the phases of :math:`\Psi'_i`.

    :Step 3: We fix the phases of :math:`\Psi'_k, k \geq 1`, relative to :math:`\Psi'_0`, by asking

        .. math::
            T[\Psi_0 + \Psi_k] = [\Psi'_0 + \Psi'_k].

        To see why this is possible, note first that :math:`T[\Psi_0 + \Psi_k] = [\alpha \Psi'_0 + \beta \Psi'_k]`, where :math:`\alpha, \beta` are phase factors, due to :eq:`t_preserves_probability` and the basis being orthonormal. Now :math:`[\alpha \Psi'_0 + \beta \Psi'_k] = [\Psi'_0 + (\beta/\alpha) \Psi'_k]` and we can absorb the phase :math:`\beta/\alpha` into the definition of :math:`\Psi'_k`. This is indeed the best one can do, because the last one degree of freedom, which is to multiply all :math:`\Psi'_i` by a phase, cannot be fixed.

    :Step 4: We have so far specified the value of :math:`U` on all of :math:`\Psi_i, i \geq 0`, and :math:`\Psi_0 + \Psi_k, k \geq 1`. Notice that all the coefficients of :math:`\Psi` are real. It is therefore instructive to ask what :math:`\Psi_0 + \i \Psi_1` should be. By the same argument as in the previous step, we can write

        .. math::
            T[\Psi_0 + \i \Psi_1] = [\Psi'_0 + c \Psi'_1]

        where :math:`c` is a phase. Let's apply :eq:`t_preserves_probability` once again as follows

        .. math::
            \sqrt{2} &= \left| \left( [\Psi_0 + \i \Psi_1], [\Psi_0 + \Psi_1] \right) \right| \\
                &= \left| \left( T[\Psi_0 + \i \Psi_1], T[\Psi_0 + \Psi_1] \right) \right| \\
                &= \left| \left( [\Psi'_0 + c \Psi'_1], [\Psi'_0 + \Psi'_1] \right) \right| \\
                &= |1 + c|

        It follows that :math:`c = \pm\i`, which correspond to :math:`U` being (complex) linear or anti-linear, respectively.

    At this point, we can extend :math:`U` to either a linear or anti-linear map of the Hilbert space. But we'll not be bothered about any further formal argument, including showing that (anti-)linearity must be coupled with (anti-)unitarity, respectively.

.. note::
    The *adjoint* of a linear operator :math:`A` is another linear operator :math:`A^{\dagger}` such that

    .. math::
        (\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)

    for all any two state vectors :math:`\Psi` and :math:`\Phi`. On the other hand, the adjoint of an anti-linear :math:`A` is another anti-linear :math:`A^{\dagger}` such that

    .. math::
        (\Psi, A\Phi) = (A^{\dagger} \Psi, \Phi)^{\ast}

    A (anti-)unitary operator :math:`U` thus satisfies :math:`U^{\dagger} = U^{-1}`.

In general we're not interested in just one symmetry transformation, but rather a group -- whether continuous or discrete -- of symmetry transformations, or just symmetry for short. In particular, if :math:`T_1, T_2` are two symmetry transformations, then we'd like :math:`T_2 T_1` to also be a symmetry transformation. In light of the :math:`U`-realization of symmetry transformations discussed above, we can rephrase this condition as

.. math::
    :label: u_depends_on_psi

    U(T_2 T_1) \Psi = e^{\i \theta(T_1, T_2, \Psi)} U(T_2) U(T_1) \Psi

where :math:`\theta(T_1, T_2, \Psi)` is an angle, which depends a priori on :math:`T_1, T_2`, and :math:`\Psi`.

It turns out, however, the angle :math:`\theta(T_1, T_2, \Psi)` cannot depend on the state because if we apply :eq:`u_depends_on_psi` to the sum of two linearly independent state vectors :math:`\Psi_A + \Psi_B`, then we'll find

.. math::
    e^{\pm \i \theta(\Psi_A)} \Psi_A + e^{\pm \i \theta(\Psi_B)} \Psi_B = e^{\pm \i \theta(\Psi_A + \Psi_B)} (\Psi_A + \Psi_B)

where we have suppressed the dependency of :math:`\theta` on :math:`T`, and the signs correspond to the cases of :math:`U` being linear or anti-linear, respectively. In any case, it follows that

.. math::
    e^{\pm \i \theta(\Psi_A)} = e^{\pm \i \theta(\Psi_B)} = e^{\pm \i \theta(\Psi_A + \Psi_B)}

which says nothing but the independence of :math:`\theta` on :math:`\Psi`.

.. todo::
    While the argument here appears to be purely mathematical, Weinberg pointed out in the book the potential inabilities to create a state like :math:`\Psi_A + \Psi_B`. More precisely, he mentioned the general believe that it's impossible to prepare a superposition of two states, one with integer total angular momentum and the other with half-integer total angular momentum, in which case there will be a "super-selection rule" between different classes of states. After all, one Hilbert space may just not be enough to describe all states. It'd be nice to elaborate a bit more on the super-selection rules.

We can now simplify :eq:`u_depends_on_psi` to the following

.. math::
    :label: u_not_depend_on_psi

    U(T_2 T_1) = e^{\i \theta(T_1, T_2)} U(T_2) U(T_1)

which, in mathematical terms, says that :math:`U` furnishes a *projective representation* of :math:`T`, or a representation up to a phase. It becomes a genuine representation if the phase is constantly one.

.. important::
    We always assume that :math:`U` furnishes a genuine representation of :math:`T` since, as it turns out, the phase factor in :eq:`u_not_depend_on_psi` is more of a mathematical artifact than something that bears any physical significance.

Continuous symmetry
*******************

Besides a handful of important discrete symmetries such as the time, charge, and parity conjugations, most of the interesting symmetries come in a continuous family, mathematically known as *Lie groups*. Note that continuous symmetries are necessarily unitary (and linear) because they can be continuously deformed into the identity, which is obviously unitary.

In fact, it will be of great importance to just look at the symmetry up to the first order at the identity transformation, mathematically known as the *Lie algebra*. Let :math:`\theta` be an element in the Lie algebra such that :math:`T(\theta) = 1 + \theta` up to the first order. We can expand :math:`U(T(\theta))` in a power series as follows

.. math::
    :label: u_expansion

    U(T(\theta)) = 1 + \i \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots

where :math:`\theta^a` are the (real) components of :math:`\theta`, and :math:`u_a` are operators independent of :math:`\theta`, and as a convention, repeated indexes are summed up. Here we put a :math:`\i` in front of the linear term so that the unitarity of :math:`U` implies that :math:`u_a` are Hermitian.

Now let :math:`\eta` be another element of the Lie algebra, and expand both sides of :math:`U(T(\eta)) U(T(\theta)) = U(T(\eta) T(\theta))` as follows

.. math::
    :nowrap:

    \begin{eqnarray}
        U(T(\eta)) U(T(\theta)) &=& \left( 1 + \i \eta^a u_a + \tfrac{1}{2} \eta^a \eta^b u_{ab} + \cdots \right) \left( 1 + \i \theta^a u_a + \tfrac{1}{2} \theta^a \theta^b u_{ab} + \cdots \right) \\
            &=& 1 + \i (\eta^a + \theta^a) u_a \blue{- \eta^a \theta^b u_a u_b} + \cdots \\
        \\
        U(T(\eta) T(\theta)) &=& U \left( 1 + \eta + \theta + f_{ab} \eta^a \theta^b + \cdots \right) \\
            &=& 1 + \blue{\i} \left( \eta^c + \theta^c + \blue{f^c_{ab} \eta^a \theta^b} + \cdots \right) \blue{u_c} + \blue{\tfrac{1}{2}} \left( \blue{\eta^a + \theta^a} + \cdots \right) \left( \blue{\eta^b + \theta^b} + \cdots \right) \blue{u_{ab}} + \cdots
    \end{eqnarray}

where :math:`f^c_{ab}` are the coefficients of the expansion of :math:`T(f(\eta, \theta)) = T(\eta) T(\theta)`. Equating the coefficients of :math:`\eta^a \theta^b`, i.e., the terms colored in blue, we get

.. math::
    -u_a u_b = \i f^c_{ab} u_c + u_{ab} \quad \Longrightarrow \quad u_{ab} = -u_a u_b - \i f^c_{ab} u_c.

It implies that one can calculate the higher-order operator :math:`u_{ab}` from the lower-order ones, assuming of course that we know the structure of the symmetry (Lie) group/algebra. In fact, this bootstrapping procedure can be continued to all orders, but we'll not be bothered about the details.

Next, note that :math:`u_{ab} = u_{ba}` since they are just partial derivatives. It follows that

.. math::
    [u_a, u_b] \coloneqq u_a u_b - u_b u_a = \i (f^c_{ba} - f^c_{ab}) u_c \eqqcolon \i \Gamma^c_{ab} u_c

where the bracket is known as the *Lie bracket* and :math:`\Gamma^c_{ab}` are known as the *structure constants*.

We conclude the general discussion about continuous symmetry by considering a special, but important, case when :math:`T` is additive in the sense that :math:`T(\eta) T(\theta) = T(\eta + \theta)`. Notable examples of such symmetry include translations and rotations about a fixed axis. In this case :math:`f` vanishes, and it follows from :eq:`u_expansion` that

.. math::
    :label: u_abelian_as_exponential

    U(T(\theta)) = \lim_{N \to \infty} (U(T(\theta / N)))^N = \lim_{N \to \infty} (1 + \i \theta^a u_a / N)^N = \op{exp}(\i \theta^a u_a)

Lorentz symmetry
****************

A particularly prominent continuous symmetry in our physical world is the Lorentz symmetry postulated by Einstein's special relativity, which supersedes the Galilean symmetry, which is respected by the Newtonian mechanics. Lorentz symmetry is a symmetry that acts on the (flat) spacetime and preserves the so-called *proper time*

.. math::
    :label: proper_time

    d\tau^2 = dx_0^2 - dx_1^2 - dx_2^2 - dx_3^2 \eqqcolon \eta^{\mu \nu} dx_{\mu} dx_{\nu}

where

1. :math:`x_0` is also known as the time, and sometimes denoted by :math:`t`,
2. the speed of light is set to :math:`1`, and
3. :math:`\eta = \op{diag}(1, -1, -1, -1)` and the indexes :math:`\mu, \nu` run from :math:`0` to :math:`3`.

.. note::
    We will follow the common convention in physics that greek letters such as :math:`\mu, \nu, \dots` run from :math:`0` to :math:`3`, while roman letters such as :math:`i, j, \dots` run from :math:`1` to :math:`3`. Moreover, we often write :math:`x` for a spacetime point :math:`(x_0, x_1, x_2, x_3)`, and :math:`\xbf` for a spatial point :math:`(x_1, x_2, x_3)`.

.. dropdown:: Einstein's special theory of relativity
    :animate: fade-in-slide-down

    Using the notations introduced above, we can rewrite :eq:`proper_time` as :math:`d\tau^2 = dt^2 - d\xbf^2`, so that it's obvious that if a particle travels at the speed of light in one inertial frame, i.e., :math:`|d\xbf / dt| = 1`, and equivalently :math:`d\tau^2 = 0`, then it travels at the speed of light in any other inertial frame, in direct contradiction with Newtonian mechanics.

    Instead of working with the spacetime coordinates, it can sometimes be convenient to work with the "dual" energy-momentum coordinates, also known as the *four momentum*. The transition can be done by imagining a particle of mass :math:`m`, and defining :math:`p = (E, \pbf) \coloneqq m dx / d\tau`. It follows from :eq:`proper_time` that

    .. math::
        1 = (dt / d\tau)^2 - (d\xbf / d\tau)^2 ~\Longrightarrow~ m^2 = (m dt / d\tau)^2 - (m d\xbf / d\tau)^2 ~\Longrightarrow~ m^2 = E^2 - \pbf^2

    which looks just like :eq:`proper_time`, and indeed, the mass (in our convention) is invariant in all inertial frames.

    One can also recover Newtonian mechanics at the low-speed limit (i.e., :math:`|\vbf| \ll 1`) using :math:`d\tau / dt = \sqrt{1 - \vbf^2}` as follows

    .. math::
        \begin{eqnarray}
            \pbf &=& m d\xbf / d\tau = \frac{m \vbf}{\sqrt{1 - \vbf^2}} = m \vbf + O(\vbf^3) \\
            E &=& m dt / d\tau = m + \tfrac{1}{2} m \vbf^2 + O(\vbf^4)
        \end{eqnarray}
