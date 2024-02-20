The Quantum Theory of Fields (S. Weinberg)
==========================================

.. warning::
    This note is work in progress.

This note covers the three volumes [Wei95]_, [Wei96]_, and [Wei00]_, written by S. Weinberg on quantum field theory. What I like the most about these books is his attempt to make logical deductions from the most basic principles, in particular, the principle of symmetries, rather than to make analogies to experience, e.g., from classical physics (from a historical perspective). Such an endeavor may not always be possible because, after all, physics is about how we interpret Nature based on nothing but experience, and *not* about how Nature actually works. By the same token, the arguments that are considered logical here should really be interpreted as "reasonable-sounding", and have nothing to do with what mathematician would call "rigorous".

The order of the sections correspond roughly to the order of the chapters of the books.

Quantum One-Particle States
---------------------------

Quantum theory postulates that *any* physical state (of the world) can be represented by a *ray* in some complex Hilbert space. It's worth noting that it is the state, rather than the Hilbert space, that we actually care about. Let's write a state as 

.. math::
    [\Psi] \equiv \{ c\Psi ~|~ c \in \Cbb \setminus 0 \}

where :math:`\Psi` is a nonzero vector in the Hilbert space. It is, however, rather inconvenient to have to deal with :math:`[\Psi]` all the time. So instead, we will almost always pick a representative :math:`\Psi`, often out of a natural choice, and call it a *state vector*, and keep in mind that anything physically meaningful must not be sensitive to a scalar multiplication.

    Throughout this post we always assume that state vectors are normalized so that :math:`||\Psi|| = 1`.

In fact, we don't really care about the states themselves either, because they are more of an abstraction rather than something one can physically measure. What we do care about are the (Hermitian) inner products between state vectors, denoted by :math:`(\Psi, \Phi)`. According to the so-called `Copenhagen interpretation <https://en.wikipedia.org/wiki/Copenhagen_interpretation>`_ of quantum mechanics, such inner product represents an *amplitude*, i.e., a square root of the probability, of finding a state :math:`[\Psi]` in :math:`[\Phi]` if we ever perform a measurement. We can write this statement as an equation as follows:

.. math::
    P([\Psi] \to [\Phi]) = |(\Psi, \Phi)|^2

In particular, the probability of finding any state in itself is one, due to the normalization above.

Symmetry
********

By a *symmetry transformation* we mean a transformation that preserves all quantities that one can ever measure about a system. Since it is the probabilities, rather than the states themselves, that are measurable, one is led to define a quantum symmetry transformation as a transformation of states :math:`T` such that

.. math::
    :label: t_preserves_probability

    P(T[\Psi] \to T[\Phi]) = P([\Psi] \to [\Phi])

for any states :math:`[\Psi]` and :math:`[\Phi]`. Now a theorem of E. Wigner asserts that such :math:`T` can be realized either as a linear unitary or as an anti-linear anti-unitary transformation :math:`U \equiv U(T)` of state vectors in the sense that :math:`[U\Psi] = T[\Psi]` for any :math:`\Psi`. In other words, :math:`U` satisfies either

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

    :Step 3: We fix the phases of :math:`\Psi'_k, k \geq 1`, relative to :math:`\Psi'_0`, by asking :math:`T[\Psi_0 + \Psi_k] = [\Psi'_0 + \Psi'_k]`. To see why this is possible, note first that :math:`T[\Psi_0 + \Psi_k] = [\alpha \Psi'_0 + \beta \Psi'_k]`, where :math:`\alpha, \beta` are phase factors, due to :eq:`t_preserves_probability` and the basis being orthonormal. Now :math:`[\alpha \Psi'_0 + \beta \Psi'_k] = [\Psi'_0 + (\beta/\alpha) \Psi'_k]` and we can absorb the phase :math:`\beta/\alpha` into the definition of :math:`\Psi'_k`. This is indeed the best one can do, because the last one degree of freedom, which is to multiply all :math:`\Psi'_i` by a phase, cannot be fixed.
    
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
    
    At this point, we can extend :math:`U` to either a linear or anti-linear map of the Hilbert space. But we'll not be bothered with any further formal argument, including showing that (anti-)linearity must be coupled with (anti-)unitarity, respectively.

..
    Obviously, it makes no sense to describe a state of the whole world. On the contrary, it's a good start to try to describe the state(s) of possibly the simplest physical object: one particle. But what is a particle anyway? The answer to this seemingly innocent question is really not that easy. But at the very least, a particle should be a point :math:`(t, \xbb)` in the (flat) :math:`4`-dimensional spacetime :math:`\Rbb^{1, 3}`. Here :math:`t` represents time and :math:`\xbb \equiv (x_1, x_2, x_3)` represents a point in the :math:`3`-dimensional Euclidean space. Einstein's special relativity postulates that the so-called proper time

    .. math::
        d\tau^2 = -dt^2 + dx_1^2 + dx_2^2 + dx_3^2

    is invariant with respect to any inertial frame.
