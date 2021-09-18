How to do Morse theory in contact topology?
===========================================

In this article, we explain how to cut any contact manifold into simple pieces from a Morse-theoretic viewpoint.


Morse Theory in Topology
------------------------

The most basic and general idea of Morse theory is to understand global properties of a topological space by cutting it into simple pieces and keep track of the assembly process. There are many variations and complications that might sneak behind the simply looking idea. For example, the topological space itself may be `singular <https://en.wikipedia.org/wiki/Stratified_Morse_theory>`_ or `infinite dimensional <https://en.wikipedia.org/wiki/Floer_homology>`_, and the "simple pieces" and "assembly process" may take very `different forms <https://en.wikipedia.org/wiki/Triangulation_(topology)>`_.

We will be concerned with the most classical Morse theory here, where spaces are smooth finite-dimensional manifolds :math:`M` and decompositions of :math:`M` are given by the so-called Morse functions :math:`f: M \to \RR`, which are nothing but *generic* smooth functions. More precisely, the only data we care about :math:`f` is a sequence of *singular values*

.. math::

    \dots < c_{k-1} < c_k < c_{k+1} < \dots

i.e., the values at critical points where :math:`df=0`, which is bounded on both sides if :math:`M` is compact.

Any value in between two adjacent singular values is a *regular value*. Let's take, for example, two regular values :math:`a_{k-1} \in (c_{k-1}, c_k)` and :math:`a_k \in (c_k, c_{k+1})`, jamming :math:`c_k` in between. Then we cut out of :math:`M` a piece bounded between two smooth hypersurfaces :math:`\Sigma_{k-1} \coloneqq f^{-1} (a_{k-1})` and :math:`\Sigma_k \coloneqq f^{-1} (a_k)`. This piece is simple enough in the sense that it's always topologically equivalent to a *handle attachment* regardless of what manifold we're talking about: it only depends on a local invariant called the Morse index of :math:`c_k`.

Morse theory aims at recognizing manifolds by their handle decompositions. Most importantly, although one same manifold may have many different handle decompositions, as many as there are generic functions, they are all connected to each other via a sequence of handle manipulations. It is therefore fair to say the following:

    Morse theory is all about manipulations of handles: creations, cancellations, and isotopies.


First Blend of Morse theory with Contact Structures
---------------------------------------------------

A Symplectic Detour
*******************

Let's actually start by talking about the attempt of blending Morse theory into symplectic topology initiated Y. Eliashberg and M. Gromov [EG91]_. Let :math:`(W, \omega)` be a symplectic manifold and :math:`f: W \to \RR` be a Morse function. The most obvious guess of a compatibility condition between :math:`f` and :math:`\omega`, at least from a Morse theoretic viewpoint, would be that the flow of the gradient vector field :math:`\nabla f \eqqcolon X`, with respect to some metric [#gradient]_ , preserves :math:`\omega`. This is not a very good guess though since :math:`X` would then necessarily be volume preserving, and in turn :math:`f` can have neither local minima nor maxima.

It's somewhat unfortunate that no reasonable compatibility condition exists between :math:`f` and :math:`\omega` which would simply allow :math:`f` to simultaneously possess local minima and maxima as a usual smooth function does. The compromise, which is by no means obvious, is to ask :math:`X` to exponentially expand :math:`\omega`, i.e., :math:`\Lcal_X \omega = \omega`. Quite obviously such thing cannot exist on a closed manifold since :math:`X` exponentially expands the volume. In fact, with a bit more work, one realizes that the Morse index of any critical point of :math:`f` cannot exceed :math:`\tfrac{1}{2} \dim W`. In some sense, we just compromised half of the topology.

Recall the two major players in the usual Morse theory are regular level sets and critical points (or equivalently, handles). Not surprisingly, they receive additional structures from the compatibility condition with a symplectic structure.

Symplectic regular level sets
    For any regular value :math:`c`, the regular level set :math:`\Sigma \coloneqq f^{-1} (c) \subset (W, \omega)` is a contact manifold with contact form :math:`\alpha \coloneqq i_X \omega`.

Symplectic handles
    Assuming :math:`\dim W = 2n`, an index-:math:`k` handle corresponds to an index-:math:`k` critical point always takes the following standard form

    .. math::

        \begin{aligned}
            H_k &= [-1, 1]^k_{x_1, \cdots, x_k} \times [-1, 1]^{2n-k}_{x_{k+1}, \cdots, x_n, y_1, \cdots, y_n},

            \omega|_{H_k} &= \sum_{i=1}^{2n} dx_i \wedge dy_i,

            X|_{H_k} &= \sum_{i=1}^k ( -x_i \p_{x_i} + 2y_i \p_{y_i} ) + \frac{1}{2} \sum_{j=k+1}^{2n} ( x_j \p_{x_j} + y_j \p_{y_j} ).
        \end{aligned}

Such handles are known as *Weinstein handles* because A. Weinstein first wrote down these normal forms. A particularly important special case is when :math:`k = n` and we we call these handles *critical* because only these handles actually carry nontrivial symplectic information.


Back to Contact
***************

The compatibility between a contact manifold :math:`(M, \xi)` and a Morse function :math:`f: M \to \RR`, by analogy, asks the flow of :math:`X \coloneqq \nabla f` to preserve :math:`\xi`.  This turns out to be a much more flexible condition than its symplectic counterpart. Indeed, E. Giroux [Gi03]_ argued that every (closed) contact manifold admits a compatible Morse function.

The regular level sets and critical points of a contact Morse function look quite different from the symplectic case, mostly because the compatibility condition now reads :math:`\Lcal_X \alpha = g \alpha` where :math:`\xi = \ker \alpha` is a contact form and :math:`g` can be *any* function. In particular, the flow of :math:`X` doesn't have to expand the contact volume at all.

Contact regular level sets
    A regular level set :math:`\Sigma \subset (M, \xi)` can be decomposed into three pieces

    .. math::
        :label: hypersurfaceDecomposition

        \Sigma = R_+ \cup \Gamma \cup R_-,

where :math:`\Gamma = \{ \alpha (X) = 0 \}` and :math:`R_{\pm} = \{ \pm \alpha (X) > 0 \}`, respectively. It turns out that :math:`(\Gamma, \xi|_{\Gamma})` is itself a contact manifold of dimension :math:`\dim M - 2`, i.e., it's a codimension-:math:`2` contact submanifold of :math:`(M, \xi)`. Moreover :math:`( R_{\pm}, d\alpha|_{R_{\pm}} )` are symplectic manifolds. However, they don't necessarily carry a Morse structure as described above.

    .. note::

        We will rewrite the decomposition :eq:`hypersurfaceDecomposition` as :math:`\Sigma = R_+ \cup_{\Gamma} R_-` to highlight the viewpoint that :math:`\Sigma` can be obtained by gluing (closures of) :math:`R_{\pm}` along the common boundary :math:`\Gamma`.

Contact handles
    Assuming :math:`\dim M = 2n+1` and :math:`k \leq n`, an index-:math:`k` contact handle always takes the following standard form

    .. math::

        \begin{aligned}
            H_k &= [-1, 1]^k_{x_1, \cdots, x_k} \times [-1, 1]^{2n-k}_{x_{k+1}, \cdots, x_n, y_1, \cdots, y_n} \times [-1, 1]_z,

            \alpha|_{H_k} &= dz - \sum_{i=1}^n y_i dx_i,

            X|_{H_k} &= \sum_{i=1}^k ( -x_i \p_{x_i} + 2y_i \p_{y_i} ) + \frac{1}{2} \sum_{j=k+1}^{2n} ( x_j \p_{x_j} + y_j \p_{y_j} ) + z dz.
        \end{aligned}

    Note that :math:`X|_{H_k}` exponentially expands the contact volume in the above model. For :math:`k \geq n+1`, one can simply reverse the signs of :math:`X|_{H_k}` in the above equation, in which case :math:`X|_{H_k}` exponentially contracts the contact volume. These handles will just be called contact handles since nobody was interested in registering them as trademarks.

Summary
*******

The fact that every (closed) contact manifold admits a compatible Morse function means that one can build any contact manifold from the standard-looking contact handles. However, it doesn't really give us much more grip on the contact manifold itself because such handle decompositions are by no means unique. Indeed, most of the power of Morse theory lies in the ability to connect different choices of Morse functions by homotopies. Such homotopies or more generally the flexibility of contact Morse functions are unfortunately not available from [Gi03]_ due to the global nature of the argument, which is more-or-less a replica of an argument of S. Donaldson [Don96]_ for symplectic manifolds which are far more rigid.

We will follow a completely different path to build a hopefully more useful contact Morse theory. The main tools will be hypersurfaces and characteristic foliations on them. These tools are native to contact topology and were extensively used by D. Bennequin, Eliashberg, Giroux among many others for various purposes in the early days of the subject.


.. rubric:: Footnotes

.. [#gradient] Morse theory is topological in nature and doesn't care about metric very much. In particular, it's more correct and convenient but unfortunately also more cumbersome to use `gradient-like vector fields <https://en.wikipedia.org/wiki/Gradient-like_vector_field>`_ instead.


.. rubric:: References

.. [EG91] Y\. Eliashberg and M\. Gromov\. `Convex symplectic manifolds <https://www.ihes.fr/~gromov/symplecticmanifolds/163/>`_

.. [Gi03] E\. Giroux\. `Géométrie de contact: de la dimension trois vers les dimensions supérieures <https://arxiv.org/abs/math/0305129>`_

.. [Don96] S\. Donaldson\. `Symplectic submanifolds and almost-complex geometry <https://projecteuclid.org/journals/journal-of-differential-geometry/volume-44/issue-4/Symplectic-submanifolds-and-almost-complex-geometry/10.4310/jdg/1214459407.full>`_