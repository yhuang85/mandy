How to do Morse theory in contact topology?
===========================================

In this article, we explain how to cut any contact manifold into simple pieces as well as how to manipulate them from a Morse-theoretic viewpoint. The math of this article is based on many years of collaboration with K. Honda, but the interpretations are mostly mine. Moreover, the focus will be on ideas rather than technical details: one can hardly find any proofs here.


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

A Symplectic detour
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


Back to contact
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

The Game Plan
-------------

Separation of contact structure and topology
********************************************

The first key principle in the development of contact Morse theory is to separate the contact topological problem from the purely topological problem. Specifically, given a contact manifold :math:`(M, \xi)`, we always start from just any Morse function :math:`f: M \to \RR`. Suppose :math:`\dim M = 2n+1`, then :math:`f` induces a decomposition

.. math::
    :label: heegaardDecomposition

    M = H_1 \cup (\Sigma \times I) \cup H_2, \quad I \coloneqq [0,1],

where :math:`H_1` is a neighborhood of the union of the stable manifolds of all critical points of :math:`f` of index at most :math:`n`, and similarly :math:`H_2` is a neighborhood of the union of the unstable manifolds of all critical points of :math:`f` of index at least :math:`n+1`, and finally :math:`\Sigma` may be identified with either :math:`\p H_1` or :math:`\p H_2`. Such a decomposition is nothing but a higher-dimensional analogue of the `Heegaard decomposition <https://en.wikipedia.org/wiki/Heegaard_splitting>`_ for 3-manifolds.

The reason for decomposing :math:`M` as in :eq:`heegaardDecomposition` is because the contact structures :math:`\xi|_{H_1}, \xi|_{H_2}` can be made standard by Gromov's *h*-principle on isotropic/Legendrian approximations. Namely, the stable manifold of all critical points of index at most :math:`n` can be :math:`C^0`-approximated by isotropic submanifolds. It follows that :math:`H_1` becomes a neighborhood of a CW-complex built out of isotropic cells, which in turns carries a standard contact structure. The same holds for :math:`H_2` by reversing the direction of :math:`\nabla f`. Finally, the complement of :math:`H_1 \cup H_2` in :math:`M` is a product :math:`\Sigma \times I`, which, in some sense, carries all the information about :math:`\xi`.

To summarize, the decomposition :eq:`heegaardDecomposition` serves the purpose of separating contact topology from pure topology as follows:

* The handlebodies :math:`H_1, H_2` knows all about the topology of :math:`M` but nothing about :math:`\xi`.
* The middle layer :math:`\Sigma \times I` knows all about :math:`\xi` but nothing about :math:`M`.

.. note::
    The handlebodies :math:`H_1, H_2` in :eq:`heegaardDecomposition` are by no means unique, although the (restricted) contact structures are uniquely determined by their topological type.

Morsify the characteristic foliation
************************************

Given any hypersurface :math:`\Sigma \in (M, \xi)`, the *characteristic foliation* :math:`\Sigma_{\xi}` is a line field defined by

.. math::
    \Sigma_{\xi} \coloneqq \ker (d\alpha|_{T\Sigma \cap \xi}),

where :math:`\xi = \ker\alpha` is a contact form. Moreover, when both :math:`\Sigma` and :math:`\xi` are oriented, which will always be the case here, so is :math:`\Sigma_{\xi}` and it becomes a vector field (without any significance on the magnitude). Characteristic foliations play a crucial role in this story because they, to a great extent which will become clear later, uniquely determines the contact germ on :math:`\Sigma`.

Now the job of characterizing a contact structure on :math:`\Sigma \times I` boils down to characterizing the evolution of the contact germs on :math:`\Sigma \times t` for :math:`t \in I`, which, in turn, boils down to characterizing the evolution of a :math:`1`-parameter family of vector fields :math:`(\Sigma \times t)_{\xi}, t \in I`.

In the case of :math:`\dim M = 3`, the above job was successfully done by Giroux in [Gi91]_ and [Gi99]_, where he applied the results of such analysis to classify contact structures on a number of :math:`3`-manifolds. However, Giroux's argument from these two papers are rather specific to dimension :math:`3` and are of little use in higher dimensions [#falseBelieve]_. Specifically, the study of characteristic foliations on a :math:`2`-dimensional surface falls into a much bigger subject of studying dynamics of generic vector fields on surfaces. Rather mature and comprehensive theories on the later subject, such as the `Poincaré-Bendixson theorem <https://en.wikipedia.org/wiki/Poincar%C3%A9%E2%80%93Bendixson_theorem>`_, was developed long before contact topology was even recognized as an independent subject. However, it's indeed hopeless to track down every single trajectory of a generic vector field in dimensions greater than two due to the ubiquity of chaotic behavior.

The challenge is, therefore, to ensure the controllability of :math:`\Sigma_{\xi}` (e.g., as the gradient vector field of a Morse function) on sufficiently generic hypersurfaces :math:`\Sigma`. The basic idea is to wiggle :math:`\Sigma` almost everywhere to create attractors, built out of Morse critical points, which destroy any potential global dynamics of :math:`\Sigma_{\xi}`. The actual implementation of this idea is nearly perfect in dimension :math:`3` but much less so in higher dimensions. The details can be found in [HH18]_ and [HH19]_.

Summary
*******

Every (closed) contact manifold can be decomposed into three pieces: two standard contact handlebodies and a product :math:`\Sigma \times I`. The contact structure :math:`\xi|_{\Sigma \times I}` can be understood via the :math:`1`-parameter family of characteristic foliations :math:`\Sigma_t|_{\xi} \coloneqq (\Sigma \times t)_{\xi}` for :math:`t \in I`. The characteristic foliations :math:`\Sigma_t|_{\xi}, t \in I`, can be made Morse by a :math:`C^0`-small perturbation. Thus the problem is finally reduced to understanding a :math:`1`-parameter family of Morse functions on :math:`\Sigma`. More details about carrying out this game plan will be explained in the next section.


Second Blend of Morse Theory with Contact Structures
----------------------------------------------------

Recall in the first blend of Morse theory with contact structures, the result is a decomposition of :math:`(M, \xi)` into a bunch of contact handles. This approach appears to be somewhat useless since there is no way (that I know of) to connect two contact Morse functions through a family of contact Morse functions.

Instead, we'll use the ideas outlined in the :ref:`game plan<The Game Plan>` to build a contact Morse theory which works in families. To facilitate the exposition, let's use the following convention to indicate the dimension of the family of Morse functions under consideration. We say a Morse theory (of whatever flavor) is established at

* :math:`\pi_0`-level if Morse functions exist generally,
* :math:`\pi_1`-level if any two Morse functions are homotopic through Morse functions,
* :math:`\pi_2`-level if a circle-family of Morse functions can be realized as the boundary of a disk-family of Morse functions,
* and so on for :math:`\pi_k`-levels for :math:`k > 2`.

.. note::
    Critial points, among others, in families of Morse functions degenerate according to the standard `tranversality theory <https://en.wikipedia.org/wiki/Transversality_theorem>`_ on jet bundles. For example, critical points are nondegenerate at :math:`\pi_0`-level but may degenerate to birth-death type singularities at :math:`\pi_1`-level and swallowtails at :math:`\pi_2`-level and so on.

tbc...

.. rubric:: Footnotes

.. [#gradient] Morse theory is topological in nature and doesn't care about metric very much. In particular, it's more correct and convenient but unfortunately also more cumbersome to use `gradient-like vector fields <https://en.wikipedia.org/wiki/Gradient-like_vector_field>`_ instead.

.. [#falseBelieve] Ironically, the failure of Giroux's argument in dimensions greater than three went so far to even form a consensus that hypersurfaces in higher-dimensional contact manifolds are intractable and hopeless. It was at least the case when I entered the subject as a graduate student. From my own experience, there is nothing better than breaking false believes.

.. rubric:: References

.. [Don96] S\. Donaldson\. `Symplectic submanifolds and almost-complex geometry <https://projecteuclid.org/journals/journal-of-differential-geometry/volume-44/issue-4/Symplectic-submanifolds-and-almost-complex-geometry/10.4310/jdg/1214459407.full>`_

.. [EG91] Y\. Eliashberg and M\. Gromov\. `Convex symplectic manifolds <https://www.ihes.fr/~gromov/symplecticmanifolds/163/>`_

.. [Gi91] E\. Giroux\. `Convexité en topologie de contact <https://link.springer.com/article/10.1007%2FBF02566670>`_

.. [Gi99] E\. Giroux\. `Structures de contact en dimension trois et bifurcations des feuilletages de surfaces <https://arxiv.org/abs/math/9908178>`_

.. [Gi03] E\. Giroux\. `Géométrie de contact: de la dimension trois vers les dimensions supérieures <https://arxiv.org/abs/math/0305129>`_

.. [HH18] K\. Honda and Y\. Huang\. `Bypass attachments in higher-dimensional contact topology <https://arxiv.org/abs/1803.09142>`_

.. [HH19] K\. Honda and Y\. Huang\. `Convex hypersurface theory in contact topology <https://arxiv.org/abs/1907.06025>`_
