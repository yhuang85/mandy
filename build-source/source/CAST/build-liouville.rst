How to build Liouville domains?
===============================

A *Liouville domain* :math:`(W, \omega)` is by definition a compact exact symplectic manifold with convex boundary. Namely, we demand the existence of a *Liouville form* :math:`\lambda` such that :math:`\omega = d\lambda` and a *Liouville vector field* :math:`X` -- defined by :math:`\iota_X \omega = \lambda` -- such that :math:`X` is everywhere outward pointing along :math:`\p W`. The question of which (compact with boundary) manifolds can carry a Liouville structure seems difficult, but there are obvious and not-so-obvious conditions that :math:`W` must satisfy to become a Liouville domain. For example, one obvious condition is that :math:`W` must carry an almost complex structure, and one not-so-obvious condition is that :math:`W` must not be :math:`S^3 \times [0,1]`, even though it does carry an almost complex structure. The later constraint comes from some rather complex geometric analysis (e.g. pseudo-holomorphic curves or mathematical gauge theory) which works only when :math:`\dim W = 4`. In any case, there is no chance of answering the question if one cannot decompose a Liouville domain into "simple" pieces.

There is one noteworthy special case though when :math:`(W, \omega, \lambda, X)` carries one additional assumption that :math:`X` is gradient-like for a Morse function :math:`f: W \to \RR`. In this case :math:`f` induces a handle decomposition of :math:`W` where the handles come with standard symplectic structures. These manifolds are known as *Weinstein domains* -- named after A. Weinstein. On the one hand, the construction of Weinstein provides us an interesting class of Liouville domains, but on the other hand, the Morse condition severely constrains the topology of :math:`W`, e.g., the Morse index of any critical point of :math:`f` cannot exceed :math:`\tfrac{1}{2} \dim W`.

One crucial difference between Weinstein domains and Liouville domains in general is the dynamics of :math:`X`. For Weinstein domains, :math:`X` has no (broken) periodic orbits and :math:`W` can be gradually built from sub-level sets of the associated Morse function :math:`f`, which doesn't exist for a general Liouville domain. It should nevertheless be noted that it is not the presence of periodic orbits of :math:`X` that distinguishes the Liouville domains that are not Weinstein -- it is the orbits that neither close up nor leave :math:`W` in any finite time. These "trapped" orbits must be cut into pieces, and our strategy to achieve this is to, roughly speaking, "stop at the first time return". This ends up cutting :math:`W` into a collections of (partial) mapping tori. The plan is therefore to understand/classify these mapping tori and hope that they will be suffice to build just about any Liouville domain.

.. note::

    The idea of relating Liouville domains to mapping tori was first explored in [Hua19]_, but the treatment there covers only the useless case: the mapping tori that have only positive ends and no negative ends -- they are themselves Liouville domains but cannot be assembled together.

A Toy Model
-----------

It's hard to imagine, out of nowhere, how the (Liouville) mapping tori should look like, not to mention how they fit together. So we shall start by describing a concrete example, which was first studied by Y. Mitsumatsu [Mit95]_, and then examine its anatomy.

Description of the Liouville domain
***********************************

Consider the following :math:`4`-manifold

.. math::

    W \coloneqq (T^2 \times [-1, 1]) \times [0, C] / (x, y, t, 0) \sim (\phi(x, y, t), C), \quad\forall (x, y, t) \in T^2 \times [-1, 1].

where :math:`T^2 = \RR^2 / \ZZ^2`, :math:`C > 0` is a constant to be determined later, and :math:`\phi: T^2 \times [-1, 1] \to T^2 \times [-1, 1]` is a linear map defined by

.. math::

    \phi = \begin{pmatrix}
        A & 0 \\
        0 & c
    \end{pmatrix},
    \textrm{ where }
    A = \begin{pmatrix}
        2 & 1 \\
        1 & 1
    \end{pmatrix}
    \textrm{ and }
    0 < c < 1.

The key property of :math:`A` is that it has two distinct real eigenvalues :math:`a_+ > 1` and :math:`a_- < 1`. The eigenspaces of :math:`A` then define two linear :math:`1`-forms :math:`\alpha_{\pm} \in \Omega^1(T^2)` such that :math:`A^{\ast} \alpha_{\pm} = a_{\pm} \alpha_{\pm}`, respectively. We'd like to regard :math:`T^2 \times [-1, 1]` as a contact manifold with contact form :math:`\alpha_+ + t\alpha_-`, and define the Liouville form to be

.. math::

    \lambda \coloneqq e^s(t \alpha_+ + \alpha_-), \textrm{ where } s \in [0, C].

Now for :math:`\lambda` to be well-defined on :math:`W`, we need

.. math::

    \phi^{\ast} (\lambda|_{s=C}) = e^C (ct \cdot a_+ \alpha_+ + a_- \alpha_-) = \lambda|_{s=0}.

This can be achieved by setting :math:`c = a_- / a_+` and :math:`C = \log a_+` (recall :math:`a_+a_- = 1`).

.. sidebar:: A schematic picture of :math:`W`

    .. image:: static/toy_mapping_torus.svg
        :width: 400px

One can easily compute the Liouville vector field :math:`X = \p_s`. Now the quadruple :math:`(W, \omega, \lambda, X)` isn't quite a Liouville domain for two reasons: first, :math:`W` has corners, and second, :math:`X` is tangent to :math:`\p (W|_{s=0}) \times [0, C]` -- as highlighted by the shaded region in the right-handle-side picture. Both issues can be resolved at once by a corner-rounding procedure highlighted by the blue arc in the "magnified" (and straightened) picture.

Now in this example :math:`W` has a particular nice, but by no means general, property that :math:`X` is completely tangent to a :math:`3`-submanifold :math:`Y \coloneqq W|_{t = 0}` and is pointing away from :math:`Y` in the normal directions. This :math:`Y` is itself a mapping torus defined by :math:`A \in \op{Aut}(T^2)`, and is schematically drawn as the thick loop in the picture. This allows us to first decompose :math:`Y` before moving on to the (visually) more challenging case of :math:`W`.

The anatomy of :math:`Y`
************************
Strictly speaking, one can jump directly to the anatomy of the original :math:`W`. But the fact that we can only visualize things up to three dimensions makes this section a good warm-up for all the techniques we shall need in the :math:`4`-dimensional case.

Recall that :math:`Y = T^2 \times [0,C] / (x, y, 0) \sim (A(x, y), C)` which comes with the restricted Liouville vector field :math:`\p_s` where :math:`s \in [0, C]`. The plan to dissect :math:`Y` involves the following steps:

#. Find a cross section in :math:`Y` such that every trajectories of :math:`\p_s` pass through it.
#. Decompose the cross section into (topologically speaking) disks.
#. Flow each such disk in the direction of :math:`\p_s` until it hits itself again to form what we'll call a *flow tube*.
#. At this point, we have decomposed :math:`Y` into a bunch of flow tubes. But just as explained in the previous section, the boundary of a flow tube has "sides" which are completely tangent to :math:`\p_s`. We need to "tilt" the sides to make the flow tubes less degenerate from the dynamics point of view.

Let's carry out the plan. For the first step, we have a fairly obvious candidate :math:`T^2_0 \coloneqq  T^2 \times \{0\}` and we will take it as the preferred cross section.


.. rubric:: References

.. [Hua19] Y\. Huang\. `A dynamical construction of Liouville domains <https://arxiv.org/abs/1910.14132v2>`_

.. [Mit95] Y\. Mitsumatsu\. `Anosov flows and non-Stein symplectic manifolds <http://www.numdam.org/item/AIF_1995__45_5_1407_0>`_
