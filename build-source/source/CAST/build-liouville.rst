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

Description of the toy Liouville domain
***************************************

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

Now in this example :math:`W` has a particularly nice, but by no means general, property that :math:`X` is completely tangent to a :math:`3`-submanifold :math:`Y \coloneqq W|_{t = 0}` and is pointing away from :math:`Y` in the normal directions. This :math:`Y` is itself a mapping torus defined by :math:`A \in \op{Aut}(T^2)`, and is schematically drawn as the thick loop in the picture. This allows us to first decompose :math:`Y` before moving on to the (visually) more challenging case of :math:`W`.

The anatomy of :math:`Y`
************************
Strictly speaking, one can jump directly to the anatomy of the original :math:`W`. But the fact that we can only visualize things up to three dimensions makes this section a good warm-up for all the techniques we shall need in the :math:`4`-dimensional case.

Recall that :math:`Y = T^2 \times [0,C] / (x, y, 0) \sim (A(x, y), C)` which comes with the restricted Liouville vector field :math:`\p_s` where :math:`s \in [0, C]`. The plan to dissect :math:`Y` involves the following steps:

#. Find a cross section in :math:`Y` such that every trajectories of :math:`\p_s` pass through it.
#. Decompose the cross section into (topologically speaking) disks.
#. Flow each such disk in the direction of :math:`\p_s` until it hits itself again to form what we'll call a *flow tube*.
#. At this point, we have decomposed :math:`Y` into a bunch of flow tubes. But just as explained in the previous section, the boundary of a flow tube has "sides" which are completely tangent to :math:`\p_s`. We need to "tilt" the sides to make the flow tubes less degenerate from the dynamics point of view.

Let's carry out the plan. For the first step, there is a fairly obvious candidate :math:`T^2_0 \coloneqq  T^2 \times \{0\}` which we will take as the preferred cross section. For the second step we take also what we think is the simplest to work with.

.. _fig_a_decomposition_of_torus:

.. sidebar:: A decomposition of :math:`T^2_0`

    .. image:: static/torus-decomposition.svg
        :width: 400px

On the right we decompose :math:`T^2_0` into four rectangles :math:`S_y, S_b, S_p` and :math:`S_g`, distinguished by four colors yellow, brown, purple and green, respectively. As usual :math:`T^2_0` is drawn as a square such that the top and the bottom edges are identified, so are the left and the right edges. Underneath, we also draw the images under :math:`A` of all four rectangles, as well as how they overlap with the original ones, respectively. Moreover, the edges of the :math:`S`'s are colored to indicate the adjacent rectangles.

In this particular scenario, one can see that there are essentially two types of flow tubes, with different topologies. Indeed, the brown flow tube is topologically a solid torus since :math:`S_b \cap A(S_b)` is connected. On the other hand, the other three flow tubes are genus :math:`2` handlebodies since :math:`S \cap A(S)` has two components, where :math:`S \in \{ S_y, S_p, S_g \}`. In the picture, besides the obvious overlaps, pairs of regions with the same color, which are translations of each other, are also identified in the flow tube. This completes step three.

Description of a single flow tube
+++++++++++++++++++++++++++++++++

In the following, we will denote a flow tube by :math:`\tau(S)`. Then, obviously, the Liouville vector field :math:`X` is pointing out along :math:`S \setminus A(S)`, pointing in along :math:`A(S) \setminus S`, and tangential along the side :math:`\p S \times [0, C]`, which we need to tilt. As a piece of terminology, we will refer to the regions where :math:`X` is pointing out as **positive** and the regions where :math:`X` is pointing in as **negative**.

In general, we can divide the tilting procedure into two types: the uniform and the mixed. Roughly speaking, after uniform tilting, :math:`\p S \times [0, C]` becomes entirely transverse to :math:`X`, but after mixed tilting, :math:`\p S \times [0, C]` further decomposes into two regions so that :math:`X` is outward-pointing along one region, inward-pointing along the other, and remains tangential along the borderline. In this section, we will only use the uniform tilting for simplicity, but as we will see in the actual decomposition of :math:`W` in the next section, mixed tilting is inevitable.

.. sidebar:: The flow tube :math:`\tau(S_b)`

    .. image:: static/3d-mapping-torus-1.svg
        :width: 400px

To work out the details, let's start with the slightly simpler :math:`\tau(S_b)`. The picture on the right depicts how step four is carried out for :math:`\tau(S_b)`. Namely, it consists of two sub-steps -- first tilt and then round corners. The leftmost figure shows the original :math:`\tau(S_b)`, where the red region is negative, the blue region is positive, and the gray regions on the top and bottom are identified by :math:`A`. The passage to the middle figure is the process of a uniform tilting. Here we have two choices: either shrink or expand :math:`S_b` as :math:`s` runs from :math:`0` to :math:`C`. The later is chosen in the picture, which turns the entire :math:`\p S_b \times [0, C]` negative. Finally we round the corners: Most of the borderlines between the red and blue regions can be rounded to *folds*, except for four points -- highlighted as thick dots in the middle figure -- which are rounded to *cusps* as shown in the rightmost figure.

.. note::
    The names of "folds" and "cusps" are borrowed from singularity theory of smooth maps or transversality theory of R. Thom. In fact, there is a rabbit hole of singularities characterized by certain "stratification" scheme which leads to nothing but a mess. It's fair to say, at this point, that either we can manage so that all flow tubes possess only folds and cusps or the sought-after decompositions are simply useless.

To summarize, we have transformed :math:`\tau(S_b)` into a (smooth) solid torus whose boundary admits a decomposition

.. math::

    \p \tau(S_b) = R_+ \cup R_-

into positive and negative regions such that  :math:`R_+` is the disjoint union of two disks, each of which has a boundary :math:`S^1` which can be further decomposed into two semicircles along a :math:`0`-sphere as follows

.. math::
    :label: eq_two_cusps

    S^1 = U_0 \cup_{S^0} U_1

Here the two semicircles :math:`U_0, U_1` are folds and the two points :math:`S^0` are cusps. If, on the contrary, we had decided to shrink :math:`S_b` as :math:`s` runs from :math:`0` to :math:`C`, the side :math:`\p S_b \times [0, C]` would have become entirely positive. We would then end up with a different :math:`\tau(S_b)` where the descriptions of :math:`R_{\pm}` switch places.

The descriptions of the other three :math:`\tau(S_y), \tau(S_p)` and :math:`\tau(S_g)` are not so different even though they have different topologies than :math:`\tau(S_b)`. Let's go through :math:`\tau(S_y)` quickly to further familiarize ourselves with this procedure.

.. sidebar:: The flow tube  :math:`\tau(S_y)` (before tilting)

    .. image:: static/3d-mapping-torus-2.svg
        :width: 400px

On the right we have the very similar picture of :math:`\tau(S_y)` where the red region is negative, the blue region is positive, and the side is tangent to :math:`X`. Note that the blue region also consists of two pieces, although the tiny triangular piece at the lower-right corner of :math:`S_y` is not easy to see. Unfortunately :math:`A(S_y)` is drawn in two pieces, where the ":math:`\cdots`" symbols are supposed to indicate the appropriate identifications. One should compare with the :ref:`picture <fig_a_decomposition_of_torus>` of :math:`S_y` and its image under :math:`A`.

Depending on how the side is (uniformly) tilted, the resulting :math:`\tau(S_y)` will turn out to be different. For example, if we slightly shrink :math:`S_y` as :math:`s` runs from :math:`0` to :math:`C`, then in the corresponding :math:`R_- \subset \p \tau(S_y)` consists of two (topologically speaking) disk components. One of them is just like the one described by :eq:`eq_two_cusps`, and the other will have four arcs of folds and four cusps, which are marked by thick dots in the picture as before.

Fit flow tubes together
+++++++++++++++++++++++



.. rubric:: References

.. [Hua19] Y\. Huang\. `A dynamical construction of Liouville domains <https://arxiv.org/abs/1910.14132v2>`_

.. [Mit95] Y\. Mitsumatsu\. `Anosov flows and non-Stein symplectic manifolds <http://www.numdam.org/item/AIF_1995__45_5_1407_0>`_
