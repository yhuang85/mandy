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
    :label: toy_liouville_domain

    W \coloneqq (T^2 \times [-1, 1]) \times [0, C] / (x, y, t, 0) \sim (\phi(x, y, t), C), \quad\forall (x, y, t) \in T^2 \times [-1, 1].

where :math:`T^2 = \RR^2 / \ZZ^2`, :math:`C > 0` is a constant to be determined later, and :math:`\phi: T^2 \times [-1, 1] \to T^2 \times [-1, 1]` is a linear map defined by

.. math::
    :label: toy_monodromy

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

The anatomy of :math:`Y^3`
**************************

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

In the following, we will denote a flow tube by :math:`\tau(S)`. Then, obviously, the Liouville vector field :math:`X` is pointing out along :math:`S \setminus A(S)`, pointing in along :math:`A(S) \setminus S`, and tangential along the side :math:`\p S \times [0, C]`, which we need to tilt. As a piece of terminology, we will refer to the regions where :math:`X` is pointing out as *positive* and the regions where :math:`X` is pointing in as *negative*.

In general, we can divide the tilting procedure into two types: the uniform and the mixed. Roughly speaking, after uniform tilting, :math:`\p S \times [0, C]` becomes entirely transverse to :math:`X`, but after mixed tilting, :math:`\p S \times [0, C]` further decomposes into two regions so that :math:`X` is outward-pointing along one region, inward-pointing along the other, and remains tangential along the borderline. Although uniform tilting may appear to be easier to visualize, as we will see, when several flow tubes stuck together, it's not always possible to avoid mixed tilting. So to keep things as simple as possible, we shall start with uniform tilting and deal with mixed tilting when it gets in the way.

.. sidebar:: The flow tube :math:`\tau(S_b)`

    .. image:: static/3d-mapping-torus-1.svg
        :width: 400px

To work out the details, let's start with the slightly simpler :math:`\tau(S_b)`. The picture on the right depicts how step four is carried out for :math:`\tau(S_b)`. Namely, it consists of two sub-steps -- first tilt and then round corners. The leftmost figure shows the original :math:`\tau(S_b)`, where the red region is negative, the blue region is positive, and the gray regions on the top and bottom are identified by :math:`A`. The passage to the middle figure is the process of a uniform tilting. Here we have two choices: either shrink or expand :math:`S_b` as :math:`s` runs from :math:`0` to :math:`C`. The later is chosen in the picture, which turns the entire :math:`\p S_b \times [0, C]` negative. Finally we round the corners: Most of the borderlines between the red and blue regions can be rounded to *folds*, except for four points -- highlighted as thick dots in the middle figure -- which are rounded to *cusps* as shown in the rightmost figure.

.. note::
    The names of "folds" and "cusps" are borrowed from singularity theory of smooth maps or transversality theory of `R. Thom <https://en.wikipedia.org/wiki/Ren%C3%A9_Thom>`_. In fact, there is a rabbit hole of singularities/tangencies characterized by certain "stratification" scheme. Although it's not our interest to explore these structures in any systematic way, they will find us in one way or another and we shall have no choice but to follow the path in front of us. The first time we will be forced to face tangencies "deeper" than folds and cusps is when we try to decompose :math:`W^4`.

To summarize, we have transformed :math:`\tau(S_b)` into a (smooth) solid torus whose boundary admits a decomposition

.. math::

    \p \tau(S_b) = R_+ \cup R_-

into positive and negative regions such that  :math:`R_+` is the disjoint union of two disks, each of which has a boundary :math:`S^1` which can be further decomposed into two semicircles along a :math:`0`-sphere as follows

.. math::
    :label: eq_two_cusps

    S^1 = U_+ \cup_{S^0} U_-

Here the two semicircles :math:`U_+, U_-` are folds so that :math:`X` points from :math:`R_{\pm}` to :math:`R_{\mp}` along :math:`U_{\pm}`, respectively, and the two points :math:`S^0` are cusps.

If, on the contrary, we had decided to shrink :math:`S_b` as :math:`s` runs from :math:`0` to :math:`C`, the side :math:`\p S_b \times [0, C]` would have become entirely positive. We would then end up with a different :math:`\tau(S_b)` where the descriptions of :math:`R_{\pm}` switch places.

The descriptions of the other three :math:`\tau(S_y), \tau(S_p)` and :math:`\tau(S_g)` are not so different even though they have different topologies than :math:`\tau(S_b)`. Let's go through :math:`\tau(S_y)` quickly to further familiarize ourselves with this procedure.

.. sidebar:: The flow tube  :math:`\tau(S_y)` (before tilting)

    .. image:: static/3d-mapping-torus-2.svg
        :width: 400px

On the right we have the very similar picture of :math:`\tau(S_y)` where the red region is negative, the blue region is positive, and the side is tangent to :math:`X`. Note that the blue region also consists of two pieces, although the tiny triangular piece at the lower-right corner of :math:`S_y` is not easy to see. Unfortunately :math:`A(S_y)` is drawn in two pieces, where the ":math:`\cdots`" symbols are supposed to indicate the appropriate identifications. One should compare with the :ref:`picture <fig_a_decomposition_of_torus>` of :math:`S_y` and its image under :math:`A`.

Depending on how the side is (uniformly) tilted, the resulting :math:`\tau(S_y)` will turn out to be different. For example, if we slightly shrink :math:`S_y` as :math:`s` runs from :math:`0` to :math:`C`, then in the corresponding :math:`R_- \subset \p \tau(S_y)` consists of two (topologically speaking) disk components. One of them is just like the one described by :eq:`eq_two_cusps`, and the other will have four arcs of folds and four cusps, which are marked by thick dots in the picture as before.

Fit flow tubes together
+++++++++++++++++++++++

Now that we understand the structures - mostly importantly the singularities - of individual flow tubes, it's time to fit them together. In this procedure, as we will see, old singularities may disappear and new ones may emerge, but the basic principle remains the same: we will keep the species of possible singularities to only folds and cusps.

When two flow tubes stuck together, a part of their boundaries overlap and it's the most crucial to understand the behavior near the boundary of this overlapping region. To illustrate the possible scenarios, we shall take a close look at how :math:`\tau(S_b)` and :math:`\tau(S_y)` are fit together. We assume that the following uniform tilting strategy:

    :math:`S_b` is expanding and :math:`S_y` is shrinking as :math:`s` runs from :math:`0` to :math:`C`,

so that they fit each other along the overlapping boundaries.

.. _fig_overlap_between_Sb_and_Sy:

.. sidebar:: The overlap between :math:`\tau(S_b)` and :math:`\tau(S_y)`

    .. image:: static/3d-fit-brown-yellow.svg
        :width: 400px

It turns out that :math:`\p \tau(S_b) \cap \p \tau(S_y)` is a union of two annuli, whose boundary circles are depicted as the thickened line segments in the right-hand-side picture. Note that the segments parallel (at least before tilting) to the :math:`s`-direction, which are needed to join the endpoints to complete the boundary circles, are not drawn in the picture. Moreover, the two magnified regions on the left-hand-side are duplicate of each other by the identification.

Since our models are essentially made out of linear objects, e.g., polygons, it's actually easier not to round the corners just for the sake of ending up with smooth objects. There should however be no difficulty in going back-and-forth between previously encountered objects such as folds and cusps and their (piecewise) linear models.

In light of the above remark, we can see, in fact, that the components of :math:`\p (\p \tau(S_b) \cap \p \tau(S_y))`, i.e., the thickened lines in the above picture, are all piecewise linear. This allows us to investigate the gluing patterns near the line segments (*edges*) and the corners (*vertices*) separately. In the above picture, the vertices are labeled by :math:`a, b, \cdots, p` and the edges will be denoted by, for example :math:`\overline{ab}`.

.. _fig_fit_along_edges:

.. sidebar:: Fit along edges

    .. image:: static/3d-fit-edges.svg
        :width: 400px

Let's start with the simpler case of fit-along-edges. There are essentially two scenarios as depicted in right-hand-side picture, where the reference vector field :math:`X` is assumed to be vertical. Namely, either the two meeting surfaces have no tangencies (with respect to :math:`X`) along the edge but the gluing results in a fold tangency as shown on the left-hand-side or one of the meeting surface have a fold tangency and the result of gluing has no singularity as shown on the right-hand-side. We call the edge in the first scenario *obtuse*, and in the second scenario *acute*. Here the meeting surfaces are the boundaries of :math:`3`-dimensional objects which are not explicitly specified in the picture -- they are unambiguously determined as the two sides of the overlapping part.

In the case of fitting :math:`\tau(S_b)` and :math:`\tau(S_y)` together, the edges :math:`\overline{ab}` and :math:`\overline{bc}` are acute, while the edge :math:`\overline{ca}` (which is almost parallel to the :math:`s`-direction) is obtuse, for example.

.. _fig_fit_vertices:

.. sidebar:: Fit around vertices

    .. image:: static/3d-fit-vertices.svg
        :width: 400px

Now we move on to the scenarios of fit-around-vertices as illustrated in the picture to the right. The vertices are where the edges meet, so it's the easiest to organize by the types of edges involved, whether they are transverse (to :math:`X`) or folds or cusps. For example, on the upper-left corner we have a transverse sheet meeting with a folded sheet, and end up with a cusp -- hence the name TF-C, read as *Transverse-Fold-to-Cusp*. Similarly, on the upper-right corner we have two folds fit together and end up with a new fold -- hence FF-F. The bottom two look similar to each other except that the assembling pieces are different, and so are the results. Indeed, one can even make a third assemble which should be TF-C, which shall look a bit different, but turns out to be equivalent to the one on the upper-left corner. Finally, we note that other combinations of T, F, and C are either trivial (e.g., TF-T is equivalent to an acute edge) or impossible (e.g., TT-T), and therefore we do not draw pictures for them.

Let's apply our classification of the vertices to the :ref:`overlap <fig_overlap_between_Sb_and_Sy>` between :math:`\tau(S_b)` and :math:`\tau(S_y)`. For example, the vertex :math:`a` lies on a fold when viewed in :math:`\tau(S_b)`, and on a transverse face when in :math:`\tau(S_y)`. Therefore it's a TF-C type vertex, and one expects a cusp at :math:`a` after fitting :math:`\tau(S_b)` and :math:`\tau(S_y)` together, which is indeed the case. Similarly, the vertex :math:`b` lies at the interface of a fold and a cusp, and is henceforth a FC-T type vertex. It's hopefully clear at this point that the process of fitting together flow tubes, at least in the generic cases, is rather mechanical and straightforward. However, we're still missing one scenario where the so-called uniform tilting becomes inadequate in fitting more flow tubes together. We shall then wrap up this section with a discussion of the mixed tilting scenario.

.. sidebar:: Mixed tilting becomes necessary

    .. image:: static/torus-decomposition-tilted.svg
        :width: 400px

Suppose :math:`\tau(S_b)` and :math:`\tau(S_y)` have been glued together, and we will try to fit in :math:`\tau(S_g)` next. Recall that our chose to uniformly tilt the sides of :math:`\tau(S_b)` and :math:`\tau(S_y)` such that :math:`S_b` is shrinking and :math:`S_y` is expanding as :math:`s` runs from :math:`C` down to :math:`0`. In the picture to the right, we draw the decompositions of the cross section :math:`T^2` at the gluing level :math:`s=C` and at a level :math:`s=C-\epsilon` slightly below it. It's clear from the picture that :math:`S_g` (i.e., the green region) has a mixed (tilting) behavior along the boundary -- the vertical sides which are adjacent to :math:`S_b` is expanding, and the horizontal sides which are (partially) adjacent to :math:`S_y` are shrinking, as :math:`s` runs from :math:`C` to :math:`0`. Similar mixed tilting occurs when we fit the last piece :math:`\tau(S_p)` in as well.

After a moment of thoughts, it should become clear that such mixed tilting introduces nothing new. Namely, one can think of the "vertex" where the expanding side meets the shrinking side as a vertex just as in the case of fitting two flow tubes together considered above -- especially the TF-C type vertex.

Summary
+++++++

Let's summarize our knowledge so far about the decomposition of the mapping torus :math:`Y` into a collection of flow tubes as follows.

* Each flow tube :math:`\tau(S)` is a solid handlebody whose boundary can be decomposed as :math:`\p \tau(S) = R_+ \cup R_-` such that :math:`X` is outward-pointing along :math:`R_+` and inward-pointing along :math:`R_-`. Moreover, the borderline :math:`\p R_+ = \p R_-` can be thought of as a polygon whose edges are folds and vertices are cusps.

* When two flow tubes :math:`\tau(S_1)` and :math:`\tau(S_2)` fit together, they are glued along a region :math:`K \coloneqq \p \tau(S_1) \cap \p \tau(S_2)`. The boundary :math:`\p K` can, again, be thought of as a polygon such that the edges are :ref:`either acute or obtuse <fig_fit_along_edges>`, and the vertices are one of the four types: :ref:`TF-C, TC-F, FC-T and FF-F <fig_fit_vertices>`.


The anatomy of :math:`W^4`
**************************

We're now ready to take on the real objects of interests -- the building blocks of the :math:`4`-dimensional Liouville domain :math:`W`, which are nothing but the :math:`4`-dimensional flow tubes. Let's also note that all the discussions so far have been purely topological, i.e., we haven't touched on any symplectic and/or contact structures at all. Indeed, we shall continue to discuss the relative positions between :math:`X` and the flow tubes before turning into symplectic and/or contact aspects of the story.

Description of a single flow tube
+++++++++++++++++++++++++++++++++

Recall that in the :math:`3`-dimensional case, a flow tube :math:`\tau(S)` is a partial mapping tori based on a square :math:`S \subset T^2`. Now in the :math:`4`-dimensional case, we need to upgrade :math:`S` to a cube :math:`S \times [-1, 1]` equipped with the standard contact structure, where the :math:`[-1, 1]` factor comes from :eq:`toy_liouville_domain`. To keep notations simple and consistent, let's adopt the following renaming convention:

    We will from now on write :math:`S` in place of :math:`S \times [-1, 1]` for the :math:`3`-dimensional cube. Inheriting from the :math:`3`-dimensional case, we cover the transverse slice :math:`T^2 \times [-1, 1] = S_b \cup S_y \cup S_g \cup S_p` by four cubes.

.. sidebar:: The flow tube :math:`\tau(S_b)` at :math:`s=C` with uniform tilting

    .. image:: static/4d-mapping-torus-section-uniform-titling.svg
        :width: 400px

In the picture on the right, we draw the (most important) section of :math:`\tau(S_b)` at :math:`s=C`. Here we lose, unfortunately, the luxury to draw the :math:`s`-direction because it's the :math:`4`-th dimension. However, the familiarity with the :math:`3`-dimensional case should help with the imagination of the :math:`s`-direction and its tilting as well!

As before, the Liouville vector field :math:`X` is outward-pointing along a, topological speaking, solid torus :math:`S_b \setminus \phi(S_b)`, and inward-pointing along the two balls :math:`\phi(S_b) \setminus S_b`. Here :math:`\phi` is defined in :eq:`toy_monodromy`. Moreover :math:`X` is tangent to :math:`\p S_b \times [0, C]`, which needs to be tilted in some way which we now elaborate.

Let's start with the naive but simple approach where :math:`S_b` is assumed to be expanding as :math:`s` runs from :math:`0` to :math:`C`. This is the direct analog of the :math:`3`-dimensional :math:`\tau(S_b)` discussed in the previous section. As a consequence, the side :math:`\p S_b \times [0, C]` joins the negative part. The borderline between :math:`R_+` and :math:`R_-` is therefore the boundary of (the closure of) :math:`S_b \setminus \phi(S_b)`, which is a torus. Moreover, the torus itself is split into two annuli by :math:`\p S_b \cap \p \phi(S_b)`, depicted as the blue circles in the above picture, such that :math:`X` defines folds along the annuli and cusps along the circles. It might worth noting that we have so far managed to keep only folds and cusps, even though they may come in families. This turns out to be a luxury that we are about to leave behind.

The above uniform tilting, albeit simple, isn't quite a building block we need for :math:`W`. Indeed, the fact that a part of :math:`\p \tau(S_b)` lies on :math:`\p W` demands that the :math:`t`-direction is necessarily shrinking as :math:`s: 0 \to C`. So we will adopt the following mixed tilting strategy:

.. _block_mixed_tilting_instruction_4d:

    As :math:`s: 0 \to C`, :math:`S_b` is shrinking along the top and bottom sides :math:`t = \pm 1` and expanding along the rest of :math:`\p S_b`.

Of course the decomposition :math:`\p \tau(S_b) = R_+ \cup R_-` has changed by the mixed tilting. Namely, :math:`R_+` becomes the union of the following two pieces:

* :math:`S_b \setminus \phi(S_b)` at :math:`s=C`;
* :math:`\p S_b|_{t = \pm 1} \times [0, C]` where :math:`\p S_b|_{t = \pm 1}` denotes the top and bottom sides of :math:`S_b`.

Topologically speaking :math:`R_+` is now a genus :math:`3` handlebody, and therefore the borderline between :math:`R_+` and :math:`R_-` is a genus :math:`3` surface. However, the topological types of these objects are not our main interests -- we are interested in the relative position between :math:`X` and :math:`\p \tau(S_b)`. We know that :math:`X` is tangent to :math:`\p \tau(S_b)` along :math:`\p R_+ = \p R_-`, but the exact form of tangency turns out to be more complicated than just folds and cusps.

.. _fig_mixed_tilting_4d:

.. sidebar:: The flow tube :math:`\tau(S_b)` at :math:`s=C` with mixed tilting

    .. image:: static/4d-mapping-torus-section-mixed-tilting.svg
        :align: center
        :width: 400px

It follows from the :ref:`instruction <block_mixed_tilting_instruction_4d>` of the mixed tilting that there are two separating loops (actually, polygons) :math:`\gamma_1, \gamma_2 \subset \p R_+`, which are depicted blue in the picture to the right, such that :math:`\p R_+` are folds away from them. As we will see shortly, there are cusps and even one-level deeper tangencies along :math:`\gamma_1` and :math:`\gamma_2`. In what follows, we shall consider only :math:`\gamma_1` as the typical case.

To avoid getting lost in the many strata, let's paint the region bounded by :math:`\gamma_1` in green. Since we cannot visualize :math:`4`-dimensional objects, the best bet to look at :math:`3`-dimensional slices. To this end, let's pick two points :math:`a, b \in \gamma_1` around the upper-left corner. Moreover, take :math:`2`-dimensional slices at :math:`a, b` transverse to the corresponding edges of :math:`\gamma_1`, respectively. These are represented by small shaded squares in the picture.

Now we add the :math:`s`-direction to the slices as shown in the lower part of the picture, where the green arcs are the intersections between the slices and the green regions enclosed by :math:`\gamma_1`, respectively. The rest of the picture should be straightforward as it simply reduces to the :math:`3`-dimensional case. In particular, it's obvious that both :math:`a` and :math:`b` are cusps. However, they are cusps in "different directions" -- a fact that forces the corner of :math:`\gamma_1` between :math:`a` and :math:`b` (marked by a fat dot in the picture) to be a type of tangency which is neither a fold nor a cusp, namely, it's yet one-level deeper in the hierarchy.

For the moment, we shall live with the unexplained "deeper tangency" and summarize what we know about the relative position between :math:`X` and :math:`\tau(S_b)` as follows:

* :math:`X` is pointing into :math:`\tau(S_b)` along the interior of :math:`R_+`, which happens to be a genus :math:`3` handlebody, and pointing out along the interior of :math:`R_-` (whose topology we didn't bother to find out), and is tangential along :math:`\p R_+ = \p R_-`.

* There are two loops :math:`\gamma_1, \gamma_2 \subset \p R_+` away from which there are folds.

* There is a finite number (which happens to be :math:`8` in this case) of points on :math:`\gamma_1 \cup \gamma_2` away from which there are cusps.

* These points on :math:`\gamma_1 \cup \gamma_2`, which are marked by the fat dots in the above :ref:`picture <fig_mixed_tilting_4d>` are tangencies which we haven't explored yet.

It should be clear at this point that although we managed to work out most of the tangencies/transversalities between :math:`X` and :math:`\p \tau(S_b)`, it's getting quite messy already at the level of describing one single flow tube, not to say the headaches we must face when we have to fit many of them together. Therefore, in order to go further, we must pause our planned-next-step of fitting flow tubes together, and re-examine the perspective from which we look at the relative position between :math:`X` and :math:`\p \tau(S_b)`

From straight vector field to flat surface
++++++++++++++++++++++++++++++++++++++++++

So far we have been keeping the Liouville vector field straight -- it goes straight up. This appears to be a natural perspective since the flow tubes are mapping tori by definition. However, when it comes to understanding the relative position between :math:`X` and the boundary of a flow tube :math:`\p \tau(S)`, it's not really necessary to remember the fact that :math:`\tau(S)` is a mapping torus. Indeed, we shall replace :math:`\p \tau(S)` by a general (hyper)surface :math:`\Sigma` of dimension :math:`2` and :math:`3`, and recover all the local models of tangencies considered before. The main feature of this approach is that we shall keep :math:`\Sigma` flat while letting :math:`X` "draw" the tangencies.

.. sidebar:: Two viewpoints of local tangencies

    .. image:: static/3d-vector-field-local-model.svg
        :width: 400px

Let's start with the case :math:`\dim \Sigma = 2`. Here we have two types of tangencies: folds and cusps, which are depicted in the picture to the right. Specifically, on the left, we arrange so that :math:`X` is vertical, and on the right, we arrange so that :math:`\Sigma` is flat. In either perspective, we can see a stratification

.. math::

    \Sigma = \Sigma^{(0)} \cup \Sigma^{(1)} \cup \Sigma^{(2)}

such that the following hold:

* :math:`\Sigma^{(0)}` has codimension :math:`0` (i.e., surfaces) and is transverse to :math:`X`.

* :math:`\Sigma^{(1)}` has codimension :math:`1` (i.e., lines) along which :math:`X` is tangent to :math:`\Sigma^{(0)}` and transverse to :math:`\Sigma^{(1)}`. These are folds.

* :math:`\Sigma^{(2)}` has codimension :math:`2` (i.e., points) at which :math:`X` is tangent to :math:`\Sigma^{(1)}`. These are cusps.

In the flat-surface model, we can see that near a fold, :math:`\Sigma^{(0)}` are the left and right half spaces, which meet along the middle line :math:`\Sigma^{(1)}`, and there is no :math:`\Sigma^{(2)}`. Near a cusp, :math:`\Sigma^{(0)}` are again two half spaces, but the borderline consists of edges :math:`\Sigma^{(1)}` and vertices :math:`\Sigma^{(2)}`.

Next, let's look at how two surfaces fit together. We've seen both the :ref:`edge-fitting <fig_fit_along_edges>` and the :ref:`vertex-fitting <fig_fit_vertices>` before in the straight-vector-field model, and now we will redo everything in the flat-surface model.


.. rubric:: References

.. [Hua19] Y\. Huang\. `A dynamical construction of Liouville domains <https://arxiv.org/abs/1910.14132v2>`_

.. [Mit95] Y\. Mitsumatsu\. `Anosov flows and non-Stein symplectic manifolds <http://www.numdam.org/item/AIF_1995__45_5_1407_0>`_
