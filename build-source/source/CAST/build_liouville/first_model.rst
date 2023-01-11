The first building block
========================

    *I consider that I understand an equation when I can predict the properties of its solutions, without actually solving it.*

    Paul Dirac

The essential difference between a Liouville structure and a Weinstein structure is that the dynamics of the Liouville vector field :math:`X` for the former can be much more complex than then later, which is necessarily a gradient vector field. In fact, all previously known Liouville domains that are not Weinstein came as "one piece", namely, they came out of some global constructions where no understanding of the dynamics of :math:`X` was necessary. These examples, especially their topological complexity, led to an impression that Liouville structures are genuinely global and cannot be built by simple local pieces. While this impression may well turn out to be true, we shall force the issue by constructing some local building blocks that have relatively simple topology (and dynamics). This effort grew out of the general believe that if a structure cannot be understood locally, it cannot be understood at all.

To appreciate the necessary complexity of a Liouville vector field :math:`X`, it must be pointed out that the most troublesome orbits are not the closed ones, for they can easily be "Morsified", but rather those that are "trapped": They neither close up nor leave the domain in any time. It is therefore easy to imagine that the most basic building blocks of a Liouville structure should look like a mapping torus. This idea was explored `here <https://arxiv.org/abs/1910.14132>`_, but unfortunately, Liouville structures were constructed again as "one piece" and no indication was given as for how this idea may be used to construct a (Liouville) building block. The reason for this is that I didn't know how to make mapping torus simple: The section of the mapping torus was topologically complicated and so was the monodromy. It is the goal of this article to achieve what I couldn't achieve before: the mapping torus should really just be a solid torus :math:`S^1 \times D^{2n-1}`, topologically speaking. Of course, it's crucial to also specify the "attaching region" on the boundary of the solid torus where :math:`X` is inward-pointing. This attaching region should also be topologically simple and yet have sufficiently high-dimensional topology so it cannot possibly be Weinstein. Everything shall become obvious when the construction is carried out later, but as an example, when :math:`n=2`, the attaching region will be a thickened torus.

This article is (tentatively) organized as follows. We first outline the idea of the construction, then the details will be worked out in dimension :math:`4`, and finally in all dimensions.

The Main Ideas
--------------

Our little Liouville piece will look like a mapping torus, namely, it's topologically a product :math:`Y \times [0,1]` with the two ends :math:`Y_i \coloneqq Y \times \{i\}, i=0, 1`, somehow glued together. It's common in topology to identify the ends by a (possibly partially defined) map :math:`Y_1 \to Y_0`, known as the *monodromy*. This is *not* suitable to our purposes, and we shall describe our mapping tori in a way where the monodromy is not explicit as follows.

Mapping tori can in general have complicated topology, which we are not interested in. Hence we shall for the rest of the article assume the following two conditions:

* :math:`Y` is (topologically) a ball embedded in :math:`\RR^{2n-1}`.
* The mapping torus is specified by an embedding :math:`\phi: Y_1 \to \RR^{2n-1}`. Here the target :math:`\RR^{2n-1}` contains :math:`Y_0`. Note that the map :math:`\phi` is not, but can be translated into, the monodromy. [#dumb_tori]_

So far the discussions are completely topological. To bring the Liouville structure into the picture, we shall postulate another two conditions:

* The Liouville vector field :math:`X` is parallel to :math:`\p_s` where :math:`s \in [0,1]`.
* The contact structure on :math:`\RR^{2n-1}` is standard.

Since the contact form grows in the direction of :math:`X`, it's necessary that the contact form on :math:`Y_1` is "greater than" that on :math:`Y_0`, at least over the overlap. In other words, if :math:`\alpha` is the (standard) contact form on :math:`\RR^{2n-1}`, then :math:`\phi^\ast (\alpha) = g \alpha` where :math:`g \geq 1`.

As probably the simplest example, let :math:`n=2`, :math:`\alpha = dz-ydx` defining the standard contact structure on :math:`\RR^3`, and :math:`Y_0 = [-1, 1]^3` be a box around the origin. Consider :math:`\phi(x, y, z) = (2x, 2y, 4z)` which expands the box. Then :math:`\phi^\ast (\alpha) = 4\alpha > \alpha` as desired. However, this :math:`T_\phi` is nothing interesting. Indeed, if we tilt the side :math:`\p Y \times [0, 1]`, which is by construction tangent to :math:`X`, so that :math:`X` is outward-pointing, then :math:`T_\phi = S^1 \times D^3` is really just a Weinstein domain with one :math:`0`-handle and one :math:`1`-handle in the obvious way. In particular, it cannot be attached to anything because it has no negative boundary -- the region along which :math:`X` is inward-pointing. [#wrong_tilt]_

Now let's consider another :math:`\phi(x, y, z) = (4x, y/2, 2z)` which expands the box in the :math:`x,z`-directions and contracts in the :math:`y`-direction. Then :math:`\phi^\ast (\alpha) = 2\alpha` also checks out. This case is a bit more interesting because at the :math:`\RR^3`-level, :math:`X` is outward-pointing along :math:`\phi(Y_1) \setminus Y_0 \cong S^1 \times D^2` and inward-pointing along :math:`Y_0 \setminus \phi(Y_1) \cong D^3 \sqcup D^3`.

tbc

.. rubric:: Footnotes

.. [#dumb_tori] It's, strictly speaking, possible that :math:`\phi(Y_1)` completely misses :math:`Y_0`, in which case one may not wish to call the result a mapping torus. But this is just a dumb case that we'll never encounter.

.. [#wrong_tilt] Even if we tilt the side :math:`\p Y \times [0, 1]` in the other direction so that :math:`X` is inward-pointing, we still wouldn't get a handle-like :math:`T_\phi` because the tangency of :math:`X` along :math:`\p T_\phi` has fold singularities in both directions.
