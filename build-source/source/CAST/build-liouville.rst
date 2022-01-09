How to build Liouville domains?
===============================

A *Liouville domain* :math:`(W, \omega)` is by definition a compact exact symplectic manifold with convex boundary. Namely, we demand the existence of a *Liouville form* :math:`\lambda` such that :math:`\omega = d\lambda` and a *Liouville vector field* :math:`X` -- defined by :math:`\iota_X \omega = \lambda` -- such that :math:`X` is everywhere outward pointing along :math:`\p W`. The question of which (compact with boundary) manifolds can carry a Liouville structure seems difficult, but there are obvious and not-so-obvious conditions that :math:`W` must satisfy to become a Liouville domain. For example, one obvious condition is that :math:`W` must carry an almost complex structure, and one not-so-obvious condition is that :math:`W` must not be :math:`S^3 \times [0,1]`, even though it does carry an almost complex structure. The later constraint comes from some rather complex geometric analysis (e.g. pseudo-holomorphic curves or mathematical gauge theory) which works only when :math:`\dim W = 4`. In any case, there is no chance of answering the question if one cannot decompose a Liouville domain into "simple" pieces.

There is one noteworthy special case though when :math:`(W, \omega, \lambda, X)` carries one additional assumption that :math:`X` is gradient-like for a Morse function :math:`f: W \to \RR`. In this case :math:`f` induces a handle decomposition of :math:`W` where the handles come with standard symplectic structures. These manifolds are known as *Weinstein domains* -- named after A. Weinstein. On the one hand, the construction of Weinstein provides us an interesting class of Liouville domains, but on the other hand, the Morse condition severely constrains the topology of :math:`W`, e.g., the Morse index of any critical point of :math:`f` cannot exceed :math:`\tfrac{1}{2} \dim W`.

One crucial difference between Weinstein domains and Liouville domains in general is the dynamics of :math:`X`. For Weinstein domains, :math:`X` has no (broken) periodic orbits and :math:`W` can be gradually built from sub-level sets of the associated Morse function :math:`f`, which doesn't exist for a general Liouville domain. It should nevertheless be noted that it is not the presence of periodic orbits of :math:`X` that distinguishes the Liouville domains that are not Weinstein -- it is the orbits that neither close up nor leave :math:`W` in any finite time. These "trapped" orbits must be cut into pieces, and our strategy to achieve this is to, roughly speaking, "stop at the first time return". This ends up cutting :math:`W` into a collections of (partial) mapping tori. The plan is therefore to understand/classify these mapping tori and hope that they will be suffice to build just about any Liouville domain.

.. note::

   The idea of relating Liouville domains to mapping tori was first explored in [Hua19]_, but the treatment there covers only the useless case: the mapping tori that have only positive ends and no negative ends -- they are themselves Liouville domains but cannot be assembled together.

A Toy Model
-----------

It's hard to imagine, out of nowhere, how the (Liouville) mapping tori should look like, not to say how they fit together. So we shall start with a concrete example, which was first studied by Y. Mitsumatsu [Mit95]_, and examine its anatomy.

.. rubric:: References

.. [Hua19] Y\. Huang\. `A dynamical construction of Liouville domains <https://arxiv.org/abs/1910.14132v2>`_

.. [Mit95] Y\. Mitsumatsu\. `Anosov flows and non-Stein symplectic manifolds <http://www.numdam.org/item/AIF_1995__45_5_1407_0>`_
