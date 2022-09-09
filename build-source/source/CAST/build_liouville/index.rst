How to build Liouville domains?
===============================

A *Liouville domain* :math:`(W, \omega)` is by definition a compact exact symplectic manifold with convex boundary. Namely, we demand the existence of a *Liouville form* :math:`\lambda` such that :math:`\omega = d\lambda` and a *Liouville vector field* :math:`X` -- defined by :math:`\iota_X \omega = \lambda` -- such that :math:`X` is everywhere outward pointing along :math:`\p W`. The question we will try to answer here is simply

    Which (compact with boundary) manifolds carry a Liouville structure?

There are obvious and not-so-obvious conditions that :math:`W` must satisfy to become a Liouville domain. An obvious and general condition is that :math:`W` must carry an almost complex structure. On the other hand, there are not-so-obvious and rather ad hoc conditions such as that :math:`W` must not be :math:`S^3 \times [0,1]`, even though it does carry an almost complex structure. The later constraint comes from some rather complicated geometric analysis (e.g. pseudo-holomorphic curves or mathematical gauge theory) which works only in dimension :math:`4`.

There is one noteworthy special case when :math:`(W, \omega, \lambda, X)` carries one additional assumption that :math:`X` is gradient-like for a Morse function :math:`f: W \to \RR`. In this case :math:`f` induces a handle decomposition of :math:`W` where the handles come with standard symplectic structures. Such a symplectic structure on :math:`W` is called a *Weinstein structure* and such :math:`W` is called a  *Weinstein domain*. Now the good news is that the existence of Weinstein structures are completely characterized by the work of Eliashberg [Eli90]_:

    :math:`W` carries a Weinstein structure if and only if it is almost complex and has no topology above half of its dimension.

The bad news, however, is also obvious from the above statement: The severe constraint on the topology of :math:`W` is by no means necessary for the general existence of Liouville structures.

It is the goal of this series to find out what is really needed to build Liouville domains.

.. toctree::
    :titlesonly:
    :maxdepth: 2

    first_model

.. rubric:: References

.. [Eli90] Y\. Eliashberg\. Topological characterization of Stein manifolds of dimension > 2
