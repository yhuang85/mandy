.. _coLeg_5d_rst:

Rigidity and flexibility of coLegendrians in dimension :math:`5`
================================================================

The primary goal of this article is to explore the basic properties of coLegendrians in contact :math:`5`-manifolds, including how they look like, how to find them and how to manipulate them. Ideally, we'd also like to see how they interact with Legendrians, e.g., how they may play a role in the classification of Legendrians. We'll surely touch upon this point, but at this moment, I don't exactly know how far it goes.

Everything we shall talk about in this article belong to the realm of :ref:`contact Morse theory <contact_morse_theory_rst>`. So let's start by introducing the major roles and the context. Since contact Morse theory is the one stage on which all stories will be told and we'll not be talking about the usual (i.e., topological) Morse theory, we often drop "contact" from the theory to be less wordy.

Morse-theoretic Setup
---------------------

Given a contact :math:`5`-manifold :math:`(M, \xi)`, a *Legendrian* :math:`\Lambda \subset M` is a submanifold such that :math:`T_x\Lambda \subset \xi_x` is Lagrangian for all :math:`x \in \Lambda` with respect to the natural (linear) symplectic structure on :math:`\xi`. For now, let's assume :math:`\Lambda` is smooth, but as E. Murphy argued in [Mur12]_, it may become necessary, or at least beneficial, to consider also (mildly) singular submanifolds, in which case one has to, of course, make sense of :math:`T\Lambda` at singularities.

Similarly, a :math:`3`-submanifold :math:`Y \subset M` is a *coLegendrian* if :math:`T_x Y \cap \xi_x \subset \xi_x` is coisotropic for all :math:`x \in Y`. More concretely, there are two scenarios here. Namely, when the intersection is transverse, :math:`T_x Y \cap \xi_x \subset \xi_x` is Lagrangian, otherwise :math:`T_x Y \subset \xi_x` is coisotropic. A geometrically intuitive, but possibly also useless, way to think about such :math:`Y` is that it naturally comes with a singular codimension-:math:`1` foliation whose leaves are Legendrians and the singularities are precisely where the transversality fails. As in the case of Legendrians, we'll also encounter coLegendrians which are not everywhere smooth.

.. note::

    None of the submanifolds considered in this article are parametrized, i.e., they are regarded as subsets of, say, :math:`M` rather than maps into :math:`M`.

To put everything in the context of Morse theory, we shall make the following assumption.

.. _morse_assumption:

    All Legendrians :math:`\Lambda` are assumed to be contained in an ambient Morse hypersurface :math:`\Sigma` such that :math:`\Lambda_{\xi} \coloneqq \Sigma_{\xi}|_{\Lambda}` is tangent to :math:`\Lambda` and restricts to a Morse vector field on :math:`\Lambda`. The same holds for all coLegendrians.

A few remarks are in order regarding the above assumption. First of all, we don't in general care about the ambient :math:`\Sigma`, so it generally won't even be mentioned in what follows. Secondly, the above assumption fit in the context of family Morse theory by thinking of a coLegendrian :math:`Y = \Lambda \times [0, 1]` as a "discretized" isotopy from :math:`\Lambda_0` to :math:`\Lambda_1`, such that :math:`(\Lambda_0)_{\xi}` and :math:`(\Lambda_1)_{\xi}` are allowed to have birth-death type critical points, aside from the nondegenerate ones. Finally, such Morse data on (co)Legendrians amounts to additional choices which we shall always carry along. For example, a Legendrian may be contained in different hypersurfaces and hence inherits different Morse structures. They shall be considered different even through it's the very same Legendrian. On the other hand, being a Morse coLegendrian demands, in particular, trivial normal bundle, which is not necessarily satisfied by all coLegendrians. As a result, we simply don't consider those with nontrivial normal bundles.

Recall that a Morse hypersurface :math:`\Sigma` admits a decomposition :math:`\Sigma = R_+(\Sigma) \cup_{\Gamma(\Sigma)} R_-(\Sigma)` such that :math:`R_{\pm} (\Sigma)` are built out of the positive/negative handles, and :math:`\Gamma(\Sigma)` is a level set in the middle. A good thing about (co)Legendrians satisfying the above :ref:`assumption <morse_assumption>` is that they inherit similar decompositions. For example, a Legendrian :math:`\Lambda` can be decomposed as :math:`\Lambda = R_+(\Lambda) \cup_{\Gamma(\Lambda)} R_-(\Lambda)` such that :math:`R_{\pm} (\Lambda)` are built out of the positive/negative handles given by :math:`\Lambda_{\xi} = \Sigma_{\xi}|_{\Lambda}`. In particular, :math:`\Gamma(\Lambda) \subset \Gamma(\Sigma)` is in general a Legendrian link inside a contact :math:`3`-manifold. Once again, let's emphasize that such a decomposition depends on the embedding :math:`\Lambda \subset \Sigma`. Similarly, for a coLegendrian :math:`Y`, :math:`\Gamma(Y) \subset \Gamma(\Sigma)` is a (hyper)surface which we will always arrange/assume to be Morse. In other words, we will work within an environment such that the Morse theory descends to lower-dimensional submanifolds in a recursive manner. We will call any of the above :math:`\Gamma` a *dividing set*, so a dividing set can be a contact manifold, a hypersurface or a Legendrian, depending on the context.

Flexibility
-----------

From the perspective of the usual *h*-principle established by M. Gromov in [Gro86]_, one wouldn't expect to be able to approximate a smoothly embedded :math:`3`-submanifold :math:`Y \subset M` (say with trivial normal bundle) by coLegendrian, simply because :math:`\dim Y > \tfrac{1}{2} \dim M`. However, whether such an expectation is correct or not doesn't really matter because a general coLegendrian is as useless as a general (say, convex in the sense of E. Giroux) hypersurface due to intractability.

It's therefore quite a pleasing fact that by applying the folding techniques developed in [HH19]_, one can argue that any :math:`Y \subset M` with trivial normal bundle can be :math:`C^0`-approximated by a (homeomorphic) coLegendrian, which we still denote by :math:`Y`. However, the so constructed :math:`Y` is a topological approximation which is not necessarily smoothly embedded. Indeed, suppose :math:`\Sigma` is the (Morse) hypersurface containing :math:`Y`, then :math:`Y_{\xi} \coloneqq \Sigma_{\xi}|_Y` is tangent to :math:`Y` in the sense that any flow line of :math:`\Sigma_{\xi}` which intersects :math:`Y` is completely contained in :math:`Y`, and is itself a Morse vector field on :math:`Y`. One can verify that :math:`Y` indeed satisfies the coLegendrian condition where it's smooth.

.. sidebar:: An coLegendrian cone

    .. image:: static/cone.svg
        :align: center
        :width: 200px

Since :math:`Y \subset \Sigma` is tangent to :math:`\Sigma_{\xi}`, the singularities of :math:`Y` are necessarily (families of) cones. A convenient, but also coincidental, consequence of the assumption :math:`\dim M = 5` is that by further wiggling :math:`\Sigma`, one can approximate :math:`Y` by a coLegendrian with only isolated cone singularities at index :math:`0` and :math:`3` (when viewed inside :math:`Y`) critical points. More precisely, those are cones over (smooth and Morse) :math:`2`-spheres in :math:`(S^3, \xi_{\std}) = \p (B^4, \omega_{std})`. A schematic picture of such a cone is drawn on the right-hand-side, where :math:`O \in B^4` is the origin, and the shaded cone represents (part of) :math:`Y`. The procedure of simplifying singularities on :math:`Y` is explained in [Hua20]_.

.. admonition:: Todo
    :class: warning

    The structural theory of coLegendrians is so far only developed for contact :math:`5`-manifolds because (1) I'm out of time, and (2) its development may require a more thorough understanding of :math:`5`-dimensional contact topology due to the recursive nature of contact Morse theory. Nonetheless, it's expected that higher dimensional coLegendrians possess singularities more complex than just the isolated cones.

The above discussion applies to both the closed case where :math:`\p Y = \varnothing` and the case where :math:`\Lambda \coloneqq \p Y` is a Legendrian. In the later case, we also allow :math:`\Lambda` to possess cone singularities, and we'll explicitly say that :math:`\Lambda` is smooth otherwise.

Building Blocks
---------------

Now that the flexibility discussed above guarantees the existence of coLegendrians, let's move on to understand the inner structure of :math:`Y`. This is rather straightforward since :math:`Y_{\xi}` defines a handle decomposition of :math:`Y` and the handles are slices of standard Weinstein handles.

tbc...

.. rubric:: References

.. [Gro86] M\. Gromov\. `Partial differential relations <https://www.ihes.fr/~gromov/wp-content/uploads/2018/08/248.pdf>`_

.. [HH19] K\. Honda and Y\. Huang\. `Convex hypersurface theory in contact topology <https://arxiv.org/abs/1907.06025v2>`_

.. [Hua20] Y\. Huang\. `Existence of coLegendrians in contact 5-manifolds <https://arxiv.org/abs/2006.11844v1>`_

.. [Mur12] E\. Murphy\. `Loose Legendrian embeddings in high dimensional contact manifolds <https://arxiv.org/abs/1201.2245v5>`_
