Rigidity and Flexibility of coLegendrians in dimension :math:`5`
================================================================

The primary goal of this article is to explore the basic properties of coLegendrians in contact :math:`5`-manifolds, including how they look like, how to find them and how to manipulate them. Ideally, we'd also like to see how they interact with Legendrians, e.g., how they may play a role in the classification of Legendrians. We'll surely touch upon this point, but at this moment, I don't exactly know how far it goes.

Everything we shall talk about in this article belong to the realm of contact Morse theory. So let's start by introducing the major roles and the context. Since contact Morse theory is the one stage on which all stories will be told and we'll not be talking about the usual (i.e., topological) Morse theory, we often drop "contact" from the theory to be less wordy.

(Co)Legendrians in Morse theory
-------------------------------

Given a contact :math:`5`-manifold :math:`(M, \xi)`, a *Legendrian* :math:`\Lambda \subset M` is a submanifold such that :math:`T_x\Lambda \subset \xi_x` for all :math:`x \in M`. For now, let's assume :math:`\Lambda` is smooth, but as E. Murphy argued in [Mur12]_, it may become necessary, or at least beneficial, to consider also (mildly) singular submanifolds, in which case one has to, of course, make sense of :math:`T\Lambda` at singularities.

Similarly, a :math:`3`-submanifold :math:`Y \subset M` is a *coLegendrian* if :math:`T_x Y \cap \xi_x \subset \xi_x` is coisotropic for all :math:`x \in Y` with respect to the natural (linear) symplectic structure on :math:`\xi`. More concretely, there are two scenarios here. Namely, when the intersection is transverse, :math:`T_x Y \cap \xi_x \subset \xi_x` is Lagrangian, otherwise :math:`T_x Y \subset \xi_x` is coisotropic. A geometrically intuitive, but possibly also useless, way to think about such :math:`Y` is that it naturally comes with a singular codimension-:math:`1` foliation whose leaves are Legendrians and the singularities are precisely where the transversality fails. As in the case of Legendrians, we'll also encounter coLegendrians which are not everywhere smooth.

.. note::

    None of the submanifolds considered in this article are parametrized, i.e., they are regarded as subsets of, say, :math:`M` rather than maps into :math:`M`.

To put everything in the context of Morse theory, we shall make the following assumption.

.. _morse_assumption:

    All Legendrians :math:`\Lambda` are assumed to be contained in an ambient Morse hypersurface :math:`\Sigma` such that :math:`\Lambda_{\xi} \coloneqq \Sigma_{\xi}|_{\Lambda}` is tangent to :math:`\Lambda` and restricts to a Morse vector field on :math:`\Lambda`. The same holds for all coLegendrians.

A few remarks are in order regarding the above assumption. First of all, we don't in general care about the ambient :math:`\Sigma`, so it generally won't even be mentioned in what follows. Secondly, the above assumption fit in the context of family Morse theory by thinking of a coLegendrian :math:`Y = \Lambda \times [0, 1]` as a "discretized" isotopy from :math:`\Lambda_0` to :math:`\Lambda_1`, such that :math:`(\Lambda_0)_{\xi}` and :math:`(\Lambda_1)_{\xi}` are allowed to have birth-death type critical points, aside from the nondegenerate ones. Finally, such Morse data on (co)Legendrians amounts to additional choices which we shall always carry along. For example, a Legendrian may be contained in different hypersurfaces and hence inherits different Morse structures. They shall be considered different even through it's the very same Legendrian. On the other hand, being a Morse coLegendrian demands, in particular, trivial normal bundle, which is not necessarily satisfied by all coLegendrians. As a result, we simply don't consider those with nontrivial normal bundles.

Recall that a Morse hypersurface :math:`\Sigma` admits a decomposition :math:`\Sigma = R_+(\Sigma) \cup_{\Gamma(\Sigma)} R_-(\Sigma)` such that :math:`R_{\pm} (\Sigma)` are built out of the positive/negative handles, and :math:`\Gamma(\Sigma)` is a level set in the middle. A good thing about (co)Legendrians satisfying the above :ref:`assumption <morse_assumption>` is that they inherit similar decompositions. For example, a Legendrian :math:`\Lambda` can be decomposed as :math:`\Lambda = R_+(\Lambda) \cup_{\Gamma(\Lambda)} R_-(\Lambda)` such that :math:`R_{\pm} (\Lambda)` are built out of the positive/negative handles given by :math:`\Lambda_{\xi} = \Sigma_{\xi}|_{\Lambda}`. In particular, :math:`\Gamma(\Lambda) \subset \Gamma(\Sigma)` is in general a Legendrian link inside a contact :math:`3`-manifold. Once again, let's emphasize that such a decomposition depends on the embedding :math:`\Lambda \subset \Sigma`. Similarly, for a coLegendrian :math:`Y`, :math:`\Gamma(Y) \subset \Gamma(\Sigma)` is a surface which we will always arrange/assume to be Morse. In other words, we will work within an environment such that the Morse theory descends to lower-dimensional submanifolds in a recursive manner. We will call any of the above :math:`\Gamma` a *dividing set*, so a dividing set can be a contact manifold, a hypersurface or a Legendrian, depending on the context.

Flexibility of coLegendrians
----------------------------

tbc...

.. rubric:: References

.. [Mur12] E\. Murphy\. `Loose Legendrian embeddings in high dimensional contact manifolds <https://arxiv.org/abs/1201.2245v5>`_
