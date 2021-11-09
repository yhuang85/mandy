Rigidity and Flexibility of CoLegendrians in dimension :math:`5`
================================================================

The primary goal of this article is to explore the basic properties of coLegendrians in contact :math:`5`-manifolds, including how they look like, how to find them and how to manipulate them. Ideally, we'd also like to see how they interact with Legendrians, e.g., how they may play a role in the classification of Legendrians. We'll surely touch upon this point, but at this moment, I don't exactly know how far it goes.

Everything we shall talk about in this article belong to the realm of contact Morse theory. So let's start by introducing the major roles and the context. Since contact Morse theory is the one stage on which all stories will be told and we'll not be talking about the usual (i.e., topological) Morse theory, we often drop "contact" from the theory to be less wordy.

(Co)Legendrians in Morse theory
-------------------------------

Given a contact :math:`5`-manifold :math:`(M, \xi)`, a *Legendrian* :math:`\Lambda \subset M` is a submanifold such that :math:`T_x\Lambda \subset \xi_x` for all :math:`x \in M`. For now, let's assume :math:`\Lambda` is smooth, but as E. Murphy argued in [Mur12]_, it may become necessary, or at least beneficial, to consider also (mildly) singular submanifolds, in which case one has to, of course, make sense of :math:`T\Lambda` at singularities.

Similarly, a :math:`3`-submanifold :math:`Y \subset M` is a *coLegendrian* if :math:`T_x Y \cap \xi_x \subset \xi_x` is coisotropic for all :math:`x \in Y` with respect to the natural (linear) symplectic structure on :math:`\xi`. More concretely, there are two scenarios here. Namely, when the intersection is transverse, :math:`T_x Y \cap \xi_x \subset \xi_x` is Lagrangian, otherwise :math:`T_x Y = T_x Y \cap \xi_x \subset \xi_x` is coisotropic. A geometrically intuitive, but possibly also useless, way to think about such :math:`Y` is that it naturally comes with a singular codimension-:math:`1` foliation whose leaves are Legendrians and the singularities are precisely where the transversality fails. As in the case of Legendrians, we'll also encounter coLegendrians which are not everywhere smooth.

.. note::

    None of the submanifolds considered in this article are parametrized, i.e., they are regarded as subsets of, say, :math:`M` rather than maps into :math:`M`.

tbc...

.. rubric:: References

.. [Mur12] E\. Murphy\. `Loose Legendrian embeddings in high dimensional contact manifolds <https://arxiv.org/abs/1201.2245v5>`_
