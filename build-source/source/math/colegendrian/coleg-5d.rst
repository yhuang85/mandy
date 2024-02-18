.. _coleg_5d_rst:

CoLegendrians in dimension five
===============================

.. note::

    Work in progress

The primary goal of this article is to explore the basic properties of coLegendrians in contact :math:`5`-manifolds, including how they look like, how to find them and how to manipulate them. Ideally, we'd also like to see how they interact with Legendrians, e.g., how they may play a role in the classification of Legendrians. We'll surely touch upon this point, but at this moment, I don't exactly know how far it goes.

Everything we shall talk about in this article belong to the realm of :ref:`contact Morse theory <contact_morse_theory_rst>`. So let's start by introducing the major roles and the context. Since contact Morse theory is the one stage on which all stories will be told and we'll not be talking about the usual (i.e., topological) Morse theory, we often drop "contact" from the theory to be less wordy.

Morse-theoretic Setup
---------------------

Given a contact :math:`5`-manifold :math:`(M, \xi)`, a *Legendrian* :math:`\Lambda \subset M` is a submanifold such that :math:`T_x\Lambda \subset \xi_x` is Lagrangian for all :math:`x \in \Lambda` with respect to the natural (linear) symplectic structure on :math:`\xi`. For now, let's assume :math:`\Lambda` is smooth, but as argued by E. Murphy in [Mur12]_, it may become necessary, or at least beneficial, to consider also (mildly) singular submanifolds, in which case one has to, of course, make sense of :math:`T\Lambda` at singularities.

Similarly, a :math:`3`-submanifold :math:`Y \subset M` is a *coLegendrian* if :math:`T_x Y \cap \xi_x \subset \xi_x` is coisotropic for all :math:`x \in Y`. More concretely, there are two scenarios here. Namely, when the intersection is transverse, :math:`T_x Y \cap \xi_x \subset \xi_x` is Lagrangian, otherwise :math:`T_x Y \subset \xi_x` is coisotropic. A geometrically intuitive, but possibly also useless, way to think about such :math:`Y` is that it naturally comes with a singular codimension-:math:`1` foliation whose leaves are Legendrians and the singularities are precisely where the transversality fails. As in the case of Legendrians, we'll also encounter coLegendrians which are not everywhere smooth.

.. note::

    None of the submanifolds considered in this article are parametrized, i.e., they are regarded as subsets of, say, :math:`M` rather than maps into :math:`M`.

To put everything in the context of Morse theory, we shall make the following assumption.

.. _morse_assumption:

    All Legendrians :math:`\Lambda` are assumed to be contained in an ambient Morse hypersurface :math:`\Sigma` such that :math:`\Lambda_{\xi} \coloneqq \Sigma_{\xi}|_{\Lambda}` is tangent to :math:`\Lambda` and restricts to a Morse vector field on :math:`\Lambda`. The same holds for all coLegendrians.

A few remarks are in order regarding the above assumption. First of all, we don't in general care about the ambient :math:`\Sigma`, so it generally won't even be mentioned in what follows. Secondly, the above assumption fit in the context of family Morse theory by thinking of a coLegendrian :math:`Y = \Lambda \times [0, 1]` as a "discretized" isotopy from :math:`\Lambda_0` to :math:`\Lambda_1`, such that :math:`(\Lambda_0)_{\xi}` and :math:`(\Lambda_1)_{\xi}` are allowed to have birth-death type critical points, aside from the nondegenerate ones. Finally, such Morse data on (co)Legendrians amounts to additional choices which we shall always carry along. For example, a Legendrian may be contained in different hypersurfaces and hence inherits different Morse structures. They shall be considered different even through it's the very same Legendrian. On the other hand, being a Morse coLegendrian demands, in particular, trivial normal bundle, which is not necessarily satisfied by all coLegendrians. As a result, we simply don't consider those with nontrivial normal bundles.

Recall that a Morse hypersurface :math:`\Sigma` admits a decomposition :math:`\Sigma = R_+(\Sigma) \cup_{\Gamma(\Sigma)} R_-(\Sigma)` such that :math:`R_{\pm} (\Sigma)` are built out of the positive/negative handles, and :math:`\Gamma(\Sigma)` is a level set in the middle. A good thing about (co)Legendrians satisfying the above :ref:`assumption <morse_assumption>` is that they inherit similar decompositions. For example, a Legendrian :math:`\Lambda` can be decomposed as :math:`\Lambda = R_+(\Lambda) \cup_{\Gamma(\Lambda)} R_-(\Lambda)` such that :math:`R_{\pm} (\Lambda)` are built out of the positive/negative handles given by :math:`\Lambda_{\xi} = \Sigma_{\xi}|_{\Lambda}`. In particular, :math:`\Gamma(\Lambda) \subset \Gamma(\Sigma)` is in general a Legendrian link inside a contact :math:`3`-manifold. Let's emphasize once again that such a decomposition depends on the embedding :math:`\Lambda \subset \Sigma`. Similarly, for a coLegendrian :math:`Y`, :math:`\Gamma(Y) \subset \Gamma(\Sigma)` is a surface, but in general a coLegendrian itself, which we will always arrange/assume to be Morse. In other words, we will work within an environment such that the Morse theory descends to lower-dimensional submanifolds in a recursive manner. We will call any of the above :math:`\Gamma` a *dividing set*, so a dividing set can be a contact manifold, a coLegendrian or a Legendrian, depending on the context.

Flexibility
-----------

From the perspective of the usual *h*-principle established by M. Gromov in [Gro86]_, one wouldn't expect to be able to approximate a smoothly embedded :math:`3`-submanifold :math:`Y \subset M` (say with trivial normal bundle) by coLegendrian, simply because :math:`\dim Y > \tfrac{1}{2} \dim M`. However, whether such an expectation is correct or not doesn't really matter because a general coLegendrian is as useless as a general (say, convex in the sense of E. Giroux) hypersurface due to intractability.

It's therefore quite a pleasing fact that by applying the folding techniques developed in [HH19]_, one can argue that any :math:`Y \subset M` with trivial normal bundle can be :math:`C^0`-approximated by a (homeomorphic) coLegendrian, which we still denote by :math:`Y`. However, the so constructed :math:`Y` is a topological approximation which is not necessarily smoothly embedded. Indeed, suppose :math:`\Sigma` is the (Morse) hypersurface containing :math:`Y`, then :math:`Y_{\xi} \coloneqq \Sigma_{\xi}|_Y` is tangent to :math:`Y` in the sense that any flow line of :math:`\Sigma_{\xi}` which intersects :math:`Y` is completely contained in :math:`Y`, and is itself a Morse vector field on :math:`Y`. One can verify that :math:`Y` indeed satisfies the coLegendrian condition wherever it's smooth.

.. _coLeg_cone:

.. sidebar:: An coLegendrian cone

    .. figure:: static/coleg-5d/cone.svg
        :align: center
        :width: 150px

Since :math:`Y \subset \Sigma` is tangent to :math:`\Sigma_{\xi}`, the singularities of :math:`Y` are necessarily (families of) cones. A convenient, but also coincidental, consequence of the assumption :math:`\dim M = 5` is that by further wiggling :math:`\Sigma`, one can approximate :math:`Y` by a coLegendrian with only isolated cone singularities of index :math:`0` and :math:`3` (when viewed inside :math:`Y`) critical points. More precisely, those are cones over (smooth and Morse) :math:`2`-spheres in :math:`(S^3, \xi_{\std}) = \p (B^4, \omega_{std})`. A schematic picture of such a cone is drawn on the right-hand-side, where :math:`O \in B^4` is the origin, and the shaded cone represents (part of) :math:`Y`. The procedure of simplifying singularities on :math:`Y` is explained in [Hua20]_.

.. admonition:: TODO
    :class: warning

    The structural theory of coLegendrians is so far only developed for contact :math:`5`-manifolds because (1) I'm out of time, and (2) its development may require a more thorough understanding of :math:`5`-dimensional contact topology due to the recursive nature of contact Morse theory. Nonetheless, it's expected that higher dimensional coLegendrians possess singularities more complex than just the isolated cones.

The above discussion applies to both the closed case where :math:`\p Y = \varnothing` and the case where :math:`\Lambda \coloneqq \p Y` is a Legendrian. In the later case, we also allow :math:`\Lambda` to possess cone singularities, and we'll explicitly say that :math:`\Lambda` is smooth otherwise.

.. _coleg_5d_building_blocks:

Building Blocks
---------------

Now that the flexibility discussed above guarantees the existence of coLegendrians, let's move on to understand the inner structure of :math:`Y`. This is rather straightforward since :math:`Y_{\xi}` defines a handle decomposition of :math:`Y` and the handles are slices of standard Weinstein handles. However, since we're dealing with a "nested" Morse theory, it's necessary to differentiate Morse indices with respect to the ambient manifold. For example, one same critical point :math:`x \in \Lambda \subset Y \subset \Sigma` might have index :math:`0` when considered in :math:`\Lambda`, and index :math:`1` in :math:`Y`, and index :math:`2` in :math:`\Sigma`. This is reflected in the following notations: :math:`\ind_{\Lambda} (x) = 0, \ind_Y (x) = 1`, and :math:`\ind_{\Sigma} (x) = 2`.

.. note::

    In addition to showing that any :math:`Y \subset \Sigma` can be approximated by a coLegendrian, it's further explained in [Hua20]_ that such an approximation can either be made repelling in the normal direction, i.e., :math:`\ind_Y (x) = \ind_{\Sigma} (x)` for all critical points :math:`x \in Y`, or be made attracting in the normal direction, i.e., :math:`\ind_Y (x) + 1 = \ind_{\Sigma} (x)` for all critical points :math:`x \in Y`. These additional arrangements may be convenient sometimes, but we decide not to bake them into the initial :ref:`assumptions <morse_assumption>` of coLegendrians to allow more flexibility.

Since coLegendrians are the main objects of interest in the article, we'll implicitly assume unspecified indices are :math:`Y`-indices and specify :math:`\Lambda` and :math:`\Sigma`-indices as needed in what follows.

Since :math:`\dim Y = 3`, it's built out of :math:`0,1,2`, and :math:`3`-handle. Assume :math:`\partial Y = \varnothing` for the time being. The :math:`0` and :math:`3`-handles are :ref:`cones <coLeg_cone>` over (Morse) :math:`2`-spheres in :math:`(S^3, \xi_{std})`. It remains to describe the :math:`1` and :math:`2`-handles, in both signs. Let's describe the positive :math:`1`-handle :math:`H_1^+` and :math:`2`-handle :math:`H_2^+`, and note that the negative :math:`H_1^-` and :math:`H_2^-` are dual to :math:`H_2^+` and :math:`H_1^+` (by reversing the Morse vector fields), respectively.

.. sidebar:: coLegendrian handles :math:`H_1^+` and :math:`H_2^+`

    .. figure:: static/coleg-5d/one-and-two-handles-5d.svg
        :align: center
        :width: 400px

Recall that :math:`H_1^+` and :math:`H_2^+` are slices of Weinstein :math:`1` and :math:`2`-handles, respectively. While :math:`H_2^+` is essentially unambiguous, there are choices one can make when defining :math:`H_1^+` as a :math:`3`-d slice of a Weinstein :math:`1`-handle. It turns out that one only needs the ones depicted on the right to build any coLegendrian. In the picture, the coLegendrian handles are drawn as solid cylinders, and the attaching regions, i.e., where the Morse vector fields point inward, are colored in green. The characteristic foliations on :math:`\p_{\pm} H` are also drawn, where :math:`\p_- H` are the attaching regions. In particular :math:`H_1^+` is attached along a pair of disks surrounding a source and a sink, respectively, and :math:`H_2^+` is a attached along an annulus surrounding a pair of positive and negative saddles :math:`h_{\pm}` together with the two separatrices flowing from :math:`h_-` to :math:`h_+`.

By reversing the Morse vector fields, :math:`H_1^-` (dual to :math:`H_2^+`) is attached along a pair of disks with linear foliation, and :math:`H_2^-` (dual to :math:`H_1^+`) is attached along an annular neighborhood of a transverse loop.

.. important::

    Note that :math:`H_2^+` is attached long a Legendrian loop with :math:`\op{tb} = 1`. It implies, for example, that :math:`H_2^+` cannot be attached directly to a :math:`0`-handle without passing through any :math:`1`-handles since no :math:`S^2 \subset (S^3, \xi_{\std})` can contain a :math:`\op{tb} = 1` loop. This is one of the main sources of the rigidity phenomenon we shall explain below.

Finally, let's briefly comment on the case where :math:`\p Y` is a Legendrian. In this case, one has to include also "halves" of the above listed coLegendrian handles. For example, if an index :math:`0` critical point :math:`p_0 \in \p Y`, then a neighborhood of :math:`p_0` is modeled on a cone over a disk :math:`D^2 \subset (S^3, \xi_{\std})` with Legendrian boundary :math:`\p D^2`. Moreover, :math:`\Lambda = \p Y` is smooth at :math:`p_0` if and only if :math:`\p D^2` is the standard unknot.

Rigidity
--------

Topologically speaking, a coLegendrian :math:`Y` is a :math:`3`-manifold equipped with a Morse function whose critical points are signed so that the indexes of positive critical points are at most :math:`2` and the indexes of the negative ones are at least :math:`1`. Moreover, since the ambient :math:`\Sigma_{\xi}` is always assumed to be generic, we assume in addition that there is no flow lines in :math:`Y` from a negative critical point to a positive one. This assumption allows us to define the dividing set :math:`\Gamma_Y` as a level set separating the critical points of opposite signs. At this point, :math:`\Gamma_Y \subset Y` can be just about any embedded surface. This turns out to be *not* the case due to constraints imposed by the contact structure, and we call this phenomenon the rigidity of coLegendrians.

Rigidity of coLegendrians is far from being understood, letting alone applications towards better understanding contact manifolds themselves. So let's start with probably the most obvious question: Can :math:`\Gamma_Y` be empty? It turns out that the answer is always no since any Morse vector field on :math:`Y` must have at least one positive source and one negative sink (in the case :math:`\p Y \neq \varnothing`, the sources and/or the sinks may be partial). However, the nonemptiness of :math:`\Gamma_Y` is not really a rigidity phenomenon as it has nothing to do with contact structures. Our modest goal here is to show the following **connectedness rigidity**:

    If :math:`Y` is a closed coLegendrian, then :math:`\Gamma_Y` is connected.

and the **sphere rigidity**:

    If :math:`Y` is a closed coLegendrian such that :math:`\Gamma_Y \cong S^2`, then :math:`Y \cong S^3` and is standard up to Weinstein homotopy, i.e., the Morse vector field on :math:`Y` has two critical points: a (positive) maximum and a (negative) minimum.

.. note::

    The analogous statements in dimension :math:`3` are false. Namely, there is no correlation between the topology of a coLegendrian (i.e., a surface) and the topology of its dividing set (i.e., a transverse link).

The proofs of the above two rigidity statements both rely on the simple fact that :math:`\op{tb}(\gamma) < 0` for any Legendrian :math:`\gamma \subset S^2 \subset (S^3, \xi_{\std})`. [#tb_ineq]_ Let's start with the simpler statement about the connectedness of :math:`\Gamma_Y`. Suppose, on the contrary, that :math:`R_+(Y)` has more than one connected component, and let :math:`K \subset R_+(Y)` be a component. By construction, :math:`\p K` is a surface transverse to the characteristic foliation :math:`Y_{\xi}`. The idea is to extend :math:`K` following :math:`Y_{\xi}` as much as possible, i.e., we'd like to consider the closure of all the trajectories of :math:`Y_{\xi}` flowing out of :math:`\p K`, denoted by :math:`\overline{K}`, and draw a contradiction. In practice, the trajectories are packed into (coLegendrian) handles, so we shall keep attaching (partial) handles to :math:`K` until we run into trouble.

.. sidebar:: Borderline between transverse and tangent boundaries

    .. figure:: static/coleg-5d/tb-nonnegative.svg
        :align: center
        :width: 400px

Since :math:`K` is not the only component, there must be at least one :math:`H_1^-` which connects :math:`K` to another (positive) component. Then within :math:`\overline{K}`, only half of the :math:`H_1^-`, cut by the unstable disk, are attached to :math:`K`. Let :math:`K_1` be the resulting domain. Then :math:`\p K_1` can be decomposed into the transverse part :math:`\p_{\tau} K_1` and the tangent part :math:`\p_t K_1`, which is nothing but the unstable disk of :math:`H_1^-`. Now :math:`\p_{\tau} K_1` can be viewed as a surface in a contact :math:`3`-manifold with a Legendrian boundary :math:`\gamma` with :math:`\op{tb} (\gamma) = 1`. The characteristic foliation on :math:`\p_{\tau} K_1` near :math:`\gamma` is shown on the left-hand-side in the picture on the right, where :math:`\gamma` is drawn in blue. Now if no :math:`2`-handles pass through :math:`H_1^-`, then :math:`\gamma` must lie on the boundary of a :math:`3`-handle, i.e., a :math:`2`-sphere in :math:`(S^3, \xi_{\std})`, and this is impossible. On the other hand, if there are :math:`2`-handles passing through :math:`H_1^-`, the characteristic foliation on :math:`\p_{\tau} K_1` near :math:`\gamma` may become more complicated. But in any case, a piece of :math:`\gamma` will be completed to a Legendrian with :math:`\op{tb} \geq 0` after all the (partial) :math:`2`-handle attachments. An example is shown on the right-hand-side of the picture, where a :math:`\op{tb} = 0` Legendrian (depicted also in blue) must lie on the boundary of a :math:`S^2 \subset (S^3, \xi_{std})` as before. This is also impossible and we run into a contradiction where the only way out is that the very first half-attached :math:`H_1^-` shouldn't exist in the first place. This is equivalent to the statement that :math:`\Gamma_Y` is connected.

tbc...

.. rubric:: Footnotes

.. [#tb_ineq] This inequality was first discovered by D. Bennequin to argue that there exists more than one contact structures on :math:`\Rbb^3`. One can prove it using (parametrized) :ref:`contact Morse theory <contact_morse_theory_rst>` by showing that the only possible Morse configuration on :math:`S^2 \subset (\Rbb^3, \xi_{\std})` is the standard one, modulo Weinstein homotopy.
