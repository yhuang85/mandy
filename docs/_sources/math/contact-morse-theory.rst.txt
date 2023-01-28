.. _contact_morse_theory_rst:

Morse theory in contact topology
================================

In this article, we explain how to cut any contact manifold into simple pieces as well as how to manipulate them from a Morse-theoretic viewpoint. The math of this article is based on many years of collaboration with K. Honda, but the interpretations are mostly mine. Moreover, the focus will be on ideas rather than technical details: one can hardly find any proofs here.

Morse Theory in Topology
------------------------

The most basic and general idea of Morse theory is to understand global properties of a topological space by cutting it into simple pieces and keep track of the assemble process. There are many variations and complications that might sneak behind the simply looking idea. For example, the topological space itself may be `singular <https://en.wikipedia.org/wiki/Stratified_Morse_theory>`__ or `infinite dimensional <https://en.wikipedia.org/wiki/Floer_homology>`__, and the "simple pieces" and "assemble process" may take very `different forms <https://en.wikipedia.org/wiki/Triangulation_(topology)>`__.

We will be concerned with the most classical Morse theory here, where spaces are smooth finite-dimensional manifolds :math:`M` and decompositions of :math:`M` are given by the so-called Morse functions :math:`f: M \to \Rbb`, which are nothing but *generic* smooth functions. More precisely, the only data we care about :math:`f` is a sequence of *singular values*

.. math::

	\dots < c_{k-1} < c_k < c_{k+1} < \dots

i.e., the values at critical points where :math:`df=0`, which is bounded on both sides if :math:`M` is compact.

Any value in between two adjacent singular values is a *regular value*. Let's take, for example, two regular values :math:`a_{k-1} \in (c_{k-1}, c_k)` and :math:`a_k \in (c_k, c_{k+1})`, jamming :math:`c_k` in between. Then we cut out of :math:`M` a piece bounded between two smooth hypersurfaces :math:`\Sigma_{k-1} \coloneqq f^{-1} (a_{k-1})` and :math:`\Sigma_k \coloneqq f^{-1} (a_k)`. This piece is simple enough in the sense that it's always topologically equivalent to a *handle attachment* regardless of what manifold we're talking about: it only depends on a local invariant called the Morse index of :math:`c_k`.

Morse theory aims at recognizing manifolds by their handle decompositions. Most importantly, although one same manifold may have many different handle decompositions, as many as there are generic functions, they are all connected to each other via a sequence of handle manipulations. It is therefore fair to say the following:

	Morse theory is all about manipulations of handles: creations, cancellations, and isotopies.


.. _first_blend:

First Blend of Morse Theory with Contact Structures
---------------------------------------------------

A symplectic detour
*******************

Let's actually start by talking about the attempt of blending Morse theory into symplectic topology initiated by Y. Eliashberg and M. Gromov in [EG91]_. Let :math:`(W, \omega)` be a symplectic manifold and :math:`f: W \to \Rbb` be a Morse function. The most obvious guess of a compatibility condition between :math:`f` and :math:`\omega`, at least from a Morse theoretic viewpoint, would be that the flow of the gradient vector field :math:`\nabla f \eqqcolon X`, with respect to some metric [#gradient]_ , preserves :math:`\omega`. This is not a very good guess though since :math:`X` would then necessarily be volume preserving, and in turn :math:`f` can have neither local minima nor maxima.

It's somewhat unfortunate that no reasonable compatibility condition exists between :math:`f` and :math:`\omega` which would simply allow :math:`f` to simultaneously possess local minima and maxima as a usual smooth function does. The compromise, which is by no means obvious, is to ask :math:`X` to exponentially expand :math:`\omega`, i.e., :math:`\Lcal_X \omega = \omega`. Quite obviously such thing cannot exist on a closed manifold since :math:`X` exponentially expands the volume. In fact, with a bit more work, one realizes that the Morse index of any critical point of :math:`f` cannot exceed :math:`\tfrac{1}{2} \dim W`. In some sense, we just compromised half of the topology.

Recall the two major players in the usual Morse theory are regular level sets and critical points (or equivalently, handles). Not surprisingly, they receive additional structures from the compatibility condition with a symplectic structure.

Regular level sets
++++++++++++++++++

For any regular value :math:`c`, the regular level set :math:`\Sigma \coloneqq f^{-1} (c) \subset (W, \omega)` is a contact manifold with contact form :math:`\alpha \coloneqq i_X \omega`.

Symplectic handles
++++++++++++++++++

Assuming :math:`\dim W = 2n`, an index-:math:`k` handle corresponds to an index-:math:`k` critical point always takes the following standard form

.. math::

	\begin{aligned}
		H_k &= [-1, 1]^k_{x_1, \cdots, x_k} \times [-1, 1]^{2n-k}_{x_{k+1}, \cdots, x_n, y_1, \cdots, y_n},

		\omega|_{H_k} &= \sum_{i=1}^{2n} dx_i \wedge dy_i,

		X|_{H_k} &= \sum_{i=1}^k ( -x_i \p_{x_i} + 2y_i \p_{y_i} ) + \frac{1}{2} \sum_{j=k+1}^{2n} ( x_j \p_{x_j} + y_j \p_{y_j} ).
	\end{aligned}

Such handles are known as *Weinstein handles* because A. Weinstein first wrote down these normal forms. A particularly important special case is when :math:`k = n` and we we call these handles *critical* because only these handles actually carry nontrivial symplectic information.


Back to contact
***************

The compatibility between a contact manifold :math:`(M, \xi)` and a Morse function :math:`f: M \to \Rbb`, by analogy, asks the flow of :math:`X \coloneqq \nabla f` to preserve :math:`\xi`.  This turns out to be a much more flexible condition than its symplectic counterpart. Indeed, Giroux argued in [Gir02]_ that every (closed) contact manifold admits a compatible Morse function.

The regular level sets and critical points of a contact Morse function look quite different from the symplectic case, mostly because the compatibility condition now reads :math:`\Lcal_X \alpha = g \alpha` where :math:`\xi = \ker \alpha` is a contact form and :math:`g` can be *any* function. In particular, the flow of :math:`X` doesn't have to expand the contact volume at all.

Regular level sets
++++++++++++++++++

A regular level set :math:`\Sigma \subset (M, \xi)` can be decomposed into three pieces

.. math::
	:label: hypersurface-decomposition

	\Sigma = R_+ \cup \Gamma \cup R_-,

where :math:`\Gamma = \{ \alpha (X) = 0 \}` and :math:`R_{\pm} = \{ \pm \alpha (X) > 0 \}`, respectively. It turns out that :math:`(\Gamma, \xi|_{\Gamma})` is itself a contact manifold of dimension :math:`\dim M - 2`, i.e., it's a codimension-:math:`2` contact submanifold of :math:`(M, \xi)`. Moreover :math:`( R_{\pm}, d\alpha|_{R_{\pm}} )` are symplectic manifolds. However, they don't necessarily carry a Morse structure as described above.

.. note::

	* We will rewrite the decomposition :eq:`hypersurface-decomposition` as :math:`\Sigma = R_+ \cup_{\Gamma} R_-` to highlight the viewpoint that :math:`\Sigma` can be obtained by gluing (closures of) :math:`R_{\pm}` along the common boundary :math:`\Gamma`.

	* It's shown by `Giroux <https://eudml.org/doc/140253>`__ that any hypersurface transverse to a (locally defined) contact vector field admits a decomposition as in :eq:`hypersurface-decomposition`. Such hypersurfaces were named *convex* by Eliashberg and Gromov in a `paper <https://www.ihes.fr/~gromov/wp-content/uploads/2018/08/976.pdf>`__ which covers both symplectic and contact cases. However, while convexity makes perfect sense in the symplectic world (e.g. it synchronizes well with convexities in complex and Riemannian geometry wherever these subjects overlap), it doesn't make any sense in the contact world. Indeed, they're more of a "flat" kind because the contact structure is invariant in the transverse direction. This is the main reason why we don't use the term "convex hypersurface" in this article. Another reason for not considering hypersurfaces like :eq:`hypersurface-decomposition` in general is that the domains :math:`R_{\pm}` are not necessarily Morse-theory friendly, i.e., they may be Liouville and not Weinstein. More about the later point will be elaborated in the :ref:`second blend <section_2nd_blend_morse_contact>`.

Contact handles
+++++++++++++++

Assuming :math:`\dim M = 2n+1` and :math:`k \leq n`, an index-:math:`k` contact handle always takes the following standard form

.. math::

	\begin{aligned}
		H_k &= [-1, 1]^k_{x_1, \cdots, x_k} \times [-1, 1]^{2n-k}_{x_{k+1}, \cdots, x_n, y_1, \cdots, y_n} \times [-1, 1]_z,

		\alpha|_{H_k} &= dz - \sum_{i=1}^n y_i dx_i,

		X|_{H_k} &= \sum_{i=1}^k ( -x_i \p_{x_i} + 2y_i \p_{y_i} ) + \frac{1}{2} \sum_{j=k+1}^{2n} ( x_j \p_{x_j} + y_j \p_{y_j} ) + z dz.
	\end{aligned}

Note that :math:`X|_{H_k}` exponentially expands the contact volume in the above model. For :math:`k \geq n+1`, one can simply reverse the signs of :math:`X|_{H_k}` in the above equation, in which case :math:`X|_{H_k}` exponentially contracts the contact volume. These handles will just be called contact handles since nobody was interested in registering them as trademarks.

Summary
*******

The fact that every (closed) contact manifold admits a compatible Morse function means that one can build any contact manifold from the standard-looking contact handles. However, it doesn't really give us much more grip on the contact manifold itself because such handle decompositions are by no means unique. Indeed, most of the power of Morse theory lies in the ability to connect different choices of Morse functions by homotopies. Such homotopies or more generally the flexibility of contact Morse functions are unfortunately not available from `Giroux's argument <https://arxiv.org/abs/math/0305129>`__ due to the global nature of the argument, which is more-or-less a replica of an argument of `S. Donaldson <https://projecteuclid.org/journals/journal-of-differential-geometry/volume-44/issue-4/Symplectic-submanifolds-and-almost-complex-geometry/10.4310/jdg/1214459407.full>`__ for symplectic manifolds which are far more rigid.

We will follow a completely different path to build a hopefully more useful contact Morse theory. The main tools will be hypersurfaces and characteristic foliations on them. These tools are native to contact topology and were extensively used by D. Bennequin, Eliashberg, Giroux among many others for various purposes in the early days of the subject.

.. _section_main_ideas:

The Main Ideas
--------------

The main ideas in the attempt to understanding contact structures via Morse theory were introduced in [HH18]_ and [HH19]_, which we briefly recall now.

Separation of contact structure and topology
********************************************

The first key principle in the development of contact Morse theory is to separate the contact topological problem from the purely topological problem. Specifically, given a contact manifold :math:`(M, \xi)`, we always start from just any Morse function :math:`f: M \to \Rbb`. Suppose :math:`\dim M = 2n+1`, then :math:`f` induces a decomposition

.. math::
	:label: heegaard-decomposition

	M = H_1 \cup (\Sigma \times I) \cup H_2, \quad I \coloneqq [0,1],

where :math:`H_1` is a neighborhood of the union of the stable manifolds of all critical points of :math:`f` of index at most :math:`n`, and similarly :math:`H_2` is a neighborhood of the union of the unstable manifolds of all critical points of :math:`f` of index at least :math:`n+1`, and finally :math:`\Sigma` may be identified with either :math:`\p H_1` or :math:`\p H_2`. Such a decomposition is nothing but a higher-dimensional analogue of the `Heegaard decomposition <https://en.wikipedia.org/wiki/Heegaard_splitting>`__ for 3-manifolds.

The reason for decomposing :math:`M` as in :eq:`heegaard-decomposition` is because the contact structures :math:`\xi|_{H_1}, \xi|_{H_2}` can be made standard by Gromov's *h*-principle on isotropic/Legendrian approximations. Namely, the stable manifold of all critical points of index at most :math:`n` can be :math:`C^0`-approximated by isotropic submanifolds. It follows that :math:`H_1` becomes a neighborhood of a CW-complex built out of isotropic cells, which in turns carries a standard contact structure. The same holds for :math:`H_2` by reversing the direction of :math:`\nabla f`. Finally, the complement of :math:`H_1 \cup H_2` in :math:`M` is a product :math:`\Sigma \times I`, which, in some sense, carries all the information about :math:`\xi`.

To summarize, the decomposition :eq:`heegaard-decomposition` serves the purpose of separating contact topology from pure topology as follows:

* The handlebodies :math:`H_1, H_2` knows all about the topology of :math:`M` but nothing about :math:`\xi`.

* The middle layer :math:`\Sigma \times I` knows all about :math:`\xi` but nothing about :math:`M`.

.. note::

	The handlebodies :math:`H_1, H_2` in :eq:`heegaard-decomposition` are by no means unique, although the (restricted) contact structures are uniquely determined by their topological type.

Morsify the characteristic foliation
************************************

Given any hypersurface :math:`\Sigma \in (M, \xi)`, the *characteristic foliation* :math:`\Sigma_{\xi}` is a line field defined by

.. math::

	\Sigma_{\xi} \coloneqq \ker (d\alpha|_{T\Sigma \cap \xi}),

where :math:`\xi = \ker\alpha` is a contact form. Moreover, when both :math:`\Sigma` and :math:`\xi` are oriented, which will always be the case here, so is :math:`\Sigma_{\xi}` and it becomes a vector field (without any significance on the magnitude). Characteristic foliations play a crucial role in this story because they, to a great extent which will become clear later, uniquely determines the contact germ on :math:`\Sigma`.

Now the job of characterizing a contact structure on :math:`\Sigma \times I` boils down to characterizing the evolution of the contact germs on :math:`\Sigma \times t` for :math:`t \in I`, which, in turn, boils down to characterizing the evolution of a :math:`1`-parameter family of vector fields :math:`(\Sigma \times t)_{\xi}, t \in I`.

In the case of :math:`\dim M = 3`, the above job was successfully done by Giroux in [Gir99]_, where he applied the results of such analysis to classify contact structures on a number of :math:`3`-manifolds. However, Giroux's argument from these two papers are rather specific to dimension :math:`3` and are of little use in higher dimensions [#false_believe]_. Specifically, the study of characteristic foliations on a :math:`2`-dimensional surface falls into a much bigger subject of studying dynamics of generic vector fields on surfaces. Rather mature and comprehensive theories on the later subject, such as the `Poincaré-Bendixson theorem <https://en.wikipedia.org/wiki/Poincar%C3%A9%E2%80%93Bendixson_theorem>`__, was developed long before contact topology was even recognized as an independent subject. However, it's indeed hopeless to track down every single trajectory of a generic vector field in dimensions greater than two due to the ubiquity of chaotic behavior.

The challenge is, therefore, to ensure the controllability of :math:`\Sigma_{\xi}` (e.g., as the gradient vector field of a Morse function) on sufficiently generic hypersurfaces :math:`\Sigma`. The basic idea is to wiggle :math:`\Sigma` almost everywhere to create attractors, built out of Morse critical points, which destroy any potential global dynamics of :math:`\Sigma_{\xi}`. The actual implementation of this idea is nearly perfect in dimension :math:`3` but much less so in higher dimensions. The details can be found in [HH19]_.

Summary
*******

Every (closed) contact manifold can be decomposed into three pieces: two standard contact handlebodies and a product :math:`\Sigma \times I`. The contact structure :math:`\xi|_{\Sigma \times I}` can be understood via the :math:`1`-parameter family of characteristic foliations :math:`\Sigma_t|_{\xi} \coloneqq (\Sigma \times t)_{\xi}` for :math:`t \in I`. The characteristic foliations :math:`\Sigma_t|_{\xi}, t \in I`, can be made Morse by a :math:`C^0`-small perturbation. Thus the problem is finally reduced to understanding a :math:`1`-parameter family of Morse functions on :math:`\Sigma`. More details about implementing these ideas will be explained in the next section.

.. _section_2nd_blend_morse_contact:

Second Blend of Morse Theory with Contact Structures
----------------------------------------------------

Recall in the first blend of Morse theory with contact structures, the result is a decomposition of :math:`(M, \xi)` into a bunch of contact handles. This approach appears to be somewhat useless since there is no way (that I know of) to connect two contact Morse functions through a family of contact Morse functions.

Instead, we'll use the ideas outlined above to build a contact Morse theory which works in families. To facilitate the exposition, let's use the following convention to indicate the dimension of the family of Morse functions under consideration. We say a Morse theory (of whatever flavor) is established at

* :math:`\pi_0`-level if Morse functions exist generically,
* :math:`\pi_1`-level if any two Morse functions are homotopic through Morse functions,
* :math:`\pi_2`-level if a circle-family of Morse functions can be realized as the boundary of a disk-family of Morse functions,
* and so on for :math:`\pi_k`-levels for :math:`k > 2`.

.. note::

	Critical points, among others, in families of Morse functions degenerate according to the standard `transversality theory <https://en.wikipedia.org/wiki/Transversality_theorem>`__ on jet bundles. For example, critical points are nondegenerate at :math:`\pi_0`-level but may degenerate to birth-death type singularities at :math:`\pi_1`-level and swallowtails at :math:`\pi_2`-level and so on.

For example, the usual Morse theory is fully-established in the category of smooth functions and provides deep insights into the structure of smooth manifolds via `Cerf theory <https://en.wikipedia.org/wiki/Cerf_theory>`__, `h-cobordism theorem <https://en.wikipedia.org/wiki/H-cobordism>`__ and so on. In the contact category, we need to at least impose one additional compatibility condition between functions and contact structures: the gradient vector field must preserve the contact structure. However, as we'll see, this condition alone is not enough to build a useful (family) contact Morse theory.

Topological skeleta
*******************

Recall that although the existence of contact Morse functions, in abundance as a matter of fact, was established in [Gir02]_, nearly no flexibility is available for these rather abstract functions, which makes it hardly useful in practice. On the other hand, one cannot expect genericity to hold in the sense of usual transversality theory as in the smooth case because contact structures are by no means generic in that sense.

As a matter of fact, it makes little sense to even look for (generic) homotopies between contact Morse functions because it violates the first principle of separation between topological and contact topological considerations. Instead, let's emphasize once again that the decomposition :eq:`heegaard-decomposition` is always the first step when decomposing a contact manifold :math:`(M, \xi)`. Recall that the contact handlebodies :math:`H_1, H_2 \subset M` are uniquely determined by the corresponding isotropic skeleta, which also capture the topology of :math:`M`. For this reason, we introduce the following terminology:

	Up to a negligible ambiguity, either :math:`H_1, H_2` or their skeleta are called *topological skeleta* of :math:`(M, \xi)`.

Of course, one contact manifold may have many different topological skeleta, and it's far from obvious how two choices are related to each other in a Morse theoretic way. However, such difficulty doesn't bother us, at least for now, since we're not really interested in the topology of :math:`M`. Indeed, it'd already be a great success of contact Morse theory if one could get some insights into contact structures on :math:`S^{2n+1}, n \geq 2`.

.. _section_family_morse:

A family Morse theory on hypersurfaces
**************************************

Away from the topological skeleta, the contact manifold reduces to a product :math:`\Sigma \times I` as in :eq:`heegaard-decomposition`. As explained in the :ref:`main ideas <section_main_ideas>`, up to a :math:`C^0`-small perturbation, the characteristic foliations :math:`\Sigma_t|_{\xi}, t \in I` can be realized as the gradient of a :math:`1`-parameter family of Morse functions on :math:`\Sigma`. It is this Morse theory which can be made "generic" and work in families. In what follows, we'll spell out the details of this Morse theory on hypersurfaces at :math:`\pi_0, \pi_1`, and :math:`\pi_2`-levels. As a convention, all explicitly mentioned (Morse) critical points are assumed to be nondegenerate unless otherwise specified.

:math:`\pi_0`-level
+++++++++++++++++++

The :math:`\pi_0`-level Morse theory means that for any :math:`t_0 \in I`, the hypersurface :math:`\Sigma = \Sigma_{t_0}` can be :math:`C^0`-perturbed such that :math:`\Sigma_{\xi}` is Morse. Let :math:`p \in \Sigma` be a critical point. Then we say :math:`p` is *positive* if :math:`T_p \Sigma = \xi_p` as oriented vector spaces and *negative* if :math:`T_p \Sigma = -\xi_p`. It turns out that the stable manifolds of the positive critical points build up a Weinstein manifold :math:`R_+ \subset \Sigma`, i.e., a symplectic manifold built out of (finitely many) Weinstein handles explained in the :ref:`first blend <first_blend>`. Likewise, the unstable manifolds of the negative critical points build up another Weinstein manifold :math:`R_- \subset \Sigma`. Denoting the remaining borderline between :math:`R_+` and :math:`R_-` by :math:`\Gamma`, we arrive at the familiar :math:`\Sigma = R_+ \cup_{\Gamma} R_-` which appeared as the structure of a regular level set in :eq:`hypersurface-decomposition`.

	We say a hypersurface :math:`\Sigma` is *Morse* if :math:`\Sigma_{\xi}` is Morse. Moreover, genericity is always appropriately understood according to the :math:`\pi_k`-level of the Morse theory under discussion.

.. note::

	Morse hypersurfaces are not generic. They are only :math:`C^0`-dense among all hypersurfaces, which is enough for all we care. It's important to note that contact Morse theory lives on hypersurfaces rather than the contact manifold itself.

:math:`\pi_1`-level
+++++++++++++++++++

Suppose :math:`\Sigma_0, \Sigma_1` are Morse, where :math:`\Sigma_t \coloneqq \Sigma \times t, t \in I`. This is indeed the case when they are boundaries of standard neighborhoods of the isotropic skeleta :math:`H_0, H_1`. Then the :math:`\pi_1`-level Morse theory means that, up to a :math:`C^0`-small perturbation, the :math:`1`-parameter family :math:`\Sigma_t|_{\xi}` can be realized as the gradient of a :math:`1`-parameter family of Morse functions. It turns out that for most of the time :math:`t \in I`, the contact germ on :math:`\Sigma_t` doesn't change, up to isotopy.

	We say a Morse hypersurface is *invariant* if the contact germ is invariant in the transverse direction. This is equivalent to, as it turns out, the nonexistence of flow lines from negative critical points to positive critical points.

Due to genericity and the index constraint on Weinstein handles, :math:`\Sigma_t` may fail to be invariant only when there is a (unique) trajectory of :math:`\Sigma_t|_{\xi}` from a negative index-:math:`n` critical point :math:`p_n^-` to a positive index-:math:`n` critical point :math:`p_n^+`, assuming :math:`\dim \Sigma = 2n`. Moreover, such failure may happen for only finitely many :math:`t \in I`, which we call the :math:`\pi_1`-*critical moments*.

	Depending on the context, a :math:`\pi_1`-*switch* at a :math:`\pi_1`-critical moment :math:`t_0 \in I` can refer to any one of the following:

	* The (transversely cut out) trajectory from :math:`p_n^-` to :math:`p_n^+`.
	* The hypersurface :math:`\Sigma_{t_0}`.
	* The contact structure on :math:`\Sigma \times [t_0 - \epsilon, t_0 + \epsilon]` for :math:`\epsilon > 0` sufficiently small.

Topological speaking, the difference between :math:`\Sigma_{t_0 - \epsilon}` and :math:`\Sigma_{t_0 + \epsilon}` is a handle slide of a negative :math:`n`-handle over a positive :math:`n`-handle. However, not every topological handle slide of this kind can be realized as a :math:`\pi_1`-switch, even after assuming all isotopies involved in the handle slide are contact isotopies. Namely, suppose :math:`Y \subset \Sigma_{t_0}` is a regular level set between :math:`p_n^-` and :math:`p_n^+` such that the unstable manifold of :math:`p_n^-` intersects :math:`Y` along a Legendrian sphere :math:`\Lambda_-` and the stable manifold of :math:`p_n^+` intersects :math:`Y` along :math:`\Lambda_+`. Here we recall :math:`Y` is naturally a contact submanifold. [#contact_submfd]_ Then :math:`\Lambda_{\pm}` intersect :math:`\xi|_Y`-transversely at exactly one point :math:`q` (on the :math:`\pi_1`-switch), i.e.,

.. math::
	:label: xi-transverse

	T_q \Lambda_+ \oplus T_q \Lambda_- = (\xi|_Y)_q.

Extending the definitions of :math:`Y` and :math:`\Lambda_{\pm}` to all :math:`t` close to :math:`t_0`, we require that :math:`\Lambda_+` is slightly "below" :math:`\Lambda_-`, measured against the positive co-orientation of :math:`\xi|_Y`, near :math:`q` for :math:`t < t_0` and "above" for :math:`t > t_0`.

	In plain words, the handle slide corresponding to a :math:`\pi_1`-switch isotopes :math:`\Lambda_+` up across :math:`\Lambda_-` as :math:`t` passes over :math:`t_0`.

.. note::

	Historically speaking, a :math:`\pi_1`-switch is trivially a special case of "bifurcations" considered by Giroux [Gir99]_ in his dynamical convex surface theory, and less trivially a special case of the so-called "bypass attachments" considered by Honda [Hon99]_ in his combinatorial convex surface theory, both in dimension :math:`3`. The later was generalized to all dimensions in [HH18]_. In particular, the decomposition :eq:`heegaard-decomposition` indeed gives rise to a contact Morse function. However, none of these developments are relevant here and we don't even care about general contact Morse functions per se.

Besides :math:`\pi_1`-switches, there are many other :math:`\pi_1`-level Morse theoretic degenerations, such as creation and elimination of critical points, that may happen in the family :math:`\Sigma_t|_{\xi}, t \in I`. However, these phenomena may happen either within :math:`R_+` or :math:`R_-`, and they belong to the subject of Weinstein homotopies, whose general understanding is completely out of reach by the current technology.

To summarize, the :math:`\pi_1`-level contact Morse theory asserts that, modulo Weinstein homotopies, any contact structure on :math:`\Sigma \times I` can be realized as a finite sequence of :math:`\pi_1`-switches.

:math:`\pi_2`-level
+++++++++++++++++++

The :math:`\pi_2`-level contact Morse theory aims at connecting two realizations of the same :math:`(\Sigma \times I, \xi)` as :math:`1`-parameter families of Morse functions on :math:`\Sigma`. It's therefore inappropriate to ignore the :math:`C^0`-perturbation part and pretend that the realizing hypersurface foliation is just :math:`\Sigma_t = \Sigma \times t, t \in I`. For the sake of distinction, let :math:`\Sigma^0_t, t \in I`, and :math:`\Sigma^1_t, t \in I`, be two different foliations realizing :math:`\pi_1`-level Morse theories as explained above. Namely, modulo Weinstein homotopies, the :math:`1`-parameter families :math:`\Sigma^0_t|_{\xi}` and :math:`\Sigma^1_t|_{\xi}` give rise to two compositions of :math:`\pi_1`-switches. Therefore, the goal is, roughly speaking, to connect different compositions of :math:`\pi_1`-switches which define the same contact structure, in a Morse theoretic way.

By analogy with the :math:`\pi_1`-switch, here is a complete list of :math:`\pi_2`-switches which at some point breaks the invariance of the contact germ. First of all, we need to work on the :math:`2`-dimensional parameter space :math:`(s, t) \in I^2`, where :math:`t` shall always parametrize the foliations and :math:`s` parametrizes the homotopies. At a :math:`\pi_2`-critical moment :math:`(s_0, t_0) \in I^2`, one of the following scenarios may happen: [#pi2labels]_

* (:math:`\pi_2^a`-switch) There exist a negative birth-death-type index-:math:`(n+1)` critical point :math:`p_{n+1, n}^-` and a positive index-:math:`n` critical point :math:`p_n^+`, such that there is a unique transversely cut out trajectory from :math:`p_{n+1, n}^-` to :math:`p_n^+`. Here the notation :math:`p_{n+1, n}` for a birth-death-type critical point indicates the dimension of the stable manifold, which is :math:`n+1`, and the unstable manifold, which is :math:`n`. Flipping the orientation, one also has the same type of switch at a trajectory from a negative :math:`p_n^-` to a positive :math:`p_{n, n-1}^+`.

* (:math:`\pi_2^b`-switch) There exist two index-:math:`n` negative critical points :math:`p_n^-, q_n^-` and two positive :math:`p_n^+, q_n^+`, such that there are exactly two transversely cut out trajectories: one from :math:`p_n^-` to :math:`p_n^+` and the other from :math:`q_n^-` to :math:`q_n^+`.

.. _bullet_pi2c_switch:

* (:math:`\pi_2^c`-switch) There exist a negative index-:math:`(n+1)` critical point :math:`p_{n+1}^-` and a positive index-:math:`n` critical point :math:`p_n^+` and a unique trajectory from :math:`p_{n+1}^-` to :math:`p_n^+` which is transversely cut out with respect to a :math:`2`-dimensional family :math:`\Sigma_t^s|_{\xi}` for :math:`(s, t)` close to :math:`(s_0, t_0)`.

.. _bullet_pi2d_switch:

* (:math:`\pi_2^d`-switch) There exist a negative index-:math:`n` critical point :math:`p_n^-` and a positive :math:`p_n^+`, such that there exists a trajectory from :math:`p_n^-` to :math:`p_n^+` which is not transversely cut out, but rather has a first-order tangency. Namely, let :math:`Y \subset \Sigma_{t_0}^{s_0}` be a regular level set between :math:`p_n^-` and :math:`p_n^+`, and :math:`\Lambda_{\pm} \subset Y` be Legendrian spheres just as in the above discussion at the :math:`\pi_1`-level. Then the unique intersection :math:`q = \Lambda_+ \cap \Lambda_-` satisfies the following

  .. math::
	  :label: xi-degenerate

	  \dim(T_q \Lambda_+ \cap T_q \Lambda_-) = 1.

  This should be compared with the :math:`\xi|_Y`-transversality condition :eq:`xi-transverse`.


Summary
*******

We start with the definition of topological skeleta, which serve the purpose of separating topology from contact structures. Then we proceed with a description of the sought-after contact Morse theory on hypersurfaces from :math:`\pi_0` to :math:`\pi_2`-level. The :math:`\pi_0`-level is the foundation for everything that follows and technically speaking, it involves all the (good and bad) techniques established in [HH19]_. The :math:`\pi_1`-level reduces the study of contact structures to the study of finite sequences of :math:`\pi_1`-switches. Finally, the :math:`\pi_2`-level provides a complete list of moves one needs to compare two different sequences of :math:`\pi_1`-switches. In principle, one could continue to build :math:`\pi_k`-level contact Morse theory for :math:`k \geq 3`. We choose not to do that for two reasons: first, as far as the classification of contact structures is concerned, the :math:`\pi_2`-level Morse theory suffices, and second, there is no significant technical advancement already from :math:`\pi_0`-level up.


Examples
--------

So far the theory has been dry and obscure. We need examples to make it sensible but as for any other theories, there is a high risk of breaking it by testing against the reality. So let's do it.

.. _section_r_pm_picture_of_pi_1_switches:

:math:`R_{\pm}`-picture of :math:`\pi_1`-switches
*************************************************

The Morse picture of :math:`\pi_1`-switches is conceptually clear but can be difficult to use in practice. So let's introduce a slightly different approach, called the :math:`R_{\pm}`-picture, which focuses less on the (Morse) gradient vector field and more on the critical points, making it easier to manipulate, especially when combined with front projections. In a nutshell, the :math:`R_{\pm}`-picture describes the changes in :math:`R_{\pm} (\Sigma_t)`, as well as how they are glued together along :math:`\Gamma(\Sigma_t)`, as :math:`\Sigma_t|_{\xi}` goes through a :math:`\pi_1`-switch.

.. _figure_r_pm_picture_of_pi_1_switch:

.. sidebar:: :math:`R_{\pm}`-picture of a :math:`\pi_1`-switch

	.. figure:: static/contact-morse-theory/pi1-handles.svg
		:align: center
		:width: 400px

The picture on the right-hand-side illustrates a completely general :math:`\pi_1`-switch decomposed into three steps, i.e., the three dashed arrows, which we now explain. Unlike the previous discussions in the :ref:`family Morse theory <section_family_morse>`, here we need to keep track of several level sets (in :math:`\Sigma`) at once and both stable and unstable manifolds of the critical points. So the notations will unfortunately become a bit more cluttered. Note that the gradient vector field (i.e., the characteristic foliation) always flows upwards (indicating that I'm not a physicist).

The upper-left corner represents a part of :math:`\Sigma` relevant to the :math:`\pi_1`-switch. Namely, there are two index-:math:`n` critical points :math:`p_n^+` and :math:`p_n^-`, and the corresponding Legendrian spheres :math:`\Lambda_+^u, \Lambda_-^s \subset \Gamma`. Here the superscripts :math:`u` and :math:`s` denote unstable and stable, respectively. Moreover, there is a small ball in :math:`\Gamma` which intersects both :math:`\Lambda_+^u` and :math:`\Lambda_-^s` in a disk such that the :math:`\Lambda_+^u`-disk is slightly below the :math:`\Lambda_-^s`-disk, where "below" is measured against the positive co-orientation of :math:`\xi|_{\Gamma}`. This small ball is magnified in the picture, and the "below"-ness is shown as an under-crossing when the :math:`\Lambda`'s appear to be :math:`1`-dimensional.

The passage to the upper-right corner is nothing but swapping the critical values of :math:`p_n^+` and :math:`p_n^-`. Note that the small ball from above carries over to the new intermediate level set :math:`Y`, inside of which the :math:`\Lambda_+^s`-disk is slightly below the :math:`\Lambda_-^u`-disk.

The passage from the upper-right to the lower-right corner is where the :math:`\pi_1`-switch really takes place. Namely, we (contact) isotop :math:`\Lambda_+^s` up across :math:`\Lambda_-^u` within the small ball such that at exactly one moment, they :math:`\xi|_Y`-transversely intersect in a point. To keep things somewhat symmetric, we denote the resulting Legendrian spheres :math:`\Lambda_+^{s, \uparrow}` and :math:`\Lambda_-^{u, \downarrow}` as if :math:`\Lambda_-^u` is simultaneously lowered while :math:`\Lambda_+^s` is raised.

Finally, the passage from the lower-right to the lower-left corner swaps :math:`p_n^+` and :math:`p_n^-` back and leave in between a new level set :math:`\Gamma'`. We can describe the new decomposition :math:`\Sigma = R'_+ \cup_{\Gamma'} R'_-` in terms of the old one as follows.

	As a Weinstein manifold, :math:`R'_+` is obtained from :math:`R_+` by removing the handle corresponding to :math:`p_n^+` and attach a handle along :math:`(\Lambda_+^u \uplus \Lambda_-^s)^{\uparrow}`. Similarly :math:`R'_-` is obtained from :math:`R_-` by removing the handle corresponding to :math:`p_n^-` and attach a handle along :math:`(\Lambda_+^u \uplus \Lambda_-^s)^{\downarrow}`. Here :math:`\uplus`, which joins two Legendrians spheres into one, is an artifact of Legendrian handle slides and can be found in p. 17 of [HH18]_. Finally since :math:`R'_{\pm}` share the same boundary :math:`\Gamma'`, it admits two equivalent Legendrian surgery descriptions, and an explicit equivalence in terms of a contact isotopy.

.. _note_y_picture:

.. note::

	The :math:`R_{\pm}`-picture of a :math:`\pi_1`-switch put some emphasis on the evolution of the decomposition :eq:`hypersurface-decomposition` assuming :math:`\Sigma` is invariant. Such emphasis is not always necessary given the local nature of :math:`\pi_1`-switches. In this case we may simply remember the second dashed arrow in the :math:`R_{\pm}`-:ref:`picture <figure_r_pm_picture_of_pi_1_switch>` above, and call it the :math:`Y`-picture since it records what happens in the level set :math:`Y`.

Free :math:`\pi_1`-switches
***************************

Recall :math:`\pi_1`-switches are exactly where the contact germs on a :math:`1`-parameter family of Morse hypersurfaces change. They don't come for free in general for otherwise contact topology would be no different than just topology. It's generally difficult to verify the existence of a particular :math:`\pi_1`-switch inside a given contact manifold. However, there exists a class of :math:`\pi_1`-switches which can always be found at the vicinity of any (invariant) Morse hypersurface. These :math:`\pi_1`-switches are called *free*, and are the subject of discussion in this section.

The creation of free :math:`\pi_1`-switches is very much a procedure of creating something out of nothing. Not surprisingly, it involves certain :math:`\pi_2`-switches. Since :math:`\pi_2`-switches are directionless, all creations can be reversed to eliminations, which we'll omit.

.. _section_trivial_pi1_switches:

Trivial :math:`\pi_1`-switches
++++++++++++++++++++++++++++++

.. sidebar:: Morse picture of a :math:`\pi_2^a`-switch

	.. figure:: static/contact-morse-theory/pi2a-morse.svg
		:align: center
		:width: 400px

Let's start with the simplest scenario of a trivial :math:`\pi_1`-switch, which can be created by a :math:`\pi_2^a`-switch as shown in the right-hand-side picture.

.. note::

	All pictures will be drawn in dimension :math:`2`, but are supposed to illustrate the general situation in any dimension. For example, saddles usually represent index-:math:`n` critical points, assuming :math:`\dim \Sigma = 2n`.

Specifically, the square in the middle represents the parameter space :math:`I^2_{s,t}` where :math:`s` is horizontal and :math:`t` is vertical. The red dot at the center of :math:`I^2` and the corresponding Morse vector field represents the critical moment when the :math:`\pi_2^a`-switch takes place. The left-side of :math:`I^2` represents a Morse homotopy :math:`\Sigma_t^0|_{\xi}` which contains no critical moments, i.e., there are no :math:`\pi_1`-switches. However, the right-side of :math:`I^2` represents a Morse homotopy :math:`\Sigma_t^1|_{\xi}` which contains exactly one :math:`\pi_1`-switch [#pi2a_morse_sign]_. Scanning from left to right, one could say that a :math:`\pi_1`-switch is born via a :math:`\pi_2^a`-switch.

	The so-created :math:`\pi_1`-switch is said to be *trivial* since the corresponding contact structure on :math:`\Sigma \times I` is isotopic, relative to the boundaries, to the :math:`I`-invariant one (modulo Weinstein homotopies of :math:`R_{\pm}` as usual).

.. sidebar:: :math:`R_{\pm}`-picture of a trivial :math:`\pi_1`-switch

	.. figure:: static/contact-morse-theory/trivial-p1-handles.svg
		:align: center
		:width: 400px

Let's turn the Morse picture of the trivial :math:`\pi_1`-switch into the :math:`R_{\pm}`-picture as shown on the right. Specifically, the top figure illustrates the relative position between :math:`\Lambda_+` and :math:`\Lambda_-`, which is the standard Legendrian unknot, corresponding to critical points :math:`p_n^{\pm}`, respectively, in :math:`\Gamma`. In contrast to the :ref:`general picture <section_r_pm_picture_of_pi_1_switches>`, we drop the superscripts :math:`u, s` from the :math:`\Lambda`'s here because it's obvious from the context. Moreover, it's arranged so that :math:`\Lambda_+` and :math:`\Lambda_-` intersect :math:`\xi|_{\Gamma}`-transversely at a point, instead of :math:`\Lambda_+` being slightly below :math:`\Lambda_-`. This serves the sole purpose of attracting our attention to around the intersection point, and one can always go back to the other picture by pushing :math:`\Lambda_+` down (or :math:`\Lambda_-` up) slightly.

The two figures at the bottom represent the new :math:`\Gamma'` after the trivial :math:`\pi_1`-switch from the perspectives of :math:`R'_+` and :math:`R'_-`, respectively. As a sanity check, one can easily see that :math:`\Gamma'` is indeed isomorphic to the original :math:`\Gamma`. Here the :math:`(\pm 1)` beside the Legendrians are coefficients of Legendrian surgeries, and correspond to removing and adding a (index-:math:`n`) critical point, respectively.

Trivial :math:`\pi_1`-switches, as its name suggests, are quite boring. But when combined with :math:`\pi_2^b`-switches, they can produce many nontrivial :math:`\pi_1`-switches. This is our next step.

.. _section_derived_pi1_switches:

Derived :math:`\pi_1`-switches from the trivial one
+++++++++++++++++++++++++++++++++++++++++++++++++++

.. sidebar:: Morse picture of a :math:`\pi_2^b`-switch

	.. figure:: static/contact-morse-theory/pi2b-morse.svg
		:align: center
		:width: 400px

The picture on the right illustrates a general :math:`\pi_2^b`-switch, where each two adjacent ovals represent two disjoint regions on :math:`\Sigma`. The square in the middle is, as before, the parameter space :math:`I^2_{s, t}` and the red dot represents the critical moment when there exist simultaneously two flow lines from index-:math:`n` critical points :math:`p_n^-, q_n^-` to :math:`p_n^+, q_n^+`, respectively. The passage from the left side :math:`\Sigma^0_t|_{\xi}`, to the right side :math:`\Sigma^1_t|_{\xi}, t \in I`, changes the order of which two :math:`\pi_1`-switches occur. It is therefore also known as the *far commutativity* of two :math:`\pi_1`-switches, which are in a sense disjoint from each other.

The plan for generating new :math:`\pi_1`-switches from the trivial one is straightforward. Namely, we'll construct two disjoint :math:`\pi_1`-switches :math:`\Delta_1, \Delta_2` such that :math:`\Delta_1` is trivial and :math:`\Delta_2` becomes trivial after :math:`\Delta_1`. However, after swapping the order, neither :math:`\Delta_2` nor :math:`\Delta_1` is trivial anymore.

.. _figure_simple_pi_1_switch_from_pi_2b_switch:

.. sidebar:: Free :math:`\pi_1`-switch from :math:`\pi_2^b`-switch

	.. figure:: static/contact-morse-theory/trivial-rotation.svg
		:align: center
		:width: 400px

On the right is a not-particularly-good-looking illustration of the above plan. Namely, in the upper-left corner, we draw the relevant Legendrian spheres :math:`\Lambda^1_{\pm}` and :math:`\Lambda^2_{\pm}` corresponding to the two trivial :math:`\pi_1`-switches :math:`\Delta_1` (black) and :math:`\Delta_2` (blue), respectively. Although :math:`\Delta_1` is obviously trivial, it's not immediately clear that :math:`\Delta_2` is also trivial after :math:`\Delta_1`. The bottom figure shows, from the perspective of :math:`R^1_+`, that it's indeed the case where ":math:`\cong`" represents a Legendrian isotopy (via a handle slide).

By swapping the two :math:`\pi_1`-switches (and forget about :math:`\Delta_1`), we get on the upper-right corner the derived :math:`\pi_1`-switch :math:`\Delta_2` which is quite general since there is no additional restrictions on :math:`\Lambda^2_-` outside of the local picture except that it must belong to the stable manifold of a negative index-:math:`n` critical point. Note that this requirement doesn't contradict our setup in the upper-left corner since :math:`\Lambda^1_-` and :math:`\Lambda^2_-` are not linked as Legendrians.

.. note::

	1. Any free :math:`\pi_1`-switch derived from a trivial one via a :math:`\pi_2^b`-switch admits a right inverse in the sense that a :math:`\pi_1`-switch can be post-composed to produce the :math:`I`-invariant contact structure.

	2. The :math:`\pi_2^b`-switch itself assumes nothing about the involved :math:`\pi_1`-switches. Hence it's possible to derive new :math:`\pi_1`-switches from the ones that are for some reason known to exist, regardless of whether the later are trivial or not. In particular, there exist an abundance of :math:`\pi_1`-switches that are iteratively derived from the trivial one. They may all be called "free" since they can always be found within an invariant neighborhood of a hypersurface.

We shall later use this type of free :math:`\pi_1`-switches to detect flexibilities of contact structures near a hypersurface.

Free :math:`\pi_1`-switches from nowhere
++++++++++++++++++++++++++++++++++++++++

As we have learned from the discussions above, free :math:`\pi_1`-switches can be created out of the trivial one using the so-called far-commutativity. What we shall do in this section is the opposite. Namely, we shall construct free :math:`\pi_1`-switches without any reference to :math:`\pi_2`-switches, and from this we shall retrospectively argue that the trivial :math:`\pi_1`-switch exists.

Since we plan not to use any :math:`\pi_2`-switches, we need to go back to the root of contact Morse theory to see why some :math:`\pi_1`-switches always exist. Recall the :math:`Y`-:ref:`picture <note_y_picture>` of :math:`\pi_1`-switch where :math:`Y \subset \Sigma` is a level set separating the two critical points :math:`p_n^{\pm}`. Moreover, the corresponding Legendrian spheres :math:`\Lambda_+^s` from :math:`p_n^+` and :math:`\Lambda_-^u` from :math:`p_n^-` in :math:`Y` are relatively positioned as shown in this :ref:`picture <figure_r_pm_picture_of_pi_1_switch>` such that :math:`\Lambda_+^s` is slightly "below" :math:`\Lambda_-^u`. Then the :math:`\pi_1`-switch amounts to isotoping :math:`\Lambda_+^s` up across :math:`\Lambda_-^u`.

If :math:`Y` were to have an arbitrarily large (standard) neighborhood, then one can always perturb :math:`\Sigma` near -- but away from -- :math:`Y` to change the characteristic foliation :math:`\Sigma_{\xi}` so that it effectively raises :math:`\Lambda_+^s` up across :math:`\Lambda_-^u` [#bump_switch]_. However, in reality, one can only assume a finite neighborhood size of :math:`Y` depending on :math:`\Sigma` or rather the contact structure at the vicinity of it. Let :math:`\phi_t: Y \to Y`, :math:`t \in [0,1]`, be the contact isotopy induced by the perturbation such that  :math:`\phi_0 = \op{id}` and :math:`\phi_1` raises :math:`\Lambda_+^s` up across :math:`\Lambda_-^u`. Then we can define an "action" of such an isotopy by

.. math::

	A(\phi) \coloneqq \max_{x \in Y} \int_0^1 \alpha(\dot{\phi}_t(x)) dt

where :math:`\alpha` denotes the contact form on :math:`Y`. It turns out that the possible magnitude of :math:`A(\phi)` is constrained by the neighborhood size of :math:`Y`. This is to say that an isotopy :math:`\phi_t` can be arranged by a small perturbation of :math:`\Sigma` near :math:`Y` if and only if the action :math:`A(\phi)` is sufficiently small.

In dimension :math:`3`, the above discussion leads to nothing because :math:`Y` in this case is a transverse knot (or link), and hence there is a positive lower bound on :math:`A(\phi)` for any :math:`\phi_t` that lifts :math:`\Lambda_+^s` (in this case :math:`0`-dimensional) up across :math:`\Lambda_-^u`. So we get no free :math:`\pi_1`-switch without knowing *a priori* that the neighborhood size :math:`Y` is sufficiently large. This is not the case in dimensions at least :math:`5`.

.. sidebar:: Free :math:`\pi_1`-switch from a "small" isotopy

	.. figure:: static/contact-morse-theory/finger-move.svg
		:align: center
		:width: 400px

In the picture to the right, we show an example of an isotopy with small action. Namely, we draw an isotropic arc (which is Legendrian if :math:`\dim Y = 3`) from :math:`\Lambda_-^u` to a place arbitrarily close to :math:`\Lambda_+^s`. Then the contact isotopy :math:`\phi_t` basically makes a finger move that guides :math:`\Lambda_-^s` towards :math:`\Lambda_+^u`, with arbitrarily small action. In this way we have shown that there exist :math:`\pi_1`-switches between *any* pairs of :math:`\Lambda_{\pm}`. Of course, we cannot in general pre-decide the isotopy :math:`\phi_t` that actually determines tha :math:`\pi_1`-switch, for otherwise we would've arrived at the absurd conclusion that all :math:`\pi_1`-switches are free.

Finally we wrap up the discussions of free :math:`\pi_1`-switches with the following conclusion.

	Both the :ref:`trivial <section_trivial_pi1_switches>` :math:`\pi_1`-switches and the :ref:`derived <section_derived_pi1_switches>` ones can be obtained by the finger move construction without any reference to :math:`\pi_2`-switches. Nonetheless, the :math:`\pi_2^b`-switches (or the far commutativity) can be used to produce new :math:`\pi_1`-switches from existing ones even if the later are not free.

More relations between :math:`\pi_1`-switches
*********************************************

It is the nature of Morse theory that :math:`\pi_2`-switches give rise to equivalence relations between :math:`\pi_1`-switches. Indeed we have examined through the investigation of free :math:`\pi_1`-switches the equivalence relations induced by :math:`\pi_2^a` and :math:`\pi_2^b`-switches. In this section we shall go through the rest, namely the :math:`\pi_2^c` and :math:`\pi_2^d`-switches. It turns out that :math:`\pi_2^d`-switches are easier to deal with, so we shall consider them first.

Relations from :math:`\pi_2^d`-switches
+++++++++++++++++++++++++++++++++++++++

.. sidebar:: :math:`Y`-picture of a :math:`\pi_2^d`-switch

	.. figure:: static/contact-morse-theory/pi2d-morse.svg
		:align: center
		:width: 400px

To the right we have the usual :math:`2`-parameter family of :math:`\Lambda_+^s` and :math:`\Lambda_-^u` in :math:`Y` -- drawn in the front projection -- around a :math:`\pi_2^d`-:ref:`switch <bullet_pi2d_switch>`. The :math:`\pi_2^d`-switch happens exactly at the red dot where :math:`\Lambda_+^s` and :math:`\Lambda_-^u` intersect, not transversely but following :eq:`xi-degenerate`. If we let :math:`s, t` parametrize the hypersurfaces as before, then we see an equivalence relation between the :math:`I`-invariant contact structure on :math:`\Sigma^0_t, t \in I`, on the left and the composition of two :math:`\pi_1`-switches on :math:`\Sigma^1_t, t \in I`, on the right. Therefore we can say that a :math:`\pi_2^d`-switch creates or cancels a pair of :math:`\pi_1`-switches.

Note that there is nothing special about the configuration of :math:`\Lambda_+^s` and :math:`\Lambda_-^u` at either the bottom :math:`t=0` or the top :math:`t=1` of the parameter space. They can for example always be arranged using the finger moves introduced above. It means, in particular, that even more free :math:`\pi_1`-switches can be produced this way by the very existence of :math:`\pi_2^d`-switches.

Relations from :math:`\pi_2^c`-switches [#pi2c_switch_error]_
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The :math:`\pi_2^c`-:ref:`switch <bullet_pi2c_switch>` is the most special among the other :math:`\pi_2`-switches because it involves critical points *not* of the middle index. It is for this reason that we will describe it using the :math:`Y`-picture where :math:`Y` is a level set that separates :math:`p_n^+` on one side and :math:`p_{n+1}^-` together with all the :math:`p_n^-`'s it interacts with on the other side. For the sake of simplicity, let's focus on the case where the contact dimension is :math:`5`, i.e., :math:`n=2`. The general case poses in principle [#higher_dim_is_harder]_ no extra difficulty but the picture becomes much harder to draw.

.. sidebar:: :math:`Y`-picture of a :math:`\pi_2^c`-switch (schematic)

	.. figure:: static/contact-morse-theory/pi2c-morse.svg
		:align: center
		:width: 400px

In the picture to the right, we first draw a solid ball around :math:`\Lambda^u(p_{n+1}^-)` -- the red dot -- which is the intersection between the unstable manifold of :math:`p_{n+1}^-` and :math:`Y`. Note that :math:`\Lambda^u(p_{n+1}^-)` is only isotropic (and not Legendrian even though we still use the letter :math:`\Lambda`) for index reasons.

Using the usual :math:`s, t` coordinates on the parameter space, we then draw schematically the Legendrian :math:`\Lambda^s(p_n^+)` at :math:`t=0` for all :math:`s` as the thickened solid curve. As :math:`t` runs from :math:`0` to :math:`1`, the semicircle part of :math:`\Lambda^s(p_n^+)` swings upwards in the direction of the green arrows either through the left hemisphere at :math:`s=0` or through the right hemisphere at :math:`s=1`, and ends up at the dashed curve. It's not hard to imagine that at :math:`s=1/2`, the semicircle part of :math:`\Lambda^s(p_n^+)` just sweeps through the longitudinal great disk and comes across :math:`\Lambda^s(p_n^+)` at exactly :math:`t=1/2`.

Finally we draw the unstable Legendrian :math:`\Lambda_-^u`'s (the blue lines) corresponding to the :math:`p_n^-`'s that emanate flow lines towards :math:`p_{n+1}^-`. It's now obvious how a :math:`\pi_2^c`-switch gives rise to an equivalence relation between :math:`\pi_1`-switches. Namely, we have one sequence of :math:`\pi_1`-switches at :math:`s=0` on the one side, and another sequence of :math:`\pi_1`-switches at :math:`s=1` on the other side. Note, however, that it is not necessary to have a nonempty sequence on either side. In particular, when both sides are empty, we obtain a pseudo-isotopy [#pi2c_switch_error]_ between :math:`\Sigma_{t, s=0}` and :math:`\Sigma_{t, s=1}`.

The discussion so far should give us at least a conceptual understanding of the consequences of :math:`\pi_2^c`-switches. We shall wrap it up by explaining why the picture above is "schematic". It appears in the picture that :math:`\Lambda^s(p_n^+)` -- ignoring the part away from the neighborhood of :math:`\Lambda^u(p_{n+1}^-)` -- sweeps out the embedded solid ball as :math:`s, t` run through the entire parameter space. However, in reality, the region swept out by :math:`\Lambda^s(p_n^+)` is not going to be an embedded ball essentially because contact structure is not an `isotropic <https://en.wikipedia.org/wiki/Isotropic_manifold>`__ structure. Since it's neither difficult nor instructive to figure out what really should happen around a :math:`\pi_2^c`-switch, we shall leave it to a future time when it turns out to be necessary.

.. _section_flexibility:

:math:`\pi_1`-switches and flexibility
**************************************

As mentioned before, in general :math:`\pi_1`-switches do not come for free because it would otherwise equate contact topology with just smooth topology. However, there is a subclass of contact structures which carry no more information than the underlying smooth structures, or more precisely, the underlying algebraic topology if one takes into account of the almost complex structures on the contact hyperplanes. This subclass of contact structures was shown to exist and was given the name *overtwisted* contact structures by Eliashberg [Eli89]_ in dimension :math:`3` and by M\. Borman, Eliashberg, and E\. Murphy [BEM14]_ in general. Again, none of these developments are relevant to our discussion of contact Morse theory here, and we shall look for our own flexibility of contact structures from a Morse theoretic viewpoint. Of course, it won't hurt to keep in mind that it's known (to some people) that a class of contact structures are quite flexible and carry no more information than an almost complex structure on the (stable) tangent bundle.

Flexible hypersurfaces
++++++++++++++++++++++

Let's start by explaining the very first sentence of the previous paragraph. Recall that a :math:`\pi_1`-switch involves a pair of critical points :math:`p_n^{\pm}` of :math:`\Sigma|_{\xi}`, a level set :math:`\Gamma`, and respectively the unstable and stable Legendrian spheres :math:`\Lambda_{\pm} \subset \Gamma` of :math:`p_n^{\pm}`, which :math:`\xi|_{\Gamma}`-transversely intersect at one point. All these data together will be called the *initial data* of the :math:`\pi_1`-switch.

.. _flexible_hypersurface:

	 A (Morse) hypersurface :math:`\Sigma` is said to be *flexible* if the :math:`\pi_1`-switch exists in the (invariant) neighborhood of :math:`\Sigma` for any initial data, again, modulo Weinstein homotopies of :math:`R_{\pm} (\Sigma)`.

.. note::

	The above flexibility of :math:`\Sigma` doesn't quite justify the claim that the invariant contact structure on :math:`\Sigma \times I` knows nothing more than the homotopy classes of the underlying (stable) almost complex structures, i.e., overtwisted, although it should if everything works out as expected.

In dimension :math:`3`, a simple criterion due to Giroux completely characterizes when the contact germ on a (convex) surface :math:`\Sigma` is overtwisted in terms of the dividing set, but it's not particularly instructive for higher dimensional cases and has nothing to do with Morse theory. However, when interpreted in Honda's theory of bypasses, one can rephrase Giroux's criterion as saying that there exists a trivial bypass which is also overtwisted. Here a trivial bypass is the same as a trivial :math:`\pi_1`-switch but we're missing the definition of an overtwisted :math:`\pi_1`-switch. Nonetheless, this perspective provides a good motivation to our attempt at understanding flexibility of hypersurfaces in this example.

A flexible configuration
++++++++++++++++++++++++

.. sidebar:: A flexible hypersurface

	.. figure:: static/contact-morse-theory/flexible-hypersurface.svg
		:align: center
		:width: 400px

At the moment, we don't know exactly when a hypersurface :math:`\Sigma` is :ref:`flexible <flexible_hypersurface>`. But let's look at a rather special configuration on :math:`\Sigma` in terms of the decomposition :math:`\Sigma = R_+ \cup_{\Gamma} R_-` as shown on the right, which we claim to be flexible. However, this picture requires some explanation to make sense.

First of all, we need a way to describe :math:`R_{\pm}` by specifying (some of) the handles, and at the same time to describe how they are glued together along :math:`\Gamma`. The most obvious choice is to describe everything within :math:`\Gamma`. Below are the key points of such a description.

* The conditions on :math:`R_+` and :math:`R_-` are not symmetric, although their roles may be swapped. We choose to view :math:`\Gamma = \p R_-`.

* The :math:`(+1)`-labeled Legendrian unknot (red) indicates the existence of a (trivial) :math:`(n-1)`-handle in :math:`R_-`. Specifically, the unknot bounds a :math:`(n+1)`-ball (shaded) in the front projection, which actually represents a :math:`(n+1)`-sphere which is the co-core sphere of the :math:`(n-1)`-handle.

* The :math:`(-1)`-labeled Legendrian sphere (blue) indicates that :math:`\Gamma` is obtained from a sub-level set in :math:`R_-` by a :math:`n`-handle attachment along the sphere.

* The only condition on :math:`R_+` is that the (black) unknot is the co-core sphere of a :math:`n`-handle in :math:`R_+`. Note that the corresponding condition is necessarily false for :math:`R_-`. Hence :math:`R_+` and :math:`R_-` are not symmetric.

It remains to show the so configured :math:`\Sigma` is indeed flexible enough to allow for any :math:`\pi_1`-switch to exist.

.. sidebar:: Generate arbitrary :math:`\pi_1`-switches

	.. figure:: static/contact-morse-theory/flexible-rotations.svg
		:align: center
		:width: 400px

This is a fairly straightforward application of our Morse theoretic techniques for generating and manipulating :math:`\pi_1`-switches. The three major steps in showing the existence of any :math:`\pi_1`-switch is shown on the right. Let's go through them one-by-one.

The top row is the only place where our specific flexible configuration is involved. We start with the trivial :math:`\pi_1`-switch where :math:`\Lambda_-` is the standard unknot and the legitimacy of :math:`\Lambda_+` is guaranteed by our assumption above. Then we perform a :math:`\pi_2^b`-:ref:`switch <figure_simple_pi_1_switch_from_pi_2b_switch>` (aka far-commutativity) to arrive at a derived :math:`\pi_1`-switch that looks almost like the trivial one. Namely, the new :math:`\Lambda_-` remains as the unknot, but instead of sitting on top of :math:`\Lambda_+` as in the trivial case, it hangs below :math:`\Lambda_-` as shown in the rightmost figure of the first row. Such a :math:`\pi_1`-switch is known as an "overtwisted bypass" in [HH18]_.

The second row follows the first by showing, in addition, that :math:`\Lambda_+` may also be made the standard unknot. This is done by another explicit :math:`\pi_2^b`-switch such that the sought-after :math:`\pi_1`-switch associated with the blue initial data becomes trivial after the first :math:`\pi_1`-switch, which is an abstraction of the one produced in the first row.

The third row also follows from the first by another explicit :math:`\pi_2^b`-switch. Instead of keeping a "parallel copy" of :math:`\Lambda^1_-` and view the :math:`\pi_2^b`-switch in :math:`R_-` as in the second row, here we keep a "parallel copy" of :math:`\Lambda^1_+` and view it in :math:`R_+`. In this case, there are no further restrictions on :math:`\Lambda^2_-` except, of course, that it has to be the stable sphere of a negative index-:math:`n` critical point.

Combine the three moves together, one can show that :math:`\pi_1`-switches exists at the vicinity of :math:`\Sigma` for any initial data. In other words, such :math:`\Sigma` is :ref:`flexible <flexible_hypersurface>` as claimed.

Epilogue
--------

In this article, we tried to explain what contact Morse theory is about, including both the general idea on how it suppose to be used to understand contact manifolds and the main objects of interests. We then went through a limited number of examples to illustrate the most basic Morse-theoretic operations in the contact setting. Such limitation is largely due to my own limited understanding of the theory at the moment of writing. It's likely that many other mathematical articles that I'm going to write in this site will serve the purpose of improving such an understanding, and therefore they may all be called "Applications of contact Morse theory".


.. rubric:: Footnotes

.. [#gradient] Morse theory is topological in nature and doesn't care about metric very much. In particular, it's more correct and convenient but unfortunately also more cumbersome to use `gradient-like vector fields <https://en.wikipedia.org/wiki/Gradient-like_vector_field>`__ instead.

.. [#false_believe] Ironically, the failure of Giroux's argument in dimensions :math:`> 3` and other "experts' insights" went so far to even form a consensus that hypersurfaces in higher-dimensional contact manifolds are intractable and hopeless. It was at least the case when I entered the subject as a graduate student. See, for example, what my back-then-advisor had to say about `this <https://youtu.be/xuw9f4huYjk?t=2820>`__. From my own experience, there is nothing better than breaking false believes even if I was one of the believers.

.. [#contact_submfd] Regular level sets in a Morse hypersurface provide a rich source of examples of contact submanifolds. However, they don't carry in themselves much information about the original contact manifold.

.. [#pi2labels] The :math:`\pi_2`-switches are labeled by alphabetic letters at a random order because I don't have a better naming strategy.

.. [#pi2a_morse_sign] One has to make (obviously) consistent choices of signs for the critical points.

.. [#bump_switch] Details about this construction -- in fact, a discussion about the entire free :math:`\pi_1`-switch -- can be found in this `blog post <https://yhuangmath.wordpress.com/2021/05/29/trivia#bypass-everywhere>`__.

.. [#pi2c_switch_error] I once thought the :math:`\pi_2^c`-switch has only something to do with pseudo-isotopies of the Morse theory on :math:`\Sigma`. The fact that this is wrong was brought to my attention by K. Honda and J. Breen.

.. [#higher_dim_is_harder] It is likely that in practice contact topology becomes more difficult as dimension increases.
