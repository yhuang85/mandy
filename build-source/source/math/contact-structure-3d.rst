What's a contact structure in 3D?
=================================

We try to explain what a contact structure is in a :math:`3`-dimensional space from an intrinsic viewpoint, i.e., we pretend to live inside the contact 3-space without any knowledge of what an outside universe means. This article is heavily influenced by M. Gromov's paper [Gro92]_.


Space and Dimension
-------------------

Contact structure is a structure living in a space. So let's first talk about spaces.

It may be either the choice of god or the perception of our kind that we all live in a space. Since the time we looked into the sky, the earth, or the ocean, we never stopped wondering, physically or mentally, about the mysterious space around us.

When talking about spaces, one of the first things that come to is the so-called *dimension*, or in more plain terms: how big is it? There is no universally agreed definition of what a dimension is. There really shouldn't be one since every definition has to serve a particular purpose and no definition makes sense as a standalone. Here are two definitions of dimensions that are relevant to contact structures.

Visual dimension -- Degrees of Freedom
	This way of defining dimension is particularly fonded by physicists and engineers. Roughly speaking, it asks about in how many independent directions can one moves around. In this sense, we all seem to live in a :math:`3`-dimensional space since there are apparently :math:`3` degrees of freedom: up :math:`\leftrightarrow` down, left :math:`\leftrightarrow` right, and forward :math:`\leftrightarrow` backward.

Measure dimension -- Volume Growth
	The other, slightly less intuitive, way to define dimension is to ask how fast does the volume of a ball grow with the radius. Namely, if :math:`\op{Vol}(B(r))` denotes the volume of a ball of radius :math:`r`, then the asymptotic :math:`\op{Vol}(B(r)) \sim r^d` as :math:`r \to 0` defines the measure dimension to be :math:`d`.

		The faster the volume grows, the larger the dimension is.

	To make some sense of this statement, imagine the space as consisting of an immensely dense population of particles. Then the volume of a ball can be thought of as counting the number of particles in the ball. Now a fast volume growth means that with a unit increment of the radius, the ball suddenly contains much more particles. This is a good indication that the space is "large".

.. note::

	1. Another intriguing definition of dimension, named after F. Hausdorff, aims at measuring the size of `fractals <https://en.wikipedia.org/wiki/Fractal>`_. We shall make no use of that since contact structures do not live on fractals.
	2. The volume growth of balls is by itself a very interesting and counter-intuitive phenomenon, especially when the dimension is large [#three_not_big]_ . Indeed, the mass of a high dimensional ball is highly concentrated near its boundary, as opposed to the interior which would have followed from common sense. This is known as the `concentration of measure <https://en.wikipedia.org/wiki/Concentration_of_measure>`_.


Homogeneity of Space
--------------------

All the above definitions of dimensions are local in nature, i.e., The dimension measures the size of a space at the vicinity of every point. It is therefore entirely possible that the dimension of a space jumps from point to point, and it doesn't even have to be continuous. Given a space :math:`V` and a point :math:`p \in V`, let's write :math:`\dim_p(V)` for the dimension of :math:`V` at :math:`p`, whose definition may be any one of the above.

.. note::

	Whenever we write just :math:`\dim(V)`, we assume implicitly that the dimension is the same at all points, i.e., the space is in a sense homogeneous.


Agreement of Dimensions
-----------------------

For a space :math:`V`, let's write :math:`\op{Vdim}(V)` for its visual dimension and :math:`\op{Mdim}(V)` for its measure dimension. For "usual" spaces :math:`V`, one can expect :math:`\op{Vdim}(V) = \op{Mdim}(V)`. For example, we have :math:`\op{Vdim}(\RR^3)=3` obviously, and also :math:`\op{Mdim}(\RR^3)=3` since :math:`\op{Vol}(B^3(r)) = \tfrac{4}{3}\pi r^3`.

.. note::

	It's in fact irrelevant what exactly the shape one uses to compute the measure dimension. We used balls above probably because it's the most intuitive shape, but it's not the shape whose volume is the easiest to compute. Indeed, it's easier to use cubes :math:`K(r)` of side length :math:`r` instead and refer to the much simpler :math:`\op{Vol}(K(r)) = r^3`.

It's not immediately obvious what kind of spaces would have unequal visual and measure dimensions. Contact structures, especially in dimension :math:`3`, provide good examples of such "weird" spaces. Indeed, we're not even interested in the underlying topology of the space, which will always be :math:`\RR^3`.


Distribution and Conductivity
-----------------------------

In :math:`\RR^3`, one can freely walk in all :math:`3` independent directions if no further constraints are enforced. However, it's not hard to imagine that a live being, for whatever reason, is disallowed to move in certain directions. From now on this live being will be called Mandy. Then for all Mandy knows, the visual dimension, aka the degree of freedom, may well be less than :math:`3` even though from the viewpoint of an outsider, the space itself has not changed at all.

Constraints like these are called distributions. More formally, a distribution :math:`D` on :math:`\RR^3` is a choice of subspace at every point, i.e., :math:`D(x) \subset T_x \RR^3 \cong \RR^3` is a linear subspace for every :math:`x \in \RR^3`. Then dimension :math:`\dim D(x)` is called the *rank* of :math:`D` at :math:`x`. To make life a bit easier, we will always work under the following two assumptions:

- The rank of :math:`D` is constant, and will be called :math:`\op{rank}(D)`.

- :math:`D` varies smoothly from point to point.

The only potentially interesting scenarios are :math:`\op{rank}(D) = 1` or :math:`2`. The :math:`\op{rank}(D) = 1` case is also not particularly interesting for our purposes here since Mandy, in this case, will just live either on a line or on a circle, and the actual space :math:`\RR^3` becomes irrelevant. So let's assume :math:`\op{rank}(D) = 2` from now on. Moreover, we want to make sure that Mandy is able to explore the entire :math:`\RR^3` under the constraint :math:`D`, which we translate into the following *conductivity* condition.

Conductivity of :math:`D`
	For any two points :math:`p, q \in \RR^3`, there exists an integral curve :math:`\gamma`, i.e., a curve everywhere tangent to :math:`D`, connecting :math:`p` and :math:`q`.

.. note::

	Since :math:`\RR^3` is connected, the conductivity of :math:`D` is equivalent to a *local* version where the points :math:`p, q` can be required to be arbitrarily close to each other. This suits better to our discussions here since dimension is also a local quantity.

Obviously not all distributions are conductive, for example, the distribution of planes all parallel to the :math:`xy`-plane is not. The following distribution, on the other hand, is quite conductive

.. math::

	D_{\std} \coloneqq \ker(dz-ydx) = \op{span} (\p_x + y\p_z, \p_y).

Indeed, for any :math:`p, q \in \RR^3`, consider the projections :math:`\pi(p), \pi(q) \in \RR^2_{xy}`, where :math:`\pi: \RR^3 \to \RR^2` is defined by :math:`\pi(x, y, z) = (x, y)`.

Pick a (parametrized) curve :math:`\bar{\gamma}: [0, 1] \to \RR^2_{xy}` such that :math:`\bar{\gamma}(0) = \pi(p)` and :math:`\bar{\gamma}(1) = \pi(q)`. Then :math:`\bar{\gamma}` has a unique integral lift :math:`\gamma: [0, 1] \to \RR^3` given the initial condition :math:`\gamma(0) = p` in the sense that :math:`\bar{\gamma} = \pi \circ \gamma`. This is plainly because :math:`\pi_{\ast}: D_{\std} \to T\RR^2` is a pointwise isomorphism. More explicitly, we can solve for :math:`z` along :math:`\gamma` by

.. math::

	z(t) = \int_0^t y(s)x'(s)ds + z(0), \quad t \in [0, 1],

where :math:`(x(t), y(t), z(t)) \coloneqq \gamma(t)` and :math:`z(0) = z(p)`. Now to make sure that :math:`z(1) = z(q)`, one just need to choose :math:`\bar{\gamma}` such that

.. math::

	\int_0^1 y(s)x'(s)ds = z(q) - z(p),

and there are simply plenty of such curves [#stokes_1]_ .

The above distribution :math:`D_{\std}` is an example of a contact structure and it's special enough to be called the *standard contact structure*.

.. note::

	The projection :math:`\pi: \RR^3 \to \RR^2` above, known as the *Lagrangian projection*, is a quite helpful tool in understanding contact structures in general essentially because of the unique lifting property mentioned above. It's an important bridge between contact and symplectic structures.


Disagreement of Dimensions
--------------------------

Obviously :math:`\op{Vdim} (\RR^3, D_{\std}) = 2` and it's our goal here to compute :math:`\op{Mdim} (\RR^3, D_{\std})`. To this end, we need a metric on :math:`\RR^3` so that we can measure length and volume. The particular choice of metric is not super important here and we will stick to the Euclidean metric for simplicity.

We shall compute :math:`\op{Mdim}(\RR^3, D_{\std})` near the origin as follows. Consider a small cube

.. math::

	K(r) \coloneqq [0, r_x] \times [0, r_y] \times [0, r_z] \subset \RR^3

of Euclidean side lengths :math:`r_x, r_y, r_z`. We'd like it to be a square cube of side length :math:`r` when measured in :math:`(\RR^3, D_{\std})` [#measure]_ . Namely, if we write :math:`X = (r, 0, 0), Y = (0, r, 0), Z = (0, 0, r)` be the three vertices adjacent to the origin :math:`O`, then we need

.. math::

	\dist(O, X) = \dist(O, Y) = \dist(O, Z) = r,

where the distance between two points is defined to be the length of the shortest path connecting them. Once this is achieved, the measure dimension is simply the growth rate of :math:`\op{Vol} (K(r)) = r_x r_y r_z` [#cube]_ as :math:`r \to 0`.

Quite obviously :math:`\dist(O, X) = r_x` and :math:`\dist(O, Y) = r_y`. Therefore :math:`r_x = r_y = r`. But the computation of :math:`\dist(O, Z)` requires some effort since the straight interval between :math:`O` and :math:`Z` is not tangent to :math:`D_{\std}`. Let

.. math::

	\gamma: [0, 1] \to \RR^3: t \mapsto (x(t), y(t), z(t))

be an integral curve from :math:`O` to :math:`Z`. Then the length :math:`\ell(\gamma)` of :math:`\gamma` can be computed as follows

.. math::

	\ell(\gamma) = \int_0^1 \left( \sqrt{\dot{x}^2 + \dot{y}^2 + \dot{z}^2} \right) dt = \int_0^1 \left( \sqrt{(1 + y^2) \dot{x}^2 + \dot{y}^2} \right) dt

Let :math:`\bar{\gamma} \coloneqq \pi \circ \gamma: [0, 1] \to \RR^2` as before. Then there is a constant :math:`K > 0` such that

.. math::

	\tfrac{1}{K} \cdot \ell(\bar{\gamma}) < \ell(\gamma) < K \cdot \ell(\bar{\gamma}).

Hence as far as the asymptotic is concerned :math:`\ell(\gamma) \sim \ell(\bar{\gamma})` as :math:`r \to 0`.

Observe that :math:`\bar{\gamma}` is a closed loop since :math:`\pi (O) = \pi (Z) = 0`. Hence the minimization of :math:`\ell(\gamma)`, up to a bounded error, can be translated into the problem of maximizing :math:`\ell(\bar{\gamma})` while keeping the area :math:`A` enclosed by :math:`\bar{\gamma} \subset \RR^2` constant [#stokes_2]_ . This is a classical isoperimetric problem, and by forgetting an irrelevant constant, we have

.. math::

	\min_{A = r_z} \ell (\bar{\gamma}) \sim \sqrt{r_z}

It follows that :math:`r_z \sim r^2`, and we conclude that

.. math::

	\op{Vol} (K(r)) = r_x r_y r_z \sim r^4.

Hence :math:`\op{Mdim} (\RR^3, D_{\std}) = 4`, at least near the origin.


Homogeneity of Contact Structures
---------------------------------

It turns out that contact structures are far more homogeneous than just having constant measure dimension. Indeed, a `theorem <https://en.wikipedia.org/wiki/Darboux%27s_theorem>`_ of Darboux asserts that any contact structure (which we haven't defined yet) is locally isomorphic to the standard :math:`(\RR^3, D_{\std})`, hence the name. It follows that :math:`\op{Mdim} (\RR^3, D_{\std}) = 4` after all [#darboux_redundant]_ .


More Distributions
------------------

Inspired by the definition of :math:`D_{\std}`, let's consider the following sequence of distributions

.. math::

	D_k \coloneqq \ker (dz - y^k dx), \quad k \geq 1,

where :math:`D_1 = D_{\std}` of course. One can verify they are all conductive, and our goal here is compute their measure dimensions.

Let's first compute the measure dimension at the origin as before. Indeed, the only difference in such a computation between different :math:`k` is the distance :math:`\dist(O, Z)`. Specifically, we need to solve the following optimization problem:

	With the 'area'

	.. math::

		A \coloneqq \iint_{\Delta (\bar{\gamma})} k y^{k-1} dxdy = r_z

	fixed, minimize the length

	.. math::

		\ell \coloneqq \int_0^1 \left( \sqrt{\dot{x}^2 + \dot{y}^2} \right) dt,

	where :math:`\Delta (\bar{\gamma})` denotes the area enclosed by :math:`\bar{\gamma} \subset \RR^2_{xy}`.

This is not particularly easy, but we don't really need to solve it literally either. All we want to know is the growth rate of :math:`\min \ell` as :math:`r_z \to 0`. To this end, it suffices to consider a rectangular :math:`\bar{\gamma}` whose sides are parallel to the axes. Of course, the lower-left corner of :math:`\bar{\gamma}` is the origin.

.. note::

	By the same token, the earlier reference to the isoperimetric inequality was unnecessary.

Let :math:`h, w` be the height and width of the rectangle, respectively. Then :math:`A = h^k w` and :math:`\ell = 2(h+w)`. This much simplified optimization problem can be easily solved to yield

.. math::

	\ell \sim r_z^{1/(k+1)}

as :math:`r_z \to 0`. We conclude that :math:`\op{Mdim}_0 (\RR^3, D_k) = k+3`. This is interesting since the measure dimension captures the power :math:`k`.

It turns out that :math:`(\RR^3, D_k)` is not homogeneous when :math:`k>1`. In particular, the measure dimension equals :math:`k+3` at points with :math:`y = 0` and equals :math:`4` as in the case of contact structures otherwise. We omit the details here as the computations are exactly the same.


What's a Contact Structure Anyway?
----------------------------------

In fact, the distributions :math:`D_k, k > 1`, are locally isomorphic to :math:`D_{\std}` on the region where :math:`y \neq 0`. Moreover, it doesn't seem possible to have a distribution on :math:`\RR^3` whose measure dimension is *everywhere* greater than :math:`4`.

Based on these observations, we can finally propose an intrinsic characterization of contact structures as follows:

A contact structure on a :math:`3`-dimensional space is a rank-:math:`2` distribution, which makes the measure dimension homogeneous and greater than :math:`3` everywhere.

Of course, one can easily find formal definition of `contact structures <https://en.wikipedia.org/wiki/Contact_geometry>`_. But the point of the above heuristic (which may not even be correct!) is that if you, at some point in your life [#drunk]_ ,  start suspecting whether you live in a contact world because you can only move in :math:`2` directions instead of the usual :math:`3`, you can find it out simply by measuring volume of tiny little balls.


.. rubric:: Footnotes

.. [#three_not_big] Dimension 3 is not large enough.

.. [#stokes_1] The easiest way to find such a curve is to invoke Stokes' theorem.

.. [#cube] A cube of side length :math:`r` is just a convenient substitute for a ball of radius :math:`r`.

.. [#stokes_2] By Stokes' theorem :math:`A = z(C) - z(O) = r_z`.

.. [#darboux_redundant] The reference to Darboux's theorem is a perfect example of laziness as one can directly verify the constancy of :math:`\op{Mdim} (\RR^3, D_{\std})` by essentially the same computation.

.. [#measure] Here we're using the Euclidean measure.

.. [#drunk] Most likely when you are sufficiently drunk.
