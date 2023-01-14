Bypass attachment revisited
===========================

.. note::

	Work in progress

Bypass attachment is an operation introduced by K. Honda in [Hon99]_ to catch the most fundamental "phase changes" among a :math:`1`-parameter family of surfaces in a contact :math:`3`-manifold. It is in effect equivalent to E. Giroux's "bifurcations" of characteristic foliations studied in [Gir99]_. The difference between the two approaches, however, is that a bypass attachment is performed by really attaching a "bypass", which is a certain half-disk with Legendrian boundary, to the surface, while a bifurcation of vector field is an abstract phenomenon that can most likely only be observed and difficult to find if we have something specific in mind. Indeed, there exist tools such as the "imbalance principle" from [Hon99]_ and the "bypass sliding" (aka. "bypass rotation") from [Hon00]_  which are designed to find bypasses.

Although the generalization of Giroux's bifurcation to hypersurfaces in higher-dimensional contact manifolds is rather straightforward [#pi1_switch]_, at least under the assumption of characteristic foliations being Morse, how a higher-dimensional bypass should look like turns out to be less obvious. A rather ugly model of such a thing, called :math:`\Delta`, was worked out in [HH18]_ (Figure 7.2.3 on page 45): it's more-or-less a half-ball with certain Legendrian foliation inside and a partially Legendrian boundary. The reason :math:`\p \Delta` not being completely Legendrian was that I couldn't then make it smooth, even with corners. Basically you cannot have three Legendrian planes [#leg_plane]_ intersect like a corner between the walls in your room. This psychological barrier disappears when we think of a bypass as a coLegendrian. They are not meant to be smooth: they come with cone singularities. The rest of the article will be devoted to spelling out the details. It would be a bonus if the more elegant construction of higher-dimensional bypasses will lead to a better understanding of contact manifolds, but that is highly unclear at this point.

.. note::

	The rest of this article assumes familiarity with Honda's bypass theory in dimension :math:`3` and coLegendrians and contact Morse theory in general in higher dimensions.

Bypass in Dimension :math:`3`
-----------------------------

.. sidebar:: Bypass in 3D

	.. figure:: static/bypass/bypass-3d.svg
		:align: center
		:width: 400px

We start by reviewing the definition of a bypass in dimension :math:`3`. In the picture to the right, a bypass :math:`\Delta` -- topologically a half-disk -- is standing on a (generic) Morse surface :math:`\Sigma` along the diameter. Moreover, the characteristic foliation on :math:`\Delta` is drawn, together with four singularities divided into two groups with opposite signs: three of one sign (green) and the other of the opposite sign (blue). Finally, the dividing set :math:`\Gamma` on both :math:`\Delta` and :math:`\Sigma` are drawn in red.

Now the bypass :math:`\Delta` can be "attached" to :math:`\Sigma` to form a new surface :math:`\Sigma'`, which is topologically isotopic to :math:`\Sigma` but contact-topologically different in general. The original recipe for such an "attachment" as described in [Hon99]_ is to first slightly thicken :math:`\Delta` and then round the corners. This procedure is then realized to be equivalent to first attaching a (contact) :math:`1`-handle whose core-disk is the semicircle part of :math:`\p \Delta`, and then a cancelling :math:`2`-handle whose core-disk is really a copy of :math:`\Delta` itself. We might on the one hand prefer the later viewpoint because it's in line with Morse theory, but on the other hand, we don't really care how bypasses are attached as long as there is some recipe, and in any case, we already know everything about the contact structure of a bypass attachment via the contact Morse theory on :math:`\Sigma`.

To wrap up our recollection of the bypass theory in dimension :math:`3`, let's remind ourselves that in this dimension, the dividing sets play a dominant role due to Giroux's convex surface theory. In particular, it doesn't really matter how the characteristic foliation on :math:`\Delta` looks like as long as the dividing set :math:`\Gamma_{\Delta}` is correct. This is *not* the case in higher dimensions, where dividing sets contain much less information. Indeed, our choice of the characteristic foliation on :math:`\Delta` is deliberate and important: it is (generically) the characteristic foliation on the core-disk of a :math:`2`-handle.

Bypass in Higher Dimensions
---------------------------

tbc...

.. rubric:: Footnotes

.. [#pi1_switch] It is a handle slide in contact Morse theory as illustrated in this :ref:`picture <figure_r_pm_picture_of_pi_1_switch>`.

.. [#leg_plane] Thinking inside of an ambient contact :math:`5`-manifold.
