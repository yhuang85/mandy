Bypass attachment revisited
===========================

.. note::

	Work in progress

Bypass attachment is an operation introduced by K. Honda in [Hon99]_ to catch the most fundamental "phase changes" among a :math:`1`-parameter family of surfaces in a contact :math:`3`-manifold. It is in effect equivalent to E. Giroux's "bifurcations" of characteristic foliations studied in [Gir99]_. The difference between the two approaches, however, is that a bypass attachment is performed by really attaching a "bypass", which is a certain half-disk with Legendrian boundary, to the surface, while a bifurcation of vector field is an abstract phenomenon that can most likely only be observed and difficult to find if we have something specific in mind. Indeed, there exist tools such as the "imbalance principle" from [Hon99]_ and the "bypass sliding" (aka. "bypass rotation") from [Hon00]_  which are designed to find bypasses.

Although the generalization of Giroux's bifurcation to hypersurfaces in higher-dimensional contact manifolds is rather straightforward [#pi1_switch]_, at least under the assumption of characteristic foliations being Morse, how a higher-dimensional bypass should look like turns out to be less obvious. A rather ugly model of such a thing -- called :math:`\Delta` -- was worked out in [HH18]_ (Figure 7.2.3 on page 45): it's more-or-less a half-ball with certain Legendrian foliation inside and a partially Legendrian boundary. The reason :math:`\p \Delta` not being completely Legendrian was that I couldn't then make it smooth, even with corners. Basically you cannot have three Legendrian planes [#leg_plane]_ intersect like a corner between the walls in your room. This psychological barrier disappears when we think of a bypass as a coLegendrian. They are not meant to be smooth: they come with cone singularities! The rest of the article will be devoted to spelling out the details. It would be a bonus if the more elegant construction of higher-dimensional bypasses will lead to a better understanding of contact manifolds, but that is highly unclear at this point.

.. note::

	The rest of this article assumes familiarity with Honda's bypass theory in dimension :math:`3` and coLegendrians and contact Morse theory in general in higher dimensions.

Bypass in Dimension :math:`3`
-----------------------------

.. _fig_bypass_attachment_3d:

.. sidebar:: Bypass attachment in 3D

	.. figure:: static/bypass/bypass-3d.svg
		:align: center
		:width: 400px

We start by reviewing the definition of a bypass in dimension :math:`3`. In the picture to the right, a bypass :math:`\psup{3}{\Delta}` [#bypass_subscript]_ -- topologically a half-disk -- is standing on a (generic) Morse surface :math:`\Sigma` along the diameter. Moreover, the characteristic foliation on :math:`\psup{3}{\Delta}` is drawn, together with four singularities divided into two groups with opposite signs: three of one sign (green) and the other of the opposite sign (blue). Finally, the dividing set :math:`\Gamma` on both :math:`\psup{3}{\Delta}` and :math:`\Sigma` are drawn in red.

Now the bypass :math:`\psup{3}{\Delta}` can be "attached" to :math:`\Sigma` to form a new surface :math:`\Sigma'`, which is topologically isotopic to :math:`\Sigma` but contact-topologically different in general. The original recipe for such an "attachment" as described in [Hon99]_ is to first slightly thicken :math:`\psup{3}{\Delta}` and then round the corners. This procedure is then realized to be equivalent to first attaching a (contact) :math:`1`-handle whose core-disk is the semicircle part of :math:`\p (\psup{3}{\Delta})`, and then a cancelling :math:`2`-handle whose core-disk is really a copy of :math:`\psup{3}{\Delta}` itself. We might on the one hand prefer the later viewpoint because it's in line with Morse theory, but on the other hand, we don't really care how bypasses are attached as long as there is some recipe, and in any case, we already know everything about the contact structure of a bypass attachment via the contact Morse theory on :math:`\Sigma`.

To wrap up our recollection of the bypass theory in dimension :math:`3`, let's remind ourselves that in this dimension, the dividing sets play a dominant role due to Giroux's convex surface theory. In particular, it doesn't really matter how the characteristic foliation on :math:`\psup{3}{\Delta}` looks like as long as the dividing set :math:`\Gamma_{\psup{3}{\Delta}}` is correct. This is *not* the case in higher dimensions, where dividing sets contain much less information. Indeed, our choice of the characteristic foliation on :math:`\psup{3}{\Delta}` is deliberate and important: it is (generically) the characteristic foliation on the core-disk of a :math:`2`-handle.

.. important::

	We will assume the characteristic foliation on :math:`\psup{3}{\Delta}` is the standard one as shown in the :ref:`figure <fig_bypass_attachment_3d>` above.

Bypass in Dimension :math:`5`
-----------------------------

Instead of moving onto bypasses in any dimension, let's first try to understand it in dimension :math:`5` for two reasons: first, it's easier to visualize a :math:`3`-dimensional bypass and second, the even higher dimensional cases are not any harder once we understand the :math:`5`-dimensional case.

Construction as A Standalone
****************************

.. _fig_bypass_5d:

.. sidebar:: Bypass in 5D

	.. figure:: static/bypass/bypass-5d.svg
		:align: center
		:width: 400px

We start by describing a standalone bypass :math:`\ps{5}{}{\Delta}` as a coLegendrian which is not attached to anything. It turns out that :math:`\psup{5}{\Delta}` is topologically just a `suspension <https://en.wikipedia.org/wiki/Suspension_(topology)>`__ of :math:`\psup{3}{\Delta}` as shown in the picture to the right. In fact, the restricted Morse vector field on :math:`\psup{5}{\Delta}` simply has one (half) source and one (half) sink, both lying on :math:`\p (\psup{5}{\Delta})`, so that all flow lines come from the source and end at the sink. The source and the sink are depicted as the blue and the green dots in the picture. Finally, the dividing set :math:`\Gamma_{\psup{5}{\Delta}}` is precisely a copy of :math:`\psup{3}{\Delta}`, which is drawn as the red half-disk.

Recall that every coLegendrian lives inside an ambient Morse hypersurface :math:`\Sigma` by definition, although it is typically unimportant how :math:`\Sigma` looks like away from the coLegendrian. Let's now put our so-far topological description of :math:`\psup{5}{\Delta}` into this context by choosing a particularly simple hypersurface :math:`\Sigma \cong S^4` such that both :math:`R_{\pm} (\Sigma) \cong B^4` are equipped with the standard (radial) Liouville vector field. In particular :math:`\Gamma_{\Sigma} \cong (S^3, \xi_{\std})`. If we now pick a bypass :math:`\psup{3}{\Delta} \subset \Gamma_{\Sigma}`, then :math:`\psup{5}{\Delta}` is simply the totality of all flow lines of the characteristic foliation :math:`\Sigma_{\xi}` passing through :math:`\psup{3}{\Delta}`, which is topologically a suspension of :math:`\psup{3}{\Delta}`.

Attachment to A Hypersurface
****************************

Our real interest lies in how a bypass :math:`\psup{5}{\Delta}` can be attached to a hypersurface :math:`\Sigma`. The basic setup for the attachment is that the ambient hypersurface at the vicinity of :math:`\psup{5}{\Delta}` -- denoted by :math:`\psup{5}{\Delta}_h` -- stands on :math:`\Sigma`, just like in the attachment of :math:`\psup{3}{\Delta}` discussed :ref:`above <fig_bypass_attachment_3d>`. In particular, we have :math:`\psup{5}{\Delta}` standing on :math:`\Sigma` along the "flat bottom", i.e., the suspension of the diameter part of :math:`\p (\psup{3}{\Delta})`, and henceforth :math:`\psup{5}{\Delta}_h` standing along a thickening of it.

.. _fig_bypass_attachment_5d:

.. sidebar:: Bypass attachment in 5D

	.. figure:: static/bypass/bypass-5d-attach.svg
		:align: center
		:width: 400px

Let's first examine how the dividing sets :math:`\Gamma_{\psup{5}{\Delta}_h}` and :math:`\Gamma_{\Sigma}` intersect. The basic principle is that the singular loci of the Legendrian foliation :math:`\Fcal` in :math:`\psup{5}{\Delta}` -- restricted to the flat bottom -- is exactly :math:`\Gamma_{\psup{5}{\Delta}} \cap \Gamma_{\Sigma}` because it is there where the contact structure :math:`\xi` is tangent to :math:`\psup{5}{\Delta}_h` [#dividing_set_ortho]_. Since :math:`\Fcal` is the suspension of the characteristic foliation on :math:`\psup{3}{\Delta}`, it follows that :math:`\Gamma_{\psup{5}{\Delta}} \cap \Gamma_{\Sigma}` is a Legendrian :math:`\Theta`-graph in :math:`\Gamma_{\Sigma}`, and :math:`\Gamma_{\psup{5}{\Delta}_h} \cap \Gamma_{\Sigma}` is a ribbon neighborhood of it as shown in the picture to the right. Moreover, the triangles :math:`\triangle{abc}` and :math:`\triangle{abd}` -- making the flat bottom of :math:`\psup{5}{\Delta}` -- are identified with Lagrangian disks in :math:`R_{\pm} (\Sigma)` respectively. Now we know exactly how a bypass -- as a coLegendrian -- is attached to a hypersurface, at least in dimension :math:`5`. The result of such an attachment is of course well-understood in contact Morse theory, so we shall not repeat.

.. note::

	In the above analysis, there is no assumption that the Lagrangian disks in :math:`R_{\pm} (\Sigma)` are compatible with (i.e., tangent to) the Morse vector field :math:`\Sigma_{\xi}`, and this is indeed not necessary. Bypass attachments at this level of generality are the basic building blocks of a contact Morse function (i.e., a Morse function whose gradient flow preserves :math:`\xi`) and *not* of a :ref:`contact Morse theory <contact_morse_theory_rst>` which lives on (families of) hypersurfaces. However, for all we concern, there is no loss by assuming that all these Lagrangians are compatible with the Morse theory on :math:`\Sigma`.

tbc...

.. rubric:: Footnotes

.. [#bypass_subscript] The superscript in :math:`\psup{3}{\Delta}` indicates that the ambient contact dimension is :math:`3`. It could be omitted when we talk about bypasses in any dimension. In general :math:`\dim (\ps{2n+1}{}{\Delta}) = n+1` because it's a coLegendrian.

.. [#pi1_switch] It is a handle slide in contact Morse theory as illustrated in this :ref:`picture <figure_r_pm_picture_of_pi_1_switch>`.

.. [#leg_plane] Thinking inside of an ambient contact :math:`5`-manifold.

.. [#dividing_set_ortho] It is often helpful, though not necessary, to think of the dividing set as the locus where the contact hyperplanes are perpendicular to the hypersurface (wrt. an auxillary Riemannian metric).
