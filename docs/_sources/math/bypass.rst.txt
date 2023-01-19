Bypass attachment revisited
===========================

.. note::

	Work in progress

Bypass attachment is an operation introduced by K. Honda in [Hon99]_ to catch the most fundamental "phase changes" among a :math:`1`-parameter family of surfaces in a contact :math:`3`-manifold. It is in effect equivalent to E. Giroux's "bifurcations" of characteristic foliations studied in [Gir99]_. The difference between the two approaches, however, is that a bypass attachment is performed by really attaching a "bypass", which is a certain half-disk with Legendrian boundary, to the surface, while a bifurcation of vector field is an abstract phenomenon that can most likely only be observed and difficult to find if we have something specific in mind. Indeed, there exist tools such as the "imbalance principle" from [Hon99]_ and the "bypass sliding" (aka. "bypass rotation") from [Hon00]_  which are designed to find bypasses.

Although the generalization of Giroux's bifurcation to hypersurfaces in higher-dimensional contact manifolds is rather straightforward [#pi1_switch]_, at least under the assumption of characteristic foliations being Morse, how a higher-dimensional bypass should look like turns out to be less obvious. A rather ugly model of such a thing -- called :math:`\Delta` -- was worked out in [HH18]_ (Figure 7.2.3 on page 45): it's more-or-less a half-ball with certain Legendrian foliation inside and a partially Legendrian boundary. The reason :math:`\p \Delta` not being completely Legendrian was that I couldn't then make it smooth, even with corners. Basically you cannot have three Legendrian planes [#leg_plane]_ intersect like a corner between the walls in your room. This psychological barrier disappears when we think of a bypass as a coLegendrian. They are not meant to be smooth: they come with cone singularities! The rest of the article will be devoted to spelling out the details. It would be a bonus if the more elegant construction of higher-dimensional bypasses will lead to a better understanding of contact manifolds, but that is unclear at this point.

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

We start by describing a standalone bypass :math:`\ps{5}{}{\Delta}` as a coLegendrian which is not attached to anything. It turns out that :math:`\psup{5}{\Delta}` is topologically just a `suspension <https://en.wikipedia.org/wiki/Suspension_(topology)>`__ of :math:`\psup{3}{\Delta}` as shown in the picture to the right [#bypass_wordpress]_. In fact, the restricted Morse vector field on :math:`\psup{5}{\Delta}` simply has one (half) source and one (half) sink, both lying on :math:`\p (\psup{5}{\Delta})`, so that all flow lines come from the source and end at the sink. The source and the sink are depicted as the blue and the green dots in the picture. Finally, the dividing set :math:`\Gamma_{\psup{5}{\Delta}}` is precisely a copy of :math:`\psup{3}{\Delta}`, which is drawn as the red half-disk.

Recall that every coLegendrian lives inside an ambient Morse hypersurface :math:`\Sigma` by definition, although it is typically unimportant how :math:`\Sigma` looks like away from the coLegendrian. Let's now put our so-far topological description of :math:`\psup{5}{\Delta}` into this context by choosing a particularly simple hypersurface :math:`\Sigma \cong S^4` such that both :math:`R_{\pm} (\Sigma) \cong B^4` are equipped with the standard (radial) Liouville vector field. In particular :math:`\Gamma_{\Sigma} \cong (S^3, \xi_{\std})`. If we now pick a bypass :math:`\psup{3}{\Delta} \subset \Gamma_{\Sigma}`, then :math:`\psup{5}{\Delta}` is simply the totality of all flow lines of the characteristic foliation :math:`\Sigma_{\xi}` passing through :math:`\psup{3}{\Delta}`, which is topologically a suspension of :math:`\psup{3}{\Delta}`.

Attachment to A Hypersurface
****************************

Our real interest lies in how a bypass :math:`\psup{5}{\Delta}` can be attached to a hypersurface :math:`\Sigma`. The basic setup for the attachment is that the ambient hypersurface at the vicinity of :math:`\psup{5}{\Delta}` -- denoted by :math:`\psup{5}{\Delta_h}` -- stands on :math:`\Sigma`, just like in the attachment of :math:`\psup{3}{\Delta}` discussed :ref:`above <fig_bypass_attachment_3d>`. In particular, we have :math:`\psup{5}{\Delta}` standing on :math:`\Sigma` along the "flat bottom", i.e., the suspension of the diameter part of :math:`\p (\psup{3}{\Delta})`, and henceforth :math:`\psup{5}{\Delta_h}` standing along a thickening of it.

.. _fig_bypass_attachment_5d:

.. sidebar:: Bypass attachment in 5D

	.. figure:: static/bypass/bypass-5d-attach.svg
		:align: center
		:width: 400px

Let's first examine how the dividing sets :math:`\Gamma_{\psup{5}{\Delta_h}}` and :math:`\Gamma_{\Sigma}` intersect. The basic principle is that the singular loci of the Legendrian foliation :math:`\Fcal` in :math:`\psup{5}{\Delta}` -- restricted to the flat bottom -- is exactly :math:`\Gamma_{\psup{5}{\Delta}} \cap \Gamma_{\Sigma}` because it is where the contact structure :math:`\xi` is tangent to :math:`\psup{5}{\Delta}` [#dividing_set_ortho]_. Since :math:`\Fcal` is the suspension of the characteristic foliation on :math:`\psup{3}{\Delta}`, it follows that :math:`\Gamma_{\psup{5}{\Delta}} \cap \Gamma_{\Sigma}` is a Legendrian :math:`\Theta`-graph in :math:`\Gamma_{\Sigma}`, and :math:`\Gamma_{\psup{5}{\Delta_h}} \cap \Gamma_{\Sigma}` is a ribbon neighborhood of it as shown in the picture to the right. Moreover, the triangles :math:`\triangle{abc}` and :math:`\triangle{abd}` -- making the flat bottom of :math:`\psup{5}{\Delta}` -- are identified with Lagrangian disks in :math:`R_{\pm} (\Sigma)` respectively. Now we know exactly how a bypass -- as a coLegendrian -- is attached to a hypersurface, at least in dimension :math:`5`. The result of such an attachment is well-understood in :ref:`contact Morse theory <contact_morse_theory_rst>`, so we shall not repeat.

.. note::

	In the above analysis, there is no assumption that the Lagrangian disks in :math:`R_{\pm} (\Sigma)` are compatible with (i.e., tangent to) the Morse vector field :math:`\Sigma_{\xi}`, and this is indeed not necessary. Bypass attachments at this level of generality are the basic building blocks of a contact Morse function (i.e., a Morse function whose gradient flow preserves :math:`\xi`) and *not* of a contact Morse theory which lives on (families of) hypersurfaces. However, for all we concern, there is no loss by assuming that all these Lagrangians are compatible with the Morse theory on :math:`\Sigma`.

What we have described above is actually how a hypersurface :math:`\Sigma` can be attached to a bypass :math:`\psup{5}{\Delta}`, but this is not how it usually works in practice. So let's turn things around and describe how :math:`\psup{5}{\Delta}` is attached to :math:`\Sigma`. All we need are the following two pieces of data on :math:`\Sigma`:

1. A Legendrian :math:`\Theta`-graph |Theta| :math:`\subset \Gamma_{\Sigma}` such that the contact planes make a half turn from one singular point of the |Theta| to the other through every one of the three connecting paths, i.e., the |upper_circle| , the |hline| , and the  |lower_circle| . Note that this is equivalent to saying that |Theta| has a ribbon neighborhood whose dividing set looks like the one shown in the :ref:`figure <fig_bypass_attachment_5d>` above.

2. The upper and lower closed semicircles of |Theta| bound Lagrangian disks in :math:`R_{\pm}`, respectively.

With these data at hand, we can then attach :math:`\psup{5}{\Delta}` to :math:`\Sigma` so that the flat bottom of :math:`\psup{5}{\Delta}` matches exactly the union of |Theta| and the two Lagrangian disks.

A Fake Bypass
*************

If all that is desired is to define a bypass and describe how it can be attached to a hypersurface, then we are already done, at least in dimension :math:`5`. But we would like to *understand* a bit more about this construction, and the best way to do so is to make mistakes. So here is one "mistake" -- which I shall call a "fake bypass" -- that had confused me a lot. Surprisingly as it turns out, we shall also learn something new from it.

.. _fig_fake_bypass_5d:

.. sidebar:: A fake bypass in 5D

	.. figure:: static/bypass/bypass-5d-fake.svg
		:align: center
		:width: 400px

The starting point is to question if the standalone bypass as constructed :ref:`above <fig_bypass_5d>` is unique. It turns out not to be the case. On the left-hand-side of the picture to the right, we have a coLegendrian :math:`\psup{5}{\Delta}^{\ast}` -- the fake bypass [#fake_bypass]_ -- with six singularities :math:`a, b, c, d, e, f` (ignore the zoomed-in half-disk about :math:`b` for the moment). The first five are of one sign and the last of the opposite sign. Using the theory of :ref:`handle decompositions <coleg_5d_building_blocks>` of coLegendrians, we can for example build such a :math:`\psup{5}{\Delta}^{\ast}` from the following partial handles [#partial_handles]_ :

1. Two (positive) :math:`0`-handles corresponding to :math:`a` and :math:`b`.
2. One (positive) half-:math:`1`-handle corresponding to :math:`e` which connects :math:`a` and :math:`b`.
3. Two positive quarter-:math:`2`-handles corresponding to :math:`c` and :math:`d`.
4. One negative :math:`3`-handle corresponding to :math:`f`.

In contrast to the actual :ref:`bypass <fig_bypass_5d>`, we note that the half-disk in the middle of :math:`\psup{5}{\Delta}^{\ast}` -- which we denote by :math:`\psup{3}{\Delta}^{\ast}` accordingly -- no longer lies in the dividing set of the ambient hypersurface containing it. Hence it makes no sense to talk about the characteristic foliation on :math:`\psup{3}{\Delta}^{\ast}`. Rather, the used-to-be characteristic foliation on :math:`\psup{3}{\Delta}^{\ast}` as drawn in the picture above is now just the intersection loci with the Legendrian foliation :math:`\Fcal^{\ast}` on :math:`\psup{5}{\Delta}^{\ast}`, which is also a suspension. It should be noted that although :math:`\Fcal^{\ast}` looks exactly like :math:`\Fcal`, they are actually different because the normal orientations along the singular loci of the Legendrian foliations are different. It follows that the contact germs on :math:`\psup{5}{\Delta}^{\ast}` and :math:`\psup{5}{\Delta}` are also different.

Another way to see that :math:`\psup{5}{\Delta}^{\ast}` and :math:`\psup{5}{\Delta}` are different is to, on the one hand, examine a linking half-disk about :math:`b \in \psup{5}{\Delta}^{\ast}`, shown as the shaded region in the :ref:`picture <fig_fake_bypass_5d>` above. Indeed, the dividing set -- or rather, the (relative) signs of the singularities -- on this half-disk makes it the opposite of a :ref:`3D bypass <fig_bypass_attachment_3d>`. On the other hand, the same link about :math:`b \in \psup{5}{\Delta}` is obviously just the usual 3D bypass.

Now comes the interesting question: can we attach :math:`\psup{5}{\Delta}^{\ast}` to a hypersurface :math:`\Sigma` in the same way that an actual bypass is attached? The answer is, interestingly enough, yes! Indeed, just as in the case of the actual bypass, all we need is a realization of a ribbon neighborhood of the Legendrian :math:`\Theta`-graph |Theta|:math:`^{\ast} \subset \Gamma_{\Sigma}` as shown in the right-hand-side of the :ref:`picture <fig_fake_bypass_5d>` above. At first sight, it might appear that such a ribbon neighborhood cannot exist in :math:`\Gamma_{\Sigma}` if, for example, the outer-circle |circle| :math:`\subset` |Theta|:math:`^{\ast}` is topologically unknotted because its Thurston-Bennequin number :math:`\op{tb}(`\ |circle|\ :math:`) = 0` relative to the ribbon framing. But this is an illusion because the ribbon framing does not have to be the same as the `Seifert framing <https://en.wikipedia.org/wiki/Self-linking_number>`__, if the later is defined at all.

.. important::

	When we attach either a bypass or a fake bypass to a hypersurface, the choice of the ribbon neighborhood of the Legendrian |Theta| in :math:`\Gamma_{\Sigma}` is not completely arbitrary, although :math:`\op{tb}(`\ |circle|\ :math:`)` is :math:`-1` for the bypass attachment and :math:`0` for the fake bypass attachment. Namely, the :math:`\op{tb}` of both the upper and lower closed semicircles of |Theta| must be :math:`-1`. This guarantees that the ribbon framing extends to the Lagrangian disks in :math:`R_{\pm}`, respectively.

The above restriction is quite strong but there is yet another possibility for the ribbon neighborhood where :math:`\op{tb}(`\ |circle|\ :math:`) = -2`. It is not the purpose of this article to exhaust all possibilities of bypass-like objects. But the discussion of the fake bypass has hopefully deepened our understanding of the bypass itself.

Back to the Old Bypass
**********************

tbc...


.. rubric:: Footnotes

.. [#pi1_switch] It is a handle slide in contact Morse theory as illustrated in this :ref:`picture <figure_r_pm_picture_of_pi_1_switch>`.

.. [#leg_plane] Thinking inside of an ambient contact :math:`5`-manifold.

.. [#bypass_subscript] The superscript in :math:`\psup{3}{\Delta}` indicates that the ambient contact dimension is :math:`3`. It could be omitted when we talk about bypasses in any dimension. In general :math:`\dim (\psup{2n+1}{\Delta}) = n+1` because it's a coLegendrian.

.. [#bypass_wordpress] A very similar-looking -- but wrong -- picture of a bypass was drawn in `one <https://yhuangmath.wordpress.com/2021/07/20/flexibility-of-legendrian-2-spheres-in-contact-5-space-iii/>`__ of my earlier wordpress posts, where I tried to attach it to a Legendrian instead of a hypersurface. While there is nothing wrong attaching a bypass to a Legendrian, it turns out to be more fundamental to understand how it is attached to a hypersurface. Every time I found something out, I found also that I had come close to the correct answer many times without realizing. I cannot tell if it is me being stupid or it is just the way it is, but I have enjoyed this process nonetheless.

.. [#dividing_set_ortho] It is often helpful, though not necessary, to think of the dividing set as the locus where the contact hyperplanes are perpendicular to the hypersurface (wrt. an auxillary Riemannian metric).

.. [#fake_bypass] The fake bypass seems to be the more natural candidate of a bypass because it looks more like the :ref:`3D bypass <fig_bypass_attachment_3d>` :math:`\psup{3}{\Delta}` in terms of their dividing sets.

.. [#partial_handles] The handles are necessarily partial because the singularities all lie on :math:`\p (\psup{5}{\Delta}^{\ast})`. For handles of middle indexes, i.e., neither :math:`0` nor :math:`\dim (\psup{5}{\Delta}^{\ast})`, we shall indicate whether it is a half or a quarter handle.


.. Unicode substitutions
.. |Theta| unicode:: U+2296
.. |circle| unicode:: U+25CB
.. |upper_circle| unicode:: U+25E0
.. |lower_circle| unicode:: U+25E1
.. |hline| unicode:: U+23AF
