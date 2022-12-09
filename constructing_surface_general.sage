from flatsurf import *




# Model surface O=O(d_1,...,d_k), the disjoint union of O_{d_i} for i=1,...,k
class ModelSurface:
	'''Model surface O=O(d_1,...,d_k), the disjoint union of infinite translation surfaces O_{d_i}=(C,z^{d_i}dz) for i=1,...,k.  Input orders is a list [d_1,...,d_k] of orders of zeros.'''
	def __init__(self,orders):
		assert sum(orders)%2==0
		orders.sort()
		self.orders_of_zeros=orders
		self.components=[ModelSurfaceComponent(d,self) for d in self.orders_of_zeros]
		cardinalities=[]
		while len(orders)>0:
			cardinalities.append(orders.count(orders[0]))
			orders=orders[orders.count(orders[0]):]
		self.cardinalities_of_distinct_zeros=cardinalities
		# Keep track of which marked segments contribute to edges of convex bodies
		self.convex_bodies_pts=None		

		# The kernel of der: Aff^+(O) -> GL^+(2,R), represented as the product S x C, where S and C are as follows: 
		### S is the product of permutation groups S_{n_1},...,S_{n_t}, where the n_i's are cardinalities of distinct orders of zeros d_1,...,d_k
		### C is the product of cyclic permutation groups C_{d_1},...,C_{d_k}
		self.Trans=cartesian_product([cartesian_product([SymmetricGroup(n) for n in self.cardinalities_of_distinct_zeros]),cartesian_product([CyclicPermutationGroup(d+1) for d in self.orders_of_zeros])])

	def __repr__(self):
			return ''.join(['Model surface O=O(',str(self.orders_of_zeros)[1:-1],'), the disjoint union of ',''.join(['O_{}, '.format(d) for d in self.orders_of_zeros[0:-1]]+['O_{}'.format(self.orders_of_zeros[-1])])])

	def add_orientation_paired_marked_segments(self,holonomy,component_index,plane_index,inverse_component_index,inverse_plane_index):
		'''Adds a marked segment with holonomy vector 'holonomy' and its orientation-paired inverse to the model surface.  The marked segment is placed in the 'plane_index'ed plane of 'component_index'ed component of 'self', and its inverse is placed in the 'inverse_plane_index'ed plane of 'inverse_component_index'ed component of 'self'.'''
		marked_seg=marked_segment(holonomy)
		marked_seg.inverse=marked_seg.reverse()
		marked_seg.inverse.inverse=marked_seg
		self.components[component_index].add_marked_segment(marked_seg,plane_index)
		self.components[inverse_component_index].add_marked_segment(marked_seg.inverse,inverse_plane_index)

	def marked_segments_by_component_and_plane(self):
		'''List of all marked segments in the model surface, decomposed into lists corresponding to componenents and copies of the plane.'''
		return [[plane.marked_segments for plane in O_d.complex_planes] for O_d in self.components]

	def marked_segments_all(self):
		'''List of all marked segments in the model surface, not distinguishing between components or copies of plane.'''
		return [ms for O_d in self.components for plane in O_d.complex_planes for ms in plane.marked_segments]

	# def convex_body_pts_all(self):
	# 	'''List of all points contributing to convex bodies in the model surface, not distinguishing between components or copies of plane.'''
	# 	return [cbp for comp in self.convex_body_pts for plane in comp for cbp in plane]		

	def marked_segments_orientation_paired_representatives(self):
		'''List of representatives of orientation-paired marked segments in the model surface, not distinguishing between components or copies of plane.'''
		ms_all=self.marked_segments_all()
		ms_reps=[]
		for ms in ms_all:
			if ms not in ms_reps and ms.inverse not in ms_reps:
				ms_reps.append(ms)
		return ms_reps

	# def convex_body_pts_orientation_paired_representatives(self):
	# 	'''List of representatives of points contributing to convex bodies in the model surface, not distinguishing between components or copies of plane.'''
	# 	cbp_all=self.convex_body_pts_all()
	# 	cbp_reps=[]
	# 	for cbp in cbp_all:
	# 		if cbp not in cbp_reps and cbp.inverse not in cbp_reps:
	# 			cbp_reps.append(cbp)
	# 	return cbp_reps	

	def holonomies_by_component_and_plane(self):
		'''List of holonomies of marked segments in the model surface, decomposed into lists corresponding to componenents and copies of the plane.'''
		return [[[ms.holonomy for ms in plane.marked_segments] for plane in O_d.complex_planes] for O_d in self.components]

	def angles_of_holonomies_by_component_and_plane(self):
		'''List of angles of holonomies of marked segments in the model surface, decomposed into lists corresponding to componenents and copies of the plane.'''
		return [[[angle(ms.holonomy) for ms in plane.marked_segments] for plane in O_d.complex_planes] for O_d in self.components]		

	def plot(self):
		for comp in self.components:
			comp.plot()

	def update_convex_bodies_pts(self,plot=False):
		'''Find the marked segments contributing to edges of the convex bodies (set of points nearer the origin of the corresponding component than to any marked segment) of self.  If plot==True, a picture of the convex body is returned.'''
		self.convex_bodies_pts=[]
		for comp in self.components:
			self.convex_bodies_pts.append(comp.update_convex_body_pts())
		if plot==True:
			show(self.convex_bodies_plots())

	def convex_bodies_plots(self,update=False):
		'''Return plots of the convex bodies of components of self. If update==True, this runs comp.update_convex_body_pts() for each component comp of self first.  Otherwise we use data from the last time this was run.'''
		for comp in self.components:
			comp.convex_body_plot(update)

	def combine_marked_segments(self,model_surface):
		'''Adds all marked segments (with orientation-pairing data) from model_surface to self.  Raises error if orientation-pairing data is contradictory, or if two distinct marked segments lie on the same ray emanating from the origin.'''
		angles_self=self.angles_of_holonomies_by_component_and_plane()
		holonomies_self=self.holonomies_by_component_and_plane()
		for marked_seg in model_surface.marked_segments_orientation_paired_representatives():
			# Indices of component and plane to which marked_seg (and its inverse) belongs in model_surface
			component_index=marked_seg.belongs_to_component.index_in(model_surface)
			plane_index=marked_seg.belongs_to_plane.index_in(marked_seg.belongs_to_component)
			inverse_component_index=marked_seg.inverse.belongs_to_component.index_in(model_surface)
			inverse_plane_index=marked_seg.inverse.belongs_to_plane.index_in(marked_seg.inverse.belongs_to_component)			
			# Check if there are any existing marked segments of the same angle in self with the same component and plane indices
			if marked_seg.angle in angles_self[component_index][plane_index]:
				# If so, check if there is a marked segment with the same holonomy vector in the same component and plane
				if marked_seg.holonomy in holonomies_self[component_index][plane_index]:
					# Find the marked segment already existing in self that lies at the same point as marked_seg
					ms_index=self.holonomies_by_component_and_plane()[component_index][plane_index].index(marked_seg.holonomy)
					ms_self=self.marked_segments_by_component_and_plane()[component_index][plane_index][ms_index]
					# Check if orientation-pairs agree
					if ms_self.inverse.belongs_to_component.index_in(self)!=inverse_component_index or ms_self.inverse.belongs_to_plane.index_in(ms_self.inverse.belongs_to_component)!=inverse_plane_index:
						raise MarkedSegmentsError('The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')#(self,marked_seg,'The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')
					# else:
						# We need not do anything since marked_seg and its inverse are already present in self
				else:
					raise MarkedSegmentsError('The marked segment to be added lies on a ray from the origin passing through an existing marked segment, but at a different distance.')#(self,marked_seg,'The marked segment to be added lies on a ray from the origin passing through an existing marked segment, but at a different distance.')
			else:
				# If not, we need to make sure the angle of the inverse of marked segment is not already present in its corresponding component and plane: if it is, we obtain one of the previous two errors for the inverse
				if marked_seg.inverse.angle in angles_self[inverse_component_index][inverse_plane_index]:
					if marked_seg.inverse.holonomy in holonomies_self[inverse_component_index][inverse_plane_index]:
						raise MarkedSegmentsError('The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')#(self,marked_seg.inverse,'The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')
					else:
						raise MarkedSegmentsError('The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')#(self,marked_seg.inverse,'The marked segment to be added agrees with an existing marked segment, but their orientation-paired inverses do not agree.')
				else:
					# marked_seg and its inverse are new; we add them to self
					self.add_orientation_paired_marked_segments(marked_seg.holonomy,component_index,plane_index,inverse_component_index,inverse_plane_index)

	##### NEED TO THINK ABOUT WHETHER OR NOT THIS WORKS FOR `COMPLICATED' VORONOI DECOMPOSITIONS, E.G. WHEN THERE ARE TWO OR MORE `NEARBY' REGULAR VERTICES OF THE DECOMPOSITION #####
	def translation_surface(self,update_convex_bodies_pts=False):
		'''Construct the translation surface determined by the points contributing to edges of convex bodies, if possible.  The resulting translation surface is built in the flatsurf package.'''
		if update_convex_bodies_pts==True or self.convex_bodies_pts==None:
			self.update_convex_bodies_pts()
		# from flatsurf import *

		all_convex_bodies_pts=[mp for comp in self.convex_bodies_pts for eq_cls in comp for mp in eq_cls]
		convex_bodies_pts_paired_reps=[]
		for mp in all_convex_bodies_pts:
			if mp not in convex_bodies_pts_paired_reps and mp.inverse not in convex_bodies_pts_paired_reps:
				convex_bodies_pts_paired_reps.append(mp)
		convex_bodies_pts_by_component=[[] for comp in O.components]
		for mp in all_convex_bodies_pts:
			convex_bodies_pts_by_component[mp.component_index()].append(mp)
		# Sort convex_bodies_pts_by_component by cumulative angle
		for comp in convex_bodies_pts_by_component:
			comp.sort(key=lambda x:x.angle_cumulative())

		# # Will contain polygons and data for identifying edges
		# polygons_and_gluings=[]
		# while len(all_convex_bodies_pts)>0:
		# 	# List of two lists: first contains vertices of a polygon, and the second contains marked pts corresponding to traversing edges of the polygon counter-clockwise
		# 	vertices_and_marked_pts=[[(0,0)],[]]
		# 	mp=all_convex_bodies_pts[0]
		# 	all_convex_bodies_pts.remove(mp)
		# 	new_vertex=mp.holonomy
		# 	while new_vertex!=vector((0,0)):
		# 		vertices_and_marked_pts[0].append(new_vertex)
		# 		vertices_and_marked_pts[1].append(mp)
		# 		# Find next mp; this will be the marked point adjacent, clockwise, to mp.inverse
		# 		mp_inverse_index=convex_bodies_pts_by_component[mp.inverse.component_index()].index(mp.inverse)
		# 		mp=convex_bodies_pts_by_component[mp.inverse.component_index()][mp_inverse_index-1]
		# 		all_convex_bodies_pts.remove(mp)
		# 		new_vertex+=mp.holonomy
		# 	vertices_and_marked_pts[1].append(mp)
		# 	polygons_and_gluings.append([polygons(vertices=vertices_and_marked_pts[0],ring=QQbar),vertices_and_marked_pts[1]])

		# Will contain polygons and data for identifying edges
		polygons_and_gluings=[]
		while len(all_convex_bodies_pts)>0:
			# List of two lists: first contains vertices of a polygon, and the second contains marked pts corresponding to traversing edges of the polygon counter-clockwise
			# vertices_and_marked_pts=[[(0,0)],[]]
			vertices_and_marked_pts=[[],[]]
			# Choose an initial marked point in the 'gluing orbit'			
			mp_init=all_convex_bodies_pts[0]
			mp=mp_init
			first=True
			all_convex_bodies_pts.remove(mp)
			new_vertex=mp.holonomy
			while mp!=mp_init or first==True:
				first=False
				vertices_and_marked_pts[0].append(new_vertex)
				vertices_and_marked_pts[1].append(mp)
				# Find next mp; this will be the marked point adjacent, clockwise, to mp.inverse
				mp_inverse_index=convex_bodies_pts_by_component[mp.inverse.component_index()].index(mp.inverse)
				mp=convex_bodies_pts_by_component[mp.inverse.component_index()][mp_inverse_index-1]
				if mp in all_convex_bodies_pts:
					all_convex_bodies_pts.remove(mp)
				new_vertex+=mp.holonomy
			vertices_and_marked_pts[1].append(mp)
			polygons_and_gluings.append([polygons(vertices=vertices_and_marked_pts[0],ring=QQbar),vertices_and_marked_pts[1]])


					# while len(VS_all_copy)>0:
					# # Will count the interior angle of a vertex of convex bodies after identifiying edges
					# orbit_angle=0
					# # Choose an initial marked point in the 'gluing orbit'
					# mp_init=VS_all_copy[0]
					# mp0=mp_init
					# first=True
					# while mp0!=mp_init or first==True:
					# 	first=False
					# 	VS_all_copy.remove(mp0)
					# 	mp0_comp_index=mp0.belongs_to_component.index_in(mp0.belongs_to_component.component_of)
					# 	mp0_index_in_comp=VS_by_component[mp0_comp_index].index(mp0)
					# 	# Subsequent, counter-clockwise marked point from mp0
					# 	mp1=VS_by_component[mp0_comp_index][(mp0_index_in_comp+1)%len(VS_by_component[mp0_comp_index])]
					# 	# Interior angle of the perpendicular bisectors of mp0, mp1
					# 	interior_ang=angle((mp0.holonomy[1],-mp0.holonomy[0]))-angle((-mp1.holonomy[1],mp1.holonomy[0]))
					# 	if interior_ang<0:
					# 		interior_ang+=1
					# 	orbit_angle+=interior_ang
					# 	mp0=mp1.inverse
					# if orbit_angle!=1:
					# 	skip=True
					# 	break



		# Construct the resulting translation surface
		surface=Surface_list(base_ring=polygons_and_gluings[0][0].base_ring())
		for poly in polygons_and_gluings:
			surface.add_polygon(poly[0])
		for mp in convex_bodies_pts_paired_reps:
			for j in range(surface.num_polygons()):
				if mp in polygons_and_gluings[j][1]:
					mp_poly_index=j
					mp_edge_index=polygons_and_gluings[j][1].index(mp)
				if mp.inverse in polygons_and_gluings[j][1]:
					mp_inverse_poly_index=j
					mp_inverse_edge_index=polygons_and_gluings[j][1].index(mp.inverse)
			surface.change_edge_gluing(mp_poly_index,mp_edge_index,mp_inverse_poly_index,mp_inverse_edge_index)

		return TranslationSurface(surface)



		all_convex_bodies_pts=[mp for comp in self.convex_bodies_pts for eq_cls in comp for mp in eq_cls]
		convex_bodies_pts_by_component=[[] for comp in O.components]
		for mp in all_convex_bodies_pts:
			convex_bodies_pts_by_component[mp.component_index()].append(mp)
		# Sort convex_bodies_pts_by_component by cumulative angle
		for comp in convex_bodies_pts_by_component:
			comp.sort(key=lambda x:x.angle_cumulative())

		# List of triangles determined by each marked point from all_convex_bodies_pts, together with said marked point
		triangles_and_marked_pts=[[polygons(vertices=[(0,0),mp.convex_body_cw_vertex,mp.convex_body_ccw_vertex],ring=QQbar),mp] for mp in all_convex_bodies_pts]
		S=Surface_list(base_ring=triangles_and_marked_pts[0][0].base_ring())
		for tr_and_mp in triangles_and_marked_pts:
			S.add_polygon(tr_and_mp[0])
		for mp_index in range(len(all_convex_bodies_pts)):
			mp=all_convex_bodies_pts[mp_index]
			mp_inverse_index=all_convex_bodies_pts.index(mp.inverse)
			# Glue edges of triangles corresponding to mp and its orientation-paired inverse
			S.change_edge_gluing(mp_index,1,mp_inverse_index,1)
			mp_index_in_comp=convex_bodies_pts_by_component[mp.component_index()].index(mp)
			# Find marked point adjacent, clockwise, to mp
			mp_cw=convex_bodies_pts_by_component[mp.component_index()][mp_index_in_comp-1]
			mp_cw_index=all_convex_bodies_pts.index(mp_cw)
			# Glue edges of triangles correpsonding to mp and mp_cw
			S.change_edge_gluing(mp_index,0,mp_cw_index,2)

		return TranslationSurface(S)





class ModelSurfaceComponent:
	'''Infinite translation surface O_d=(C,(z^d)dz), thought of as d+1 copies of the complex plane C cut along the positive real axis and glued together cyclically.'''
	def __init__(self,d,model_surface=None):
		assert d in NN
		self.order=d
		self.num_planes=d+1
		self.component_of=model_surface
		self.complex_planes=[complex_plane(self) for j in range(self.order+1)]
		self.Trans=CyclicPermutationGroup(self.num_planes)
		# We keep track of which marked segments do not contribute to the convex body (points nearer the origin than any marked segment, constructed as intersection of (generalizations of) 'half-planes')
		self.not_convex_body_pts=None
		# Also keep track of which marked segments do contribute to convex body.  As more marked segments are considered, some points here may move to self.not_convex_body_pts
		self.convex_body_pts=None
	
	def __repr__(self):
		return 'Model surface O_{}=(C,(z^{})dz) as {} copies of the complex plane C'.format(self.order,self.order,self.num_planes)

	def add_marked_segment(self,marked_segment,plane_index):
		'''Adds a marked segment to 'plane_index'ed copy of plane.'''
		self.complex_planes[plane_index].add_marked_segment(marked_segment)

	def marked_segments_by_plane(self):
		'''List of all marked segments in the surface, decomposed into lists corresponding to copies of the plane.'''
		return [plane.marked_segments for plane in self.complex_planes]

	def marked_segments_all(self):
		'''List of all marked segments in the surface, not distinguishing between copies of plane.'''
		return [ms for plane in self.complex_planes for ms in plane.marked_segments]

	def holonomies_by_plane(self):
		'''List of holonomies of marked segments in the surface, decomposed into lists corresponding to copies of the plane.'''
		return [[ms.holonomy for ms in plane.marked_segments] for plane in self.complex_planes]

	def sort_marked_segments(self):
		'''Sort marked segments in each plane by counter-clockwise angle from positive real axis.'''
		for plane in self.complex_planes:
			plane.sort_marked_segments()

	# def sort_marked_segments_by_cumulative_angle(self,marked_segments):

	def index_in(self,O):
		'''Current index of self in the list O.components, where O is the surface of which self is a component.'''
		return O.components.index(self)

	def plot(self):#,diameter=None):
		holonomies=self.holonomies_by_plane()
		for plane in holonomies:
			show(point(plane,aspect_ratio=1))
			# if diameter==None:
			# 	diam=max([hol[0] for hol in plane]+[hol[1] for hol in plane])+1
			# else:
			# 	diam=diameter
			# show(point(plane),xmin=diam,xmax=diam,ymin=diam,ymax=diam)
		# for plane in self.complex_planes:
		# 	picture=point(holonomies)
		# 	for ms in plane.marked_segments:
		# 		picture+=text('{}'.format(ms.angle_to_inverse.n(digits=1)),ms.holonomy,horizontal_alignment='left')
		# 	show(picture,xmin=x_min,xmax=x_max,ymin=y_min,ymax=y_max)

	def update_convex_body_pts(self,plot=False):
		'''Find the marked segments contributing to edges of the convex body (set of points nearer the origin than to any marked segment) of self.  If plot==True, a picture of the convex body is returned.'''
		if self.not_convex_body_pts==None:
			self.not_convex_body_pts=[]
		self.convex_body_pts=[]		
		possible_convex_body_pts=[]
		for marked_seg in self.marked_segments_all():
			if marked_seg not in self.not_convex_body_pts:
				possible_convex_body_pts.append(marked_seg)
		# We partition the set of possible_convex_body_pts into maximal sets of marked segments whose minimal subsequent angles are less than pi
		# Let x, y be two marked segments in possible_convex_body_pts.  Define a relation x~y iff the (cumulative) angles of x, y differ by less than pi.  We extend this to an equivalance relation by forcing transitivity and find the resulting equivalence classes.
		equiv_classes=[]
		# Sort possible_convex_body_pts by cumulative angle
		possible_convex_body_pts.sort(key=lambda x:x.angle_cumulative())
		# We construct the equivalence classes one at a time and add them to equiv_classes
		if len(possible_convex_body_pts)>0:
			eq_cls=[possible_convex_body_pts[0]]
			for j in range(1,len(possible_convex_body_pts)):
				# If the cumulative angle of the next marked segment in possible_convex_body_pts is within pi of the previous one, we add them to the same equivalence class
				if possible_convex_body_pts[j].angle_cumulative()-eq_cls[-1].angle_cumulative()<1/2:
					eq_cls.append(possible_convex_body_pts[j])
				# Otherwise we begin a new equivalence class
				else:
					equiv_classes.append(eq_cls)
					eq_cls=[possible_convex_body_pts[j]]
			equiv_classes.append(eq_cls)
			# Need to check if the first and last equivalence classes (if there are multiple) should really be one.  This is the case if the first and last marked segments in possible_convex_body_pts are within pi
			if len(equiv_classes)>1 and possible_convex_body_pts[0].angle_cumulative()-possible_convex_body_pts[-1].angle_cumulative()+self.num_planes<1/2:
				equiv_classes[0]=equiv_classes.pop()+equiv_classes[0]

			# We construct edges of the convex body for each equivalence class
			for j in range(len(equiv_classes)):
				eq_cls=equiv_classes[j]
				# If there is only one point in the equivalence class, then it contributes to the convex body and we can skip the remaining steps
				if len(eq_cls)==1:
					self.convex_body_pts.append(eq_cls)
					continue
				eq_cls_convex_body_pts=[]
				convex_body_compact=False
				# If there is only one equivalence class and its first and last elements are within pi (i.e. the convex body is compact), then the shortest representatives are guaranteed to contribute to the sides of the convex body.  We find one of these and mark its index in its equivalence class as k0_short.
				# (Note: modulo, '%', can act funny and give negatives, so 'diff' here is a workaround)
				diff=(eq_cls[0].angle_cumulative()-eq_cls[-1].angle_cumulative()).n()%self.num_planes
				if diff<0:
					diff=diff+self.num_planes
				if len(equiv_classes)==1 and diff<1/2:
					convex_body_compact=True
					lengths=[ms.length_squared() for ms in eq_cls]
					min_length=min(lengths)
					k0_short=lengths.index(min_length)
					k0=k0_short				
				# If there are multiple equivalence classes, then the clockwise-most and counter-clockwise-most marked segments in the equivalence class are guaranteed to contribute to the sides of the convex body.  We begin with the clockwise-most
				else:
					k0=0
				keep_finding=True
				while keep_finding==True:
					ms0=eq_cls[k0]
					eq_cls_convex_body_pts.append(ms0)
					# Find all elements of eq_cls within pi counter-clockwise of ms0
					within_pi_of_ms0=[]
					l=1
					ms1=eq_cls[(k0+l)%len(eq_cls)]
					# (Another workaround for '%'; see note above.)
					diff=(ms1.angle_cumulative()-ms0.angle_cumulative()).n()%self.num_planes
					if diff<0:
						diff+=self.num_planes
					while diff<1/2 and ms0!=ms1:
						within_pi_of_ms0.append(ms1)
						l+=1
						ms1=eq_cls[(k0+l)%len(eq_cls)]
						diff=(ms1.angle_cumulative()-ms0.angle_cumulative()).n()%self.num_planes
						if diff<0:
							diff+=self.num_planes					
					# Find perpendicular bisectors of each line segment from the origin to each element of within_pi_of_ms0, and the intersections of these bisectors with the perp. bisector corresponding to the ms0.  The marked segment whose intersection is nearest the orgin contributes to a side of the convex body. (Note: There may be multiple.  We deal with this below.)
					perp_bis_intersections=[]
					x0=ms0.holonomy[0]
					y0=ms0.holonomy[1]
					for ms1 in within_pi_of_ms0:
						x1=ms1.holonomy[0]
						y1=ms1.holonomy[1]
						t=(x1*(x1-x0)+y1*(y1-y0))/(2*(x0*y1-y0*x1))
						perp_bis_intersections.append((-y0*t+x0/2,x0*t+y0/2))
					# Find the intersection nearest the origin
					perp_bis_intersections_lengths=[z[0]^2+z[1]^2 for z in perp_bis_intersections]
					shortest_int_len=min(perp_bis_intersections_lengths)
					# There could be multiple marked segment giving this shortest intersection.  The one whose angle is furthest from that of ms0 will 'cut' into the convex body most and thus contributes to the convex body
					shortest_int_len_indices=[]
					for k in range(len(perp_bis_intersections_lengths)):
						if perp_bis_intersections_lengths[k]==shortest_int_len:
							shortest_int_len_indices.append(k)
					shortest_int_marked_segments=[within_pi_of_ms0[k] for k in shortest_int_len_indices]
					# Sort by angle from ms0
					shortest_int_marked_segments_ang_from_ms0=[(x.angle_cumulative()-ms0.angle_cumulative()).n()%self.num_planes for x in shortest_int_marked_segments]
					# (Another workaround for '%'; see note above.)
					for k in range(len(shortest_int_marked_segments_ang_from_ms0)):
						if shortest_int_marked_segments_ang_from_ms0[k]<0:
							shortest_int_marked_segments_ang_from_ms0[k]+=self.num_planes
					k1=shortest_int_marked_segments_ang_from_ms0.index(max(shortest_int_marked_segments_ang_from_ms0))
					# shortest_int_marked_segments.sort(key=lambda x:(x.angle_cumulative()-ms0.angle_cumulative())%self.num_planes)
					ms1=shortest_int_marked_segments[k1]
					# Update the vertices of the edge of the convex body determined by ms0 (its counter-clockwise vertex) and ms1 (its clock-wise vertex)
					x1=ms1.holonomy[0]
					y1=ms1.holonomy[1]
					t=(x1*(x1-x0)+y1*(y1-y0))/(2*(x0*y1-y0*x1))
					ms0.convex_body_ccw_vertex=(-y0*t+x0/2,x0*t+y0/2)
					ms1.convex_body_cw_vertex=(-y0*t+x0/2,x0*t+y0/2)
					# if ms0.convex_body_edge==None:
					# 	ms0.convex_body_edge=[infinity,shortest_int]
					# else:
					# 	ms0.convex_body_edge[1]=shortest_int
					# if ms1.convex_body_edge==None:
					# 	ms1.convex_body_edge=[shortest_int,infinity]
					# else:
					# 	ms1.convex_body_edge[0]=shortest_int
					k0=eq_cls.index(ms1)
					# If convex body is compact and we've returned to the original marked segment, we're done
					if convex_body_compact==True and k0==k0_short:
						keep_finding=False
					# Else, if we've reached the final (counter-clockwise-most) segment in eq_cls, we add this to eq_cls_convex_body_pts and are done
					elif convex_body_compact==False and k0==len(eq_cls)-1:
						eq_cls_convex_body_pts.append(ms1)
						keep_finding=False
				# Add to self.not_convex_body_pts those marked segments not in eq_cls_convex_body_pts
				for ms in possible_convex_body_pts:
					if ms not in eq_cls_convex_body_pts:
						self.not_convex_body_pts.append(ms)
				self.convex_body_pts.append(eq_cls_convex_body_pts)
		if plot==True:
			show(self.convex_body_plot())
		return self.convex_body_pts





	def convex_body_plot(self,update=False):
		'''Return a plot of the convex body of self using marked segments in (equivalence classes in) self.convex_body_pts. If update==True, this runs self.update_convex_body_pts() first.  Otherwise we use data from the last time this was run.'''
		if update==True:
			self.update_convex_body_pts()
		pictures_in_each_plane=[[] for d in range(self.num_planes)]
		# To be used to place a label in a reasonable spot in each copy of the plane
		max_image_len_squared=0
		if self.convex_body_pts!=None:
			for eq_cls in self.convex_body_pts:
				# If the edge corresponding to marked_seg is of infinite length, then we artificially choose a length for which to plot (see also below)
				if len(eq_cls)<=2:
					len_of_infinite_edge_to_plot=2
				else:				
					edge_lengths_squared=[]
					for ms in eq_cls:
						if ms.convex_body_cw_vertex!=None and ms.convex_body_ccw_vertex!=None:
							edge_lengths_squared.append((ms.convex_body_cw_vertex[0]-ms.convex_body_ccw_vertex[0])^2+(ms.convex_body_cw_vertex[1]-ms.convex_body_ccw_vertex[1])^2)			
					len_of_infinite_edge_to_plot=sqrt(max(edge_lengths_squared))

				# List of lines to plot in each copy of the plane
				lines_to_plot=[[] for d in range(self.num_planes)]
				new_line=True
				for j in range(len(eq_cls)):
					# A marked segment in eq_cls and the index of the plane it belongs to
					marked_seg=eq_cls[j]
					marked_seg_plane_index=marked_seg.belongs_to_plane.index_in(self)
					# Add each marked segment and data (index of component: copy of plane) pointing to its inverse to plot
					marked_seg_inverse_component_index=marked_seg.inverse.belongs_to_component.index_in(marked_seg.inverse.belongs_to_component.component_of)
					marked_seg_inverse_plane_index=marked_seg.inverse.belongs_to_plane.index_in(marked_seg.inverse.belongs_to_component)
					pictures_in_each_plane[marked_seg_plane_index].append(point(marked_seg.holonomy,aspect_ratio=1))
					pictures_in_each_plane[marked_seg_plane_index].append(text('  {}: c{}'.format(marked_seg_inverse_component_index,marked_seg_inverse_plane_index),marked_seg.holonomy,horizontal_alignment='left',aspect_ratio=1))

					# Coordinates of the holonomy of marked_seg
					x=marked_seg.holonomy[0]
					y=marked_seg.holonomy[1]

					# Vertices of edge of the convex body corresponding to marked_seg
					cw_vertex=marked_seg.convex_body_cw_vertex
					ccw_vertex=marked_seg.convex_body_ccw_vertex

					# If the edge corresponding to marked_seg is of infinite length, then we artificially choose a length for which to plot
					if cw_vertex==None and ccw_vertex!=None:
						cw_vertex=(y*len_of_infinite_edge_to_plot/sqrt(x^2+y^2)+ccw_vertex[0],-x*len_of_infinite_edge_to_plot/sqrt(x^2+y^2)+ccw_vertex[1])
					elif cw_vertex!=None and ccw_vertex==None:
						ccw_vertex=(-y*len_of_infinite_edge_to_plot/sqrt(x^2+y^2)+cw_vertex[0],x*len_of_infinite_edge_to_plot/sqrt(x^2+y^2)+cw_vertex[1])
					elif cw_vertex==None and ccw_vertex==None:
						cw_vertex=(y*len_of_infinite_edge_to_plot/(2*sqrt(x^2+y^2))+x/2,-x*len_of_infinite_edge_to_plot/(2*sqrt(x^2+y^2))+y/2)
						ccw_vertex=(-y*len_of_infinite_edge_to_plot/(2*sqrt(x^2+y^2))+x/2,x*len_of_infinite_edge_to_plot/(2*sqrt(x^2+y^2))+y/2)				

					# If marked_seg belongs to the (closure of the) second or third quadrant, then the corresponding edge must belong to the same plane as marked_seg.  This is also the case if marked_seg belongs to the first (or fourth) quadrant and the endpoints of the edge are non-negatative (resp. negative)
					if x<=0 or (y>=0 and cw_vertex[1]>=0) or (y<0 and ccw_vertex[1]<0):
						if marked_seg.convex_body_cw_vertex==None:
							lines_to_plot[marked_seg_plane_index].append(arrow(ccw_vertex,cw_vertex,aspect_ratio=1))
						if marked_seg.convex_body_ccw_vertex==None:
							lines_to_plot[marked_seg_plane_index].append(arrow(cw_vertex,ccw_vertex,aspect_ratio=1))
						else:
							lines_to_plot[marked_seg_plane_index].append(arrow(cw_vertex,ccw_vertex,arrowsize=.1,aspect_ratio=1))
					# If marked_seg belongs to the first quadrant, the clockwise vertex belongs to the lower half-plane and the counter-clockwise vertex belongs to the upper half-plane, then we split the edge where it meets the (positive) real axis and place the two segments in their corresponding copies of the plane
					elif y>=0 and cw_vertex[1]<0 and ccw_vertex[1]>=0:
						# Intersection of perpendicular bisector with the (positive) real-axis
						int_with_real_axis=(y^2/(2*x)+x/2,0)
						if marked_seg.convex_body_cw_vertex==None:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(int_with_real_axis,cw_vertex,aspect_ratio=1))
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,ccw_vertex,arrowsize=.1,aspect_ratio=1))
						if marked_seg.convex_body_ccw_vertex==None:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(int_with_real_axis,cw_vertex,arrowsize=.1,aspect_ratio=1))
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,ccw_vertex,aspect_ratio=1))
						else:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(int_with_real_axis,cw_vertex,arrowsize=.1,aspect_ratio=1))
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,ccw_vertex,arrowsize=.1,aspect_ratio=1))
					# If marked_seg belongs to the fourth quadrant, the clockwise vertex belongs to the lower half-plane and the counter-clockwise vertex belongs to the upper half-plane, then we split the edge where it meets the (positive) real axis and place the two segments in their corresponding copies of the plane
					elif y<=0 and cw_vertex[1]<0 and ccw_vertex[1]>=0:
						# Intersection of perpendicular bisector with the (positive) real-axis
						int_with_real_axis=(y^2/(2*x)+x/2,0)
						if marked_seg.convex_body_cw_vertex==None:
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,cw_vertex,aspect_ratio=1))
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(int_with_real_axis,ccw_vertex,arrowsize=.1,aspect_ratio=1))
						if marked_seg.convex_body_ccw_vertex==None:
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,cw_vertex,arrowsize=.1,aspect_ratio=1))
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(int_with_real_axis,ccw_vertex,aspect_ratio=1))
						else:
							lines_to_plot[marked_seg_plane_index].append(arrow(int_with_real_axis,cw_vertex,arrowsize=.1,aspect_ratio=1))
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(int_with_real_axis,ccw_vertex,arrowsize=.1,aspect_ratio=1))
					# If marked_seg belongs to the first quadrant, but the edge belongs to the fourth quadrant, then we place the edge in the previous plane
					elif y>=0 and ccw_vertex[1]<=0:
						if marked_seg.convex_body_cw_vertex==None:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(ccw_vertex,cw_vertex,aspect_ratio=1))
						if marked_seg.convex_body_ccw_vertex==None:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(cw_vertex,ccw_vertex,aspect_ratio=1))
						else:
							lines_to_plot[(marked_seg_plane_index-1)%self.num_planes].append(arrow(cw_vertex,ccw_vertex,arrowsize=.1,aspect_ratio=1))
					# If marked_seg belongs to the fourth quadrant, but the edge belongs to the first quadrant, then we place the edge in the next plane
					elif y<=0 and cw_vertex[1]>=0:
						if marked_seg.convex_body_cw_vertex==None:
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(ccw_vertex,cw_vertex,aspect_ratio=1))
						if marked_seg.convex_body_ccw_vertex==None:
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(cw_vertex,ccw_vertex,aspect_ratio=1))
						else:
							lines_to_plot[(marked_seg_plane_index+1)%self.num_planes].append(arrow(cw_vertex,ccw_vertex,arrowsize=.1,aspect_ratio=1))

					# Update max_image_len_squared (used for placement of a label of each copy of plane) as needed
					if cw_vertex[0]^2+cw_vertex[1]^2>max_image_len_squared:
						max_image_len_squared=cw_vertex[0]^2+cw_vertex[1]^2
					if ccw_vertex[0]^2+ccw_vertex[1]^2>max_image_len_squared:
						max_image_len_squared=ccw_vertex[0]^2+ccw_vertex[1]^2
					if x^2+y^2>max_image_len_squared:
						max_image_len_squared=x^2+y^2						

				for d in range(self.num_planes):
					pictures_in_each_plane[d]+=lines_to_plot[d]
			for d in range(self.num_planes):
				pictures_in_each_plane[d].append(text('{}: c{}'.format(self.index_in(self.component_of),d),[1.2*sqrt(max_image_len_squared/2),1.2*sqrt(max_image_len_squared/2)],fontsize='large',color='black',aspect_ratio=1))
				show(sum(pictures_in_each_plane[d]))			


					


def angle(holonomy):
	'''Counter-clockwise angle of vector 'holonomy' from the positive real axis, modulo 2pi.'''
	x=holonomy[0]
	y=holonomy[1]
	if x>0:
		if y>=0:
			theta=arctan(y/x)/(2*pi)
		else:
			theta=arctan(y/x)/(2*pi)+1
	elif x<0:
			theta=arctan(y/x)/(2*pi)+1/2
	else:
		if y>0:
			theta=1/4
		else: theta=3/4
	return theta





class complex_plane:
	'''A copy of the complex plane, cut along the positive real axis, embedded in the infinite translation surface O_d.'''
	def __init__(self,O_d):
		self.embedded_in=O_d
		# List of marked segments in self
		self.marked_segments=[]

	def __repr__(self):
		return 'Copy of the complex plane, cut along the positive real axis, embedded in O_{}'.format(self.embedded_in.order)

	def add_marked_segment(self,marked_segment):
		'''Adds a marked segment to this copy of the plane.'''
		self.marked_segments.append(marked_segment)
		marked_segment.belongs_to_component=self.embedded_in
		marked_segment.belongs_to_plane=self


	def sort_marked_segments(self):
		'''Sorts all marked segments in this copy of the plane by counter-clockwise angle from positive real axis.'''
		self.marked_segments.sort(key=lambda x:x.angle)

	def index_in(self,O_d):
		'''Current index of self in the list O_d.complex_planes.'''
		return O_d.complex_planes.index(self)





class marked_segment:
	'''A marked segment defined by a list or tuple 'holonomy'.  Optional inputs are a surface 'O_d' and 'plane_index', which specify where to place the marked segment.'''
	def __init__(self,holonomy,O_d=None,plane_index=None):#,O_d,plane_index):
		assert len(holonomy)==2
		self.holonomy=vector(holonomy)
		self.angle=angle(self.holonomy)
		if O_d!=None and plane_index!=None:
			O_d.add_marked_segment(self,plane_index)
		# When self is assigned to a particular plane of a particular component of a model surface, we record that data as attributes
		self.belongs_to_component=None
		self.belongs_to_plane=None
		# # When we begin constructing convex bodies in each component (set of points nearer the origin than any marked segment), we record whether or not self determines a side of the convex body.  Note: self.contributes_to_convex_body may be True for some collection of marked segments, then False as more are added, but once it is False it is always False
		# self.contributes_to_convex_body=None
		# # If self.contributes_to_convex_body==True, self.convex_body_edge will record the clockwise and counterclockwise vertices of the edge of the convex body determined by self
		self.convex_body_ccw_vertex=None
		self.convex_body_cw_vertex=None

	def __repr__(self):
		return 'Marked segment with holonomy {}'.format(self.holonomy)

	def reverse(self):
		'''Create a marked segment with reversed holonomy vector.'''
		return marked_segment((-self.holonomy[0],-self.holonomy[1]))

	def length_squared(self):
		'''Return the length of the holonomy of self, squared.'''
		return self.holonomy[0]^2+self.holonomy[1]^2

	def angle_cumulative(self):
		'''Return the cumulative counter-clockwise angle from the positive real axis in the initial plane of the component in which self belongs to self.  Note, this depends on the ordering of planes in the component.'''
		return self.angle+self.belongs_to_plane.index_in(self.belongs_to_component)

	def component_index(self):
		return self.belongs_to_component.index_in(self.belongs_to_component.component_of)

	def plane_index(self):
		return self.belongs_to_plane.index_in(self.belongs_to_component)





def affine_diffeo(model_surface,A,tau=None):
	'''Apply the affine diffeomorphism tau(f_A) to 'model_surface', where f_A, tau are themselves affine diffeomorphisms satisfying der(f_A)=A, der(tau)=Id.  Here 'A' belongs to GL^+(2,R) and 'tau' is an element of self.Trans, if specified.  By default, tau=None, which we interpret as the identity.  (Note: the map f_A is ambigous as stated, but we make a canonical choice depending on the image of e_1 under A in the plane.)'''
	assert det(A)>0
	# First apply the affine diffeomorphism f_A
	image_O=ModelSurface(model_surface.orders_of_zeros)
	angle_of_A_e1=angle(A*vector([1,0]))
	for ms in model_surface.marked_segments_orientation_paired_representatives():
		# Image of ms under A
		A_ms=marked_segment(A*ms.holonomy)
		# Components of image_O in which A_ms and its inverse will be placed
		A_ms_component_index=ms.belongs_to_component.index_in(model_surface)
		A_ms_inverse_component_index=ms.inverse.belongs_to_component.index_in(model_surface)
		# If the angle of A_ms is greater than or equal to angle_of_A_e1, then we assign this marked segment to the same plane as ms; otherwise, it gets assigned to the next plane
		if A_ms.angle>=angle_of_A_e1:
			A_ms_plane_index=ms.belongs_to_plane.index_in(ms.belongs_to_component)
		else:
			A_ms_plane_index=(ms.belongs_to_plane.index_in(ms.belongs_to_component)+1)%ms.belongs_to_component.num_planes
		# Similarly for the inverse of A_ms
		if A_ms.reverse().angle>=angle_of_A_e1:
			A_ms_inverse_plane_index=ms.inverse.belongs_to_plane.index_in(ms.inverse.belongs_to_component)
		else:
			A_ms_inverse_plane_index=(ms.inverse.belongs_to_plane.index_in(ms.inverse.belongs_to_component)+1)%ms.inverse.belongs_to_component.num_planes
		image_O.add_orientation_paired_marked_segments(A_ms.holonomy,A_ms_component_index,A_ms_plane_index,A_ms_inverse_component_index,A_ms_inverse_plane_index)

	# Next apply tau
	if tau!=model_surface.Trans.one() and tau!=None:
		permutation_of_components=tau[0]
		permutations_of_planes=tau[1]
		# First we cyclically permute copies of the plane in each component ModelSurfaceComponent(d_i)
		for j in range(len(image_O.orders_of_zeros)):
			cyclic_perm=permutations_of_planes[j]
			image_O.components[j].complex_planes=cyclic_perm(image_O.components[j].complex_planes)
		# Then we permute components of same order
		place=0
		for j in range(len(image_O.cardinalities_of_distinct_zeros)):
			n=image_O.cardinalities_of_distinct_zeros[j]
			image_O.components[place:place+n]=permutation_of_components[j](image_O.components[place:place+n])
			place+=n
	return image_O





def affine_diffeo_inverse(model_surface,A,tau=None):
	'''Apply the inverse of the affine diffeomorphism tau(f_A) to 'model_surface', where f_A, tau are themselves affine diffeomorphisms satisfying der(f_A)=A, der(tau)=Id.  Here 'A' belongs to GL^+(2,R) and 'tau' is an element of self.Trans, if specified.  By default, tau=None, which we interpret as the identity.  (Note: the map f_A is ambigous as stated, but we make a canonical choice depending on the image of e_1 under A in the plane).'''
	assert det(A)>0
	if tau==None:
		tau=model_surface.Trans.one()
	# We want to find tau_prime for which (tau*f_A)^{-1}=tau_prime*f_{A^{-1}}.  
	# Note the left-hand side is (f_A)^{-1}*tau^{-1}, and these commute, so (tau*f_A)^{-1}=tau^{-1}*(f_A)^{-1}.  
	# One can show that sigma:=f_{A^{-1}}*f_A is the identity if e_1 is an eigenvector of A with positive eigenvalue,
	# and otherwise acts as a counter-clockwise rotation by 2pi on each component.
	# Setting tau_prime=tau^{-1}*sigma^{-1}, we have tau_prime*f_{A^{-1}}=tau^{-1}*sigma^{-1}*f_{A^{-1}}=tau^{-1}*(f_A)^{-1}=(f_A*tau)^{-1}, as desired.

	# If e_1 is an eigenvector of A with positive eigenvalue, set sigma equal to the identity
	if angle(A*vector((1,0)))==0:
		sigma_inverse=model_surface.Trans.one()
	else: 
		# The element of the cartesian product of cyclic permutations corresponding to each component, each entry of which corresponds to a clockwise rotation by 2pi
		two_pi_clockwise_perms=model_surface.Trans.cartesian_factors()[1](model_surface.Trans.cartesian_factors()[1].cartesian_factors()[j].gen() for j in range(len(model_surface.orders_of_zeros)))
		# The element of model_surface.Trans with the identity permutation on components of the same order, and acting on each individual component as a clockwise rotation by 2pi
		sigma_inverse=model_surface.Trans((model_surface.Trans.cartesian_factors()[0].one(),two_pi_clockwise_perms))
	tau_inverse=Trans_inverse(model_surface,tau)
	# Construct tau_prime in model_surface.Trans as tau_inverse*sigma_inverse
	tau_prime_permutations_of_planes=[]
	for k in range(len(model_surface.orders_of_zeros)):
		tau_prime_permutations_of_planes.append(tau_inverse[1][k]*sigma_inverse[1][k])
	tau_prime=model_surface.Trans([tau_inverse[0],tau_prime_permutations_of_planes])

	return affine_diffeo(model_surface,A.inverse(),tau_prime)





def Trans_inverse(model_surface,tau):
	'''Return tau^{-1}, where tau is an affine diffeomorphism satisfying der(tau)=Id in model_surface.Trans'''

	# tau is given as (permutation_of_components,permutations_of_planes); its inverse is given as (new_permutation_of_components,new_permutations_of_planes), where new_permutation_of_components consists of inverses of permutation_of_components (in the same order), and new_permutations_of_planes consists of the inverses of the cyclic permutations comprising permutations_of_planes, permuted by permutation_of_components
	permutation_of_components=tau[0]
	permutations_of_planes=tau[1]
	new_permutation_of_components=[]
	new_permutations_of_planes=[]
	place=0
	for j in range(len(model_surface.cardinalities_of_distinct_zeros)):
		n=model_surface.cardinalities_of_distinct_zeros[j]
		# Permute the cyclic permutations permutations_of_planes by permutations in permutations_of_components and add to new_permutations_of_planes (we take their inverses below)
		perms=permutation_of_components[j](permutations_of_planes[place:place+n])
		for perm in perms:
			new_permutations_of_planes.append(perm)
		place+=n
		# Add inverses of permutations in permutation_of_components to new_permutation_of_components
		new_permutation_of_components.append(permutation_of_components[j].inverse())
	# Replace cyclic permutations in new_permutations_of_planes by their inverses
	for d in range(len(new_permutations_of_planes)):
		new_permutations_of_planes[d]=new_permutations_of_planes[d].inverse()

	return model_surface.Trans([new_permutation_of_components,new_permutations_of_planes])





class Error(Exception):
	'''Base class for exceptions.'''
	pass

class MarkedSegmentsError(Error):
	'''Raised when attempting to add marked segments to a model surface that are not allowed.'''
	pass
	



########## LOOK INTO THIS FOR M_m, M_n, m=4, n=8 ##########
def normalized_short_staple_representatives(generators,return_plot=False):
	'''Given a list of generators (in SL(2,R)) of a lattice Fuchsian group Gamma with parabolic directions dense in S^1, return a finite list of vectors containing the finite set S^1-(union M*B_1(0)) where the union is over all M in Gamma.  The finite set in this union corresponds to directions of shortest representatives of Gamma-orbits of saddle connections of a translation surface with Veech group (containing) Gamma.'''
	# We represent S^1 as the interval [0,1) with 0~1, where theta in [0,1) corresponds to e^{2pi*i*theta} on S^1
	S1_minus_images_of_ball=RealSet.closed_open(0,1)
	matrices_applied=[matrix.identity(2)]
	while S1_minus_images_of_ball.cardinality()==infinity:
		for mat in generators:
			products=[mat*N for N in matrices_applied]+[mat.inverse()*N for N in matrices_applied]+[N*mat for N in matrices_applied]+[N*mat.inverse() for N in matrices_applied]
			for M in products:
				if M not in matrices_applied:
					matrices_applied.append(M)
					a=M[0][0]
					b=M[0][1]
					c=M[1][0]
					d=M[1][1]
					# S^1-M*B_1(0) is empty if and only if M acts as a rotation (iff M has Frobenius norm sqrt(2))
					if a^2+b^2+c^2+d^2!=2:
						# The vectors (x,y) in S^1 with M*(x,y) in S^1 are those satisfying 1=(ax+by)^2+(cx+dy)^2=Ax^2+-Bx*sqrt(1-x^2)+C, where A, B, C are given by:
						A=a^2-b^2+c^2-d^2
						B=2*(a*b+c*d)
						C=b^2+d^2
						# We solve for the different solutions (x,y), and their images M*(x,y) on S^1
						X=[sqrt((2*A-2*A*C+B^2+sqrt((2*A-2*A*C+B^2)^2-4*(A^2+B^2)*(C^2-2*C+1)))/(2*(A^2+B^2))),sqrt((2*A-2*A*C+B^2-sqrt((2*A-2*A*C+B^2)^2-4*(A^2+B^2)*(C^2-2*C+1)))/(2*(A^2+B^2))),-sqrt((2*A-2*A*C+B^2+sqrt((2*A-2*A*C+B^2)^2-4*(A^2+B^2)*(C^2-2*C+1)))/(2*(A^2+B^2))),-sqrt((2*A-2*A*C+B^2-sqrt((2*A-2*A*C+B^2)^2-4*(A^2+B^2)*(C^2-2*C+1)))/(2*(A^2+B^2)))]
						images=[]
						for x in X:
							y=sqrt(1-x^2)
							image_x_plus=a*x+b*y
							image_y_plus=c*x+d*y
							image_x_minus=a*x-b*y
							image_y_minus=c*x-d*y
							if (image_x_plus)^2+(image_y_plus)^2==1 and (image_x_plus,image_y_plus) not in images:
								images.append((image_x_plus,image_y_plus))
							elif (image_x_minus)^2+(image_y_minus)^2==1 and (image_x_minus,image_y_minus) not in images:
								images.append((image_x_minus,image_y_minus))
						# Sort images by counter-clockwise angle from horizontal
						images.sort(key=lambda x:angle(x))
						# Find differences between subsequent angles.  Since the ellipse M*B_1(0) has area 1, S^1 intersected with M*B_1(0) consists of the two open arcs between points in images of minimal subsequent angles
						angles=[angle(v) for v in images]
						differences=[]
						for j in range(len(angles)):
							diff=angles[j]-angles[j-1]
							if diff>0:
								differences.append(diff)
							else:
								differences.append(diff+1)
						min_diff=min(differences)
						intervals_to_take_away=[]
						for j in range(len(differences)):
							if differences[j]==min_diff:
								if angles[j]>angles[j-1]:
									intervals_to_take_away.append(RealSet(angles[j-1],angles[j]))
								else:
									intervals_to_take_away.append(RealSet.closed_open(0,angles[j]))
									intervals_to_take_away.append(RealSet(angles[j-1],1))

						for interval in intervals_to_take_away:
							S1_minus_images_of_ball=S1_minus_images_of_ball.difference(interval)
	thetas=[angle.lower() for angle in S1_minus_images_of_ball]
	if return_plot==True:
		show(ellipse_plot(matrices_applied))
	return [(cos(2*pi*theta),sin(2*pi*theta)) for theta in thetas]							





def ellipse_plot(matrices):
	'''Plot S^1 and its image under each matrix in list matrices.'''
	t=var('t')
	return sum([parametric_plot((M[0][0]*cos(t)+M[0][1]*sin(t),M[1][0]*cos(t)+M[1][1]*sin(t)),(t,0,2*pi)) for M in matrices]+[parametric_plot((cos(t),sin(t)),(t,0,2*pi),color='orange')])





# Will print progress of various for loops
import sys
def print_percent_complete(msg, i, n_i):
    i, n_i = int(i)+1, int(n_i)
    sys.stdout.write('\r')
    sys.stdout.write("{} {}/{}" .format(msg, i, n_i)) #, (100/(n_i)*i).n()))
    sys.stdout.flush()
    if i/n_i == 1:
        print()



def possible_voronoi_staples(model_surfaces):
	'''Given a list of model surfaces, return lists of marked segments from these which could possibly be (normalized) Voronoi staples for a translation surface.'''
	





def construct(generators,stratum,iteration_limit=3,short_saddle_connection_directions_representatives=None):
	'''Constructs translation surfaces in stratum with Veech group generated by generators.  The input generators is a finite list of matrices in SL(2,R), and stratum is a list of the orders of zeros.  The optional argument short_saddle_connection_directions_representatives can be given (as a list of vectors in S^1) to bypass the construction of a finite list (determined by generators) containing the short saddle connection directions of the translation surfaces in question.'''

	##### Find directions of short saddle connections #####
	if short_saddle_connection_directions_representatives==None:
		short_saddle_connection_directions=normalized_short_staple_representatives(generators)
		# If v belongs to short_saddle_connection_directions, then so does -v.  We need only consider one of these
		short_saddle_connection_directions_representatives=[]
		for v in short_saddle_connection_directions:
			if v not in short_saddle_connection_directions_representatives and (-v[0],-v[1]) not in short_saddle_connection_directions_representatives:
				short_saddle_connection_directions_representatives.append(v)
	else:
		for v in short_saddle_connection_directions_representatives:
			assert v[0]^2+v[1]^2==1

	##### Create model surfaces with marked segments placed in various copies of the plane, and simulations from the model surfaces with different choices of Trans(O).  These Trans(O) elements will be composed with affine diffeomorphisms corresponding to elements of generators #####
	O_blank=ModelSurface(stratum)
	cart_prod_of_Trans=cartesian_product([O_blank.Trans for A in generators])
	simulations=[]

	# We minimize the number of model surfaces needed by normalizing with Trans(O); we 'un-normalize' below by applying each element of Trans(O) to each of the simulations
	for v in short_saddle_connection_directions_representatives:
		v_component_index=-O_blank.cardinalities_of_distinct_zeros[0]
		# Create model surfaces where v belongs to the (first plane of the) first component of any given order; we can normalize by Trans(O) so that this is always the case
		for n in O_blank.cardinalities_of_distinct_zeros:
			v_component_index+=n
			# Consider where to place the orientation-paired inverse of v
			v_inverse_component_index_holder=-O_blank.cardinalities_of_distinct_zeros[0]
			for m in O_blank.cardinalities_of_distinct_zeros:				
				v_inverse_component_index_holder+=m
				# Here we place the inverse of v in a component of the same order as that in which v is placed
				if v_inverse_component_index_holder==v_component_index:
					# For each plane of the component in which v is placed, we create a model surface with the inverse of v in that plane
					for v_inverse_plane_index in range(O_blank.orders_of_zeros[v_component_index]+1):
						for taus in cart_prod_of_Trans:
							O=ModelSurface(stratum)
							O.add_orientation_paired_marked_segments(v,v_component_index,0,v_component_index,v_inverse_plane_index)
							# model_surfaces.append(O)
							simulations.append([O,taus])
					# If there are multiple components of the same order as that in which v is place, we also create a model surface with the inverse of v placed in one of these other components (we can normalize by Trans(O) so that the inverse of v is in the first plane of the second component of this order)
					if n>1:
						for taus in cart_prod_of_Trans:
							O=ModelSurface(stratum)
							O.add_orientation_paired_marked_segments(v,v_component_index,0,v_component_index+1,0)
							# model_surfaces.append(O)
							simulations.append([O,taus])
				# Here we place the inverse of v in a component of a different order as that in which v is place (we can normalize by Trans(O) so that the inverse of v is placed in the first plane of the first component of this different order)
				else:
					for taus in cart_prod_of_Trans:
						O=ModelSurface(stratum)
						O.add_orientation_paired_marked_segments(v,v_component_index,0,v_inverse_component_index_holder,0)
						# model_surfaces.append(O)
						simulations.append([O,taus])



	##### Apply affine diffeormorphisms and their inverses to each simulation and update marked segments #####
	for iteration in range(iteration_limit):
		print_percent_complete('Iteration',iteration,iteration_limit)
		remove_these_sims=[]
		for j in range(len(simulations)):
			print_percent_complete('Simulation',j,len(simulations))
			sim=simulations[j]
			O=sim[0]
			for k in range(len(generators)):
				A=generators[k]
				tau=sim[1][k]
				A_O=affine_diffeo(O,A,tau)
				try:
					O.combine_marked_segments(A_O)
				except:
					remove_these_sims.append(sim)
					break
				A_inv_O=affine_diffeo_inverse(O,A,tau)
				try:
					O.combine_marked_segments(A_inv_O)
				except:
					remove_these_sims.append(sim)
					break
		for sim in remove_these_sims:
			simulations.remove(sim)			

	# 'Un-normalize' by Trans(O); see note above
	for sim in copy(simulations):
		O=sim[0]
		taus=sim[1]
		for tau in O_blank.Trans:
			if tau!=O_blank.Trans.one():
				new_sim=[affine_diffeo(O,matrix.identity(2),tau),taus]
				##### UPDATE CONVEX BODIES HERE??? #####
				simulations.append(new_sim)

	##### Try to construct surfaces from remaining simulations #####
	# Sort simulations by taus
	simulations.sort(key=lambda x:x[1])		
	# Group together (into lists) simulations with common taus
	simulations_by_taus=[]
	sims_with_same_taus=[]
	taus=simulations[0][1]
	for sim in simulations:
		if sim[1]==taus:
			sims_with_same_taus.append(sim)
		else:
			simulations_by_taus.append(sims_with_same_taus)
			taus=sim[1]
			sims_with_same_taus=[sim]
	simulations_by_taus.append(sims_with_same_taus)

	surfaces=[]
	# for sims_with_same_taus in simulations_by_taus:
	for j in range(len(simulations_by_taus)):
		sims_with_same_taus=simulations_by_taus[j]
		print_percent_complete('simulations_by_taus',j,len(simulations_by_taus))
		# sims_with_same_taus is a list of simulations constructed with the same taus, and sims is a subset of these
		# for sims in list(powerset(sims_with_same_taus)):
		for k in range(len(list(powerset(sims_with_same_taus)))):
			sims=list(powerset(sims_with_same_taus))[k]
			print_percent_complete('sims_with_same_taus',k,len(list(powerset(sims_with_same_taus))))
			### MAYBE 'UN-NORMALIZE' HERE RATHER THAN ABOVE TO REDUCE REDUNDANCY ###
			if sims!=[]:
				surfaces+=find_surfaces(sims)
	return(surfaces)

			# ### (I) FIND POSSIBLE (normalized) VORONOI STAPLES ###
			# # Will hold convex body points from each model in models
			# convex_body_pts=[[[] for j in range(d+1)] for d in stratum]
			# for model in models
			# 	model.update_convex_body_pts()

			# # O=ModelSurface(stratum)
			# # for model in models:
			# # 	model.update_convex_body_pts()
			# # 	# model_copy will be a model surface with just the convex body points from model.  We will combine all of the model_copy's to form O
			# # 	model_copy=ModelSurface(stratum)
			# # 	for cbp in model.convex_body_pts_orientation_paired_representatives():

			# ### (II) GIVEN SUCH, FIND POSSIBLE SCALARS ###



def find_surfaces(sims):
	'''Given a list of simulations, search for scalars to apply to each simulation so that the marked points of each scaled simulation determine (convex bodies giving rise to) a translation surface.'''
	return_these=[]
	# Place all marked points contributing to edges of convex bodies in their respective simulations in a list (with sublists corresponding to components)
	stratum=sims[0][0].orders_of_zeros
	convex_bodies_pts=[[] for d in stratum]
	# Flag to end construction if needed
	end=False
	for sim in sims:
		o=sim[0]
		o.update_convex_bodies_pts()
		for comp_index in range(len(o.convex_bodies_pts)):
			cumulative_angles_of_pts=[mp.angle_cumulative() for mp in convex_bodies_pts[comp_index]]
			for eq_cls in o.convex_bodies_pts[comp_index]:
				# Check if any marked point in eq_cls lies on the same ray from the origin as any existing marked point in o.convex_bodies_pts[comp_index]; if so, either the simulations cannot give rise to a translation surface, or two simulations really correspond to the same Aff^+_O(X,omega)-orbit of marked points.  In either case, we end this construction
				for mp in eq_cls:
					if mp.angle_cumulative() in cumulative_angles_of_pts:
						end=True
						break
				if end==True:
					break
				else:
					convex_bodies_pts[comp_index]+=eq_cls
			if end==True:
				break
		if end==True:
			break

	if end==False:
		# Of the marked points in convex_bodies_pts, we must find a subset VS ('Voronoi Staples') for which:
		## (i) Elements come in orientation-pairs
		## (ii) The convex bodies formed by elements of VS are compact
		## (iii) Vertices of the convex bodies are regular points of the resulting translation surface (after appropriately scaling simulations and identifying edges)

		# List of all points in convex_bodies_pts
		convex_bodies_pts_all=[]
		for comp in convex_bodies_pts:
			convex_bodies_pts_all+=comp

		## (i) Elements come in orientation-pairs
		# List of orientation-paired representatives from convex_bodies_pts_all whose inverse is also in convex_bodies_pts_all
		convex_bodies_pts_paired_reps=[]
		for mp in convex_bodies_pts_all:
			if mp.inverse in convex_bodies_pts_all and mp.inverse not in convex_bodies_pts_paired_reps:
				convex_bodies_pts_paired_reps.append(mp)

		## (ii) The convex bodies formed by elements of VS are compact
		# We consider subsets of convex_bodies_pts_paired_reps (and together with their inverses) as the possible set VS
		# For (ii), a necessary condition is that each open half-plane of each component contain an element of VS; hence the minimum size of VS is \sum_{i=1}^\kappa(2d_i+3)
		##### THINK ABOUT UPPER BOUND AS WELL #####
		min_size_VS=sum([2*d+3 for d in stratum])
		for VS_reps in list(powerset(convex_bodies_pts_paired_reps)):
			# Flag to go to next element of power set
			skip=False
			if 2*len(VS_reps)>=min_size_VS:
				VS_all=VS_reps+[mp.inverse for mp in VS_reps]
				VS_by_component=[[] for d in stratum]
				for mp in VS_all:
					VS_by_component[mp.belongs_to_component.index_in(mp.belongs_to_component.component_of)].append(mp)
				# Sort components by cumulative angle, and check if subsequent angles are less than pi
				for comp in VS_by_component:
					comp.sort(key=lambda x:x.angle_cumulative())
					if (comp[0].angle_cumulative()+comp[0].belongs_to_component.num_planes)-comp[-1].angle_cumulative()>=1/2:
						skip=True
						break
					# assert((comp[0].angle_cumulative()+comp[0].belongs_to_component.num_planes)-comp[-1].angle_cumulative()<1/2)
					for j in range(1,len(comp)):
						mp0=comp[j]
						mp1=comp[j-1]
						if mp0.angle_cumulative()-mp1.angle_cumulative()>=1/2:
							skip=True
							break
					if skip==True:
						break
						# assert(mp0.angle_cumulative()-mp1.angle_cumulative()<1/2)
			else:
				skip=True

		## (iii) Vertices of the convex bodies are regular points of the resulting translation surface (after appropriately scaling simulations and identifying edges)
			if skip==False:
				VS_all_copy=copy(VS_all)
				while len(VS_all_copy)>0:
					# Will count the interior angle of a vertex of convex bodies after identifiying edges
					orbit_angle=0
					# Choose an initial marked point in the 'gluing orbit'
					mp_init=VS_all_copy[0]
					mp0=mp_init
					first=True
					while mp0!=mp_init or first==True:
						first=False
						VS_all_copy.remove(mp0)
						mp0_comp_index=mp0.belongs_to_component.index_in(mp0.belongs_to_component.component_of)
						mp0_index_in_comp=VS_by_component[mp0_comp_index].index(mp0)
						# Subsequent, counter-clockwise marked point from mp0
						mp1=VS_by_component[mp0_comp_index][(mp0_index_in_comp+1)%len(VS_by_component[mp0_comp_index])]
						# Interior angle of the perpendicular bisectors of mp0, mp1
						interior_ang=angle((mp0.holonomy[1],-mp0.holonomy[0]))-angle((-mp1.holonomy[1],mp1.holonomy[0]))
						if interior_ang<0:
							interior_ang+=1
						orbit_angle+=interior_ang
						mp0=mp1.inverse
					if orbit_angle!=1:
						skip=True
						break

			if skip==False:
				# We have a potential set VS; now we need positive scalars such that 
				## (iv) scalars are constant on simulations
				VS_by_sim=[[] for j in range(len(sims))]
				model_surfaces=[sim[0] for sim in sims]
				for vs in VS_all:
					o=vs.belongs_to_component.component_of
					VS_by_sim[model_surfaces.index(o)].append(vs)
				# We will store the necessary relations in a matrix M
				M=matrix(1,len(VS_all))
				first=True	
				for VS_in_sim in VS_by_sim:
					if len(VS_in_sim)>0:
						vs0=VS_in_sim[QQ(0)]
						vs0_index=VS_all.index(vs0)
						for vs1 in VS_in_sim[1:]:
							vs1_index=VS_all.index(vs1)
							new_row=[0]*len(VS_all)
							new_row[vs0_index]=1
							new_row[vs1_index]=-1
							if first==True:
								M.set_row(0,new_row)
								first=False
							else:						
								M=M.insert_row(M.nrows(),new_row)


			
				## (v) the marked points from VS are the only marked points from convex_bodies_pts contributing nontrivially to edgs of convex bodies

				## (vi) identified edges of convex bodies are of equal length
				for vs in VS_reps:
					vs_comp_index=vs.belongs_to_component.index_in(vs.belongs_to_component.component_of)
					vs_inverse_comp_index=vs.inverse.belongs_to_component.index_in(vs.inverse.belongs_to_component.component_of)
					# Adjacent Voronoi staples clockwise and counterclockwise, respectively, of vs and vs.inverse, resp.
					vs_cw=VS_by_component[vs_comp_index][(VS_by_component[vs_comp_index].index(vs)-1)%len(VS_by_component[vs_comp_index])]
					vs_ccw=VS_by_component[vs_comp_index][(VS_by_component[vs_comp_index].index(vs)+1)%len(VS_by_component[vs_comp_index])]
					vs_inverse_cw=VS_by_component[vs_inverse_comp_index][(VS_by_component[vs_inverse_comp_index].index(vs.inverse)-1)%len(VS_by_component[vs_inverse_comp_index])]
					vs_inverse_ccw=VS_by_component[vs_inverse_comp_index][(VS_by_component[vs_inverse_comp_index].index(vs.inverse)+1)%len(VS_by_component[vs_inverse_comp_index])]				
					# Indices of vs, vs_cw, vs_ccw, vs_inverse_cw, vs_inverse_ccw in VS_all
					vs_index=VS_all.index(vs)
					vs_cw_index=VS_all.index(vs_cw)
					vs_ccw_index=VS_all.index(vs_ccw)
					vs_inverse_cw_index=VS_all.index(vs_inverse_cw)
					vs_inverse_ccw_index=VS_all.index(vs_inverse_ccw)

					new_row=vector([QQbar(0)]*len(VS_all))
					new_row[vs_index]+=determ(vs.inverse.holonomy,vs_inverse_ccw.holonomy)*determ(vs.inverse.holonomy,vs_inverse_cw.holonomy)*(determ(vs.holonomy,vs_ccw.holonomy)*dot(vs.holonomy,vs_cw.holonomy)-determ(vs.holonomy,vs_cw.holonomy)*dot(vs.holonomy,vs_ccw.holonomy))-determ(vs.holonomy,vs_ccw.holonomy)*determ(vs.holonomy,vs_cw.holonomy)*(determ(vs.inverse.holonomy,vs_inverse_ccw.holonomy)*dot(vs.inverse.holonomy,vs_inverse_cw.holonomy)-determ(vs.inverse.holonomy,vs_inverse_cw.holonomy)*dot(vs.inverse.holonomy,vs_inverse_ccw.holonomy))
					new_row[vs_cw_index]+=-determ(vs.inverse.holonomy,vs_inverse_ccw.holonomy)*determ(vs.inverse.holonomy,vs_inverse_cw.holonomy)*determ(vs.holonomy,vs_ccw.holonomy)*dot(vs_cw.holonomy,vs_cw.holonomy)
					new_row[vs_ccw_index]+=determ(vs.inverse.holonomy,vs_inverse_ccw.holonomy)*determ(vs.inverse.holonomy,vs_inverse_cw.holonomy)*determ(vs.holonomy,vs_cw.holonomy)*dot(vs_ccw.holonomy,vs_ccw.holonomy)
					new_row[vs_inverse_cw_index]+=determ(vs.holonomy,vs_ccw.holonomy)*determ(vs.holonomy,vs_cw.holonomy)*determ(vs.inverse.holonomy,vs_inverse_ccw.holonomy)*dot(vs_inverse_cw.holonomy,vs_inverse_cw.holonomy)
					new_row[vs_inverse_ccw_index]+=-determ(vs.holonomy,vs_ccw.holonomy)*determ(vs.holonomy,vs_cw.holonomy)*determ(vs.inverse.holonomy,vs_inverse_cw.holonomy)*dot(vs_inverse_ccw.holonomy,vs_inverse_ccw.holonomy)
					N=M.rows()
					N.append(new_row)
					M=matrix(N)
				
				##### IF M.right_nullity()==0, THESE VS_reps DON'T GIVE RISE TO A TRANSLATION SURFACE #####
				##### IF M.right_nullity()==1, CONSTRUCT THE RESULTING TRANSLATION SURFACE #####
				if M.right_nullity()==1:
					O=ModelSurface(stratum)
					for staple in VS_reps:
						# scalar=scalars[model_surfaces.index(staple.belongs_to_component.component_of)]
						# scaled_holonomy=[scalar*staple.holonomy[0],scalar*staple.holonomy[1]]
						# O.add_orientation_paired_marked_segments(scaled_holonomy,staple.component_index(),staple.plane_index(),staple.inverse.component_index(),staple.inverse.plane_index())
						O.add_orientation_paired_marked_segments(staple.holonomy,staple.belongs_to_component.index_in(staple.belongs_to_component.component_of),staple.belongs_to_plane.index_in(staple.belongs_to_component),staple.inverse.belongs_to_component.index_in(staple.inverse.belongs_to_component.component_of),staple.inverse.belongs_to_plane.index_in(staple.inverse.belongs_to_component))
					try:
						X=O.translation_surface()
						X=X.canonicalize()
						if X not in return_these:
							# Check if Veech group of X contains generators
							contains_generators=True
							for A in generators:
								if X!=(A*X).canonicalize():
									contains_generators=False
									break
							if contains_generators==True:
								return_these.append(X)

					except:
						continue
				##### IF M.right_nullity()>1, USE DELAUNAY TRIANGLES TO DETERMINE SCALARS #####
				elif M.right_nullity()>1:
					scalars=[None for sim in sims]
					scalars[0]=1
					potential_scalars=[scalars]
					construct_with_these_scalars=[]
					while len(potential_scalars)>0:
						scalars=potential_scalars.pop()
						# We will find new scalars inductively
						for new_scalar_index in range(len(scalars)):
							# We seek new scalars (whose values are currently 'None') which correspond to simulations containing an element of VS_reps
							if scalars[new_scalar_index]!=None or new_scalar_index not in [model_surfaces.index(vs.belongs_to_component.component_of) for vs in VS_reps]:
								continue
							use_this_vs=False 
							j=-1
							# Find a Voronoi staple whose correpsonding scalar has yet to be determined, and whose adjacent, clockwise staple has been determined
							while use_this_vs==False:
								j+=1
								vs=VS_reps[j]
								# Check if this vs corresponds to the simulation whose scalar we are after
								if model_surfaces.index(vs.belongs_to_component.component_of)==new_scalar_index:
									vs_comp_index=vs.belongs_to_component.index_in(vs.belongs_to_component.component_of)
									vs_cw=VS_by_component[vs_comp_index][(VS_by_component[vs_comp_index].index(vs)-1)%len(VS_by_component[vs_comp_index])]
									# Check if the scalar of vs_cw has already been determined
									vs_cw_scalar=scalars[model_surfaces.index(vs_cw.belongs_to_component.component_of)]
									if vs_cw_scalar!=None:
										use_this_vs=True 
							angle_from_vs_cw_to_vs=vs.angle-vs_cw.angle
							if angle_from_vs_cw_to_vs<0:
								angle_from_vs_cw_to_vs+=1
							# Find all marked points mp from current simulations which could give a triangle (vs_cw.inverse,0,mp) isometric to the triangle (0,vs_cw,vs) after appropriate scalings of simulations; both of these triangles are isometric to a Delaunay triangle on the translation surface we're constructing
							potential_third_points=[]
							vs_cw_inverse_comp_index=vs_cw.inverse.belongs_to_component.index_in(vs_cw.inverse.belongs_to_component.component_of)
							for o in model_surfaces:
								for segments_in_plane in o.marked_segments_by_component_and_plane()[vs_cw_inverse_comp_index]:
									potential_third_points+=segments_in_plane
							# mp must be within pi-angle_from_vs_cw_to_vs clockwise of vs_cw.inverse							
							min_possible_angle=vs_cw.inverse.angle_cumulative()-(1/2-angle_from_vs_cw_to_vs)
							for mp in copy(potential_third_points):
								if min_possible_angle>=0 and (mp.angle_cumulative()<=min_possible_angle or mp.angle_cumulative()>=vs_cw.inverse.angle_cumulative()):
									potential_third_points.remove(mp)
									continue
								elif min_possible_angle<0 and (mp.angle_cumulative()<=min_possible_angle+stratum[vs_cw_inverse_comp_index]+1 and mp.angle_cumulative()>=vs_cw.inverse.angle_cumulative()):
									potential_third_points.remove(mp)
									continue
								# mp.inverse must belong to the same component as vs.inverse
								if mp.inverse.belongs_to_component.index_in(mp.inverse.belongs_to_component.component_of)!=vs.inverse.belongs_to_component.index_in(vs.inverse.belongs_to_component.component_of):
									potential_third_points.remove(mp)
									continue 
								# mp must satisfy that the angles from (i) vs_cw to vs, (ii) mp to vs_cw.inverse and (iii) vs.inverse to mp.inverse sum to pi
								angle_from_mp_to_vs_cw_inverse=vs_cw.inverse.angle_cumulative()-mp.angle_cumulative()
								if angle_from_mp_to_vs_cw_inverse<0:
									angle_from_mp_to_vs_cw_inverse+=stratum[mp.belongs_to_component.index_in(mp.belongs_to_component.component_of)]
								angle_from_vs_inverse_to_mp_inverse=mp.inverse.angle_cumulative()-vs.inverse.angle_cumulative()
								if angle_from_vs_inverse_to_mp_inverse<0:
									angle_from_vs_inverse_to_mp_inverse+=stratum[vs.inverse.belongs_to_component.index_in(vs.inverse.belongs_to_component.component_of)]
								if angle_from_vs_cw_to_vs+angle_from_mp_to_vs_cw_inverse+angle_from_vs_inverse_to_mp_inverse!=1/2:
									potential_third_points.remove(mp)
							# Update scalars
							for mp in potential_third_points:
								scalars_copy=copy(scalars)
								vs_scalar=vs_cw_scalar*determ(mp.holonomy,vs_cw.holonomy)/determ(mp.holonomy,vs.holonomy)
								mp_scalar=vs_cw_scalar*determ(vs.holonomy,vs_cw.holonomy)/determ(mp.holonomy,vs.holonomy)
								# Update scalar corresponding to simulation of vs
								scalars_copy[new_scalar_index]=vs_scalar
								# Update scalar corresponding to simulation of ms if possible. Note that this scalar may already have been found, and may disagree with mp_scalar above.  In this case, mp cannot be a vertex of the desired triangle
								if scalars_copy[model_surfaces.index(mp.belongs_to_component.component_of)]==None or scalars_copy[model_surfaces.index(mp.belongs_to_component.component_of)]==mp_scalar:
									scalars_copy[model_surfaces.index(mp.belongs_to_component.component_of)]=mp_scalar
									if None in scalars_copy:
										construct=True
										for k in range(len(scalars)):
											if scalars_copy[k]==None and k in [model_surfaces.index(x.belongs_to_component.component_of) for x in VS_reps]:
												construct=False 
												break 
										if construct==True:
											# Make sure the scalars from scalar_copy that correspond to simulations with Voronoi staples belong to the kernel of M
											v=[0]*len(VS_all)
											for k in range(len(VS_all)):
												staple=VS_all[k]
												staple_sim_index=model_surfaces.index(staple.belongs_to_component.component_of)
												v[k]=scalars_copy[staple_sim_index]
											if M*matrix(v).transpose()==matrix.zero(M.nrows(),1):
												construct_with_these_scalars.append(scalars_copy)
										else:
											potential_scalars.append(scalars_copy)
										# potential_scalars.append(scalars_copy)
									else:
										# Make sure the scalars from scalar_copy that correspond to simulations with Voronoi staples belong to the kernel of M
										v=[0]*len(VS_all)
										for k in range(len(VS_all)):
											staple=VS_all[k]
											staple_sim_index=model_surfaces.index(staple.belongs_to_component.component_of)
											v[k]=scalars_copy[staple_sim_index]
										if M*matrix(v).transpose()==matrix.zero(M.nrows(),1):
											construct_with_these_scalars.append(scalars_copy)
					# Scale and combine the simulations, and construct the resulting surfaces
					for scalars in construct_with_these_scalars:
						O=ModelSurface(stratum)
						for staple in VS_reps:
							scalar=scalars[model_surfaces.index(staple.belongs_to_component.component_of)]
							scaled_holonomy=[scalar*staple.holonomy[0],scalar*staple.holonomy[1]]
							# O.add_orientation_paired_marked_segments(scaled_holonomy,staple.component_index(),staple.plane_index(),staple.inverse.component_index(),staple.inverse.plane_index())
							O.add_orientation_paired_marked_segments(scaled_holonomy,staple.belongs_to_component.index_in(staple.belongs_to_component.component_of),staple.belongs_to_plane.index_in(staple.belongs_to_component),staple.inverse.belongs_to_component.index_in(staple.inverse.belongs_to_component.component_of),staple.inverse.belongs_to_plane.index_in(staple.inverse.belongs_to_component))

						try:
							X=O.translation_surface()
							X=X.canonicalize()
							if X not in return_these:
								# Check if Veech group of X contains generators
								contains_generators=True
								for A in generators:
									if X!=(A*X).canonicalize():
										contains_generators=False
										break
								if contains_generators==True:
									return_these.append(X)							
						except:
							continue
	return(return_these)


				# for vs in VS_reps:
				# 	vs_comp_index=vs.belongs_to_component.index_in(vs.belongs_to_component.component_of)
				# 	vs_inverse_comp_index=vs.inverse.belongs_to_component.index_in(vs.inverse.belongs_to_component.component_of)
				# 	# Adjacent Voronoi staples clockwise and counterclockwise, respectively, of vs and vs.inverse, resp.
				# 	vs_cw=VS_by_component[vs_comp_index][(VS_by_component[vs_comp_index].index(vs)-1)%len(VS_by_component[vs_comp_index])]
				# 	vs_ccw=VS_by_component[vs_comp_index][(VS_by_component[vs_comp_index].index(vs)+1)%len(VS_by_component[vs_comp_index])]
				# 	vs_inverse_cw=VS_by_component[vs_inverse_comp_index][(VS_by_component[vs_inverse_comp_index].index(vs.inverse)-1)%len(VS_by_component[vs_inverse_comp_index])]
				# 	vs_inverse_ccw=VS_by_component[vs_inverse_comp_index][(VS_by_component[vs_inverse_comp_index].index(vs.inverse)+1)%len(VS_by_component[vs_inverse_comp_index])]				
				# 	# Indices of vs, vs_cw, vs_ccw, vs_inverse_cw, vs_inverse_ccw in VS_all
				# 	vs_index=VS_all.index(vs)
				# 	vs_cw_index=VS_all.index(vs_cw)
				# 	vs_ccw_index=VS_all.index(vs_ccw)
				# 	vs_inverse_cw_index=VS_all.index(vs_inverse_cw)
				# 	vs_inverse_ccw_index=VS_all.index(vs_inverse_ccw)








def determ(u,v):
	u_x=u[0]
	u_y=u[1]
	v_x=v[0]
	v_y=v[1]
	return u_x*v_y-u_y*v_x

def dot(u,v):
	u_x=u[0]
	u_y=u[1]
	v_x=v[0]
	v_y=v[1]
	return u_x*v_x+u_y*v_y	


			




				


	# 			convex_bodies_pts[comp_index]+=eq_cls


	# # Sort points in convex_body_pts by cumulative angle
	# for comp in convex_bodies_pts:
	# 	comp.sort(key=lambda x:x.angle_cumulative())

	# # Check if any marked points lie on the same ray from the origin





