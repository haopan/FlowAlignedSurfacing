#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <cmath> 
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <queue>
#include <cstdlib>


#include "Geex\mathematics\glsl_linear.h"

#include "CurvNet.h"

#include "WildMagic5\Include\Wm5Vector2.h"
#include "WildMagic5\Include\Wm5Vector3.h"


using namespace Wm5;


//declare classes for the compiler
class HE_vert;
class HE_face;
class HE_edge;
class Mesh3D;

//boundary data classes.
class BD_Corner{
public:
	HE_vert* p_corner_pos;
	HE_edge* p_corner_drv;
	int step;
	BD_Corner(){
		p_corner_pos = NULL;
		p_corner_drv = NULL;
		step = 0;
	}
	BD_Corner(HE_vert* cptr, HE_edge* eptr, int u){
		p_corner_pos = cptr;
		p_corner_drv = eptr;
		step = u;
	}
};

class BD_Seg{
public:
	BD_Corner two_ends[2];

	HE_edge* p_bdseg_drv;
	HE_vert* p_bdseg_pos;
	int step;

private:
	int counter;
public:
	BD_Seg(){
		counter = 0;
		p_bdseg_drv = NULL;
		p_bdseg_pos = NULL;
		step = 0;
	}

	BD_Seg(HE_vert* bd_v, HE_edge* bd_dir){
		counter = 0;
		p_bdseg_pos = bd_v;
		p_bdseg_drv = bd_dir;
		step = 0;
	}

	void push_corner(BD_Corner& cn){
		assert(counter<2);
		two_ends[counter++] = cn;
	}

	bool intersect_at(BD_Seg& other, bool isFront){
		if (other.two_ends[0].p_corner_pos == two_ends[isFront?0:1].p_corner_pos)
		{
			return true;
		}
		else if (other.two_ends[1].p_corner_pos == two_ends[isFront?0:1].p_corner_pos)
		{
			std::swap(other.two_ends[0], other.two_ends[1]);
			return true;
		}
		return false;
	}
};

enum PatchType {REGULAR_PATCH, SPHERICAL_PATCH, UNDEFINED_PATCH} ;

enum VertexType {FREE_VTX, FIXED_VTX, APPROX_BOUNDARY_VTX};


/*!
*	The basic vertex class for half-edge structure.
*/

class HE_vert
{
public:
	vec3 pos;		//!< 3D coordinate
	HE_edge* edge;			//!< one of the half-edges_list emanating from the vertex
	vec3 normal;	    //!< vertex normal
	ptrdiff_t id;					//!< index 
	unsigned int degree;			//!< the degree of vertex 
	bool tag;				//!< tag for programming easily

	vec3 pdir_min;		//principal direction, the minor axis.
	vec3 pdir_max;		//principal direction, the major axis.

	int singular_idx;

	bool  frame_fixed;// !< indicate whether the frame of cyclide patches are fixed , bob

	//bool on_boundary;
	bool fixed;	

	bool fairing;

	bool selected;

	int var_id; //store variable id, hpan.

	double quality; //holder for various quality measures.

	//double m_funcvalue_for_colorcode[COLORCODEING_COUNT];// for color coding show , see VERTEX_CONST_TYPE

	std::vector<HE_edge*> pg_edges; //used for polygraph construction only.

	double	m_fCurvePara; //parameter of corresponding point on curve. hpan
	SplineCurv*	m_pParaCurve; //pointer to curve with corresponding point. hpan
	int	m_dOnCurveParaID; //id of parameter variable, used for indexing the array of variables. hpan

	enum VertexType m_vtype; //type of the vertex. hpan

	std::vector<BD_Seg> hermite_data;	//to keep the boundary Hermite data. hpan

public:
	//! constructor
	HE_vert(const vec3 &v)
		:pos(v), edge(0), id(-1), degree(0), tag(false), frame_fixed(false), fixed(false), selected(false), fairing(true), pdir_min(vec3(0,0,0)),
		pdir_max(vec3(0,0,0)), m_fCurvePara(0), m_pParaCurve(NULL), m_dOnCurveParaID(-1), m_vtype(VertexType::FREE_VTX), singular_idx(0), var_id(-1)
	{

		quality =  0.0;
	}
	//destructor
	~HE_vert()
	{
	}
};


/*!
*	The basic edge class for half-edge structure.
*/

class HE_edge
{ 
public:
	HE_vert* vert;	//!< vertex at the end of the half-edge
	HE_edge* pair;	//!< oppositely oriented adjacent half-edge 
	HE_face* face;	//!< face the half-edge borders
	HE_edge* next;	//!< next half-edge around the face
	HE_edge* prev;	//!< prev half-edge around the face
	ptrdiff_t id;			//!< index
	bool tag;				//!< tag for programming easily

	bool selected;

	//SplineCurv* spln_curv; //pointer to the corresponding spline curve
	std::vector<HE_vert*> sample_vh; //array of sample vertices. for merging subdivided polygraph patches.

	int pij; //period jump of cross directions.
	double kij; //chain curvature along dual edge.

	int var_id; //variable id.

	vec3 ortho_dual; //the orthogonal dual. hpan

	vec3 debug_vec; //displaying some quantities for debug. hpan
	vec3 debug_target_vec;
	vec3 debug_cur_vec;

	double target_sin_hdi;
	double dual_len;

	PolyGraph_Edge* cycle_curv;

	//cyclide computation.
	bool fixed;
	//vec3 contactpoint;
	//double m_funcvalue_for_colorcode[COLORCODEING_COUNT];
	vec3 tangent;
	//bool m_bCutOff;
	//int m_interpolation_id;
	//int m_free_id;
	//double m_tangent_parameter[2];
	//std::vector<Wm5::Vector2d>  m_curvatureline_pointuv;
	//std::vector<int >  m_sample_id;
	//std::vector<SURF_POINT>  m_curvatureline_sample_points;
	//std::vector<std::vector<Vector2d>>   m_face_curvatureline_list;//curvature lines on surface, with starting point on this edge

public:

	//!constructor
	HE_edge()
		:vert(0), pair(0), face(0), next(0), prev(0), id(-1), tag(false), selected(false), /*spln_curv(NULL),*/ fixed(false), tangent(vec3(0,0,0)),
		pij(0), kij(0), var_id(-1), target_sin_hdi(0), dual_len(0.0)
	{

		ortho_dual = vec3(0,0,0);
		debug_vec = vec3(0,0,0);
		debug_target_vec = debug_cur_vec = vec3(0,0,0);

		cycle_curv = NULL;

		//for(int i = 0; i < COLORCODEING_COUNT; i++)
		//	m_funcvalue_for_colorcode[i] =  0;

		//m_tangent_parameter[0] = m_tangent_parameter[1] = 0;
	}
	//!destructor
	~HE_edge()
	{
	}

	//! compute the middle point
	inline vec3 GetMidPoint()
	{
		return 0.5*(vert->pos+pair->vert->pos);
	}
};



/*!
*	The basic face class for half-edge structure.
*/

class HE_face
{ 
public:
	HE_edge* edge;		//!< one of the half-edges_list bordering the face	
	unsigned int valence;		//!< the number of edges_list 
	vec3 normal;	//!< face normal
	int id;					    //!< index
	bool tag;				    //!< tag for programming easily
	std::vector<ptrdiff_t>  texture_indices; //! texture indices
	std::vector<ptrdiff_t>  normal_indices; //! texture indices
	bool fixed;

	bool selected;

	double theta[2]; //cross field angle. use two for conjugate field. hpan.
	int var_id[2]; //var id. hpan.

	vec3 tensor_axis [2]; //axis of estimated curvature tensor. i.e. the principal directions.
	double tensor_diag [2]; //the estimated principal curvature values.

	//debug. hpan
	vec3 unrestricted_tensor [2]; //principal curvature directions.
	double unrestricted_pcv [2]; //principal curvature values.

	PolyGraph_Cycle* curve_cycle; //curve cycle this facet belongs to.

	double quality; //holder for quality measurements.
	bool is_umbilic;

private:
	//FLAG_TYPE m_flags; //!< face-flags


public:
	//!constructor
	HE_face()
		:edge(0), id(-1), tag(false), fixed(false), selected(false)
	{
		var_id[0] = var_id[1] = -1;
		theta[0] = 0.0;
		theta[1] = M_PI*0.5;

		tensor_axis[0] = tensor_axis[1] = vec3(0,0,0);
		tensor_diag[0] = tensor_diag[1] = 0;

		unrestricted_tensor[0] = unrestricted_tensor[1] = vec3(0,0,0);

		curve_cycle = NULL;

		quality = 0.0;
		is_umbilic = false;

		unrestricted_pcv[0] = unrestricted_pcv[1] = 0.0;
	}
	//!destructor
	~HE_face()
	{
	}

	//! compute the barycenter
	inline vec3 GetCentroid()
	{
		vec3 V(0, 0, 0);
		HE_edge* he = edge;
		int i = 0;
		do 
		{
			V+= he->vert->pos;
			he = he->next;
			i++;
		} while(he!=edge);
		return V/double(i);
	}

	//! get facet area
	double get_area(){
		vec3 x[3];
		HE_edge* cur_edge = edge;
		for (int i = 0; i < 3; i ++)
		{
			x[i] = cur_edge->vert->pos;
			cur_edge = cur_edge->next;
		}

		return 0.5*cross(x[2]-x[0],x[1]-x[0]).length();
	}

	//judge if the facet is degenerate, by computing the facet area over edge length.
	bool is_degenerate(){
			vec3 x[3];
			HE_edge* cur_edge = edge;
			for (int i = 0; i < 3; i ++)
			{
				x[i] = cur_edge->vert->pos;
				cur_edge = cur_edge->next;
			}

			double edge_len = ((x[1]-x[0]).length() + (x[2]-x[1]).length() + (x[0]-x[2]).length())/3.0;

			double ratio = cross(x[2]-x[0],x[1]-x[0]).length()/(edge_len*edge_len);
			const double threshold = std::sin(5.0/180.0*M_PI);
			if (ratio < threshold)
			{
				return true;
			}

			return false;
	}

	//compute the axis of the cross at the facet.
	void get_cross_axis(vec3& u_axis, vec3& v_axis){
		vec3 faxis = normalize(edge->vert->pos - edge->prev->vert->pos);
		vec3 binormal = cross(normal, faxis);

		u_axis = vec3(0,0,0);
		v_axis = vec3(0,0,0);

		double cos_theta = std::cos(theta[0]), sin_theta = std::sin(theta[0]);

		u_axis = cos_theta*faxis + sin_theta*binormal;
		v_axis = -sin_theta*faxis + cos_theta*binormal;
	}

	//! compute index of the vertex in the face. hpan
	int index(const HE_vert* vh){
		int idx = 0;
		HE_edge* eh = edge;
		do 
		{
			if (eh->vert == vh)
			{
				return idx;
			}
			idx++;
			eh = eh->next;
		} while (eh != edge);
		return -1;
	}
	//! whether texture_indices exists
	bool has_texture_map()
	{
		return (!texture_indices.empty()) && (texture_indices.size() == (size_t)valence);
	}
	//! whether normal_indices exists
	bool has_normal_map()
	{
		return (!normal_indices.empty()) && (normal_indices.size() == (size_t)valence);
	}

	void copy_regular_tri_data(Mesh3D* p_mesh, HE_face* rhs);

	HE_edge* has_edge(int v0, int v1) {
		HE_edge* eh = edge;
		do
		{
			if (eh->prev->vert->id == v0 && eh->vert->id == v1)
			{
				return eh;
			}
			eh = eh->next;
		} while (eh != edge);

		return NULL;
	}
};

//////////////////////////////////////////////////////////////////////////

bool inline CompareEdgeID (HE_edge* he1, HE_edge* he2 )
{
	return he1->id < he2->id;
}

bool inline CompareVertexID (HE_vert* hv1, HE_vert* hv2 )
{
	return hv1->id < hv2->id;
}

bool inline CompareFaceID (HE_face* hf1, HE_face* hf2 )
{
	return hf1->id < hf2->id;
}


//////////////////////////////////////////////////////////////////////////

//! Mesh3D class: Half edge data structure \ingroup MeshCore
/*!
* a half-edge based mesh data structure
* For understanding half-edge structure, 
* please read the article in http://www.flipcode.com/articles/article_halfedge.shtml
*/

class Mesh3D
{

public:

	// type definition
	typedef std::vector<HE_vert* > VERTEX_LIST;
	typedef std::vector<HE_face* > FACE_LIST;
	typedef std::vector<HE_edge* > EDGE_LIST;

	typedef VERTEX_LIST* PTR_VERTEX_LIST;
	typedef FACE_LIST* PTR_FACE_LIST;
	typedef EDGE_LIST* PTR_EDGE_LIST;

	typedef  VERTEX_LIST::iterator VERTEX_ITER;
	typedef  FACE_LIST::iterator FACE_ITER;
	typedef  EDGE_LIST::iterator EDGE_ITER;

	typedef  VERTEX_LIST::reverse_iterator VERTEX_RITER;
	typedef  FACE_LIST::reverse_iterator FACE_RITER;
	typedef  EDGE_LIST::reverse_iterator EDGE_RITER;
	typedef std::pair<HE_vert*, HE_vert* > PAIR_VERTEX;

protected:

	// mesh data

	PTR_VERTEX_LIST vertices_list;		//!< store vertices
	PTR_EDGE_LIST edges_list;			//!< store edges
	PTR_FACE_LIST faces_list;			//!< store faces

	// mesh type

	bool m_closed;						//!< indicate whether the mesh is closed
	bool m_quad;						//!< indicate whether the mesh is quadrilateral
	bool m_tri;							//!< indicate whether the mesh is triangular
	bool m_hex;							//!< indicate whether the mesh is hexagonal
	bool m_pentagon;					//!< indicate whether the mesh is pentagonal

	//! associate two end vertices with its edge: only useful in creating mesh
	std::map<PAIR_VERTEX, HE_edge*> m_edgemap;	

	// mesh info

	int m_num_components;				//!< number of components
	int m_num_boundaries;				//!< number of boundaries
	int m_genus;						//!< the genus value

	bool m_encounter_non_manifold;


public:

	//! values for the bounding box
	double xmax, xmin, ymax, ymin, zmax, zmin;

	//! store all the boundary vertices, each vector corresponds to one boundary
	std::vector<std::vector<HE_vert*> > boundaryvertices;

public:
	//! constructor
	Mesh3D(void);

	//! destructor
	~Mesh3D(void);

	//! get the pointer of vertices list
	inline PTR_VERTEX_LIST get_vertices_list()
	{
		return vertices_list;
	}

	//! get the pointer of edges list
	inline PTR_EDGE_LIST get_edges_list()
	{
		return edges_list;
	}

	//! get the pointer of faces list
	inline PTR_FACE_LIST get_faces_list()
	{
		return faces_list;
	}

	//! get the total number of vertices
	inline ptrdiff_t get_num_of_vertices()
	{	
		return vertices_list?(ptrdiff_t)vertices_list->size():0;
	}

	//! get the total number of faces
	inline ptrdiff_t get_num_of_faces()
	{
		return faces_list?(ptrdiff_t)faces_list->size():0;
	}

	//! get the total number of half-edges
	inline ptrdiff_t get_num_of_edges()
	{
		return edges_list?(ptrdiff_t)edges_list->size():0;
	}

	//! get the pointer of the id-th vertex
	inline HE_vert* get_vertex(ptrdiff_t id) 
	{
		return id >= get_num_of_vertices()||id<0? NULL:(*vertices_list)[id];
	}

	//! get the pointer of the id-th edge
	inline HE_edge* get_edge(ptrdiff_t id) 
	{
		return id >= get_num_of_edges()||id<0? NULL: (*edges_list)[id];
	}

	//! get the pointer of the id-th face
	inline HE_face* get_face(ptrdiff_t id) 
	{
		return id >= get_num_of_faces()||id<0? NULL: (*faces_list)[id];
	}

	//! get the number of components
	inline int get_num_of_components()
	{
		return m_num_components;
	}

	//! get the number of boundaries
	inline int get_num_of_boundaries()
	{
		return m_num_boundaries;
	}

	//! get the genus
	inline int genus()
	{
		return m_genus;
	}

	//! check whether the mesh is valid
	inline bool is_valid() 
	{
		if( get_num_of_vertices()==0 || get_num_of_faces() == 0 )
			return false;
		return true;
	}

	//! check whether the mesh is closed
	inline bool is_closed() 
	{
		return m_closed;
	}

	//! check whether the mesh is triangular
	inline bool is_tri()
	{
		return m_tri;
	}

	//! check whether the mesh is quadrilateral
	inline bool is_quad()
	{
		return m_quad;
	}

	//! check whether the mesh is hexgaonal
	inline bool is_hex()
	{
		return m_hex;
	}

	//! check whether the mesh is pentagonal
	inline bool is_pentagon()
	{
		return m_pentagon;
	}

	//! insert a vertex 
	/*!
	*	\param v a 3d point
	*	\return a pointer to the created vertex
	*/
	HE_vert* insert_vertex(vec3 &v);

	//! insert a face
	/*!
	*	\param vec_hv the vertices list of a face
	*	\param texture the pointer of texture vector
	*	\param normal the normal of texture vector
	*	\return a pointer to the created face
	*/
	HE_face* insert_face(VERTEX_LIST& vec_hv, std::vector<ptrdiff_t> * texture = 0, std::vector<ptrdiff_t> * normal = 0);

	//! check whether the vertex is on border
	bool is_on_boundary(HE_vert* hv); 
	//! check whether the face is on border
	bool is_on_boundary(HE_face* hf); 
	//! check whether the edge is on border
	bool is_on_boundary(HE_edge* he); 

	//FILE IO

	//! load a 3D mesh from an OFF format file
	bool load_off(const char* fins);
	//! export the current mesh to an OFF format file
	void write_off(const char* fouts);
	//! load a 3D mesh from an OBJ format file
	bool load_obj(const char* fins);
	//! export the current mesh to an OBJ format file
	void write_obj(const char* fouts);
	//! export to a VTK format file
	void write_vtk(const char* fouts);

	//! update mesh:
	/*! 
	*	call it when you have created the mesh
	*/
	void update_mesh(); 

	//! update normal
	/*!
	*	compute all the normals of vertices and faces
	*/
	void update_normal(bool onlyupdate_facenormal = false);

	//! compute the bounding box
	void compute_boundingbox();

	//! get a copy of the current
	Mesh3D* make_copy();

	//get a copy of the current. no copying of edge id.
	Mesh3D* make_copy_simple();

	//! return a face-orientation-changed mesh
	Mesh3D* reverse_orientation();

	//! init edge tags
	/*!
	*	for a pair of edges, only one of them is tagged to be true.
	*/
	void init_edge_tag();

	//! reset all the vertices' tag
	void reset_vertices_tag(bool tag_status);
	//! reset all the faces' tag
	void reset_faces_tag(bool tag_status);
	//! reset all the edges' tag
	void reset_edges_tag(bool tag_status);
	//! reset all tag exclude edges' tag2
	void reset_all_tag(bool tag_status);

	//! translate the mesh with tran_V
	void translate(const vec3& tran_V);
	//! scale the mesh
	void scale(double factorx, double factory, double factorz);

	//! check whether there is any non-manifold case
	inline bool is_encounter_nonmanifold()
	{
		return m_encounter_non_manifold;
	}

	//! make the mesh Delaunay through edge-flipping.
	void make_delaunay( bool& topo_updated);

	//! check if edge is Delaunay.
	bool delaunay_test( HE_edge* triedge );

	//! swap edge in trimesh.
	void swap_edge( HE_edge* triedge );

	//! insert a vertex in triface.
	void split_face(HE_face* face, HE_vert* vtx);

	//build kd-tree for fast projection.
	void build_kdtree();

	//project pt to mesh.
	vec3 project_pt(vec3 pt, HE_face* & fh);


	HE_edge* pg_insert_edge(HE_vert* vstart, HE_vert* vend);
	HE_face* pg_insert_face(std::vector<HE_edge*>& edges);

	PolyGraph_Edge* vertex_on_curve(HE_vert* vh);
	int vertex_curves(HE_vert* vh);

	Mesh3D* make_copy_regular();

	//measure chain curvature of given edge eij
	double edge_chain_crvt(HE_edge* eij);

	//measure length of barycentric dual edge
	double bc_dual_edge_length(HE_edge* edge);

	//measure length of orthogonal dual edge, and the normal difference of two facets. return the weight used.
	double ortho_dual_edge_data(HE_edge* edge, vec3& dir, vec3& ndiff, bool use_cot_w = false) const;

	//get function and gradient of cvt_cmc energy.
	void get_cvt_fg( double& f, std::vector<double>& g );

	//get cvt centroid.
	vec3 get_cvt_centroid(HE_vert* vh);

	//initialize list of vertices, edges, faces.
	void init_vef_list();

	

private:

	//! insert an edge
	HE_edge* insert_edge(HE_vert* vstart, HE_vert* vend); 

	//! clear all the data
	void clear_data();

	//! clear vertices
	void clear_vertices();
	//! clear edges
	void clear_edges();
	//! clear faces
	void clear_faces();

	//! check whether the mesh is closed
	void check_closed();
	//! check the mesh type
	void check_meshtype();

	//! compute all the normals of faces
	void compute_faces_list_normal(); 
	//! compute the normal of a face
	void compute_perface_normal(HE_face* hf);
	//! compute all the normals of vertices
	void compute_vertices_list_normal();
	//! compute the normal of a vertex
	void compute_pervertex_normal(HE_vert* hv); 

	//! compute the number of components
	void compute_num_components();
	//! compute the number of boundaries
	void compute_num_boundaries();
	//! compute the genus
	void compute_genus();

	//! handle the boundary half edges specially
	void set_nextedge_for_border_vertices(); 

	//! remove the vertices which have no connection to others.
	void remove_hanged_vertices();

	//! align edges's id
	/*!
	set mesh's edge(vertex(startid), vertex(endid)->id = edgeid.
	only used in make_copy
	*/
	void copy_edge_id(ptrdiff_t edgeid, ptrdiff_t startid, ptrdiff_t endid, Mesh3D* mesh);


	//kd-tree for fast projection.
	ANNkd_tree* m_kdtree_;

	static const int facet_sample_num_ = 30;
};

	