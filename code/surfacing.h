///////////////////////////////////////////////////////////////////////////////
//  Demo code for the core iterative process of flow aligned surfacing of curve networks.
//  Author: Hao Pan
///////////////////////////////////////////////////////////////////////////////

#include "Mesh3D.h"


class Surfacing{

public:

    //adapt the tri mesh by iteratively computing cross field and solving for new mesh
	//whose principal directions better match the cross field.
	void iterative_surfacing();

protected:

    //build the curvature curve aligned smooth cross field on mesh.
    void curve_aligned_cross_field(Mesh3D* p_mesh);

    //adapt the mesh to given cross field on the mesh.
    //assume the cross field data is already stored on facets and edges of the mesh.
    void adapt_surface_to_df(Mesh3D* p_mesh);


    //mark facets with degenerate area, by setting tag = false.
	void mark_degenerate_facets(Mesh3D* p_mesh);

	//fit curvature tensor for each facet.
	void fit_curvature_tensor(Mesh3D* p_mesh);

	//smooth principal curvatures.
    void smooth_curvature_variation(Mesh3D* p_mesh);

	//compute target curvatures and update surface.
	void match_target_curvature(Mesh3D* p_mesh);


    //compute per edge dihedral diff gradient wrt. vertices.
    double dihedral_diff_gradient(HE_edge* edge, double target_sinhd, std::vector<HE_vert*>& vts_list, std::vector<vec3>& grad_list);
    //fill in linear system for the edge's target curvature normal.
    void edgecurv_fillmatrix(HE_edge* edge, vec3 targt_curvnorm, Sparse_Matrix& mat, double coeff = 1.0);
    //fill matrix with RMF normal constraint.
    void fillmatrix_normalConstraint(HE_edge* edge, Sparse_Matrix& mat);
    //get facet based differential data.
    void facet_differential(const Mesh3D* pmesh, const HE_face* face, const vec3 u_axis, const vec3 v_axis, 
                std::vector<vec2>& dx_list, std::vector<vec2>& dn_list, const HE_face* prev_face = NULL);
    


    //mesh to evolve
    Mesh3D* m_mesh_;

}

	


    