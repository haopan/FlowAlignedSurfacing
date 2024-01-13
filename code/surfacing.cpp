#include "surfacing.h"


void Surfacing::iterative_surfacing(){

    if (!m_mesh_)
    {
        return;
    }

    if (!m_mesh_->is_tri())
    {
        std::cout << "The mesh is not triangular!!" << std::endl;
        return;
    }

    const int max_itr = 10;
    for (int itr_counter = 0; itr_counter < max_itr; itr_counter ++)
    {
        //remesh
        //make mesh delaunay to improve quality.
        {
            std::cout << "make mesh delaunay..." << std::endl;
            bool topo;
            m_mesh_->make_delaunay(topo);
        }
        m_mesh_->update_normal();

        //solve the smooth cross field.
        curve_aligned_cross_field(m_mesh_);
        
        //fit the mesh to given principal direction data.
        adapt_surface_to_df(m_mesh_, b_cross_conjugate_switch_);
    }

    m_mesh_->update_normal();
}





void Surfacing::curve_aligned_cross_field(Mesh3D* p_mesh){

    class FacetDist{
        public:
            HE_edge* coming_edge_;
            HE_face* fh_;
            double dist_;
    
            FacetDist(HE_edge* edge, HE_face* face, double dist){
                coming_edge_ = edge;
                fh_ = face;
                dist_ = dist;
            }

            bool operator< (const FacetDist& rhs) const{
                if (dist_ > rhs.dist_)
                {
                    return true;
                }
                return false;
            }
    };

    //solution vector x
    std::vector< double > x;

    //ids of variables to round to integer
    std::vector< int > ids_to_round;

    std::vector<HE_face*> map_var_fh;
    std::vector<HE_edge*> map_var_eh;

    //initialize num of variables and integer variables

    //std::queue<HE_face*> q_fixed_jumps;

    std::priority_queue<FacetDist> q_fixed_jumps;

    p_mesh->reset_edges_tag(false);
    p_mesh->reset_faces_tag(false);

    //initialize edge curvature.
    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
    {
        HE_edge* eh = p_mesh->get_edge(eitr);
        if (!eh->tag)
        {
            eh->tag = eh->pair->tag = true;
            
            eh->kij = p_mesh->edge_chain_crvt(eh);
            eh->pair->kij = -eh->kij;
        }
    }

    //bool debug_oneface = false;
    
    //get constrained and free facets.
    for (int fitr = 0; fitr < p_mesh->get_num_of_faces(); fitr++)
    {
        HE_face* fh = p_mesh->get_face(fitr);
        HE_edge* eh = fh->edge;

        bool constrained = false;
        fh->var_id[0] = -1;

        do
        {

            if ( (delaunay_->m_b_align_bdry_curv_ && !eh->pair->face) //constrain all boundary facets.
                ||
            (!delaunay_->m_b_align_bdry_curv_ && (eh->cycle_curv && eh->cycle_curv->curve_type == CURVE_TYPE::CRVT_LINE)) )
            {
                constrained = true;

                if (eh == fh->edge)
                {
                    fh->theta[0] = 0;
                }
                else
                {
                    vec3 target_dir = normalize(eh->vert->pos - eh->prev->vert->pos);
                    vec3 axis = normalize(fh->edge->vert->pos - fh->edge->prev->vert->pos);
                    
                    double dp = dot(axis, target_dir);

                    fh->theta[0] = std::acos(dp);

                    if (dot(cross(axis,target_dir), fh->normal) < 0)
                    {
                        fh->theta[0] *= -1;
                    }

                    if (Numeric::is_nan(fh->theta[0]))
                    {
                        fh->theta[0] = 0;
                    }

                }


                break;
            }

            eh = eh->next;
        } while (eh != fh->edge);
        
        if (!constrained)
        {
            fh->var_id[0] = x.size();
            x.push_back(fh->theta[0]);
            map_var_fh.push_back(fh);
        }
        else
        {
            q_fixed_jumps.push(FacetDist(NULL,fh,0));
            //fh->tag = true;
        }
    }

    p_mesh->reset_edges_tag(false);
    //reduce free period jump numbers.
    while (q_fixed_jumps.size()>0)
    {
        FacetDist fdist = q_fixed_jumps.top();
        q_fixed_jumps.pop();

        HE_face* fh = fdist.fh_;

        if (fh->tag)
        {
            continue;
        }

        fh->tag = true;

        if (fdist.coming_edge_)
        {
            fdist.coming_edge_->var_id = fdist.coming_edge_->pair->var_id = -1;
            fdist.coming_edge_->tag = fdist.coming_edge_->pair->tag = true; //mark zero jump for this edge.
            fdist.coming_edge_->pij = fdist.coming_edge_->pair->pij = 0;
        }

        HE_edge* eh = fh->edge;
        do
        {
            if (eh->pair->face)
            {
                HE_face* adj_fh = eh->pair->face;
                if (!adj_fh->tag)
                {
                    //adj_fh->tag = true;
                    q_fixed_jumps.push(FacetDist(eh, adj_fh, p_mesh->bc_dual_edge_length(eh) + fdist.dist_ ));

                    //eh->var_id = eh->pair->var_id = -1;
                    //eh->tag = eh->pair->tag = true; //mark zero jump for this edge.
                    //eh->pij = eh->pair->pij = 0;
                }
                
                if (fh->var_id[0] < 0 && adj_fh->var_id[0] < 0) //two facets are both constrained.
                {
                    eh->var_id = eh->pair->var_id = -1;
                    eh->tag = eh->pair->tag = true;
                    eh->pij = MY_ROUND(2*(adj_fh->theta - fh->theta - eh->kij)/M_PI);
                    eh->pair->pij = - eh->pij;
                }

            }
            eh = eh->next;
        } while (eh != fh->edge);
    }
    
    //setup free jump edges.
    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
    {
        HE_edge* eh = p_mesh->get_edge(eitr);

        eh->var_id = -1;

        if (!eh->tag)
        {
            eh->tag = eh->pair->tag = true;
            if (eh->face && eh->pair->face)
            {
                eh->var_id = x.size();
                x.push_back(eh->pij);
                map_var_eh.push_back(eh);
                ids_to_round.push_back(eh->var_id);
            }
        }
    }

    //setup matrices.

    // setup BtB
    gmm::col_matrix< gmm::wsvector< double > > BtB( x.size(), x.size() );

    // setup rhs
    std::vector< double > rhs(x.size(), 0);

    //fill matrices up.
    {
        p_mesh->reset_edges_tag(false);
        for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
        {
            HE_edge* eh = p_mesh->get_edge(eitr);
            
            if (eh->tag)
            {
                continue;
            }

            if (!eh->face || !eh->pair->face)
            {
                continue;
            }

            eh->tag = eh->pair->tag = true;

            //measure the curvature of cross field.
            double kij = eh->kij;
            double theta_i = eh->face->theta[0];
            double theta_j = eh->pair->face->theta[0];

            //double fij = theta_i + kij + M_PI*0.5*eh->pij - theta_j;
            double bij = ((eh->face->var_id[0] > -1)?0:theta_i) + kij 
                            + ((eh->var_id > -1 || eh->pair->var_id > -1)?0:(M_PI*0.5*eh->pij))
                            - ((eh->pair->face->var_id[0] > -1)?0:theta_j)
                            ;

            int varid [3] = {eh->face->var_id[0], eh->pair->face->var_id[0], eh->var_id};
            double diff [3] = {1, -1, M_PI*0.5};

            if (eh->var_id < 0 && eh->pair->var_id > -1)
            {
                varid[2] = eh->pair->var_id;
                diff[2] *= -1;
            }

            //fill up matrix.
            for (int var_itr = 0; var_itr < 3; var_itr++)
            {
                if (varid[var_itr] > -1)
                {
                    rhs[varid[var_itr]] -= diff[var_itr]*bij;//- 

                    for (int var_itr_1 = 0; var_itr_1 < 3; var_itr_1++)
                    {
                        if (varid[var_itr_1] > -1){
                            BtB(varid[var_itr], varid[var_itr_1]) += diff[var_itr]*diff[var_itr_1];		
                        }
                    }
                }
                
            }
        }
    }

    // BtB -> CSC
    gmm::csc_matrix<double> BtBCSC;
    BtBCSC.init_with_good_format( BtB);

    //std::cout << "                      the linear system now looks like..." << std::endl;
    //std::cout << "                      Matrix A\n " << BtBCSC << std::endl;
    //std::cout << "                      Right hand side b\n" << rhs << std::endl << std::endl;

    //std::cout << "---------- ---------- 3.4) solve the system using the mixed-integer solver..." << std::endl;
    // create solver
    COMISO::MISolver miso;
    //miso.set_direct_rounding(true);
    //miso.set_no_rounding(true);
    // miso solve
    miso.solve(BtBCSC, x, rhs, ids_to_round);
    //std::cout << "                      solution vector x is\n" << x << std::endl << std::endl;

    //read solution.
    for (int fitr = 0; fitr < map_var_fh.size(); fitr++)
    {
        HE_face* fh = map_var_fh[fitr];
        fh->theta[0] = x[fitr];
    }
    for (int eitr = 0; eitr < map_var_eh.size(); eitr++)
    {
        HE_edge* eh = map_var_eh[eitr];
        eh->pij = (int)x[eitr + map_var_fh.size()]; //
        eh->pair->pij = - eh->pij;
    }


    //compute vertex singularity index.
    for (int vitr = 0; vitr < p_mesh->get_num_of_vertices(); vitr++)
    {
        HE_vert* vh = p_mesh->get_vertex(vitr);
        double sum_angle = 0;
        double sum_kij = 0;
        double sum_pij = 0;

        if (!p_mesh->is_on_boundary(vh))
        {
            HE_edge* eh = vh->edge;

            do
            {
                HE_edge* prev_eh = eh->prev;

                sum_angle += acos(dot(normalize(eh->vert->pos-eh->prev->vert->pos),-normalize(prev_eh->vert->pos-prev_eh->prev->vert->pos)));

                sum_kij += eh->kij;
                sum_pij += eh->pij;

                eh = eh->pair->next;
            } while (eh != vh->edge);

            vh->singular_idx = (int) (4*((2*M_PI - sum_angle + sum_kij)/(2*M_PI) + sum_pij/4));

            //std::cout << vh->singular_idx << std::endl;
            
        }
    }

}


void Surfacing::adapt_surface_to_df(Mesh3D* p_mesh){

    mark_degenerate_facets(p_mesh);

    fit_curvature_tensor(p_mesh);

    smooth_curvature_variation(p_mesh);

    match_target_curvature(p_mesh);
}


void Surfacing::mark_degenerate_facets(Mesh3D* p_mesh){
    int valid_fnum = 0;
    for (int fitr = 0; fitr < p_mesh->get_num_of_faces(); fitr++)
    {
        HE_face* fh = p_mesh->get_face(fitr);
        fh->tag = true;

        if (fh->is_degenerate())
        {
            fh->tag = false;
        }
        else
        {
            valid_fnum ++;
        }
    }

    std::cout << "Valid facet num: " << valid_fnum << " out of " << p_mesh->get_num_of_faces() <<std::endl;
}


void Surfacing::fit_curvature_tensor(Mesh3D* p_mesh){

    //compute best fit curvature tensor per face.
    //std::cout << "recover curvature tensors..." << std::endl;
    for (int fitr = 0; fitr < p_mesh->get_num_of_faces(); fitr++)
    {
        HE_face* fh = p_mesh->get_face(fitr);

        if (!fh->tag)
        {
            continue;
        }

        //solve II such that II(ddir.u, ddir.v) = (ndiff.u, ndiff.v) + northo * (-ddir.v, ddir.u) for three edges.

        //initialize.
        fh->tensor_axis[0] = vec3(0,0,0);
        fh->tensor_axis[1] = vec3(0,0,0);
        fh->tensor_diag[0] = 0;
        fh->tensor_diag[1] = 0;

        vec3 faxis = normalize(fh->edge->vert->pos - fh->edge->prev->vert->pos);
        vec3 binormal = cross(fh->normal, faxis);

        vec3 u_axis(0,0,0), v_axis(0,0,0);

        {
            double cos_theta = cos(fh->theta[0]), sin_theta = sin(fh->theta[0]);

            u_axis = cos_theta*faxis + sin_theta*binormal;
            v_axis = -sin_theta*faxis + cos_theta*binormal;
        }

        std::vector<vec2> dx_list, dn_list;

        double lambda_1 = 0, lambda_2 = 0;

        facet_differential(p_mesh, fh, u_axis, v_axis, dx_list, dn_list);
        
        

        if (dx_list.size() >= 2)
        {
        
            //given cross field, the unkowns are at diagonal of tensor.

            TNT::Array2D<double> mat_A (dx_list.size(), 2, 0.0);
            TNT::Array1D<double> rhs_b (dx_list.size(), 0.0);

            for (int dx_itr = 0; dx_itr < dx_list.size(); dx_itr++)
            {
                int rid = dx_itr;
                vec2 u_dx = normalize(dx_list[dx_itr]);
                
                mat_A[rid][0] = u_dx.x * dx_list[dx_itr].x;
                mat_A[rid][1] = u_dx.y * dx_list[dx_itr].y;

                rhs_b[rid] = dot(u_dx, dn_list[dx_itr]);
            }

            TNT::Array2D<double> mat_AtA = trans_mult(mat_A, mat_A);
            TNT::Array1D<double> vec_Atb = trans_mult(mat_A, rhs_b);

            double mat_AtA_array [4];
            mat_AtA_array[0] = mat_AtA[0][0];
            mat_AtA_array[1] = mat_AtA[1][0];
            mat_AtA_array[2] = mat_AtA[0][1];
            mat_AtA_array[3] = mat_AtA[1][1];

            double mat_AtA_inv [4] = {0.0};

            double det = invert2x2(mat_AtA_array, mat_AtA_inv);
            if(abs(det) < 1e-15)
            {
                
            }
            else{
                lambda_1 = mat_AtA_inv[0]*vec_Atb[0] + mat_AtA_inv[2]*vec_Atb[1];
                lambda_2 = mat_AtA_inv[1]*vec_Atb[0] + mat_AtA_inv[3]*vec_Atb[1];
            }

            fh->tensor_axis[0] = u_axis;
            fh->tensor_axis[1] = v_axis;
            fh->tensor_diag[0] = lambda_1;
            fh->tensor_diag[1] = lambda_2;
        
        }
    }
}



void Surfacing::smooth_curvature_variation(Mesh3D* p_mesh){

    //smooth the curvature tensor for reducing curvature variation.
    //std::cout << "smooth curvature tensors..." << std::endl;

    if (curvature_smooth_weight_ == 0.0)
        return;


    double smooth_weight = curvature_smooth_weight_;

    double smooth_coeff = sqrt(smooth_weight);
    double preserve_coeff = (1.0-smooth_weight);

    //std::cout << "smooth_coeff: " << smooth_coeff << std::endl;

    const int face_num = p_mesh->get_num_of_faces();

    int var_num = 0;
    std::vector<MyPair<int>> map_face_varid (face_num);

    for (int fitr = 0; fitr < face_num; fitr++)
    {
        HE_face* fh = p_mesh->get_face(fitr);

        fh->id = fitr;
        {
            if (!fh->tag)
            {
                map_face_varid[fitr] = (MyPair<int>(-1,-1));
            }
            else
            {
                var_num += 2;
                map_face_varid[fitr] = (MyPair<int>(var_num-2,var_num-1));
            }
        }
    }

    std::vector<double> lambda_weights (var_num, 0);

    Sparse_Matrix normal_mat (var_num, var_num, SYM_LOWER);

    //the smooth part.
    p_mesh->reset_edges_tag(false);
    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
    {
        HE_edge* eh = p_mesh->get_edge(eitr);
        if (!eh->tag && eh->face && eh->pair->face && eh->face->tag && eh->pair->face->tag)
        {
            eh->tag = eh->pair->tag = true;

            double wij [2] = {0.0,0.0};

            HE_face* fh[2];
            fh[0] = eh->face;
            fh[1] = eh->pair->face;

            int adj_axis_map[2];
            adj_axis_map[0] = mod(eh->pij,2);
            adj_axis_map[1] = (adj_axis_map[0]+1)%2;

            vec3 dual_edge, ndiff;// = eh->ortho_dual;
            p_mesh->ortho_dual_edge_data(eh, dual_edge, ndiff, true);
            
            double dual_len = dual_edge.length();

            if (dual_len < 1e-5) //cocircular. 1e-7
            {
                wij[0] = wij[1] = 1.0; //give the maximum value.
            }
            else
            {
                //double scaling = 1.0/(dual_len*dual_len);//std::exp(-dual_len)/dual_len;
                //
                //dual_edge *= scaling;

                dual_edge = normalize(dual_edge);

                vec3 pair_dual_edge, pair_ndiff;// = eh->pair->ortho_dual * scaling;
                p_mesh->ortho_dual_edge_data(eh->pair, pair_dual_edge, pair_ndiff, true);
                //pair_dual_edge *= scaling;
                pair_dual_edge = normalize(pair_dual_edge);
                
                for (int axis_itr = 0; axis_itr < 2; axis_itr++)
                {
                    wij[axis_itr] = 0.5*( abs(dot(dual_edge,fh[0]->tensor_axis[axis_itr])) + abs(dot(pair_dual_edge,fh[1]->tensor_axis[adj_axis_map[axis_itr]])) );
                }
            }

            //fill two quadratic terms.
            for (int axis_itr = 0; axis_itr < 2; axis_itr++)
            {
                double coeff = smooth_coeff*wij[axis_itr];
                double diff[2] = {coeff, -coeff};

                int varid[2] = {map_face_varid[fh[0]->id].id[axis_itr], map_face_varid[fh[1]->id].id[adj_axis_map[axis_itr]]};
                
                double cur_val[2] = { fh[0]->tensor_diag[axis_itr], fh[1]->tensor_diag[adj_axis_map[axis_itr]] };

                //collect weights for preservence term.
                for (int var_itr = 0; var_itr < 2; var_itr++){
                    if (varid[var_itr] > -1) 
                    {
                        //std::cout << "lambda size: " << lambda_weights.size() << " index: " << varid[var_itr] << std::endl;
                        lambda_weights[varid[var_itr]] += wij[axis_itr];
                    }

                    ////debug
                    //if (varid[var_itr] >= lambda_weights.size())
                    //{
                    //	std::cout << "varid[2]: " << varid[0] << " " << varid[1] << std::endl;
                    //	std::cout << "adj_axis_map[2]: " << adj_axis_map[0] << " " << adj_axis_map[1] << std::endl;
                    //	int junk;
                    //	std::cin >> junk;
                    //}
                }
                
                

                //std::cout << eitr << ": " << cur_val[0] << " " << cur_val[1] << std::endl;

                for (int var_itr = 0; var_itr < 2; var_itr++)
                {
                    if (varid[var_itr] > -1)
                    {
                        for (int var_itr_1 = 0; var_itr_1 < 2; var_itr_1++)
                        {
                            if (varid[var_itr_1] > -1)
                            {
                                if (varid[var_itr_1] <= varid[var_itr])
                                {
                                    normal_mat.fill_entry(varid[var_itr], varid[var_itr_1], diff[var_itr]*diff[var_itr_1]);
                                }
                            }
                            else
                            {
                                normal_mat.fill_rhs_entry(varid[var_itr], -diff[var_itr]*diff[var_itr_1]*cur_val[var_itr_1]);
                            }
                        }
                        
                    }
                }
            }
        }
    }

    //the preservence part.
    for (int fitr = 0; fitr < face_num; fitr++)
    {
        MyPair<int>& var_ids = map_face_varid[fitr];

        if (p_mesh->is_on_boundary(p_mesh->get_face(fitr)))
        {
            continue;
        }

        for (int axis_itr = 0; axis_itr < 2; axis_itr++)
        {
            if (var_ids.id[axis_itr] > -1)
            {
                double lweight = lambda_weights[var_ids.id[axis_itr]] * 0.5;
                lweight *= lweight;

                normal_mat.fill_entry(var_ids.id[axis_itr], var_ids.id[axis_itr], preserve_coeff*lweight);
                normal_mat.fill_rhs_entry(var_ids.id[axis_itr], preserve_coeff*lweight*(p_mesh->get_face(fitr)->tensor_diag[axis_itr]));
            }
        }
    }


    //solve.
    solve_by_CHOLMOD(&normal_mat);

    //std::cout << normal_mat;

    //read solution.
    for (int fitr = 0; fitr < face_num; fitr++)
    {
        MyPair<int>& var_ids = map_face_varid[fitr];

        //if (p_mesh->is_on_boundary(p_mesh->get_face(fitr)))
        //{
        //	continue;
        //}

        for (int axis_itr = 0; axis_itr < 2; axis_itr++)
        {
            if (var_ids.id[axis_itr] > -1)
            {
                double smoothed_curv = normal_mat.get_solution()[var_ids.id[axis_itr]];
                if (!Numeric::is_nan(smoothed_curv))
                {
                    p_mesh->get_face(fitr)->tensor_diag[axis_itr] = smoothed_curv;
                }
            }
        }
    }
}

void Surfacing::match_target_curvature(Mesh3D* p_mesh){


    std::vector<HE_vert*> map_var_vtx;
    for (int vitr = 0; vitr < p_mesh->get_num_of_vertices(); vitr++)
    {
        HE_vert* vh = p_mesh->get_vertex(vitr);
        vh->var_id = -1;
        if (!p_mesh->is_on_boundary(vh) && !p_mesh->vertex_on_curve(vh))
        {
            map_var_vtx.push_back(vh);
            vh->var_id = map_var_vtx.size()-1;
        }
    }
    Sparse_Matrix Lvert (map_var_vtx.size()*3, map_var_vtx.size()*3, SYM_LOWER);

    //compute target normal differences per edge, according to curvature tensor, and build linear system for edge-based curvatures.
    //solve the system for updated vertices.
    //std::cout << "compute target dihedral angles..." << std::endl;
    p_mesh->reset_edges_tag(false);

    bool complained = false;//

    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
    {
        HE_edge* eh = p_mesh->get_edge(eitr);

        if (eh->tag || eh->pair->pair != eh)
        {
            continue;
        }

        eh->tag = eh->pair->tag = true;

        if (eh->face && eh->pair->face && eh->face->tag && eh->pair->face->tag)
        {
            
            vec2 dx0 ( dot(eh->ortho_dual, eh->face->tensor_axis[0]), dot(eh->ortho_dual, eh->face->tensor_axis[1]) );

            vec2 dn0 ( eh->face->tensor_diag[0]*dx0[0], eh->face->tensor_diag[1]*dx0[1] );

            vec2 dx0_dir (0,0);
            if (!Numeric::is_nan(dx0.length()) && dx0.length() > 1e-10)
            {
                dx0_dir = normalize(dx0);
            }

            double dn0_len = dot(dn0, dx0_dir);

            vec2 dx1 ( dot(eh->pair->ortho_dual, eh->pair->face->tensor_axis[0]), dot(eh->pair->ortho_dual, eh->pair->face->tensor_axis[1]) );

            vec2 dn1 ( eh->pair->face->tensor_diag[0]*dx1[0],	eh->pair->face->tensor_diag[1]*dx1[1] );

            vec2 dx1_dir (0,0);
            if (!Numeric::is_nan(dx1.length()) && dx1.length() > 1e-10)
            {
                dx1_dir = normalize(dx1);
            }

            double dn1_len = dot(dn1, dx1_dir);

            double sin_dihedral = (dn0_len + dn1_len)/2.0;
            //debug
            
            if (!complained && (Numeric::is_nan(sin_dihedral) || abs(sin_dihedral)>1))
            {
                std::cout << "WARNING: sin(theta) = "<<sin_dihedral<<std::endl;
                //sin_dihedral = sign(sin_dihedral);
                sin_dihedral = 0;

                complained = true;
            }
            
    
            /////////
            double sin_half_dihedral = sign(sin_dihedral)*sqrt( abs(0.5*( 1 - sqrt(abs(1-sin_dihedral*sin_dihedral)) )) );
            
            eh->target_sin_hdi = eh->pair->target_sin_hdi = sin_half_dihedral;

            //compute target mean curvature normal at edges and solve the linear equations.
            {
                vec3 curv_norm = normalize(eh->face->normal + eh->pair->face->normal);
                curv_norm *= sin_half_dihedral*2*(eh->vert->pos-eh->pair->vert->pos).length();

                if (eh->cycle_curv || eh->pair->cycle_curv)
                {
                    curv_norm = vec3(0,0,0);
                }

                edgecurv_fillmatrix(eh, curv_norm, Lvert);
                //std::cout << curv_norm << std::endl;
            }
            
        }
        
        if ( (!eh->face || !eh->pair->face) 
            && ((eh->cycle_curv && eh->cycle_curv->curve_type==CURVE_TYPE::CRVT_LINE) || 
            (eh->pair->cycle_curv && eh->pair->cycle_curv->curve_type==CURVE_TYPE::CRVT_LINE) )
            )
        {
            //the case of boundary curve, where RMF normal should be used.
            //std::cout << "RMF normal constraint"<<std::endl;

            fillmatrix_normalConstraint(eh, Lvert);
            fillmatrix_normalConstraint(eh->pair, Lvert);
        }
    }

    //std::cout << "evolve vertices..." << std::endl;

    //solve the system of equations and update vertices.
    solve_by_CHOLMOD(&Lvert);

    double avg_delta = 0;

    for (int var_itr = 0; var_itr < map_var_vtx.size(); var_itr++)
    {
        HE_vert* vh = map_var_vtx[var_itr];
        vec3 new_pos (Lvert.get_solution()[3*var_itr], Lvert.get_solution()[3*var_itr+1], Lvert.get_solution()[3*var_itr+2]);

        if (Numeric::is_nan(new_pos.x) || Numeric::is_nan(new_pos.y) || Numeric::is_nan(new_pos.z))
        {
            continue;
        }

        avg_delta += (new_pos - vh->pos).length();

        vh->pos = new_pos;
    }

    avg_delta /= map_var_vtx.size();

    vec3 bb_ll (p_mesh->xmin, p_mesh->ymin, p_mesh->zmin);
    vec3 bb_ur (p_mesh->xmax, p_mesh->ymax, p_mesh->zmax);
    double diagonal_len = (bb_ll-bb_ur).length();

    if (avg_delta/diagonal_len < 0.005)
    {
        std::cout << "Movement is " << avg_delta/diagonal_len*100 <<"% (< 0.5%), should stop iterations."<<std::endl;
    }
}



double Surfacing::dihedral_diff_gradient(HE_edge* edge, double target_sinhd, std::vector<HE_vert*>& vts_list, std::vector<vec3>& grad_list){

    vec3 current_curvnormal = cross(edge->vert->pos - edge->prev->vert->pos, edge->face->normal) 
                - cross(edge->vert->pos - edge->prev->vert->pos, edge->pair->face->normal);

    double cur_sign = sign(dot(current_curvnormal, edge->face->normal+edge->pair->face->normal));

    double cur_curvnormal_len = cur_sign * current_curvnormal.length();

    double edge_len = (edge->vert->pos - edge->prev->vert->pos).length();

    double f_cdiff = cur_curvnormal_len/edge_len - 2*target_sinhd;

    vts_list[0] = edge->next->vert;
    vts_list[1] = edge->prev->vert;
    vts_list[2] = edge->pair->next->vert;
    vts_list[3] = edge->vert;
    
    double cot013 = dot(vts_list[3]->pos-vts_list[1]->pos, vts_list[0]->pos-vts_list[1]->pos)/cross(vts_list[3]->pos-vts_list[1]->pos, vts_list[0]->pos-vts_list[1]->pos).length() ; 
    double cot031 = dot(vts_list[0]->pos-vts_list[3]->pos, vts_list[1]->pos-vts_list[3]->pos)/cross(vts_list[0]->pos-vts_list[3]->pos, vts_list[1]->pos-vts_list[3]->pos).length() ; 

    double cot213 = dot(vts_list[2]->pos - vts_list[1]->pos, vts_list[3]->pos - vts_list[1]->pos) / cross(vts_list[2]->pos - vts_list[1]->pos, vts_list[3]->pos - vts_list[1]->pos).length();
    double cot231 = dot(vts_list[2]->pos - vts_list[3]->pos, vts_list[1]->pos - vts_list[3]->pos) / cross(vts_list[2]->pos - vts_list[3]->pos, vts_list[1]->pos - vts_list[3]->pos).length();

    double diff[4];
    diff[0] = -cot013 - cot031;
    diff[1] = cot031 + cot231;
    diff[2] = -cot213 - cot231;
    diff[3] = cot213 + cot013;

    for (int i = 0; i < 4; i++)
    {
        diff[i] *= cur_sign;
    }

    vec3 unit_crvnorm (0,0,0);
    if (abs(cur_curvnormal_len) > 1e-8)
    {
        unit_crvnorm = normalize(current_curvnormal);
    }

    vec3 grad_He [4];
    vec3 grad_elen[4];

    for (int i = 0; i < 4; i++)
    {
        grad_He[i] = unit_crvnorm * diff[i];
    }

    grad_elen[0] = vec3(0,0,0);
    grad_elen[1] = normalize(vts_list[1]->pos - vts_list[3]->pos);
    grad_elen[2] = vec3(0,0,0);
    grad_elen[3] = normalize(vts_list[3]->pos - vts_list[1]->pos);

    for (int i = 0; i < 4; i++)
    {
        grad_list[i] = f_cdiff*(grad_He[i]*edge_len - cur_curvnormal_len*grad_elen[i])/(edge_len*edge_len);
    }

    return f_cdiff;
}


void Surfacing::edgecurv_fillmatrix(HE_edge* edge, vec3 targt_curvnorm, Sparse_Matrix& mat, double coeff){

    HE_vert* vts_list[4];

    vts_list[0] = edge->next->vert;
    vts_list[1] = edge->prev->vert;
    vts_list[2] = edge->pair->next->vert;
    vts_list[3] = edge->vert;
    
    double cot013 = dot(vts_list[3]->pos-vts_list[1]->pos, vts_list[0]->pos-vts_list[1]->pos)/cross(vts_list[3]->pos-vts_list[1]->pos, vts_list[0]->pos-vts_list[1]->pos).length() ; 
    double cot031 = dot(vts_list[0]->pos-vts_list[3]->pos, vts_list[1]->pos-vts_list[3]->pos)/cross(vts_list[0]->pos-vts_list[3]->pos, vts_list[1]->pos-vts_list[3]->pos).length() ; 

    double cot213 = dot(vts_list[2]->pos - vts_list[1]->pos, vts_list[3]->pos - vts_list[1]->pos) / cross(vts_list[2]->pos - vts_list[1]->pos, vts_list[3]->pos - vts_list[1]->pos).length();
    double cot231 = dot(vts_list[2]->pos - vts_list[3]->pos, vts_list[1]->pos - vts_list[3]->pos) / cross(vts_list[2]->pos - vts_list[3]->pos, vts_list[1]->pos - vts_list[3]->pos).length();
    
    ////debug.
    //vec3 cot_crvnorm = (vts_list[3]->pos - vts_list[0]->pos)*cot013 + (vts_list[1]->pos - vts_list[0]->pos)*cot031 + 
    //					(vts_list[3]->pos - vts_list[2]->pos)*cot213 + (vts_list[1]->pos - vts_list[2]->pos)*cot231;
    //std::cout <<"cur normal diff: " << (current_curvnormal - cot_crvnorm).length() << std::endl;

    double diff[4];
    diff[0] = -cot013 - cot031;
    diff[1] = cot031 + cot231;
    diff[2] = -cot213 - cot231;
    diff[3] = cot213 + cot013;

    //normalize edge length.
    //coeff *= 1.0/(vts_list[1]->pos - vts_list[3]->pos).length2();

    ////final check.

    //if (Numeric::is_nan(diff[0]) || Numeric::is_nan(diff[1]) || Numeric::is_nan(diff[2]) || Numeric::is_nan(diff[3]) || Numeric::is_nan(coeff))
    //{
    //	return;
    //}

    ////use cotangent weight?
    //double cot103 = dot(vts_list[1]->pos - vts_list[0]->pos, vts_list[3]->pos - vts_list[0]->pos) / cross(vts_list[1]->pos - vts_list[0]->pos, vts_list[3]->pos - vts_list[0]->pos).length();
    //double cot123 = dot(vts_list[1]->pos - vts_list[2]->pos, vts_list[3]->pos - vts_list[2]->pos) / cross(vts_list[1]->pos - vts_list[2]->pos, vts_list[3]->pos - vts_list[2]->pos).length();

    //double coeff = 0.5*(cot103 + cot123);
    //if (Numeric::is_nan(coeff))
    //{
    //	coeff = 0.0;
    //	return;
    //}
    ////coeff *= coeff;

    for (int coord_itr = 0; coord_itr < 3; coord_itr++)
    {
        for (int v_itr0 = 0; v_itr0 < 4; v_itr0++)
        {
            if (vts_list[v_itr0]->var_id > -1)
            {
                int ridx = vts_list[v_itr0]->var_id*3 + coord_itr;

                mat.fill_rhs_entry(ridx, diff[v_itr0]*targt_curvnorm[coord_itr]*coeff);

                for (int v_itr1 = 0; v_itr1 < 4; v_itr1++)
                {
                    if (vts_list[v_itr1]->var_id > -1)
                    {
                        int cidx = vts_list[v_itr1]->var_id*3 + coord_itr;
                        if (ridx >= cidx)
                        {
                            ////debug.
                            //if (Numeric::is_nan(diff[v_itr0]*diff[v_itr1]*coeff))
                            //{
                            //	std::cout << "Invalid numeric: " << diff[v_itr0] << " " << diff[v_itr1] << " " << coeff << std::endl;
                            //	char junk;
                            //	std::cin >> junk;
                            //}

                            mat.fill_entry(ridx, cidx, diff[v_itr0]*diff[v_itr1]*coeff);
                        }
                    }
                    else
                    {
                        mat.fill_rhs_entry(ridx, -diff[v_itr0]*diff[v_itr1]*vts_list[v_itr1]->pos[coord_itr]*coeff);
                    }
                }
            }

        }
        
    }
}

void Surfacing::fillmatrix_normalConstraint(HE_edge* edge, Sparse_Matrix& mat){
    if (edge && edge->face && edge->cycle_curv && edge->cycle_curv->curve_type==CURVE_TYPE::CRVT_LINE)
    {
        HE_vert* vh = edge->next->vert;
        double coeff = 1.0;

        if (vh->var_id < 0){
            return;
        }

        vec3 normal = edge->tangent;
        vec3 rhs = 2.0 * normal * (dot(normal,edge->vert->pos + edge->pair->vert->pos));
        vec3 lhs [3] = {4.0 * normal.x * normal, 4.0 * normal.y * normal, 4.0 * normal.z * normal};

        int mat_idx [3] = {vh->var_id*3, vh->var_id*3+1, vh->var_id*3+2};

        
        for (int coord_itr = 0; coord_itr < 3; coord_itr++)
        {
            int rid = mat_idx[coord_itr];
            mat.fill_rhs_entry(rid, rhs[coord_itr]*coeff);
            for (int coord_itr_inner = 0; coord_itr_inner < 3; coord_itr_inner++)
            {
                int cid = mat_idx[coord_itr_inner];
                if (rid >= cid)
                {
                    mat.fill_entry(rid, cid, lhs[coord_itr][coord_itr_inner]*coeff);
                }
            }
        }			
    }
}


void Surfacing::facet_differential(const Mesh3D* pmesh, const HE_face* face, 
    const vec3 u_axis, const vec3 v_axis, std::vector<vec2>& dx_list, std::vector<vec2>& dn_list, const HE_face* prev_face){

    HE_edge* eh = face->edge;
    do
    {
        if (eh->pair->face && eh->pair->face != prev_face)
        {
            vec3 ddir;
            vec3 ndiff;

            double cot_weight = pmesh->ortho_dual_edge_data(eh, ddir, ndiff, true);
            eh->ortho_dual = ddir;

            vec2 dx = vec2(dot(ddir,u_axis), dot(ddir,v_axis));
            vec2 dn = vec2(dot(ndiff,u_axis), dot(ndiff,v_axis));

            dx_list.push_back(dx);
            dn_list.push_back(dn);
        }

        eh = eh->next;
    } while (eh != face->edge);
}