

#include "Mesh3D.h"





//////////////////////////////////////////////////////////////////////////

Mesh3D::Mesh3D(void)
{
	//initialization
	vertices_list = NULL;
	faces_list = NULL;
	edges_list = NULL;

	xmax = ymax = zmax = (double)1.0;
	xmin = ymin = zmin = (double)-1.0;

	m_closed = true;
	m_quad = false;
	m_tri = false;
	m_hex = false;
	m_pentagon = false;

	m_num_components = 0;
	m_num_boundaries = 0;
	m_genus = 0;
	m_encounter_non_manifold = false;

	m_kdtree_=NULL;
}
//////////////////////////////////////////////////////////////////////////

Mesh3D::~Mesh3D(void)
{
	clear_data();

	if (m_kdtree_)
	{
		delete m_kdtree_;
	}
}
//////////////////////////////////////////////////////////////////////////
void Mesh3D::clear_data()
{
	clear_vertices();
	clear_edges();
	clear_faces();
}
//////////////////////////////////////////////////////////////////////////
void Mesh3D::clear_vertices()
{
	if (vertices_list==NULL) 
		return;

	for (VERTEX_ITER viter = vertices_list->begin(); viter != vertices_list->end(); viter++) 
		delete *viter;

	delete vertices_list;
	vertices_list = NULL;
}
//////////////////////////////////////////////////////////////////////////
void Mesh3D::clear_edges()
{
	m_edgemap.clear();
	if (edges_list==NULL) 
		return;

	for (EDGE_ITER eiter = edges_list->begin(); eiter != edges_list->end(); eiter++) 
		delete *eiter;

	delete edges_list;
	edges_list = NULL;
}
//////////////////////////////////////////////////////////////////////////
void Mesh3D::clear_faces()
{
	if (faces_list==NULL) 
		return;

	for (FACE_ITER fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++) 
		delete *fiter;

	delete faces_list;
	faces_list = NULL;
}
//////////////////////////////////////////////////////////////////////////
HE_vert* Mesh3D::insert_vertex(vec3 &v)
{
	HE_vert* hv = new HE_vert(v);
	if (vertices_list==NULL) 
		vertices_list = new VERTEX_LIST;

	hv->id = (ptrdiff_t)vertices_list->size();
	vertices_list->push_back(hv);

	return hv;
}
//////////////////////////////////////////////////////////////////////////
HE_face* Mesh3D::insert_face(VERTEX_LIST & vec_hv, std::vector<ptrdiff_t> * texture, std::vector<ptrdiff_t> * normal)
{
	int vsize = (int)vec_hv.size();
	if (vsize < 3) 
		return NULL;

	//////////////////////////////////////////////////////////////////////////
	//detect non-manifold
	bool b_find = false;
	for (int i = 0; i < vsize; i++)
	{
		HE_edge* c_he = m_edgemap[PAIR_VERTEX(vec_hv[i], vec_hv[(i+1)%vsize])];
		if (c_he && c_he->face)
		{
			//detect nonmanifold
			
			b_find = true;
			break;
		}
	}

	if (b_find)
	{
		m_encounter_non_manifold = true;

		//guess there are faces with reverse orientation, try to reverse them.
		//If these faces are inserted after their neighbor faces which have correct orientation,
		//probably we can obtain correct meshes.
		for (int i = 0; i < vsize; i++)
		{
			HE_edge* c_he = m_edgemap[PAIR_VERTEX(vec_hv[(i+1)%vsize], vec_hv[i])];
			if (c_he && c_he->face)
			{
				return NULL; 
			}
		}
		//it seems fine, reverse the vector.
		std::reverse(vec_hv.begin(), vec_hv.end());
		if (texture)
		{
			std::reverse(texture->begin(), texture->end());
		}
		if (normal)
		{
			std::reverse(normal->begin(), normal->end());
		}
	}

	//////////////////////////////////////////////////////////////////////////

	if (faces_list==NULL) 
		faces_list = new FACE_LIST;

	HE_face* hf = new HE_face;
	hf->valence = vsize;
	VERTEX_ITER viter = vec_hv.begin();
	VERTEX_ITER nviter = vec_hv.begin();
	nviter++;

	HE_edge *he1, *he2;
	std::vector<HE_edge* > v_he;
	int i;
	for (i=0; i<vsize-1; i++) 
	{
		he1 = insert_edge( *viter, *nviter);
		he2 = insert_edge( *nviter, *viter);

		if (hf->edge==NULL) 
			hf->edge = he1;

		he1->face = hf;
		he1->pair = he2;
		he2->pair = he1;
		v_he.push_back(he1);
		viter++, nviter++;
	}

	nviter = vec_hv.begin();

	he1 = insert_edge(*viter, *nviter);
	he2 = insert_edge(*nviter , *viter);
	he1->face = hf;
	if (hf->edge==NULL) 
		hf->edge = he1;

	he1->pair = he2;
	he2->pair = he1;
	v_he.push_back(he1);

	for (i=0; i<vsize-1; i++) 
	{
		v_he[i]->next = v_he[i+1];
		v_he[i+1]->prev = v_he[i];
	}
	v_he[i]->next = v_he[0];
	v_he[0]->prev = v_he[i];

	hf->id = (int)faces_list->size();
	faces_list->push_back(hf);


	if (texture)
	{
		hf->texture_indices.resize(texture->size());
		std::copy(texture->begin(), texture->end(), hf->texture_indices.begin());
	}
	if (normal)
	{
		hf->normal_indices.resize(normal->size());
		std::copy(normal->begin(), normal->end(), hf->normal_indices.begin());
	}

	return hf;
}
//////////////////////////////////////////////////////////////////////////

HE_edge* Mesh3D::insert_edge(HE_vert* vstart, HE_vert* vend)
{
	if (vstart == NULL || vend == NULL)
	{
		return NULL;
	}

	if (edges_list==NULL) 
		edges_list = new EDGE_LIST;

	if( m_edgemap[PAIR_VERTEX(vstart, vend)] != NULL )
		return m_edgemap[PAIR_VERTEX(vstart, vend)];

	HE_edge* he = new HE_edge;

	he->vert = vend;
	he->vert->degree++;
	vstart->edge = he;
	m_edgemap[PAIR_VERTEX(vstart, vend)] = he;

	he->id = (int)edges_list->size();
	edges_list->push_back(he);

	return he;
}

/////////////////////////////////////////////////////////////////////////
HE_edge* Mesh3D::pg_insert_edge(HE_vert* vstart, HE_vert* vend){
	if (vstart == NULL || vend == NULL)
	{
		return NULL;
	}

	if (edges_list==NULL) 
		edges_list = new EDGE_LIST;

	//if( m_edgemap[PAIR_VERTEX(vstart, vend)] != NULL )
	//	return m_edgemap[PAIR_VERTEX(vstart, vend)];

	HE_edge* he = new HE_edge;

	he->vert = vend;
	he->vert->degree++;
	vstart->pg_edges.push_back(he); //hpan.
	vstart->edge = he;
	//m_edgemap[PAIR_VERTEX(vstart, vend)] = he;

	he->id = (int)edges_list->size();
	edges_list->push_back(he);

	return he;
}

//////////////////////////////////////////////////////////////////////////
HE_face* Mesh3D::pg_insert_face(std::vector<HE_edge*>& edges){
	if (faces_list==NULL) 
		faces_list = new FACE_LIST;

	HE_face* hf = new HE_face;
	hf->valence = edges.size();

	for(int i=0; i<edges.size(); i++){
		edges[i]->face = hf;

		HE_edge* eh0 = edges[i], *eh1 = edges[(i+1)%edges.size()];
		if(eh0->vert == eh1->pair->vert)
		{
			eh0->next = eh1;
			eh1->prev = eh0;
		}
		else if (eh0->pair->vert == eh1->vert)
		{
			eh0->prev = eh1;
			eh1->next = eh0;
		}
		else{
			std::cerr<<"ERROR pg_insert_face()."<<std::endl;
		}
	}

	hf->edge = edges[0];
	hf->id = (int)faces_list->size();
	faces_list->push_back(hf);
	compute_perface_normal(hf);

	return hf;
}

/////////////////////////////////////////////////////////////////////////
void Mesh3D::init_vef_list(){
	if(!vertices_list) vertices_list = new VERTEX_LIST;
	if(!edges_list) edges_list = new EDGE_LIST;
	if(!faces_list) faces_list = new FACE_LIST;
}

//////////////////////////////////////////////////////////////////////////

void Mesh3D::set_nextedge_for_border_vertices()
{
	if (!is_valid()) 
		return;

	EDGE_ITER eiter = edges_list->begin();

	for (; eiter!=edges_list->end(); eiter++) 
	{
		if ((*eiter)->next==NULL&& (*eiter)->face==NULL)
			(*eiter)->pair->vert->edge = *eiter;

	}

	for (eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++) 
	{
		if ( (*eiter)->next ==NULL ) 
		{
			HE_vert* hv = (*eiter)->vert;
			if (hv->edge != (*eiter)->pair)
			{
				(*eiter)->next = hv->edge;
				hv->edge->prev = *eiter;
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////

bool Mesh3D::is_on_boundary(HE_vert* hv)
{
	HE_edge* edge = hv->edge; 
	do {
		if (edge==NULL||edge->pair->face==NULL||edge->face==NULL) 
			return true;

		edge = edge->pair->next; 

	} while (edge != hv->edge);

	return false;
}
//////////////////////////////////////////////////////////////////////////

bool Mesh3D::is_on_boundary(HE_face* hf)
{
	HE_edge* edge = hf->edge; 

	do {
		if (is_on_boundary(edge)) 
			return true;
		edge = edge->next;
	} while (edge != hf->edge);

	return false;
}
//////////////////////////////////////////////////////////////////////////

bool Mesh3D::is_on_boundary(HE_edge* he)
{
	if(he->face==NULL||he->pair->face==NULL)
		return true;
	return false;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::check_closed()
{
	m_closed = true; 
	VERTEX_ITER viter = vertices_list->begin();
	for (;viter!=vertices_list->end(); viter++) 
	{
		if (is_on_boundary(*viter)) {
			m_closed = false;
			return;
		}
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::check_meshtype()
{
	if (get_num_of_faces() == 0)
	{
		m_tri = false;
		m_quad = false;
		m_hex = false;
		return;
	}

	FACE_ITER fiter = faces_list->begin();
	m_tri = true;
	m_quad = true;
	m_hex = true;
	m_pentagon = true;
	for (;fiter!=faces_list->end(); fiter++)
	{
		int d = (*fiter)->valence;
		if (d != 3)
		{
			m_tri = false;
		}
		if (d != 4)
		{
			m_quad = false;
		}
		if (d != 5)
		{
			m_pentagon = false;
		}
		if (d !=6)
		{
			m_hex = false;
		}
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::update_mesh()
{
	if (!is_valid()) 
		return;

	remove_hanged_vertices(); 
	set_nextedge_for_border_vertices();
	check_closed();
	check_meshtype();
	m_edgemap.clear();
	update_normal();
	compute_boundingbox();
	compute_genus();

}
//////////////////////////////////////////////////////////////////////////
//FILE IO
//////////////////////////////////////////////////////////////////////////

bool Mesh3D::load_off(const char* fins)
{
	std::ifstream fin(fins);

	try{

		clear_data();

		int vsize, fsize, esize;

		std::string head;

		fin>>head;

		if (head == "OFF") {
		}
		else
			return false;

		fin>>vsize>>fsize>>esize;

		double x, y, z;
		for (int i = 0; i<vsize; i++) {
			fin>>x>>y>>z;
			vec3 nvv(x, y, z) ;
			insert_vertex(nvv);
		}

		for (int i = 0; i<fsize; i++) {

			VERTEX_LIST v_list;
			int valence;
			fin>>valence;

			for (int j=0; j<valence; j++) {
				int id;
				fin>>id;
				HE_vert* hv = get_vertex(id);

				bool findit = false;
				for (int i = 0; i <(int) v_list.size(); i++)
				{
					if (hv == v_list[i])
					{
						findit = true;
						break;
					}
				}
				if (findit == false && hv != NULL)
				{
					v_list.push_back(hv);
				}

			}

			if ((int)v_list.size() >= 3)
			{
				insert_face(v_list);
			}
		}

		update_mesh();
	}
	catch (...) {
		//catch any error
		clear_data();
		xmax = ymax = zmax = (double)1.0;
		xmin = ymin = zmin = (double)-1.0;
		fin.close();
		return false;
	}
	fin.close();
	return is_valid();
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::write_off(const char* fouts)
{

	std::ofstream fout(fouts);
	fout.precision(16);
	fout<<"OFF\n";
	//output the number of vertices_list, faces_list-> edges_list
	fout<<(int)vertices_list->size()<<" "<<(int)faces_list->size()<<" "<<(int)edges_list->size()/2<<"\n";

	//output coordinates of each vertex
	VERTEX_ITER viter = vertices_list->begin();
	for (;viter!=vertices_list->end(); viter++) {

		fout<<std::scientific<<(*viter)->pos<<"\n";
	}

	//output the valence of each face and its vertices_list' id

	FACE_ITER fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) 
	{
		fout<<(*fiter)->valence;

		HE_edge* edge = (*fiter)->edge; 

		do {
			fout<<" "<<edge->pair->vert->id;
			edge = edge->next;

		} while (edge != (*fiter)->edge);
		fout<<"\n";
	}

	fout.close();
}
//////////////////////////////////////////////////////////////////////////

bool Mesh3D::load_obj(const char* fins)
{
	FILE *m_pFile = fopen(fins, "r");

	if (!m_pFile)
	{
		return false;
	}

	char *tok;
	char temp[128];

	try
	{
		clear_data();
		//read vertices
		fseek(m_pFile, 0, SEEK_SET);
		char pLine[512];
		while(fgets(pLine, 512, m_pFile))
		{
			if(pLine[0] == 'v' && pLine[1] == ' ')
			{
				vec3 nvv;
				tok = strtok(pLine," ");
				for (int i=0; i<3; i++) 
				{
					tok = strtok(NULL," ");
					strcpy(temp, tok);
					temp[strcspn(temp," ")] = 0;
					nvv[i] = (double)atof(temp);
				}
				insert_vertex(vec3(nvv.x,nvv.y,nvv.z));
			}
			else if(pLine[0] == 'v' && pLine[1] == 't')
			{
// 					vec3g<float> nvv;
// 					tok = strtok(pLine," ");
// 					for (int i=0; i<2; i++) 
// 					{
// 						tok = strtok(NULL," ");
// 						strcpy(temp, tok);
// 						temp[strcspn(temp," ")] = 0;
// 						nvv[i] = (float)atof(temp);
// 					}
// 					texture_array.push_back(nvv);
			}
			else if(pLine[0] == 'v' && pLine[1] == 'n')
			{
// 					TinyVector<float, 3> nvv;
// 					tok = strtok(pLine," ");
// 					for (int i=0; i<3; i++) 
// 					{
// 						tok = strtok(NULL," ");
// 						strcpy(temp, tok);
// 						temp[strcspn(temp," ")] = 0;
// 						nvv[i] = (float)atof(temp);
// 					}
// 					normal_array.push_back(nvv);
			}
		}

		//read facets

		fseek(m_pFile, 0, SEEK_SET);


		while(fgets(pLine, 512, m_pFile))
		{
			char *pTmp = pLine;
			if(pTmp[0] == 'f')
			{
				VERTEX_LIST s_faceid;
				std::vector<ptrdiff_t> normal_ind, texture_ind;
				tok = strtok(pLine," ");
				while ((tok = strtok(NULL," ")) != NULL)
				{
					strcpy(temp, tok);
					size_t len = strlen(temp);
					int start_pos = 0;
					int pos[2] = {-1,-1};
					for (int k = 0; k < (int)len; k++)
					{
						if (temp[k] == '/')
						{
							pos[start_pos] = k;
							start_pos++;
						}
					}
					size_t end_pos = len;
					if (start_pos == 1)
					{
						end_pos = pos[0];
						std::string mstr(&temp[pos[0]+1], len-pos[0]-1);
						int id = (int)strtol(mstr.c_str(), NULL, 10) - 1;
						texture_ind.push_back(id);
					}
					else if (start_pos == 2)
					{
						end_pos = pos[0];
						if (pos[0]+1!=pos[1])
						{
							std::string mstr(&temp[pos[0]+1], pos[1]-pos[0]-1);
							int id = (int)strtol(mstr.c_str(), NULL, 10) - 1;
							texture_ind.push_back(id);
						}
						std::string mstr2(&temp[pos[1]+1], len-pos[1]-1);
						int id2 = (int)strtol(mstr2.c_str(), NULL, 10) - 1;
						normal_ind.push_back(id2);

					}
					std::string mstr(&temp[0], end_pos);
					int id = (int)strtol(mstr.c_str(), NULL, 10) - 1;
					HE_vert* hv = get_vertex(id);
					bool findit = false;
					for (int i = 0; i <(int) s_faceid.size(); i++)
					{
						if (hv == s_faceid[i])	//remove the redundant vertex id if it exists
						{
							findit = true;
							break;
						}
					}
					if (findit == false && hv != NULL)
					{
						s_faceid.push_back(hv);
					}
				}
				if (s_faceid.size() >= 3)
				{
					insert_face(s_faceid, texture_ind.size()==s_faceid.size()?&texture_ind:0, normal_ind.size() == s_faceid.size()?&normal_ind:0);
				}
			}
		}
		update_mesh();
	}
	catch (...)
	{
		clear_data();
		xmax = ymax = zmax = (double)1.0;
		xmin = ymin = zmin = (double)-1.0;

		fclose(m_pFile);
		return false;
	}

	fclose(m_pFile);
	return is_valid();
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::write_obj(const char* fouts)
{
	std::ofstream fout(fouts);

	fout<<"g object\n";
	fout.precision(16);
	//output coordinates of each vertex
	VERTEX_ITER viter = vertices_list->begin();
	for (;viter!=vertices_list->end(); viter++) {

		fout<<"v "<< std::scientific <<(*viter)->pos<<"\n";
	}


	FACE_ITER fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) {

		fout<<"f";

		HE_edge* edge = (*fiter)->edge; 
		unsigned int count = 0;
		do {
			fout<<" "<<edge->pair->vert->id+1;
			if ((*fiter)->has_texture_map())
			{
				fout<<"/" <<(*fiter)->texture_indices[count]+1;
			}
			else
				fout <<"/";

			if ((*fiter)->has_normal_map())
			{
				fout<<"/" <<(*fiter)->normal_indices[count]+1;
			}
			else
				fout<<"/" <<edge->pair->vert->id+1;


			edge = edge->next;

		} while (edge != (*fiter)->edge);
		fout<<"\n";
	}
	fout.close();
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::write_vtk(const char* fouts)
{
	std::ofstream fout(fouts);
	fout.precision(16);
	fout<<"# vtk DataFile Version 3.0\n"
		<<"Polygonal Mesh\n"
		<<"ASCII\n"
		<<"DATASET POLYDATA\n"
		<<"POINTS " << vertices_list->size() <<" double\n";

	//output coordinates of each vertex
	VERTEX_ITER viter = vertices_list->begin();
	for (;viter!=vertices_list->end(); viter++) {
		fout<<std::scientific<<(*viter)->pos<<"\n";
	}

	FACE_ITER fiter = faces_list->begin();
	size_t count = 0;
	for (;fiter!=faces_list->end(); fiter++) 
	{
		count += (*fiter)->valence + 1;
	}

	fout <<"POLYGONS " << faces_list->size() << " " << count << "\n";

	//output the valence of each face and its vertices_list' id

	fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) 
	{
		fout<<(*fiter)->valence;

		HE_edge* edge = (*fiter)->edge; 

		do {
			fout<<" "<<edge->pair->vert->id;
			edge = edge->next;

		} while (edge != (*fiter)->edge);
		fout<<"\n";
	}

	fout.close();
}
//////////////////////////////////////////////////////////////////////////
//For rendering

void Mesh3D::compute_boundingbox()
{
	if (vertices_list->size()<3) 
		return;

	xmax = ymax = zmax = (double)-10e10;
	xmin = ymin = zmin = (double)10e10;

	//zmin = zmax = 0;

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter!=vertices_list->end(); viter++) 
	{
		xmin = (*viter)->pos[0]<xmin? (*viter)->pos[0]:xmin;
		ymin = (*viter)->pos[1]<ymin? (*viter)->pos[1]:ymin;
		zmin = (*viter)->pos[2]<zmin? (*viter)->pos[2]:zmin;
		xmax = (*viter)->pos[0]>xmax? (*viter)->pos[0]:xmax;
		ymax = (*viter)->pos[1]>ymax? (*viter)->pos[1]:ymax;
		zmax = (*viter)->pos[2]>zmax? (*viter)->pos[2]:zmax;
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_faces_list_normal()
{
	for (FACE_ITER fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++) 
	{
		compute_perface_normal(*fiter);
	}
}
////////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_perface_normal(HE_face* hf)
{
	size_t i = 0;
	HE_edge* pedge = hf->edge;
	HE_edge* nedge = hf->edge->next;

	hf->normal = vec3(0, 0, 0);
	for (i=0; i<hf->valence; i++) 
	{
		//cross product
		HE_vert* p = pedge->vert;
		HE_vert* c = pedge->next->vert;
		HE_vert* n = nedge->next->vert;
		vec3 pc, cn;
		pc = c->pos - p->pos; 
		cn = n->pos - c->pos; 

		hf->normal += cross(pc,cn);
		pedge = nedge;
		nedge = nedge->next;
		if (hf->valence==3) 
		{
			break;
		}
	}

	//if (hf->normal.length() > 1e-7)
	{
		hf->normal = normalize(hf->normal);
	}
	//else
	//{
	//	hf->normal = vec3(0,0,0);
	//}
	
}
////////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_vertices_list_normal()
{
	VERTEX_ITER viter = vertices_list->begin();

	for (; viter!=vertices_list->end(); viter++) 
	{
		compute_pervertex_normal(*viter);
	}
}
////////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_pervertex_normal(HE_vert* hv)
{

	HE_edge* edge = hv->edge; 

	if (edge==NULL) {
		hv->normal = vec3(0, 0, 0);
		return;
	}
	hv->normal = vec3(0, 0, 0);

	do {
		if (edge->face!=NULL)
		{
			if (edge->face->get_area() > 1e-10) //check degenerate facets. hpan.
			{
				hv->normal += edge->face->normal;
			}
			
		}
		edge = edge->pair->next; 

	} while (edge && edge != hv->edge);

	hv->normal = normalize(hv->normal);
}
////////////////////////////////////////////////////////////////////////////

void Mesh3D::update_normal(bool onlyupdate_facenormal)
{
	compute_faces_list_normal();
	if (onlyupdate_facenormal == false)
	{
		compute_vertices_list_normal();
	}

}
////////////////////////////////////////////////////////////////////////////

void Mesh3D::copy_edge_id(ptrdiff_t edgeid, ptrdiff_t startid, ptrdiff_t endid, Mesh3D* mesh)
{
	mesh->m_edgemap[PAIR_VERTEX( mesh->get_vertex(startid), mesh->get_vertex(endid))]->id = edgeid;
}

////////////////////////////////////////////////////////////////////////////

Mesh3D* Mesh3D::make_copy_simple(){
	Mesh3D* new_mesh = new Mesh3D;

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter!=vertices_list->end(); viter++)
	{
		HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

		//hpan
		if ((*viter)->fixed)
		{
			vh->fixed = true;
		}
		vh->m_vtype = (*viter)->m_vtype;
	}

	FACE_ITER fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) 
	{
		HE_face* hf = *fiter;
		HE_edge* edge = hf->edge;
		VERTEX_LIST mvlist;
		do {
			mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
			edge = edge->next;
		} while(edge!=hf->edge);

		HE_face* rfh = new_mesh->insert_face(mvlist, hf->texture_indices.empty()?0:&hf->texture_indices, hf->normal_indices.empty()?0:&hf->normal_indices);
	}

	new_mesh->update_mesh();

	//hpan. Copy pointer to curves.
	PTR_EDGE_LIST edges_list = new_mesh->get_edges_list();
	for (int i = 0; i < edges_list->size(); i++)
	{
		//(*edgelist)[i]->spln_curv = (*edges_list)[i]->spln_curv;
		(*edges_list)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
	}


	return new_mesh;
}

//////////////////////////////////////////////////////////////////////////

Mesh3D* Mesh3D::make_copy()
{
	Mesh3D* new_mesh = new Mesh3D;

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter!=vertices_list->end(); viter++)
	{
		HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

		//hpan
		if ((*viter)->fixed)
		{
			vh->fixed = true;
		}
		vh->m_vtype = (*viter)->m_vtype;
	}

	FACE_ITER fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) 
	{
		HE_face* hf = *fiter;
		HE_edge* edge = hf->edge;
		VERTEX_LIST mvlist;
		do {
			mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
			edge = edge->next;
		} while(edge!=hf->edge);

		HE_face* rfh = new_mesh->insert_face(mvlist, hf->texture_indices.empty()?0:&hf->texture_indices, hf->normal_indices.empty()?0:&hf->normal_indices);

		rfh->curve_cycle = hf->curve_cycle;
	}

	for (EDGE_ITER eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++)
	{
		HE_edge* he = *eiter;
		copy_edge_id(he->id, he->pair->vert->id, he->vert->id, new_mesh);
	}

	PTR_EDGE_LIST edgelist = new_mesh->get_edges_list();

	std::sort(edgelist->begin(), edgelist->end(), CompareEdgeID);

	new_mesh->update_mesh();

	VERTEX_ITER cviter = new_mesh->get_vertices_list()->begin();
	for (viter = vertices_list->begin(); viter!=vertices_list->end(); viter++, cviter++)
	{
		if ((*viter)->edge != NULL)
		{
			(*cviter)->edge = new_mesh->get_edge((*viter)->edge->id);
		}
	}

	FACE_ITER cfiter = new_mesh->get_faces_list()->begin();
	for (fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++ , cfiter++)
	{
		(*cfiter)->edge = new_mesh->get_edge((*fiter)->edge->id);
	}

	//hpan. Copy pointer to curves.
	for (int i = 0; i < edgelist->size(); i++)
	{
		//(*edgelist)[i]->spln_curv = (*edges_list)[i]->spln_curv;
		(*edgelist)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
	}



	return new_mesh;
}
//////////////////////////////////////////////////////////////////////////

Mesh3D* Mesh3D::reverse_orientation()
{
	Mesh3D* new_mesh = new Mesh3D;

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter!=vertices_list->end(); viter++)
	{
		new_mesh->insert_vertex((*viter)->pos);
	}

	FACE_ITER fiter = faces_list->begin();

	for (;fiter!=faces_list->end(); fiter++) 
	{
		HE_face* hf = *fiter;
		HE_edge* edge = hf->edge;
		VERTEX_LIST mvlist;
		do {
			mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
			edge = edge->next;
		} while(edge!=hf->edge);
		std::reverse(mvlist.begin(), mvlist.end());
		new_mesh->insert_face(mvlist);
	}

	new_mesh->update_mesh();
	return new_mesh;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_num_components()
{
	FACE_ITER fiter = faces_list->begin();
	for (; fiter != faces_list->end(); fiter++ )
	{
		(*fiter)->tag = false;
	}

	m_num_components = 0;

	for (fiter = faces_list->begin(); fiter != faces_list->end(); fiter++ )
	{
		HE_face* hf = *fiter;
		if (hf->tag == false)
		{
			m_num_components++;

			std::queue<HE_face* > facets;
			facets.push(hf);
			hf->tag = true;

			while (!facets.empty())
			{
				HE_face* pFacet = facets.front();
				facets.pop();
				pFacet->tag = true;

				HE_edge* he = pFacet->edge;
				do 
				{
					if(he->pair->face!=NULL && he->pair->face->tag==false)
					{
						facets.push(he->pair->face);
						he->pair->face->tag = true;

						HE_edge* vhe = he->vert->edge;

						do 
						{
							if (vhe->face && vhe->face->tag == false)
							{
								facets.push(vhe->face);
								vhe->face->tag = true;

								HE_edge* mvhe = vhe->vert->edge;
								do 
								{
									if (mvhe->face && mvhe->face->tag == false)
									{
										facets.push(mvhe->face);
										mvhe->face->tag = true;
									}
									mvhe = mvhe->pair->next;
								} while(mvhe != vhe->vert->edge);

							}
							vhe = vhe->pair->next;
						} while(vhe != he->vert->edge);

					}
					he = he->next;
				} while(he != pFacet->edge);

				he = pFacet->edge;
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_num_boundaries()
{
	//if the mesh is not manifold, this function may report the wrong number

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter != vertices_list->end(); viter++ )
	{
		(*viter)->tag = false;
	}

	boundaryvertices.clear();

	m_num_boundaries = 0;

	for(viter = vertices_list->begin(); viter != vertices_list->end(); viter++)
	{

		HE_vert* hv = *viter;

		if (is_on_boundary(hv) && hv->tag == false)
		{

			std::vector<HE_vert*> vedges;
			vedges.push_back(hv);

			m_num_boundaries++;
			std::queue<HE_vert* > vertices;
			hv->tag = true;
			vertices.push(hv);

			while (!vertices.empty())
			{
				HE_vert* pVertex = vertices.front();
				pVertex->tag = true;
				vertices.pop();

				HE_edge* he = pVertex->edge;

				do 
				{
					if (is_on_boundary(he) && is_on_boundary(he->vert) && he->vert->tag==false)
					{
						he->vert->tag = true;
						vertices.push(he->vert);
						vedges.push_back(he->vert);

					}
					he = he->pair->next;
				} while(he!=pVertex->edge);

			}

			HE_vert* hvtmp = vedges[0];
			vedges[0] = vedges[1];
			vedges[1] = hvtmp;
			boundaryvertices.push_back(vedges);
		}
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::compute_genus()
{
	compute_num_components();
	compute_num_boundaries();
	int c = m_num_components;
	int b = m_num_boundaries;
	ptrdiff_t v = get_num_of_vertices();
	ptrdiff_t e = get_num_of_edges()/2;
	ptrdiff_t f = get_num_of_faces();

	m_genus = (int)(2*c+e-b-f-v)/2;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::remove_hanged_vertices()
{

	bool find = false;
	VERTEX_ITER viter = vertices_list->begin();
	for (;viter!=vertices_list->end(); viter++)
	{
		if ((*viter)->edge==NULL)
		{
			find = true;
			break;
		}
	}

	if (find)
	{
		VERTEX_LIST* new_vertices_list = new VERTEX_LIST;
		int i = 0;
		for (viter = vertices_list->begin();viter!=vertices_list->end(); viter++)
		{
			if ((*viter)->edge!=NULL)
			{
				new_vertices_list->push_back(*viter);
				(*viter)->id = i;
				i++;
			}
			else
				delete *viter;
		}
		delete vertices_list;
		vertices_list = new_vertices_list;

		EDGE_LIST* new_edges_list = new EDGE_LIST;
		i = 0;
		for (EDGE_ITER eiter = edges_list->begin();eiter!=edges_list->end(); eiter++)
		{
			if ((*eiter)->face!=NULL||(*eiter)->pair->face!=NULL)
			{
				new_edges_list->push_back(*eiter);
				(*eiter)->id = i;
				i++;
			}
			else
				delete *eiter;
		}
		delete edges_list;
		edges_list = new_edges_list;
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::reset_vertices_tag(bool tag_status)
{
	VERTEX_ITER viter = vertices_list->begin();
	for (; viter != vertices_list->end(); viter++)
		(*viter)->tag = tag_status;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::reset_faces_tag(bool tag_status)
{
	FACE_ITER fiter = faces_list->begin();
	for (; fiter !=faces_list->end(); fiter++)
		(*fiter)->tag = tag_status;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::reset_edges_tag(bool tag_status)
{
	EDGE_ITER eiter = edges_list->begin();
	for (; eiter != edges_list->end(); eiter++)
		(*eiter)->tag = tag_status;
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::reset_all_tag(bool tag_status)
{
	reset_edges_tag(tag_status);
	reset_faces_tag(tag_status);
	reset_vertices_tag(tag_status);
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::translate(const vec3& tran_V)
{
	VERTEX_ITER viter = vertices_list->begin();
	for (; viter != vertices_list->end(); viter++)
	{
		HE_vert* hv = *viter;
		hv->pos += tran_V;
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::scale(double factorx, double factory, double factorz)
{
	VERTEX_ITER viter = vertices_list->begin();
	for (; viter != vertices_list->end(); viter++)
	{
		HE_vert* hv = *viter;
		hv->pos[0] *= factorx;
		hv->pos[1] *= factory;
		//hv->pos[2] *= factorz;
	}
}
//////////////////////////////////////////////////////////////////////////

void Mesh3D::init_edge_tag()
{
	reset_edges_tag(false);
	EDGE_ITER eiter = edges_list->begin();
	for (; eiter != edges_list->end(); eiter++)
	{
		HE_edge* he = *eiter;
		if (he->tag==false && he->pair->tag == false)
		{
			he->tag = true;
		}
	}
}


void Mesh3D::swap_edge( HE_edge* triedge )
{
	if (triedge==NULL || triedge->face==NULL || triedge->pair->face==NULL || triedge->face->valence!=3 || triedge->pair->face->valence != 3)
	{
		return;
	}

	HE_vert* hv2 = triedge->vert;
	HE_vert* hv1 = triedge->pair->vert;
	HE_vert* hv3 = triedge->next->vert;
	HE_vert* hv4 = triedge->pair->next->vert;


	//check whether the swap causes non_manifold
	HE_edge* che = hv3->edge;
	do 
	{
		if (che->vert == hv4)
		{
			return;
		}
		che = che->pair->next;
	} while (che != hv3->edge);


	HE_face* hf1 = triedge->face;
	HE_face* hf2 = triedge->pair->face;


	HE_edge* he1 = triedge->next->next;
	HE_edge* he2 = triedge->next;
	HE_edge* he3 = triedge->pair->next->next;
	HE_edge* he4 = triedge->pair->next;




	triedge->vert = hv3;
	triedge->next = he1;
	triedge->prev = he4;

	triedge->pair->vert = hv4;
	triedge->pair->next = he3;
	triedge->pair->prev = he2;

	hf1->edge = triedge;
	hf2->edge = triedge->pair;

	he1->face = hf1;
	he2->face = hf2;
	he3->face = hf2;
	he4->face = hf1;

	he1->next = he4;
	he1->prev = triedge;
	he2->next = triedge->pair;
	he2->prev = he3;
	he3->next = he2;
	he3->prev = triedge->pair;
	he4->next = triedge;
	he4->prev = he1;

	if (hv1->edge == triedge)
	{
		hv1->edge = he1->pair;
	}
	if (hv2->edge == triedge->pair)
	{
		hv2->edge = he3->pair;
	}

	hv1->degree--;
	hv2->degree--;
	hv3->degree++;
	hv4->degree++;

}

void Mesh3D::build_kdtree(){
	if (m_kdtree_)
	{
		delete m_kdtree_;
	}

	std::vector<vec3> samples;
	for (int fitr = 0; fitr < get_num_of_faces(); fitr++)
	{
		HE_face* fh = get_face(fitr);
		vec3 fvpts [3];
		fvpts[0] = fh->edge->vert->pos, fvpts[1] = fh->edge->next->vert->pos, fvpts[2] = fh->edge->next->next->vert->pos;

		for (int sitr = 0; sitr < facet_sample_num_; sitr++)
		{
			double alpha = ((double) rand() / (RAND_MAX));
			double beta = (1.0-alpha)*((double) rand() / (RAND_MAX));
			double gamma = 1-alpha-beta;

			samples.push_back(alpha*fvpts[0]+beta*fvpts[1]+gamma*fvpts[2]);
		}
	}

	ANNpointArray dataPts = annAllocPts(samples.size(), 3);

	for (int sitr = 0; sitr < samples.size(); sitr++)
	{
		dataPts[sitr][0] = samples[sitr][0];
		dataPts[sitr][1] = samples[sitr][1];
		dataPts[sitr][2] = samples[sitr][2];
	}

	m_kdtree_ = new ANNkd_tree(		// build search structure
		dataPts, 					// the data points
		samples.size(), 						// number of points
		3);
}

inline vec3 project_to_segment(vec3 pt, vec3 s0, vec3 s1){
	vec3 dir = normalize(s1-s0);
	double t = dot(pt-s0, dir)/(s1-s0).length();

	if (t<0)
	{
		t=0;
	}
	else if (t>1)
	{
		t=1;
	}

	return (1-t)*s0 + t*s1;
}

vec3 Mesh3D::project_pt(vec3 pt, HE_face* & fh){
	if (!m_kdtree_)
	{
		std::cout << "Error: no kdtree built" << std::endl;
		return pt;
	}

	ANNpoint		queryPt;				// query point
	int nnIdx ;								// allocate near neighbor indices
	double dists ;							// allocate near neighbor dists

	queryPt = annAllocPt(3);
	queryPt[0] = pt[0];
	queryPt[1] = pt[1];
	queryPt[2] = pt[2];

	m_kdtree_->annkSearch(					// search
		queryPt, 						// query point
		1, 								// number of near neighbors
		&nnIdx, 						// nearest neighbors (returned)
		&dists, 						// distance (returned)
		0.0);							// error bound

	annDeallocPt(queryPt);	

	HE_face* pj_fh = get_face(nnIdx/facet_sample_num_);
	
	//find projection point.
	HE_edge* pj_eh = pj_fh->edge;
	vec3 face_vts[3];
	face_vts[0] = pj_eh->vert->pos;
	pj_eh=pj_eh->next;
	face_vts[1] = pj_eh->vert->pos;
	pj_eh=pj_eh->next;
	face_vts[2] = pj_eh->vert->pos;
	
	vec3 fnormal = pj_fh->normal;
	vec3 pj_pt = dot(face_vts[0]-pt, fnormal)*fnormal + pt;

	vec3 area_element = cross(face_vts[1]-face_vts[0], face_vts[2]-face_vts[0]);
	for (int eitr = 0; eitr < 3; eitr++)
	{
		vec3 ev0 = face_vts[(eitr+1)%3], ev1 = face_vts[(eitr+2)%3];
		vec3 sub_element = cross(ev0-pj_pt, ev1-pj_pt);
		if (dot(sub_element, area_element)<0)
		{
			pj_pt = project_to_segment(pj_pt, ev0, ev1);
		}
	}

	fh = pj_fh;
	return pj_pt;
}

void Mesh3D::split_face(HE_face* face, HE_vert* vtx)
{
	int nv = (int)vertices_list->size(); 
	int ne = (int)edges_list->size();
	int nf = (int)faces_list->size();
	int face_valence = face->valence;

	assert(face_valence == 3);

	//vec3 V = face->GetCentroid();
	HE_vert* new_vertex = vtx;
	//new_vertex->id = nv;
	new_vertex->degree = face_valence;

	std::vector<HE_face* > new_faces(face_valence);
	std::vector<HE_edge* > new_edges(2*face_valence);

	new_edges.resize(2*face_valence);
	new_faces.resize(face_valence);

	std::vector<HE_edge*> face_edges(face_valence);
	std::vector<HE_vert*> face_vertices(face_valence);
	HE_edge* he = face->edge;
	for (int i = 0; i < face_valence; i++)
	{
		face_vertices[i] = he->vert;
		face_vertices[i]->degree++;
		face_edges[i] = he;
		he = he->next;
		if (i == 0)
		{
			new_faces[i] = face;
		}
		else
		{
			new_faces[i] = new HE_face;
			new_faces[i]->id = nf;
			nf++;
		}
		new_faces[i]->valence = 3;
		new_faces[i]->edge = face_edges[i];
		face_edges[i]->face = new_faces[i];
	}

	for (int i = 0; i < face_valence; i++)
	{
		new_edges[2*i] = new HE_edge;
		new_edges[2*i+1] = new HE_edge;
		new_edges[2*i]->pair = new_edges[2*i+1];
		new_edges[2*i+1]->pair = new_edges[2*i];
		new_edges[2*i]->face = new_faces[i];
		new_edges[2*i+1]->face = new_faces[(i+1)%face_valence];
		new_edges[2*i]->vert = new_vertex;
		new_edges[2*i+1]->vert = face_vertices[i];
		new_edges[2*i]->id = ne;
		new_edges[2*i+1]->id = ne+1;
		ne += 2;
	}

	new_vertex->edge = new_edges[1];

	for (int i = 0; i < face_valence; i++)
	{
		new_edges[2*i]->prev = face_edges[i];
		new_edges[2*i]->next = i==0? new_edges[2*face_valence-1] : new_edges[2*i-1];
		new_edges[2*i+1]->prev = i==face_valence-1? new_edges[0] : new_edges[2*i+2];
		new_edges[2*i+1]->next = face_edges[(i+1)%face_valence];

		new_edges[2*i]->prev->next = new_edges[2*i];
		new_edges[2*i]->next->prev = new_edges[2*i];
		new_edges[2*i+1]->prev->next = new_edges[2*i+1];
		new_edges[2*i+1]->next->prev = new_edges[2*i+1];
	}
	//vertices_list->push_back(new_vertex);
	edges_list->insert(edges_list->end(), new_edges.begin(), new_edges.end());
	std::vector<HE_face*>::iterator fiter = new_faces.begin();
	fiter++;
	faces_list->insert(faces_list->end(), fiter, new_faces.end());

}


void Mesh3D::make_delaunay( bool& topo_updated)
{


	topo_updated = false;

	std::vector<HE_edge*> edge_v;

	//update the regular triangulation by edge flipping
	reset_edges_tag(false);
	int edge_counter = 0;
	for (EDGE_ITER eitr = edges_list->begin(); eitr != edges_list->end(); ++eitr)
	{
		HE_edge* eh = *eitr;
		if (eh->tag || eh->pair->tag)
		{
			continue;
		}
		edge_v.push_back(eh);
		eh->tag = true;
		eh->pair->tag = true;
	}

	//std::set<HE_edge*> set_edges;
	//std::set<HE_face*> set_faces;
	//set_edges.insert(edges_list->begin(),edges_list->end());
	//set_faces.insert(faces_list->begin(),faces_list->end());

	const int max_itr = edges_list->size()*5;
	int itr_count = 0;

	while (edge_v.size()>0 && itr_count < max_itr)
	{
		itr_count ++;

		HE_edge* eh = edge_v.back();
		edge_v.pop_back();
		eh->tag = false;
		eh->pair->tag = false;

		if (eh->face == NULL || eh->pair->face == NULL)
		{
			continue;
		}

		//if (eh->spln_curv)
		//{
		//	continue;
		//}
		if (eh->cycle_curv)
		{
			continue;
		}

		if (!delaunay_test(eh))
		{

			if (0)
			{//unflippable


			} else{
				topo_updated = true;

				//flip
				swap_edge(eh);
				topo_updated = true;

				if (!eh->next->tag && !eh->next->pair->tag)
				{
					eh->next->tag = eh->next->pair->tag = true;
					edge_v.push_back(eh->next);
				}
				if (!eh->prev->tag && !eh->prev->pair->tag)
				{
					eh->prev->tag = eh->prev->pair->tag = true;
					edge_v.push_back(eh->prev);
				}
				if (!eh->pair->next->tag && !eh->pair->next->pair->tag)
				{
					eh->pair->next->tag = eh->pair->next->pair->tag = true;
					edge_v.push_back(eh->pair->next);
				}
				if (!eh->pair->prev->tag && !eh->pair->prev->pair->tag)
				{
					eh->pair->prev->tag = eh->pair->prev->pair->tag = true;
					edge_v.push_back(eh->pair->prev);
				}
			}
		}
	}

	//update_mesh();
	//write_obj("2d.obj");

	std::cout <<"Topo updated: "<<(topo_updated?"Yes":"No")<<std::endl;
}


bool Mesh3D::delaunay_test( HE_edge* triedge )
{
	vec3 v0, v1, e0, e1;
	v0 = triedge->next->vert->pos;
	v1 = triedge->pair->next->vert->pos;
	e0 = triedge->vert->pos;
	e1 = triedge->pair->vert->pos;

	double angle1 = acos(dot(normalize(e0-v0),normalize(e1-v0)));
	double angle2 = acos(dot(normalize(e0-v1),normalize(e1-v1)));

	if (angle1 + angle2 < M_PI)
	{
		return true;
	}

	return false;
}

PolyGraph_Edge* Mesh3D::vertex_on_curve(HE_vert* vh){
	HE_edge* eh = vh->edge;
	do
	{
		//if(eh->spln_curv)
		//	return eh->spln_curv;

		//if (eh->pair->spln_curv)
		//	return eh->pair->spln_curv;

		if(eh->cycle_curv)
			return eh->cycle_curv;

		if (eh->pair->cycle_curv)
			return eh->pair->cycle_curv;

		eh=eh->pair->next;
	} while (eh!=vh->edge);

	return NULL;
}

int Mesh3D::vertex_curves(HE_vert* vh){
	int num = 0;
	HE_edge* eh = vh->edge;
	do
	{
		//if(eh->spln_curv || eh->pair->spln_curv)
		//	num++;

		if(eh->cycle_curv || eh->pair->cycle_curv)
			num++;

		eh=eh->pair->next;
	} while (eh!=vh->edge);

	return num;
}

double Mesh3D::edge_chain_crvt(HE_edge* eij){
	HE_face* f0 = eij->face, *f1 = eij->pair->face;
	if (f0 && f1)
	{
		vec3 f0_axis = normalize(f0->edge->vert->pos - f0->edge->prev->vert->pos);
		vec3 f1_axis = normalize(f1->edge->vert->pos - f1->edge->prev->vert->pos);
		vec3 vec_ij = normalize(eij->vert->pos - eij->pair->vert->pos);
		vec3 vec_ji = -vec_ij;

		vec3 rdir0 = cross(f0_axis, vec_ij);
		vec3 rdir1 = cross(vec_ji, f1_axis);

		double angle0 , angle1;

		if (f0->edge == eij)
		{
			angle0 = 0;
		}
		else
		{
			angle0 = acos(dot(f0_axis, vec_ij));
			if (dot(rdir0, f0->normal) < 0)
			{
				angle0 *= -1;
			}
		}

		if (f1->edge == eij->pair)
		{
			angle1 = 0;
		}
		else
		{
			angle1 = acos(dot(f1_axis, vec_ji));
			if (dot(rdir1, f1->normal) < 0)
			{
				angle1 *= -1;
			}
		}

		double kij = angle0 + angle1 - M_PI;
		if (kij > M_PI)
		{
			kij -= 2*M_PI;
		}
		else if (kij <= - M_PI)
		{
			kij += 2*M_PI;
		}

		if (!Numeric::is_nan(kij))
		{
			return -kij;
		}
		else
		{
			return 0;
		}

		
	}

	return 0;
}

double Mesh3D::bc_dual_edge_length(HE_edge* edge){
	//for now a barycentric dual is used.
	if (edge->face && edge->pair->face)
	{
		double l0 = (edge->face->GetCentroid() - edge->GetMidPoint()).length();
		double l1 = (edge->pair->face->GetCentroid() - edge->GetMidPoint()).length();

		return l0+l1;
	}
	return (edge->face->GetCentroid() - edge->GetMidPoint()).length();
}

double Mesh3D::ortho_dual_edge_data(HE_edge* edge, vec3& dir, vec3& ndiff, bool use_cot_w /*= false*/) const{

	double edge_len = (edge->vert->pos - edge->prev->vert->pos).length();
	
	//if (edge_len < 1e-10)
	//{
	//	dir = vec3(0,0,0);
	//	ndiff = vec3(0,0,0);
	//	return;
	//}

	dir = -cross(edge->face->normal, edge->vert->pos - edge->prev->vert->pos); //normalize
	if (edge->face && edge->pair->face)
	{
		ndiff = (edge->pair->face->normal - edge->face->normal);//*edge_len

		if (use_cot_w)
		{
			//compute ortho dual length through cot. formula.
			vec3 pt[4];
			pt[0] = edge->prev->vert->pos;
			pt[1] = edge->vert->pos;
			pt[2] = edge->next->vert->pos;
			pt[3] = edge->pair->next->vert->pos;

			vec3 e20 = normalize(pt[0] - pt[2]);
			vec3 e21 = normalize(pt[1] - pt[2]);
			vec3 e30 = normalize(pt[0] - pt[3]);
			vec3 e31 = normalize(pt[1] - pt[3]);

			double cos_2 = dot(e20, e21);
			double cos_3 = dot(e30, e31);

			double sin_2 = cross(e20, e21).length();
			double sin_3 = cross(e30, e31).length();

			double cot_23 = cos_2/sin_2 + cos_3/sin_3;

			double dual_ratio = cot_23*0.5;
			//double dual_len = (pt[1]-pt[0]).length();

			if (Numeric::is_nan(dual_ratio)) // || dual_len > 1e6
			{
				dual_ratio = 0;
			}

			dir *= dual_ratio ;

			return dual_ratio;
		}
	}

	return 1;
}

vec3 Mesh3D::get_cvt_centroid(HE_vert* vh){

	vec3 cent (0,0,0);

	if (vh->fixed)
	{
		return vh->pos;
	}

	vec3 p0 = vh->pos;
	double V = 0.0 ;

	vec3 gradi;
	double Vi = 0;
	vec3 cvt_grad (0,0,0);

	/*For each facet*/
	HE_edge* eh = vh->edge;
	do
	{
		if (!eh->face)
		{
			eh = eh->pair->next;
			continue;
		}


		/*compute the mid points, the centroid of the facet, that defines the voronoi region on this facet*/
		double facet_area;
		vec3 pts[3], mid01, mid02, fcentroid, n1, n2;
		pts[0] = (p0);
		pts[1] = (eh->vert->pos);
		pts[2] = (eh->next->vert->pos);
		mid01 = 0.5*(pts[0]+pts[1]);
		mid02 = 0.5*(pts[0]+pts[2]);
		n1 = normalize( (dot( pts[2] - pts[1] , pts[0] - pts[1] ) / (pts[0] - pts[1]).length2()) * (pts[0] - pts[1]) - (pts[2] - pts[1]) );
		n2 = normalize( (dot( pts[1] - pts[2] , pts[0] - pts[2] ) / (pts[0] - pts[2]).length2()) * (pts[0] - pts[2]) - (pts[1] - pts[2]) );

		facet_area = tri_area(pts[0],pts[1],pts[2]);

		double alpha = (pts[1]-pts[2]).length2()*dot(pts[1]-pts[0],pts[2]-pts[0]),
			beta = (pts[0]-pts[2]).length2()*dot(pts[0]-pts[1],pts[2]-pts[1]),
			gamma = (pts[0]-pts[1]).length2()*dot(pts[0]-pts[2],pts[1]-pts[2]);
		fcentroid = 1.0/(8*facet_area*facet_area)*(alpha*pts[0]+beta*pts[1]+gamma*pts[2]);

		/*compute the energy and gradient*/


		//compute energy and traditional cvt gradient.
		//discuss the shape of facet (acute or not).
		if (alpha < 0)
		{
			//obtuse at p0
			vec3 p01, p02;
			double tan_210 = 2*(pts[0]-pts[2]).length2()*facet_area/beta, tan_120 = 2*(pts[0]-pts[1]).length2()*facet_area/gamma;
			p01 = mid01 + tan_210 * (mid01-pts[0]).length() * normalize(fcentroid - mid01);
			p02 = mid02 + tan_120 * (mid02-pts[0]).length() * normalize(fcentroid - mid02);

			Lloyd_energy(pts[0], pts[0], mid01, p01, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;

			Lloyd_energy(pts[0], pts[0], mid02, p02, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;

			Lloyd_energy(pts[0], pts[0], p01, p02, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;
		}
		else if (beta < 0)
		{
			//obtuse at p1
			vec3 p01;
			double tan_201 = 2*(pts[1]-pts[2]).length2()*facet_area/alpha;
			p01 = mid01 + tan_201 * (mid01-pts[0]).length() * normalize(fcentroid - mid01);

			Lloyd_energy(pts[0], pts[0], mid01, p01, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;
		}
		else if (gamma < 0)
		{
			//obtuse at p2
			vec3 p02;
			double tan_201 = 2*(pts[1]-pts[2]).length2()*facet_area/alpha;
			p02 = mid02 + tan_201 * (mid02-pts[0]).length() * normalize(fcentroid - mid02);

			Lloyd_energy(pts[0], pts[0], mid02, p02, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;
		}
		else{
			Lloyd_energy(pts[0], pts[0], mid01, fcentroid, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;

			Lloyd_energy(pts[0], pts[0], mid02, fcentroid, gradi, Vi);
			V += Vi;
			cvt_grad += gradi;
		}

		////cvt gradient second term
		//double elength_1 = (pts[1] - pts[0]).length(), elength_2 = (pts[2] - pts[0]).length();
		//vec3 cvt_g_2nd_term = 1.0/24.0 * ( elength_1*elength_1*elength_1*n1 + elength_2*elength_2*elength_2*n2 );
		//cvt_grad += cvt_g_2nd_term;

		eh = eh->pair->next;
	} while(eh != vh->edge);

	cent = vh->pos - cvt_grad / (2.0*V);

	return cent;
}

void Mesh3D::get_cvt_fg( double& f, std::vector<double>& g )
{
	f = 0;
	g.resize(g.size(),0);

	for (int i = 0; i < vertices_list->size(); i ++)
	{
		HE_vert* vh = get_vertex(i);
		if (vh->fixed)
		{
			continue;
		}

		vec3 p0 = vh->pos;
		double V = 0.0 ;

		vec3 gradi;
		double Vi = 0;
		vec3 cvt_grad (0,0,0);

		/*For each facet*/
		HE_edge* eh = vh->edge;
		do
		{
			if (!eh->face)
			{
				eh = eh->pair->next;
				continue;
			}


			/*compute the mid points, the centroid of the facet, that defines the voronoi region on this facet*/
			double facet_area;
			vec3 pts[3], mid01, mid02, fcentroid, n1, n2;
			pts[0] = (p0);
			pts[1] = (eh->vert->pos);
			pts[2] = (eh->next->vert->pos);
			mid01 = 0.5*(pts[0]+pts[1]);
			mid02 = 0.5*(pts[0]+pts[2]);
			n1 = normalize( (dot( pts[2] - pts[1] , pts[0] - pts[1] ) / (pts[0] - pts[1]).length2()) * (pts[0] - pts[1]) - (pts[2] - pts[1]) );
			n2 = normalize( (dot( pts[1] - pts[2] , pts[0] - pts[2] ) / (pts[0] - pts[2]).length2()) * (pts[0] - pts[2]) - (pts[1] - pts[2]) );

			facet_area = tri_area(pts[0],pts[1],pts[2]);

			double alpha = (pts[1]-pts[2]).length2()*dot(pts[1]-pts[0],pts[2]-pts[0]),
				beta = (pts[0]-pts[2]).length2()*dot(pts[0]-pts[1],pts[2]-pts[1]),
				gamma = (pts[0]-pts[1]).length2()*dot(pts[0]-pts[2],pts[1]-pts[2]);
			fcentroid = 1.0/(8*facet_area*facet_area)*(alpha*pts[0]+beta*pts[1]+gamma*pts[2]);

			/*compute the energy and gradient*/


			//compute energy and traditional cvt gradient.
			//discuss the shape of facet (acute or not).
			if (alpha < 0)
			{
				//obtuse at p0
				vec3 p01, p02;
				double tan_210 = 2*(pts[0]-pts[2]).length2()*facet_area/beta, tan_120 = 2*(pts[0]-pts[1]).length2()*facet_area/gamma;
				p01 = mid01 + tan_210 * (mid01-pts[0]).length() * normalize(fcentroid - mid01);
				p02 = mid02 + tan_120 * (mid02-pts[0]).length() * normalize(fcentroid - mid02);

				f += Lloyd_energy(pts[0], pts[0], mid01, p01, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;

				f += Lloyd_energy(pts[0], pts[0], mid02, p02, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;

				f += Lloyd_energy(pts[0], pts[0], p01, p02, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;
			}
			else if (beta < 0)
			{
				//obtuse at p1
				vec3 p01;
				double tan_201 = 2*(pts[1]-pts[2]).length2()*facet_area/alpha;
				p01 = mid01 + tan_201 * (mid01-pts[0]).length() * normalize(fcentroid - mid01);

				f += Lloyd_energy(pts[0], pts[0], mid01, p01, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;
			}
			else if (gamma < 0)
			{
				//obtuse at p2
				vec3 p02;
				double tan_201 = 2*(pts[1]-pts[2]).length2()*facet_area/alpha;
				p02 = mid02 + tan_201 * (mid02-pts[0]).length() * normalize(fcentroid - mid02);

				f += Lloyd_energy(pts[0], pts[0], mid02, p02, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;
			}
			else{
				f += Lloyd_energy(pts[0], pts[0], mid01, fcentroid, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;

				f += Lloyd_energy(pts[0], pts[0], mid02, fcentroid, gradi, Vi);
				V += Vi;
				cvt_grad += gradi;
			}

			//cvt gradient second term
			double elength_1 = (pts[1] - pts[0]).length(), elength_2 = (pts[2] - pts[0]).length();
			vec3 cvt_g_2nd_term = 1.0/24.0 * ( elength_1*elength_1*elength_1*n1 + elength_2*elength_2*elength_2*n2 );
			cvt_grad += cvt_g_2nd_term;

			eh = eh->pair->next;
		} while(eh != vh->edge);


		g[vh->id*3] = cvt_grad.x;
		g[vh->id*3+1] = cvt_grad.y;
		g[vh->id*3+2] = cvt_grad.z;
	}
}

Mesh3D* Mesh3D::make_copy_regular(){
	Mesh3D* new_mesh = new Mesh3D;

	VERTEX_ITER viter = vertices_list->begin();
	for (; viter!=vertices_list->end(); viter++)
	{
		HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

		//hpan
		if ((*viter)->fixed)
		{
			vh->fixed = true;
		}
	}

	new_mesh->edges_list = new EDGE_LIST;

	for (EDGE_ITER eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++)
	{
		HE_edge* he = *eiter;
		HE_edge* nedge = new HE_edge;
		nedge->vert = new_mesh->get_vertex(he->vert->id);
		nedge->id = he->id;
		new_mesh->edges_list->push_back(nedge);
	}

	PTR_EDGE_LIST edgelist = new_mesh->get_edges_list();

	std::sort(edgelist->begin(), edgelist->end(), CompareEdgeID);

	
	VERTEX_ITER cviter = new_mesh->get_vertices_list()->begin();
	for (viter = vertices_list->begin(); viter!=vertices_list->end(); viter++, cviter++)
	{
		if ((*viter)->edge != NULL)
		{
			(*cviter)->edge = new_mesh->get_edge((*viter)->edge->id);
		}
	}

	for (EDGE_ITER ceitr = new_mesh->get_edges_list()->begin(), eitr = edges_list->begin(); ceitr != new_mesh->get_edges_list()->end(); ceitr++, eitr++)
	{
		HE_edge* neh = *ceitr, *eh = *eitr;
		neh->next = new_mesh->get_edge(eh->next->id);
		neh->prev = new_mesh->get_edge(eh->prev->id);
		neh->pair = new_mesh->get_edge(eh->pair->id);
	}

	FACE_ITER fiter = faces_list->begin();
	new_mesh->faces_list = new FACE_LIST;

	for (;fiter!=faces_list->end(); fiter++) 
	{
		HE_face* hf = *fiter;
		HE_edge* edge = hf->edge;

		HE_face* nface = new HE_face;
		nface->edge = new_mesh->get_edge(edge->id);
		nface->valence = hf->valence;

		nface->id = hf->id;
		new_mesh->faces_list->push_back(nface);
		
		do {
			
			new_mesh->get_edge(edge->id)->face = nface;

			edge = edge->next;
		} while(edge!=hf->edge);

	}



	new_mesh->update_mesh();



	//hpan. Copy pointer to curves.
	for (int i = 0; i < edgelist->size(); i++)
	{
		//(*edgelist)[i]->spln_curv = (*edges_list)[i]->spln_curv;
		(*edgelist)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
	}


	return new_mesh;
}

