#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

/*Fonciton du tp */

/* get the vertexhandle of current edge */
void MainWindow::getPointsEdge(MyMesh *_mesh,EdgeHandle *edge,VertexHandle *vh1,
                               VertexHandle *vh2)
{
    auto v1= _mesh->to_vertex_handle(_mesh->halfedge_handle(*edge, 0));
    auto v2= _mesh->from_vertex_handle(_mesh->halfedge_handle(*edge, 0));
    *vh1 = v1;
    *vh2 = v2;
    //qDebug()<<"sommet de l'arete :"<<vh1->idx()<<" "<<vh2->idx();
}

double MainWindow::getLenghtPoint1Point2(MyMesh *_mesh,VertexHandle v1,VertexHandle v2)
{
    //Euler
    double a,b,c,edgeLenght;
    a = _mesh->point(v2)[0] - _mesh->point(v1)[0];
    b = _mesh->point(v2)[1] - _mesh->point(v1)[1];
    c = _mesh->point(v2)[2] - _mesh->point(v1)[2];

    edgeLenght = sqrt(a*a+b*b+c*c);
    return edgeLenght;
}

/*return lenght of current edge*/
double MainWindow::getLenghtEdge(MyMesh *_mesh,EdgeHandle edge)
{
    VertexHandle v1,v2;
    double edgeLenght;
    getPointsEdge(_mesh,&edge,&v1,&v2);

    //Euler
    double a,b,c;
    a = _mesh->point(v2)[0] - _mesh->point(v1)[0];
    b = _mesh->point(v2)[1] - _mesh->point(v1)[1];
    c = _mesh->point(v2)[2] - _mesh->point(v1)[2];

    edgeLenght = sqrt(a*a+b*b+c*c);
    return edgeLenght;
}

/* get the commun faces shared by current edge */
void MainWindow::getFacesSharedByEdge(MyMesh *_mesh,FaceHandle *Fh0,FaceHandle *Fh1
                                        ,EdgeHandle *edge)
{
    VertexHandle vh1,vh2;
    getPointsEdge(&mesh,edge,&vh1,&vh2);
    std::vector<int> vertexFace_idx;

    //On parcours toutes les faces adjacentes au sommet vh1
    for(MyMesh::VertexFaceIter curFace = _mesh->vf_iter(vh1); curFace.is_valid();
        curFace ++){
       vertexFace_idx.push_back(curFace->idx());
    }

    //On parcours toutes les faces adjacentes au sommet vh2
    for(MyMesh::VertexFaceIter curFace = _mesh->vf_iter(vh2);
        curFace.is_valid(); curFace ++){
       vertexFace_idx.push_back(curFace->idx());
    }

    int element1,element2;
    std::vector<int> communfaces_idx;

    //On ne garde que les points commun aux deux faces
    for(unsigned int i=0; i<vertexFace_idx.size(); i++){
        element1 = vertexFace_idx[i];
        for(unsigned int j=0; j<vertexFace_idx.size(); j++){
            element2 = vertexFace_idx[j];
            if((i != j) &&(element1 == element2)){
                communfaces_idx.push_back(element1);
            }
        }
    }

    //suppression des doublons
    std::sort(communfaces_idx.begin(), communfaces_idx.end());
    communfaces_idx.erase(std::unique(communfaces_idx.begin(), communfaces_idx.end()), communfaces_idx.end());

    if(communfaces_idx.size() > 1){
        *Fh0 = _mesh->face_handle(communfaces_idx[0]);
        *Fh1 = _mesh->face_handle(communfaces_idx[1]);
    }
    qDebug()<<"face adjacentes : "<<Fh0->idx()<<Fh1->idx();
}

/* get the opposite vertex of neighboor faces */
void MainWindow::getOppositeVertexOfNeighoorFaces(MyMesh *_mesh,FaceHandle fh0,FaceHandle fh1
                                        ,VertexHandle *vh0,VertexHandle *vh1)
{
    std::vector<int> faceVertex_idx;

    //On parcours tous les sommets de la face fh0
    for(MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh0); curVert.is_valid();
        curVert ++){
       faceVertex_idx.push_back(curVert->idx());
    }

    //On parcours tous les sommets de la face fh1
    for(MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh1); curVert.is_valid();
        curVert ++){
       faceVertex_idx.push_back(curVert->idx());
    }

    int element1,element2,count=0;
    std::vector<int> oppositeVertex_idx;

    //On ne garde que les points non commun au deux faces
    for(unsigned int i=0; i<faceVertex_idx.size();i++){
        element1 = faceVertex_idx[i];
        count = 0;
        for(unsigned int j=0; j<faceVertex_idx.size(); j++){
            element2 = faceVertex_idx[j];
            if((i!=j) && (element1 == element2)){
                count ++;
            }
        }
        if(count == 0) //unique element
            oppositeVertex_idx.push_back(element1);
    }

    if(faceVertex_idx.size() > 1){
        *vh0 = _mesh->vertex_handle(oppositeVertex_idx[0]);
        *vh1 = _mesh->vertex_handle(oppositeVertex_idx[1]);
    }

}

void MainWindow::edgesSplit(MyMesh *_mesh, double maxLength)
{
    //called befor a vertex/face/edge can be deleted, it grants to the status attribute
    _mesh->request_face_status();
    _mesh->request_edge_status();
    _mesh->request_vertex_status();

    for(MyMesh::EdgeIter eh = _mesh->edges_begin(); eh!= _mesh->edges_end(); ++eh)
    {
           EdgeHandle edge = *eh;
           edgeSplit(_mesh,&edge,maxLength);
           qDebug()<<"id arete courante"<<eh->idx();
           qDebug()<<"nb sommet"<<_mesh->n_vertices(); // nombre de sommets
           qDebug()<<"nb arete"<<_mesh->n_edges(); // nombre d'arêtes
           qDebug()<<"nb demis aret"<<_mesh->n_halfedges(); // nombre de demi-arêtes, (=_mesh->n_edges()*2)
           qDebug()<<"nb face"<<_mesh->n_faces(); // nombre de faces

        break;
    }
   //delate all elements that are marked as deleted from memory
   _mesh->garbage_collection();
}

void MainWindow::edgeSplit(MyMesh *_mesh,EdgeHandle *edge, double maxLength)
{

    VertexHandle vh1,vh2,vh3,vh4;
    FaceHandle fh1,fh2;

    getPointsEdge(&mesh,edge,&vh1,&vh2);

    int oldNumberFace = _mesh->n_faces();

    double a,b,c;
    a = _mesh->point(vh2)[0]-_mesh->point(vh1)[0];
    b = _mesh->point(vh2)[1]-_mesh->point(vh1)[1];
    c = _mesh->point(vh2)[2]-_mesh->point(vh1)[2];

    double edgeLength = sqrt(a*a+b*b+c*c);
    qDebug()<<"edgeLength : "<<edgeLength;

    if(edgeLength > maxLength)
    {
        //on sélectionne le(s) face(s) partagé(es) par l'arête courante
        getFacesSharedByEdge(_mesh, &fh1, &fh2, edge);

        //on récupère les sommets opposés à l'arête
        getOppositeVertexOfNeighoorFaces(_mesh,fh1,fh2,&vh3,&vh4);
        qDebug()<<"vertex opposé à l'aret:"<<vh3.idx()<<vh4.idx();
        qDebug()<<"vertex de l'arete :"<<vh1.idx()<<vh2.idx();

        //on détruit les faces adjacentes à l'arête ainsi que l'arête
        _mesh->delete_face(fh1,false);
        _mesh->delete_face(fh2,false);
        _mesh->delete_edge(*edge,false);

        //on ajoute le sommet au centre de l'arete détruite
        double x_vh1=_mesh->point(vh1)[0],y_vh1=_mesh->point(vh1)[1],z_vh1=_mesh->point(vh1)[2];
        MyMesh::Point middlePoint(x_vh1+a/2, y_vh1+b/2, z_vh1 + c/2);
        VertexHandle middleVh = mesh.add_vertex(middlePoint);

        //on reconstruit les faces
        std::vector<MyMesh::VertexHandle> uneNouvelleFace;
        //face 0
        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(middleVh);
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh2.idx()));
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh4.idx()));
        _mesh->add_face(uneNouvelleFace);
        //face 1
        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(middleVh);
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh4.idx()));
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh1.idx()));
        _mesh->add_face(uneNouvelleFace);
        //face 2
        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(middleVh);
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh1.idx()));
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh3.idx()));
        _mesh->add_face(uneNouvelleFace);
        //face 3
        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(middleVh);
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh3.idx()));
        uneNouvelleFace.push_back(_mesh->vertex_handle(vh2.idx()));
        _mesh->add_face(uneNouvelleFace);

        _mesh->request_face_normals();
        //delate all elements that are marked as deleted from memory
        _mesh->garbage_collection();

        qDebug()<<"old number face :"<<oldNumberFace<<"new :"<<_mesh->n_faces();
        if(oldNumberFace >= _mesh->n_faces())
        {
            qDebug()<<"is inferieur";
            //face 0
            uneNouvelleFace.clear();
            uneNouvelleFace.push_back(middleVh);
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh1.idx()));
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh4.idx()));
            _mesh->add_face(uneNouvelleFace);
            //face 1
            uneNouvelleFace.clear();
            uneNouvelleFace.push_back(middleVh);
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh4.idx()));
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh2.idx()));
            _mesh->add_face(uneNouvelleFace);
            //face 2
            uneNouvelleFace.clear();
            uneNouvelleFace.push_back(middleVh);
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh2.idx()));
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh3.idx()));
            _mesh->add_face(uneNouvelleFace);
            //face 3
            uneNouvelleFace.clear();
            uneNouvelleFace.push_back(middleVh);
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh3.idx()));
            uneNouvelleFace.push_back(_mesh->vertex_handle(vh1.idx()));
            _mesh->add_face(uneNouvelleFace);

            _mesh->request_face_normals();
            //delate all elements that are marked as deleted from memory
            _mesh->garbage_collection();
        }

    }
}

/* Get barycenter of neighbors 1 ring vertices  */
void MainWindow::oneRingNeighbors_Barycenter(MyMesh * _mesh, std::vector<double> *barycenter,
                                             VertexHandle v)
{
    double x=0.0,y=0.0,z=0.0;

    std::vector<int> ring1Neigbors;

    //On parcours toutes les faces adjacentes à v
    for(MyMesh::VertexFaceIter curFace = _mesh->vf_iter(v); curFace.is_valid(); curFace ++){
        for(MyMesh::FaceVertexIter curVert = _mesh->fv_iter(curFace); curVert.is_valid();
            curVert ++){
            if(curVert->idx() != v.idx())
            {
               ring1Neigbors.push_back(curVert->idx());
            }
        }
    }

    //On elimine les doublons
    std::sort(ring1Neigbors.begin(), ring1Neigbors.end());
    ring1Neigbors.erase(std::unique(ring1Neigbors.begin(), ring1Neigbors.end()),
                          ring1Neigbors.end());

    for(int i=0; i<ring1Neigbors.size(); i++){
        VertexHandle curVert = _mesh->vertex_handle(ring1Neigbors[i]);
        x += _mesh->point(curVert)[0];
        y += _mesh->point(curVert)[1];
        z += _mesh->point(curVert)[2];
    }
    barycenter->push_back(x/(double)ring1Neigbors.size());
    barycenter->push_back(y/(double)ring1Neigbors.size());
    barycenter->push_back(z/(double)ring1Neigbors.size());
}

void MainWindow::getVertexSurfaceApproximation(MyMesh *_mesh,VertexHandle v,
                                               std::vector<double> barycenterPoint)
{
    /*To do*/

    //double x,y,z;
    //x = barycenterPoint[0] - mesh.normal(v)[0];
    /*newVertex = _mesh->vertex_handle(v.idx());
    MyMesh::Point NewCoordonatesCurrentVertex(x, y, z);
    _mesh->set_point(newVertex, NewCoordonatesCurrentVertex);*/
    //+mesh.normal(*v_it)
}

/* set coordinate of current vertex */
void MainWindow::setVertexCoordonate(MyMesh *_mesh,VertexHandle v,double x,double y, double z)
{
    VertexHandle newVertex = _mesh->vertex_handle(v.idx());
    MyMesh::Point NewCoordonatesCurrentVertex(x, y, z);
    _mesh->set_point(newVertex, NewCoordonatesCurrentVertex);
}

/* Shift all vertices of Mesh */
void MainWindow::verticesShift(MyMesh *_mesh)
{
    for(MyMesh::VertexIter currentVertex = _mesh->vertices_begin(); currentVertex!= _mesh->vertices_end(); ++currentVertex){
        std::vector<double> barycenterPoint;
        oneRingNeighbors_Barycenter(_mesh, &barycenterPoint, currentVertex);
        setVertexCoordonate(_mesh,currentVertex,barycenterPoint[0],barycenterPoint[1],barycenterPoint[2]);
    }
}

/* get angle create by the two vectors of the differents edges shared by vertexID on the faceID */
double MainWindow::getAngle(MyMesh* _mesh, int vertexID,  int faceID)
{
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    FaceHandle fh = _mesh->face_handle(faceID);
    QVector<VertexHandle> listePoints;
    QVector<VertexHandle> listePointsOnFace;
    QVector<float> vectors;

    //On identifie les point voisins de vertexID qui appartiennent à faceID
    //tout les points voisins de vertexID
    for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh); curVertex.is_valid(); curVertex ++)
    {
        VertexHandle v = *curVertex;
        listePoints.append(v);
    }
    //parmis ces points ceux qui appartiennent à la faceID

    for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid()
        ; curVertex ++){
        VertexHandle v = *curVertex;
        if(listePoints.contains(v) && v.idx() != vertexID)
        {
            listePointsOnFace.append(v);
        }
    }

    //On créer des vecteurs a partir des point obtenu
    for(int i=0; i<listePointsOnFace.size();i++){
        vectors.append((_mesh->point(listePointsOnFace[i])[0])-(_mesh->point(vh)[0]));
        vectors.append((_mesh->point(listePointsOnFace[i])[1])-(_mesh->point(vh)[1]));
        vectors.append((_mesh->point(listePointsOnFace[i])[2])-(_mesh->point(vh)[2]));
    }

    //on normalise les vecteurs obtenu
    QVector<float> norme;
    float tmp = vectors[0]*vectors[0]+vectors[1]*vectors[1]+vectors[2]*vectors[2];
    norme.append(sqrt(tmp)); //sqrt(x^2+y^2+z^2) : C'est la norme du premier vecteur
    tmp = vectors[3]*vectors[3]+vectors[4]*vectors[4]+vectors[5]*vectors[5];
    norme.append(sqrt(tmp)); //C'est la norme du second vecteurs

    //normalisation des vecteur :  chacun des vecteur / par la norme
    for(int i=0; i<listePointsOnFace.size(); i++){
        vectors[i*3] = vectors[i*3]/norme[i];
        vectors[i*3+1] = vectors[i*3+1]/norme[i];
        vectors[i*3+2] = vectors[i*3+2]/norme[i];
    }

    float a = vectors[0]*vectors[0+3];
    float b = vectors[1]*vectors[1+3];
    float c = vectors[2]*vectors[2+3];
    float prodScal = a+b+c;
    abs(prodScal); //garantie que le produit scalaire soit positif

    float angle = acos(prodScal)/**180/M_PI*/;
    return angle; //en radians
}

void MainWindow::collapseEdge(MyMesh* _mesh, int edgeID)
{
    /*
     * cette fonction utilise l'opérateur collapse pour supprimer l'arête d'index edgeID
     * Attention à ne pas oublier garbage_collection() !
     */

     // Collapse edge
    _mesh->edge_handle(edgeID);
    _mesh->edge_handle(edgeID);
    MyMesh::HalfedgeHandle heh= _mesh->halfedge_handle(_mesh->edge_handle(edgeID),0);
    MyMesh::VertexHandle v1 = _mesh->to_vertex_handle(heh);
    MyMesh::VertexHandle v0 = _mesh->from_vertex_handle(heh);
    MyMesh::Point point0 = _mesh->point(v0);
    MyMesh::Point point1 = _mesh->point(v1);

    MyMesh::Point point = (point1 + point0)/2;
    _mesh->point(v1) = point;
    mesh.collapse(heh);

    _mesh->garbage_collection();
}

int MainWindow::getVertexValence(MyMesh *_mesh, VertexHandle v)
{
    int valence=0;
    for(MyMesh::VertexVertexIter curVert = _mesh->vv_iter(v); curVert.is_valid(); curVert ++){
        valence ++;
    }
    return valence;
}

void MainWindow::valenceEgalisation(MyMesh *_mesh)
{
    // Request required status flags

    int test = 0;
    for(MyMesh::EdgeIter eh = _mesh->edges_begin(); eh!=_mesh->edges_end(); ++eh)
    {
        if(test < 2)
        {
         test ++;
         EdgeHandle curent_edge = *eh;
         VertexHandle a,b,c,d;
         FaceHandle fh0,fh1;

         getPointsEdge(_mesh,&curent_edge,&c,&d);
         getFacesSharedByEdge(_mesh,&fh0,&fh1,&curent_edge);
         getOppositeVertexOfNeighoorFaces(_mesh,fh0,fh1,&c,&d);

         int valenceObjective = 6; //by default
         int valenceA = getVertexValence(_mesh,a);
         int valenceB = getVertexValence(_mesh,b);
         int valenceC = getVertexValence(_mesh,c);
         int valenceD = getVertexValence(_mesh,d);

         if(_mesh->is_boundary(*eh))
         {
            valenceObjective = 4;
         }

         int deviationAvant = abs(valenceA-valenceObjective)+abs(valenceA-valenceObjective)
                 +abs(valenceA-valenceObjective)+abs(valenceA-valenceObjective);

         // Find this edge and then flip it
         for(MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
                 if(!mesh.is_boundary(*it)) {
                         // Flip edge
                         mesh.flip(*it);
                 }
         }

         getPointsEdge(_mesh,&curent_edge,&c,&d);
         getFacesSharedByEdge(_mesh,&fh0,&fh1,&curent_edge);
         getOppositeVertexOfNeighoorFaces(_mesh,fh0,fh1,&c,&d);

         valenceObjective = 6; //by default
         valenceA = getVertexValence(_mesh,a);
         valenceB = getVertexValence(_mesh,b);
         valenceC = getVertexValence(_mesh,c);
         valenceD = getVertexValence(_mesh,d);

         if(_mesh->is_boundary(*eh))
         {
            valenceObjective = 4;
         }

         int deviationApres = abs(valenceA-valenceObjective)+abs(valenceA-valenceObjective)
                 +abs(valenceA-valenceObjective)+abs(valenceA-valenceObjective);

         if(deviationAvant <= deviationApres)
         {
             for(MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
                     if(!mesh.is_boundary(*it)) {
                             // Flip edge
                             mesh.flip(*it);
                     }
             }
         }
        }
    }
}

void MainWindow::collapsShortEdge(MyMesh *_mesh, double minLength, double maxLength)
{
    int count =0;
    bool hasMinimalEdge=true,fusionner;
    // Request required status flags
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();

    //do{
    //while(hasMinimalEdge == true)
    //{
        for(MyMesh::EdgeIter eh = _mesh->edges_begin(); eh!=_mesh->edges_end() /*&&*/ /*count<*/; ++eh)
        {
            qDebug()<<"=====================";
            hasMinimalEdge = false;
            VertexHandle vh1,vh2;
            EdgeHandle curent_edge = *eh;
            if(_mesh->is_collapse_ok(_mesh->halfedge_handle(curent_edge.idx())))
            {
                qDebug()<<"longueur arete"<<getLenghtEdge(_mesh, eh);
                qDebug()<<"longueur min"<<minLength;

                if(getLenghtEdge(_mesh, curent_edge) < minLength && !(_mesh->is_boundary(curent_edge)))
                {
                    hasMinimalEdge = true;
                    getPointsEdge(&mesh,&curent_edge,&vh1,&vh2);
                    fusionner = true;
                    for(MyMesh::VertexVertexIter curVert = _mesh->vv_iter(vh1); curVert.is_valid(); curVert ++)
                    {
                        qDebug()<<"longueur p1p2"<<getLenghtPoint1Point2(_mesh,vh2,*curVert);
                        if(getLenghtPoint1Point2(_mesh,vh2,*curVert) > maxLength)
                        {
                            fusionner = false;
                        }
                    }
                    if(fusionner)
                    {
                        for(MyMesh::HalfedgeIter it = _mesh->halfedges_begin(); it != _mesh->halfedges_end(); ++it)
                        {
                          if( _mesh->to_vertex_handle(*it) ==  vh1 &&
                              _mesh->from_vertex_handle(*it) == vh2
                              )

                          {
                              //collapseEdge(_mesh,curent_edge.idx());
                            // Collapse edge
                              count ++;
                             qDebug()<<"arete : ";
                            _mesh->collapse(*it);
                            _mesh->garbage_collection();

                            //_mesh->garbage_collection();
                            // permet de nettoyer le maillage et de garder la cohérence des indices après un collapse
                            //test = true;
                               //qDebug()<<"la";
                            break;

                          }
                        }

                    }
                }
            }
        }
    //}
    //}
    //}
    //while(hasMinimalEdge == true && test == false);
}

double MainWindow::getmeshAngleQuality(MyMesh *_mesh)
{
    double critere  = 0.0;

    for(MyMesh::FaceIter fh = _mesh->faces_begin(); fh!= _mesh->faces_end(); ++fh)
    {
        std::vector<int> vertexIndex;
        std::vector<double> angles;
        for(MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh); curVert.is_valid(); curVert ++)
        {
            vertexIndex.push_back(curVert->idx());
        }
        angles.push_back( getAngle(_mesh,vertexIndex[0],fh->idx()));
        angles.push_back( getAngle(_mesh,vertexIndex[1],fh->idx()));
        angles.push_back(getAngle(_mesh,vertexIndex[2],fh->idx()));

        const auto [min, max] = std::minmax_element(begin(angles), end(angles));
        critere += *max - *min;
        //qDebug()<<"ecart max angle"<<*max - *min;
    }
    return critere;
}

/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

/* **** fin de la partie boutons et IHM **** */

// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_remaillage_incremental_clicked()
{

    EdgeHandle e = mesh.edge_handle(3);
    VertexHandle v = mesh.vertex_handle(1);
    double lao = 2.83;//largeur arête objectif
    FaceHandle f0,f1;
    VertexHandle v0,v1;

    edgesSplit(&mesh,1.5);
    //verticesShift(&mesh);
    //getmeshAngleQuality(&mesh);
    //collapsShortEdge(&mesh, (4.0/5.0)*lao, (4.0/3.0)*lao);
    //valenceEgalisation(&mesh);
    mesh.update_normals();
    resetAllColorsAndThickness(&mesh);
    displayMesh(&mesh);
}
