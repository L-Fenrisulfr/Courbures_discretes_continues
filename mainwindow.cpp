#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "courbures.h"


/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj
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

// exemple pour construire un mesh face par face
void MainWindow::on_pushButton_generer_clicked()
{
    MyMesh mesh;

    // on construit une liste de sommets
    MyMesh::VertexHandle sommets[8];
    sommets[0] = mesh.add_vertex(MyMesh::Point(-1, -1,  1));
    sommets[1] = mesh.add_vertex(MyMesh::Point( 1, -1,  1));
    sommets[2] = mesh.add_vertex(MyMesh::Point( 1,  1,  1));
    sommets[3] = mesh.add_vertex(MyMesh::Point(-1,  1,  1));
    sommets[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
    sommets[5] = mesh.add_vertex(MyMesh::Point( 1, -1, -1));
    sommets[6] = mesh.add_vertex(MyMesh::Point( 1,  1, -1));
    sommets[7] = mesh.add_vertex(MyMesh::Point(-1,  1, -1));


    // on construit des faces à partir des sommets

    std::vector<MyMesh::VertexHandle> uneNouvelleFace;

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[3]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[4]);
    uneNouvelleFace.push_back(sommets[5]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[6]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[7]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);

}

// Courbures

void MainWindow::on_pushButton_courbures_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();
    mesh.request_vertex_colors() ;

    Courbures courb(mesh) ;
    resetAllColorsAndThickness(&mesh);

    courb.set_fixed_colors();
    displayMesh(&mesh, false, -1, ColorBy::CBVertex);

}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(255, 0, 0));
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

MyMesh::Point MainWindow::getFaceNormal(MyMesh* _mesh,VertexHandle v0, VertexHandle v1, VertexHandle v2)
{
     MyMesh::Point resultat;
     QVector<float> vecteur_U, vecteur_V, produitVect, normal;

     //on créer les vecteurs vecteur_U :v0v1 et vecteur_V : v0v2
     for(int i=0; i<3; i++)
     {
         vecteur_U.push_back(_mesh->point(v1)[i] -_mesh->point(v0)[i]);
         vecteur_V.push_back(_mesh->point(v2)[i] -_mesh->point(v0)[i]);
     }

     //On fait le produit vectoriel de U et V

     produitVect.push_back(vecteur_U[1]*vecteur_V[2] - vecteur_U[2]*vecteur_V[1]);
     produitVect.push_back(vecteur_U[0]*vecteur_V[2] - vecteur_U[2]*vecteur_V[0]);
     produitVect.push_back(vecteur_U[0]*vecteur_V[1] - vecteur_U[1]*vecteur_V[0]);

     //On calcul sa norme:
     float normeVect = sqrt(produitVect[0]*produitVect[0]+produitVect[1]*produitVect[1]+produitVect[2]*produitVect[2]);

     //On determine la normal
     for(int i=0; i<3; i++)
     {
         normal.push_back(produitVect[i]/normeVect);
     }

     resultat = MyMesh::Point(normal[0],normal[1],normal[2]);

     return resultat;

}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    qDebug() << __FUNCTION__;

    QVector<MyMesh::Point> * listPoints = new QVector<MyMesh::Point>();
    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(_mesh->face_handle(static_cast<unsigned int>(faceID))); curVert.is_valid(); curVert ++) {
        if (*curVert != _mesh->vertex_handle(static_cast<unsigned int>(vertexID))) {
            listPoints->push_back(_mesh->point(*curVert));
        }
    }

    MyMesh::Point originPoint = _mesh->point(_mesh->vertex_handle(static_cast<unsigned int>(vertexID)));
    //qDebug() << __FUNCTION__ << "norm :" << static_cast<double>(((listPoints->at(0) - originPoint) | (listPoints->at(1) - originPoint)).norm());

    // float angle = static_cast<float>( acos( static_cast<double>(((listPoints->at(0) - originPoint) | (listPoints->at(1) - originPoint)).norm()) ) );

    float norm_v1 = (listPoints->at(0) - originPoint).norm();
    float norm_v2 = (listPoints->at(1) - originPoint).norm();

    float angle = static_cast<float>( acos( ((listPoints->at(0) - originPoint)/norm_v1) | ((listPoints->at(1) - originPoint)/norm_v2) ) );


    //qDebug() << __FUNCTION__ << "angle :" << angle;
    listPoints->clear();
    delete listPoints;

    return angle;
}

float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    qDebug() << __FUNCTION__;

    QVector<MyMesh::Point> * listPoints = new QVector<MyMesh::Point>();
    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(_mesh->face_handle(static_cast<unsigned int>(faceID))); curVert.is_valid(); curVert ++) {
        //qDebug() << "    vertID :" << (*curVert).idx();
        listPoints->push_back(_mesh->point(*curVert));
    }

    float aire = ( ((listPoints->at(1) - listPoints->at(0)) % (listPoints->at(2) - listPoints->at(0))) / 2 ).norm();

    listPoints->clear();
    delete listPoints;

    return aire;
}

float MainWindow::barycentricArea(MyMesh* _mesh, int vertexID)
{
    float barycentric_area = 0.0;

    for(MyMesh::VertexFaceIter curFace = _mesh->vf_iter(_mesh->vertex_handle(static_cast<unsigned int>(vertexID))); curFace.is_valid(); curFace++)
        barycentric_area += faceArea(_mesh, curFace->idx());

    return barycentric_area/3;
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    float somme_teta = 0.0;
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        somme_teta = 0;
        for (MyMesh::VertexFaceIter curFace = _mesh->vf_iter(*curVert); curFace.is_valid(); curFace++) {
            somme_teta += angleEE (_mesh, curVert->idx(), curFace->idx());
        }
        qDebug() << "somme_teta : " << somme_teta;
        float aire_barycentric = barycentricArea(_mesh, curVert->idx());
        qDebug() << "aire barycentric : " << aire_barycentric;
        float courbure_gausienne = (static_cast<float>(2*M_PI) - somme_teta) / aire_barycentric;
        qDebug() << "courbure_gausienne : " << courbure_gausienne;
        mesh.data(*curVert).value = courbure_gausienne;
    }
}

double MainWindow::getNeighboringFacesNormale_Angle(MyMesh* _mesh, MyMesh::Point p1, MyMesh::Point p2)
{

    float a = p1[0]*p2[0];
    float b = p1[1]*p2[1];
    float c = p1[2]*p2[2];
    float prodScal = a+b+c;
    abs(prodScal); //garantie que le produit scalaire soit positif

    double angle = acos(prodScal)/**180/M_PI*/;
    return angle; //en radians
}

double MainWindow::getLenghtPoint1Point2(MyMesh *_mesh,MyMesh::Point p1, MyMesh::Point p2)
{
    //Euler
    double a,b,c,edgeLenght;
    a = p2[0] - p1[0];
    b = p2[1] - p1[1];
    c = p2[2] - p1[2];

    edgeLenght = sqrt(a*a+b*b+c*c);
    return edgeLenght;
}

void MainWindow::getPointSharedByFaces(MyMesh *_mesh,FaceHandle f1, FaceHandle f2,MyMesh::Point *p1,MyMesh::Point *p2)
{
    std::vector<int> faceVertices_idx;
    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(f1); curVert.is_valid(); curVert++)
    {
        faceVertices_idx.push_back(curVert->idx());
    }
    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(f2); curVert.is_valid(); curVert++)
    {
        faceVertices_idx.push_back(curVert->idx());
    }

    int element1,element2;
    std::vector<int> communedge_idx;

    //On ne garde que les points commun aux deux faces
    for(unsigned int i=0; i<faceVertices_idx.size(); i++){
        element1 = faceVertices_idx[i];
        for(unsigned int j=0; j<faceVertices_idx.size(); j++){
            element2 = faceVertices_idx[j];
            if((i != j) &&(element1 == element2)){
                communedge_idx.push_back(element1);
            }
        }
    }

    //suppression des doublons
    std::sort(communedge_idx.begin(), communedge_idx.end());
    communedge_idx.erase(std::unique(communedge_idx.begin(), communedge_idx.end()), communedge_idx.end());

    VertexHandle v0,v1;
    if(communedge_idx.size() > 1){
        v0= _mesh->vertex_handle(communedge_idx[0]);
        v1= _mesh->vertex_handle(communedge_idx[1]);
    }
    *p1 = _mesh->point(v0);
    *p2 = _mesh->point(v1);

}

void MainWindow::H_Curv(MyMesh* _mesh)
{

    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        std::vector<MyMesh::Point> vertexNormals;
        std::vector<FaceHandle> vertexFaces;
        qDebug()<<"=========== sommet courant : ============"<<curVert->idx();

        float somme_teta = 0.0;
        for (MyMesh::VertexFaceIter curFace = _mesh->vf_iter(*curVert); curFace.is_valid(); curFace++) {
            std::vector<VertexHandle> opositeVertex;
            for (MyMesh::FaceVertexIter curVert2 = _mesh->fv_iter(*curFace); curVert2.is_valid(); curVert2++) {
                if(*curVert2 != *curVert)
                    opositeVertex.push_back(*curVert2);

            }
            vertexNormals.push_back(getFaceNormal(&mesh,*curVert, opositeVertex[0],opositeVertex[1]));
            vertexFaces.push_back(curFace);
        }

        for(int i=0; i<vertexNormals.size(); i++)
        {
            MyMesh::Point n1,n2;
            FaceHandle f1,f2;

            if(i==vertexNormals.size()-1)
            {
                n1 = vertexNormals[i];
                n2 = vertexNormals[0];
                f1 = vertexFaces[i];
                f2 = vertexFaces[0];
            }
            else
            {
                n1 = vertexNormals[i];
                n2 = vertexNormals[i+1];
                f1 = vertexFaces[i];
                f2 = vertexFaces[i+1];
            }
            MyMesh::Point p1_shared,p2_shared;

            getPointSharedByFaces(&mesh,f1,f2, &p1_shared,&p2_shared);
            float length_edge = getLenghtPoint1Point2(&mesh,p1_shared,p2_shared);
            somme_teta += getNeighboringFacesNormale_Angle(&mesh,n1,n2)*length_edge;
            qDebug()<<"angle :"<<getNeighboringFacesNormale_Angle(&mesh,n1,n2);
            qDebug()<<"longueur arête"<<length_edge;
        }

        float aire_barycentric = barycentricArea(_mesh, curVert->idx());

        float courbure_moyenne = (somme_teta) / (4*aire_barycentric);

        mesh.data(*curVert).value = courbure_moyenne;
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange, ColorBy cb)
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
            if(_mesh->data(*fvIt).value > 0)
            {
                triCols[3*i+0] = 255;
                triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            else
            {
                triCols[3*i+2] = 255;
                triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            triVerts[3*i+0] = _mesh->point(*fvIt)[0];
            triVerts[3*i+1] = _mesh->point(*fvIt)[1];
            triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0)
            {
                triCols[3*i+0] = 255;
                triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            else
            {
                triCols[3*i+2] = 255;
                triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            triVerts[3*i+0] = _mesh->point(*fvIt)[0];
            triVerts[3*i+1] = _mesh->point(*fvIt)[1];
            triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0)
            {
                triCols[3*i+0] = 255;
                triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            else
            {
                triCols[3*i+2] = 255;
                triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            triVerts[3*i+0] = _mesh->point(*fvIt)[0];
            triVerts[3*i+1] = _mesh->point(*fvIt)[1];
            triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else if (cb != ColorBy::CBVertex)
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
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
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
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true);
}
