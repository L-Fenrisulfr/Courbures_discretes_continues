#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void getPointsEdge(MyMesh *_mesh,EdgeHandle *edge,VertexHandle *vh1,VertexHandle *vh2);
    void getFacesSharedByEdge(MyMesh *_mesh, FaceHandle *Fh0, FaceHandle *Fh1,
                                EdgeHandle *edge);
    void getOppositeVertexOfNeighoorFaces(MyMesh *_mesh, FaceHandle fh0,
                                                     FaceHandle fh1, VertexHandle *vh0, VertexHandle *vh1);
    void edgeSplit(MyMesh *_mesh, EdgeHandle *edge, double maxLength);
    void verticesShift(MyMesh *_mesh);
    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);
    void oneRingNeighbors_Barycenter(MyMesh * _mesh, std::vector<double> *barycenter,
                                                 VertexHandle v);
    void getVertexSurfaceApproximation(MyMesh *_mesh,VertexHandle v,
                                                   std::vector<double> barycenterPoint);
    void setVertexCoordonate(MyMesh *_mesh,VertexHandle v,double x,double y, double z);
    double getAngle(MyMesh* _mesh, int vertexID,  int faceID);
    double getmeshAngleQuality(MyMesh *_mesh);
    double getLenghtEdge(MyMesh *_mesh,EdgeHandle edge);
    double getLenghtPoint1Point2(MyMesh *_mesh,VertexHandle v1,VertexHandle v2);
    void collapsShortEdge(MyMesh *_mesh, double minLength, double maxLength);
    void collapseEdge(MyMesh* _mesh, int edgeID);
    void edgesSplit(MyMesh *_mesh, double maxLength);
    void valenceEgalisation(MyMesh *_mesh);
    int getVertexValence(MyMesh *_mesh, VertexHandle v);

private slots:

    void on_pushButton_chargement_clicked();
    void on_remaillage_incremental_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
