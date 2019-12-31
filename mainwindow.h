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
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

enum class ColorBy { CBFace , CBEdge , CBVertex } ;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1, ColorBy cb = ColorBy::CBFace);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    float faceArea(MyMesh* _mesh, int faceID);
    float barycentricArea(MyMesh* _mesh, int vertexID);
    double getNeighboringFacesNormale_Angle(MyMesh* _mesh, MyMesh::Point p1, MyMesh::Point p2);
    MyMesh::Point getFaceNormal(MyMesh* _mesh,VertexHandle v0, VertexHandle v1, VertexHandle v2);
    void getPointSharedByFaces(MyMesh *_mesh, FaceHandle f1, FaceHandle f2, MyMesh::Point *p1, MyMesh::Point *p2);
    double getLenghtPoint1Point2(MyMesh *_mesh,MyMesh::Point p1, MyMesh::Point p2);
    void K_Curv(MyMesh* _mesh);
    void H_Curv(MyMesh* _mesh);

private slots:
    void on_pushButton_chargement_clicked();

    void on_pushButton_generer_clicked();

    void on_pushButton_courbures_clicked();

    void on_pushButton_K_clicked();

    //void on_pushButton_H_clicked();

    void on_pushButton_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
