////////////////////////////////////////////////////////////////////////////////
// Code for ECS 289H Assignment 0
////////////////////////////////////////////////////////////////////////////////
#include "Viewer.hh"

int main(int argc, char *argv[]) {
    // The following two matrices `V` and `F` represent a simple triangle mesh
    // that is the boundary of a (rotation and translation of) the following
    // regular tetrahedron:
    //       3
    //       *       
    //      / \`.    
    //     /   \ `* 2
    //    / _.--\ /  
    //  0*-------* 1 
    // Each *row* of the `V` matrix gives the coordinates of the corresponding vertex.
    Eigen::MatrixX3d V(4, 3);
    V <<  1.268725608389637438, -1.268204617906207199, -1.038364355233982517, // v0
         -0.572727684923425295, -2.019727444873399058, -1.490504660004000392, // v1
          0.419861124573105204, -2.679782825950963066,  0.164566938228722592, // v2
         -0.340943355162440620, -0.787898531890284914,  0.118545121612435245; // v3

    // Each row of the `F` matrix describes a triangle by listing the indices
    // of its corner vertices. This is the "indexed face set" triangle mesh
    // representation that we will learn about later on.
    // Vertices of each face are listed in counter-clockwise order
    // (when the tetrahedron is viewed from outside).
    Eigen::MatrixX3i F(4, 3);
    F << 0, 2, 1,
         0, 3, 2,
         0, 1, 3,
         1, 2, 3;

    // TODO: Compute outgoing edge vectors for the three edges incident vertex 0
    Eigen::Vector3d e1, e2, e3;
    e1 = V.row(1) - V.row(0);
    e2 = V.row(2) - V.row(0);
    e3 = V.row(3) - V.row(0);
    // TODO: Compute the tetrahedron's volume using the cross and dot products
    double vol = 0.0;
     vol = (1.0/6)*(e1.cross(e2)).dot(e3);
    // TODO: Compute the tetrahedron's volume using the determinant
    Eigen::Matrix3d E;
    double vol_det = 0.0;
    E << e1, e2, e3;
    vol_det = ((1.0)/6)*E.determinant();

    // TODO: Compute translation bringing the tetrahedron's barycenter
    // (average coordinate vector) to the origin.
    Eigen::Vector3d t = Eigen::Vector3d::Zero();
    Eigen::Vector3d barycenter = V.colwise().mean();
    t -= barycenter;
    // TODO :Compute a rotation matrix that, when applied to the points in V,
    // orients face 0 parallel to the XZ plane with edge vector `v1 - v0`
    // pointing in the +x direction and ensures the opposite vertex `3` is in
    // the +y (not -z) direction.
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d unit_e1,unit_n;
    
    unit_e1 = e1/e1.norm();
    unit_n = e1.cross(e2).normalized();
    R << unit_e1, unit_n, unit_e1.cross(unit_n);
    R = R.transpose().eval();

    // TODO: Transform all of the points
    V = V;
    V = (V.rowwise()+t.transpose())*R.transpose();

    std::cout << "Volume computed with cross and dot product: " << vol     << std::endl;
    std::cout << "Volume computed with determinant: "           << vol_det << std::endl;
    std::cout << "translation: " << t.transpose() << std::endl << std::endl;

    std::cout << "Rotation: " << std::endl;
    std::cout << R << std::endl << std::endl;

    std::cout << "Transformed vertices: " << std::endl;
    std::cout << V << std::endl << std::endl;

    std::cout << "Transformed barycenter: " << V.colwise().mean() << std::endl << std::endl;

    // Plot the resulting mesh
    Viewer viewer("289H Homework 0 - Eigen Demo");
    viewer.set_mesh(V, F);
    viewer.run();
}
