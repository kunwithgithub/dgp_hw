////////////////////////////////////////////////////////////////////////////////
// Code for ECS 289H Assignment 1
////////////////////////////////////////////////////////////////////////////////
#include "Viewer.hh"
#include <igl/readOBJ.h>
#include "extract_contours.hh"
#include <cctype>
#include <cmath>

using PtRef = Eigen::Ref<const Eigen::Vector3d>;
using ScalarField = std::function<double(PtRef)>;

// TODO: Implement your additional scalar fields here.
std::vector<ScalarField> scalarFields = {{
    [](PtRef p) { return p.norm() - 1.0; },
},{

    [](PtRef p) { return std::pow(p[1],2.0) - std::sin(std::pow(p[0],2.0)); }
},{

    [](PtRef p) { return std::sin(2.0*p[0]+2.0*p[1]) - std::cos(4.0*p[0]*p[1])+1.0; }
}};

// Global variables to be configured by GUI/keyboard
int scalarFieldIdx = 0;
bool useLinearInterpolation = false;
float snapEpsilon = 0.0;

Eigen::Vector3d sampling_offset = Eigen::Vector3d::Zero();

void update(IGLViewer &v) {
    if (scalarFieldIdx >= scalarFields.size()) {
        std::cerr << "Field " << scalarFieldIdx << " does not exist" << std::endl;
        return;
    }

    const auto &V = v.data().V;
    const auto &F = v.data().F;

    Eigen::VectorXd D(V.rows());
    ScalarField shifted_sf = [&](PtRef p) { return scalarFields[scalarFieldIdx](p + sampling_offset); };
    for (int i = 0; i < V.rows(); ++i)
        D[i] = shifted_sf(V.row(i));

    Eigen::MatrixX3d P;
    Eigen::MatrixX2i E;
    extract_contours(V, F, shifted_sf, P, E, useLinearInterpolation, snapEpsilon);

    v.data().set_data(D);
    v.data().set_edges(P, E, Eigen::RowVector3d(1.0, 1.0, 1.0));
}

// Use the WASD keys to shift the scalar field with respect to the grid.
// This can be helpful to debug linear interpolation/snapping.
bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) {
    key = std::toupper(key);
    bool handled = false;
    if (key == 'A') { sampling_offset[0] += 0.01f; handled = true; }
    if (key == 'D') { sampling_offset[0] -= 0.01f; handled = true; }
    if (key == 'W') { sampling_offset[1] -= 0.01f; handled = true; }
    if (key == 'S') { sampling_offset[1] += 0.01f; handled = true; }
    if (handled) update(viewer);
    return handled;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " mesh.obj" << std::endl;
        exit(-1);
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    bool success = igl::readOBJ(argv[1], V, F);
    if (!success) throw std::runtime_error(std::string("Failed to open input file ") + argv[1]);

    Viewer viewer("289H Homework 1 - Contour Extraction");
    viewer.set_mesh(V, F, /* update_base_camera= */ true);
    update(viewer);

    viewer.menu.callback_draw_viewer_menu = [&]() {
        bool changed = false;
        std::string scalarFieldOptions;
        for (size_t i = 0; i < scalarFields.size(); ++i) { scalarFieldOptions += "Field " + std::to_string(i); scalarFieldOptions.push_back('\0'); }
        scalarFieldOptions.push_back('\0');
        changed |= ImGui::Combo("Scalar Field", &scalarFieldIdx, scalarFieldOptions.c_str());
        changed |= ImGui::Checkbox("Linear Interpolation", &useLinearInterpolation);
        changed |= ImGui::DragFloat("Snapping Threshold", &snapEpsilon, /* vspeed */ 0.05f, /* vmin */ 0.0f, /* vmax */ 0.25f);
        if (changed) update(viewer);
    };

    viewer.callback_key_down    = callback_key_pressed;
    viewer.callback_key_pressed = callback_key_pressed;

    viewer.run();
}
