////////////////////////////////////////////////////////////////////////////////
// extract_contours.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Contour extraction framework code for ECS 289H HW 1.
*/
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
////////////////////////////////////////////////////////////////////////////////
#include "extract_contours.hh"
#include <array>
#include <vector>
#include <map>
#include <iostream>

struct UnorderedPair {
    UnorderedPair(size_t a, size_t b) {
        m_data[0] = std::min(a, b);
        m_data[1] = std::max(a, b);
    }

    // Lexicographic comparison
    bool operator<(const UnorderedPair &other) const {
        if (m_data[0] < other.m_data[0]) return true;
        if (m_data[0] > other.m_data[0]) return false;
        return m_data[1] < other.m_data[1]; // break ties with second entry
    }
private:
    Eigen::Vector2i m_data;
};

// Extract the contours for the scalar field `sf` defined over the triangle
// mesh (V, F).
// The contour points are returned in `P` and the edges connecting them in `E`.
// If `linearInterpolation` is false, the contour points should be generated
// at the edge midpoint. If it is true, it should be placed at the exact
// zero-crossing of the linearly interpolated scalar field along the edge.
// Bonus:
// If `snapEpsilon` > 0, contour points whose relative distance along the edge
// from one of the edge endpoints is less than `snapEpsilon` will be snapped
// to that endpoint. Contour edges that degenerate to zero length due to
// this snapping should be neglected.
void extract_contours(Eigen::Ref<const Eigen::MatrixX3d> V,
                      Eigen::Ref<const Eigen::MatrixX3i> F,
                      const ScalarField &sf,
                      Eigen::MatrixX3d &P,
                      Eigen::MatrixX2i &E,
                      bool linearInterpolation,
                      double snapEpsilon) {
    // Map used to look up the index of the contour points already generated
    // for certain undirected edges of the triangle mesh. (The undirected edge
    // is represented as an unordered pair of indices into the `V` array). This
    // map will be used to ensure the two triangles sharing an edge do not
    // generate two copies of the contour point for that edge.
    // Hint: you can record a point index `ptIdx` to be shared by edges (0, 1)
    // and (1, 0) by doing
    //      pointOnTriEdge[UnorderedPair(0, 1)] = ptIdx;
    // After this, pointOnTriEdge.count(UnorderedPair(0, 1)) should return 1.
    std::map<UnorderedPair, size_t> pointOnTriEdge;
    // Map used to prevent duplication of snapped contour points when `snapEpsilon > 0`
    // (Only needed for the bonus snapping task.)
    std::map<size_t, size_t> pointOnTriVertex;

    std::vector<Eigen::Vector3d> contourPts;
    std::vector<Eigen::Vector2i> contourEdges;

    // TODO: implement your contour extraction algorithm here.
    // You can add a point `p` with `contourPts.push_back(p);` and
    // add an edge connecting the contour points indexed `i0` and `i1`
    // with `contourEdges.emplace_back(i0, i1);`.
    //
    // Suggestion: first ignore `linearInterpolation` (place
    // contour points at the edge midpoint), and don't worry about
    // generating duplicate contour points (ignore `pointOnTriEdge`).
    // Then implement the `linearInterpolation` mode, and finally
    // add code to avoid contour point duplication using `pointOnTriEdge`.

        std::vector<int> in;
        std::vector<int> out;
        for(size_t i = 0; i < F.rows();i++){

            for(size_t j = 0; j < F.row(i).cols(); j++){

                if((sf(V.row(F(i,j)))<=0.0? true:false)){  in.push_back(F(i,j));   }
                else{   in.push_back(F(i,j));   }

            }

            std::vector<int> edge_index;

            if(in.empty()||out.empty()){    continue;   }
            else if(in.size()>out.size()){
                
                for(auto it = in.begin();it!=in.end();it++){
                    int snapped_vertex = -1;

                    if(!pointOnTriEdge.count(UnorderedPair(*it, out[0]))){
                        if(!linearInterpolation){
                            Eigen::Vector3d midpoint( (V.row(*it)[0]+V.row(out[0])[0])/2.0 , (V.row(*it)[1]+V.row(out[0])[1])/2.0 , (V.row(*it)[2]+V.row(out[0])[2])/2.0 );
                            contourPts.push_back(midpoint);
                        }else{
                            double s1,s2,t;
                            s1 = sf(V.row(*it));
                            s2 = sf(V.row(out[0]));
                            t = s1/(s1-s2);
                            Eigen::Vector3d p_t = (1-t)*V.row(*it) + t*V.row(out[0]);
                            if(snapEpsilon>0){
                                if((V.row(*it).transpose()-p_t).norm()<(V.row(out[0]).transpose()-p_t).norm()){
                                    double relative_distance = (V.row(*it).transpose()-p_t).norm();
                                    if(relative_distance<snapEpsilon){
                                        if(!pointOnTriVertex.count(*it)){
                                            contourPts.push_back(V.row(*it));
                                            pointOnTriVertex[*it]=contourPts.size()-1;
                                        }else{
                                            snapped_vertex = *it;
                                        }

                                    }else{
                                        contourPts.push_back(p_t);
                                    }
                                }else{
                                    double relative_distance = (V.row(out[0]).transpose()-p_t).norm();
                                    if(relative_distance<snapEpsilon){
                                        if(!pointOnTriVertex.count(out[0])){
                                            contourPts.push_back(V.row(out[0]));
                                            pointOnTriVertex[out[0]]=contourPts.size()-1;
                                        }else{
                                            snapped_vertex = out[0];
                                        }

                                    }else{
                                        contourPts.push_back(p_t);
                                    }

                                }
                            }else{
                                contourPts.push_back(p_t);
                            }
                            

                        }
                        if(snapped_vertex < 0){
                            edge_index.push_back(contourPts.size()-1);
                            pointOnTriEdge[UnorderedPair(*it, out[0])] = (contourPts.size()-1);
                        }else{
                            edge_index.push_back(pointOnTriVertex[snapped_vertex]);
                            pointOnTriEdge[UnorderedPair(*it, out[0])] = pointOnTriVertex[snapped_vertex];
                        }



                    }else{
                        edge_index.push_back(pointOnTriEdge[UnorderedPair(*it, out[0])]);
                    }


                }
                if(edge_index[0]!=edge_index[1]){
                    contourEdges.emplace_back(edge_index[0], edge_index[1]);
                }

            }else{
                
                for(auto it = out.begin();it!=out.end();it++){
                    int snapped_vertex = -1;

                    if(!pointOnTriEdge.count(UnorderedPair(*it, in[0]))){
                        if(!linearInterpolation){
                            Eigen::Vector3d midpoint( (V.row(*it)[0]+V.row(in[0])[0])/2.0 , (V.row(*it)[1]+V.row(in[0])[1])/2.0 , (V.row(*it)[2]+V.row(in[0])[2])/2.0 );
                            contourPts.push_back(midpoint);
                        }else{
                            double s1,s2,t;
                            s1 = sf(V.row(*it));
                            s2 = sf(V.row(in[0]));
                            t = s1/(s1-s2);
                            Eigen::Vector3d p_t = (1-t)*V.row(*it) + t*V.row(in[0]);
                            if(snapEpsilon>0){
                                if((V.row(*it).transpose()-p_t).norm()<(V.row(in[0]).transpose()-p_t).norm()){
                                    double relative_distance = (V.row(*it).transpose()-p_t).norm();
                                    if(relative_distance<snapEpsilon){
                                        if(!pointOnTriVertex.count(*it)){
                                            contourPts.push_back(V.row(*it));
                                            pointOnTriVertex[*it]=contourPts.size()-1;
                                        }else{
                                            snapped_vertex = *it;
                                        }

                                    }else{
                                        contourPts.push_back(p_t);
                                    }
                                }else{
                                    double relative_distance = (V.row(in[0]).transpose()-p_t).norm();
                                    if(relative_distance<snapEpsilon){
                                        if(!pointOnTriVertex.count(in[0])){
                                            contourPts.push_back(V.row(in[0]));
                                            pointOnTriVertex[in[0]]=contourPts.size()-1;
                                        }else{
                                            snapped_vertex = in[0];
                                        }

                                    }else{
                                        contourPts.push_back(p_t);
                                    }

                                }
                            }else{
                                contourPts.push_back(p_t);
                            }
                            

                        }
                        if(snapped_vertex < 0){
                            edge_index.push_back(contourPts.size()-1);
                            pointOnTriEdge[UnorderedPair(*it, in[0])] = (contourPts.size()-1);
                        }else{
                            edge_index.push_back(pointOnTriVertex[snapped_vertex]);
                            pointOnTriEdge[UnorderedPair(*it, in[0])] = pointOnTriVertex[snapped_vertex];
                        }



                    }else{
                        edge_index.push_back(pointOnTriEdge[UnorderedPair(*it, in[0])]);
                    }


                }

                if(edge_index[0]!=edge_index[1]){
                    contourEdges.emplace_back(edge_index[0], edge_index[1]);
                }

            }
            in.clear();
            out.clear();
            edge_index.clear();

        }
    // Convert to Eigen format
    P.resize(  contourPts.size(), 3);
    E.resize(contourEdges.size(), 2);
    for (size_t i = 0; i < contourPts  .size(); ++i) P.row(i) = contourPts[i];
    for (size_t i = 0; i < contourEdges.size(); ++i) E.row(i) = contourEdges[i];
}
