#include "MPI_GridSlicer.h"
#include "MPI_Topology.h"

namespace mpi
{

// Slices a min gird into sections and outputs grid parameters to a sub-grid associated with the specified rank
// Grid divisions are specified in an array of std::vectors, divisions[0] - divisions along the X axis, [1] - Y, [2] - Z
void GridSlicer::getSubGridParameters(pfc::FP3& out_origin, pfc::Int3& out_size, 
                                      const pfc::FP3& main_origin, const pfc::Int3& main_size, const pfc::FP3& step,
                                      std::vector<int> divisions[3],
                                      int rank, const Topology& topology)
{
    pfc::Int3 node = topology.getNodeSection(rank);
    const pfc::Int3& sections = topology.getSections();

    Topology::LoopType loopType = topology.getLoopType();

    pfc::Int3 l_border_overlap;
    l_border_overlap.x = topology.doesLoopOverX(loopType) ? (-1) : 0;
    l_border_overlap.y = topology.doesLoopOverY(loopType) ? (-1) : 0;
    l_border_overlap.z = topology.doesLoopOverZ(loopType) ? (-1) : 0;

    pfc::Int3 r_border_overlap;
    r_border_overlap.x = topology.doesLoopOverX(loopType) ? (sections.x) : sections.x - 1;
    r_border_overlap.y = topology.doesLoopOverY(loopType) ? (sections.y) : sections.y - 1;
    r_border_overlap.z = topology.doesLoopOverZ(loopType) ? (sections.z) : sections.z - 1;

    pfc::Int3 lb;
    pfc::Int3 rb;

    lb.x = ( (node.x - 1 < 0) ? l_border_overlap.x : divisions[0][node.x - 1] - 1 );
    rb.x = ( (node.x == r_border_overlap.x) ? main_size.x : divisions[0][node.x] + 1 );

    lb.y = ( (node.y - 1 < 0) ? l_border_overlap.y : divisions[1][node.y - 1] - 1 );
    rb.y = ( (node.y == r_border_overlap.y) ? main_size.y : divisions[1][node.y] + 1 );

    lb.z = ( (node.z - 1 < 0) ? l_border_overlap.z : divisions[2][node.z - 1] - 1 );
    rb.z = ( (node.z == r_border_overlap.z) ? main_size.z : divisions[2][node.z] + 1 );

    // Outputs
    out_origin = main_origin + step * lb;
    out_size = rb - lb;
}

}