#include "SignalMesh.hpp"
using namespace std ;


class CommonNeighbor
{
public:
    int cellIDCommonNeighbor ;
    int nghbrIDCell_CommonNeighbor ;
    int nghbrIDNeighboringCell_CommonNeighbor ;
    int nghbrIDCommonNeighbor_Cell ;
    int nghbrIDCommonNeighbor_NeighboringCell ;
};




class Neighbor
{
public:
    int CellID_Neighbor ;       // stores cell ID number of the neighboring cell
    int NeighborID_Cell ;       // stores the neighbor ID number that points to our primary cell ID in the neighboring cell struct
    vector<vector<double > > nodeDistoNghbr ;
    vector<double> minDistNode ;
    vector< int > indexMinDistNode ;
    vector<int> indexCellNode ;
    vector<int> indexNeighborNode ;
    bool common_edge_ID ;
    vector <double> intfX ;
    vector <double> intfY ;
    vector <double> intersectX ;
    vector <double> intersectY ;
    vector<bool> realIntersect ;
    bool curvedInterface = false ;
    vector< CommonNeighbor > commonNeighbors ;

    
};


class SignalCell
{   public:
    SignalCell() ;
    vector<SignalMesh> meshes ;
    vector <Neighbor> neighbors ;
    int cellID ;
    int layer ;
    vector<double> nodesX ;
    vector<double> nodesY ;
    vector<double> nodesXNew ;
    vector<double> nodesYNew ;
    vector<double> noNeighboringNodesX ;
    vector<double> noNeighboringNodesY ;
    vector<double> verticesX ;
    vector<double> verticesY ;
    vector<double> newVertX ;
    vector<double> newVertY ;
    vector<bool> newVertUpdateStatus ;
    vector<int> connections ;
    bool boundary = false ;
    vector<int> cyclic4 ;
    
    double productionW ;
    double productionC ;
    double productionCk ;
    double productionCkR ;
    double productionPMad ;
    
    vector <double> centroid ;
    vector<double> cntrToCntr ;
    vector <int > nghbrCandidate ;
    vector<vector<vector<double> > > nodeDistoNghbrsCandidate ;
    vector<double> cellConcentration ;
    double cellU ;
    double areaCell ;
    
    void Cal_Centroid() ;
    vector <vector <double> > Cal_NodeToNodeDist (vector<double> , vector<double>) ;
    void Find_NghbrCandidate () ;
    void Find_NghbrProperties () ;
    void DetectCommonEdge () ;
    void NewEdge () ;
    void SortCCW (); //vector<double> &vectX , vector<double> &vectY) ;
    void Cal_Vertices () ;
    void Cal_Connections () ;
    void Add_BoundaryVertice () ;
    void Add_BoundaryVertice2 () ;
    void Add_BoundaryVertice3 () ;
    
    
    void Find_Mesh () ;
    void Self_Diffusion () ;
    void Refine_NoBoundary () ;
    void Refine_NodeXNew () ;
    void FullModel_SelfDiffusion (bool type) ;
    void FullModel_ProductionCell () ;
    void CellLevelConcentration (bool type) ;
    void CellLevelConcentration2(bool type) ;
    
    void FullModel_WingDisc_SelfDiffusion () ;
    
    
    
};





